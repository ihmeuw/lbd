#####################################################################
## Generic parallel script for running MBG models                  ##
## Roy Burstein, Nick Graetz, Aaron Osgood-Zimmerman, Jon Mosser   ##
#####################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## grab arguments
## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

# check for existence of core repo
# if it was not passed, default to cwd
core_repo <- commandArgs()[11]
if(is.na(core_repo)) {
  core_repo <- getwd()
  message(sprintf("Core repo was not passed, defaulting to current directory: %s.", core_repo))
  if(!file.exists(file.path(core_repo, 'mbg_central', 'setup.R'))) {
    stop("Could not locate mbg setup. You must either pass the core repo as an
          argument or launch this script from the `lbd_core` directory. For
          more assistance, contact the Core Code team.")
  }
}

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
commondir <- paste(core_repo, 'mbg_central/share_scripts/common_inputs', sep = '/')
package_list <- c(t(read.csv(paste(commondir, 'package_list.csv', sep = '/'), header = FALSE)))
mbg_setup(package_list = package_list, repos = core_repo)

# save package list
old_package_list <- package_list

# get path builder for pre run image
rpb <- RunPathBuilder$from_globals(c('indicator',
                                    'indicator_group',
                                    'run_date',
                                    'reg',
                                    'age',
                                    'holdout'))

## load an image of the main environment
load(rpb$get_pre_run_image_file())

# keep path builder
pb <- rpb

# mbg_setup again with package_list from image
if(!setequal(old_package_list, package_list)) {
  mbg_setup(package_list = package_list, repos = core_repo)
}

# get arguments after loading pre run image in case they conflict
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

pathaddin <- pb$get_add_in()
outputdir <- pb$get_output_dir()
dir.create(outputdir, showWarnings = FALSE)

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

# print out session info so we have it on record
sessionInfo()

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
check_config()

# We need to be in the singularity image, and specifically the LBD one if using TMB
if(!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if(as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Print the core_repo hash and check it
message("Printing git hash for 'core_repo' and checking against LBD Core Code master repo")
record_git_status(core_repo = core_repo, check_core_repo = TRUE)

## Make sure this inla patch is implemented if running on geos
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

## cores to use
cores_to_use <- Sys.getenv("SGE_HGR_fthread")
message(paste("Model set to use", cores_to_use, "cores"))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PID <- Sys.getpid()
tic("Entire script") # Start master timer

## skip a large chunk if requested in config
if(as.logical(skiptoinla) == FALSE){
  message('You have chosen to not skip directly to inla.')

  ## some set up
  if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
  if (class(z_list)    == "character") z_list    <- eval(parse(text=z_list))

  ## Load simple polygon template to model over
  gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
                                             shapefile_version = modeling_shapefile_version)
  subset_shape        <- simple_polygon_list[[1]]
  simple_polygon      <- simple_polygon_list[[2]]

  ## Load list of raster inputs (pop and simple)
  raster_list        <- build_simple_raster_pop(subset_shape)
  simple_raster      <- raster_list[['simple_raster']]
  pop_raster         <- raster_list[['pop_raster']]

  ## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
  ## For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
  if(holdout==0) {
    message('Holdout == 0 so loading in full dataset using load_input_data()')

    if(!exists("loadinputdata_from_rundate")) {
      df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                            agebin      = age,
                            removeyemen = FALSE,
                            pathaddin   = pb$get_add_in(),
                            years       = yearload,
                            withtag     = as.logical(withtag),
                            datatag     = datatag,
                            use_share   = as.logical(use_share),
                            yl          = year_list,
                            region      = reg)

      ## Create df with pixel ids (and polygon centroids and age group weights if aggregation)
      df <- process_input_data(df,
                               pop_release                = pop_release,
                               interval_mo                = interval_mo,
                               modeling_shapefile_version = modeling_shapefile_version,
                               poly_ag                    = as.logical(poly_ag),
                               zcol                       = zcol,
                               zcol_ag                    = zcol_ag,
                               zcol_ag_id                 = zcol_ag_id,
                               z_map                      = z_map,
                               z_ag_mat                   = z_ag_mat)
      save(df, file = pb$get_processed_input_file())
    } else {
      load(pb$with_updated_nodes(run_date = loadinputdata_from_rundate)
           $get_processed_input_file())
      save(df, file = pb$get_processed_input_file())
    }
  } else {
    message(paste0('Holdout != 0 so loading holdout data only from holdout ',holdout))
    message('Please be sure you have a list object called stratum_ho in your environment.')
    ## if stratifies by age then make sure loads correctly
    if(age!=0) df <- as.data.table(stratum_ho[[paste('region',reg,'_age',age,sep='__')]])
    if(age==0) df <- as.data.table(stratum_ho[[paste('region',reg,sep='__')]])
    df <- df[fold != holdout, ]
    df$first_entry <- 1
    df$agg_weight <- 1
  }

  if(as.logical(use_subnat_res)) {

    #get adm0 codes for countries to get subnat REs
    if("all" %in% subnat_country_to_get){
      countries_to_get_subnat_res <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
    } else {
      countries_to_get_subnat_res <- get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version)[get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version) %in% get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)]
    }

    #load and subset standard admin1 shape to countries with subnational REs
    subnat_full_shp    <- readRDS(get_admin_shapefile( admin_level = 1, raking = F, suffix = '.rds', version = modeling_shapefile_version ))
    subnat_shapefile <- raster::subset(subnat_full_shp,
                                       ADM0_CODE %in% countries_to_get_subnat_res)

    simple_polygon_list2 <- load_simple_polygon(gaul_list = NULL,
                                                buffer = 1, tolerance = 0.4,
                                                custom_shapefile = subnat_shapefile)
    subset_shape2        <- simple_polygon_list2[[1]]

    ## Load list of raster inputs (pop and simple)
    raster_list2        <- build_simple_raster_pop(subset_shape2, field = "ADM1_CODE")
    simple_raster2      <- raster_list2[['simple_raster']]

    #simple_raster2 has to be the same size as the simple raster for predict to work correctly
    simple_raster2 <- raster::extend(simple_raster2, extent(simple_raster))
    simple_raster2 <- raster::crop(simple_raster2, extent(simple_raster))

    ## Merge ADM0/1 codes to df
    adm1_subset_lox <- over(SpatialPoints(df[,.(long = longitude, lat = latitude)],
                                          CRS(proj4string(subnat_shapefile))), subnat_shapefile)
    df[, subnat_re_ADM1_CODE := as.numeric(as.character(adm1_subset_lox$ADM1_CODE))]
    df[, subnat_re_ADM0_CODE := as.numeric(as.character(adm1_subset_lox$ADM0_CODE))]

    #create new ADM1 columns for each country in subnat_country_to_get so data can be fit separately
    for(i in 1:length(unique(na.omit(df$subnat_re_ADM0_CODE)))) {
      df[subnat_re_ADM0_CODE == unique(na.omit(df$subnat_re_ADM0_CODE))[i], (paste0("SUBNAT", i)) := subnat_re_ADM1_CODE]
    }
  } else {
    simple_raster2 <- NULL
  }

  ## if testing, we only keep 1000 or so observations
  if(as.logical(test)){
    test_pct <- as.numeric(test_pct)

    message(paste0('Test option was set on and the test_pct argument was found at ',test_pct,'% \n
                 ... keeping only ', round(nrow(df)*(test_pct/100),0),' random rows of data.'))

    set.seed(seed)
    seed <- increment_seed(seed)
    df <- df[sample(nrow(df),round(nrow(df)*(test_pct/100),0)),]

    message('Also, making it so we only take 100 draws')
    samples <- 100
  }

  ## if there is another weight column, multiply it with weight now
  if(exists('other_weight')) if(other_weight!='') {
    message(paste0('Multiplying weight and ',other_weight))
    df[['weight']] <- df[['weight']]*df[[other_weight]]
  }

  ## Some built in data checks that cause known problems later on
  if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
    stop('You have binomial data where k > N. Check your data before proceeding')
  if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
    stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')

  ## Save distribution of data for this region
  png(file.path(pb$get_output_dir(), sprintf('%s.png', reg)))
  if(indicator_family=='binomial') hist(df[df$first_entry==1, get(indicator)]/df$N[df$first_entry==1]) else hist(df[df$first_entry==1, get(indicator)])
  dev.off()

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Define modeling space. In years only for now.
  if(yearload=='annual') period_map <- make_period_map(modeling_periods = c(min(year_list):max(year_list)))
  if(yearload=='five-year') period_map <- make_period_map(modeling_periods = seq(min(year_list), max(year_list), by = 5))

  ## Make placeholders for covariates
  cov_layers <- gbd_cov_layers <- NULL

  ## Pull all covariate bricks/layers
  if (nrow(fixed_effects_config) > 0) {
    message('Grabbing raster covariate layers')
    loader <- MbgStandardCovariateLoader$new(start_year = min(year_list),
                                             end_year = max(year_list),
                                             interval = as.numeric(interval_mo),
                                             covariate_config = fixed_effects_config)
    cov_layers <- loader$get_covariates(simple_polygon)
  }

  ## Pull country level gbd covariates
  if (nchar(gbd_fixed_effects) > 0) {
    message('Grabbing GBD covariates')

    effects <- trim(strsplit(gbd_fixed_effects, "\\+")[[1]])
    measures <- trim(strsplit(gbd_fixed_effects_measures, "\\+")[[1]])
    gbd_cov_layers <- load_gbd_covariates(covs     = effects,
                                          measures = measures,
                                          year_ids = year_list,
                                          age_ids  = gbd_fixed_effects_age,
                                          template = cov_layers[[1]][[1]],
                                          modeling_shapefile_version = modeling_shapefile_version,
                                          simple_polygon = simple_polygon,
                                          interval_mo = interval_mo)
  }

  ## Combine all covariates
  all_cov_layers <- c(cov_layers, gbd_cov_layers)

  ## regenerate all fixed effects equation from the cov layers
  all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")

  ## Make stacker-specific formulas where applicable
  all_fixed_effects_brt <- all_fixed_effects

  ## Set Up Country Fixed Effects
  if(use_child_country_fes == TRUE | use_inla_country_fes == TRUE) {
    message('Setting up country fixed effects')
    fe_gaul_list <- unique(c(gaul_convert(unique(df[, country]),
                                          shapefile_version =
                                            modeling_shapefile_version),
                             gaul_list))
    fe_template  <- cov_layers[[1]][[1]]
    simple_polygon_list <- load_simple_polygon(gaul_list   = fe_gaul_list,
                                               buffer      = 0.4,
                                               subset_only = TRUE,
                                               shapefile_version = modeling_shapefile_version)
    fe_subset_shape     <- simple_polygon_list[[1]]
    gaul_code <- rasterize_check_coverage(fe_subset_shape, fe_template, field = 'ADM0_CODE')
    gaul_code <- setNames(gaul_code,'gaul_code')
    gaul_code <- create_categorical_raster(gaul_code)

    ## update covlayers and add country fixed effects to the
    all_cov_layers = update_cov_layers(all_cov_layers, gaul_code)
    all_fixed_effects_cfes = paste(all_fixed_effects,
                                   paste(names(gaul_code)[1:length(names(gaul_code))],
                                         collapse = " + "), sep=" + ")

    ## update specific stacker formulas (for now we just want country effects in BRT)
    all_fixed_effects_brt <- all_fixed_effects_cfes
  }

  ## Add these to the fixed effects if we want them in stacking
  if(use_child_country_fes == TRUE) {
    gaul_fes <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
    all_fixed_effects = paste(all_fixed_effects, gaul_fes, sep = " + ")
  }


  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Stacking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  tic("Stacking - all") ## Start stacking master timer

  ## Figure out which models we're going to use
  child_model_names <- stacked_fixed_effects        %>%
    gsub(" ", "", .)          %>%
    strsplit(., "+", fixed=T) %>%
    unlist
  message(paste0('Child stackers included are: ',paste(child_model_names,collapse=' // ')))

  the_covs <- format_covariates(all_fixed_effects)

  ## copy the dataset to avoid unintended namespace conflicts
  the_data <- copy(df)

  ## only use data where we know what age group or point
  ag_data <- the_data[the_data$agg_weight!=1, ]
  the_data <- the_data[the_data$agg_weight==1, ]

  ## shuffle the data into six folds
  set.seed(seed)
  seed <- increment_seed(seed)
  the_data <- the_data[sample(nrow(the_data)),]
  the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=as.numeric(n_stack_folds),labels=FALSE)]

  ## add a row id column
  the_data[, a_rowid := seq(1:nrow(the_data))]

  ## extract covariates to the points and subset data where its missing covariate values
  cs_covs <- extract_covariates(the_data,
                                all_cov_layers,
                                id_col              = "a_rowid",
                                return_only_results = TRUE,
                                centre_scale        = TRUE,
                                period_var          = 'year',
                                period_map          = period_map)

  # A check to see if any of the variables do not vary across the data. This could break model later so we check and update some objects
  covchecklist <- check_for_cov_issues(check_pixelcount = check_cov_pixelcount,
                                       check_pixelcount_thresh = ifelse(exists("pixelcount_thresh"), as.numeric(pixelcount_thresh), 0.95))
  for(n in names(covchecklist)){
    assign(n, covchecklist[[n]])
  }

  # plot covariates as a simple diagnostic here
  pdf(file.path(pb$get_output_dir(),
                sprintf('raw_covariates_%s.pdf', pb$get_add_in())), height = 12, width = 12)
  for(covname in names(all_cov_layers)){
    plot(all_cov_layers[[covname]],main=covname,maxpixel=1e6)
  }
  dev.off()

  ## Check for data where covariate extraction failed
  rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
  if (rows_missing_covs > 0) {
    pct_missing_covs <- round((rows_missing_covs/nrow(the_data))*100, 2)
    warning(paste0(rows_missing_covs, " out of ", nrow(the_data), " rows of data ",
                   "(", pct_missing_covs, "%) do not have corresponding ",
                   "covariate values and will be dropped from child models..."))
    if (rows_missing_covs/nrow(the_data) > 0.1) {
      stop(paste0("Something has gone quite wrong: more than 10% of your data does not have ",
                  "corresponding covariates.  You should investigate this before proceeding."))
    }
  }

  the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)

  ## store the centre scaling mapping
  covs_cs_df  <-  cs_covs[[2]]

  if(as.logical(use_raw_covs) | as.logical(use_stacking_covs)) {
    ## this will drop rows with NA covariate values
    the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))
  }

  ## stop if this na omit demolished the whole dataset
  if(nrow(the_data) == 0) stop('You have an empty df, make sure one of your covariates was not NA everywhere.')

  ## seperated out into a different script
  if(as.logical(use_stacking_covs)){
    message('Fitting Stackers')

    # Run the child stacker models
    child_model_run <- run_child_stackers(models = child_model_names, input_data = the_data)

    # Bind the list of predictions into a data frame
    child_mods_df <- do.call(cbind, lapply(child_model_run, function(x) x[[1]]))

    ## combine the children models with the_data
    the_data  <- cbind(the_data, child_mods_df)

    ## Rename the child model objects into a named list
    child_model_objs <- setNames(lapply(child_model_run, function(x) x[[2]]), child_model_names)



    ## return the stacked rasters
    stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                          period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                          child_models     = child_model_objs,
                                          indicator_family = indicator_family,
                                          centre_scale_df  = covs_cs_df)

    ## plot stackers
    pdf(file.path(pb$get_output_dir(), sprintf('stacker_rasters%s.pdf', pb$get_add_in())))
    for(i in 1:length(stacked_rasters))
      plot(stacked_rasters[[i]],main=names(stacked_rasters[[i]]),maxpixel=ncell(stacked_rasters[[i]]))
    dev.off()

    message('Stacking is complete')
  } ## if(use_stacking_covs)

  ## add aggregate data back in, with stacking predictions from the full model
  if (nrow(ag_data) > 0) {
    ag_data[, a_rowid := 1:.N + max(the_data$a_rowid)]
    if(as.logical(use_stacking_covs)) {
      ag_stackers <- extract_covariates(ag_data,
                                        stacked_rasters,
                                        id_col              = "a_rowid",
                                        return_only_results = TRUE,
                                        centre_scale        = FALSE,
                                        period_var          = "year",
                                        period_map          = period_map)
      ag_stackers <- ag_stackers[, c("a_rowid", child_model_names, child_model_names), with = F]

      stacker_names <- c(paste0(child_model_names, "_full_pred"), paste0(child_model_names, "_cv_pred"))
      setnames(ag_stackers, c("a_rowid", stacker_names))

      ag_data <- merge(ag_data, ag_stackers, by = "a_rowid")

      if(any(is.na(ag_data[, ..stacker_names]))) {
        stop("There are NAs in predictions from stackers for aggregated data.
               Please contact the core code team if you encounter this problem.")
      }
    } else {
      ag_covs <- extract_covariates(ag_data,
                                    all_cov_layers,
                                    id_col              = "a_rowid",
                                    return_only_results = TRUE,
                                    centre_scale        = FALSE,
                                    period_var          = 'year',
                                    period_map          = period_map)
      cov_names <- names(all_cov_layers)
      cs_ag_covs <- centreScale(ag_covs[, ..cov_names], df = covs_cs_df)
      ag_covs <- cbind(ag_covs[,- ..cov_names], cs_ag_covs)

      ag_data <- merge(ag_data, ag_covs, by = "a_rowid")

      if(as.logical(use_raw_covs)) {
        if(any(is.na(ag_data[, ..cov_names]))) {
          stop("There are NAs in covariates for aggregated data.
               Please contact the core code team if you encounter this problem.")
        }
      }
    }

    the_data <- rbind(the_data, ag_data, fill = T)
  }




  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## set the fixed effects to use in INLA based on config args
  if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- ''
  }
  if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- stacked_fixed_effects
  }
  if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(fixed_effects, gbd_fixed_effects, sep = " + ") ## from config
  }
  if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
  }
  if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, sep = " + ")
  }
  if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }
  if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }
  if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
    all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  }

  ## copy things back over to df
  df <- copy(the_data)

  ## remove the covariate columns so that there are no name conflicts when they get added back in
  df <- df[,paste0(the_covs) := rep(NULL, length(the_covs))]

  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
  df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]

  ## create a full raster list to carry though to the shiny/next steps
  if(as.logical(use_stacking_covs)){
    cov_list      <- c(unlist(stacked_rasters),unlist(all_cov_layers))
    child_mod_ras <- cov_list[child_model_names]
  }else{
    cov_list <- unlist(all_cov_layers)
    child_model_names <- ''
  }

  toc(log = T) ## End stacking master timer

  ## make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

  set.seed(seed)
  seed <- increment_seed(seed)
  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d           = df,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset,
                             s2mesh = as.logical(use_s2_mesh),
                             s2params = s2_mesh_params)

  ## Build temporal mesh (standard for now)
  if (length(unique(year_list)) == 1) {
    mesh_t <- NULL
  } else {
    mesh_t <- build_time_mesh(periods = eval(parse(text = mesh_t_knots)))
  }


  ## ## For raw covs, don't want to center-scale (as that will happen in `build_mbg_data_stack()`)
  ##
  ## ## This is a bit weird, but when stacking covs are used the oos-stackers (used in `fit_mbg()`)
  ## ## do not get center scaled in save_mbg_input() - this just harmonizes the measures.  If this
  ## ## step isn't done, then the covs get double-center-scaled with odd results.
  ##
  ## ## For predict_mbg, non-center-scaled covs are pulled from cov_list (either stackes or raw) and
  ## ## center-scaled within the function.  So both fit and predict take place on center-scaled covs
  ##
  ## ## TODO: move all center-scaling to a single location to avoid these crazy acrobatics.
  ## ## But for now, we just do this
  if (as.logical(use_raw_covs) == TRUE) {
    centre_scale_covs <- FALSE
  } else {
    centre_scale_covs <- TRUE
  }

  ## Save all inputs for MBG model into correct location on /share
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 mesh_s            = mesh_s,
                 mesh_t            = mesh_t,
                 cov_list          = cov_list,
                 pathaddin         = pb$get_add_in(),
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = centre_scale_covs)

} else { ## END !SKIPTOINLA
  message(paste0('You have chosen to skip directly to INLA. Picking up objects from run_date ',skiptoinla_from_rundate))
  message('Now copying saved MBG inputs from that chosen run_date.')

  file.copy(from = pb$with_updated_nodes(run_date = skiptoinla_from_rundate)
            $get_image_file(),
            to = pb$get_image_file())
}

## reload data an prepare for MBG
load(file.path(pb$get_image_file()))

# Bound GBM to 0-1 if desired
if (exists("gbm_bounded_0_1")) {
  if (as.logical(gbm_bounded_0_1) == T) {
    message("Truncating GBM values > 1 to 0.999")
    values(cov_list[["gbm"]])[values(cov_list[["gbm"]]) >= 1 & !is.na(values(cov_list[["gbm"]]))] <- 0.999
    gbm_cols <- grep(paste0("(gbm)(.*_pred)"), names(df), value=T)
    replace_one <- function(x) {
      x[x>=1 & !is.na(x)] <- 0.999
      return(x)
    }
    df[, (gbm_cols) := lapply(.SD, replace_one), .SDcols = gbm_cols]
  }
}

## convert stackers to transform space, if desired
## NOTE: we do this here to ensure that the stacker rasters are saved in prevalence/untransformed space
## this is useful for diagnostics and other code that was built expecting the untransformed rasters
if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)){
  message('Converting stackers to logit space')

  ## transform the rasters
  for (ii in child_model_names) {

    ## Preserve variable names in the raster first
    tmp_rastvar <- names(cov_list[[ii]])

    ## Logit
    cov_list[[ii]] <- logit(cov_list[[ii]])

    ## Reassign names
    names(cov_list[[ii]]) <- tmp_rastvar
    rm(tmp_rastvar)
  }

  ## transform the stacker values that are in df
  stacker_cols <- grep(paste0("(", paste(child_model_names, collapse="|"), ")(.*_pred)"), names(df), value=T)
  df[, (stacker_cols) := lapply(.SD, logit), .SDcols = stacker_cols]

}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


tic("MBG - all") ## Start MBG master timer

## for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
if(as.logical(use_stacking_covs)){
  df[, paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]
}

## Generate MBG formula for INLA call (will run but not used by TMB)
mbg_formula <- build_mbg_formula_with_priors(fixed_effects               = all_fixed_effects,
                                             add_nugget                  = use_inla_nugget,
                                             nugget_prior                = nugget_prior,
                                             add_ctry_res                = use_country_res,
                                             ctry_re_prior               = ctry_re_prior,
                                             temporal_model_theta1_prior = rho_prior,
                                             no_gp                       = !as.logical(use_gp),
                                             use_space_only_gp           = as.logical(use_space_only_gp),
                                             use_time_only_gmrf          = as.logical(use_time_only_gmrf),
                                             time_only_gmrf_type         = time_only_gmrf_type,
                                             coefs.sum1                  = coefs_sum1,
                                             subnat_RE                   = use_subnat_res,
                                             subnat_country_to_get       = subnat_country_to_get,
                                             subnat_re_prior             = subnat_re_prior,
                                             timebycountry_RE            = use_timebyctry_res,
                                             adm0_list                    = gaul_list)

## For INLA we need to add data for missing time points to ensure we get predictions
##  for all relevant time points. The 0 observations do not contribute to the
##  model fitting but they prevent INLA from auto-removing
##  random effects that (conditionally) have no data impacting their fit
## This also seems to improve TMB (??)
  if(use_timebyctry_res) {
    ## If we are using a time only effect by country then we need to make sure
    ##  all year effects are estimated for each country.
    df$adm0code <- gaul_convert(df$country)
    for(adm0_code in gaul_list) {
      dfsub <- df[df$adm0code == adm0_code, ]
      missing_years <- setdiff(year_list, dfsub$year)
      if (length(missing_years) > 0) {
        fake_data <- dfsub[1:length(missing_years), ]
        fake_data[, year := missing_years]
        fake_data[, c(indicator, 'N', 'weight') := 0]
        fake_data[, period := NULL]
        fake_data <- merge(fake_data, period_map)
        df <- rbind(df, fake_data)
      }
    }
  } else {
    ## If not, we only need to make sure we have an observation for each missing
    ##  year (country irrelevant)
    missing_years <- setdiff(year_list, df$year)
    if (length(missing_years) > 0) {
      fake_data <- df[1:length(missing_years), ]
      fake_data[, year := missing_years]
      fake_data[, c(indicator, 'N', 'weight') := 0]
      fake_data[, period := NULL]
      fake_data <- merge(fake_data, period_map)
      df <- rbind(df, fake_data)
    }
  }

# get covariate constraints
cov_constraints <- covariate_constraint_vectorize(fixed_effects,
                                                  gbd_fixed_effects,
                                                  fixed_effects_constraints,
                                                  gbd_fixed_effects_constraints)

## Create SPDE INLA stack
input_data <- build_mbg_data_stack(df            = df, # note that merge (if using TMB) will return data in a different (but internally consistent) order, just different than df
                                   fixed_effects = all_fixed_effects,
                                   mesh_s        = mesh_s,
                                   mesh_t        = mesh_t, # not currently implemented with tmb
                                   spde_prior    = spde_prior,
                                   use_ctry_res  = use_country_res,
                                   use_subnat_res  = use_subnat_res,
                                   use_nid_res   = use_nid_res,
                                   use_gp = as.logical(use_gp),
                                   st_gp_int_zero = as.logical(st_gp_int_zero),
                                   use_space_only_gp = as.logical(use_space_only_gp),
                                   s_gp_int_zero = as.logical(s_gp_int_zero),
                                   use_time_only_gmrf = as.logical(use_time_only_gmrf),
                                   use_timebyctry_res = use_timebyctry_res,
                                   adm0_list = gaul_list,
                                   use_age_only_gmrf = as.logical(use_age_only_gmrf),
                                   use_nugget    = use_inla_nugget, # implemented with tmb
                                   exclude_cs    = child_model_names, # raw covs will get center scaled here though (see notes above)
                                   coefs.sum1    = coefs_sum1, # not currenlty implemented tmb
                                   tmb           = fit_with_tmb,
                                   scale_gaussian_variance_N = scale_gaussian_variance_N,
                                   shapefile_version = modeling_shapefile_version,
                                   zl            = z_list, # if this is not zero and tmb==TRUE, it will trigger 3rd kronecker and fixed effects
                                   zcol          = zcol, # must not be null if z_list is present
                                   cov_constraints = cov_constraints)

## combine all the inputs, other than cs_df these are not used if you are using TMB
stacked_input  <- input_data[[1]]
spde           <- input_data[[2]] ## used for space-time gp
cs_df          <- input_data[[3]]
spde.sp        <- input_data[[4]] ## used for space only (time stationary) gp

## Generate other inputs necessary
outcome <- df[[indicator]] # N+_i - event obs in cluster
N       <- df$N                  # N_i - total obs in cluster
weights <- df$weight

## catch in case there is no weight column
if(is.null(weights)){
  weights = rep(1,nrow(df))
}

tic("MBG - fit model") ## Start MBG - model fit timer

## Set the number of cores to be equal to input;
## If missing, then revert to cores_to_use value
if(Sys.getenv("OMP_NUM_THREADS") != "") {
  setompthreads(Sys.getenv("OMP_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input OpenMP threads.")
  setompthreads(cores_to_use)
}

if(Sys.getenv("MKL_NUM_THREADS") != "") {
  setmklthreads(Sys.getenv("MKL_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input MKL threads.")
  setmklthreads(cores_to_use)
}

## Fit MBG model
if(!as.logical(skipinla)) {
  if(fit_with_tmb == FALSE) {
    message('Fitting model with R-INLA')

    model_fit <- fit_mbg(indicator_family = indicator_family,
                         stack.obs        = stacked_input,
                         spde             = spde,
                         cov              = outcome,
                         N                = N,
                         int_prior_mn     = intercept_prior,
                         f_mbg            = mbg_formula,
                         run_date         = run_date,
                         keep_inla_files  = keep_inla_files,
                         cores            = cores_to_use,
                         wgts             = weights,
                         intstrat         = intstrat,
                         fe_sd_prior      = 1 / 9, ## this actually sets precision!. prec=1/9 -> sd=3
                         sparse_ordering  = as.logical(sparse_ordering))

  } else {
    message('Fitting model with TMB')
    message(sprintf('%i Data points and %i mesh nodes',nrow(df),length(input_data$Parameters$Epsilon_stz)))

    # save RDS file of input data for replication
    saveRDS(object = input_data, ## save this here in case predict dies
            file = file.path(pb$get_output_dir(),
                             sprintf('%s_TMB_data_input_list_%s_holdout_%i_agebin_%i.RDS',
                                     ifelse(fit_with_tmb, 'tmb', 'inla'),
                                     reg,
                                     holdout,
                                     age)))

   # run the model
    system.time(
    model_fit <- fit_mbg_tmb( lbdcorerepo     = core_repo,
                              cpp_template    = 'mbg_tmb_model',
                              tmb_input_stack = input_data,
                              control_list    = list(trace=1, eval.max=500, iter.max=300, abs.tol=1e-20),
                              optimizer       = 'nlminb', # TODO add optimx
                              ADmap_list      = NULL,
                              sparse_ordering = as.logical(sparse_ordering),
                              seed            = seed)
    )
  }

  saveRDS(object = model_fit, ## save this here in case predict dies
          file = file.path(pb$get_output_dir(),
                           sprintf('%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                                   ifelse(fit_with_tmb, 'tmb', 'inla'),
                                   reg,
                                   holdout,
                                   age)))
}else{
  ## skipped fitting INLA so just load model and move to predict
  model_fit <- readRDS(file = file.path(pb$get_output_dir(),
                                        sprintf('%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                                                ifelse(fit_with_tmb, 'tmb', 'inla'),
                                                reg,
                                                holdout,
                                                age)))
}

toc(log = T) ## End MBG - model fit timer

tic("MBG - predict model") ## Start MBG - model predict timer

## Run predict_mbg on chunks of 50 samples (to avoid memory issues)
message('Making predictions in 50 draw chunks.')

max_chunk <- 50
samples   <- as.numeric(samples)

## Create vector of chunk sizes
chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
pm <- lapply(chunks, function(samp) {
  if(fit_with_tmb == FALSE){
    predict_mbg(res_fit       = model_fit,
                cs_df         = cs_df,
                mesh_s        = mesh_s,
                mesh_t        = mesh_t,
                cov_list      = cov_list,
                samples       = samp,
                simple_raster = simple_raster,
                transform     = transform,
                coefs.sum1    = coefs_sum1,
                pred_gp       = as.logical(use_gp),
                yl            = year_list,
                use_space_only_gp = as.logical(use_space_only_gp),
                use_time_only_gmrf = as.logical(use_time_only_gmrf),
                use_timebyctry_res = as.logical(use_timebyctry_res),
                shapefile_version = modeling_shapefile_version,
                simple_raster_subnats = simple_raster2,
                subnat_country_to_get = subnat_country_to_get,
                seed = seed)[[3]]
  } else {
    predict_mbg_tmb(samples              = samp,
                    seed                 = seed,
                    tmb_input_stack      = input_data,
                    model_fit_object     = model_fit,
                    fes                  = all_fixed_effects, # TODO use input_data or model_fit object for this (in case its changed due to checks)
                    sr                   = simple_raster,
                    yl                   = year_list,
                    zl                   = z_list,
                    transform            = transform,
                    covs_list            = cov_list,
                    clamp_covs           = clamp_covs,
                    cov_constraints = cov_constraints,
                    use_full_interacting_effect = as.logical(use_gp),
                    use_space_only_gp = as.logical(use_space_only_gp),
                    use_time_only_gmrf = as.logical(use_time_only_gmrf),
                    use_age_only_gmrf = as.logical(use_age_only_gmrf),
                    use_timebyctry_res = as.logical(use_timebyctry_res),
                    coefs.sum1           = coefs_sum)
  }
})


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# if z dimension has more than one level, then save each z as a different indicator
if(length(z_list) > 1){

  # reorder pm list, right now its z within each chunk. rbind all z's together
  for(z in z_list) # z_list must be integers starting with 1
    if(length(chunks) > 1)
      for(ch in 2:length(chunks))
        pm[[1]][[z]] <- cbind(pm[[1]][[z]], pm[[ch]][[z]])
  pm <- pm[[1]] # pm is now a list of cell_preds by z

  toc(log = T) # Stop MBG - model predict timer

  # loop over z and save as an indicator each one
  orig_indic  <- indicator

  message('Wrapping up')

 for(z in z_list) {
    cptmp <- pm[[z]]

    indicator <- sprintf('%s_%s%i',orig_indic,zcol,z)  # new indicator name
    zpb <- pb$with_updated_nodes(indicator = indicator)
    dir.create(zpb$get_output_dir())
    message(sprintf('New indicator: %s',indicator))

    # make a mean raster
    library(matrixStats)
    mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cptmp),ncol = max(period_map$period)))
    sd_ras    <- insertRaster(simple_raster,matrix(  rowSds(cptmp),ncol = max(period_map$period)))

    # save z specific objects
    writeRaster(
      mean_ras,
      file = file.path(zpb$get_output_dir(),
                       sprintf('%s_prediction_eb%s', indicator, zpb$get_add_in())),
      overwrite = TRUE
    )

    save(
      cptmp,
      file = file.path(zpb$get_output_dir(),
                       sprintf('%s_cell_draws_eb%s.RData', indicator, zpb$get_add_in())),
      compress = TRUE
    )

    pdf(file.path(zpb$get_output_dir(), sprintf('mean_raster%s.pdf', zpb$get_add_in())))
    plot(mean_ras,main='mean',maxpixel=1e6)
    plot(sd_ras,main='sd',maxpixel=1e6)
    dev.off()

    rm(cptmp)
  }

  indicator <- orig_indic

  # save training data
  write.csv(
    df,
    file = file.path(pb$get_output_dir(),
                     sprintf('%s_trainingdata%s.csv', indicator, pb$get_add_in())),
    row.names = FALSE
  )

  message('done saving indicator-specific outputs by z')

}  else { # if no z colums (most peoples cases)


  ## Make cell preds and a mean raster
  cell_pred <- do.call(cbind, pm)
  mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))
  toc(log = T) # Stop MBG - model predict timer



  message('Wrapping up')
  save_mbg_preds(config     = config,
                 time_stamp = time_stamp,
                 run_date   = run_date,
                 mean_ras   = mean_ras,
                 sd_ras     = NULL,
                 res_fit    = model_fit,
                 cell_pred  = cell_pred,
                 df         = df,
                 pathaddin  = pb$get_add_in())


  # plot the mean raster
  pdf(file.path(pb$get_output_dir(),
                sprintf('mean_rasterXX%s.pdf', pb$get_add_in())))
  plot(mean_ras,maxpixel=1e6)
  dev.off()
}

## timer stuff
toc(log = T) # End master timer

## Format timer
ticlog   <- tic.log(format = F)
df_timer <- generate_time_log(ticlog)
df_timer[, region := reg]
df_timer[, holdout := holdout]
setcolorder(df_timer, c("region", "holdout", "step", "time"))

## Pull run time for this run
run_time_all <- df_timer[step == "Entire script", time]

## Write to a run summary csv file in the output directory
output_file <- file.path(pb$get_output_dir(),
                         sprintf('run_summary_%s_%s.csv', indicator, run_date))

## Pull in file contents from other region/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}


# Write a an empty file to indicate done with this parallel script
write(NULL, file = file.path(pb$get_output_dir(), sprintf('fin_%s', pb$get_add_in())))

## Write CSV
write.csv(df_timer, file = output_file, row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~