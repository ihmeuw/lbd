#####################################################################
## Generic parallel script for running MBG models                  ##
#####################################################################
# ---SETUP-------------------------------------------------------------------------------------------------------------

## clear environment
rm(list=ls())

#load external packages
pacman::p_load(assertthat, ccaPP, fasterize, fst, mgsub, randtoolbox, magrittr, mgcViz)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

## if running interactively, set arguments
if (interactive) {
  warning('interactive is set to TRUE - if you did not mean to run MBG interactively then kill the model and set interactive to FALSE in parallel script')

  #REDACTED
  
  ## set output directory
  outputdir <- file.path('#REDACTED',indicator_group,indicator,'output',run_date,'/')
  
  ## load an image of the main environment
  load(paste0('#REDACTED/', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))
  
  ## Set BRT parameters from optimizer sheet
  if (any(grepl('gbm', stacked_fixed_effects))) {
    
    gbm_params <- data.table::fread(paste0(core_repo, indicator_group, '/model/configs/gbm_params.csv'), stringsAsFactors=F)
    gbm_tc <- gbm_params[indi == indicator & region == Regions, gbm_tc]
    gbm_lr <- gbm_params[indi == indicator & region == Regions, gbm_lr]
    gbm_bf <- gbm_params[indi == indicator & region == Regions, gbm_bf]
    gbm_nminobs <- gbm_params[indi == indicator & region == Regions, gbm_nminobs]
    gbm_ntrees <- gbm_params[indi == indicator & region == Regions, gbm_ntrees]
    gbm_cv <- gbm_params[indi == indicator & region == Regions, gbm_cv]
    
  } else {
    
    gbm_cv <- NA
    gbm_tc <- NA
    gbm_bf <- NA
    
  }
  
  ## Override skip options
  skiptoinla <- T
  skipinla <- T
  
} else {
  
  ## otherwise, grab arguments from qsub
  ## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
  reg                      <- as.character(commandArgs()[4])
  age                      <- as.numeric(commandArgs()[5])
  run_date                 <- as.character(commandArgs()[6])
  test                     <- as.numeric(commandArgs()[7])
  holdout                  <- as.numeric(commandArgs()[8])
  indicator                <- as.character(commandArgs()[9])
  indicator_group          <- as.character(commandArgs()[10])
  
  ## make a pathaddin that gets used widely
  pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
  
  ## load an image of the main environment
  load(paste0('#REDACTED/', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))
  
  ## In case anything got overwritten in the load, reload args
  reg                      <- as.character(commandArgs()[4])
  age                      <- as.numeric(commandArgs()[5])
  run_date                 <- as.character(commandArgs()[6])
  test                     <- as.numeric(commandArgs()[7])
  holdout                  <- as.numeric(commandArgs()[8])
  indicator                <- as.character(commandArgs()[9])
  indicator_group          <- as.character(commandArgs()[10])
  
  pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
  outputdir <- file.path('#REDACTED',indicator_group,indicator,'output',run_date,'/')
  
} 

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

# print out session info so we have it on record
sessionInfo()
#***********************************************************************************************************************
 
# ---FUNCTIONS----------------------------------------------------------------------------------------------------------
## Load MBG packages
package_list <- c(t(read.csv(paste0(core_repo, '#REDACTED/package_list.csv'), header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

##Be explicit about namespace conflicts

#use your own diacritics fx, due to inscrutable error
#note: requires mgsub pkg
fix_diacritics <<- function(x) {
  
  require(mgsub)
  
  #first define replacement patterns as a named list
  defs <-
    list('??'='S', '??'='s', '??'='Z', '??'='z', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', 
         '??'='C', '??'='E', '??'='E','??'='E', '??'='E', '??'='I', '??'='I', '??'='I', '??'='I', '??'='N', '??'='O', 
         '??'='O', '??'='O', '??'='O', '??'='O', '??'='O', '??'='U','??'='U', '??'='U', '??'='U', '??'='Y', '??'='B', 
         '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='c','??'='e', '??'='e', '??'='e', 
         '??'='e', '??'='i', '??'='i', '??'='i', '??'='i', '??'='o', '??'='n', '??'='o', '??'='o', '??'='o', '??'='o',
         '??'='o', '??'='o', '??'='u', '??'='u', '??'='u', '??'='y', '??'='y', '??'='b', '??'='y', '??'='Ss')
  
  #then force conversion to UTF-8 and replace with non-diacritic character
  enc2utf8(x) %>% 
    mgsub(., pattern=enc2utf8(names(defs)), replacement = defs) %>% 
    return
  
}
#***********************************************************************************************************************
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

## cores to use
cores_to_use <- ifelse(Sys.getenv("OMP_NUM_THREADS")=='', 6, Sys.getenv("OMP_NUM_THREADS"))
message('Running calculations with, ', cores_to_use, ' cores.')
#***********************************************************************************************************************

# ---PREP DATA----------------------------------------------------------------------------------------------------------
PID <- Sys.getpid()
tic("Entire script") # Start master timer

## Set seed for reproducibility
message('Setting seed for reproducibility')
set.seed(98118)

## skip a large chunk if requested in config
if (as.logical(skiptoinla) == FALSE) {
  
  message('You have chosen to not skip directly to inla.')
  
  ## some set up
  if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
  if (class(z_list)    == "character") z_list    <- eval(parse(text=z_list))
  
  ## Load simple polygon template to model over
  gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
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
                            pathaddin   = pathaddin,
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
      save(df, file=paste0(outputdir,"processed_input_data", pathaddin,".rda"))
    } else {
      load('<<< FILEPATH REDACTED >>>')
      save(df, file=paste0(outputdir,"processed_input_data", pathaddin,".rda")) # save to current rundate dir
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
    message('Prepping for subnational REs')
    
    ## Prep subnat location info
    # subnat_country_to_get <- "IND"
    gaul_list_a1        <- get_adm_codes_subnat(gaul_convert(subnat_country_to_get), 
                                                admin_level = 1, shapefile_version = modeling_shapefile_version)
    location_names <- get_location_code_mapping(modeling_shapefile_version)
    
    ## Get India subnational polygon
    subnat_gaul        <- get_adm_codes_subnat(get_adm0_codes(subnat_country_to_get), admin_level = 1, shapefile_version = modeling_shapefile_version)
    subnat_full_shp    <- readRDS(get_admin_shapefile( admin_level = 1, raking = F, suffix = '.rds', version = modeling_shapefile_version ))
    
    subnat_shapefile <- raster::subset(subnat_full_shp, 
                                       ADM0_NAME == location_names[ihme_lc_id %in% subnat_country_to_get, loc_name])
    simple_polygon_list2 <- load_simple_polygon(gaul_list = subnat_gaul, 
                                                buffer = 1, tolerance = 0.4, 
                                                custom_shapefile = subnat_shapefile)
    subset_shape2        <- simple_polygon_list2[[1]]
    simple_polygon2      <- simple_polygon_list2[[2]]
    
    ## Load list of raster inputs (pop and simple)
    raster_list2        <- build_simple_raster_pop(subset_shape2, field = "ADM1_CODE")
    simple_raster2      <- raster_list2[['simple_raster']]
    
    ## Merge ADM1 names to df
    adm1_subset_lox <- over(SpatialPoints(df[,.(long = longitude, lat = latitude)], 
                                          CRS(proj4string(subnat_shapefile))), subnat_shapefile)
    df[, ADM1_CODE:= as.numeric(adm1_subset_lox$ADM1_CODE)]
    
    ## Assign same names to all missings
    df[is.na(ADM1_CODE), ADM1_CODE:= 0]
    
  }  else {
    subset_shape2        <- NULL
    simple_polygon2      <- NULL
    raster_list2        <- NULL
    simple_raster2      <- NULL
  }
  
  ## if testing, we only keep 1000 or so observations
  if (test){
    test_pct <- as.numeric(test_pct)
    
    message(paste0('Test option was set on and the test_pct argument was found at ',test_pct,'% \n
                   ... keeping only ', round(nrow(df)*(test_pct/100),0),' random rows of data.'))
    df <- df[sample(nrow(df),round(nrow(df)*(test_pct/100),0)),]
    
    message('Also, making it so we only take 100 draws')
    samples <- 100
  }
  
  ## Some built in data checks that cause known problems later on
  if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
    stop('You have binomial data where k > N. Check your data before proceeding')
  if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
    stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')
  
  ## Save distribution of data for this region
  png(paste0(outputdir, reg, '.png'))
  if(indicator_family=='binomial') hist(df[, get(indicator)]/df$N) else hist(df[, get(indicator)])
  dev.off()
  
#***********************************************************************************************************************

# ---PREP COVS----------------------------------------------------------------------------------------------------------
  
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
  
  
#***********************************************************************************************************************

# ---STACK--------------------------------------------------------------------------------------------------------------
  
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
  pdf(sprintf('%s/raw_covariates_%s.pdf',outputdir,pathaddin), height=12, width=12)
  for(covname in names(all_cov_layers)){
    plot(all_cov_layers[[covname]],main=covname,maxpixel=1e8)
  }
  dev.off()
  
  ## Check for data where covariate extraction failed
  rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]] %>% na.omit(., the_covs))
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
    # sourcing this script will run the child stackers:
    source(paste0(my_repo, '/#REDACTED/run_child_stackers.R'))

    ## combine the children models
    the_data  <- cbind(the_data, do.call(cbind, lapply(lapply(child_model_names, 'get'), function(x) x[['dataset']])))
    child_model_objs <- lapply(child_model_names, function(x) get(x) %>% .[[x]]) %>% setNames(., child_model_names) 

    ## return the stacked rasters
    stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                          period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                          child_models     = child_model_objs,
                                          indicator_family = indicator_family,
                                          centre_scale_df  = covs_cs_df)

    ## plot stackers
    pdf(paste0(outputdir, "/stackers/", 'stacker_rasters', pathaddin, '.pdf'))
    for(i in 1:length(stacked_rasters))
      plot(stacked_rasters[[i]],main=names(stacked_rasters[[i]]),maxpixel=ncell(stacked_rasters[[i]]))
    dev.off()
    
    message('Stacking is complete')
  }
  
#***********************************************************************************************************************

# ---PREP FOR MBG--------------------------------------------------------------------------------------------------------
  
  ## set the fixed effects to use in INLA based on config args
  if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- ''
  }
  if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- stacked_fixed_effects
  }
  if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
    all_fixed_effects <- fixed_effects ## from config
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
  
  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d           = df,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset)
  saveRDS(mesh_s, file=file.path(outputdir, paste0(reg, '_mesh_s.RDS')))
  
  ## Build temporal mesh (standard for now)
  mesh_t <- build_time_mesh(periods=eval(parse(text=mesh_t_knots)))
  saveRDS(mesh_t, file=file.path(outputdir, paste0(reg, '_mesh_t.RDS')))
  
  ## ## For raw covs, don't want to center-scale (as that will happen in `build_mbg_data_stack()`)
  ##
  ## ## This is a bit weird, but when stacking covs are used the oos-stackers (used in `fit_mbg()`)
  ## ## do not get center scaled in save_mbg_input() - this just harmonizes the measures.  If this
  ## ## step isn't done, then the covs get double-center-scaled with odd results.
  ##
  ## ## For predict_mbg, non-center-scaled covs are pulled from cov_list (either stackes or raw) and
  ## ## center-scaled within the function.  So both fit and predict take place on center-scaled covs

  if (as.logical(use_raw_covs) == TRUE) {
    centre_scale_covs <- FALSE
  } else {
    centre_scale_covs <- TRUE
  }
  
  #produce covariate interrogation plots for 5 random pixels
  #sourced from cov_interrogation_functions.R
  covInterrogation(pixel_count=5) 

  ## Save all inputs for MBG model into correct location on /share
  cov_list <- lapply(cov_list, readAll)
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 simple_raster2    = simple_raster2,
                 mesh_s            = mesh_s,
                 mesh_t            = mesh_t,
                 cov_list          = cov_list,
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = centre_scale_covs)
  
} else { ## END !SKIPTOINLA
  message(paste0('You have chosen to skip directly to INLA. Picking up objects from run_date ',skiptoinla_from_rundate))
  message('Now copying saved MBG inputs from that chosen run_date.')
  
  file.copy(from = paste0('#REDACTED/', indicator_group, '/', indicator, '/model_image_history/', skiptoinla_from_rundate, pathaddin, '.RData'),
            to = paste0('#REDACTED/', indicator_group, '/', indicator, '/model_image_history/', run_date, pathaddin, '.RData'))
  
  ## also need to load the simple_raster2 for subnational RE models
  message('Prepping for subnational REs')
  ## Prep subnat location info
  if (as.logical(use_subnat_res)) {
    # subnat_country_to_get <- "IND"
    gaul_list_a1        <- get_adm_codes_subnat(gaul_convert(subnat_country_to_get), 
                                                admin_level = 1, shapefile_version = modeling_shapefile_version)
    location_names <- get_location_code_mapping(modeling_shapefile_version)
    
    ## Get India subnational polygon
    subnat_gaul        <- get_adm_codes_subnat(get_adm0_codes(subnat_country_to_get), admin_level = 1, shapefile_version = modeling_shapefile_version)
    subnat_full_shp    <- readRDS(get_admin_shapefile( admin_level = 1, raking = F, suffix = '.rds', version = modeling_shapefile_version ))
    
    subnat_shapefile <- raster::subset(subnat_full_shp, 
                                       ADM0_NAME == location_names[ihme_lc_id %in% subnat_country_to_get, loc_name])
    simple_polygon_list2 <- load_simple_polygon(gaul_list = subnat_gaul, 
                                                buffer = 1, tolerance = 0.4, 
                                                custom_shapefile = subnat_shapefile)
    subset_shape2        <- simple_polygon_list2[[1]]
    simple_polygon2      <- simple_polygon_list2[[2]]
    
    ## Load list of raster inputs (pop and simple)
    raster_list2        <- build_simple_raster_pop(subset_shape2, field = "ADM1_CODE")
    simple_raster2      <- raster_list2[['simple_raster']]
    
  }  else {
    simple_raster2      <- NULL
  }
}

## reload data an prepare for MBG
load(paste0('#REDACTED/', indicator_group, '/', indicator, '/model_image_history/', run_date, pathaddin, '.RData'))

## convert stackers to transform space, if desired
## NOTE: we do this here to ensure that the stacker rasters are saved in prevalence/untransformed space
## this is useful for diagnostics and other code that was built expecting the untransformed rasters
if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)){
  message('Converting stackers to logit space')
  if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
  
  ## transform the rasters (and offset extreme values)
  rclmat <- matrix(c(1,Inf,0.9999999, -Inf,0,0.0000001), ncol = 3, byrow = TRUE)
  for (i in child_model_names) {
    cov_list[[i]] <- reclassify(cov_list[[i]], rcl = rclmat, include.lowest = TRUE, right = TRUE)
    cov_list[[i]] <- logit(cov_list[[i]])
    names(cov_list[[i]]) <- paste0(i, '.', 1:length(year_list))
  }
  
  ## transform the stacker values that are in df (and offset extreme values)
  stacker_cols <- grep(paste0("(", paste(child_model_names, collapse="|"), ")(.*_pred)"), names(df), value=T)
  for (i in stacker_cols) {
    df[get(i) <= 0, (i) := 0.0000001]
    df[get(i) >= 1, (i) := 0.9999999]
  }
  df[, (stacker_cols) := lapply(.SD, logit), .SDcols = stacker_cols]
  
}

#***********************************************************************************************************************

# ---RUN MBG------------------------------------------------------------------------------------------------------------
tic("MBG - all") ## Start MBG master timer

## for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
if(as.logical(use_stacking_covs)){
  df[, paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]
}

## Generate MBG formula for INLA call (will run but not used by TMB)
mbg_formula <- build_mbg_formula_with_priors(fixed_effects = all_fixed_effects,
                                             add_nugget    = use_inla_nugget,
                                             nugget_prior  = nugget_prior,
                                             add_ctry_res  = use_inla_country_res,
                                             ctry_re_prior = ctry_re_prior,
                                             ctry_re_sum0  = ctry_re_sum0,
                                             temporal_model_theta1_prior = rho_prior,
                                             temporal_model_theta_prior = theta_prior,
                                             no_gp         = !as.logical(use_gp),
                                             coefs.sum1    = coefs_sum1,
                                             subnat_RE     = as.logical(use_subnat_res),
                                             subnat_re_prior = subnat_re_prior,
                                             use_space_only_gp = as.logical(use_space_only_gp),
                                             use_time_only_gmrf = as.logical(use_time_only_gmrf),
                                             time_only_gmrf_type = time_only_gmrf_type,
                                             nid_RE        = use_nid_res)


## If needed, add fake data to make sure INLA estimates all years
missing_years <- setdiff(year_list, df$year)
if (length(missing_years) > 0 & !as.logical(fit_with_tmb)) {
  fake_data <- df[1:length(missing_years), ]
  fake_data[, year := missing_years]
  fake_data[, c(indicator, 'N', 'weight') := 0]
  fake_data[, period := NULL]
  fake_data <- merge(fake_data, period_map)
  df <- rbind(df, fake_data)
}

## Create SPDE INLA stack

time_fe=FALSE 
input_data <- build_mbg_data_stack(df            = df, # note that merge (if using TMB) will return data in a different (but internally consistent) order, just different than df
                                   fixed_effects = all_fixed_effects,
                                   mesh_s        = mesh_s,
                                   mesh_t        = mesh_t, # not currently implemented with tmb
                                   spde_prior    = spde_prior,
                                   use_ctry_res  = use_inla_country_res,
                                   use_subnat_res  = as.logical(use_subnat_res),
                                   use_gp = as.logical(use_gp),
                                   st_gp_int_zero = as.logical(st_gp_int_zero), 
                                   use_space_only_gp = as.logical(use_space_only_gp),
                                   s_gp_int_zero = as.logical(s_gp_int_zero),
                                   use_time_only_gmrf = as.logical(use_time_only_gmrf),
                                   remove_non_subnats = TRUE,
                                   use_nugget    = use_inla_nugget, # implemented with tmb
                                   exclude_cs    = child_model_names, # raw covs will get center scaled here though (see notes above)
                                   coefs.sum1    = coefs_sum1, # not currenlty implemented tmb
                                   tmb           = fit_with_tmb,
                                   scale_gaussian_variance_N = scale_gaussian_variance_N,
                                   shapefile_version = modeling_shapefile_version, 
                                   zl            = z_list, # if this is not zero and tmb==TRUE, it will trigger 3rd kronecker and fixed effects
                                   zcol          = zcol,   # must not be null if z_list is present
                                   nid_RE        = use_nid_res)

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
                         blas_cores       = cores_to_use,
                         wgts             = weights,
                         intstrat         = intstrat,
                         verbose_output   = TRUE,
                         fe_sd_prior      = 1 / 9, ## this actually sets precision!. prec=1/9 -> sd=3
                         pardiso_license  = '/#REDACTED',
                         omp_strat        = 'pardiso.parallel'
                         ) 
  }
  
  if(fit_with_tmb == TRUE) {
    message('Fitting model with TMB')
    message(sprintf('%i Data points and %i mesh nodes',nrow(df),length(input_data$Parameters$Epsilon_stz)))
    
    # save RDS file of input data for replication
    saveRDS(object = input_data, ## save this here in case predict dies
            file = sprintf('#REDACTED',
                           indicator_group, indicator, run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
    # run the model
    system.time(
      model_fit <- fit_mbg_tmb( lbdcorerepo     = core_repo,
                                cpp_template    = 'mbg_tmb_model',
                                tmb_input_stack = input_data,
                                control_list    = list(trace=1, eval.max=500, iter.max=300, abs.tol=1e-20),
                                optimizer       = 'nlminb', # add optimx
                                ADmap_list      =  NULL ) #list(log_gauss_sigma = factor(NA))  ) # list(  zrho = 0.80))
    )
    
    # clamping
    clamp_covs <- TRUE 
    
    
  }
  
  saveRDS(object = model_fit, ## save this here in case predict dies
          file = sprintf('#REDACTED',
                         indicator_group, indicator, run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
  
  
}else{
  ## skipped fitting INLA so just load model and move to predict
  model_fit <- readRDS( file = sprintf('#REDACTED',
                                       indicator_group, indicator, run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
  
}

toc(log = T) ## End MBG - model fit timer

#***********************************************************************************************************************

# ---PREDICT MBG--------------------------------------------------------------------------------------------------------

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
                shapefile_version = modeling_shapefile_version,
                yl            = year_list,
                use_space_only_gp = as.logical(use_space_only_gp),
                use_time_only_gmrf = as.logical(use_time_only_gmrf),
                simple_raster_subnats = simple_raster2,
                subnat_country_to_get = subnat_country_to_get
    )[[3]]
  } else {
    predict_mbg_tmb(samples              = samp,
                    seed                 = NULL,
                    tmb_input_stack      = input_data,
                    model_fit_object     = model_fit,
                    fes                  = all_fixed_effects, 
                    sr                   = simple_raster,
                    yl                   = year_list,
                    zl                   = z_list,
                    use_space_only_gp    = as.logical(use_space_only_gp),
                    covs_list            = cov_list,
                    clamp_covs           = clamp_covs)
  }
})



#***********************************************************************************************************************

# ---FINISH-------------------------------------------------------------------------------------------------------------
## Make cell preds and a mean, ci raster, and cfb variability raster
cell_pred <- do.call(cbind, pm)
mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))
cirange_ras  <- insertRaster(simple_raster,matrix(apply(cell_pred, 1, cirange),ncol = max(period_map$period)))
lower_ras <- insertRaster(simple_raster,matrix(apply(cell_pred, 1, lower),ncol = max(period_map$period)))
upper_ras <- insertRaster(simple_raster,matrix(apply(cell_pred, 1, upper),ncol = max(period_map$period)))
cfb_ras <- insertRaster(simple_raster,matrix(apply(cell_pred, 1, cfb),ncol = max(period_map$period)))
toc(log = T) # Stop MBG - model predict timer

message('Wrapping up')
message('saving preds')
save_mbg_preds(config     = config,
               time_stamp = time_stamp,
               run_date   = run_date,
               mean_ras   = mean_ras,
               sd_ras     = NULL,
               res_fit    = model_fit,
               cell_pred  = cell_pred,
               df         = df,
               pathaddin  = pathaddin)

# save lower
message('writing upper/lower')
writeRaster(
  lower_ras,
  file = paste0(outputdir, '/', indicator,'_lower_eb', pathaddin),
  overwrite = TRUE
)
#save upper
writeRaster(
  upper_ras,
  file = paste0(outputdir, '/', indicator,'_upper_eb', pathaddin),
  overwrite = TRUE
)

# plot the mean, cirange, and cfb variability rasters in a pdf
message('plotting summmary rasters')
dir.create(paste0(outputdir,'results_maps/'))
pdf(paste0(outputdir,'results_maps/summary_rastersXX', pathaddin, '.pdf'))
plot(mean_ras, main=paste0('mean.', names(mean_ras)), maxpixel=ncell(mean_ras))
plot(lower_ras, main=paste0('lower.', names(lower_ras)), maxpixel=ncell(lower_ras))
plot(upper_ras, main=paste0('upper.', names(upper_ras)), maxpixel=ncell(upper_ras))
plot(cirange_ras, main=paste0('ci.', names(cirange_ras)), maxpixel=ncell(cirange_ras))
plot(cfb_ras, main=paste0('cfb.', names(cfb_ras)), maxpixel=ncell(cfb_ras))
dev.off()

# plot the mean rasters as png
message('plotting mean')
png(paste0(outputdir,'results_maps/mean_rasterXX', pathaddin, '.png'), 1800, 1200)
print(spplot(mean_ras, main=list(label=paste0(toupper(indicator), ': mean'),cex=2)))
dev.off()

## timer stuff
message('final mbg exit steps')
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
output_file <- paste0(outputdir, "run_summary_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}

# Write a an empty file to indicate done with this parallel script
write(NULL, file = paste0(outputdir, "/fin_", pathaddin))

## Write CSV
write.csv(df_timer,file = output_file, row.names = FALSE)

#***********************************************************************************************************************

# ---LAUNCH_AGG---------------------------------------------------------------------------------------------------------
#save necessary objects


#dont launch aggregation if we are modelling kerosene
if(indicator%in%c('cooking_fuel_solid', 'cooking_fuel_dirty')) {

  # set specific arguments
  measure         <- 'prev'
  jname           <- paste('eDL', reg, indicator, sep = '_')
  
  # set memory by region
  individual_countries <- ifelse(nchar(reg) == 3, TRUE, FALSE)
  if (reg %in% c('dia_chn_mng', 'dia_s_america-GUY', 'dia_s_america-BRA')) { mymem <- '900G'
  } else if (reg %in% c('wssa-CPV-NGA', 'trsa-GUF', 'CHN', 'soas', 'ansa-VEN', 'ocea-MYS')) { mymem <- '500G'
  } else mymem <- '200G'
  
  # set up qsub
  sys.sub <- paste0('qsub -e ', outputdir, '/errors -o ', outputdir, '/output ', 
                    '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                    '-l fthread=2 -l h_rt=00:24:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
  r_shell <- file.path(core_repo, 'mbg_central/share_scripts/shell_sing.sh')
  script <- file.path(my_repo, indicator_group, 'post/2_entry.R')
  args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, run_date, measure, holdout, my_repo)
                      
  # submit qsub
  paste(sys.sub, r_shell, script, args) %>% 
    system
  
}
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~