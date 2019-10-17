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
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

## make a pathaddin that get used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

## load an image of the main environment
load('<<<< FILEPATH REDACTED >>>>')

## In case anything got overwritten in the load, reload args
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
outputdir <- file.path('<<<< FILEPATH REDACTED >>>>')
dir.create(outputdir, showWarnings = FALSE)
# Create a model fit directory, where a small file will be saved to indicate
#  that model FITTING has finished
model_fit_dir <- paste0(outputdir,'/model_fit_progress')
dir.create(model_fit_dir, showWarnings=FALSE)

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

# print out session info so we have it on record
sessionInfo()

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

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
cores_to_use <- round(as.numeric(slots)*.5)

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
  pop_raster         <- raster_list[['pop_raster']] # TODO DELETE LINE

  ## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
  ## For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
  message('Holdout == 0 so loading in full dataset using load_input_data()')
  df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                        simple      = simple_polygon,
                        agebin      = age,
                        removeyemen = FALSE,
                        pathaddin   = pathaddin,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share),
                        yl          = year_list)
  if(holdout!=0) {
    message(paste0('Holdout != 0 so loading holdout data only from holdout ',holdout))
    nids_to_drop <- oos_fold_table[fold==holdout,]$nid
    df <- subset(df,! nid %in% nids_to_drop)
  }
  write.csv(df,sprintf('%s/input_data_%s.csv',outputdir,pathaddin))
  
  ## if testing, we only keep 1000 or so observations
  if(test == 1){
    test_pct <- as.numeric(test_pct)

    message(paste0('Test option was set on and the test_pct argument was found at ',test_pct,'% \n
                 ... keeping only ', round(nrow(df)*(test_pct/100),0),' random rows of data.'))
    df <- df[sample(nrow(df),round(nrow(df)*(test_pct/100),0)),]

    message('Also, making it so we only take 100 draws')
    samples <- 100
  }

  ## for u5m in particular,  make sure indicator is properly named here, wont affect others
  df[[indicator]] <- df[[gsub(paste0('_age',age),'',indicator)]]

  ##### TEST FOR u5m only
  if(indicator=='died'){
  #  message('DOWNWEIGHTING SBH OBSERVATIONS')
  #  df[data_type=='sbh', died := died*0.10]
  #  df[data_type=='sbh', N := N*0.10]
    df <- df[nid != 281548,]
  }
  
  ## if there is another weight column, multiply it with weight now
  # NOTE: For u5m sbh_weight is already incorporated into N and D by this point
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
  png(paste0(outputdir, reg, '.png'))
  if(indicator_family=='binomial') hist(df[, get(indicator)]/df$N) else hist(df[, get(indicator)])
  dev.off()

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Define modeling space. In years only for now.
  if(yearload=='annual') period_map <-
                         make_period_map(modeling_periods = c(min(year_list):max(year_list)))
  if(yearload=='five-year') period_map <-
                            make_period_map(modeling_periods = seq(min(year_list),max(year_list),by=5))

  ## Make placeholders for covariates
  cov_layers <- gbd_cov_layers <- NULL

  ## Pull all covariate bricks/layers
  if (nchar(fixed_effects) > 0) {
    message('Grabbing raster covariate layers')

    effects <- trim(strsplit(fixed_effects, "\\+")[[1]])
    measures <- trim(strsplit(fixed_effects_measures, "\\+")[[1]])
    cov_layers <- load_and_crop_covariates_annual(covs            = effects,
                                                  measures        = measures,
                                                  simple_polygon  = simple_polygon,
                                                  start_year      = min(year_list),
                                                  end_year        = max(year_list),
                                                  interval_mo     = as.numeric(interval_mo))
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
    gaul_code <- rasterize(fe_subset_shape, fe_template, field = 'ADM0_CODE')
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
    plot(all_cov_layers[[covname]],main=covname,maxpixel=1e6)
  }
  dev.off()

  ## Check for data where covariate extraction failed
  # TODO CHECK THIS CHECK
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

  ## this will drop rows with NA covariate values
  the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))

  ## stop if this na omit demolished the whole dataset
  if(nrow(the_data) == 0) stop('You have an empty df, make sure one of your covariates was not NA everywhere.')

  ## seperated out into a different script
  if(as.logical(use_stacking_covs)){
    message('Fitting Stackers')
    # sourcing this script will run the child stackers:
    source(paste0(core_repo, '/mbg_central/share_scripts/run_child_stackers.R'))

    ## combine the children models
    the_data  <- cbind(the_data, do.call(cbind, lapply(lapply(child_model_names, 'get'), function(x) x[[1]])))
    child_model_objs <- setNames(lapply(lapply(child_model_names, 'get'), function(x) x[[2]]), child_model_names)

    ## fit GAM stacker -- SOON TO BE DEPRECIATED, right Daniel?
    stacked_results <- gam_stacker(the_data,
                                   model_names      = child_model_names,
                                   indicator        = indicator,
                                   indicator_family = indicator_family)

    ## return the stacked rasters
    stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                          period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                          child_models     = child_model_objs,
                                          stacker_model    = stacked_results[[2]],
                                          indicator_family = indicator_family,
                                          return_children  = TRUE,
                                          centre_scale_df  = covs_cs_df)

    ## plot stackers
    pdf(paste0(outputdir, 'stacker_rasters', pathaddin, '.pdf'))
    for(i in 1:length(stacked_rasters))
      plot(stacked_rasters[[i]],main=names(stacked_rasters[[i]]),maxpixel=ncell(stacked_rasters[[i]]))
    dev.off()

    message('Stacking is complete')
  } ## if(use_stacking_covs)


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

  ## make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d           = df,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset,
                             s2mesh = as.logical(use_s2_mesh),
                             s2params = s2_mesh_params)

  ## Build temporal mesh (standard for now)
  if (length(unique(year_list)) == 1)
    mesh_t <- NULL
  else
    mesh_t <- build_time_mesh(periods=eval(parse(text=mesh_t_knots)))


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
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = centre_scale_covs)

} else { ## END !SKIPTOINLA
  message(paste0('You have chosen to skip directly to INLA. Picking up objects from run_date ',skiptoinla_from_rundate))
  message('Now copying saved MBG inputs from that chosen run_date.')

  file.copy(from = '<<<< FILEPATH REDACTED >>>>',
            to = '<<<< FILEPATH REDACTED >>>>')
}

## reload data an prepare for MBG
load('<<<< FILEPATH REDACTED >>>>')

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
  for (ii in child_model_names) cov_list[[ii]] <- logit(cov_list[[ii]])

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
                                             add_ctry_res                = use_inla_country_res,
                                             ctry_re_prior               = ctry_re_prior,
                                             temporal_model_theta1_prior = rho_prior,
                                             no_gp                       = !as.logical(use_gp),
                                             coefs.sum1                  = coefs_sum1)

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

############ TESTING! SET UP SBH TO BE MODELLED AS GAUSSIAN
if(testgauss==TRUE){
  df[,sd_i:=0]
#  df[data_type == 'sbh', sd_i:= (qlogis(haz_upper)-qlogis(haz))/1.96] # jenky for now, later use draws for this MAKE SURE TO DO IN LOGIT SPACE BECAUSE THATS HOW ITLL GET MODELLED
  df[data_type == 'sbh', sd_i := sqrt(var_logit_lss_w)] # TODO add in the (small) model uncertainty
  # QUESTION: WHAT DO ABOUT SAMPLE SIZES HERE
  df[data_type == 'sbh', died := qlogis(haz)]
  df[,lik_fam_binom := ifelse(data_type=='sbh',0,1)]
  df[,lik_fam_gauss := ifelse(data_type=='sbh',1,0)]
scale_gaussian_variance_N<- FALSE
}

## Create SPDE INLA stack
input_data <- build_mbg_data_stack(df            = df, # note that merge (if using TMB) will return data in a different (but internally consistent) order, just different than df
                                   fixed_effects = all_fixed_effects,
                                   mesh_s        = mesh_s,
                                   mesh_t        = mesh_t, # not currently implemented with tmb
                                   use_ctry_res  = use_inla_country_res,
                                   use_nugget    = use_inla_nugget, # implemented with tmb
                                   exclude_cs    = child_model_names, # raw covs will get center scaled here though (see notes above)
                                   coefs.sum1    = coefs_sum1, # not currenlty implemented tmb
                                   tmb           = fit_with_tmb,
                                   scale_gaussian_variance_N = scale_gaussian_variance_N,
                                   shapefile_version = modeling_shapefile_version, 
                                   zl            = z_list, # if this is not zero and tmb==TRUE, it will trigger 3rd kronecker and fixed effects
                                   zcol          = zcol)   # must not be null if z_list is present

## combine all the inputs, other than cs_df these are not used if you are using TMB
stacked_input  <- input_data[[1]]
spde           <- input_data[[2]]
cs_df          <- input_data[[3]]

## Generate other inputs necessary
outcome <- df[[indicator]] # N+_i - event obs in cluster
N       <- df$N                  # N_i - total obs in cluster
weights <- df$weight

## catch in case there is no weight column
if(is.null(weights)){
  weights = rep(1,nrow(df))
}

tic("MBG - fit model") ## Start MBG - model fit timer

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
                         fe_sd_prior      = 1 / 9) ## this actually sets precision!. prec=1/9 -> sd=3
  }

  if(fit_with_tmb == TRUE) {
    message('Fitting model with TMB')
    message(sprintf('%i Data points and %i mesh nodes',nrow(df),length(input_data$Parameters$Epsilon_stz)))

    # save RDS file of input data for replication
    saveRDS(object = input_data, ## save this here in case predict dies
            file = '<<<< FILEPATH REDACTED >>>>')
    # run the model
    system.time(
    model_fit <- fit_mbg_tmb( lbdcorerepo     = core_repo,
                              cpp_template    = 'mbg_tmb_model', #_testing
                              tmb_input_stack = input_data,
                              control_list    = list(trace=1, eval.max=500, iter.max=300, abs.tol=1e-20),
                              optimizer       = 'nlminb', # TODO add optimx
                              ADmap_list      = NULL,
                              sparse_ordering = as.logical(sparse_ordering))
    )

    # clamping
    clamp_covs <- TRUE # TODO CONFIG THIS
    
    
  }

  saveRDS(object = model_fit, ## save this here in case predict dies
          file = '<<<< FILEPATH REDACTED >>>>')


}else{
  ## skipped fitting INLA so just load model and move to predict
  model_fit <- readRDS( file = '<<<< FILEPATH REDACTED >>>>')

}

toc(log = T) ## End MBG - model fit timer


## *** MODEL PREDICTION HAPPENS IN ANOTHER SCRIPT ***
## Run a Cholesky decomposition to make sure that the joint precision matrix
##  is positive definite (meaning that the model converged correctly)
dummy_data <- Cholesky(model_fit$sdrep$jointPrecision, super = TRUE)


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# save training data
write.csv(
  df,
  file = (paste0(outputdir, '/', indicator,'_trainingdata',pathaddin,'.csv')),
  row.names = FALSE
)

message('Done with model fitting.')

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
output_file <- paste0(outputdir, "run_summary_fitting_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}


# Write a an empty file to indicate done with this parallel FITTING script
write(NULL, file = paste0(model_fit_dir, "/fin_", pathaddin))

## Write CSV
write.csv(df_timer,file = output_file, row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
