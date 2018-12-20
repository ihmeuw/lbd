## Process arguments for this model job from submission.
  reg=as.character(commandArgs()[4])
  age=as.numeric(commandArgs()[5])
  run_date=as.character(commandArgs()[6])
  test=as.character(commandArgs()[7])
  holdout=as.character(commandArgs()[8])
  indicator=as.character(commandArgs()[9])
  indicator_group=as.character(commandArgs()[10])
  
  pathaddin = paste0('_bin',age,'_',reg,'_',holdout)
  message(pathaddin)
  message(run_date)
  message(test)
  
## Load an image of the main environment
  load(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))
  slots = as.numeric(slots)

## Reload functions
  setwd(repo)
  root <- ifelse(Sys.info()[1]=="Windows", "<<<<< FILEPATH REDACTED >>>>>", "<<<<< FILEPATH REDACTED >>>>>")
  package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                        <<<< FILEPATH REDACTED >>>>,
                        <<<< FILEPATH REDACTED >>>>)
  .libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library(). Necessary for seeg libraries.
  source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
  source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
  source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
  source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
  source('mbg_central/post_estimation_functions.R')
  source('mbg_central/gbd_functions.R')
  source('mbg_central/shiny_functions.R')
  source('mbg_central/stacking_functions.R')
  source('mbg_central/categorical_variable_functions.R')
  source('mbg_central/seegMBG_transform_functions.R')
  source('mbg_central/validation_functions.R') 

## Functions come from the saved imaged
  package_list <- c('rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr','foreign', 'magrittr', 'tictoc')
  for(package in package_list) {
    library(package, lib.loc = package_lib, character.only=TRUE)
  }

tic("Entire script") # Start master timer

## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cores_to_use <- round(slots*.5)

## Load simple polygon template to model over
  gaul_list <- get_gaul_codes(reg)
  suppressMessages(suppressWarnings(simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                                                               buffer = 0.4)))
  subset_shape     <- simple_polygon_list[[1]]
  simple_polygon   <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
  raster_list <- suppressMessages(suppressWarnings(build_simple_raster_pop(subset_shape)))
  simple_raster <- raster_list[['simple_raster']]
  pop_raster <- raster_list[['pop_raster']]

if(skip_stacking == FALSE) {

  ## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
    if(holdout!=0) {
      df <- as.data.table(stratum_qt[[paste0('region__',reg)]])
      oos_df <- df[fold == holdout, ]
      i <- length(unique(oos_df$year))
      periods <- data.frame(group = rep(1:i,5),years = rep(sort(unique(oos_df$year)),5))
      oos_df$period <- match(oos_df$year, periods$years) # add these to df
      df <- df[fold != holdout, ]
    }
    if(holdout==0) {
      df <- load_input_data(indicator = indicator,
                            simple = simple_polygon,
                            removeyemen = TRUE,
                            pathaddin = pathaddin,
                            years = 'annual',
                            use_share = TRUE)
    }
  
  ## Save distribution of data for this region
    png(paste0('<<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/', reg, '.png'))
    hist(df[, get(indicator)])
    dev.off()

  ## Define modeling space (right now, just in years)
    period_map <- make_period_map(modeling_periods = c(min(year_list):max(year_list))) 

  ## Make covariates conditional
    cov_layers = NULL
    gbd_cov_layers = NULL
    mbg_cov_layers = NULL

  ## Pull all covariate bricks/layers
    if(nchar(fixed_effects)> 0){
      message(fixed_effects)
      selected_fixed_effects <- strsplit(fixed_effects," ")
      selected_fixed_effects <- selected_fixed_effects[[1]][selected_fixed_effects[[1]] != "+"]
      selected_measures <- strsplit(fixed_effects_measures," ")
      selected_measures <- selected_measures[[1]][selected_measures[[1]] != "+"]
      cov_layers <- load_and_crop_covariates_annual(covs = selected_fixed_effects,                
                                                 measures = selected_measures,          
                                                 simple_polygon = simple_polygon,
                                                 start_year  = min(year_list),
                                                 end_year    = max(year_list),
                                                 interval_mo = 12,
                                                 agebin=1)
    }
    if(nchar(gbd_fixed_effects)>0){
      selected_gbd_fixed_effects <- strsplit(gbd_fixed_effects," ")
      selected_gbd_fixed_effects <- selected_gbd_fixed_effects[[1]][selected_gbd_fixed_effects[[1]] != "+"]
      ## Create list of gaul codes for region + any countries in the data that aren't in region (buffer zone)
      gbd_gaul_list <- unique(c(gaul_convert(unique(df[, country])), gaul_list))
      ## Use a layer from the geospatial cov_layers list as a template (raster of simple_polygon) to rasterize GBD cov spdfs.
      gbd_cov_layers <- suppressMessages(suppressWarnings(load_gbd_covariates(gbd_fixed_effects = selected_gbd_fixed_effects, 
                                                                              year_ids = year_list, 
                                                                              gaul_list = gbd_gaul_list,
                                                                              template = cov_layers[[1]][[1]])))
    }
    if(nchar(mbg_fixed_effects)>0){
      mbg_cov_layers <- suppressMessages(suppressWarnings(load_mbg_covariates(mbg_fixed_effects = mbg_fixed_effects, 
                                                                              simple_polygon = simple_polygon)))
    }

  ## Update all cov layers with an indicator variable on country
    all_cov_layers <- c(cov_layers, gbd_cov_layers)
    all_fixed_effects = paste(names(all_cov_layers), collapse = " + ")
  ## Make stacker-specific formulas where applicable 
    all_fixed_effects_brt <- all_fixed_effects
    source('mbg_central/stacking_functions.R')
  
  ## Process country-level effects, if included in model
  if(use_child_country_fes == TRUE | use_inla_country_fes == TRUE) {
    fe_gaul_list <- unique(c(gaul_convert(unique(df[, country])), gaul_list))
    fe_template = cov_layers[[1]][[1]]
    suppressMessages(suppressWarnings(simple_polygon_list <- load_simple_polygon(gaul_list = fe_gaul_list,
                                                                                 buffer = 0.4,
                                                                                 subset_only = TRUE)))
    fe_subset_shape     <- simple_polygon_list[[1]]
    gaul_code <- rasterize(fe_subset_shape, fe_template, field = 'GAUL_CODE')
    gaul_code = setNames(gaul_code,'gaul_code')
    gaul_code = create_categorical_raster(gaul_code)
    
    ## Update covlayers and add country fixed effects to the
    all_cov_layers = update_cov_layers(all_cov_layers, gaul_code)
    all_fixed_effects_cfes = paste(all_fixed_effects, paste(names(gaul_code)[1:length(names(gaul_code))], collapse = " + "), sep=" + ")
    
    ## Update specific stacker formulas (for now we just want country effects in BRT)
    all_fixed_effects_brt <- all_fixed_effects_cfes
  }

  ## Estimate stacking procedure (GAM, BRT, penalized regressions)
  tic("Stacking - all") # Start stacking master timer
  
  ## Figure out which models we're going to use
    child_model_names <- stacked_fixed_effects %>% 
                            gsub(" ", "", .) %>%
                            strsplit(., "+", fixed=T) %>%
                            unlist
    
    the_covs = format_covariates(all_fixed_effects)
  
  ## Copy the dataset to avoid unintended namespace conflicts
    the_data = copy(df)
  
  ## Shuffle the data into six folds
    n_stack_folds <- 5
    the_data = the_data[sample(nrow(the_data)),]
    the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=as.numeric(n_stack_folds),labels=FALSE)] #make folds
  
  ## Extract covariates to the points and subset data where its missing covariate values
    if(just_covs_model == FALSE) {
      save_mbg_input_cs <- TRUE
    }
  
  ## If modeling with just linear raw covariates, don't centerscale here. We don't need the CS'd values for stacking because we
  ##   aren't using the stackers, and if we CS here they get CS'd again in save_mbg_input().
    if(just_covs_model == TRUE) {
      save_mbg_input_cs <- FALSE
    }
    cs_covs = extract_covariates(the_data, all_cov_layers, return_only_results = T, centre_scale = T, period_var = 'year', period_map = period_map)
    the_data = cbind(the_data, cs_covs[[1]])
    covs_cs_df = cs_covs[[2]]
    the_data = na.omit(the_data, c(indicator, 'N', the_covs)) #this will drop rows with NA covariate values
  
  ## Also check to see if we have enough data in all countries to do country level fixed effects
    ct.ct  <- numeric()
    for(ii in 1:n_stack_folds){
      ct.ct[ii] <- length(unique(subset(the_data, fold_id == ii)$country))
    }
    
  ## If country count (ct.ct) has different numbers, remove gaul codes from formulas
    if(length(unique(ct.ct)) > 1){
      sep.fixed.effects <- strsplit(all_fixed_effects, " ")[[1]]
      sep.fixed.effects <- sep.fixed.effects[seq(1, length(sep.fixed.effects), by = 2)] ## remove +s
      sep.fixed.effects <- sep.fixed.effects[!apply(t(sep.fixed.effects), 1, FUN = function(ch){grepl(x = ch, pattern = "gaul_code_")})[, 1]]
      all_fixed_effects <- paste0(paste0(sep.fixed.effects[1:length(sep.fixed.effects) - 1], " + ", collapse = ""), sep.fixed.effects[length(sep.fixed.effects)])
    }
  
  ## Fit a GAM model
    if ('gam' %in% child_model_names) {
      tic("Stacking - GAM")     # Start stacking timer - GAM
      gam = fit_gam_child_model(df = the_data, #data frame
                                model_name = 'gam', #identifier for the child model-- needs to be consistent through the steps
                                fold_id_col = 'fold_id',
                                covariates = all_fixed_effects, #rhs formula
                                additional_terms = 'year', #column(s) in df that should be included in the fit. Ideally, there is a raster companion to the column. These columns are not splined
                                weight_column = 'weight', #column in the data frame specifying weights
                                bam =F, #should mgcv::bam be called rather than gam?
                                spline_args = list(bs = 'ts', k = 3), #spline arguments to be applied to the columns specified by the covariates argument
                                auto_model_select =T, #should the function override bam calls in low data situations (boosts probability of convergence)
                                indicator = indicator, #outcome variable column
                                indicator_family = indicator_family, #family of the outcome. Binomial and Gaussian are implemented.
                                cores = 10) #number of compute cores available
      toc(log = T)              # End stacking timer - GAM
    }
  
  ## Fit a GBM/BRT model
    if ('gbm' %in% child_model_names) {
      tic("Stacking - GBM")     # Start stacking timer - GBM
      gbm = fit_gbm_child_model(df = the_data,
                                model_name = 'gbm',
                                fold_id_col = 'fold_id',
                                covariates = all_fixed_effects_brt,
                                weight_column = 'weight',
                                tc = 3, #tree complexity, change back to 4 for real runs
                                lr = 0.005, #learning rate
                                bf = 0.75, #bag fraction
                                indicator = indicator,
                                indicator_family = indicator_family,
                                cores = cores_to_use)
      toc(log = T)              # End stacking timer - GBM
    }
  
  ## Fit some nets
  ## Lasso
    if ('lasso' %in% child_model_names) {
      tic("Stacking - lasso")        # Start stacking timer - lasso
      lasso = fit_glmnet_child_model(df = the_data,
                                     model_name = 'lasso',
                                     covariates =all_fixed_effects,
                                     fold_id_col = 'fold_id',
                                     additional_terms = NULL,
                                     indicator_family = indicator_family,
                                     indicator = indicator,
                                     cores = cores_to_use,
                                     alpha = 0,
                                     weight_column = 'weight')
      toc(log = T)                   # End stacking timer - lasso
    }
  
  ## Ridge
    if ('ridge' %in% child_model_names) {
      tic("Stacking - ridge")        # Start stacking timer - ridge
      ridge = fit_glmnet_child_model(df = the_data,
                                     model_name = 'ridge',
                                     covariates = all_fixed_effects,
                                     fold_id_col = 'fold_id',
                                     additional_terms = NULL,
                                     indicator_family = indicator_family,
                                     indicator = indicator,
                                     cores = cores_to_use,
                                     alpha = 1,
                                     weight_column = 'weight')
      toc(log = T)                   # End stacking timer - ridge
    }
  
  ## Enet
    if ('enet' %in% child_model_names) {
      tic("Stacking - enet")         # Start stacking timer - enet
      enet = fit_glmnet_child_model(df = the_data,
                                    model_name = 'enet',
                                    covariates = all_fixed_effects,
                                    fold_id_col = 'fold_id',
                                    additional_terms = NULL,
                                    indicator_family = indicator_family,
                                    indicator = indicator,
                                    cores = cores_to_use,
                                    alpha = .5,
                                    weight_column = 'weight')
      toc(log = T)                   # End stacking timer - enet
    }
  
  ## Combine the children models
    the_data = cbind(the_data, do.call(cbind, lapply(lapply(child_model_names, 'get'), function(x) x[[1]])))
    child_model_objs = setNames(lapply(lapply(child_model_names, 'get'), function(x) x[[2]]), child_model_names)
  
  ## Fit stacker
    stacked_results = gam_stacker(the_data, #the dataset in data table format
                                  model_names= child_model_names, #prefixes of the models to be stacked
                                  indicator = indicator, #the indicator of analysis
                                  indicator_family = indicator_family) #indicator family (e.g. binomial)
  
  ## Return the stacked rasters
    stacked_rasters = make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                         period = min(period_map[, period_id]):max(period_map[, period_id]), #period of analysis, NULL returns all periods
                                         child_models = child_model_objs, #model objects fitted to full data
                                         stacker_model = stacked_results[[2]],
                                         indicator_family = indicator_family,
                                         return_children = T,
                                         centre_scale_df = covs_cs_df)
  
    all_fixed_effects <- stacked_fixed_effects
    if(use_inla_country_fes) all_fixed_effects = paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep=" + ")
    if(use_year_in_inla) all_fixed_effects = paste(all_fixed_effects, 'year_cov', sep=" + ")
  
  ## Copy things back over to df
    df = copy(the_data)
  
  ## Remove the covariate columns so that there are no name conflicts when they get added back in
    df = df[,paste0(the_covs) := rep(NULL, length(the_covs))]
  
  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
    df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]
  
  ## Build spatial mesh over modeling area
    if(grepl("geos", Sys.info()[4])){
      INLA:::inla.dynload.workaround()
    }
    mesh_s <- build_space_mesh(d = df,
                               simple = simple_polygon,
                               max_edge = mesh_s_max_edge,
                               mesh_offset = mesh_s_offset)
  
  ## Build temporal mesh using knots from config
    mesh_t <- build_time_mesh(periods=eval(parse(text=mesh_t_knots)))
  
  ## Create a full raster list to carry though to the shiny/next steps
    full_raster_list = c(unlist(stacked_rasters),unlist(all_cov_layers))
    child_mod_ras = full_raster_list[child_model_names]
  
  toc(log = T) # End stacking master timer
  
  ## Save all inputs for MBG model into correct location on 
    save_mbg_input(indicator = indicator, 
                   indicator_group = indicator_group,
                   df = df,
                   simple_raster = simple_raster,
                   mesh_s = mesh_s,
                   mesh_t = mesh_t, 
                   cov_list = full_raster_list,
                   pathaddin = pathaddin,
                   run_date = run_date,
                   child_model_names = child_model_names,
                   all_fixed_effects = all_fixed_effects,
                   period_map = period_map,
                   centre_scale = save_mbg_input_cs) #specify by region

} # Finish loop for stacking 

## Reload data an prepare for MBG
  load(paste0('<<< FILEPATH REDACTED >>>', indicator_group, '/', indicator, '/model_image_history/', run_date, pathaddin, '.RData'))

## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tic("MBG - all") # Start MBG master timer

## For stacking, overwrite the columns matching the model_names 
  df = df[,paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]
  message(paste0("Dropping ", length(df[is.na(get(grep(names(df), pattern = "gaul_code_", value = T)[1])), period]), " rows missing country random effects indicator."))
  df <- df[!is.na(get(grep(names(df), pattern = "gaul_code_", value = T)[1])), ]            

## Generate MBG formula for INLA call
  mbg_formula <- build_mbg_formula_with_priors(fixed_effects = all_fixed_effects,
                                               add_nugget = use_inla_nugget,
                                               add_ctry_res = use_inla_country_res,
                                               no_gp = drop_gp_from_formula,
                                               coefs.sum1 = use_coefs.sum1)
  if(grepl("geos", Sys.info()[4])){
    INLA:::inla.dynload.workaround()
  }
  
## Create SPDE INLA stack
  input_data <- build_mbg_data_stack(df = df,
                                     fixed_effects = all_fixed_effects,
                                     mesh_s = mesh_s,
                                     mesh_t = mesh_t,
                                     use_ctry_res = use_inla_country_res,
                                     use_nugget = use_inla_nugget,
                                     coefs.sum1 = use_coefs.sum1,
                                     exclude_cs = exclude_cs_option)

## Combine all the inputs
  stacked_input <- input_data[[1]]
  spde <- input_data[[2]]
  cs_df <- input_data[[3]]

## Generate other inputs necessary
  outcome=df[[indicator]] # N+_i - event obs in cluster
  N=df$N                  # N_i - total obs in cluster
  weights=df$weight

## Catch in case there is no weight column
  if(is.null(weights)){
    weights = rep(1,nrow(df))
  }

## Fit MBG model in INLA
  if(skip_fit == FALSE) {
  tic("MBG - fit model") # Start MBG - model fit timer
  
  ## Fit MBG model
    source('mbg_central/mbg_functions.R')
    if(grepl("geos", Sys.info()[4])){
      INLA:::inla.dynload.workaround()
    }
    model_fit <- fit_mbg(indicator_family = indicator_family,
                         stack.obs = stacked_input,
                         spde = spde,
                         cov = outcome,
                         N = N,
                         int_prior_mn = intercept_prior,
                         f_mbg = mbg_formula,
                         run_date = run_date,
                         keep_inla_files = keep_inla_files,
                         cores = cores_to_use,
                         wgts = weights)
    toc(log = T) # End MBG - model fit timer
    
  }
  
  ## If loading a previously fit model instead...
  if(skip_fit == TRUE) {
    load(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_model_eb_bin0_', reg, '_0.RData'))
    model_fit <- res_fit
  }

tic("MBG - predict model") # Start MBG - model predict timer

## Predict cell-wise draws object using INLA fit to sample from full posterior.
  source('mbg_central/mbg_functions.R')

## Create vector of chunk sizes
  max_chunk <- 100
  samples <- as.numeric(samples)
  chunks <- rep(max_chunk, samples %/% max_chunk)
  if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
  pm <- lapply(chunks, function(samp) {
    predict_mbg(res_fit       = model_fit,
                cs_df         = cs_df,
                mesh_s        = mesh_s,
                mesh_t        = mesh_t,
                cov_list      = cov_list,
                samples       = samp,
                simple_raster = simple_raster,
                simple_polygon = simple_polygon,
                transform     = transform,
                pred_gp       = use_gp_in_predict,
                max_period    = 17,
                full_df       = df,
                ho_id         = holdout,
                coefs.sum1    = use_coefs.sum1)[[3]]
  })

## Make cell preds and a mean raster
  cell_pred <- do.call(cbind, pm)

toc(log = T) # Stop MBG - model predict timer


##################################################################################################
############################### SAVE ALL OUTPUT ##################################################
##################################################################################################
## Save MBG outputs in standard outputs folder structure
  save_mbg_preds(config     = config,
                 time_stamp = time_stamp,
                 run_date   = run_date,
                 mean_ras   = NULL,
                 sd_ras     = NULL,
                 res_fit    = model_fit,
                 cell_pred  = cell_pred,
                 df         = df,
                 pathaddin  = pathaddin)

toc(log = T) # End master timer

## Format timer
  ticlog <- tic.log(format = F)
  df_timer <- generate_time_log(ticlog)
  df_timer[, region := reg]
  df_timer[, holdout := holdout]
  setcolorder(df_timer, c("region", "holdout", "step", "time"))

## Pull run time for this run
  run_time_all <- df_timer[step == "Entire script", time]

## Write to a run summary csv file in the output directory
  output_dir <- paste0('<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, "/")
  output_file <- paste0(output_dir, "run_summary_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
  if (file.exists(output_file)) {
    file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
    df_timer <- rbind(file_contents, df_timer)
  }

## Write CSV
  write.csv(df_timer, 
            file      = output_file, 
            row.names = F)
