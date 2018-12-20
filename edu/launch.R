###########################################################################################
###########################################################################################
## Run MBG model educational attainment.
## Nick Graetz
###########################################################################################
###########################################################################################

## Set repo location and indicator group
  repo <- <<<< FILEPATH REDACTED >>>>
  indicator_group <- 'education'
  indicator <- as.character(commandArgs()[5])
  
## Load libraries and miscellaneous MBG project functions.
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
  source('mbg_central/holdout_functions.R')
  source('mbg_central/categorical_variable_functions.R')
  source('mbg_central/validation_functions.R') 
  source('mbg_central/validation_report_functions.R') 
  source('mbg_central/seegMBG_transform_functions.R')     
  package_list <- c('foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr')
  for(package in package_list) {
    library(package, lib.loc = package_lib, character.only=TRUE)
  }
  
  crossval <- TRUE
  rerun_run_date <- NA
  
## Read config file and save all parameters in memory
  config <- load_config(repo = repo,
                        indicator_group = indicator_group,
                        indicator = indicator,
                        config_name = as.character(commandArgs()[4]))
  Regions <- strsplit(region_list," ")
  Regions <- Regions[[1]][Regions[[1]] != "+"]
  message(fixed_effects)
  year_list <- eval(parse(text=year_list))
  
## Create run date in correct format
  run_date <- make_time_stamp(time_stamp)
  
## Create directory structure for this model run
  create_dirs(indicator_group = indicator_group,
              indicator = indicator)
  
## Load simple polygon template to model over
  gaul_list <- get_gaul_codes('africa')
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                             buffer = 0.4)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
## Load input data and make sure it is cropped to modeling area
  df_list <- load_input_data(indicator = indicator,
                             simple = simple_polygon,
                             removeyemen = TRUE,
                             update_run_date = TRUE,
                             use_share = TRUE)
  df <- df_list[[1]]
  if(is.na(rerun_run_date)) run_date <- df_list[[2]]

## Add GAUL_CODE and region to df given the lat/longs. We need this to stratify in holdout function.
  df <- add_gauls_regions(df = df,
                          simple_raster = simple_raster)

## Run function to create holdouts (returns list of data.frames with an additional "folds" column)
  if(crossval==TRUE) {
    table(df$region, df$year)
    long_col = 'longitude'
    lat_col = 'latitude'
    n_folds = as.numeric(n_folds)
    stratum_qt <- make_folds(data = df, n_folds = n_folds, spat_strat = spat_strat,
                             temp_strat = temp_strat, strat_cols = "region",
                             ts = 20, mb = 10)
  }
  
## Set strata as character vector of each strata (in my case, just stratifying by region)
  strata <- Regions

## ~~~~~~~~~~~~~~~~~~~~~~~~  Parallel MBG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~ Submit job by strata/holdout  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Get GBD estimates for calibrating after predicting.
  if(indicator=='edu_mean_male' | indicator == 'edu_mean_20_24_male') {
    raking_sex <- 1
    raking_pop <- 'total'
  }
  if(indicator=='edu_mean' | indicator=='edu_mean_logit18' | indicator=='edu_mean_log' | indicator=='edu_mean_20_24') {
    raking_sex <- 2
    raking_pop <- 'wocba'
  }
  
  if(indicator=='edu_mean' | indicator=='edu_mean_logit18' | indicator=='edu_mean_log' | indicator == 'edu_mean_male') {
    year_ids  <- c(min(year_list):max(year_list))
    gaul_list <- get_gaul_codes('africa')
    metadata <- get_covariate_metadata()
    source(<<<< FILEPATH REDACTED >>>>)
    gbd_estimates <- get_covariate_estimates(covariate_id = metadata[covariate_name_short == 'education_yrs_pc', covariate_id], gbd_round_id = 4)
    
    # Merge pops to collapse age-specific to 10-54 women
    source(<<<< FILEPATH REDACTED >>>>)
    gbd_pops <- get_population(age_group_id = paste(as.character(unique(gbd_estimates[, age_group_id])), collapse = " "),
                               sex_id = "1 2",
                               location_id = paste(as.character(unique(gbd_estimates[, location_id])), collapse = " "),
                               year_id = paste(as.character(unique(gbd_estimates[, year_id])), collapse = " "))
    gbd_pop_estimates <- merge(gbd_estimates, gbd_pops, by=c('age_group_id','sex_id','location_id','year_id'))
    # Subset 
    gbd_pop_estimates <- gbd_pop_estimates[age_group_id %in% c(7:15) & sex_id == raking_sex, ]
    # Pop-weighted collapse to mean 
    gbd_means <- gbd_pop_estimates[,.(gbd_mean=weighted.mean(x=mean_value,w=population)), by=.(location_id,year_id)]
    # Convert to GAUL_CODE
    gbd_estimates <- gbd_means
    gaul_to_loc_id <- fread(<<<< FILEPATH REDACTED >>>>)
    names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"
    gbd_estimates <- gbd_estimates[year_id %in% year_ids,]
    gbd_estimates <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
    gbd_estimates <- gbd_estimates[GAUL_CODE %in% gaul_list,]
    names(gbd_estimates)[names(gbd_estimates)=="GAUL_CODE"] <- "name"
    names(gbd_estimates)[names(gbd_estimates)=="year_id"] <- "year"
    names(gbd_estimates)[names(gbd_estimates)=="gbd_mean"] <- "mean"
    gbd_estimates <- gbd_estimates[, c('name', 'year', 'mean'), with = FALSE]
    gbd <- gbd_estimates
  }
  if(indicator=='edu_mean_20_24' | indicator=='edu_mean_20_24_male') {
    year_ids  <- c(min(year_list):max(year_list))
    gaul_list <- get_gaul_codes('africa')
        metadata <- get_covariate_metadata()
    source(<<<< FILEPATH REDACTED >>>>)
    gbd_estimates <- get_covariate_estimates(covariate_id = metadata[covariate_name_short == 'education_yrs_pc', covariate_id], gbd_round_id = 4)
    # Merge pops to collapse age-specific to this subset
    source(<<<< FILEPATH REDACTED >>>>)
    gbd_pops <- get_population(age_group_id = paste(as.character(unique(gbd_estimates[, age_group_id])), collapse = " "),
                               sex_id = "1 2",
                               location_id = paste(as.character(unique(gbd_estimates[, location_id])), collapse = " "),
                               year_id = paste(as.character(unique(gbd_estimates[, year_id])), collapse = " "))
    gbd_pop_estimates <- merge(gbd_estimates, gbd_pops, by=c('age_group_id','sex_id','location_id','year_id'))
    # Subset 
    gbd_pop_estimates <- gbd_pop_estimates[age_group_id == 9 & sex_id == raking_sex, ]
    # Pop-weighted collapse to mean for 10-54
    gbd_means <- gbd_pop_estimates[,.(gbd_mean=weighted.mean(x=mean_value,w=population)), by=.(location_id,year_id)]
    # Convert to GAUL_CODE
    gbd_estimates <- gbd_means
    gaul_to_loc_id <- fread(<<<< FILEPATH REDACTED >>>>)
    names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"
    gbd_estimates <- gbd_estimates[year_id %in% year_ids,]
    gbd_estimates <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
    gbd_estimates <- gbd_estimates[GAUL_CODE %in% gaul_list,]
    names(gbd_estimates)[names(gbd_estimates)=="GAUL_CODE"] <- "name"
    names(gbd_estimates)[names(gbd_estimates)=="year_id"] <- "year"
    names(gbd_estimates)[names(gbd_estimates)=="gbd_mean"] <- "mean"
    gbd_estimates <- gbd_estimates[, c('name', 'year', 'mean'), with = FALSE]
    gbd <- gbd_estimates
  }

## Submit all parallel MBG models.
  parallel_script <- 'parallel_model_full'
  slots <- as.numeric(slots)
  sharedir       <- sprintf('<<<< FILEPATH REDACTED >>>>/%s/%s',indicator_group,indicator)
  loopvars <- NULL
  for(r in strata){
    
    fin <- 'FALSE'
    if(!is.na(rerun_run_date)) {
      # Only qsub if that strata has any failed holdout jobs ("fin" file doesn't exist)
      fins <- 1
      for(holdout in c(0:n_folds)) {
        holdout_fin <- file.exists(paste0('<<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/fin__bin0_', r, '_', holdout))
        fins <- c(fins, holdout_fin)
      }
      if(sum(fins)==length(fins)) {
        # If this strata was completed finished, mark to not rerun
        fin <- 'TRUE'
      }
      # If we are rerunning this strata because one or more holdouts broke, delete fin files for all holdouts.
      # This is because if one broke, we have to relaunch all because we have to reorganize holdout ids.
      if(fin==FALSE) {
        delete_fin_files <- list.files(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date), 
                                       pattern = paste0('fin__bin0_', r, '_'), 
                                       full.names = TRUE)
        lapply(delete_fin_files, unlink)
      }
    }
    
    ## Submit models
    for(holdout in c(0:n_folds)) {

      if(fin==FALSE) {
        
        qsub <- make_qsub(code = parallel_script,
                          reg = r,
                          saveimage = TRUE,
                          test = TRUE,
                          holdout = holdout,
                          log_location = 'sharedir',
                          proj = 'proj_geo_nodes',
                          geo_nodes = TRUE)
        system(qsub)
        
      }
      
      loopvars <- rbind(loopvars, c(r,holdout))
      
    }
  }

## ~~~~~~~~~~~~~~~~~~~~~~~~ Post-Estimation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Check to make sure models are done before continuing
  waitformodelstofinish()

## Save strata for Shiny to use in producing aggregated fit statistics
  dir.create(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/fit_stats'))
  save(strata, file = paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/fit_stats/strata.RData'))
  
## Loop over regions, calibrate to GBD national estimates, and save relevant objects.
for(reg in strata){
  
  message(reg)

  load(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))

  ## Get aggregated estimates for all admin0. Aggregate to level you rake to.
    simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = FALSE)
    subset_shape   <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
  
  ## Pull 2000-2015 annual population brick 
    pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',                
                                                         measures = raking_pop,          
                                                         simple_polygon = simple_polygon,
                                                         start_year  = min(year_list),
                                                         end_year    = max(year_list),
                                                         interval_mo = 12,
                                                         agebin=1)
    pop_raster_annual  <- pop_raster_annual[[1]]
    pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
    pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
    pop_raster_annual  <- mask(pop_raster_annual, simple_raster)
  
  ## Create population weights using the annual brick and feed custom year argument to aggregation function
  pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                          simple_raster = simple_raster,
                                          pop_raster    = pop_raster_annual,
                                          gaul_list     = get_gaul_codes(reg))
  
  cond_sim_raw_adm0 <- make_condSim(pop_wts_object = pop_wts_adm0,
                                    gaul_list      = get_gaul_codes(reg),
                                    admin_level    = 0,
                                    cell_pred      = cell_pred,
                                    summarize      = TRUE,
                                    years          = c(min(year_list):max(year_list)))

  rf <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
                              gaul_list   = gaul_list,
                              rake_to     = gbd)
  
  raked_cell_pred <- rake_predictions(raking_factors = rf,
                                      pop_wts_object = pop_wts_adm0,
                                      cell_pred      = cell_pred,
                                      logit_rake     = FALSE)

  ## Summarize raked predictions for each cell
  mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
                                         mask                 = simple_raster,
                                         return_as_raster     = TRUE,
                                         summary_stat         = 'mean')

  raked_mean_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                                                     mask                 = simple_raster,
                                                     return_as_raster     = TRUE,
                                                     summary_stat         = 'mean')

  assign(sprintf('%s_rf',reg),rf)
  assign(sprintf('%s_mean_raster',reg),mean_raster)
  assign(sprintf('%s_raked_mean_raster',reg),raked_mean_raster)

  save(
    raked_cell_pred,
    file = (paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/', indicator,'_raked_cell_draws_eb_bin0_',reg,'_0.RData')),
    compress = TRUE
  )
  
  rm(mean_raster); rm(raked_mean_raster); rm(raked_cell_pred); rm(cell_pred); rm(pop_wts)
  
}

## Combine regions raster to all Africa and save objects.
  rf <- do.call(rbind.fill, lapply(paste0(Regions, '_rf'), get))
  if(length(paste0(Regions, '_mean_raster'))!=1) m = do.call(raster::merge, lapply(paste0(Regions, '_mean_raster'), get))
  if(length(paste0(Regions, '_mean_raster'))==1) m = get(paste0(Regions, '_mean_raster'))
  if(length(paste0(Regions, '_raked_mean_raster'))!=1) m_raked = do.call(raster::merge, lapply(paste0(Regions, '_raked_mean_raster'), get))
  if(length(paste0(Regions, '_raked_mean_raster'))==1) m_raked = get(paste0(Regions, '_raked_mean_raster'))
  save_post_est(rf,'csv','rf')
  save_post_est(m,'raster','mean_raster')
  save_post_est(m_raked,'raster','mean_raked_raster')
  
## Submit in- and out-of-sample compile for cross-valiation.
  config <- load_config(repo = repo,
                        indicator_group = indicator_group,
                        indicator = indicator,
                        post_est_only = TRUE,
                        run_date = run_date)
  year_list <- eval(parse(text=year_list))
  results_dir <- paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date)
  pull_region_data <- function(reg) {
    d <- fread(paste0(results_dir, '/input_data_bin0_', reg, '_0.csv'))
    return(d)
  }
  df <- rbindlist(lapply(get_output_regions(in_dir = results_dir), pull_region_data))
  
## Submit in- and out-of-sample compile jobs.
  if(crossval==FALSE) {
    run_in_oos <- get_is_oos_draws(ind_gp = indicator_group,
                                   ind = indicator,
                                   rd = run_date,
                                   ind_fm = indicator_family,
                                   model_domain = 'africa',
                                   age = 0,
                                   nperiod = length(year_list),
                                   yrs = min(year_list):max(year_list),
                                   get.oos = FALSE, 
                                   write.to.file = TRUE,
                                   input_data = input_df)
  }
  if(crossval==TRUE) {
    run_in_oos <- get_is_oos_draws(ind_gp = indicator_group,
                                   ind = indicator,
                                   rd = run_date,
                                   ind_fm = indicator_family,
                                   model_domain = 'africa',
                                   age = 0,
                                   nperiod = length(year_list),
                                   yrs = min(year_list):max(year_list),
                                   get.oos = TRUE, 
                                   write.to.file = TRUE)
  }
  d <- fread(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/output_draws_data.csv'))
  d <- d[year >= 2000 & year <= 2002, year := 2000]
  d <- d[year >= 2003 & year <= 2007, year := 2005]
  d <- d[year >= 2008 & year <= 2012, year := 2010]
  d <- d[year >= 2013 & year <= 2017, year := 2015]
  df <- d
  
  for(admin in c('country','ad1','ad2','ho_id')) {
    pv_table <- get_pv_table(d = df,
                             indicator = indicator,
                             indicator_group = indicator_group,
                             rd = run_date,
                             aggregate_on = admin, # column of spatial aggregation (admin1 or admin2?)
                             draws = as.numeric(samples), # number of draws
                             coverage_probs = c(95), # coverage percentage
                             result_agg_over =  c('year','oos'), # final table aggregates over
                             weighted = TRUE, # weight PV metrics on SS?
                             family = ind_family, # distribution
                             plot = TRUE,
                             plot_ci = TRUE,
                             out.dir = paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/'))
    write.csv(pv_table, paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/pv_table_', admin, '.csv'))
  }
  
################ Calculate all admin-level aggregates, for both attainment and probabilities regarding targets ##########################
  submit_aggregation_script(indicator       = indicator,
                            indicator_group = indicator_group, 
                            run_date        = run_date, 
                            raked           = c(TRUE), 
                            pop_measure     = pop_measure, 
                            overwrite       = T, 
                            ages            = 0, 
                            holdouts        = 0,
                            regions         = strata, 
                            corerepo        = repo, 
                            log_dir         = log_dir, 
                            geo_node        = TRUE, 
                            slots           = 8)
  
  waitforaggregation(rd = run_date, indic = indicator, ig = indicator_group,
                     ages = 0,
                     regions = strata,
                     holdouts = 0,
                     raked = c(T))
  
  combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                      ages = 0, 
                      regions = strata,
                      holdouts = 0,
                      raked = c(T))
  
  load(paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_raked_admin_draws_eb_bin0_0.RData'))
  
  sum_admin_probs <- function(admin_level, target_type, goal_threshold, indicator, indicator_group, run_date) {
    
    admin <- get(paste0('admin_', admin_level))
    
    classify_pixel <- function(x,...) {
      if(target_type == 'less') {
        x[!is.na(x) & x > goal_threshold] <- 0
        x[!is.na(x) & x != 0 & x < goal_threshold] <- 1
      }
      if(target_type == 'greater') {
        x[!is.na(x) & x < goal_threshold] <- 0
        x[!is.na(x) & x != 0 & x > goal_threshold] <- 1
      }
      return(x)
    }
    
    adm_code <- paste0('ADM', admin_level, '_CODE')
    
    ## Admin probabilities
    message(paste0('Calculating ', adm_code, ' probability that ', indicator, ' is ', target_type, ' than ', goal_threshold, '...'))
    admin_probs <- admin[, lapply(.SD, classify_pixel), by=c(adm_code,'year'), .SDcols=grep("^V", names(admin))]
    admin_probs <- admin_probs[, p_goal := apply(.SD, 1, mean), by=c(adm_code,'year'), .SDcols=grep("^V", names(admin_probs))]
    admin_probs <- admin_probs[, c(grep("^ADM", names(admin_probs), value = TRUE),'year','p_goal'), with=FALSE]
    write.csv(admin_probs, paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/admin', admin_level, '_probs.csv'))
    
    ## Admin summaries
    message(paste0('Calculating ', adm_code, ' mean/upper/lower...'))
    admin_summary <- admin[, lower := apply(.SD, 1, quantile, c(.025), na.rm=T), .SDcols=grep("^V", names(admin))]
    admin_summary <- admin_summary[, mean := apply(.SD, 1, mean), .SDcols=grep("^V", names(admin_summary))]
    admin_summary <- admin_summary[, upper := apply(.SD, 1, quantile, c(.975), na.rm=T), .SDcols=grep("^V", names(admin_summary))]
    admin_summary <- admin_summary[, c(grep("^ADM", names(admin_summary), value = TRUE),'year','mean','upper','lower'), with=FALSE]
    write.csv(admin_summary, paste0('<<<< FILEPATH REDACTED >>>>', indicator_group, '/', indicator, '/output/', run_date, '/admin', admin_level, '_summary.csv'))
    
  }
  
  lapply(c(0,1,2), sum_admin_probs,
         goal_threshold = 6,
         target_type = 'greater',
         indicator_group = indicator_group,
         indicator = indicator,
         run_date = run_date)
  
