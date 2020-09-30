###########################################################################################
###########################################################################################
## Run MBG model
##    1. Set up and create holdouts
##    2. Submit script to run strata stacking, fitting, and predicting in parallel
##          - For this example, I'm only stratifying by region
##    3. Combine and rake
##    4. Save results to Shiny
###########################################################################################
###########################################################################################

## Set repo location and indicator group
  repo <- [REPO_PATH]
  indicator_group <- [INDCATOR_GROUP]
  indicator <- [INDICATOR]
  Regions=c('cssa','essa','name','sssa','wssa') #[STRATA]

## Load libraries and miscellaneous MBG project functions.
  setwd(repo)
  package_lib <- '<<< FILEPATH REDACTED >>>' # Library for all MBG versioned packages. Ensures that none of this code is
  #    dependent on the machine where the user runs the code.
  .libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library(). 
  #    Necessary for seeg libraries.
  source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
  source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
  source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
  source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
  source('mbg_central/post_estimation_functions.R')
  source('mbg_central/gbd_functions.R')
  source('mbg_central/shiny_functions.R')
  source('mbg_central/holdout_functions.R')
  source('mbg_central/seegMBG_transform_functions.R')     # Using Roy's edit for now that can take temporally varying covariates,
  #   TODO: will need to send pull request to seegMBG of github
  package_list <- c('foreign', 'rgeos', 'data.table','raster','rgdal','INLA','seegSDM','seegMBG','plyr','dplyr')
  for(package in package_list) {
    library(package, lib.loc = package_lib, character.only=TRUE)
  }

## Read config file and save all parameters in memory
  config <- load_config(repo = repo,
                        indicator_group = indicator_group,
                        indicator = indicator)

## Create run date in correct format
  run_date <- make_time_stamp(time_stamp)

## Create directory structure for this model run
  create_dirs(indicator_group = indicator_group,
              indicator = indicator)

## Load simple polygon template to model over
  gaul_list <- get_gaul_codes('africa')
  simple_polygon <- load_simple_polygon(gaul_list = gaul_list,
                                        buffer = 0.4)
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

## Load input data and make sure it is cropped to modeling area
  df <- load_input_data(indicator = indicator,
                        simple = simple_polygon,
                        removeyemen = TRUE)

## Add GAUL_CODE and region to df given the lat/longs. We need this to stratify in holdout function.
  df <- add_gauls_regions(df = df,
                          simple_raster = simple_raster)

## Run function to create holdouts (returns list of data.frames with an additional "folds" column)
table(df$region, df$year)
long_col = 'longitude'
lat_col = 'latitude'
n_folds = as.numeric(n_folds)
stratum_qt <- make_folds(data = df, n_folds = n_folds, spat_strat = spat_strat,
                         temp_strat = temp_strat, strat_cols = "region",
                         ts = 20, mb = 10)

## Set strata as character vector of each strata (in my case, just stratifying by region whereas U5M stratifies by region/age)
strata <- Regions

## ~~~~~~~~~~~~~~~~~~~~~~~~  Parallel MBG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~ Submit job by strata/holdout  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# get gbd estimates for raking
gbd <- load_gbd_data(gbd_type = "output",
                     gbd_name = 322,
                     gaul_list = get_gaul_codes('africa'),
                     measure_id = 6,
                     age_group_id = 1,
                     metric_id = 3)

## This is where we begin to do things by region
# if(stacking_method=='brt') parallel_script <- 'mbg_byregion_brt'
# if(stacking_method=='gam') parallel_script <- 'mbg_byregion'
parallel_script <- 'stacking_mbg'
slots=45

for(r in strata){
  for(holdout in 0:n_folds) {
    
    qsub <- make_qsub(code = parallel_script,
                      reg = r,
                      saveimage = TRUE,
                      test = TRUE,
                      holdout = holdout)
    
    system(qsub)
    
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~ Post-Estimation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## check to make sure models are done before continuing
waitformodelstofinish()

## Save strata for Shiny to use in producing aggregated fit statistics
save(strata, file = '<<< FILEPATH REDACTED >>>')

# Post estimation to be done by strata (in my case, this is just region)
# To-Do: parallelize this loop
for(reg in strata){
  message(reg)
  
  load('<<< FILEPATH REDACTED >>>')
  
  ## get aggregated estimates for all admin0. Aggregate to level you rake to
  load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = TRUE)
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                          simple_raster = simple_raster,
                                          pop_raster    = pop_raster,
                                          gaul_list     = get_gaul_codes(reg))
  
  cond_sim_raw_adm0 <- make_condSim(pop_wts_object = pop_wts_adm0,
                                    gaul_list      = get_gaul_codes(reg),
                                    admin_level    = 0,
                                    cell_pred      = cell_pred,
                                    summarize      = TRUE)
  
  
  ## Get raking factors
  rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
                              gaul_list   = gaul_list,
                              rake_to     = gbd)
  
  raked_cell_pred <- rake_predictions(raking_factors = rf,
                                      pop_wts_object = pop_wts_adm0,
                                      cell_pred      = cell_pred)
  
  ## summarize raked predictions for each cell
  mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
                                         mask                 = simple_raster,
                                         return_as_raster     = TRUE,
                                         summary_stat         = 'mean')
  
  raked_mean_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                                               mask                 = simple_raster,
                                               return_as_raster     = TRUE,
                                               summary_stat         = 'mean')
  
  for(ad in 1:2){
    pop_wts <- make_population_weights(admin_level   = ad,
                                       simple_raster = simple_raster,
                                       pop_raster    = pop_raster,
                                       gaul_list     = get_gaul_codes(reg))
    
    condsim <-    make_condSim(pop_wts_object = pop_wts,
                               cell_pred      = raked_cell_pred,
                               gaul_list      = get_gaul_codes(reg),
                               admin_level    = ad,
                               summarize      = TRUE)
    
    condsim <- data.table(cbind(split_geo_names(as.matrix(condsim)),mean=unname(condsim) ))
    assign(paste0('cond_sim_raked_adm',ad,'_',reg),condsim)
    rm(condsim)
    
  }
  
  assign(sprintf('%s_rf',reg),rf)
  assign(sprintf('%s_mean_raster',reg),mean_raster)
  assign(sprintf('%s_raked_mean_raster',reg),raked_mean_raster)
  rm(mean_raster); rm(raked_mean_raster); rm(raked_cell_pred); rm(cell_pred); rm(pop_wts)
  
}


# combine regions raster (child only for now)
rf <- do.call(rbind.fill, list(name_rf,
                               essa_rf,
                               sssa_rf,
                               wssa_rf,
                               cssa_rf))
m = do.call(raster::merge,list(name_mean_raster,
                               essa_mean_raster,
                               sssa_mean_raster,
                               wssa_mean_raster,
                               cssa_mean_raster))
m_raked = do.call(raster::merge,list(name_raked_mean_raster,
                                     essa_raked_mean_raster,
                                     sssa_raked_mean_raster,
                                     wssa_raked_mean_raster,
                                     cssa_raked_mean_raster))
save_post_est(rf,'csv','rf')
save_post_est(m,'raster','mean_raster')
save_post_est(m_raked,'raster','mean_raked_raster')


# put adm1 and adm2 means in a csv
for(r in Regions){
  c=get(paste0('cond_sim_raked_adm1_',r))
  c=c[as.numeric(name)%in%get_gaul_codes(r),]
}

ad1 = do.call(rbind,
              mget(grep('cond_sim_raked_adm1_*', ls(), value = TRUE)))
ad2 = do.call(rbind,
              mget(grep('cond_sim_raked_adm2_*', ls(), value = TRUE)))

save_post_est(ad1,'csv','child_adm1')
save_post_est(ad2,'csv','child_adm2')


## ~~~~~~~~~~~~~~~~~~~~~~~~ Shiny Diagnostics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~http://mbg-viz.duckdns.org:3456/~~~~~~~~~~~~~~~~~~~~~~~

# Combine model image history from each reason (these Shiny functions assume a single model image)
all_fixed_effects <- paste(fixed_effects, gbd_fixed_effects, mbg_fixed_effects, sep=" + ")
combine_region_image_history(indicator = indicator,
                             indicator_group = indicator_group,
                             run_date = run_date,
                             fixed_effects = all_fixed_effects)

# Load combined image history
load('<<< FILEPATH REDACTED >>>')

# Load all Africa polygon for plotting
load_simple_polygon(gaul_list = get_gaul_codes('africa'), buffer = 0.4, subset_only = TRUE)

# Plot Shiny stuff
shiny_data_and_preds(gaul_list = get_gaul_codes('africa'),
                     run_date = run_date,
                     indicator = indicator,
                     
                     indicator_group = indicator_group,
                     pred_file = 'has_lri_mean_raster.tif',
                     layer_name = 'has_lri_mean_raster.') 

shiny_raked(gaul_list = get_gaul_codes('africa'),
            run_date = run_date,
            indicator = indicator,
            indicator_group = indicator_group,
            pred_file = 'has_lri_mean_raked_raster.tif',
            layer_name = 'has_lri_mean_raked_raster.') 

shiny_cov_layers(fixed_effects = all_fixed_effects,
                 gaul_list = gaul_list,
                 run_date = run_date,
                 indicator = indicator,
                 indicator_group = indicator_group)

shiny_raking_map()

shiny_data_scatter(df = df,
                   run_date = run_date,
                   indicator = indicator,
                   indicator_group = indicator_group,
                   year_var = 'original_year')


## ~~~~~~~~~~~~~~~~~~~~~~~~ Mark model best ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~ Also would need to run Shiny functions on run_date = "best" ~~~

# save_best_model(indicator_group = indicator_group,
#                 indicator = indicator,
#                 run_date = '2016_10_16_19_27_10',
#                 pred_file = 'has_lri_mean_raster.tif')

##############################################################################
##############################################################################
##############################################################################
