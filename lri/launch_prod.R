###############################################################################
## MBG launch script
##
## Example inital launch script for LRI
##
###############################################################################
#comment

## Setup -------------------------------------------------------------------------

## clear environment
rm(list=ls())

## Set repo location, indicator group, and some arguments
user            <- Sys.info()['user']
use_run_date    <- NULL
indicator_group <- 'lri'
indicator       <- 'has_lri'
core_repo       <- '<<<< FILEPATH REDACTED >>>>' # use the core repo so you don't have to worry about keeping it up to date
indicator_repo  <- '<<<< FILEPATH REDACTED >>>>' # source the configs from your indicator repo
config_par      <- '06' # you can change these for your custom configs (for example, you may want change to the GADM shapefiles using modeling_shapefile_version 2018_08_01)
cov_par         <- '03' # you can change these for your custom covariates (for example, you may want to add other MBG covariates and add GBD covariates)
gbm_par         <- '02'
region_specific_cov <- T
Regions         <- c('essa', 'cssa', 'wssa', 'sssa', 'name') # can be all regions for a global model, or can be an individual country, e.g. 'PER'
stage_1         <- TRUE #sets region for post-processing to 'africa'
stage_2         <- FALSE #sets region for post-processing to all stage 2 locations
rake_measures   <- c('incidence','prevalence','mortality') #what measures to rake to
agg_measures    <- rake_measures #what measures to aggregate over
rakecores       <- 8 # how many slots to request for raking, for each region

message(indicator)


## drive locations
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#source any functions with custom edits from sharedir
source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/misc_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/post_estimation_functions.R')

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')

## Read config file and save all parameters in memory
config <- load_config(repo            = indicator_repo,
                      indicator_group = '',
                      indicator       = '',
                      config_name     = paste0('config_', config_par),
                      covs_name       = paste0('covs_', cov_par))

## Ensure you have defined all necessary settings in your config
check_config()

##set pop
pop_measure <- 'a0004t'

## Set project
proj <- ifelse(as.logical(use_geos_nodes), 'proj_geo_nodes', 'proj_geospatial')

## Create run date in correct format
if (is.null(use_run_date)){
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- use_run_date
}

## Create output folder with the run_date
outputdir      <- '<<<< FILEPATH REDACTED >>>>'
dir.create(outputdir)

## Make sure year object is in the correct format
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

## If running individual countries make sure all country FEs and REs off
if (individual_countries) {
  use_child_country_fes <- FALSE
  use_inla_country_fes  <- FALSE
  use_inla_country_res  <- FALSE
}

## Set optimized BRT parameters
if (stacked_fixed_effects %like% 'gbm') {
  gbm_params <- fread('<<<< FILEPATH REDACTED >>>>', stringsAsFactors = F)
} else {
  gbm_cv <- NA
  gbm_tc <- NA
  gbm_bf <- NA
}

## Record model parameters in a google sheet model tracker
library('googlesheets')
region_as_char <- paste(Regions, sep=" ", collapse=" + ")
gbm_bf <- paste(gbm_params$gbm_bf, sep=" ", collapse=" + ")
gbm_cv <- paste(gbm_params$gbm_cv, sep=" ", collapse=" + ")
gbm_tc <- paste(gbm_params$gbm_tc, sep=" ", collapse=" + ")
model_params <- c(indicator_group, indicator, run_date, outputdir, 
                  individual_countries, region_as_char, 
                  config_par, cov_par, gbm_par, region_specific_cov, fixed_effects, gbd_fixed_effects,
                  stacked_fixed_effects, makeholdouts, use_inla_country_fes, gbm_cv, gbm_bf, gbm_tc,
                  mesh_t_knots, rho_prior, nugget_prior, ctry_re_prior, modeling_shapefile_version)

ttt <- readRDS('~/repos/reference/gs_authorization_token.RData')
gs_auth(token = ttt)
lri_tracker <- gs_title('LRI model tracker')
gs_add_row(lri_tracker, ws = 'model runs', input = model_params)
rm(lri_tracker)

# ^^^ I find having a model tracker very helpful. I can help you set this up.


## Make holdouts -------------------------------------------------------------------------

if(makeholdouts){
  message('Making holdouts')
  
  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        simple      = NULL,
                        removeyemen = TRUE,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share))
  
  # add in location information
  df <- merge_with_ihme_loc(df)
  
  # make a list of dfs for each region, with 5 qt folds identified in each
  stratum_ho <- make_folds(data       = df,
                           n_folds    = as.numeric(n_ho_folds),
                           spat_strat = 'qt',
                           temp_strat = 'prop',
                           strat_cols = 'region',
                           ts         = as.numeric(ho_ts),
                           mb         = as.numeric(ho_mb))
}


## Launch parallel script -------------------------------------------------------------------------

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(makeholdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)

## loop over them, save images and submit qsubs

# set region-specific covariates if needed
for(i in 1:nrow(loopvars)){
  
  if (region_specific_cov){
    reg <- loopvars[i,1]
    cov_par_reg <- paste0('covs_', cov_par, '_', reg)
    
    ## Read config file and save all parameters in memory
    config <- load_config(repo            = indicator_repo,
                          indicator_group = '',
                          indicator       = '',
                          config_name     = paste0('config_', config_par),
                          covs_name       = cov_par_reg)
    
    ## Ensure you have defined all necessary settings in your config
    check_config()
    
    ##set pop
    pop_measure <- 'a0004t'
    
    message(paste('Used region-specific covariate file', cov_par_reg, 'for region', reg))
  }
  
  if (!is.null(gbm_par)){
    reg <- loopvars[i,1]
    gbm_tc <- gbm_params[region == reg, gbm_tc]
    gbm_lr <- gbm_params[region == reg, gbm_lr]
    gbm_bf <- gbm_params[region == reg, gbm_bf]
    gbm_nminobs <- gbm_params[region == reg, gbm_nminobs]
    gbm_ntrees <- gbm_params[region == reg, gbm_ntrees]
    gbm_cv <- gbm_params[region == reg, gbm_cv]
    
    message(paste0('Used region-specific gbm parameters ', 'gbm_params_', gbm_par, '.csv for region ', reg))
  }
  message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))
  # make a qsub string
  qsub <- make_qsub_share(age           = loopvars[i,2],
                          reg           = as.character(loopvars[i,1]),
                          holdout       = loopvars[i,3],
                          test          = F,
                          indic         = indicator,
                          saveimage     = TRUE,
                          use_c2_nodes  = TRUE,
                          memory        = ifelse(individual_countries, 10, 50), # this is conservative, will need to profile your jobs by region
                          cores         = ifelse(as.logical(use_geos_nodes), 10, as.numeric(slots)),
                          proj          = proj,
                          geo_nodes     = as.logical(use_geos_nodes),
                          corerepo      = core_repo,
                          code          = NULL, # eventually you may want to specifiy your own custom parallel script (for example, if you want to output means maps or directly go to post-estimate and aggregation)
                          addl_job_name = paste0(indicator, '_parallel'),
                          singularity   = 'default')
  # submit job
  system(qsub)
  
}

## check to make sure models are done before continuing
waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

##############################################################################
## Summarize model results
##############################################################################

clean_model_results_table()

# Stop if individual countries (no need for post-est)
if(as.logical(individual_countries) == F) {
  
  ###############################################################################
  ## Post-Estimation
  ###############################################################################
  
  ## Save strata for Shiny to use in producing aggregated fit statistics
  strata <- unique(as.character(loopvars[,1]))
  dir.create(paste0(sharedir, '/fit_stats'))
  save(strata, file = paste0(sharedir, '/fit_stats/strata.RData'))
  
  
  ## Load GBD Estimates for this indicator which will be used in raking
  ## set region
  if (stage_1){
    region <- 'africa'
  } else if (stage_2){
    region <- 'stage_2'
  }
  
  gbd <- get_gbd_estimates(gbd_name = 322,
                           region = 'africa',
                           measure_id = 5,
                           age_group_id = 1,
                           metric_id = 3,
                           year_ids = year_list,
                           shapefile_version = raking_shapefile_version,
                           rake_subnational = subnational_raking, 
                           gbd_round_id = 4) #TODO should be 5
  
  ## Prepare for parallel post-estimation - save file with objs to re-load in child processes
  prep_postest(indicator = indicator,
               indicator_group = indicator_group,
               run_date = run_date,
               save_objs = c("core_repo", "gbd", "year_list", "summstats",
                             "rake_transform", "pop_measure"))
  
  ## Parallelized post-estimation over region
  postest_script <- "postest_script"
  
  for (s in strata) {
    qsub <- make_qsub_postest(code = postest_script,
                              stratum = s,
                              log_location = 'sharedir',
                              memory        = 10,
                              singularity   = "default",
                              subnat_raking = subnational_raking,
                              modeling_shapefile_version = modeling_shapefile_version,
                              raking_shapefile_version = raking_shapefile_version, 
                              geo_nodes     = as.logical(use_geos_nodes),
                              cores         = 10)
    system(qsub)
  }
  
  ## check to make sure post-est done before continuing
  waitformodelstofinish(lv = cbind(strata, 0), sleeptime=60)
  
  ## Combine post est stuff across regions and save needed outputs
  post_load_combine_save(summstats = summstats)
  
  # Clean up / delete unnecessary files
  clean_after_postest(indicator             = indicator,
                      indicator_group       = indicator_group,
                      run_date              = run_date,
                      strata                = strata,
                      delete_region_rasters = F)
  
  ###############################################################################
  ## Launch model diagnostics script for shiny tool
  ###############################################################################
  make_model_diagnostics(indic      = indicator,
                         ig         = indicator_group,
                         rd         = run_date,
                         geo_nodes  = TRUE,
                         cores      = 10)
  ###############################################################################
  ## Rake to GBD
  ###############################################################################
  ## Launcher script for rake_lri_cell_pred
  #for each region, submit a separate job:
  # set variables
  nodes <- 'geos'
  
  # set repo
  setwd('<<<< FILEPATH REDACTED >>>>')
  
  # Define nodes
  if (nodes == 'geos') {
    proj <- '-P proj_geo_nodes -l gn=TRUE'
    r_shell <- '<<<< FILEPATH REDACTED >>>>'
  } else {
    proj <- '-P proj_geospatial'
    r_shell <- '<<<< FILEPATH REDACTED >>>>'
  }
  
  
  # Launch child scripts
  for (reg in Regions) {
    for (raketo in rake_measures){
      jname <- paste0(reg, '_', raketo, '_rake')
      sys.sub <- paste0('qsub ',proj,paste0(' -e <<<< FILEPATH REDACTED >>>>',
                                            ' -o <<<< FILEPATH REDACTED >>>> '),
                        '-cwd -N ', jname, ' ', '-pe multi_slot ', rakecores, ' ')
      script <- 'rake_lri_cell_pred.R'
      args <- paste(raketo, reg, outputdir, makeholdouts, 2000, indicator, modeling_shapefile_version, raking_shapefile_version)
      system(paste(sys.sub, r_shell, 1, script, args))
    }
  }
  
  print('recall that year_list is set in rake_lri_cell_pred.R to c(2000:2017)')
  ###############################################################################
  ## Aggregate to admin2, admin1, and national levels
  ###############################################################################
  #aggregate raked
  submit_aggregation_script (indicator       = indicator,
                             indicator_group = indicator_group,
                             run_date        = run_date,
                             raked           = T,
                             pop_measure     = pop_measure,
                             overwrite       = T,
                             ages            = 0, # Note: can take vector of ages
                             holdouts        = 0,
                             regions         = strata,
                             corerepo        = '<<<< FILEPATH REDACTED >>>>',
                             log_dir         = paste0(sharedir, "/output/", run_date),
                             geo_nodes       = T,
                             singularity     = "default",
                             slots           = 10,
                             modeling_shapefile_version = modeling_shapefile_version,
                             raking_shapefile_version = raking_shapefile_version,
                             measures        = agg_measures)
  
  #aggregate unraked
  submit_aggregation_script (indicator       = indicator,
                             indicator_group = indicator_group,
                             run_date        = run_date,
                             raked           = F,
                             pop_measure     = pop_measure,
                             overwrite       = T,
                             ages            = 0, # Note: can take vector of ages
                             holdouts        = 0,
                             regions         = strata,
                             corerepo        = '<<<< FILEPATH REDACTED >>>>',
                             log_dir         = paste0(sharedir, "/output/", run_date),
                             geo_nodes       = T,
                             singularity     = "default",
                             slots           = 10,
                             modeling_shapefile_version = modeling_shapefile_version,
                             raking_shapefile_version = raking_shapefile_version,
                             measures        = 'prevalence')
  
  waitforaggregation(rd = run_date, indic = indicator, ig = indicator_group,
                     ages     = 0,
                     regions  = strata,
                     holdouts = 0,
                     raked    = c(T, F),
                     measure  = agg_measures)
  
  #unraked combine and summarize
  combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                      ages     = 0,
                      Regions  = strata,
                      holdouts = 0,
                      raked    = F,
                      measures = 'prevalence')
  
  summarize_admins(summstats = c("mean", "upper", "lower", "cirange"),
                   ad_levels = c(0,1,2),
                   raked     = F,
                   measures  = 'prevalence')
  
  #raked combine and summarize
  combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                      ages     = 0,
                      Regions  = strata,
                      holdouts = 0,
                      raked    = T,
                      measures = agg_measures)
  
  summarize_admins(summstats = c("mean", "upper", "lower", "cirange"),
                   ad_levels = c(0,1,2),
                   raked     = T,
                   measures  = agg_measures)
  
  
  
  # make scatter of old vs. new model
  source('<<<< FILEPATH REDACTED >>>>/lri/adm_compare.R')
  
  #unraked
  adm_compare(indicator <- indicator,
              indicator_group <- indicator_group,
              rake_list <- FALSE,
              admin_lvls <- c(0:2),
              measure_list <- 'prevalence',
              age <- 0,
              holdout <- 0,
              run_date_1 <- '2019_04_01_15_28_31',
              run_date_2 <- run_date,
              region_list <- strata,
              collapse_years <- T)
  
  #raked
  adm_compare(indicator <- indicator,
              indicator_group <- indicator_group,
              rake_list <- TRUE,
              admin_lvls <- c(0:2),
              measure_list <- agg_measures,
              age <- 0,
              holdout <- 0,
              run_date_1 <- '2019_04_01_15_28_31',
              run_date_2 <- run_date,
              region_list <- strata,
              collapse_years <- T)
  
  # make line plots
  source('<<<< FILEPATH REDACTED >>>>/lri/stacker_admin_time_series.R')
  
  plot_stackers_by_adm01(indicator = indicator,
                         indicator_group = indicator_group,
                         run_date = run_date,
                         regions = strata,
                         draws = T,
                         raked = T,
                         credible_interval = 0.95,
                         N_breaks = c(0, 10, 50, 100, 500, 1000, 2000, 4000),
                         admin_data = NULL,
                         data_tag = 'lriprevalence',
                         measure <- 'prevalence',
                         rm_yemen <- T)
  
  message("After summarize_admins - about to save CSV")
  
  # Combine csv files
  csvs <- list.files(paste0(sharedir, '/output/', run_date, '/'),
                     pattern = "input_data(.*).csv",
                     full.names = T)
  
  csv_master <- lapply(csvs, fread) %>%
    rbindlist %>%
    subset(., select = names(.) != "V1")
  write.csv(csv_master, file=paste0(sharedir, '/output/', run_date, '/input_data.csv'))
  
  ###############################################################################
  ## Create AROC objects & do projections
  ###############################################################################
  
  # note! you need to use raking_shapefile_version is using raked results
  make_aroc(ind_gp           = indicator_group,
            ind              = indicator,
            rd               = run_date,
            matrix_pred_name = NULL,
            type             = c("cell", "admin"),
            measure          = "prevalence",
            year_list        = year_list,
            uselogit         = FALSE,
            raked            = TRUE,
            weighting_res    = 'domain',
            weighting_type   = 'exponential',
            pow              = 1,
            input_data = read.csv('<<<< FILEPATH REDACTED >>>>'),
            mult_emp_exp     = FALSE,
            extra_file_tag = "_exp_domain",
            shapefile_version = raking_shapefile_version)
  
  # note! you need to use raking_shapefile_version is using raked results
  make_proj(ind_gp     = indicator_group,
            ind        = indicator,
            rd         = run_date,
            type       = c("cell", "admin"),
            proj_years = c(2020, 2025, 2030),
            measure    = "prevalence",
            skip_cols  = NULL,
            year_list  = c(2000:2015),
            uselogit   = FALSE,
            raked = TRUE, 
            extra_file_tag = "_exp_domain", #TODO
            shapefile_version = raking_shapefile_version)
  
  ###############################################################################
  # Look at performance against goals
  ###############################################################################
  
  # Define goals: start by initializing goal object
  goals <- add_goal(target_year = 2025,
                    target = 0.03,
                    target_type = "less",
                    abs_rel = "absolute",
                    pred_type = c("cell", "admin"))
  
  # Add goals to existing goal object by specifying goal_obj
  goals <- add_goal(goal_obj = goals,
                    target_year = 2025,
                    baseline_year = 2010,
                    target = 0.75,
                    target_type = "less",
                    abs_rel = "relative",
                    pred_type = c("cell", "admin"))
  
  # Run comparisons
  #TODO: this is wrong, write in options for correct!
  compare_to_target(ind_gp = indicator_group,
                    ind = indicator,
                    rd = run_date,
                    goal_obj = goals,
                    measure = "prevalence",
                    year_list = c(2000:2015),
                    uselogit = FALSE,
                    shapefile_version = modeling_shapefile_version)
  
  ###############################################################################
  # Make summary metrics
  ###############################################################################
  
  # Get in and out of sample draws
  run_in_oos <- get_is_oos_draws(ind_gp = indicator_group,
                                 ind = indicator,
                                 rd = run_date,
                                 ind_fm = 'binomial',
                                 model_domain = region,
                                 age = 0,
                                 nperiod = length(year_list),
                                 yrs = year_list,
                                 get.oos = as.logical(makeholdouts),
                                 write.to.file = TRUE,
                                 shapefile_version = modeling_shapefile_version)
  
  ## set out_dir
  out_dir <- paste0(sharedir, "/output/", run_date, "/summary_metrics/")
  dir.create(out_dir, recursive = T, showWarnings = F)
  
  ## for admin0
  draws.df <- fread('<<<< FILEPATH REDACTED >>>>')
  
  country.pvtable <- get_pv_table(d = draws.df,
                                  indicator_group = indicator_group,
                                  rd = run_date,
                                  indicator=indicator,
                                  aggregate_on='country',
                                  draws = as.numeric(samples),
                                  out.dir = out_dir)
  
  write.csv(country.pvtable,
            file = '<<<< FILEPATH REDACTED >>>>')
  
  ad1.pvtable <- get_pv_table(d = draws.df,
                              indicator_group = indicator_group,
                              rd = run_date,
                              indicator=indicator,
                              aggregate_on='ad1',
                              draws = as.numeric(samples),
                              out.dir = out_dir)
  write.csv(ad1.pvtable,
            file = '<<<< FILEPATH REDACTED >>>>')
  
  ad2.pvtable <- get_pv_table(d = draws.df,
                              indicator_group = indicator_group,
                              rd = run_date,
                              indicator=indicator,
                              aggregate_on='ad2',
                              draws = as.numeric(samples),
                              out.dir = out_dir)
  write.csv(ad2.pvtable,
            file = '<<<< FILEPATH REDACTED >>>>')
  
  
}