###############################################################################
## MBG launch script
## for fair cluster
##
###############################################################################

## Setup -------------------------------------------------------------------------

## clear environment
rm(list=ls())

## Set repo location, indicator group, and other arguments
user            <- Sys.info()['user']
use_run_date    <- '2019_06_25_12_40_32'
indicator_group <- 'lri'
indicator       <- 'has_lri'
core_repo       <- '<<<< FILEPATH REDACTED >>>>'
indicator_repo  <- '<<<< FILEPATH REDACTED >>>>'
config_par      <- '11'
cov_par         <- '04'
gbm_par         <- '02'
region_specific_cov <- FALSE
Regions         <- c('wssa','essa','cssa','sssa','name')
go_to_inla      <- FALSE  #skips to INLA, reading in stackers from saved model image history
stage_1         <- FALSE #sets region for post-processing to 'africa'
stage_2         <- TRUE #sets region for post-processing to all stage 2 locations
rake_measures   <- c('incidence','prevalence','mortality') #what measures to rake to
agg_measures    <- rake_measures #what measures to aggregate over

message(indicator)


## drive locations
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#source any functions with custom edits from lbd_custom folder
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/misc_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/post_estimation_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/stacking_functions.R')

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')

## Read config file and save all parameters in memory
source(paste0(core_repo, '/mbg_central/LBDCore/R/set_up_config.R'))
load(paste0(core_repo, '/mbg_central/LBDCore/data/config_values.rda'))
load(paste0(core_repo, '/mbg_central/LBDCore/data/config_tests.rda'))
setwd(paste0(core_repo, '/mbg_central/LBDCore'))

config <- set_up_config(repo = '<<<< FILEPATH REDACTED >>>>',
                        core_repo = core_repo,
                        indicator_group,
                        indicator,
                        config_name = paste0('config_', config_par),
                        covs_name = paste0('covs_', cov_par))

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
if (is.null(use_run_date)){
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
  
  ttt <- readRDS('<<<< FILEPATH REDACTED >>>>')
  gs_auth(token = ttt)
  lri_tracker <- gs_title('LRI model tracker')
  gs_add_row(lri_tracker, ws = 'model runs', input = model_params)
  rm(lri_tracker)
}

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

  ## Save strata for Shiny to use in producing aggregated fit statistics
  strata <- unique(as.character(loopvars[,1]))

  
  #unraked combine and summarize
  combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                      ages     = 0,
                      Regions  = strata,
                      holdouts = 0,
                      raked    = F,
                      measures = 'prevalence',
                      check_for_dupes = T)
  
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
                      measures = agg_measures,
                      check_for_dupes = T)
  
  summarize_admins(summstats = c("mean", "upper", "lower", "cirange"),
                   ad_levels = c(0,1,2),
                   raked     = T,
                   measures  = agg_measures)

