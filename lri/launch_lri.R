#####################################################################

# MBG launch script: 1) launch MBG models
                   # 2) rake to incidence, prevalence, mortality

#####################################################################

#--------------------------------------------------------------------
# (1) Launch MBG models
#--------------------------------------------------------------------

# Setup -------------------------------------------------------------

## set arguments
user                <- commandArgs()[4]
use_run_date        <- commandArgs()[5]
indicator_group     <- commandArgs()[6]
indicator           <- commandArgs()[7]
core_repo           <- commandArgs()[8]
indicator_repo      <- commandArgs()[9]
config_par          <- commandArgs()[10]
cov_par             <- commandArgs()[11]
gbm_par             <- commandArgs()[12]
region_specific_cov <- commandArgs()[13]
Regions             <- unlist(strsplit(commandArgs()[14], '~'))
go_to_inla          <- commandArgs()[15]
rake_measures       <- unlist(strsplit(commandArgs()[16], '~'))
agg_measures        <- unlist(strsplit(commandArgs()[17], '~'))
skip_to_rake        <- commandArgs()[18]

## drive locations
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## source any functions with custom edits from lbd_custom folder
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
if (is.na(use_run_date) | use_run_date == 'NA'){
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- use_run_date
}

## Create output folder with the run_date
outputdir <- '<<<< FILEPATH REDACTED >>>>'

dir.create(outputdir)

## Make sure year object is in the correct format
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

## If running individual countries make sure all country FEs and REs off
if (individual_countries) {
  use_child_country_fes <- FALSE
  use_inla_country_fes  <- FALSE
  use_inla_country_res  <- FALSE
}

## Set optimized BRT parameters if using GBM
if (stacked_fixed_effects %like% 'gbm') {
  gbm_params <- fread('<<<< FILEPATH REDACTED >>>>')
} else {
  gbm_cv <- NA
  gbm_tc <- NA
  gbm_bf <- NA
}

## Record model parameters in a google sheet model tracker
if (is.na(use_run_date) | use_run_date == 'NA'){
  library('googlesheets')
  region_as_char <- paste(Regions, sep=" ", collapse=" + ")
  if (stacked_fixed_effects %like% 'gbm') {
    gbm_bf <- paste(gbm_params$gbm_bf, sep=" ", collapse=" + ")
    gbm_cv <- paste(gbm_params$gbm_cv, sep=" ", collapse=" + ")
    gbm_tc <- paste(gbm_params$gbm_tc, sep=" ", collapse=" + ")
  } else {
    gbm_bf <- 'NA'
    gbm_cv <- 'NA'
    gbm_tc <- 'NA'
    gbm_par <- 'NA'
  }
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
}else if(is.null(gbm_par)){
  gbm_par <- 'NA'
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

if (skip_to_rake == FALSE){
  ## loop over strata, save images and submit qsubs
  
  ## set region-specific covariates if needed
  for(i in 1:nrow(loopvars)){
    
    if (region_specific_cov){
      reg <- loopvars[i,1]
      cov_par_reg <- paste0('covs_', cov_par, '_', reg)
      
      config <- set_up_config(repo = '<<<< FILEPATH REDACTED >>>>',
                              core_repo = core_repo,
                              indicator_group,
                              indicator,
                              config_name = paste0('config_', config_par),
                              covs_name = paste0('covs_', cov_par))
      
      ##double check pop
      if (pop_measure != 'a0004t'){
        print('check population measure!')
      }
      
      message(paste('Used region-specific covariate file', cov_par_reg, 'for region', reg))
    }
    
    if (gbm_par != 'NA'){
      reg <- loopvars[i,1]
      gbm_tc <- gbm_params[region == reg, gbm_tc]
      gbm_lr <- gbm_params[region == reg, gbm_lr]
      gbm_bf <- gbm_params[region == reg, gbm_bf]
      gbm_nminobs <- gbm_params[region == reg, gbm_nminobs]
      gbm_ntrees <- gbm_params[region == reg, gbm_ntrees]
      gbm_cv <- gbm_params[region == reg, gbm_cv]
      
      message(paste0('Used region-specific gbm parameters ', 'gbm_params_', gbm_par, '.csv for region ', reg))
    }
    
    # set cores by region
    region_cores <- 8
    reg <- loopvars[i,1]
    if(reg == 'dia_central_asia') region_cores <- 4
    if(reg == 'sssa' | reg == 'dia_se_asia' | reg == 'dia_sssa') region_cores <- 6
    if(reg == 'dia_malay' | reg == 'dia_name' | reg == 'name' | reg == 'dia_afr_horn') region_cores <- 10
    if(reg == 'dia_chn_mng' | reg == 'dia_wssa' | reg =='dia_south_asia' | reg == 'wssa') region_cores <- 12
    if(reg == 'dia_s_america') region_cores <- 16
    
    # convert cores approximately to memory and run time
    region_rt <- '06:00:00:00'
    if (region_cores < 9) region_rt <- '03:00:00:00'
    if (region_cores < 6) region_rt <- '01:12:00:00'
    if (region_cores > 9) region_rt <- '16:00:00:00'
    region_mem <- region_cores*15
    
    if (go_to_inla == TRUE){
      skiptoinla <- TRUE
      skiptoinla_from_rundate <- run_date
    }
    
    #message age, region, holdout
    message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))
    
    # make a qsub string
    qsub <- make_qsub_share(age           = loopvars[i,2],
                            reg           = as.character(loopvars[i,1]),
                            holdout       = loopvars[i,3],
                            test          = F,
                            indic         = indicator,
                            user          = indicator_group,
                            rd            = run_date,
                            saveimage     = TRUE,
                            use_c2_nodes  = TRUE,
                            memory        = region_mem,
                            cores         = region_cores,
                            run_time      = region_rt,
                            proj          = proj,
                            geo_nodes     = as.logical(use_geos_nodes),
                            corerepo      = core_repo,
                            code_path          = '<<<< FILEPATH REDACTED >>>>', # custom parallel model
                            addl_job_name = paste0(indicator, '_parallel'),
                            singularity   = 'default')
    # submit job
    system(qsub)
  }
  
  waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)
  
} else if (skip_to_rake == TRUE){
  message(paste0('Skipping to raking without launching new models for run_date ', run_date))
}

#--------------------------------------------------------------------
# (2) Rake to GBD
#--------------------------------------------------------------------

## Launcher script for rake_lri_cell_pred

# set variables
nodes <- 'geos'
user <- '<<<< USERNAME REDACTED >>>>'

rake_rt = '10:00:00'

# Define nodes
if (nodes == 'geos') {
  proj <- '-P proj_geo_nodes'
  r_shell <- '<<<< FILEPATH REDACTED >>>>'
} else {
  proj <- '-P proj_geospatial'
  r_shell <- '<<<< FILEPATH REDACTED >>>>'
}

# Launch child scripts
for (reg in Regions) {
  for (raketo in rake_measures){
    rake_mem = '100G'
    if (reg == 'dia_s_america'){
      rake_mem = '150G'
    }
    
    jname <- paste0(reg, '_', raketo, '_rake')
    sys.sub <- paste0('qsub -l m_mem_free=', rake_mem, ' -l fthread=1 -v sing_image=default -q geospatial.q -l h_rt=', rake_rt, ' ', proj,' -e <<<< FILEPATH REDACTED >>>>',
                      ' -o <<<< FILEPATH REDACTED >>>> -N ', jname)
    script <- '<<<< FILEPATH REDACTED >>>>/lri/rake_lri_cell_pred.R'
    args <- paste(raketo, reg, outputdir, makeholdouts, 2000, indicator, modeling_shapefile_version, raking_shapefile_version)
    system(paste(sys.sub, r_shell, script, args))
  }
}

print('recall that year_list is set in rake_lri_cell_pred.R to c(2000:2017)')