#######################################################################################
## MBG core diagnostics launcher script
##
## Launches diagnostics.R:
## - cleans up and summarizes aggregations, outputs csvs in pred_derivatives 
## - aggregates data and stackers
## - plots stacker line plots
## - plots INLA hyperparameters
## - creates in-sample and out-of-sample fit statistics
#######################################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'ort'

# set cluster arguments
use_geos_nodes  <- TRUE
proj_arg        <- ifelse(use_geos_nodes, 'proj_geo_nodes_dia', 'proj_geospatial_dia')

# list indicators
indicators <- c('ors_or_rhf', 'rhf', 'ors')

# indicate run date
run_date <- '2019_11_07_13_50_11'

# set measure and configs
measure <- 'prevalence'
config_par <- 'config_ort_best'
config_file <- 'ors/3_modeling/'
cov_par <- 'covs_ort_standard'
cov_file <- config_file


## Run diagnostic functions -------------------------------------------------------------------------

# run diagnostics for each indicator
for (indicator in indicators) {
  message(indicator)
    
  # set indicator arguments
  jname <- paste0(indicator, '_diagnostics')
  mymem <- '100G'
  
  # set up qsub
  sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
                    '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                    '-l fthread=1 -l h_rt=00:12:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
  r_shell <- paste0('<<<< FILEPATH REDACTED >>>>/ort/mbg_central/share_scripts/shell_sing.sh')
  script <- paste0('<<<< FILEPATH REDACTED >>>>/ort/post_estimation/diagnostics.R')
  args <- paste(user, repo, indicator_group, indicator, config_par, config_file, cov_par, cov_file,
                run_date, measure)
  
  # run launch script
  system(paste(sys.sub, r_shell, script, args))
  
}
