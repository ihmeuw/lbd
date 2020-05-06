#######################################################################################
## MBG core diagnostics launcher script for diarrhea
##
## Launches diagnostics.R:
## - cleans up and summarizes aggregations as csvs
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
indicators <- c('had_diarrhea')

# indicate run date
run_date <- '2019_09_17_14_12_53'

# covariate and config settings
config_par   <- 'config_dia_best'
config_file <- 'had_diarrhea/'
cov_par      <- 'covs_dia_standard_inc_mbg'
cov_file <- config_file


## Run diagnostic functions -------------------------------------------------------------------------

# set measures
measures <- c('prevalence', 'deaths', 'incidence')

# run diagnostics for each measure
for (measure in measures) {
  message(measure)
  
  # set specific arguments
  jname           <- paste0(indicator, '_diagnostics')
  mymem           <- ifelse(measure == 'prevalence', '100G', '30G')
  
  # set up qsub
  sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
                    '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                    '-l fthread=1 -l h_rt=00:12:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
  r_shell <- paste0('<<<< FILEPATH REDACTED >>>>/mbg_central/share_scripts/shell_sing.sh')
  script <- paste0('/<<<< FILEPATH REDACTED >>>>/post_estimation/diagnostics.R')
  args <- paste(user, repo, indicator_group, indicator, config_par, config_file, cov_par, cov_file,
                run_date, measure)
  
  # run launch script
  system(paste(sys.sub, r_shell, script, args))
  
}
