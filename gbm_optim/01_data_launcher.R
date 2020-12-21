##############################################################################
## Launch script to crop covariates and data together for BRT optimization
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- paste0('<<<< FILEPATH REDACTED >>>>')
core_repo       <- repo
indicator_group <- 'ort'
indicators      <- c('any_ors', 'rhf_only', 'no_ort')

# set cluster arguments
use_geos_nodes  <- TRUE
proj_arg        <- 'proj_geo_nodes_dia'

# set whether running for individual countries
individual_countries <- FALSE


## Run covariate data scripts -------------------------------------------------------------------------

for (indicator in indicators) {
  
  # set config and covariate files
  config_par   <- 'config_ort_gam3'
  config_file  <- 'ors/3_modeling/'
  covar_par      <- 'region_specific'
  cov_file     <- config_file
  file_addin <- FALSE
  
  # list all regions or countries
  if (indicator == 'any_ors' | indicator == 'rhf_only' | indicator == 'no_ort') {
    regions <- c('SEN', 'SLE', 'MLI')
  }
  
  for (r in regions) {
    
    # set covariates
    cov_par <- paste0('covs_', indicator, '/covs_ort_', r)
    
    # set specific arguments
    region          <- r
    jname           <- paste0(r, '_brt_data')
    mymem           <- '20G'
    
    # set up qsub
    sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
                      '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' '),
                      '-l fthread=1 -l h_rt=00:04:00:00 -v sing_image=default -N ', jname, ' ')
    r_shell <- paste0('<<<< FILEPATH REDACTED >>>>/shell_sing.sh')
    script <- paste0('<<<< FILEPATH REDACTED >>>>/ort/gbm_optim/save_cov_data.R')
    args <- paste(region, user, repo, indicator_group, indicator, config_par, config_file, 
                  cov_par, cov_file, individual_countries, file_addin)
    
    # run launch script
    system(paste(sys.sub, r_shell, script, args))
    
  }
  
}
