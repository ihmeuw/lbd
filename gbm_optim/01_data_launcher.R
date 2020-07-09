##############################################################################
## Launch script to crop covariates and data together for BRT optimization
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
core_repo       <- repo
indicator_group <- 'ort'
indicators      <- c('ors', 'rhf', 'ors_or_rhf')

# set cluster arguments
use_geos_nodes  <- TRUE
proj_arg        <- 'proj_geo_nodes_dia'

# set cov version
cov_version <- 'standard_inc_mbg'

# set whether running for individual countries
individual_countries <- FALSE


## Run covariate data scripts -------------------------------------------------------------------------

for (indicator in indicators) {
  
  # set config and covariate files
  config_par   <- 'config_ort_best'
  config_file  <- 'ors/3_modeling/'
  covar_par      <- 'region_specific'
  cov_file     <- config_file
  file_addin <- FALSE
  
  # list all regions or countries
  regions <- c('dia_afr_horn-eth-yem', 'ETH', 'YEM', 'dia_name', 'dia_sssa', 
               'dia_mcaca', 'dia_central_asia', 'MNG', 
               'dia_se_asia', 'dia_malay', 'dia_mid_east',
               'ZWE', 'KEN', 'NGA', 'COD', 'IND', 'PAK',
               'dia_essa-zwe-ken', 'dia_cssa-cod', 'dia_south_asia-ind-pak',
               'dia_s_america_n', 'dia_s_america_s')
  
  for (r in regions) {
    
    # set region specific covariates, if desired
    if (covar_par == 'region_specific') cov_par <- paste0('covs_', indicator, '/covs_ort_', r)

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
