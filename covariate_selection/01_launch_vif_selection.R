##############################################################################
## Variance inflation factor (VIF) launcher script
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'ort'

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = repo)

# set cluster arguments
use_geos_nodes  <- TRUE
proj_arg        <- 'proj_geo_nodes_dia'

# indicate whether to use old run date
use_old_run_date <- FALSE
old_run_date_input <- ''

# set whether running individual countries
individual_countries <- FALSE

# indicate threshold parameters
threshold_min <- 2
threshold_max <- 5
threshold_step <- 1

# indicate whether or not to crop covariates (must give .RData file with cropped covariate file if crop_covs = FALSE)
crop_covs <- TRUE

# indicate whether there are any covariates you'd like to make sure we keep, otherwise set to NULL
# note: if you pick more than one, you may not fully remove collinearity
covs_to_keep <- 'edu_mean_stage2'

# list indicators
indics <- c('ors', 'rhf', 'ors_or_rhf')

# list either all regions or one region/country (argument cannot be a vector/list of regions)
regions <- 'all'


## Run vif selection script -------------------------------------------------------------------------

for (region in regions ) {

  for (i in indics) {
    
    # set run date
    if (use_old_run_date == FALSE) {
      run_date <- make_time_stamp(TRUE)
    } else {
      run_date <- old_run_date_input
    }
      
    # set specific arguments
    indicator       <- i
    jname           <- paste0(indicator, '_vif_selection')
    mymem           <- ifelse(crop_covs, '50G', '10G')
    
    # set config and covariate files
    config_par   <- 'config_ort_best'
    config_file  <- 'ors/3_modeling/'
    cov_par      <- 'covs_ort_standard_inc_mbg'
    cov_file     <- config_file
    
    # set cropped covs file (will only be used if crop_covs is FALSE)
    cropped_covs_file <- paste0('<<<< FILEPATH REDACTED >>>>', i, '/covariate_selection_2019_07_16_12_11_48/covs_cropped_to_data_2019-07-16.RData')
    
    # set up qsub
    sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
                      '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                      '-l fthread=1 -l h_rt=01:00:00:00 -v sing_image=default -N ', jname, ' ')
    r_shell <- paste0('<<<< FILEPATH REDACTED >>>>/ort/mbg_central/share_scripts/shell_sing.sh')
    script <- paste0('<<<< FILEPATH REDACTED >>>>/ort/covariate_selection/02_vif_selection_script.R')
    args <- paste(user, repo, indicator_group, indicator, config_par, config_file, cov_par, cov_file, 
                  run_date, region, threshold_min, threshold_max, threshold_step, individual_countries,
                  crop_covs, cropped_covs_file, covs_to_keep)
    
    # run launch script
    system(paste(sys.sub, r_shell, script, args))
    
  }
  
}
