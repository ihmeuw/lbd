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
indicator_group <- 'lri'
indicator       <- 'has_lri'
file_addin      <- FALSE #append a suffix to the end of the saved csv

# set cluster arguments
proj_arg <- 'proj_geospatial'
memory <- '10G'
run_time <- '02:00:00'

# set config and covariate files
config_version   <- '11'
cov_version <- '04'
region_specific_cov <- T

# set whether running for individual countries
individual_countries <- FALSE

# list all regions or countries
regions <- c('sssa','wssa','essa','name','cssa')


## Run covariate data scripts -------------------------------------------------------------------------

for (r in regions) {
  
  # set covariate file name
  if (region_specific_cov){
    cov_file <- paste0('covs_', cov_version, '_', r)
  } else {
    cov_file <- paste0('covs_', cov_version)
  }
  
  # set specific arguments
  region          <- r
  jname           <- paste0(r, '_brt_data')
  
  # set up qsub
  sys.sub <- paste0('qsub -l m_mem_free=', memory, 
                    ' -l fthread=1 -l h_rt=', run_time,  
                    ' -v sing_image=default -q all.q -l archive=TRUE -P ', proj_arg, 
                    ' -N ', jname, 
                    ' <<<< FILEPATH REDACTED >>>> <<<< FILEPATH REDACTED >>>>')
  args <- paste(region, user, repo, indicator_group, indicator, config_version, 
                cov_file, individual_countries, file_addin)
  
  # run launch script
  system(paste(sys.sub, args))
  
}

