##############################################################################
## Launch script to run BRT optimization
##############################################################################


## Setup -------------------------------------------------------------------------

# set node preference
use_geos_nodes <- TRUE
user <- Sys.info()['user']
if (use_geos_nodes) {
  proj_arg <- 'proj_geo_nodes'
} else {
  proj_arg <- 'proj_geospatial'
}

# Use hf_inla
r_shell <- paste0('<<<< FILEPATH REDACTED >>>>/ort/mbg_central/share_scripts/shell_sing.sh')

# indicators to be launched
indicators <- c('ors', 'rhf', 'ors_or_rhf')

# set cov version
cov_version <- 'standard_inc_mbg'

# bounds versions to test
bounds_versions <- c(1:16)

# train and bag fractions to test
cv_fold <- 3
bag_fraction <- 0.5


## BRT optimizer script -------------------------------------------------------------------------

for (bv in bounds_versions) {

  for (ind in indis) {
    
    # set file addin
    file_addin <- FALSE
    
    # set regions
    regions <- c('dia_afr_horn-eth-yem', 'ETH', 'YEM', 'dia_name', 'dia_sssa', 
                 'dia_mcaca', 'dia_central_asia', 'MNG', 
                 'dia_se_asia', 'dia_malay', 'dia_mid_east',
                 'ZWE', 'KEN', 'NGA', 'COD', 'IND', 'PAK',
                 'dia_essa-zwe-ken', 'dia_cssa-cod', 'dia_south_asia-ind-pak',
                 'dia_s_america_n', 'dia_s_america_s')
    
    for (r in regions) {
      
      # get arguments
      jname <- paste0(ind, '_', r, '_optim')
      mymem <- '2G'
      sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
                        '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' '),
                        '-l fthread=1 -l h_rt=02:00:00:00 -v sing_image=default -N ', jname, ' ')
      # launch script name to qsub
      script <- '<<<< FILEPATH REDACTED >>>>/ort/gbm_optim/runGBM.R'
      jobnum <- paste0(ind, '_', r)
      opt_type <- 'gp'
      lrnr_type <- 'brt'
      bounds_version <- bv
      experiment_version <- paste0('test', bounds_version, '_', Sys.Date(), ifelse(file_addin == FALSE, '', paste0('_', file_addin)))
      args <- paste0(jobnum, ' ', opt_type, ' ', lrnr_type, ' ', bounds_version, ' ', 
                     experiment_version, ' ', ind, ' ', cv_fold, ' ', bag_fraction, ' ', file_addin)
      system(paste(sys.sub, r_shell, script, args))
  
    }
    
    Sys.sleep(1)
  }

}