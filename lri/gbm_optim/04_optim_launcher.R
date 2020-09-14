##############################################################################
## Launch script to run BRT optimization
##############################################################################


## Setup -------------------------------------------------------------------------

# set repo
setwd('<<<< FILEPATH REDACTED >>>>')

# set node preference
proj <- 'proj_geospatial'
memory <- '10G'
run_time <- '30:00:00'

# Use hf_inla
hf_inla_shell <- '<<<< FILEPATH REDACTED >>>>'
r_shell <- hf_inla_shell
mkl <- 1

# indicators to be launched
indis <- c('has_lri')

# regions to be launched
regions <- c('wssa','essa','name','sssa','cssa')

# set cov version
cov_version <- '04'
file_addin <- FALSE

# bounds versions to test
bounds_versions <- c(2,4,7,11,14,19,14,29,30,31,32,33,34,35)

# train and bag fractions to test
cv_fold <- 3
bag_fraction <- 0.5


## BRT optimizer script -------------------------------------------------------------------------

for (bv in bounds_versions) {

  for (ind in indis) {
    
    for (r in regions) {
      
      # get arguments
      jname <- paste0(ind, '_', r, '_optim')
      sys.sub <- paste0('qsub -l m_mem_free=', memory, 
                        ' -l fthread=1 -l h_rt=', run_time,  
                        ' -v sing_image=default -q all.q -l archive=TRUE -P ', proj, 
                        ' -N ', jname, 
                        ' ', r_shell,
                        ' ', mkl,
                        '<<<< FILEPATH REDACTED >>>>')
      # launch script name to qsub
      jobnum <- paste0(ind, '_', r)
      opt_type <- 'gp'
      lrnr_type <- 'brt'
      bounds_version <- bv
      experiment_version <- paste0('test', bounds_version, '_', Sys.Date(), ifelse(file_addin == FALSE, '', paste0('_', file_addin)))
      args <- paste0(jobnum, ' ', opt_type, ' ', lrnr_type, ' ', bounds_version, ' ', 
                     experiment_version, ' ', ind, ' ', cv_fold, ' ', bag_fraction, ' ', file_addin)
      system(paste(sys.sub, args))
  
    }
    
    Sys.sleep(1)
  }

}

