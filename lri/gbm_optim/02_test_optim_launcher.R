##############################################################################
## Launch script to run BRT optimization
##############################################################################


## Setup -------------------------------------------------------------------------

# set repo
setwd('<<<< FILEPATH REDACTED >>>>')

# set node preference
nodes <- 'geos'
if (nodes == 'geos') {
  proj <- '-P proj_geo_nodes_dia -l gn=TRUE'
  r_shell <- 'shell_geos.sh'
} else {
  proj <- '-P proj_geospatial'
  r_shell <- 'shell_prod.sh'
}

# Use hf_inla
hf_inla_shell <- '<<<< FILEPATH REDACTED >>>>'

r_shell <- hf_inla_shell
mkl <- 1

# region(s) to be launched
regions <- c('dia_se_asia')

# indicator(s) to be launched
indis <- c('had_diarrhea')

# cv_folds and bag_fractions to run pairwise tests on
cv_folds <- c(2, 3, 4, 5)
bag_fractions <- c(1-(1/2), 1-(1/3), 1-(1/4), 1-(1/5))


## Launch BRT optimizer scripts -------------------------------------------------------------------------

i <- 1

for (fold in cv_folds) {
  
  for (frac in bag_fractions) {
    
    # get test arguments
    cv_fold <- fold
    bag_fraction <- frac
    experiment <- i
    i <- i + 1
    
    # get other arguments
    r <- regions[1]
    ind <- indis[1]
    jname <- paste0(ind, '_', r, '_exp', experiment)
    mycores <- 10
    sys.sub <- paste0('qsub ',proj,paste0(' -e <<<< FILEPATH REDACTED >>>>',' -o <<<< FILEPATH REDACTED >>>>'),
                      '-cwd -N ', jname, ' ', '-pe multi_slot ', mycores)
    
    # launch script name to qsub
    script <- 'test_runGBM.R'
    jobnum <- paste0(ind, '_', r)
    opt_type <- 'gp'
    lrnr_type <- 'brt'
    bounds_version <- 18
    experiment_version <- paste0('test', bounds_version, '_', Sys.Date())
    args <- paste0(jobnum, ' ', opt_type, ' ', lrnr_type, ' ', bounds_version, ' ', experiment_version, ' ', ind, ' ', cv_fold, ' ', bag_fraction, ' ', experiment)
    system(paste(sys.sub, r_shell, mkl, script, args))
    
  }
  
}
  
