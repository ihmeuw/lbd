##############################################################################
## MBG launch, aggregate results, and diagnostics launcher script for ORT
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'ort'
parallel_script <- 'ors/3_modeling/parallel_ort'
launch_script   <- 'launch_ort.R' # MBG launch script: 'launch_ort.R'; Raking launch script: 'launch_ort_aggregate.R'

# Load MBG packages and functions
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = repo)

# set cluster arguments
use_geos_nodes  <- TRUE
proj_arg        <- 'proj_geo_nodes_dia'

# set covariate arguments
plot_covariates <- TRUE
covariate_plotting_only <- FALSE

# indicate whether to use old run date
use_old_run_date <- FALSE
old_run_date_input <- ''

# indicate holdout (also need to do so in config)
holdout <- TRUE # only matters if running aggregation, set to TRUE for holdouts

# set config and covariate files
config_pars <- c('ort_best_ho_ad_1000draws')
covar_par <- 'region_specific'

# list all regions or countries
regions <- c('SEN', 'SLE', 'MLI')

# list indicators
indics <- c('any_ors', 'no_ort', 'rhf_only')


## Run launch scripts -------------------------------------------------------------------------

for (config_par in config_pars) {
  
  # set run date
  if (use_old_run_date == FALSE) {
    run_date <- make_time_stamp(TRUE)
  } else {
    run_date <- old_run_date_input
  }

	for (i in indics) {
	  
	  # make sure that only selecting a previous run_date intentionally
	  if (use_old_run_date == TRUE) {
	    prev <- readline('Are you sure you want to use a previous run date? Y or N: ')
	    if (prev != 'Y') stop('Set use_old_run_date to FALSE.')
	  }
	  
	  for (r in regions) {
	    
	    # set specific arguments
	    Regions         <- r
	    indicator       <- i
	    jname           <- paste(indicator, Regions, sep = '_')
	    mymem           <- '20G'
	    
	    # set region specific covariates, if desired
	    if (covar_par == 'region_specific') {
	      cov_par <- paste0('covs_', indicator, '/covs_ort_', Regions)
	    } else {
	      cov_par <- covar_par
	    }
	    
	    # some quick checks for the arguments
	    if(use_old_run_date == TRUE & old_run_date_input == '') stop('You indicated using an old run date; please provide an old run date')
	    
	    # set up qsub
	    sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
	                      '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
	                      '-l fthread=1 -l h_rt=00:02:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE ')
	    r_shell <- paste0(repo, 'mbg_central/share_scripts/shell_sing.sh')
	    script <- paste0('<<<< FILEPATH REDACTED >>>>', launch_script)
	    args <- paste(user, repo, indicator_group, indicator, config_par, cov_par, Regions, parallel_script,
	                  plot_covariates, covariate_plotting_only, proj_arg, use_geos_nodes, run_date, holdout)
	    
	    # run launch script
	    system(paste(sys.sub, r_shell, script, args))
	    
	  }
	  
	}

}