###############################################################################
## MBG launch, aggregate results, and diagnostics launcher script for diarrhea
###############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'ort'
parallel_script <- 'had_diarrhea/parallel_dia'
launch_script   <- 'launch_dia.R' # MBG launch script: 'launch_dia.R'; Raking launch script: 'launch_dia_rake.R'

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

# set config and covariate files
config_par   <- 'dia_best'
covar_par      <- 'region_specific'

# set cov version
cov_version <- 'standard_inc_mbg'

# set whether making holdouts (only matters if launching raking script)
makeholdouts <- TRUE

# list all regions or countries
regions <- c('dia_south_asia-ind', 'dia_malay', 'dia_wssa', 
             'dia_s_america', 'dia_mcaca', 'dia_name',
             'dia_mid_east', 'dia_essa', 'dia_afr_horn', 'dia_cssa',
             'dia_se_asia', 'dia_sssa', 
             'dia_central_asia', 'IND', 'MNG')

# list indicators
indics <- 'had_diarrhea'

# set cores
mymem <- ifelse(launch_script == 'launch_dia.R', '20G', '10G')


## Run launch scripts -------------------------------------------------------------------------
  
# set run date
if (use_old_run_date == FALSE) {
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- old_run_date_input
}

# make sure that only selecting a previous run_date intentionally
if (use_old_run_date == TRUE) {
  prev <- readline('Are you sure you want to use a previous run date? Y or N: ')
  if (prev != 'Y') stop('Set use_old_run_date to FALSE.')
}

for (r in regions) {
  
  # set specific arguments
  Regions         <- r
  indicator       <- indics
  jname           <- paste(indicator, Regions, sep = '_')
  
  # set region specific covariates, if desired
  if (covar_par == 'region_specific') {
    cov_par <- paste0('covs_', cov_version, '/covs_ort_', r)
  } else {
    cov_par <- covar_par
  }
  message(cov_par)
  
  # some quick checks for the arguments
  if(use_old_run_date == TRUE & old_run_date_input == '') stop('You indicated using an old run date; please provide an old run date')
  
  # set up qsub
  sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
                    '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                    '-l fthread=1 -l h_rt=00:02:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE')
  r_shell <- paste0(repo, 'mbg_central/share_scripts/shell_sing.sh')
  script <- paste0('<<<< FILEPATH REDACTED >>>>/had_diarrhea/', launch_script)
  args <- paste(user, repo, indicator_group, indicator, config_par, cov_par, Regions, parallel_script,
                plot_covariates, covariate_plotting_only, proj_arg, use_geos_nodes, run_date,
                makeholdouts, cov_version)
  
  # run launch script
  system(paste(sys.sub, r_shell, script, args))
  
}

