##############################################################################
## MBG launch, aggregate results, and diagnostics launcher script 
## Created 2018/02/23
##############################################################################
## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

#REDACTED

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
package_list <- c(t(read.csv(paste0(core_repo, '/#REDACTED/package_list.csv'), header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# set cluster arguments
use_geos_nodes  <- T
proj_arg        <- ifelse(use_geos_nodes, 'proj_geo_nodes', 'proj_geospatial_dia')
proj            <- ifelse(use_geos_nodes, paste0(' -P ', proj_arg, ' -l gn=TRUE '), paste0(' -P ', proj_arg, ' '))

# set covariate arguments
plot_covariates <- TRUE
covariate_plotting_only <- FALSE

# indicate whether to use old run date
use_old_run_date <- F
old_run_date_input <- '2020_09_01_11_42_52'

# set run date
if (use_old_run_date == FALSE) {
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- old_run_date_input
}

# set config and covariate files
config_par   <- 'hap_sp_fine'
covar_par      <- 'region_specific'

# set whether running for individual countries
individual_countries <- FALSE


# indicate holdout (also need to do so in config)
holdout <- TRUE # only matters if running aggregation, set to TRUE for holdouts

# list all regions or countries
# standard regions
regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_sssa', 
             'dia_mcaca', 'dia_s_america-GUF', 'dia_central_asia', 'dia_chn_mng', 
             'dia_se_asia', 'dia_malay', 'dia_south_asia', 'dia_mid_east', 'dia_essa')

# custom country-specifics
regions <- c('essa-ERI-DJI-YEM', "ERI+DJI+YEM",
             'sssa-ZAF', 'ZAF',
             'cssa-AGO-GNQ', 'AGO',
             'wssa-CPV-NGA', 'NGA',
             'noaf-ESH',
             'caca-CUB',
             'ansa-VEN', 'trsa-GUF',
             'stan-TKM',
             'CHN', 'MNG',
             'ocea-MYS',
             'seas',
             'mide+TKM', 'soas')

# list indicators
indics <- 'cooking_fuel_solid'

## Run launch scripts -------------------------------------------------------------------------

for (i in indics) {
   
  # make sure that only selecting a previous run_date intentionally
  if (use_old_run_date == TRUE) {
    prev <- readline('Are you sure you want to use a previous run date? Y or N: ')
    if (prev != 'Y') stop('Set use_old_run_date to FALSE.')
  }
  
  for (reg in regions) {
    
    # set specific arguments
    indicator       <- i
    jname           <- paste('rocket', indicator_group, reg, sep = '_')
    mymem           <- '20G'
    
    # set region specific covariates, if desired
    if (covar_par == 'region_specific') cov_par <- paste0('cooking_', reg)
    else cov_par <- covar_par
    
    # some quick checks for the arguments
    if(use_old_run_date == TRUE & old_run_date_input == '') stop('You indicated using an old run date; please provide an old run date')
    
    # set up qsub
    sys.sub <- paste0('qsub -e /#REDACTED/', user,'/errors -o /#REDACTED/', user, '/output ', 
                      '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                      '-l fthread=1 -l h_rt=', ifelse(use_geos_nodes, '16:00:00:00', '3:00:00:00'),
                      ' -v sing_image=default -N ', jname, ' -l archive=TRUE ')
    r_shell <- file.path(core_repo, 'mbg_central/share_scripts/shell_sing.sh')
    script <- '<<< FILEPATH REDACTED >>>'
    args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, parallel_script,
                  plot_covariates, covariate_plotting_only, proj_arg, use_geos_nodes, run_date, my_repo)

    # run launch script
    paste(sys.sub, r_shell, script, args) %>% 
      system
    
  }
  
}