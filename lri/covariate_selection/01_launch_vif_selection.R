##############################################################################
## Variance inflation factor (VIF) launcher script
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
config_par   <- 'config_62'
cov_par      <- 'covs_27'
use_old_run_date <- FALSE
old_run_date_input <- ''
individual_countries <- FALSE

user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'lri'
indicator <- 'has_lri'
config_file <- '/lri/'
cov_file <- config_file

# indicate threshold parameters
threshold_min <- 2
threshold_max <- 5
threshold_step <- 1

# indicate whether or not to crop covariates (must give .RData file with cropped covariate file if crop_covs = FALSE)
crop_covs <- TRUE
cropped_covs_file <- ''

# list indicators
indics <- c('has_lri')

# specify 'all' for parent regions or 'countries' for country-specific regions
regions <- 'all'

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>',header=FALSE)))
source(paste0(repo, '/lbd_core/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos ='<<<< FILEPATH REDACTED >>>>')

# set run date
if (use_old_run_date == FALSE) {
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- old_run_date_input
}


# Make qsub string and submit -----------------------------------------
mem <- '20G'
rt <- '8:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
name <- 'has_lri_vif_selection'
shell <- '<<<< FILEPATH REDACTED >>>>'
code <- '<<<< FILEPATH REDACTED >>>>/lri/covariate_selection/02_vif_selection_script.R'

args <- paste(user, 
              repo, 
              indicator_group, 
              indicator, 
              config_par, 
              config_file, 
              cov_par, 
              cov_file, 
              run_date,
              regions, 
              threshold_min, 
              threshold_max, 
              threshold_step, 
              individual_countries,
              crop_covs, 
              cropped_covs_file)

qsub <- paste0('qsub -l m_mem_free=', mem,
               ' -l fthread=1 -l h_rt=', rt, 
               ' -v sing_image=default -q ', queue, 
               ' -P ', proj, 
               ' -N ', name, 
               ' ', shell, 
               ' ', code, 
               ' ', args)

system(qsub)
