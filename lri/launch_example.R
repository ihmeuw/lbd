###############################################################################
## MBG launch script
##
## Example inital launch script for LRI
##
###############################################################################
#comment

## Setup -------------------------------------------------------------------------

## clear environment
rm(list=ls())

## Set repo location, indicator group, and some arguments
user            <- Sys.info()['user']
indicator_group <- 'lri'
indicator       <- 'has_lri'
core_repo       <- '<<<< FILEPATH REDACTED >>>>' # use the core repo so you don't have to worry about keeping it up to date
indicator_repo  <- '<<<< FILEPATH REDACTED >>>>' # source the configs from your indicator repo
config_par      <- 'example' # you can change these for your custom configs (for example, you may want change to the GADM shapefiles using modeling_shapefile_version 2018_08_01)
cov_par         <- 'example' # you can change these for your custom covariates (for example, you may want to add other MBG covariates and add GBD covariates)
Regions         <- c('essa') # can be all regions for a global model, or can be an individual country, e.g. 'PER'
message(indicator)

## drive locations
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#load my custom functions
souce(paste)

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')

## Read config file and save all parameters in memory
config <- load_config(repo            = indicator_repo,
                      indicator_group = '',
                      indicator       = '',
                      config_name     = paste0('config_', config_par),
                      covs_name       = paste0('covs_', cov_par))

## Ensure you have defined all necessary settings in your config
check_config()

## Set project
proj <- ifelse(as.logical(use_geos_nodes), 'proj_geo_nodes', 'proj_geospatial')

## Create ru'n date in correct format
run_date <- '2019_02_20_11_44_56'

## Create output folder with the run_date
outputdir      <- '<<<< FILEPATH REDACTED >>>>'
dir.create(outputdir)

## Make sure year object is in the correct format
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

## If running individual countries make sure all country FEs and REs off
if (individual_countries) {
  use_child_country_fes <- FALSE
  use_inla_country_fes  <- FALSE
  use_inla_country_res  <- FALSE
}

## Make holdouts -------------------------------------------------------------------------

if(makeholdouts){
  message('Making holdouts')
  
  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        simple      = NULL,
                        removeyemen = TRUE,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share))
  
  # add in location information
  df <- merge_with_ihme_loc(df)
  
  # make a list of dfs for each region, with 5 qt folds identified in each
  stratum_ho <- make_folds(data       = df,
                           n_folds    = as.numeric(n_ho_folds),
                           spat_strat = 'qt',
                           temp_strat = 'prop',
                           strat_cols = 'region',
                           ts         = as.numeric(ho_ts),
                           mb         = as.numeric(ho_mb))
}


## Launch parallel script -------------------------------------------------------------------------

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(makeholdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)

## loop over them, save images and submit qsubs
for(i in 1:nrow(loopvars)){
  
  message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))
  
  # make a qsub string
  qsub <- make_qsub_share(age           = loopvars[i,2],
                          reg           = as.character(loopvars[i,1]),
                          holdout       = loopvars[i,3],
                          test          = F,
                          indic         = indicator,
                          saveimage     = TRUE,
                          memory        = ifelse(individual_countries, 10, 90), # this is conservative, will need to profile your jobs by region
                          cores         = ifelse(as.logical(use_geos_nodes), 10, as.numeric(slots)),
                          proj          = proj,
                          geo_nodes     = as.logical(use_geos_nodes),
                          corerepo      = core_repo,
                          code          = NULL, # eventually you may want to specifiy your own custom parallel script (for example, if you want to output means maps or directly go to post-estimate and aggregation)
                          addl_job_name = paste0(indicator, '_parallel'),
                          singularity   = 'default')
  # submit job
  system(qsub)
  
}