##############################################################################
## MBG aggregate results launcher script for ORT
##############################################################################


## Setup -------------------------------------------------------------------------

## clear environment
rm(list=ls())

## Set repo location, indicator group, and some arguments
user            <- commandArgs()[4]
core_repo       <- commandArgs()[5]
indicator_group <- commandArgs()[6]
indicator       <- commandArgs()[7]
config_par      <- commandArgs()[8]
cov_par         <- commandArgs()[9]
Regions         <- commandArgs()[10]
holdout         <- as.numeric(commandArgs()[15])

## Load MBG packages
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')

## Read config file and save all parameters in memory
config <- set_up_config(repo            = core_repo,
                        indicator_group = '',
                        indicator       = '',
                        config_name     = paste0('ors/3_modeling/config_', config_par),
                        covs_name       = paste0('ors/3_modeling/', cov_par))

## Set to prod or geos nodes
proj <- commandArgs()[11]
use_geos_nodes <- commandArgs()[12]

## Set run date
run_date <- commandArgs()[13]

## Create output folder with the run_date
outputdir      <- '<<<< FILEPATH REDACTED >>>>'

## Create proper year list object
if (class(year_list) == 'character') year_list <- eval(parse(text=year_list))

## Get measure to aggregate
measure <- commandArgs()[14]

# set cores by region
individual_countries <- ifelse(nchar(Regions) == 3, TRUE, FALSE)
r <- Regions
region_cores <- 6
if(r == 'dia_malay' | r == 'dia_name') region_cores <- 8
if(r == 'dia_chn_mng' | r == 'dia_wssa' | r =='dia_south_asia') region_cores <- 10
if(r == 'dia_s_america') region_cores <- 12
if(makeholdouts) region_cores <- round(region_cores*0.8)

# convert cores approximately to memory and run time
region_rt <- '01:00:00:00'
region_mem <- region_cores*15


## Aggregate to admin2, admin1 and admin0 -------------------------------------------------------------------------

message('Submitting aggregation script')


submit_aggregation_script(indicator       = indicator, 
                          indicator_group = indicator_group,
                          run_date        = run_date,
                          raked           = FALSE,
                          pop_measure     = pop_measure,
                          overwrite       = T,
                          ages            = 0,
                          holdouts        = holdout,
                          regions         = Regions,
                          corerepo        = core_repo,
                          log_dir         = outputdir,
                          cores           = 4,
                          memory          = region_mem,
                          run_time        = region_rt,
                          proj            = proj,
                          geo_nodes       = as.logical(use_geos_nodes),
                          singularity     = 'default',
                          measure         = measure,
                          modeling_shapefile_version = modeling_shapefile_version,
                          raking_shapefile_version = raking_shapefile_version)


