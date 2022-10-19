#####################################################################################################################################
### Script objective: Parent script to run MBG models for generating spatio-temporal estimates
#            - tracks progress through stages (polygon resampling and full MBG) and launches child scripts (launch & resampling scripts
#            - specific script is defined in config) when necessary inputs have been generated.
#
### Modeling rationale:
#   The modeling process involves use of point and polygon data. Polygon data are converted to point data through polygon resampling.
#   This resampling is done using a population-based approach (probability of underlying survey location is proportional to population density).
#####################################################################################################################################

#### Setup for qsub commands ########################################################################################################
## Input settings via control files #################################################################################################
user <- Sys.info()[["user"]] ## Get current user name

## Settings for printing qsub command
config_file <- "config_file"
lbd_core_branch <- "develop"

indicator_group <- "oncho"
indicator <- "had_oncho_poly"
parent_script <- "parent_master"
singularity_image <- <<<< FILEPATH REDACTED >>>>

# #Read arguments from qsub command
config_file <- as.character(commandArgs(trailingOnly = T)[1])
indicator_group <- as.character(commandArgs(trailingOnly = T)[2])
indicator <- as.character(commandArgs(trailingOnly = T)[3])
lbd_core_branch <- as.character(commandArgs(trailingOnly = T)[4])
message(as.character(commandArgs(trailingOnly = T)))

#####################################################################################################################################


#### Setup for job launch ###########################################################################################################
## Set location of LBD core repo, based on qsub argument
if (lbd_core_branch == "master") { # version of lbd_core which is updated (manually) from master branch of lbd_core
  core_repo <- <<<< FILEPATH REDACTED >>>>
} else if (lbd_core_branch == "develop") { # version of lbd_core which is updated (manually) from develop branch of lbd_core
  core_repo <- <<<< FILEPATH REDACTED >>>>
} else if (lbd_core_branch == "legacy") { # legacy Focal 3 version of lbd_core
  core_repo <- <<<< FILEPATH REDACTED >>>>
}

message(core_repo)

## Set location of Focal 3 repo (contains Focal 3-specific code)
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))
                  
message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

## Focal 3-specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

## Load config; set_up_config_focal_3() is a wrapper for set_up_config() which accommodates custom covariate options
config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator, 
                        config_file = <<<< FILEPATH REDACTED >>>>, run_tests = FALSE)
message(paste0("covs_name = ", covs_name))

## Focal 3-specific workflow: to restart from or overwrite an existing run_date; add from_run_date argument to config
if (from_run_date != "NULL") {
  run_date <- from_run_date
} else {
  run_date <- make_time_stamp(time_stamp)
}
message(paste0("run_date: ", run_date))


### Module 1: Resample polygons  ####################################################################################################
if (run_poly_resamp) {
  message("Running polygon resampling as specified in your config file. You have specified that '", poly_resamp_script, "' should be used")
  indicator <- indicator_poly
  sharedir <- <<<< FILEPATH REDACTED >>>>
  
  ## Create directory structure for this model run
  create_dirs(indicator_group = indicator_group, indicator = indicator_poly)
  output_dir <- <<<< FILEPATH REDACTED >>>>
  print(output_dir)
  dir.create(paste0(output_dir), showWarnings = TRUE)
  dir.create(<<<< FILEPATH REDACTED >>>>, showWarnings = TRUE)
  dir.create(<<<< FILEPATH REDACTED >>>>, showWarnings = TRUE)
  
  ## Save log of config file
  write.csv(config, <<<< FILEPATH REDACTED >>>>, row.names = FALSE)
  
  ## Launch polygon resampling (parallelized over regions)
  polygon_resampling_job <- parallelize(slots            = fthread_launch,
                                        memory           = m_mem_free_launch,
                                        script           = poly_resamp_script,
                                        geo_nodes        = as.logical(use_geos_nodes),
                                        expand_vars      = list(region = region_list),
                                        save_objs        = c("indicator", "indicator_group", "run_date", "config_file", "core_repo"),
                                        prefix           = "resample",
                                        log_location     = 'sharedir',
                                        script_dir       = <<<< FILEPATH REDACTED >>>>,
                                        run_time         = runtime_launch,
                                        threads          = fthread_launch,
                                        singularity_opts = list(SET_OMP_THREADS = fthread_launch, SET_MKL_THREADS = fthread_launch))
  
  monitor_jobs(polygon_resampling_job)
} else {
  message("Skipping polygon resampling as specified in config")
}

message("Resampling complete")

### Module 2. Running an MBG model on the full dataset (points and resampled polygons)  #############################################
if (run_full == "TRUE") {
  message("Running a full MBG model as specified in your config file. You have specified that '", mbg_model_script, "' should be used")
  indicator <- indicator_pts_poly
  sharedir <- <<<< FILEPATH REDACTED >>>>
  
  ## Create directory structure for this model run
  create_dirs(indicator_group = indicator_group, indicator = indicator_pts_poly)
  output_dir <- <<<< FILEPATH REDACTED >>>>
  print(output_dir)
  dir.create(paste0(output_dir), showWarnings = FALSE)
  dir.create(<<<< FILEPATH REDACTED >>>>, showWarnings = FALSE)
  dir.create(<<<< FILEPATH REDACTED >>>>, showWarnings = FALSE)
  
  ## Save log of config file
  write.csv(config, <<<< FILEPATH REDACTED >>>>, row.names = FALSE)
  
  ## Launch parallel MBG model
  launch_MBG_job <- parallelize(slots            = fthread_launch,
                                memory           = m_mem_free_launch,
                                script           = mbg_model_script,
                                geo_nodes        = as.logical(use_geos_nodes),
                                expand_vars      = list(region = region_list),
                                save_objs        = c("indicator", "indicator_group", "run_date", "config_file", "core_repo"),
                                prefix           = "launch",
                                log_location     = 'sharedir',
                                script_dir       = paste0(indic_repo, "scripts/"),
                                run_time         = runtime_launch,
                                threads          = fthread_launch,
                                singularity_opts = list(SET_OMP_THREADS = fthread_launch, SET_MKL_THREADS = fthread_launch))
  
} else {
  message("Skipping full model run as specified in config")
}
