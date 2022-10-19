##########################################################################################
###########  Launch bootstraps of frequentist crosswalk models, in parallel  #############
##########################################################################################

##################################################################################
###### I. Setup
##################################################################################

diagnostic_test <- commandArgs(trailingOnly= T)[1]
if (diagnostic_test == "nod_ss") {
  diagnostics <- c("nod", "ss")
}
message(diagnostics)

initial_iter <- commandArgs(trailingOnly= T)[2]
message(initial_iter)

burnin_iter <- commandArgs(trailingOnly= T)[3]
message(burnin_iter)

final_iter <- as.integer(commandArgs(trailingOnly= T)[4]); message(final_iter)
save_name <- commandArgs(trailingOnly= T)[5]; message(save_name)
load_and_continue <- commandArgs(trailingOnly= T)[6]; message(load_and_continue)
use_freq <- commandArgs(trailingOnly= T)[7]; message(use_freq)
replicates <- commandArgs(trailingOnly= T)[8]; message(replicates)
save_folder <- commandArgs(trailingOnly= T)[9]; message(save_folder)

script <- commandArgs(trailingOnly= T)[10]; message(script)

perform_bootstrap <- commandArgs(trailingOnly= T)[11]; message(script)

user <- Sys.info()[["user"]] ## Get current user name
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")
path <- <<<< FILEPATH REDACTED >>>>

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

library(jsonlite)
mbg_setup(package_list = package_list, repos = core_repo)

library(ggthemes)
library(fasterize)
library(tidyverse)
library(fda)
library(subplex)
library(faraway)
library(BayesianTools)
library(spdep)
library(boot)

singularity_version <- <<<< FILEPATH REDACTED >>>>
indicator_group <- "oncho"
indicator <- "had_oncho_w_resamp"

run_date <- make_time_stamp(TRUE)
message(paste0("run_date: ", run_date))
## Create directory structure for this model run
create_dirs(indicator_group = indicator_group, indicator = indicator)
output_dir <- <<<< FILEPATH REDACTED >>>>
print(output_dir)
dir.create(paste0(output_dir), showWarnings = TRUE)
dir.create(paste0(output_dir, "/errors"), showWarnings = TRUE)
dir.create(paste0(output_dir, "/output"), showWarnings = TRUE)

##################################################################################
###### II. Load and prep data
##################################################################################

## Load crosswalk training data set from file
data_final <- fread(<<<< FILEPATH REDACTED >>>>)

# Post hoc edits
data_final <- data_final[data_final$diagnostic %in% diagnostics]

## Subset to studies reporting > 1 age group for either diagnostic (ss or nod)

if (perform_bootstrap) {
  ##################################################################################
  ###### III. Generate bootstrap samples
  ##################################################################################
  bootstrap_samples <- list()
  
  ## Bootstrap data_final by cohort (not by individual data rows)
  cohorts <- unique(data_final$cohort)
  n_cohorts <- length(cohorts)
  cohorts_sampled <- as.data.table(t(replicate(replicates, sample(cohorts, replace = TRUE))))
  
  for (i in 1:replicates) {
    temp_dt <- data.table()
    for (j in 1:n_cohorts) {
      temp_dt <- rbind(temp_dt, data_final[cohort == as.data.frame(cohorts_sampled)[i, j]])
    }
    bootstrap_samples[[i]] <- copy(temp_dt)
  }
  
  ##################################################################################
  ###### IV. Launch parallel bootstrap jobs
  ##################################################################################
  bootstrap_job <-          parallelize(slots            = 1,
                                        memory           = 1,
                                        script           = script,
                                        geo_nodes        = TRUE,
                                        expand_vars      = list(replicate = 1:replicates),
                                        save_objs        = c("diagnostic_test", "initial_iter", "burnin_iter", "final_iter", "save_name", "load_and_continue", "use_freq", "replicates", "save_folder", "bootstrap_samples", "core_repo", "perform_bootstrap"),
                                        prefix           = paste0("bootstrap_", diagnostic_test),
                                        log_location     = 'sharedir',
                                        script_dir       = paste0(indic_repo, "scripts/"),
                                        run_time         = "360:00:00",
                                        threads          = 1,
                                        priority         = -1)
  
} else { # run analysis with full data set
  replicates <- 0
  bootstrap_job <-          parallelize(slots            = 1,
                                        memory           = 1,
                                        script           = script,
                                        geo_nodes        = TRUE,
                                        expand_vars      = list(replicate = 0),
                                        save_objs        = c("diagnostic_test", "initial_iter", "burnin_iter", "final_iter", "save_name", "load_and_continue", "use_freq", "replicates", "save_folder", "data_final", "core_repo", "perform_bootstrap"),
                                        prefix           = paste0("oncho_crosswalk_", diagnostic_test),
                                        log_location     = 'sharedir',
                                        script_dir       = paste0(indic_repo, "scripts/"),
                                        run_time         = "360:00:00",
                                        threads          = 1,
                                        priority         = -1)
}

