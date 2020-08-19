##########################################################################################
###########  Launch bootstraps of frequentist crosswalk models, in parallel  #############
##########################################################################################

##################################################################################
###### I. Setup
##################################################################################

diagnostic_test <- commandArgs(trailingOnly= T)[1]
if (diagnostic_test == "ict_fts") {
  diagnostics <- c("ict", "fts")
} else if (diagnostic_test == "ict_fts_mf") {
  diagnostics <- c("ict", "fts", "mf")
} else {
  diagnostics <- diagnostic_test
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

user <- Sys.info()[["user"]] ## Get current user name
core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")
path <- paste0(<<<< FILEPATH REDACTED >>>>)

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

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
indicator_group <- "lf"
indicator <- "had_lf_w_resamp"
run_date <- NULL

##################################################################################
###### II. Load and prep data
##################################################################################

## Load crosswalk training data set from file
data_final <- fread(file=<<<< FILEPATH REDACTED >>>>)

# Post hoc edits
data_final <- as.data.table(data_final)
data_final <- data_final[!(nid == 403329 & age_start == 2 & age_end == 94)]
data_final <- data_final[!(nid == 409067 & age_start == 5 & age_end == 18)]
data_final[nid == 143009, nid := 387657]
data_final <- rbind(data_final, data_final[nid == 387657 & age_start == 61 & age_end == 70])
data_final[nrow(data_final), c("age_start", "age_end", "cases", "sample_size") := list(71, 94, 6, 52)]
data_final <- data_final[!(nid == "403316" & ((site_memo == "Lewomada, Sikka District" & year > 2012) | (site_memo == "Paga, Sikka District" & year > 2011) | (site_memo == "Pruda, Sikka District" & year > 2011)))]

data_final <- data_final[data_final$diagnostic %in% diagnostics]

### If both ICT and FTS are reported for a given study, prefer FTS
ag_cohorts <- unique(data_final[diagnostic %in% c("fts", "ict"), cohort])
dx_by_group <- as.data.table(aggregate(diagnostic ~ cohort, data = data_final[cohort %in% ag_cohorts,], unique))
fts_ict <- data_final[cohort %in% dx_by_group[diagnostic %in% c("c(\"ict\", \"mf\")"), cohort]]
data_final <- data_final[!((cohort %in% fts_ict$cohort) & diagnostic == "ict")]

### Subset to studies reporting > 1 age group for diagnostic
age_groups_by_nid <- as.data.table(aggregate(uid_crosswalk ~ cohort, data = data_final, length))
data_final <- data_final[cohort %in% age_groups_by_nid[uid_crosswalk > 1, cohort]]

##################################################################################
###### III. Generate bootstrap samples
##################################################################################
bootstrap_samples <- list()

if (script %in% c("dx_crosswalk_cluster_freq_parallel")) { # dx crosswalk
  ## Retrieve cohorts with multiple dx types (skipping br and fts)
  dx_by_group <- as.data.table(aggregate(diagnostic ~ cohort, data = data_final[order(diagnostic), ], unique))
  multiple_dx <- data_final[cohort %in% dx_by_group[diagnostic %in% c("c(\"ict\", \"mf\")"), cohort]]
  
  ## Convert to wide format
  multiple_dx_wide <- reshape(multiple_dx, idvar = c("cohort", "age_start", "age_end"), timevar = "diagnostic", direction = "wide")
  multiple_dx_wide$geom_mean_sample_size <- round(sqrt(multiple_dx_wide$sample_size.ict * multiple_dx_wide$sample_size.mf))
  multiple_dx_wide$prev.ict <- multiple_dx_wide$cases.ict / multiple_dx_wide$sample_size.ict
  multiple_dx_wide$prev.mf <- multiple_dx_wide$cases.mf / multiple_dx_wide$sample_size.mf
  multiple_dx_wide$geom_mean_sample_size.cases.ict <- round(multiple_dx_wide$prev.ict * multiple_dx_wide$geom_mean_sample_size)
  multiple_dx_wide$geom_mean_sample_size.cases.mf <- round(multiple_dx_wide$prev.mf * multiple_dx_wide$geom_mean_sample_size)
  multiple_dx_wide <- multiple_dx_wide[!is.na(cases.ict) & !is.na(cases.mf)]
  data_final <- copy(multiple_dx_wide)
  
  data_final <- data_final[(cases.mf > 0) & (cases.ict > 0)]
  
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
} else { # age crosswalk
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
}

##################################################################################
###### IV. Launch parallel bootstrap jobs
##################################################################################
bootstrap_job <-          parallelize(slots            = 1,
                                      memory           = 1,
                                      script           = script,
                                      geo_nodes        = TRUE,
                                      expand_vars      = list(replicate = 1:replicates),
                                      save_objs        = c("diagnostic_test", "initial_iter", "burnin_iter", "final_iter", "save_name", "load_and_continue", "use_freq", "replicates", "save_folder", "bootstrap_samples", "core_repo"),
                                      prefix           = paste0("bootstrap_", diagnostic_test),
                                      log_location     = 'sharedir',
                                      script_dir       = paste0(indic_repo, "scripts/"),
                                      run_time         = "144:00:00",
                                      threads          = 1,
                                      priority         = -1)


