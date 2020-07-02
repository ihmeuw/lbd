####################################################################################################
## Description: Checks Circumcision pre-collapsed/geomatched datasets
## Indicator:    Male Circumcision
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this indicator/topic and dataset dates
indicator <- "male_circumcision"
topic <- "male_circumcision"

# set arguments for this run and user.
core_repo  <- "/lbd_core/"
indic_repo <- "/lbd_hiv/"

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))

source(paste0(indic_repo, '/mbg/covariates/functions/subset_geomatch_collapse_functions.r'))
source(paste0(indic_repo, 'data/covariate_surveys/functions/geomatch_functions.r'))

# Set collapse and precollapse dates
######
geomatched_old <- "2020_02_07"
geomatched_new <- "2020_03_05"

old <- readRDS("<<<< FILEPATH REDACTED >>>>")
old <- subset_geomatched_male_circumcision(old)
old[, location_code := as.numeric(location_code)]

new <- readRDS("<<<< FILEPATH REDACTED >>>>")
new <- subset_geomatched_male_circumcision(new)

X <- new[new$nid == '411301',]

# Change to generic indicator column
setnames(old, indicator, "indicator")
setnames(new, indicator, "indicator")

# Create directory for saving
save_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
dir.create(save_dir, showWarnings = F)

sink(paste0(save_dir, "/data_changes.txt"))
## tells us what new nids we have in the pre-collapsed dataset
cat(paste0("\n\n Checking differences between geomatched data from current date (", geomatched_new, ") to data from ", geomatched_old)) 
cat(paste0("\n\n*****NIDs that have been added:\n", paste(sort(setdiff(new$nid, old$nid)), collapse=","), "\n\n"))
cat(paste0("\n\n*****NIDs that have been dropped:\n", paste(sort(setdiff(old$nid, new$nid)), collapse=","), "\n\n"))

# See data changes
geomatched_data_changes(old, new)

# Make time trend plots
geomatch_time_plots(new, old, save_plot = save_dir)
sink()
