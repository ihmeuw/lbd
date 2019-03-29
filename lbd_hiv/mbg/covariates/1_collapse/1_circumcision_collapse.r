####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   Circumcision in men, ages 15-49
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- '<<<< FILEPATH REDACTED >>>>/lbd_hiv/'
cores <- 10

# set arguments for this indicator/topic
indicator <- "male_circumcision"
topic <- "male_circumcision"

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))

# find most recent versions
geomatched_version <- most_recent_date(paste0("<<<< FILEPATH REDACTED >>>>", topic, "/"))
prev_version <- most_recent_date("<<<< FILEPATH REDACTED >>>>", file_pattern = indicator)
collapse_version <- format(Sys.time(), "%Y_%m_%d")
## Load and subset data ----------------------------------------------------------------------------
# load geomatched data
data <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>", geomatched_version, "/", topic, "_", geomatched_version, ".RDS"))
data <- data.table(data)

# Keep data for later crosswalk
data_for_crosswalk <- copy(data)

# Drop unneeded variables and remove annoying stata attributes
data <- data[, list(nid, country, survey_series, survey_name, year,
                    strata, psu, point, shapefile, location_code, latitude, longitude,
                    sex_id, age_year, int_year, int_month, male_circumcision, pweight)]
for (ii in 1:ncol(data)) attributes(data[[ ii]]) <- NULL

# rename variables
setnames(data, c("survey_series"), c("source"))

orig_nids <- unique(data$nid)

# fix variable class
data[, year := as.numeric(year)]
data[, point := as.numeric(point)]
data[, latitude := as.numeric(latitude)]
data[, longitude := as.numeric(longitude)]

# Drop NID's related to the country specific Nigeria survey
drop_nigeria_nid <- c(325046, 151719, 324443)
data <- data[!nid %in% drop_nigeria_nid,]

# subset to ages 15-49 (and only surveys with this full range)
data <- data[between(age_year, 15, 49),]
drop <- data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
data <- data[!nid %in% drop,]

# drop observations with missing circumcision response or individual weight
data <- data[!is.na(male_circumcision) & !is.na(pweight),]

# drop point data with missing latitude/longitude
data <- data[(!is.na(latitude) & !is.na(longitude)) | point == 0,]

# drop polygon data with missing shapefile/location code
data <- data[(!is.na(shapefile) & !is.na(location_code)) | point == 1,]

# Drop countries outside Africa
source("<<<< FILEPATH REDACTED >>>>/get_location_metadata.R")
loc <- data.table(get_location_metadata(location_set_id=2, gbd_round_id=4))
loc <- loc[grepl("Africa", region_name) & location_type == "admin0", ihme_loc_id]
data <- data[country %in% loc,]

# store NIDs and sample size to compare at the end
nid_list <- unique(data$nid)
start_n <- data[, list(start = .N), by='country,nid']

# save microdata for indicator after all processing but before collapse
saveRDS(data, file = paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_pre_collapse_", collapse_version, ".rds"))


## Collapse all data (polygon and point) -----------------------------------------------------------
data <- data[, .(int_year = floor(median(int_year, na.rm = T)),
                 male_circumcision = weighted.mean(male_circumcision, pweight),
                 N = sum(pweight) ^ 2 / sum(pweight ^ 2),
                 sum_of_sample_weights = sum(pweight),
                 N_obs = .N),
             by = .(nid, country, source, year, point, shapefile, location_code, latitude, longitude)]

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")


## Add survey report data ---------------------------------------------------------------------------
report <- fread("<<<< FILEPATH REDACTED >>>>/circumcision_extraction.csv")
report[, male_circumcision := male_circumcision / 100]

# Remove survey reports that don't overlap 15-49
report <- report[(start_age <= 15 & end_age >= 15) | (start_age >= 15 & start_age <= 49), ]

# Add survey data that is not in standard format
# identify all age ranges other than 15-49, subset report to covering those age_ranges
age_list <- unique(report[,c("start_age", "end_age")])

#Only add age lists that cover at least 20 years for crosswalk
age_list <-
  age_list %>%
  filter(start_age != 15 | end_age != 49)

crosswalk_report_list <- apply(age_list, 1, function(x) return(report[start_age == x[1] & end_age == x[2]]))

# Prepare microdata for crosswalk
data_for_crosswalk <- crosswalk_data_prepare(data_for_crosswalk, indicator, "pweight")

# Apply crosswalk
crosswalked <-
  rbindlist(lapply(crosswalk_report_list, function(x) {
    crosswalk_reports(data_for_crosswalk, x, indicator, "pweight")
  }))

setnames(crosswalked, c("prev_old", "prev_new"), c("male_circumcision_old", "male_circumcision_new"))
rm(data_for_crosswalk)

# save cross-walk data for future reference
write.csv(crosswalked, file = paste0('<<<< FILEPATH REDACTED >>>>/',
                                     collapse_version, '_crosswalked_data.csv'), row.names = F)

# recombine cross-walked report data with data for 15-49
crosswalked[, c("male_circumcision", "N") := list(male_circumcision_new, N_eff)]
report <- rbind(report[start_age == 15 & end_age == 49,], crosswalked, fill = T)

# drop national estimates
report <- report[location_type != "National",]

# drop countries outside Africa
report <- report[country %in% loc,]

# drop any NIDs we have microdata for
report <- report[!nid %in% data$nid,]

# approximate the effective sample size using the design effects from the polygon microdata
# Note: this is *very* rough; the purpose is to make sure the report data doesn't have an unfair
# advantage compared to microdata; the calculation below excludes surveys where there is no apparent
# design effect -- this is not super plausible, and we don't want these to skew the distribution.
de <- data[point == 0, list(de = sum(N)/sum(N_obs)), by='nid,shapefile,location_code'][de < 0.99999, median(de)]
report[, N_obs := N]
report[, N := N_obs * de]

# update the survey list
nid_list <- c(nid_list, unique(report$nid))

# combine with collapsed microdata
report <- report[, list(nid, country, survey_series, int_year, point, shapefile, location_code, latitude, longitude, male_circumcision, N, N_obs)]
setnames(report, c("int_year", "survey_series"), c("year", "source"))

# Mark report data, and make sum of sample weights equal to N to approximate later in visualization plots
report[, sum_of_sample_weights := N]
report[, report := 1]
data <- rbind(data, report, fill=T)

# Now that all data is compiled, assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]

## Convert to counts and resample the polygon data -------------------------------------------------
data[, male_circumcision := male_circumcision * N]
data <- resample_polygons(data = data, indic = "male_circumcision", cores = cores, pull_poly_method = "fast", seed = 98121)

# use resampling weights to down-weight N and male_circumcision, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, male_circumcision := male_circumcision * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, male_circumcision, weight, sum_of_sample_weights, report, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# report differences compared to the most recent version
compare_most_recent(data, indicator, prev_version)

## save input data for model runs & log dates-------------------------------------------
save_data(data, indicator, geomatched_version, collapse_version)

