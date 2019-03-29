####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   sti symptoms (genital discharge & sore/ulcer) in men & women, ages 15-49
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- '<<<< FILEPATH REDACTED >>>>/lbd_hiv/'
cores <- 10

# set arguments for this indicator/topic
indicator <- "sti_symptoms"
topic <- "sti_symptoms"

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

# drop unneeded variables and remove annoying stata attributes
data <- data[, list(nid, country, survey_series, survey_name, year_n,
                    strata, psu, point, shapefile, location_code, latitude, longitude,
                    sex_id, age_year, int_year, had_intercourse, discharge_or_sore, pweight)]
for (ii in 1:ncol(data)) attributes(data[[ ii]]) <- NULL

# rename variables
setnames(data, c("survey_series", "year_n"), c("source", "year"))

# fix variable class
data[, point := as.numeric(point)]
data[, latitude := as.numeric(latitude)]
data[, longitude := as.numeric(longitude)]

# only want surveys 1998 or later for now
data <- data[year >= 1998,]

# only want surveys that asked both men and women
svys = data[,list(sex_split = mean(sex_id, na.rm=T)),by='nid,country']
data <- data[nid %in% svys[sex_split > 1 & sex_split < 2,nid]]

# only want respondents who have had intercourse
data <- data[had_intercourse == 1]

# subset to ages 15-49 (and only surveys with this full range)
data <- data[between(age_year, 15, 49),]
drop <- data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
data <- data[!nid %in% drop,]

# add a new variable combining the sti_symptoms: discharge_OR_sore
data[,sti_symptoms := discharge_or_sore]

# drop observations with missing sti_symptoms or pweight
data <- data[!is.na(pweight) & !is.na(sti_symptoms),]

# drop point data with missing latitude/longitude
data <- data[(!is.na(latitude) & !is.na(longitude)) | point == 0,]

# drop polygon data with missing shapefile/location code
data <- data[(!is.na(shapefile) & !is.na(location_code)) | point == 1,]

# drop countries outside Africa
loc <- source("<<<< FILEPATH REDACTED >>>>/get_location_metadata.R")
loc <- data.table(get_location_metadata(location_set_id=2, gbd_round_id=4))
loc <- loc[grepl("Africa", region_name) & location_type == "admin0", ihme_loc_id]
data <- data[country %in% loc,]

# store NIDs and sample size to compare at the end
nid_list <- unique(data$nid)
start_n <- data[, list(start = .N), by='country,nid']

# save microdata for indicator after all processing but before collapse
saveRDS(data, file = paste0("<<<< FILEPATH REDACTED >>>>", topic, "/",
                            "pre_collapse/", indicator, "_pre_collapse_", collapse_version, ".rds"))

## Collapse all data (polygon and point) -----------------------------------------------------------
data <- data[, list(int_year = floor(median(int_year, na.rm=T)),
                    #int_month = mean(int_month, na.rm=T),
                    sti_symptoms = weighted.mean(sti_symptoms, pweight),
                    N = sum(pweight)^2/sum(pweight^2),
                    sum_of_sample_weights = sum(pweight),
                    N_obs = .N),
             by='nid,country,source,year,point,shapefile,location_code,latitude,longitude']

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")

# assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]

## Convert to counts and resample the polygon data -------------------------------------------------
data[, sti_symptoms := sti_symptoms * N]
data <- resample_polygons(data = data, indic = "sti_symptoms", cores = cores, pull_poly_method = "fast", seed = 98121)

# use resampling weights to down-weight N and sti_symptoms, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, sti_symptoms := sti_symptoms * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, sti_symptoms, weight, sum_of_sample_weights, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# report differences compared to the most recent version
compare_most_recent(data, indicator, prev_version)

## save input data for model runs & log dates-------------------------------------------
save_data(data, indicator, geomatched_version, collapse_version)

