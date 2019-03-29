####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   Condom use at the last sexual encounter in the last 12 months
##              among men & women who had sex in the last 12 months, ages 15-49
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- '<<<< FILEPATH REDACTED >>>>/lbd_hiv/'
cores <- 10

# set arguments for this indicator/topic
indicator <- "condom_last_time_BOTH"
topic <- "rsp"

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
data <- data[, list(nid, country, survey_series, survey_name, year,
                    strata, psu, point, shapefile, location_code, latitude, longitude,
                    sex_id, age_year, int_year, pweight, sex_in_last_12_months, condom_last_time)]
for (ii in 1:ncol(data)) attributes(data[[ ii]]) <- NULL

# rename variables
setnames(data, c("survey_series"), c("source"))

# fix variable class
data[, year := as.numeric(year)]
data[, point := as.numeric(point)]
data[, latitude := as.numeric(latitude)]
data[, longitude := as.numeric(longitude)]

# only want surveys 1998 or later for now
data <- data[year >= 1998,]

orig_nids <- unique(data$nid)

# subset to ages 15-49 (and only surveys with this full range)
# doing pre sex in last year restriction because one survey has
# all 15-17 year-old-years who have not had sex in last year
data <- data[between(age_year, 15, 49),]
drop <- data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
data <- data[!nid %in% drop,]


# Subset to those who've had sex in the last 12 months
data <- data[sex_in_last_12_months == 1, ]

# drop observations with missing condom_last_time or pweight
data <- data[!is.na(pweight) & !is.na(condom_last_time),]

# Only keep surveys that asked both men and women about condom last time
svys = data[,list(sex_split = mean(sex_id, na.rm=T)),by='nid,country']
data <- data[nid %in% svys[sex_split > 1 & sex_split < 2,nid]]

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
saveRDS(data, file = paste0("<<<< FILEPATH REDACTED >>>>", topic, "/","pre_collapse/", indicator, "_pre_collapse_", collapse_version, ".rds"))

## Collapse all data (polygon and point) -----------------------------------------------------------
# collapse data by survey and location
data <- data[, list(int_year = floor(median(int_year, na.rm=T)),
                    condom_last_time = weighted.mean(condom_last_time, pweight),
                    N = sum(pweight)^2/sum(pweight^2),
                    sum_of_sample_weights = sum(pweight),
                    N_obs = .N),
             by='nid,country,source,year,point,shapefile,location_code,latitude,longitude']

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")

# assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]

# Convert to counts and resample the polygon data -------------------------------------------------
data[, condom_last_time := condom_last_time * N]
data <- resample_polygons(data = data, indic = "condom_last_time", cores = cores, pull_poly_method = "fast", seed = 98121)

# use resampling weights to down-weight N and condom_last_time, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, condom_last_time := condom_last_time * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, condom_last_time, weight, sum_of_sample_weights, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# report differences compared to the most recent version
compare_most_recent(data, indicator, prev_version)

## save input data for model runs & log dates-------------------------------------------
save_data(data, indicator, geomatched_version, collapse_version)
