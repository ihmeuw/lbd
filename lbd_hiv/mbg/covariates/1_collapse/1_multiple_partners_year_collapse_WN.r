####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   multiple sexual partners in the last year among women, ages 15-49
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- '<<<< FILEPATH REDACTED >>>>/lbd_hiv/'
cores <- 10

# set arguments for this indicator/topic
indicator <- "multiple_partners_year_WN"
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
                    sex_id, age_year, int_year, pweight, num_partners_year)]
for (ii in 1:ncol(data)) attributes(data[[ ii]]) <- NULL

# rename variables
setnames(data, c("survey_series"), c("source"))

# fix variable class
data[, year := as.numeric(year)]
data[, point := as.numeric(point)]
data[, latitude := as.numeric(latitude)]
data[, longitude := as.numeric(longitude)]

# only want surveys 1998 or later for now
data <- data[year >= 1998 & sex_id == 2,]

# drop all DHS3 surveys (plus DHS4 GIN & DHS4 GHA which also required a spousal adjust)
# for differing question construction (see docs)
dhs3_nids <- c(19076, 18531, 20537, 20132, 19198, 20909, 20796, 19305, 21090, 20382,
               20780, 20212, 20301, 19493, 19546, 20852, 18938, 19370, 19292, 20976,
               19604, 21139, 18519, 20120, 19670, 19614)
data <- data[!(nid %in% dhs3_nids),]

# drop GAB DHS 2000 because of data oddities
data <- data[nid != 19579]

data[, multiple_partners_year := ifelse(num_partners_year <= 1, 0, 1)]

# subset to ages 15-49 (and only surveys with this full range)
data <- data[between(age_year, 15, 49),]
drop <- data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
data <- data[!nid %in% drop,]

# drop observations with missing multiple_partners_year or pweight
data <- data[!is.na(pweight) & !is.na(multiple_partners_year),]

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
# collapse data by survey and location
data <- data[, list(int_year = floor(median(int_year, na.rm=T)),
                    multiple_partners_year = weighted.mean(multiple_partners_year, pweight),
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
data[, multiple_partners_year := multiple_partners_year * N]
data <- resample_polygons(data = data, indic = "multiple_partners_year", cores = cores, pull_poly_method = "fast", seed = 98121)

# use resampling weights to down-weight N and multiple_partners_year, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, multiple_partners_year := multiple_partners_year * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, multiple_partners_year, weight, sum_of_sample_weights, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# report differences compared to the most recent version
compare_most_recent(data, indicator, prev_version)


## save input data for model runs & log dates-------------------------------------------
save_data(data, indicator, geomatched_version, collapse_version)

# changed
#353526,  287639,  77384,  59339, 270627

# dropped: 3935,91506
