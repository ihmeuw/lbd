####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   married or living as married among men and women, ages 15-49
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- '<<<< FILEPATH REDACTED >>>>/lbd_hiv/'
cores <- 10

# set arguments for this indicator/topic
indicator <- "in_union"
topic <- "rsp"

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>/package_list.csv"))
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
                    sex_id, age_year, int_year, pweight, marital_status, survey_module)]
for (ii in 1:ncol(data)) attributes(data[[ ii]]) <- NULL

# rename variables
setnames(data, c("survey_series"), c("source"))

# fix variable class
data[, year := as.numeric(year)]
data[, point := as.numeric(point)]
data[, latitude := as.numeric(latitude)]
data[, longitude := as.numeric(longitude)]

# only want surveys 1998 and later for now
data <- data[year >= 1998,]

# drop PMA WN mods to prevent duplication (see docs)
data <- data[!(grepl("*PMA2020*", survey_name) & survey_module == "WN"),]

# drop DHS WN mods where HHM mod was extracted to prevent duplication (see docs)
data <- data[!((nid == 19539 | nid == 111432) & survey_module == "WN")]

data[, survey_module := NULL]

# drop COD PMA 2013-14, 2014, 2015, 2015-2016, 2016 surveys because of data oddities (see docs)
data <- data[!nid %in% c(257822, 257823, 257826, 286019, 286054, 286020),]

# drop MLI LSMS 2014-2015, KEN WHO STEPS_NCD 2015,
# NAM HH INCOME & EXPENDITURE 2009-2010, UGA LSMS ISA 2009-2010,
# BEN HH Living survey 2009 for data oddities (see docs)
data <- data[!nid %in% c(260407, 165492, 134371, 81004, 151768)]

# drop ZAF ISSP surveys for data oddities (see docs)
data <- data[!(source == "ISSP" & country == "ZAF")]

# create a new indicator in_union that's true for respondents
# who are currently married/living with partner, false otherwise
data[, in_union := ifelse(marital_status < 3, 1, 0)]

# drop observations with missing in_union or pweight
data <- data[!is.na(pweight) & !is.na(in_union),]

# Only keep surveys that asked both men and women about marital status
svys = data[,list(sex_split = mean(sex_id, na.rm=T)),by='nid,country']
data <- data[nid %in% svys[sex_split > 1 & sex_split < 2,nid]]

# subset to ages 15-49 (and only surveys with this full range)
data <- data[between(age_year, 15, 49),]
drop <- data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
data <- data[!nid %in% drop,]

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

data <- data[, .(int_year = floor(median(int_year, na.rm = T)),
                    in_union = weighted.mean(in_union, pweight),
                    N = sum(pweight) ^ 2 / sum(pweight ^ 2),
                    N_obs = .N,
                    sum_of_sample_weights = sum(pweight)),
             by=.(nid, country, source, year, point, shapefile, location_code, latitude, longitude)]

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")

# assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]

# Convert to counts and resample the polygon data -------------------------------------------------
data[, in_union := in_union * N]
data <- resample_polygons(data = data, indic = "in_union", cores = cores, pull_poly_method = "fast", seed = 98121)

# use resampling weights to down-weight N and in_union, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, in_union := in_union * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, in_union, weight, sum_of_sample_weights, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# # report differences compared to the most recent version
# compare_most_recent(data, paste0(indicator, "_nonconservative"), prev_version)
#
# ## save input data for model runs & log dates-------------------------------------------
# save_data(data, paste0(indicator, "_nonconservative"), geomatched_version, collapse_version)

# drop nids that are borderline outliers for conservative version
data <- data[!nid %in% c(1912, 22950, 7375, 224223, 151797, 1981, 2209)]

# report differences compared to the most recent version
#compare_most_recent(data, paste0(indicator, "_conservative"), prev_version)
old <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>",
                      indicator, "_conservative_", prev_version, ".rds"))
both_names <- intersect(names(old), names(data))
cat(paste0("\n\n*****Columns have been added: \n", paste(setdiff(names(data),
                                                                 names(old)), collapse = ", "), "\n\n"))
old <- old[, both_names, with = F]
setkeyv(old, key(data))
cat(paste0("\n\n*****NIDs that have been added:\n", paste(setdiff(data$nid,
                                                                  old$nid), collapse = ","), "\n\n"))
cat(paste0("\n\n*****NIDs that have been dropped:\n", paste(setdiff(old$nid,
                                                                    data$nid), collapse = ","), "\n\n"))
indicator <- gsub("_WN|_MN|_BOTH", "", indicator)
cat("\n\n*****NIDs that have changed:\n")
for (nn in intersect(data$nid, old$nid)) {
  temp_old <- old[nid == nn, mget(c("nid", "country",
                                    "source", "year", "latitude", "longitude", "N",
                                    indicator, "weight"))]
  temp_new <- data[nid == nn, mget(c("nid", "country",
                                     "source", "year", "latitude", "longitude", "N",
                                     indicator, "weight"))]
  if (all.equal(temp_old, temp_new)[1] != "TRUE") {
    cat(paste("\n\n", nn, "- data changed\n"))
    print(all.equal(data.frame(temp_old), data.frame(temp_new)))
  }
}

save_data(data, paste0(indicator, "_conservative"), geomatched_version, collapse_version)

