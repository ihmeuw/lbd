####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, and resample
##              polygons.
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list = ls())

# set arguments for this run and user
core_repo <- paste0('<<<< FILEPATH REDACTED >>>>', '/lbd_core/')
indic_repo <- paste0('<<<< FILEPATH REDACTED >>>>', '/lbd_hiv/')
cores <- 10

# find most recent versions
most_recent_date <- function(dir, date_format = "%Y_%m_%d", out_format = "%Y_%m_%d", file_pattern = NULL) {
  date_pattern <- gsub("%y|%m|%d", "[[:digit:]]{2}", gsub("%Y", "[[:digit:]]{4}", date_format))
  dates <- dir(dir, pattern = date_pattern)
  if (!is.null(file_pattern)) dates <- grep(file_pattern, dates, value = T)
  dates <- gsub(paste0("(.*)(", date_pattern, ")(.*)"), "\\2", dates)
  dates <- as.Date(dates, date_format)
  format(max(dates), out_format)
}

geomatch_date <- most_recent_date("<<<< FILEPATH REDACTED >>>>")
anc_date <- most_recent_date("<<<< FILEPATH REDACTED >>>>")
collapse_date <- format(Sys.time(), "%Y_%m_%d")
report_commit <- system(paste0('cd ', indic_repo, ';git log -1 --format="%h"'), intern = T)

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = c(package_list, 'openxlsx'), repos = c(core_repo, indic_repo))

setompthreads(cores)

## Load and subset data ----------------------------------------------------------------------------
# load geomatched data
data <- readRDS("<<<< FILEPATH REDACTED >>>>")
data <- data.table(data)
data_for_crosswalk <- copy(data) # used for cross-walking function

# drop unneeded variables and remove annoying stata attributes
data <- data[, list(nid, country, survey_series, survey_name, year,
                    strata, psu, point, shapefile, location_code, latitude, longitude,
                    sex_id, age_year, int_year, hiv_test, hiv_weight)]
for (ii in 1:ncol(data)) attributes(data[[ ii]]) <- NULL

# rename variables
setnames(data, c("survey_series"), c("source"))

# fix variable class
data[, point := as.numeric(point)]
data[, latitude := as.numeric(latitude)]
data[, longitude := as.numeric(longitude)]

# subset to ages 15-49 (and only surveys with this full range)
data <- data[between(age_year, 15, 49),]
drop <- data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
data <- data[!nid %in% drop,]

# drop any surveys that are sex-specific
drop <- data[, list(sum(sex_id == 1, na.rm = T), sum(sex_id == 2, na.rm = T)), nid][V1 == 0 | V2 == 0, nid]
data <- data[!nid %in% drop,]

# drop observations with missing hiv test result or hiv weight
data <- data[!is.na(hiv_test) & !is.na(hiv_weight),]

# drop PSUs where all hiv weights are 0 (otherwise this causes NaNs)
drop <- data[, sum(hiv_weight), by = 'nid,psu'][V1 == 0, list(nid, psu)]
data <- data[!drop, on = c('nid', 'psu')]

# drop point data with missing latitude/longitude
data <- data[(!is.na(latitude) & !is.na(longitude)) | point == 0,]

# drop polygon data with missing shapefile/location code
data <- data[(!is.na(shapefile) & !is.na(location_code)) | point == 1,]

# drop countries outside Africa
loc <- source("<<<< FILEPATH REDACTED >>>>/get_location_metadata.R")
loc <- data.table(get_location_metadata(location_set_id = 2, gbd_round_id = 4))
loc <- loc[grepl("Africa", region_name) & location_type == "admin0", ihme_loc_id]
data <- data[country %in% loc,]

## Collapse all data (polygon and point) -----------------------------------------------------------
# collapse data by survey and location
data <- data[, list(int_year = floor(median(int_year, na.rm = T)),
                    hiv_test = weighted.mean(hiv_test, hiv_weight),
                    N = sum(hiv_weight) ^ 2 / sum(hiv_weight ^ 2),
                    N_obs = .N,
                    sum_of_sample_weights = sum(hiv_weight)),
             by = 'nid,country,source,year,point,shapefile,location_code,latitude,longitude']

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")

# note the data type
data[, type := "Survey microdata"]

## Load and subset survey report data --------------------------------------------------------------
# load report data
report <- fread("<<<< FILEPATH REDACTED >>>>")
report <- report[, list(country, nid, location_type, survey_series, survey_name, int_year,
                        start_age, end_age, sex_id, point, shapefile, location_code, latitude,
                        longitude, hiv_test, hiv_test_lb, hiv_test_ub, N)]

# drop national estimates, countries outside Africa, and NIDs we have microdata for
report <- report[location_type != "National",]
report <- report[country %in% loc]
report <- report[!nid %in% data$nid]

# keep only for both sexes combined
report[, id := paste(nid, country, int_year)]
by_sex_only <- setdiff(intersect(report[sex_id == 1, id], report[sex_id == 2, id]), report[sex_id == 3, id])
if (length(by_sex_only) > 0) warning(paste("Report data dropped because only sex-specific estimates were available in", length(by_sex_only), "source-country-years"))
report <- report[sex_id == 3]

# change hiv_test from a percentage to a proportion
report[, `:=` (hiv_test = hiv_test/100, hiv_test_lb = hiv_test_lb/100, hiv_test_ub = hiv_test_ub/100)]

# estimate N if only CIs are available (using a normal approximation, assuming 95% CIs)
report[is.na(N) & !is.na(hiv_test_lb) & !is.na(hiv_test_ub),
       N := (hiv_test * (1 - hiv_test)) / (((hiv_test_ub - hiv_test_lb) / (2 * 1.96))^2)]
ci0 <- report[is.na(N) & hiv_test_lb == 0 & hiv_test_ub == 0, .N]
if (ci0 > 0) warning(paste("Report data dropped because no sample size was provided and lower and upper bound are both 0 for", ci0, "source-country-year-locations"))
report <- report[!(is.na(N) & hiv_test_lb == 0 & hiv_test_ub == 0),]
report[, c("hiv_test_lb", "hiv_test_ub") := NULL]

# drop alternate ages if 15 to 49 is available
age_15_49 <- report[start_age == 15 & end_age == 49, id]
report <- report[(start_age == 15 & end_age == 49) | !id %in% age_15_49,]

# if there are multiple non-15-49 age groups, keep only the one with the largest overlap
mult_ages <- report[, uniqueN(paste(start_age, end_age)), id][V1 > 1, id]
if (length(mult_ages) > 0) warning(paste("Report data with multiple age groups identified in", length(mult_ages), "source-country-years; all but one age group has been dropped."))
report[, overlap := length(intersect(15:49, start_age:end_age)), by = seq_len(nrow(report))]
report[, max_overlap := max(overlap), id]
report <- report[overlap == max_overlap,]

## Load and subset lit review data -----------------------------------------------------------------
# load lit data
lit <- read.xlsx('<<<< FILEPATH REDACTED >>>>')
lit <- data.table(lit[-1, ])
lit <- lit[!is.na(nid), ]

# subset to rows marked for inclusion
lit <- lit[use_for_adult_prev == 1, ]

# drop national estimates, countries outside Africa, and NIDs we have microdata for
lit <- lit[substr(ihme_loc_id, 1, 3) %in% loc,]
lit <- lit[!nid %in% data$nid,]

# reformat variables
lit <- lit[, list(nid, country = ihme_loc_id, shapefile = poly_reference, location_code = poly_id,
                  latitude = lat, longitude = long, sex, year_start, year_end, age_start, age_end,
                  mean, lower, upper, cases, sample_size)]
lit <- mutate_if(lit, is.factor, as.character)
lit <- mutate_at(lit, c('nid', 'location_code', 'age_start', 'age_end', 'year_start', 'year_end'), as.integer)
lit <- mutate_at(lit, c('latitude', 'longitude', 'mean', 'lower', 'upper', 'cases', 'sample_size'), as.numeric)
lit <- data.table(lit)

lit[, shapefile := gsub('.shp$', '', shapefile)]
lit[, point := as.numeric(!is.na(latitude))]

# exclude if there is missing location info
lit <- lit[(!is.na(latitude) & !is.na(longitude)) | point == 0,]
lit <- lit[(!is.na(shapefile) & !is.na(location_code)) | point == 1,]

# fill in 'sample size' from uncertainty intervals, or drop if both are unavailable
lit[is.na(sample_size), sample_size := (mean * (1 - mean)) / (((upper - lower) / (2 * 1.96))^2)]
if (lit[, sum(is.na(sample_size))] > 0) warnings('lit data dropped because no sample size or lower and upper bounds were available')
lit <- lit[!is.na(sample_size), ]

# if a mean is provided, use this to repopulate cases (in case sample weights were incorporated)
lit[!is.na(mean), cases := mean * sample_size]

# use the midpoint for the time localization
lit[, year := floor((year_start + year_end)/2)]

# combine ages, sexes, and HIV infection types as needed
if (lit[, .N, by = 'country,nid,shapefile,location_code,latitude,longitude,year,sex,age_start'][N > 1, .N] > 0) stop("Unexplained duplicate rows in lit review")
lit <- lit[, list(hiv_test = sum(cases) / sum(sample_size), N = sum(sample_size), age_start = min(age_start), age_end = max(age_end)),
           by = 'country,nid,point,shapefile,location_code,latitude,longitude,year']

## Age crosswalk report and lit extraction data ----------------------------------------------------
# combine report and lit data
setnames(report, 'int_year', 'year')
report[, type := "Survey report"]

setnames(lit, c('age_start', 'age_end'), c('start_age', 'end_age'))
lit[, c('survey_name', 'survey_series') := "COUNTRY_SPECIFIC"]
lit[, type := "Literature review"]

if (length(intersect(report$nid, lit$nid)) > 0) stop('duplicated NIDs in report and lit review extractions')
add_data <- rbind(lit, report[, names(lit), with = F])
rm(lit, report)

# remove survey report and lit extractions that don't overlap 15-49
add_data <- add_data[(start_age <= 15 & end_age >= 15) | (start_age >= 15 & start_age <= 49), ]

# identify all age ranges other than 15-49, subset report and lit extraction to covering those age ranges
age_list <- unique(add_data[!(start_age == 15 & end_age == 49), c("start_age", "end_age")])
crosswalk_report_list <- apply(age_list, 1, function(x) return(add_data[start_age == x[1] & end_age == x[2]]))

# prepare microdata for crosswalk
data_for_crosswalk <- crosswalk_data_prepare(data_for_crosswalk)

# apply crosswalk
crosswalked <- lapply(crosswalk_report_list, function(x) crosswalk_reports(data_for_crosswalk, x))
crosswalked <- data.table(bind_rows(crosswalked))
setnames(crosswalked, c("prev_old", "prev_new"), c("hiv_old", "hiv_new"))
rm(data_for_crosswalk)

# drop r_square column
crosswalked$r_square <- NULL

# recombine cross-walked report data with data for 15-49
crosswalked[, c("hiv_test", "N") := list(hiv_new, N_eff)]
add_data <- rbind(add_data[start_age == 15 & end_age == 49,], crosswalked, fill = T)

## Combine report & lit review data with collapsed microdata ---------------------------------------
# approximate the effective sample size using the design effects from the polygon microdata
# Note: this is *very* rough; the purpose is to make sure the report/lit data doesn't have an unfair
# advantage compared to microdata; the calculation below excludes surveys where there is no apparent
# design effect -- this is not super plausible, and we don't want these to skew the distribution.
de <- data[point == 0, list(de = sum(N)/sum(N_obs)), by = 'nid,shapefile,location_code'][de < 0.99999, median(de)]
add_data[, N_obs := N]
add_data[, N := N_obs * de]

# update the survey list
nid_list <- c(nid_list, unique(add_data$nid))

# combine with collapsed microdata
add_data <- add_data[, list(nid, country, survey_series, type, point, shapefile, location_code, latitude, longitude, year, hiv_test, N, N_obs)]
setnames(add_data, "survey_series", "source")
add_data[, sum_of_sample_weights := N] # this is a poor approximation...really only a placeholder
data <- rbind(data, add_data, fill = T)
rm(add_data)

## Add ANC data ------------------------------------------------------------------------------------
# load anc data
anc <- readRDS("<<<< FILEPATH REDACTED >>>>")
anc <- data.table(mutate_if(anc, is.factor, as.character))
setnames(anc, "anc_hiv_test", "hiv_test")

# remove potentially problematic non-standard characters in site names
for (i in 1:nrow(anc)) {
  if (class(try(nchar(anc[i, site]), silent = T)) == 'try-error') {
    prob <- which(sapply(1:50, function(x) class(try(substr(anc[i, site], x, x), silent = T))) == 'try-error')
    anc[i, site := paste0(substr(site, 1, prob - 1), 'XX')]
  }
}

# combine with survey data
anc[, sum_of_sample_weights := N] # sample weights are not a thing for ANC, but for completeness...
anc[, type := "ANC"]
data <- rbind(data, anc[, c(names(data), "site"), with = F], fill = T)

# drop north Africa data
north_africa <- c("ESH","LBY","EGY","DZA","TUN", "MAR")
data <- subset(data, !(country %in% north_africa))

# now that all data are compiled, assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]

## Convert to counts and resample the polygon data -------------------------------------------------
# run polygon resampling
data[, hiv_test := hiv_test * N]
data <- resample_polygons(data = data, indic = "hiv_test", cores = cores, pull_poly_method = "fast", seed = 98121)

# use resampling weights to down-weight N and hiv_test, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, hiv_test := hiv_test * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

# check for missing lat/long or weights
if (sum(is.na(data$latitude)) > 0 | sum(is.na(data$longitude)) > 0 | sum(is.na(data$weight)) > 0 | sum(is.na(data$point)) > 0) stop("Issue with resample_polygons.")
if (sum(is.na(data$hiv_test)) > 0 | sum(is.na(data$nid)) > 0 | sum(is.na(data$cluster_id)) > 0 | sum(is.na(data$year)) > 0) stop("There is missingness in hiv_test, nid, cluster_id, or year.")

## Save collapsed data -----------------------------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, site, year, latitude, longitude, N, N_obs, hiv_test, weight, type, point, cluster_id, sum_of_sample_weights)]
setkey(data, nid, country, source, site, year, latitude, longitude, N, hiv_test)

# add dates of geomatched data (geomatched_version) & collapsed data (collapse_version)
data$geomatch_date <- geomatch_date
data$report_commit <- report_commit
data$anc_date <- anc_date
data$collapse_date <- collapse_date

# write files
for (save_type in c("survey", "anc", "all", "no_poly")) {
  if (save_type == "survey") {
    all_collapsed <- data[type != "ANC", ]
    file_tag <- "_survey"
  }
  if (save_type == "anc") {
    all_collapsed <- data[type == "ANC", ]
    file_tag <- "_anc"
  }
  if (save_type == "all") {
    all_collapsed <- copy(data)
    file_tag <- NULL
  }
  if (save_type == "no_poly") {
    all_collapsed <- data[point == 1, ]
    file_tag <- "_no_poly"
  }

  write.csv(all_collapsed, file = paste0("<<<< FILEPATH REDACTED >>>>/hiv_test", file_tag, ".csv"), row.names = FALSE)
  saveRDS(all_collapsed, file = paste0("<<<< FILEPATH REDACTED >>>>/hiv_test", file_tag, "_", format(Sys.Date(), "%Y_%m_%d"), ".rds"))
}
