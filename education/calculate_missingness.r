########### FOR DROPPING BY MISSINGNESS ######################

cores <- 15
collapse_new <- TRUE # whether to collapse all extracted data for diagnostics (TRUE) or use most recent saved version (FALSE)

package_list <- c("plyr", "dplyr", "openxlsx", "data.table")
source("<<<< FILEPATH REDACTED >>>>")
load_R_packages(package_list)

if (collapse_new == TRUE) {
  source("<<<< FILEPATH REDACTED >>>>")
  data <- collapse_for_diagnostics(cores, return_data=T, drop_age_missing=T, age_15_to_49=T, sex_separate=T)
  data_version <- format(Sys.time(), "%Y_%m_%d")
} else {
  # read in the most recent data for diagnostics
  diag_files <- sort(grep("*collapsed_for_diagnostics*", list.files("<<<< FILEPATH REDACTED >>>>")))
  file <- diag_files[length(diag_files)]
  data <- readRDS(file)
  data_version <- substr(file, nchar(file)-13, nchar(file)-4)
}

data <- data.table(data)

data$point <- as.numeric(data$point)
data$year <- as.numeric(data$year)
data$nid <- as.numeric(data$nid)

data$location_code <- as.numeric(data$location_code)
data$lat <- as.numeric(data$lat)
data$long <- as.numeric(data$long)
data$point <- as.numeric(data$point)

data$point <- with(data, ifelse(!is.na(lat) & !is.na(long), 1, ifelse(!is.na(shapefile) & shapefile != "", 0, NA)))

data <- data[, list(nid, iso3, year, num_persons, n_missing, survey_series, point, lat, long, shapefile, location_code, sex)]

# iso3 as we've extracted it is actually ihme_loc_id (e.g., could be CHN_3476)
data[, loc_id := iso3]
data[, iso3 := substr(loc_id, 1, 3)]

# subset data to only keep post 1998 and stage 1 & 2 surveys
data <- data[year >= 1998,]
stage_dta <- data.table(read.xlsx(paste0(j, "<<<< FILEPATH REDACTED >>>>")))[,c("iso3", "Stage")]

# only looking at 1998 and later, stage 1 and 2 surveys
data <- merge(data, stage_dta, by="iso3", all.x=T)
data <- data[Stage != "3",]

data <- data[!is.na(sex)]

survey_data <- data[,.(N=sum(num_persons),prop_missing=sum(n_missing, na.rm=T)/sum(num_persons, na.rm=T),
                       survey_series=unique(survey_series[!is.na(survey_series)])[1],
                       point=mean(point, na.rm=T)),
                    by=.(year,iso3,nid)]


missing_threshold <- 0.5

svy_sex_missing <- data[,.(N=sum(num_persons),prop_missing=sum(n_missing, na.rm=T)/sum(num_persons, na.rm=T),
                           survey_series=unique(survey_series[!is.na(survey_series)])[1],
                           point=mean(point, na.rm=T)),
                        by=.(year,iso3,nid, sex)]
svy_sex_missing[, drop := ifelse(prop_missing >= missing_threshold, 1, 0)]
svy_sex_missing[, sex := ifelse(sex == 1, "mn", "wn")]



drop_tracker <- svy_sex_missing[,list(year, iso3, nid, survey_series, drop, sex)]
drop_tracker <- reshape(drop_tracker, idvar=c("year", "nid", "iso3", "survey_series"), v.names="drop", timevar="sex", direction="wide")
over_threshold_svys <- survey_data[prop_missing >= missing_threshold, nid]
drop_tracker[, drop.both := ifelse(nid %in% over_threshold_svys, 1, 0)]

write.csv(drop_tracker, "<<<< FILEPATH REDACTED >>>>")
