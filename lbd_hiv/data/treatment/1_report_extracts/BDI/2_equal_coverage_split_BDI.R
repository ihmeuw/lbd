####################################################################################################
## Obtain people alive on ART based on PLHIV for each GADM location in BDI

####################################################################################################


rm(list = ls())

libs <- c("data.table", "raster", "foreign", "rgdal", "geosphere", "colorspace", "readxl",
          "dplyr", "rgeos", "car","plyr", "ggplot2", "scales", "sf", "gridExtra")
sapply(libs, require, character.only = T)

## Read in extraction sheet ------------------------------------------------------------------------
## Load new & old data -----------------------------------------------------------------------------
iso3 <- c("BDI")
country <- c("Burundi")
treat_file <- data.table()
treat_file <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
treat_file <- na.omit(treat_file, cols = c('end_year', "other_id"))
treat_file$V1 <- NULL
treat_file$occ_id = NULL
## Add PLHIV data and merge to treatment data ------------------------------------------------------
pop_data <- fread("<<<< FILEPATH REDACTED >>>>")

## Subset to only data years and data country ------------------------------------------------------
data_year <- as.vector(unique(treat_file$end_year))
subset_PLHIV <- subset(pop_data, subset = pop_data$year %in% data_year)
subset_PLHIV <- subset(subset_PLHIV, subset = subset_PLHIV$ADM0_NAME %in% country)

treat_file[, ADM0_NAME := country]

## Merge PLHIV data onto ART extraction sheet ------------------------------------------------------
final_data <- treat_file[age_group_category == "adults",]
final_data <- merge(subset_PLHIV, final_data, by.x = c("ADM2_CODE", "year", "ADM0_NAME"), by.y = c("location_code", "end_year", "ADM0_NAME"), all.y = T, allow.cartesian = T)
## Aggregate PLHIV based on pepfar location and year -----------------------------------------------
final_data[, pep_hiv_pop := sum(mean, na.rm = T), by = c('year', "other_id")]
final_data[, plhiv_frac := (mean/pep_hiv_pop)]

final_data[, lbd_hiv_treat_count := (final_data$plhiv_frac * final_data$hiv_treat_count)]

#final_data[year == 2012 | year == 2013, lbd_hiv_treat_count:= hiv_treat_count] ## for 2012 and 2013 the counts don't change (already matched to admin 2)
setnames(final_data, old = c("lbd_hiv_treat_count", "hiv_treat_count", "year", "ADM2_CODE"), new = c("hiv_treat_count", "hiv_treat_count_report", "end_year", "location_code"))
## Check column names and same the file ------------------------------------------------------------
colnames(final_data)
write.csv(final_data, paste0("<<<< FILEPATH REDACTED >>>>"))



