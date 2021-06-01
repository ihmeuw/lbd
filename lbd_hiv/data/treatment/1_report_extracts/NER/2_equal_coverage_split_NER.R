####################################################################################################
## Obtain people alive on ART based on PLHIV for each GADM location in NER
####################################################################################################

rm(list = ls())

libs <- c("data.table", "raster", "foreign", "rgdal", "geosphere", "colorspace", "fossil",
          "dplyr", "rgeos", "car","plyr", "ggplot2", "scales", "sf", "gridExtra")
sapply(libs, require, character.only = T)

## Read in extraction sheet ------------------------------------------------------------------------
## Load new & old data -----------------------------------------------------------------------------
iso3 <- c("NER")
country <- c("Niger")
user <- "USER"
treat_file <- data.table()
treat_file <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
treat_file <- na.omit(treat_file, cols = c('end_year', "location_code"))
col_names <- c("end_year", "location_code", "hiv_treat_count", "age_group_category")
data_adm1 <- treat_file
data_adm1$location_code <- as.numeric(as.character(data_adm1$location_code))

## Add GADM location information and PLHIV columns -------------------------------------------------
new_cols <- c("lbd_locs", "lbd_hiv_pop", "pep_hiv_pop", "plhiv_frac", "lbd_hiv_treat_count")
data_adm1[, (new_cols) := NA]
data_adm1$V1 <- NULL
look_up_final <- fread("<<<< FILEPATH REDACTED >>>>")
look_up_final$V1 <- NULL
data_adm1 <- merge(look_up_final, data_adm1, by.x = "other_id", by.y = "location_code", all.x = T, allow.cartesian = T)

## Add PLHIV data and merge to treatment data ------------------------------------------------------
pop_data <- fread("<<<< FILEPATH REDACTED >>>>")

## Subset to only data years and data country ------------------------------------------------------
data_year <- as.vector(unique(data_adm1$end_year))
subset_PLHIV <- subset(pop_data, subset = pop_data$year %in% data_year)
subset_PLHIV <- subset(subset_PLHIV, subset = subset_PLHIV$ADM0_NAME %in% country)

data_adm1[, ADM0_NAME := country]

## Merge PLHIV data onto ART extraction sheet ------------------------------------------------------
final_data <- data_adm1[age_group_category == "adults",]
final_data$hiv_treat_count <- as.numeric(as.character(final_data$hiv_treat_count))
setnames(final_data, old = c("end_year", "other_id"), new = c("year", "location_code"))
required_columns <- c("ADM2_CODE", "year")
for (col in required_columns){
  final_data <- final_data[ !is.na(get(col)), ]
}
final_data_with_pop <- merge(final_data, subset_PLHIV, by = c("ADM2_CODE", "year"), all.x = T)
final_data_with_pop[, lbd_locs := ADM2_CODE]
final_data_with_pop[, lbd_hiv_pop := mean]

## Aggregate PLHIV based on adm1 location and year -----------------------------------------------
final_data_with_pop$pep_hiv_pop <- NULL
setnames(final_data_with_pop, old = "location_code", new = "location_code_adm1")
final_data_with_pop[, adm1_hiv_pop := sum(lbd_hiv_pop, na.rm = T), by = c('year', "location_code_adm1", "Sites")]
final_data_with_pop[, plhiv_frac := (lbd_hiv_pop/adm1_hiv_pop)]
if ("hiv_treat_count" %in% colnames(final_data_with_pop)) {
  setnames(final_data_with_pop, old = "hiv_treat_count", new = "hiv_treat_count_adm1")
}

final_data_with_pop[, lbd_hiv_treat_count := (plhiv_frac * hiv_treat_count_adm1)]

## Fill in the rest of the extraction sheet and save the file ------------------------------------------------------------
sort(colnames(final_data_with_pop))
extract_temp <- as.data.frame(read_excel("<<<< FILEPATH REDACTED >>>>"))
final_data_with_pop[, site_memo:= paste0(Sites, " (Health Facility) |",
                                         ADM2_NAME, " (Department) |",
                                         ADM1_NAME.x, " (Region) |", " Niger (Country)")]
final_data_with_pop[, shapefile_ref:= "lbd_standard_admin_2_stage_1"]
final_data_with_pop[, location_name:= "Niger|NER"]
final_data_with_pop[, location_id:= "213"]
final_data_with_pop[, ihme_loc_id:= "NER"]
final_data_with_pop[, admin_level:= 2]
final_data_with_pop[, "shape_type (point or poly)":= 1]
final_data_with_pop[, poly_id_field_name:= "ADM2_CODE"]
setnames(final_data_with_pop, old = c("lbd_hiv_treat_count", "ADM2_CODE", "year"), new = c("hiv_treat_count", "location_code", "end_year"))
coldiffs <- setdiff(colnames(extract_temp), colnames(final_data_with_pop))
final_data_with_pop[, (coldiffs):= "NA"]

write.csv(final_data_with_pop, paste0("<<<< FILEPATH REDACTED >>>>"))



