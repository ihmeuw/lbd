####################################################################################################
## Obtain people alive on ART based on PLHIV for each GADM location in UGA
####################################################################################################


rm(list = ls())

libs <- c("data.table", "raster", "foreign", "rgdal", "geosphere", "colorspace", 'readxl',
          "dplyr", "rgeos", "car","plyr", "ggplot2", "scales", "sf", "gridExtra")
sapply(libs, require, character.only = T)

## Read in extraction sheet ------------------------------------------------------------------------
## Load new & old data -----------------------------------------------------------------------------
iso3 <- c("NAM")
country <- c("Namibia")
user <- "USER"
treat_file <- data.table()
treat_file <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
treat_file <- na.omit(treat_file, cols = c('end_year', "location_code"))
col_names <- c("end_year", "shapefile_ref", "location_code", "site_memo", "hiv_treat_count", "age_group_category")
data_pepfar <- treat_file[, col_names, with = FALSE]
data_pepfar$location_code <- as.numeric(as.character(data_pepfar$location_code))

## Add GADM location information and PLHIV columns -------------------------------------------------
new_cols <- c("lbd_locs", "lbd_hiv_pop", "pep_hiv_pop", "plhiv_frac", "lbd_hiv_treat_count")
data_pepfar[, (new_cols) := NA]

look_up_final <- fread("<<<< FILEPATH REDACTED >>>>")
look_up_final$V1 <- NULL # added 2019_11_27
look_up_final$V1 <- NULL # added 2019_11_27
data_pepfar <- merge(look_up_final, data_pepfar, by.x = "other_id", by.y = "location_code", all = T, allow.cartesian = T)

## Add PLHIV data and merge to treatment data ------------------------------------------------------
pop_data <- fread("<<<< FILEPATH REDACTED >>>>")


## Subset to only data years and data country ------------------------------------------------------
data_year <- as.vector(unique(data_pepfar$end_year))
subset_PLHIV <- subset(pop_data, subset = pop_data$year %in% data_year)
subset_PLHIV <- subset(subset_PLHIV, subset = subset_PLHIV$ADM0_NAME %in% country)

data_pepfar[, ADM0_NAME := country]
#setnames(data_pepfar, old = "site_memo", new = "ADM2_NAME")

lbd_loc_only <- distinct(select(subset_PLHIV, ADM2_NAME, ADM2_CODE))

## Merge PLHIV data onto ART extraction sheet ------------------------------------------------------
final_data <- data_pepfar[age_group_category == "adults",]
final_data$hiv_treat_count <- as.numeric(as.character(final_data$hiv_treat_count))
setnames(final_data, old = c("end_year", "other_id"), new = c("year", "location_code"))
required_columns <- c("ADM2_CODE", "year")
for (col in required_columns){
  final_data <- final_data[ !is.na(get(col)), ]
}
final_data_with_pop <- merge(final_data, subset_PLHIV, by = c("ADM2_CODE", "year", "ADM0_NAME"), all = T)
final_data_with_pop[, lbd_locs := ADM2_CODE]
final_data_with_pop[, lbd_hiv_pop := mean]

## Aggregate PLHIV based on pepfar location and year -----------------------------------------------
final_data_with_pop$pep_hiv_pop <- NULL
setnames(final_data_with_pop, old = "location_code", new = "location_code_pep")
final_data_with_pop[, pep_hiv_pop := sum(lbd_hiv_pop, na.rm = T), by = c('year', "location_code_pep", "site_memo")] # 2019_11_27 added site_memo to ensure input sum equals output sum
final_data_with_pop[, plhiv_frac := (lbd_hiv_pop/pep_hiv_pop)]
if ("hiv_treat_count" %in% colnames(final_data_with_pop)) {
  setnames(final_data_with_pop, old = "hiv_treat_count", new = "hiv_treat_count_pep_report")
}

final_data_with_pop[, lbd_hiv_treat_count := (plhiv_frac * hiv_treat_count_pep_report)]


sort(colnames(final_data_with_pop))
extract_temp <- as.data.frame(read_excel("<<<< FILEPATH REDACTED >>>>"))
final_data_with_pop[, site_memo:= paste0(ADM2_NAME, " (Constituency) |",
                                         ADM1_NAME, " (Region) |", " Namibia (Country)")]
final_data_with_pop[, shapefile_ref:= "lbd_standard_admin_2_stage_1"]
final_data_with_pop[, location_name:= "Namibia|NAM"]
final_data_with_pop[, location_id:= "195"]
final_data_with_pop[, ihme_loc_id:= "NAM"]
final_data_with_pop[, admin_level:= 2]
final_data_with_pop[, "shape_type (point or poly)":= 1]
final_data_with_pop[, poly_id_field_name:= "ADM2_CODE"]
setnames(final_data_with_pop, old = c("lbd_hiv_treat_count", "ADM2_CODE", "year"), new = c("hiv_treat_count", "location_code", "end_year"))
coldiffs <- setdiff(colnames(extract_temp), colnames(final_data_with_pop))
final_data_with_pop[, (coldiffs):= "NA"]



## Check column names and same the file ------------------------------------------------------------
colnames(final_data_with_pop)
write.csv(final_data_with_pop, paste0("<<<< FILEPATH REDACTED >>>>"))


