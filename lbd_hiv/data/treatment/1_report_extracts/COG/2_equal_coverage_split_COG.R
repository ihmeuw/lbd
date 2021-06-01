####################################################################################################
## Obtain people alive on ART based on PLHIV for each GADM location in COG
####################################################################################################

rm(list = ls())
libs <- c("data.table", "raster", "foreign", "rgdal", "geosphere", "colorspace", "readxl",
          "dplyr", "rgeos", "car","plyr", "ggplot2", "scales", "sf", "gridExtra")
sapply(libs, require, character.only = T)

###################################################################################################
## 0. NEED TO CHANGE THESE VARIABLES BEFORE THIS CAN BE RUN
iso3 <- c("COG")
country <- c("Republic of Congo")
username <- "USER"
treat_input <- paste0("<<<< FILEPATH REDACTED >>>>")
lookup_input <- paste0("<<<< FILEPATH REDACTED >>>>")
# In-script things to change ##
# - Columns to remove in section 2 (set to NULL)
# - PLHIV model run (if applicable)
# - Country-specific information in section 7
# - IF input != output, check aggregation groups in section 6 first (default  = )
###################################################################################################

equal_coverage_split <- function(country, iso3, user, treat_input, lookup_input){
  ## 1.Read in extraction sheet, load new and old data, create new data.table -----------------------
  treat_file <- data.table()
  treat_file <- fread(treat_input)
  treat_file <- na.omit(treat_file, cols = c('end_year', "location_code"))
  col_names <- c("end_year", "location_code", "hiv_treat_count", "age_group_category")
  data_adm1 <- treat_file ###You can change the new variable name here based on the granularity of the initial data
  data_adm1$location_code <- as.numeric(as.character(data_adm1$location_code))
  data_adm1$occ_id = NULL
  if ("V1" %in% colnames(data_adm1)){
    data_adm1$V1 = NULL
  }
  if ("ADM2_NAME"  %in% colnames(data_adm1)){
    data_adm1$ADM2_NAME = NULL
  }

  ## 2.Add GADM location information and PLHIV columns from the lookup table to the data.table ------
  new_cols <- c("lbd_locs", "lbd_hiv_pop", "plhiv_frac", "lbd_hiv_treat_count")
  data_adm1[, (new_cols) := NA]
  look_up_final <- fread(lookup_input)
  look_up_final = subset(look_up_final, select = c("ADM2_CODE", "other_id"))
  data_adm1 <- merge(look_up_final, data_adm1, by.x = "other_id", by.y = "location_code", all.x = T, allow.cartesian = T)


  ## 3.load PLHIV data from most recent model run ---------------------------------------------------
  pop_data <- fread("<<<< FILEPATH REDACTED >>>>")

  ## 4.Subset PLHIV to only data years and data country ---------------------------------------------
  data_year <- as.vector(unique(data_adm1$end_year))
  subset_PLHIV <- subset(pop_data, subset = pop_data$year %in% data_year)
  subset_PLHIV <- subset(subset_PLHIV, subset = subset_PLHIV$ADM0_NAME %in% country)
  #data_adm1[, ADM0_NAME := country]

  ## 5.Merge PLHIV data onto ART extraction sheet ---------------------------------------------------
  final_data <- data_adm1[age_group_category == "adults",]
  final_data$hiv_treat_count <- as.numeric(as.character(final_data$hiv_treat_count))
  setnames(final_data, old = c("end_year", "other_id"), new = c("year", "location_code"))
  required_columns <- c("ADM2_CODE", "year")
  for (col in required_columns){
    final_data <- final_data[ !is.na(get(col)), ]
  }
  final_data <- merge(final_data, look_up_final, by = c("ADM2_CODE"), allow.cartesian = T)
  final_data_with_pop <- merge(final_data, subset_PLHIV, by = c("ADM2_CODE", "year"), all.x = T)
  final_data_with_pop[, lbd_locs := ADM2_CODE]
  final_data_with_pop[, lbd_hiv_pop := mean]

  ## 6.Aggregate PLHIV based on pepfar location and year --------------------------------------------
  setnames(final_data_with_pop, old = "location_code", new = "location_code_adm1")
  final_data_with_pop[, adm1_hiv_pop := sum(lbd_hiv_pop, na.rm = T), by = c('year', "location_code_adm1", 'site_memo')] # If under- or over-counting in final output, check this line first
  final_data_with_pop[, plhiv_frac := (lbd_hiv_pop/adm1_hiv_pop)]
  if ("hiv_treat_count" %in% colnames(final_data_with_pop)) {
    setnames(final_data_with_pop, old = "hiv_treat_count", new = "hiv_treat_count_adm1")
  }
  final_data_with_pop[, lbd_hiv_treat_count := (plhiv_frac * hiv_treat_count_adm1)]

  ## 7.Fill in the rest of the extraction sheet and save the file -----------------------------------
  sort(colnames(final_data_with_pop))
  extract_temp <- as.data.frame(read_excel(paste0("<<<< FILEPATH REDACTED >>>>")))
  final_data_with_pop[, site_memo:= paste0(ADM2_NAME, " (District) |",
                                           ADM1_NAME, " (Department) |", " Congo (Country)")]
  final_data_with_pop[, shapefile_ref:= "lbd_standard_admin_2_stage_1"]
  final_data_with_pop[, location_name:= "Congo|COG"]
  final_data_with_pop[, location_id:= 170]
  final_data_with_pop[, ihme_loc_id:= "COG"]
  final_data_with_pop[, admin_level:= 2]
  final_data_with_pop[, "shape_type (point or poly)":= 1]
  final_data_with_pop[, poly_id_field_name:= "ADM2_CODE"]
  setnames(final_data_with_pop, old = c("lbd_hiv_treat_count", "ADM2_CODE", "year"), new = c("hiv_treat_count", "location_code", "end_year"))
  coldiffs <- setdiff(colnames(extract_temp), colnames(final_data_with_pop))

  ## 8.Initial input vs final output check ---------------------------------------------------------
  treat_file_check = na.omit(treat_file, cols = c('hiv_treat_count'))
  treat_file_check$hiv_treat_count = as.numeric(as.character(treat_file_check$hiv_treat_count))
  final_data_with_pop_check = na.omit(final_data_with_pop, cols = c('hiv_treat_count'))
  if (sum(treat_file_check[treat_file_check$age_group_category == "adults", "hiv_treat_count"]) == sum(as.numeric(final_data_with_pop_check$hiv_treat_count))){
    print("Input treatment count equals output treatment count")
  }  else {
    stop("Equal Coverage did not work! Input treatment count does not equal the output treatment count!")
  }

  ## 9.Save the file -------------------------------------------------------------------------------
  write.csv(final_data_with_pop, paste0("<<<< FILEPATH REDACTED >>>>"))

}


equal_coverage_split(country, iso3, user, treat_input, lookup_input)
