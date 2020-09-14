###### Drop microdata that are outside the identified start and end ages, by nid
source(paste0(repo, "/collapse_functions/crosswalk/age_crosswalk.R"))

#### Initial setup
temp_lri_data <- lri_data
age_crosswalk_setup()
tracking <- get_tracking_sheet(sheet_name="LRI Vetting Sheet")
surveys_non_standard <- unique(tracking)

temp_lri_data <- merge(temp_lri_data, surveys_non_standard, by.x="nid", by.y="nid", all.x=TRUE, allow.cartesian = T)
temp_lri_data <- temp_lri_data[is.na(Youngest_age) | (!is.na(Youngest_age) & age_year >= Youngest_age)]
temp_lri_data <- temp_lri_data[is.na(Oldest_age) | (!is.na(Oldest_age) & age_year < Oldest_age)]
temp_lri_data[, Youngest_age := NULL]
temp_lri_data[, Oldest_age := NULL]

lri_data <- temp_lri_data
