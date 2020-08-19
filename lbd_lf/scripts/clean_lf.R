### Clean LF data
### Purpose: Lymphatic Filariasis (LF) data downloaded from the online ESPEN data portal needs to be combined
###          with LF data from other sources (systematic review, WHO, GAHI). This involves:
###               - Standardizing fields
###               - Checking for, and removing duplicates
###               - Applying necessary data adjustments (age and diagnostic)
#############################################################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Script Parameters ######################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

produce_crosswalk_data_set <- FALSE ### produce crosswalk data set? If false, the final modeling data set will be produced instead
perform_crosswalk_new <- TRUE ### turn on or off age/dx crosswalk code

save_outputs <- TRUE ### save cleaned data files?
drop_outliered_data <- TRUE ### drop outliered data?
output_crosswalk_quantiles <- FALSE ### save data files with had_lf_poly set at lower and upper 95% of crosswalked values?
save_by_region <- TRUE ### save separate data files by region? (prevents overwriting during resampling)

save_nonobvious_dups <- FALSE ### should identified nonobvious duplicates be written to disk? FALSE after initial run
nonobvious_dup_filepath <- <<<< FILEPATH REDACTED >>>>

diag <- "ICT"

pathadd <- <<<< FILEPATH REDACTED >>>>
shapefile_version <- "2020_02_20"

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Script Prep #############################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Set repo locations
## Get current user name
user <- Sys.info()[["user"]] ## Get current user name

## Set repo location and indicator group
core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

library(stringi)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Load data ##############################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load non-espen data (processed and combined in cleaning_WHO_data.R) and espen data
df <- fread(<<<< FILEPATH REDACTED >>>>)
espen_old <- fread(<<<< FILEPATH REDACTED >>>>)
espen <- fread(<<<< FILEPATH REDACTED >>>>)

df <- df[is.na(uniqueID)] # Drop rows from crosswalk data set

espen_old$Latitude <- round(espen_old$Latitude, 6)
espen_old$Longitude <- round(espen_old$Longitude, 6)
espen$Latitude <- round(espen$Latitude, 6)
espen$Longitude <- round(espen$Longitude, 6)

espen_old$temp_id <- 1:nrow(espen_old)
espen$temp_id <- 1:nrow(espen)

### Create Master_UIDs for ESPEN
# Recreate older Master_UIDs to retain later cleaning steps
espen_old$Master_UID <- paste0("temp_", c(1:nrow(espen_old)))

espen_merged2 <- merge(espen_old[, c("Country", "SurveyYear", "Longitude", "Latitude", "Examined", "Positive", "Master_UID")], espen, by.x = c("Country", "SurveyYear", "Longitude", "Latitude", "Examined", "Positive"), by.y = c("Country", "SurveyYear", "Longitude", "Latitude", "Examined", "Positive"), all = TRUE)
Master_UID_link_table <- espen_merged2[!is.na(Master_UID) & !is.na(ID), c("ID", "Master_UID")]
Master_UID_link_table <- Master_UID_link_table[, .SD[1], ID]

espen_merged <- merge(espen_old, espen, by.x = c("Country", "SurveyYear", "Longitude", "Latitude", "Examined", "Positive"), by.y = c("Country", "SurveyYear", "Longitude", "Latitude", "Examined", "Positive"), all = TRUE)
espen_dropped <- espen_merged[!is.na(temp_id.x) & is.na(temp_id.y)]

espen_working <- rbind(espen_old[temp_id %in% espen_dropped$temp_id.x], espen, fill = TRUE)
espen_working <- merge(espen_working, Master_UID_link_table, by = "ID", all.x = TRUE)
espen_working[!is.na(Master_UID.x), Master_UID := Master_UID.x]
espen_working[is.na(Master_UID) & !is.na(Master_UID.y), Master_UID := Master_UID.y]

# Create Master_UIDs for new ESPEN data
espen_working[is.na(Master_UID), Master_UID := paste0("espen_ID_", ID)]

espen <- copy(espen_working)

### Reconcile country field for ESPEN data
lookup_table <- fread(<<<< FILEPATH REDACTED >>>>)
espen_working <- merge(espen[is.na(ISO3)], lookup_table[, c("location_name", "iso3")], by.x = "Country", by.y = "location_name", all.x = TRUE)
espen_working[Country == "Congo (Kinshasa)", iso3 := "COD"]
espen_working[Country == "Sao Tome & Principe", iso3 := "STP"]
espen_working[Country == "Tanzania (Mainland)", iso3 := "TZA"]
espen_working$ISO3 <- espen_working$iso3
espen_working[, iso3 := NULL]
espen <- rbind(espen[!is.na(ISO3)], espen_working)

setnames(espen, c("ISO3", "Longitude", "Latitude", "Method_0", "Examined", "Positive", "Prevalence", "Method_2"), c("iso3", "longitude", "latitude", "data_collect_method", "N", "had_lf_poly", "lf_prev", "diagnostic"))
espen <- espen[!is.na(Year_start) | !is.na(Year_end) | !is.na(SurveyYear), ]
espen[, original_year := SurveyYear]
espen[, year := SurveyYear]
espen[, shapefile := NA]
espen[, location_code := NA]
espen[, point := 1]
espen[, source := "ESPEN"]
espen[, nid := NA]
espen$latitude <- as.double(espen$latitude)
espen$longitude <- as.double(espen$longitude)
espen$latitude <- round(espen$latitude, 5)
espen$longitude <- round(espen$longitude, 5)

# drop espen data with missing prevalence information
espen <- espen[!is.na(lf_prev) & !is.na(had_lf_poly) & !is.na(N)]

### Identify and fix data with incorrectly assigned country var ############################################
espen <- espen[!(is.na(latitude) | is.na(longitude))]

espen[iso3 == "MDG" & latitude > 0, latitude := -latitude] # Fix incorrect latitudes
setnames(espen, "iso3", "country")

setnames(df, "note_modler", "note_modeler")

df[is.na(point) == T] %>% nrow()
df[is.na(point) & !is.na(shapefile), point := 0]
df <- df[is.na(point) == F]

df[data_source == "" & extractor == "donkers", data_source := "sys_rev"]
df[data_source == "" & nid == "293159", data_source := "gahi"]
df[data_source == "" & nid == "136419", data_source := "gahi"]
df[nid == "136489", data_source := "sys_rev"]
df[nid == "222103", data_source := "sys_rev"]

df$data_collect_method <- tolower(df$data_collect_method)
table(df$data_collect_method)
df[grep("check", df$data_collect_method), data_collect_method := "SS/Spot"]
df[grep("site", df$data_collect_method), data_collect_method := "SS/Spot"]
df[unique(c(grep("lqas", df$data_collect_method), grep("mapping", df$data_collect_method))), data_collect_method := "Mapping"]
df[grep("tas", df$data_collect_method), data_collect_method := "TAS"]
df[data_collect_method %in% c("mapping", "baseline"), data_collect_method := "Mapping"]
df[data_collect_method %in% c("ss", "sentinal site", "ss/spot", "ss/sc", "sc"), data_collect_method := "SS/Spot"]
df[data_collect_method %in% c("c survey", "village surveillance", "other", "cl_tas", "operational research", "surveillance", "census", "post-endemic", "special survey", "post survey", "research"), data_collect_method := "Other"]
df[data_collect_method %in% c(
  "", "5-year programme to control bancroftian filariasis",
  "analyze 2007-2008 data to also conform to the 2011 antigen guidelines for the 6 to 7 yr old age group"
), data_collect_method := NA]
table(df$data_collect_method)

# clean/impute year_survey
df[is.na(year_survey) == T & year_end == 0 & is.na(year_start) == F & year_start != 0, year_survey := year_start]
df[is.na(year_survey) == T & is.na(year_end) == F & is.na(year_start) == F & year_start != 0 & year_end != 0, year_survey := round((year_start + year_end) / 2)]
df[is.na(year_survey) == T & is.na(year_start) == F & year_start != 0, year_survey := year_start]

df[(is.na(year_survey) == T | year_survey == 0 | year_survey == "") & reference_year != "" & is.na(reference_year) == F, year_survey := as.integer(reference_year) - 2]

### number being dropped due to missing year
table(is.na(df$year_survey) | df$year_survey == 0 | df$year_survey < 0)
dropped <- df[is.na(year_survey) == T | year_survey == 0] %>% .[, drop_reason := "year unspecified"]
df <- df[is.na(year_survey) == F & year_survey > 0]

## drop if WHO_DUP or EPI_REF_DUP
df[(WHO_DUP %in% c(1, 2)) | (EPI_REF_DUP %in% c(1)), ] %>% nrow()
dropped <- rbind(dropped, df[(WHO_DUP %in% c(1, 2)) | (EPI_REF_DUP %in% c(1)), ] %>% .[, drop_reason := "duplicated"])
df <- df[!(WHO_DUP %in% c(1, 2)) & !(EPI_REF_DUP %in% c(1)), ]

# drop out of range coordinates
df$latitude <- as.numeric(df$latitude)
df$longitude <- as.numeric(df$longitude)
df[(longitude > 180 | longitude < -180 | latitude > 180 | latitude < -180)] %>% nrow()
dropped <- rbind(dropped, df[point == 1 & is.na(longitude) == F & is.na(latitude) == F &
                                    (longitude > 180 | longitude < -180 | latitude > 180 | latitude < -180), ] %>% .[, ihme_loc_id := NULL] %>%
                   .[, drop_reason := "coordinates out of bounds"], fill = T)
df <- df[!(point == 1 & is.na(longitude) == F & is.na(latitude) == F & (longitude > 180 | longitude < -180 | latitude > 180 | latitude < -180))]

# standardize data source field
df[grep("gahi", df$data_source), data_source := "GAHI"]
df[data_source == "sys_rev", data_source := "sys_rev"]
df[grep("WHO", df$data_source), data_source := "WHO"]
df[grep("Dossier", df$data_source), data_source := "WHO"]
df[data_source %in% c("Ben_TZ", "Brazil_MOH"), data_source := "country"]

# geo_review field
df[, geo_review := 1]

# sort out prevalences
df$pop_mf <- as.numeric(gsub(",", "", df$pop_mf)) %>% round()
df$pop_ict <- as.numeric(gsub(",", "", df$pop_ict)) %>% round()
df$np_ict <- as.numeric(gsub(",", "", df$np_ict)) %>% round()
df$np_mf <- as.numeric(gsub(",", "", df$np_mf)) %>% round()
df$prev_mf <- as.double(gsub(",", ".", df$prev_mf))
df$prev_ict <- as.double(gsub(",", ".", df$prev_ict))
df$pop_test_strip <- as.numeric(df$pop_test_strip)

df[is.na(prev_ict) == T & pop_ict != 0 & is.na(np_ict) == F, prev_ict := np_ict / pop_ict]
df[is.na(prev_mf) == T & pop_mf != 0 & is.na(np_mf) == F, prev_mf := np_mf / pop_mf]
df[is.na(prev_test_strip) == T & pop_test_strip != 0 & is.na(np_test_strip) == F, prev_test_strip := np_test_strip / pop_test_strip]

df[is.na(np_ict) == T & is.na(prev_ict) == F & prev_ict < 1 & is.na(pop_ict) == F, np_ict := prev_ict * pop_ict]
df[is.na(np_mf) == T & is.na(prev_mf) == F & prev_mf < 1 & is.na(pop_mf) == F, np_mf := prev_mf * pop_mf]
df[is.na(np_test_strip) == T & is.na(prev_test_strip) == F & prev_test_strip < 1 & is.na(pop_test_strip) == F, np_test_strip := prev_test_strip * pop_test_strip]

df[is.na(np_ict) == T & is.na(prev_ict) == F & prev_ict >= 1 & is.na(pop_ict) == F, np_ict := (prev_ict / 100) * pop_ict]
df[is.na(np_mf) == T & is.na(prev_mf) == F & prev_mf >= 1 & is.na(pop_mf) == F, np_mf := (prev_mf / 100) * pop_mf]
df[is.na(np_test_strip) == T & is.na(prev_test_strip) == F & prev_test_strip >= 1 & is.na(pop_test_strip) == F, np_test_strip := (prev_test_strip / 100) * pop_test_strip]

## generate analysis fields using hierarchy for diagnosis preference (ICT > MF > Test Strip > BR > ELISA > PCR)
df[
  is.na(prev_ict) == F & pop_ict != 0 & is.na(pop_ict) == F & pop_ict != "" &
    ((pop_ict > 20 & pop_ict > pop_mf) | (pop_mf != 0 & is.na(pop_mf) == F & pop_mf / pop_ict < 2) | is.na(pop_mf) == T),
  c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_ict, np_ict, prev_ict, "ict")
  ]
df[is.na(lf_prev) == T & pop_mf != 0 & is.na(pop_mf) == F & pop_mf != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_mf, np_mf, prev_mf, "mf")]
df[is.na(lf_prev) == T & pop_test_strip != 0 & is.na(pop_test_strip) == F & pop_test_strip != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_test_strip, np_test_strip, prev_test_strip, "test_strip")]
df[is.na(lf_prev) == T & pop_br != 0 & is.na(pop_br) == F & pop_br != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_br, np_br, prev_br, "br")]
df[is.na(lf_prev) == T & pop_elisa != 0 & is.na(pop_elisa) == F & pop_elisa != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_elisa, np_elisa, prev_elisa, "elisa")]
df[is.na(lf_prev) == T & pop_pcr != 0 & is.na(pop_pcr) == F & pop_pcr != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_pcr, np_pcr, prev_pcr, "pcr")]

dropped <- rbind(dropped, df[is.na(N) == T, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "no valid sample size"], fill = T)

df <- df[is.na(N) == F]
dropped <- rbind(dropped, df[N < 20, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "N < 20"], fill = T)
df <- df[N > 20]
df$lf_pos <- df$lf_pos %>% as.integer()
dropped <- rbind(dropped, df[is.na(lf_pos) == T & is.na(lf_prev) == T, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "No numerator/prevalence"], fill = T)
df <- df[!(is.na(lf_pos) == T & is.na(lf_prev) == T)]

dropped <- rbind(dropped, df[lf_pos > N, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Illogical sample size"], fill = T)
df <- df[N >= lf_pos, ]
df[, lf_prev := lf_pos / N]

# change one polygon to point (known issue, check this if running again)
# fixing other discovered issues
df[point == 0 & (is.na(shapefile) == F | shapefile == "") & (location_code == "" | is.na(location_code) == T), point := 1]
df[location_code == "BINAH", location_code := 1]
df[shapefile == "lf_g2015_2001_0", shapefile := "lf_g2015_2010_0"]
df[shapefile == "lf_g2015_2010_0" & adm0 %in% c("togo", "vietnam") & adm1 != "", shapefile := "lf_g2015_2011_2"]
dropped <- rbind(dropped, df[grep(123.6693, df$latitude)] %>% .[, drop_reason := "coordinates out of bounds"], fill = T)
df <- df[!(grep(123.6693, df$latitude))]

### drop non-georeferenced data ###############################################
dropped <- rbind(dropped, df[point == 0 & (is.na(shapefile) == T | shapefile == "")] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Polygon not geolocated"], fill = T)
df <- df[!(point == 0 & (is.na(shapefile) == T | shapefile == ""))]

dropped <- rbind(dropped, df[point == 1 & (is.na(latitude) == T | is.na(longitude) == T | latitude == "" | longitude == "")] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Point not geolocated"], fill = T)
df <- df[!(point == 1 & (is.na(latitude) == T | is.na(longitude) == T | latitude == "" | longitude == ""))]

# drop outliered data
dropped <- rbind(dropped, df[outlier == 1] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Outliered"], fill = T)
df <- df[is.na(outlier) == T]

# tag test polygons
df[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Kenema", "Kailahun", "Kono"), location_code := 2654]
df[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Bombali", "Koinadugu", "Port Loko", "Tonkolili", "Kambia"), location_code := 2655]
df[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Pujehun", "Bonthe", "Moyamba", "Bo"), location_code := 2656]
df[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Western Rural", "Freetown"), location_code := 2657]
df[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & location_code == 2655, resamp_test_case := 3]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Step 5: Final subsetting and export ################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# currently imputing population-based survey design for unknowns (since these would generally be TAS)
df[sampling_prod == 1 | sampling_prod == 9, sampling := 1]
df[sampling_prod == 0, sampling := 0]
df[point == 0 & is.na(sampling_prod) == T, sampling := 1]
df$sampling_prod <- NULL

### Harmonize country fields
lookup_table <- fread(<<<< FILEPATH REDACTED >>>>)
df <- merge(df, lookup_table[, c("location_name", "iso3")], by.x = "adm0", by.y = "location_name", all.x = TRUE)
df[(is.na(ihme_loc_id) | (ihme_loc_id == "")) & !is.na(iso3), ihme_loc_id := iso3]

df[nid == "136477", ihme_loc_id := "IND"]
df[nid == "288829", ihme_loc_id := "IDN"]
df[nid == "136419", ihme_loc_id := "GHA"]
df[nid == "222121", ihme_loc_id := "BRA"]
df[adm0 == "CAR", ihme_loc_id := "CAR"]
df[adm0 == "Congo", ihme_loc_id := "COG"]
df[adm0 == "Gambia", ihme_loc_id := "GMB"]
df[adm0 == "Micronesia (Federated States of)", ihme_loc_id := "FSM"]
df[adm0 == "Nauru", ihme_loc_id := "NRU"]
df[adm0 == "Sao Tome Principe", ihme_loc_id := "STP"]
df[adm0 == "Wallis and Futuna", ihme_loc_id := "WLF"]
df[location_name == "Papua New Guinea", ihme_loc_id := "PNG"]
df[adm0 == "Mayotte", ihme_loc_id := "MYT"]

df$country <- df$ihme_loc_id

## Calculate sample sizes
df[is.na(N) & !is.na(lf_pos) & !is.na(mean), N := round(lf_pos / mean, 0)]
df <- df[!is.na(N) & !is.nan(N) & !(N == "")] # drop rows without sample size

#  Drop India TAS data in df
df <- df[!(country == "IND" & data_collect_method == "TAS")]
india_tas <- fread(<<<< FILEPATH REDACTED >>>>)
setnames(india_tas, "adm0", "country")
setnames(india_tas, "year_survey", "year")
india_tas[, V1 := NULL]
india_tas[, original_year := year]
india_tas$country <- "IND"

setnames(df, "year_survey", "year")
setnames(df, "lf_pos", "had_lf_poly")

df <- rbindlist(list(df, india_tas), fill = TRUE)

df$weight <- 1
df$cluster_id <- 1
df$sampling <- 1

# load master shape
if (!("master_shape_all" %in% ls())) {
  master_shape_all <- readRDS(<<<< FILEPATH REDACTED >>>>)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Standardize fields #####################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make coordinates uniform
df$latitude <- round(df$latitude, 5)
df$longitude <- round(df$longitude, 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Manually update nids missing from GAHI data #######################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df[nid == "Eigege 2017", nid := "378953"]
df[nid == "", nid := NA]
df$nid <- as.integer(df$nid)

setnames(df, "data_source", "source")

# Update GAHI rows with missing nid, when possible
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(1028), nid := 147755]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(1549, 1550, 1564), nid := 136492]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(1589, 1590), nid := 147743]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(7101, 7102, 7638, 7639, 7640, 7642, 7731, 7781, 7861), nid := 136409]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(10901, 10902, 10903, 10904, 10905, 10906, 10907, 10908, 10909, 10910, 10911, 10912), nid := 147687]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(10854, 10853, 10852, 10848, 10849, 10850, 10871, 10872, 10870), nid := 147859]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(10851, 10855), nid := 147811]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(10897, 10898, 10899, 10900), nid := 147837]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(10919), nid := 148400]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(10873), nid := 147819]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(676), nid := 136521]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(745), nid := 136552]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(747), nid := 147783]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9065, 9066, 9067, 9070, 9079, 9080, 9085, 9086, 9088), nid := 147799]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9194), nid := 147696]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9188, 9189), nid := 136508]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(8900), nid := 147855]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9143), nid := 147705]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9114, 9116, 9121, 9130), nid := 136516]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(8898), nid := 136587]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9046, 9054), nid := 147830]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(8861, 8862, 8864), nid := 147865]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(8886), nid := 147849]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9102, 9104, 9120, 9125, 9127, 9129, 9131, 9132, 9133, 9137, 9144), nid := 147833]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9083), nid := 147789]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9162), nid := 137360]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(8867), nid := 136583]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(8899, 8901), nid := 136458]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9099, 9101, 9105, 9107, 9108, 9109, 9110, 9111, 9112, 9115, 9117, 9118, 9119, 9128, 9135), nid := 136506]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9056), nid := 147839]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(9002, 9023, 9044), nid := 136447]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(11353, 11356, 11357, 11358, 11359, 11361, 11363, 11368, 11369, 11370, 11372, 11375), nid := 147790]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(11350, 11373), nid := 136547]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(11414, 11415), nid := 147851]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(11362, 11380), nid := 136463]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(11535, 11536, 11538, 11540, 11541, 11543, 11545), nid := 137364]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(13693), nid := 136466]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(13711), nid := 147726]
df[source == "GAHI" & is.na(nid) & Master_UID %in% c(13717), nid := 147681]

df[Master_UID %in% c(19035, 19036, 19037, 19038, 19039, 19040, 19041, 19042, 19043, 19044, 19045, 19046, 19047, 19048, 19049), nid := 414527]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Split out data by geographies ##########################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This is to chunk out the data cleaning process some of which is necessarily manual

### Africa ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# subset df to africa data. Clean Africa geography first. Re-combine data from other geographies later.
africa_iso <- c(
  "DZA", "AGO", "SHN", "BEN", "BWA",
  "BFA", "BDI", "CMR", "CPV", "CAF",
  "TCD", "COM", "COG", "COD", "DJI",
  "EGY", "GNQ", "ERI", "ETH", "GAB",
  "GMB", "GHA", "GIN", "GNB", "CIV",
  "KEN", "LSO", "LBR", "LBY", "MDG",
  "MWI", "MLI", "MRT", "MUS", "MYT",
  "MAR", "MOZ", "NAM", "NER", "NGA",
  "STP", "REU", "RWA", "STP", "SEN",
  "SYC", "SLE", "SOM", "ZAF", "SSD",
  "SHN", "SDN", "SWZ", "TZA", "TGO",
  "TUN", "UGA", "COD", "ZMB", "TZA",
  "ZWE"
)

# split out non-africa data
df_outside_afr <- df[!(country %in% africa_iso)]
df <- df[country %in% africa_iso]
df_polys <- df[point == 0]
df <- df[point == 1 & is.na(latitude) == F & is.na(longitude) == F]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Bulk cleaning of identified data quality issues ############################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Ethiopia 2013 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: collapse split out same loc + year surveys into multiple rows
# Fix: Sum up unique locations

eth_locs <- df[country == "ETH" & year == 2013, .N, by = .(latitude, longitude)] %>% .[N > 1]
eth_locs$N <- NULL
eth_locs$had_lf_poly <- apply(eth_locs, 1, FUN = function(x) sum(df[latitude == x[1] & longitude == x[2], had_lf_poly]))
eth_locs$N <- apply(eth_locs, 1, FUN = function(x) sum(df[latitude == x[1] & longitude == x[2], N]))
eth_locs$Master_UID <- apply(eth_locs, 1, FUN = function(x) as.numeric(paste0(df[latitude == x[1] & longitude == x[2], Master_UID])))
eth_locs[, lf_prev := had_lf_poly / N]

eth_locs <- cbind(eth_locs, df[country == "ETH" & year == 2013, names(df)[!(names(df) %in% names(eth_locs))], with = F][1:nrow(eth_locs)])
df <- df[!(latitude %in% eth_locs$latitude & country == "ETH" & year == 2013)]
df <- rbind(df, eth_locs)

### Gabon 2014 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Duplication between non-ESPEN and ESPEN datasets
# Fix: Keep ESPEN rows

gabon_drop <- df[country == "GAB" & year == 2014, Master_UID]
df <- df[!(Master_UID %in% gabon_drop)]

### Togo mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Duplication in ESPEN dataset
# Fix: De-duplication (N = 7 duplicates)

tgo_check <- espen[Country == "Togo" & data_collect_method == "Mapping survey", ]
tgo_ulocs <- unique(tgo_check$LocationName) %>% as.data.table()
tgo_ulocs$uloc_id <- c(1:nrow(tgo_ulocs))
setnames(tgo_ulocs, ".", "LocationName")
tgo_check <- merge(tgo_check, tgo_ulocs, by = "LocationName")
tgo_dups <- table(tgo_check$uloc_id) %>%
  as.data.table() %>%
  .[N > 1]
drop_tgo_uids <- c()
for (i in tgo_dups$V1) {
  drop_tgo_uids <- c(drop_tgo_uids, tgo_check[uloc_id == i, Master_UID][2:length(tgo_check[uloc_id == i, Master_UID])])
}
espen <- espen[!(Master_UID %in% drop_tgo_uids), ]


### ESPEN outliers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Some country survey series have been identified as outliers
# Fix: drop outliers according to script parameter

if (drop_outliered_data) {
  drop_espen_countries <- c("AGO", "MRT") # (N = 340)
} else {
  drop_espen_countries <- NULL
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Automated & manual de-duplication: Step 1 ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 1: Split out by type of duplication:
# Obvious duplication = same location, year, and sample size; these can be de-duplicated easily by keeping only one copy of each duplicate set
# Nonobvious duplication = generally the same location, year, and sample size; requires manual evaluation. Could be result of multiple reporting mechanisms & data processing artifacts

# split_dups() in data_cleaning_functions.R

### ESPEN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Split obvious and nonobvious duplicates
espen_check2 <- unique(espen[, .(latitude, longitude, year, country)])
espen_check2$espen.loc.yr.id <- c(1:nrow(espen_check2))
espen <- merge(espen, espen_check2, by = c("latitude", "longitude", "year", "country"))
espen_loc.yr.tab <- espen$espen.loc.yr.id %>%
  table() %>%
  as.data.table()
setnames(espen_loc.yr.tab, ".", "espen.loc.yr.id")
espen_check_dups <- espen_loc.yr.tab[N > 1]
espen_check_dups <- espen[espen.loc.yr.id %in% espen_check_dups$espen.loc.yr.id]

espen_split_dups <- split_dups(data = espen_check_dups, group_id = "espen.loc.yr.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = c("N", "had_lf_poly", "lf_prev"), get_unique = c("diagnostic", "LocationName", "N", "had_lf_poly", "lf_prev", "data_collect_method", "Age_start", "Age_end", "Sex"))

espen_obvi_drop <- espen_split_dups[[1]]
espen_nonobvi <- espen_split_dups[[2]]
espen_nonobvi_summary <- espen_split_dups[[3]]

# Resolve nonobvious duplication due to different reported/retained diagnostic for the same survey
espen_diff_diagnostics <- espen_nonobvi_summary[grep(";", espen_nonobvi_summary$unique_diagnostic)] # could output this and add to diagnostic crosswalk
espen_diff_diagnostics[grep("ICT", espen_diff_diagnostics$unique_diagnostic), ict := 1]
espen_diff_diagnostics[grep("Blood smear", espen_diff_diagnostics$unique_diagnostic), mf := 1]
espen_diff_diagnostics[grep("rapid", espen_diff_diagnostics$unique_diagnostic), rapid_test := 1]
espen_diff_diagnostics[grep("Filtration", espen_diff_diagnostics$unique_diagnostic), filtration := 1]

# Remove obvious duplicates and non-obvious duplicates due to different reported/retained diagnostics
espen_obvi_drop <- c(espen_obvi_drop, espen[espen.loc.yr.id %in% espen_diff_diagnostics[ict == 1, group_id] & diagnostic != "ICT", Master_UID])
espen_obvi_drop <- c(espen_obvi_drop, espen[espen.loc.yr.id %in% espen_diff_diagnostics[is.na(ict) == T & mf == 1, group_id] & diagnostic != "Blood smear", Master_UID])
espen_obvi_drop <- c(espen_obvi_drop, espen[espen.loc.yr.id %in% espen_diff_diagnostics[is.na(ict) == T & is.na(mf) == T & rapid_test == 1, group_id] & diagnostic != "New rapid test", Master_UID])
espen_obvi_drop <- c(espen_obvi_drop, espen[espen.loc.yr.id %in% espen_diff_diagnostics[is.na(ict) == T & is.na(mf) == T & is.na(rapid_test) == T & filtration == 1, group_id] & diagnostic != "Filtration", Master_UID])
espen_nonobvi <- espen_nonobvi[!(espen_nonobvi %in% espen_obvi_drop)]
espen_nonobvi_summary <- espen_nonobvi_summary[!(group_id %in% espen_diff_diagnostics$group_id)]

### Non-ESPEN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Split obvious and nonobvious duplicates
check2 <- unique(df[point == 1, .(latitude, longitude, year, country)])
check2$loc.yr.id <- c(1:nrow(check2))
df <- merge(df, check2, by = c("latitude", "longitude", "year", "country"))
loc.yr.tab <- df$loc.yr.id %>%
  table() %>%
  as.data.table()
setnames(loc.yr.tab, ".", "loc.yr.id")
check_dups <- loc.yr.tab[N > 1]
check_dups <- df[loc.yr.id %in% check_dups$loc.yr.id]

df_split_dups <- split_dups(data = check_dups, group_id = "loc.yr.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = c("N", "had_lf_poly", "lf_prev"), get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "source"))
df_obvi_drop <- df_split_dups[[1]]
df_nonobvi <- df_split_dups[[2]]
df_nonobvi_summary <- df_split_dups[[3]]

# Resolve nonobvious duplication due to different reported/retained diagnostic for the same survey
df_diff_diagnostics <- df_nonobvi_summary[grep(",", df_nonobvi_summary$unique_diagnostic)] # could output this and add to diagnostic crosswalk
df_diff_diagnostics[grep("ict", df_diff_diagnostics$unique_diagnostic), ict := 1]
df_diff_diagnostics[grep("mf", df_diff_diagnostics$unique_diagnostic), mf := 1]
df_diff_diagnostics[grep("elisa", df_diff_diagnostics$unique_diagnostic), elisa := 1]

# Remove obvious duplicates and non-obvious duplicates due to different reported/retained diagnostics
df_obvi_drop <- c(df_obvi_drop, df[loc.yr.id %in% df_diff_diagnostics[ict == 1, group_id] & diagnostic != "ict", Master_UID])
df_obvi_drop <- c(df_obvi_drop, df[loc.yr.id %in% df_diff_diagnostics[is.na(ict) == T & mf == 1, group_id] & diagnostic != "mf", Master_UID])
df_obvi_drop <- c(df_obvi_drop, df[loc.yr.id %in% df_diff_diagnostics[is.na(ict) == T & is.na(mf) == T & elisa == 1, group_id] & diagnostic != "elisa", Master_UID])
df_nonobvi <- df_nonobvi[!(df_nonobvi %in% df_obvi_drop)]
df_nonobvi_summary <- df_nonobvi_summary[!(group_id %in% df_diff_diagnostics$group_id)]

if (save_nonobvious_dups) {
  write.csv(df_nonobvi_summary, paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
  write.csv(espen_nonobvi_summary, paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Automated & manual de-duplication: Step 2 ##############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 2: Resolve both types of duplication identified in step 1

### Obvious duplication ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Drop (df N = 138, espen N = 577)
df <- df[!(Master_UID %in% unlist(df_obvi_drop))]
espen <- espen[!(Master_UID %in% unlist(espen_obvi_drop))]
setnames(espen, c("Age_start", "Age_end"), c("age_start", "age_end"))

### Non-obvious duplication ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# These are manually reviewed.

# read in manual checks of nonobvious duplicates (487 total)
espen_nonobvi_manual <- fread(<<<< FILEPATH REDACTED >>>>)
df_nonobvi_manual <- fread(<<<< FILEPATH REDACTED >>>>)
espen_nonobvi_checked <- fread(<<<< FILEPATH REDACTED >>>>)
df_nonobvi_checked <- fread(<<<< FILEPATH REDACTED >>>>)

# extract rows that should be retained
no_drop <- c(espen_nonobvi_checked[drop == 0, Master_UID], df_nonobvi_checked[drop == 0, Master_UID])

# extract rows that are not clearly
unsure_uids <- c(espen_nonobvi_manual[ask_liz == 1 | ignore == 1, uids], df_nonobvi_manual[ask_liz == 1 | ignore == 1, uids])
unsure_uids <- gsub(";", ",", unsure_uids)
unsure_uids <- strsplit(unsure_uids, ",") %>% unlist()
unsure_uids <- unsure_uids[!(unsure_uids %in% no_drop)]
message(paste0("Dropping ", length(unsure_uids), " data rows due to unclear duplicates at location-year. Please resolve if able."))

drop_uids <- c(espen_nonobvi_manual[is.na(drop) == F, drop], df_nonobvi_manual[is.na(drop) == F, drop])
drop_uids <- gsub(";", ",", drop_uids)
drop_uids <- strsplit(drop_uids, ",") %>% unlist()
message(paste0("Dropping ", length(drop_uids), " data rows to resolve location-year duplicates."))

### Fixing some incorrect coordinates identified during de-duplication ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Fixing some coordinates per manual checks")
df_fix_coords <- df[Master_UID %in% df_nonobvi_checked[drop == 0 & replace_coords == 1, Master_UID], ]
df_fix_coords$latitude <- NULL
df_fix_coords$longitude <- NULL
new_coords <- df_nonobvi_checked[drop == 0 & replace_coords == 1, .(Master_UID, new_lat, new_long)]
setnames(new_coords, c("new_lat", "new_long"), c("latitude", "longitude"))
df_fix_coords$Master_UID <- as.integer(df_fix_coords$Master_UID)
df_fix_coords <- merge(df_fix_coords, new_coords, by = "Master_UID")
df <- rbind(df[!(Master_UID %in% df_fix_coords)], df_fix_coords)

espen_fix_coords <- espen[Master_UID %in% espen_nonobvi_checked[drop == 0 & replace_coords == 1, Master_UID], ]
espen_fix_coords$latitude <- NULL
espen_fix_coords$longitude <- NULL
new_coords <- espen_nonobvi_checked[drop == 0 & replace_coords == 1, .(Master_UID, new_lat, new_long)]
setnames(new_coords, c("new_lat", "new_long"), c("latitude", "longitude"))
espen_fix_coords <- merge(espen_fix_coords, new_coords, by = "Master_UID")
espen <- rbind(espen[!(Master_UID %in% espen_fix_coords)], espen_fix_coords, fill = TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Automated & manual de-duplication: Step 3 #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 3: Do the same as step 1 & 2 but between the ESPEN & non-ESPEN datasets

# Split obvious (N = 3571) and nonobvious duplicates
combo_ref <- merge(check2, espen_check2, by = c("latitude", "longitude", "year", "country"), all = T)
combo_ref$combo.loc.yr.id <- c(1:nrow(combo_ref))
combo_ref[is.na(loc.yr.id) == F & is.na(espen.loc.yr.id) == F, shared := 1]
shared_cols <- c("data_collect_method", "source", "latitude", "longitude", "year", "N", "had_lf_poly", "lf_prev", "Master_UID", "country", "diagnostic", "age_start", "age_end", "original_year", "point", "nid")
combo <- rbind(df[!(Master_UID %in% c(unsure_uids, drop_uids)), c(shared_cols, "age_group"), with = F], espen[!(Master_UID %in% c(unsure_uids, drop_uids)), shared_cols, with = F], fill = T)
combo <- merge(combo, combo_ref[, .(latitude, longitude, year, combo.loc.yr.id, shared, country)], by = c("latitude", "longitude", "year", "country"), all.x = T)
combo[, had_lf_poly := as.integer(had_lf_poly)]
combo[, N := as.integer(N)]
combo[, lf_prev := had_lf_poly / N]
combo[diagnostic == "Blood smear", diagnostic := "mf"]
combo[diagnostic == "ICT", diagnostic := "ict"]
combo[data_collect_method == "Mapping survey", data_collect_method := "Mapping"]
combo[data_collect_method == "Surveillance", data_collect_method := "SS/Spot"]

combo_check_dups <- combo[shared == 1, ]
if (nrow(combo_check_dups) > 0) {
  combo_split_dups <- split_dups(data = combo_check_dups, group_id = "combo.loc.yr.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = c("N", "had_lf_poly", "lf_prev"), get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "source"))
  combo_obvi_drop <- combo_split_dups[[1]]
  combo_nonobvi <- combo_split_dups[[2]]
  combo_nonobvi_summary <- combo_split_dups[[3]]
}

# Resolve nonobvious duplication due to different reported/retained diagnostic for the same survey
combo_diff_diagnostics <- combo_nonobvi_summary[grep(",", combo_nonobvi_summary$unique_diagnostic)] # could output this and add to diagnostic crosswalk
combo_diff_diagnostics[grep("ict", combo_diff_diagnostics$unique_diagnostic), ict := 1]
combo_diff_diagnostics[grep("mf", combo_diff_diagnostics$unique_diagnostic), mf := 1]
combo_diff_diagnostics[grep("elisa", combo_diff_diagnostics$unique_diagnostic), elisa := 1]

# Remove obvious duplicates and non-obvious duplicates due to different reported/retained diagnostics
combo_obvi_drop <- c(combo_obvi_drop, combo[combo.loc.yr.id %in% combo_diff_diagnostics[ict == 1, group_id] & diagnostic != "ict", Master_UID])
combo_obvi_drop <- c(combo_obvi_drop, combo[combo.loc.yr.id %in% combo_diff_diagnostics[is.na(ict) == T & mf == 1, group_id] & diagnostic != "mf", Master_UID])
combo_nonobvi <- combo_nonobvi[!(combo_nonobvi %in% combo_obvi_drop)]
combo_nonobvi_summary <- combo_nonobvi_summary[!(group_id %in% combo_diff_diagnostics$group_id)]

# read in and sort out nonobvious duplicates (checked manually) (108 total)
combo_nonobvi_manual <- fread(<<<< FILEPATH REDACTED >>>>)

unsure_combo <- combo_nonobvi_manual[ask_liz == 1 | ignore == 1, uids]
unsure_combo <- gsub(";", ",", unsure_combo)
unsure_combo <- strsplit(unsure_combo, ",") %>% unlist()
message(paste0("Dropping ", length(unsure_combo), " data rows due to unclear duplicates at location-year. Please resolve if able."))
drop_combo <- combo_nonobvi_manual[is.na(drop) == F, drop]
drop_combo <- gsub(";", ",", drop_combo)
drop_combo <- strsplit(drop_combo, ",") %>% unlist()
message(paste0("Dropping ", length(drop_combo), " data rows to resolve location-year duplicates."))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Automated & manual de-duplication: Step 4 #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 4: Apply drops accumulated from previous steps & combine ESPEN & non-ESPEN datasets

all_drop_uids <- c(unsure_uids, drop_uids, unsure_combo, drop_combo, df_split_dups[[1]], espen_split_dups[[1]], combo_split_dups[[1]])

potential_full <- rbind(df[!(Master_UID %in% all_drop_uids), c(shared_cols, "age_group"), with = F], espen[!(Master_UID %in% all_drop_uids) & !(country %in% drop_espen_countries), shared_cols, with = F], fill = T)
potential_full <- potential_full[!(Master_UID %in% espen_split_dups[[2]])]
potential_full <- potential_full[!(Master_UID %in% df_split_dups[[2]])]
potential_full <- potential_full[is.na(latitude) == F]
potential_full[data_collect_method == "Mapping survey", data_collect_method := "Mapping"]
potential_full[data_collect_method == "Surveillence" | data_collect_method == "Surveillance", data_collect_method := "SS/Spot"]
potential_full[, c("weight", "shapefile", "location_code", "resamp_test_case", "cluster_id", "sampling") := list(1, NA, NA, NA, 1, NA)]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Automated & manual de-duplication: Step 5 #############################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Step 5: Do the same for outside Africa geographies

c <- unique(df_outside_afr[, .(latitude, longitude, year, country, shapefile, location_code)])
c$loc.yr.id <- c(1:nrow(c))
df_outside_afr <- merge(df_outside_afr, c, by = c("latitude", "longitude", "year", "country", "shapefile", "location_code"))
loc.yr.tab <- df_outside_afr$loc.yr.id %>%
  table() %>%
  as.data.table()
setnames(loc.yr.tab, ".", "loc.yr.id")
check_dups <- loc.yr.tab[N > 1]
check_dups <- df_outside_afr[loc.yr.id %in% check_dups$loc.yr.id]

oafr_split_dups <- split_dups(data = check_dups, group_id = "loc.yr.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = c("N", "had_lf_poly", "lf_prev"), get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "source"))
oafr_obvi_drop <- oafr_split_dups[[1]]
oafr_nonobvi <- oafr_split_dups[[2]]
oafr_nonobvi_summary <- oafr_split_dups[[3]]

# drop obvi_drops
df_outside_afr <- df_outside_afr[!(Master_UID %in% oafr_obvi_drop)]

# drop manually excluded duplicates
manual_checks <- fread(<<<< FILEPATH REDACTED >>>>)
keeps <- manual_checks[is.na(keep_uid) == F & keep_uid != "", keep_uid]
keeps <- c(keeps, oafr_nonobvi_summary[country %in% c("BRA", "CHN", "FJI", "FSM", "GUY", "HTI", "PYF", "TUV", "VUT", "WSM", "YEM"), uids])
keeps <- c(keeps, unlist(strsplit(grep(",", keeps, value = T), ",")))
keeps <- c(keeps, unlist(strsplit(grep(";", keeps, value = T), ";")))
keeps <- keeps[!(keeps %in% grep(",", keeps, value = T))] %>% as.integer()
drops <- oafr_nonobvi[!(oafr_nonobvi %in% keeps)]
df_outside_afr <- df_outside_afr[!(Master_UID %in% drops), ]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### De-duplication and outliering ##########################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Manually identified duplicates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

temp1 <- fread(<<<< FILEPATH REDACTED >>>>)
temp2 <- fread(<<<< FILEPATH REDACTED >>>>)
temp3 <- fread(<<<< FILEPATH REDACTED >>>>)
temp4 <- fread(<<<< FILEPATH REDACTED >>>>)

drop_mids <- c(temp1[leaflet_drop == 1, Master_UID], temp2[leaflet_drop == 1, Master_UID], temp3[leaflet_drop == 1, Master_UID], temp4[leaflet_drop == 1 & country != "GHA", Master_UID])
potential_full <- potential_full[!(Master_UID %in% drop_mids)]
potential_full <- as.data.table(potential_full)

### Various outliers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (drop_outliered_data) {
  # Cameroon 2003 mapping data
  cmr_outlier <- potential_full[country == "CMR" & year == 2003 & data_collect_method == "Mapping" & source == "ESPEN", Master_UID]
  potential_full <- potential_full[!(Master_UID %in% cmr_outlier)]
  
  # CIV 2000-2001 mapping data
  civ_outlier <- potential_full[country == "CIV" & year %in% c(2000, 2001), Master_UID]
  potential_full <- potential_full[!(Master_UID %in% civ_outlier)]
  
  # data points in northern Mali, Niger and Chad
  north_MNC <- potential_full[country %in% c("MLI", "NER", "TCD") & (latitude > 15.9), Master_UID]
  potential_full <- potential_full[!(Master_UID %in% north_MNC)]
  
  # large SS Egypt datapoint
  potential_full <- potential_full[!(Master_UID %in% 4782)]
}


### Nigeria duplicates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Duplicates surveys different by year
# Fix: De-duplication

nga <- potential_full[country == "NGA"]
nga_locs <- nga[, .(latitude, longitude, N, had_lf_poly)] %>% unique()
nga_locs$nga_id <- c(1:nrow(nga_locs))
nga <- merge(nga, nga_locs, by = c("latitude", "longitude", "N", "had_lf_poly"))
nga_dup_ids <- table(nga$nga_id) %>%
  as.data.table() %>%
  .[N > 1, V1]
nga <- nga[nga_id %in% nga_dup_ids]
nga_drop <- c()

for (id in unique(nga$nga_id)) {
  print(id)
  y <- nga[nga_id == id, original_year]
  if (!is.na(y)) {
    if (max(y) - min(y) <= 2) {
      nga_drop <- c(nga_drop, grep("temp", nga[nga_id == id, Master_UID], value = T))
    } else {
      if (max(y) == 2008) nga_drop <- c(nga_drop, nga[nga_id == id & year == 2008, Master_UID])
    }
  }
}
nga_drop <- unlist(nga_drop)
nga_drop <- c(nga_drop, "14316", "14296", "14080", "14291", "14286") # (N = 108)
potential_full <- potential_full[!(Master_UID %in% nga_drop)]

### Ghana duplicates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Duplicates surveys different by diagnostics
# Fix: De-duplication

gha <- potential_full[country == "GHA" & original_year <= 2000]
gha$latitude <- round(gha$latitude, 5)
gha$longitude <- round(gha$longitude, 5)
gha_locs <- gha[, .(latitude, longitude, N)] %>% unique()
gha_locs$gha_id <- c(1:nrow(gha_locs))
gha <- merge(gha, gha_locs, by = c("latitude", "longitude", "N"))
gha_dup_ids <- table(gha$gha_id) %>%
  as.data.table() %>%
  .[N > 1, V1]
gha <- gha[gha_id %in% gha_dup_ids]
gha_drop <- c()

for (id in unique(gha$gha_id)) {
  gha_drop <- c(gha_drop, gha[gha_id == id & source == "GAHI", Master_UID])
}
gha_drop <- unlist(gha_drop)
potential_full <- potential_full[!(Master_UID %in% gha_drop)]

### Burkina Faso duplicates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Duplicates surveys different by year
# Fix: De-duplication

bfa <- potential_full[country == "BFA" & year %in% c(2000, 2001)]
bfa$latitude <- round(bfa$latitude, 5)
bfa$longitude <- round(bfa$longitude, 5)
bfa_locs <- bfa[, .(latitude, longitude, N)] %>% unique()
bfa_locs$bfa_id <- c(1:nrow(bfa_locs))
bfa <- merge(bfa, bfa_locs, by = c("latitude", "longitude", "N"))
bfa_dup_ids <- table(bfa$bfa_id) %>%
  as.data.table() %>%
  .[N > 1, V1]
bfa <- bfa[bfa_id %in% bfa_dup_ids]
bfa_drop <- c()

for (id in unique(bfa$bfa_id)) {
  ids <- bfa[bfa_id == id, Master_UID]
  diag_pref <- bfa[bfa_id == id & diagnostic %in% c("ict", "ICT"), Master_UID]
  keep <- bfa[Master_UID %in% diag_pref & source != "WHO", Master_UID][1]
  
  bfa_drop <- c(bfa_drop, ids[!(ids %in% keep)])
}
bfa_drop <- unlist(bfa_drop)
potential_full <- potential_full[!(Master_UID %in% bfa_drop)]

### DRC duplicates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Duplicates surveys across datasets & years (GAHI 2010, ESPEN 2011, WHO 2015)
# Fix: De-duplication. ESPEN preferred
drc <- potential_full[country == "COD" & year %in% c(2010, 2011, 2015)]
drc$latitude <- round(drc$latitude, 5)
drc$longitude <- round(drc$longitude, 5)
drc_locs <- drc[, .(latitude, longitude, N)] %>% unique()
drc_locs$drc_id <- c(1:nrow(drc_locs))
drc <- merge(drc, drc_locs, by = c("latitude", "longitude", "N"))
drc_dup_ids <- table(drc$drc_id) %>%
  as.data.table() %>%
  .[N > 1, V1]
drc <- drc[drc_id %in% drc_dup_ids]
drc_drop <- c()

for (id in unique(drc$drc_id)) {
  ids <- drc[drc_id == id, Master_UID]
  if (nrow(drc[drc_id == id & source == "GAHI"]) > 0) {
    keep <- drc[drc_id == id & source == "GAHI", Master_UID][1]
  } else {
    keep <- drc[drc_id == id & source == "ESPEN", Master_UID][1]
  }
  
  drc_drop <- c(drc_drop, ids[!(ids %in% keep)])
}
drc_drop <- unlist(drc_drop) # (N = 120)
potential_full <- potential_full[!(Master_UID %in% drc_drop)]

### India nonrepresentative data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: Extracted data is not representative of the polygon unit referenced
# Fix: Drop
nonrep_nids <- c(136407, 136407, 136477, 136505, 136540, 28840, 28830, 288862, 288865, 288920, 288923, 136462, 136503)
nonrep_muids <- c(9498, 9343, 9330, 9146, 9344, 9637)

### Nepal duplicate TAS polygons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: same numerator & denominator & shapefile/location_code
# Fix: Drop duplicates
dup_tas_muids <- c(13904, 13906, 13908, 13910, 13912, 13915, 23924, 23926, 23927, 23928, 23929)

### Myanmar duplicate TAS polygons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: same numerator & denominator & shapefile/location_code
# Fix: Drop duplicates
dup_tas_muids <- c(dup_tas_muids, 21473, 21480)

### Bangladesh duplicate TAS polygons ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: same numerator & denominator & shapefile/location_code
# Fix: Drop duplicates
dup_tas_muids <- c(dup_tas_muids, 108, 109, 164, 104, 241, 117, 116, 114, 105, 112, 118, 113, 149, 154, 155)

### Thailand problematic polygons  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Issue: same numerator & denominator & shapefile/location_code
# Fix: Drop duplicates
dup_tas_muids <- c(dup_tas_muids, 17650, 23960, 17647, 23958, 23957, 23956, 23955)

# Issue: small N, unable to reference original source, pre-2000
# Fix: Drop
nonrep_nids <- c(nonrep_nids, 285490, 147658, 147884, 285491, 136490)

### TAS manually identified duplicates (N = 25) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fix: Drop duplicates
dup_tas_muids <- c(
  dup_tas_muids, 21959, 21896, 21899, 23751, 23755, 23669, 23625, 23784, 11137, 23628, 23681, 23616, 23762, 23495, 23429, 23445, 10760, 23351, 10744,
  23459, 23460, 23461, 23398, 23399, 23400
)

### Apply Asia drops ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df_outside_afr <- df_outside_afr[!(nid %in% c(nonrep_nids)) & !(Master_UID %in% c(nonrep_muids, dup_tas_muids)), ] # (N = 57)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Combine all & add new data #############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
full_output <- rbind(potential_full, df_polys, fill = TRUE) %>% as.data.table()
full_output <- rbind(full_output, df_outside_afr[, names(full_output), with = F])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Collapse new lit extraction #################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# load new lit extraction datasets
extraction <- fread(<<<< FILEPATH REDACTED >>>>)
dossier_extraction <- fread(<<<< FILEPATH REDACTED >>>>)
extraction3 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction4 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction5 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction7 <- fread(<<<< FILEPATH REDACTED >>>>)
extraction8 <- fread(<<<< FILEPATH REDACTED >>>>)

setnames(extraction, c("lat", "long"), c("latitude", "longitude"))
setnames(extraction3, c("lat", "long"), c("latitude", "longitude"))

extraction$row_num <- paste0("ext1_", 1:nrow(extraction))
extraction_raw <- copy(extraction)

dossier_extraction[, location_name := Admn_0]
dossier_extraction[, diagnostic := case_diagnostics]
dossier_extraction[diagnostic == "ICT", diagnostic := "ict"]
dossier_extraction[diagnostic == "no information", diagnostic := NA]
setnames(dossier_extraction, c("prevalence"), c("mean prevalence (%)"))
colnames(extraction)[39] <- "poly_id"
dossier_extraction$row_num <- paste0("dos_", 1:nrow(dossier_extraction))
dossier_extraction$nid <- as.character(dossier_extraction$nid)
dossier_extraction$nid <- "dossier_temp" # replace with final nid
extraction <- rbind(extraction, dossier_extraction, all.x = T, fill = T, use.names = T)

setnames(extraction3, c("mean prevalence (decimal)"), c("mean prevalence (%)"))
extraction3 <- extraction3[c(2:nrow(extraction3)), ]
extraction3$latitude <- as.numeric(extraction3$latitude)
extraction3$longitude <- as.numeric(extraction3$longitude)
extraction3$row_num <- paste0("ext3_", 1:nrow(extraction3))

setnames(extraction4, c("mean prevalence (decimal)", "lat", "long"), c("mean prevalence (%)", "latitude", "longitude"))
extraction4 <- extraction4[c(2:nrow(extraction4)), ]
extraction4$latitude <- as.numeric(extraction4$latitude)
extraction4$longitude <- as.numeric(extraction4$longitude)
extraction4$row_num <- paste0("ext4_", 1:nrow(extraction4))

extraction <- rbind(extraction, extraction3, use.names = T, fill = T)
extraction[, diagnostic := NA]
setnames(extraction4, "location_id", "iso3")
extraction <- rbind(extraction, extraction4, use.names = T, fill = T)

setnames(extraction5, c("mean prevalence (decimal)", "lat", "long"), c("mean prevalence (%)", "latitude", "longitude"))
extraction5$latitude <- as.numeric(extraction5$latitude)
extraction5$longitude <- as.numeric(extraction5$longitude)
extraction5$row_num <- paste0("ext5_", 1:nrow(extraction5))
extraction5 <- extraction5[!is.na(nid)]
setnames(extraction5, "location_id", "iso3")
extraction <- rbind(extraction, extraction5, use.names = T, fill = T)

setnames(extraction7, c("mean prevalence (decimal)", "lat", "long"), c("mean prevalence (%)", "latitude", "longitude"))
extraction7$latitude <- as.numeric(extraction7$latitude)
extraction7$longitude <- as.numeric(extraction7$longitude)
extraction7$row_num <- paste0("ext7_", 1:nrow(extraction7))
setnames(extraction7, c("location_id"), c("iso3"))
extraction7 <- extraction7[!(title == "")]
extraction <- rbind(extraction, extraction7, use.names = T, fill = T)

setnames(extraction8, c("mean prevalence (decimal)", "lat", "long"), c("mean prevalence (%)", "latitude", "longitude"))
extraction8$latitude <- as.numeric(extraction8$latitude)
extraction8$longitude <- as.numeric(extraction8$longitude)
extraction8$row_num <- paste0("ext8_", 1:nrow(extraction8))
setnames(extraction8, c("location_id"), c("iso3"))
extraction <- rbind(extraction, extraction8, use.names = T, fill = T)

### Standardize some fields ###################################################
# delete extraction template explanatory row
extraction <- extraction[c(2:nrow(extraction))]
extraction <- extraction[!is.na(row_num)]

extraction[group == 104, "MDA start" := 2007] # fix likely extraction error
extraction[specificity == "age/sex", specificity := "sex/age"] # standardize on sex/age
extraction[specificity %in% c("sex/diagnostic/age", "diagnostic/age/sex"), specificity := "age/sex/diagnostic"] # standardize on age/sex/diagnostic

extraction <- extraction[is.na(location_name) | !(location_name == "location name from locations tab")]

# Convert location names as necessary
extraction[location_name == "Palau", iso3 := "PLW"]
extraction[location_name == "Niue", iso3 := "NIU"]
extraction[location_name == "Cook Islands", iso3 := "COK"]
extraction[location_name == "Brazil - Pernambuco", location_name := "Brazil"]
extraction[location_name == "India - Gujarat", location_name := "India"]
extraction[location_name == "India - Puducherry", location_name := "India"]
extraction[location_name == "Kenya - Kwale", location_name := "Kenya"]
extraction[location_name == "New Zealand - Cook Islands", location_name := "Cook Islands"]
extraction[location_name == "Nigeria - Kano", location_name := "Nigeria"]
extraction[location_name == "Nigeria - Nasarawa", location_name := "Nigeria"]
extraction[location_name == "Nigeria - Osun", location_name := "Nigeria"]
extraction[location_name == "Nigeria - Plateau", location_name := "Nigeria"]
extraction[location_name == "United States - American Samoa", location_name := "American Samoa"]
extraction[location_name == "Nigeria - Ondo", location_name := "Nigeria"]
extraction[ihme_loc_id == "Fiji|FJI", location_name := "Fiji"]
extraction[ihme_loc_id == "Haiti|HTI", location_name := "Haiti"]
extraction[ihme_loc_id == "India|IND", location_name := "India"]
extraction[ihme_loc_id == "Indonesia|IDN", location_name := "Indonesia"]
extraction[ihme_loc_id == "Papua New Guinea|PNG", location_name := "Papua New Guinea"]
extraction[location_name == "Federal States of Micronesia", location_name := "Federated States of Micronesia"]

adm0s <- sort(unique(extraction$location_name))
if (length(adm0s) != 33) warning("Number of unique countries differ from when I specified the below string matching. Please check and revise as needed.")
source(<<<< FILEPATH REDACTED >>>>)
location_hierarchy <- get_location_metadata(location_set_id = 35, gbd_round_id = 6) # round = 2019

message(paste0(length(adm0s[!(adm0s %in% location_hierarchy$location_name)]), " location_names were unmatched to location_hierarchy. Implementing fix"))
unmatched <- adm0s[!(adm0s %in% location_hierarchy$location_name)]
iso_link <- data.table(location_name=c("Federated States of Micronesia", "New Caledonia", "Wallis and Futuna", "Tanzania", "The Gambia"), ihme_loc_id=c("FSM", "NCL", "WLF", "TZA", "GMB"))
location_link <- location_hierarchy[location_name %in% adm0s, .(location_name, ihme_loc_id)]
location_link <- rbind(location_link, iso_link) %>% as.data.table()
setnames(location_link, "ihme_loc_id", "country")
extraction <- merge(extraction, location_link, by = "location_name")

# standardize some fields
extraction[, source := "sys_rev"]
extraction[shape_type %in% c("Point", "point"), point := 1]
extraction[shape_type %in% c("Polygon", "polygon"), point := 0]
extraction[is.na(point) & tolower(recall_type) == "point", point := 1]
extraction[, latitude := round(as.double(latitude), 5)]
extraction[, longitude := round(as.double(longitude), 5)]

# process polygon reference fields
extraction[grep("_shp", extraction$poly_reference), poly_reference := sub("_shp", ".shp", poly_reference)]
extraction[(is.na(poly_reference) == F | poly_reference == ""), shapefile := sub("\\..*", "", basename(poly_reference))]
extraction[, location_code := poly_id]
extraction[shapefile == "", shapefile := NA]
extraction[location_code == "", location_code := NA]

# process year
extraction[is.na(year_end) == F & is.na(year_start) == F & year_start != 0 & year_end != 0, year_survey := floor((as.integer(year_start) + as.integer(year_end)) / 2)]
message(paste0(nrow(extraction[is.na(year_survey) == T]), " rows missing year_survey after step 1 of year cleaning"))
extraction[is.na(year_survey) == T & is.na(year_start) == F, year_survey := year_start]
message(paste0(nrow(extraction[is.na(year_survey) == T]), " rows missing year_survey after step 2 of year cleaning"))
extraction[is.na(year_survey) == T & is.na(year_end) == F, year_survey := year_end]
message(paste0(nrow(extraction[is.na(year_survey) == T]), " rows missing year_survey after step 3 of year cleaning"))
extraction[year_end == 0 & is.na(year_start) == F & year_start != 0, year_survey := year_start]
extraction[year_start == 0 & is.na(year_end) == F & year_end != 0, year_survey := year_end]

# process diagnostics
cv_fields <- sort(grep("diag_", names(extraction), value = T))
cv_names <- c("br", "mf", "elisa_ab", "elisa_ag", "fts", "ict", "other")
cv_table <- data.table(cv_field = cv_fields, cv_name = cv_names)

for (i in cv_fields) {
  extraction[get(i) == 1, diagnostic := cv_table[which(cv_fields == i), cv_name]]
}
extraction[is.na(diagnostic) == T & (is.na(case_diagnostics) == F & case_diagnostics != ""), diagnostic := case_diagnostics]
extraction[diagnostic == "ELISA", diagnostic := "elisa"]
extraction[grep("blood smear", extraction$diagnostic), diagnostic := "mf"]
extraction[diagnostic == "SD Bioline Onchocerciasis/LF Ig G4 Rapid Test", diagnostic := "ab"] # new IgG4 test for both oncho and LF
extraction[grep("membrane filtration", extraction$diagnostic), diagnostic := "mf"]
extraction[diagnostic == "FTS", diagnostic := "fts"]
extraction[grep("ICT", extraction$diagnostic), diagnostic := "ict"]
extraction[diagnostic == "finger prick blood sample, collected between 20:00 and 24:00, three thick smears of 20 cmm, examined for mf", diagnostic := "mf"]
extraction[diagnostic == "no information" & case_definition == "Mf rate", diagnostic := "mf"]

# process programmatic stage
extraction[(cv_SS == 1 | cv_spotcheck == 1), data_collect_method := "SS/Spot"]
extraction[cv_mapping_survey == 1, data_collect_method := "Mapping"]
extraction[cv_LQAS == 1, data_collect_method := "TAS"]
extraction[is.na(data_collect_method) == T, data_collect_method := "Unknown"]

# numerator & denominator
setnames(extraction, "mean prevalence (%)", "mean")
extraction$cases <- as.integer(extraction$cases)
extraction$sample_size <- as.integer(extraction$sample_size)
extraction$mean <- as.numeric(extraction$mean)
extraction$upper <- as.numeric(extraction$upper)
extraction$lower <- as.numeric(extraction$lower)

check_swap <- extraction[(is.na(sample_size) == T | sample_size == "" | sample_size < 10) & (is.na(cases) == F & cases > sample_size), .(row_num, sample_size, mean, upper, lower, cases)] # check for swapped sample size & cases
if (nrow(check_swap) > 0) { # review and edit based on above output
  extraction[(is.na(sample_size) == T | sample_size == "" | sample_size < 10) & (is.na(cases) == F & cases > sample_size), c("N", "sample_size", "cases") := list(cases, cases, sample_size)]
} else {
  extraction[, N := as.integer(NA)]
}
extraction[(is.na(cases) == T & is.na(mean) == F & is.na(sample_size) == F), cases := round((mean / 100) * sample_size, 0)]
extraction[(is.na(sample_size) == F & sample_size != "" & sample_size >= cases & is.na(N) == T), N := sample_size]
table(is.na(extraction$N))

extraction <- extraction[!is.na(N)]

### Check "issue" fields #####################################################
# View in RStudio is preferred for reviewing these
message(paste0(nrow(extraction[year_issue == 1]), " rows with year_issues"))
table(extraction[year_issue == 1, note_extracter])
message(paste0(nrow(extraction[age_issue == 1]), " rows with age_issues"))
table(extraction[age_issue == 1, note_extracter])
message(paste0(nrow(extraction[sex_issue == 1]), " rows with sex_issues"))
table(extraction[sex_issue == 1, note_extracter])

other_issues <- c("cv_Slum", "cv_Pregnant", "cv_Indigenous")
for (issue in other_issues) {
  if (nrow(extraction[get(issue) == 1]) > 1) message(paste0(nrow(extraction[get(issue) == 1]), " rows marked as ", issue, ". Check on whether these should be excluded"))
}

extraction[stri_detect_fixed(row_num, "ext8_"), group := paste0(group, "_ext8")]

### Resolve multiple specificity within the same group ################################
extraction <- data.table(extraction)
groups <- extraction[!(group %in% c("", NA)), group] %>% unique()
spec_num <- c()
for (i in groups) {
  spec_num <- c(spec_num, length(unique(extraction[group == i, specificity])))
}
groups <- data.table(group = groups, spec = spec_num)
combo_look <- extraction[group %in% groups[spec > 1, group], .(group, specificity)] %>%
  table() %>%
  as.data.table() %>%
  .[N > 1]
groups <- groups[spec > 1, ] %>% as.data.table()
groups$spec_keep <- lapply(unique(groups$group), FUN = function(x) {
  if ("location" %in% combo_look[group == x, specificity]) { # preferentially keep location specificity
    keep <- "location"
  } else { # otherwise no preference
    keep <- combo_look[group == x, specificity][1]
    if (keep == "diagnostic") keep <- combo_look[group == x, specificity][2]
  }
  return(keep)
}) %>% unlist()

keep_specs <- lapply(1:nrow(groups), FUN = function(x) {
  g <- groups[x, group]
  k <- groups[x, spec_keep]
  extraction[group == g & specificity == k, ]
}) %>% rbindlist()

extraction <- rbind(extraction[!(group %in% unique(groups$group))], keep_specs, use.names = T)

### Subset rows by specificity ###############################################
# identify and weed out unused/duplicate fields
unused_fields <- c()
for (i in names(extraction)) {
  if ((nrow(extraction[(get(i) == "" | is.na(get(i)) == T), ]) / nrow(extraction)) == 1) unused_fields <- c(unused_fields, i)
}

extraction <- subset(extraction, select = names(extraction)[!(names(extraction) %in% c(unused_fields))])

extraction[group == 240 & latitude == 2.850616, specificity := "none"] # fix location specificity error
extraction[group == 240 & latitude == 2.850616, group := ""] # fix location specificity error

processing <- list()
processing[["incorrect"]] <- extraction[(specificity %in% c(0:100) | specificity %in% as.character(c(0:100))), ]
warning(paste0(nrow(processing[["incorrect"]]), " rows with irregular specificity (numeric). Please check for extraction error and fix"))
extraction[(specificity == "" | is.na(specificity)), specificity := "none"]
for (i in unique(extraction[!(specificity %in% c(0:100) | specificity %in% as.character(c(0:100))), specificity])) {
  processing[[as.character(i)]] <- extraction[specificity == i & !(row_num %in% processing[["incorrect"]]$row_num), ]
}
message("Specificity in extraction:")
lapply(processing, nrow) %>%
  unlist() %>%
  print()
processing_copy <- copy(processing)

# ### age specificity ##########################################################
temp <- copy(processing[["age"]])
if (nrow(temp[(is.na(group) == T | group == ""), ]) > 0) warning("Error: Group field was incorrectly filled")
varying_fields <- c("row_num", "mean", "lower", "upper", "cases", "N", "unique", "age_start", "age_end", "Master_UID")
nonvarying_fields <- c(
  "nid", "group", "source", "country", "data_collect_method", "point", "shapefile", "location_code",
  "latitude", "longitude", "geo_accuracy", "geo_precision", "year_survey", "MDA", "MDA start", "diagnostic", "sex"
)

if (nrow(unique(temp[, nonvarying_fields, with = F])) != length(unique(temp$group))) {
  warning("Fields that should be non-varying are varying in some groups. Check and fix")
} else {
  processing[["age"]] <- unique(temp[, nonvarying_fields, with = F])
  for (g in unique(temp$group)) {
    min_age <- min(as.integer(temp[group == g, age_start]), na.rm = T)
    max_age <- max(as.integer(temp[group == g, age_end]), na.rm = T)
    total_n <- sum(as.integer(temp[group == g, N]))
    total_cases <- sum(as.integer(temp[group == g, cases]))
    row_nums <- paste(temp[group == g, row_num], collapse = ";")
    Master_UIDs <- paste(temp[group == g, Master_UID], collapse = ";")
    processing[["age"]][group == g, c("N", "cases", "age_start", "age_end", "row_num", "Master_UID") := list(total_n, total_cases, min_age, max_age, row_nums, Master_UIDs)]
  }
}

### sex specificity ##########################################################
temp <- copy(processing[["sex"]])
if (nrow(temp[(is.na(group) == T | group == ""), ]) > 0) warning("Error: Group field was incorrectly filled")
varying_fields <- c("row_num", "mean", "lower", "upper", "cases", "N", "unique", "sex", "Master_UID")
nonvarying_fields <- c(
  "nid", "group", "source", "country", "data_collect_method", "point", "shapefile", "location_code",
  "latitude", "longitude", "geo_accuracy", "geo_precision", "year_survey", "MDA", "MDA start", "diagnostic", "age_start", "age_end"
)

if (nrow(unique(temp[, nonvarying_fields, with = F])) != length(unique(temp$group))) {
  warning("Fields that should be non-varying are varying in some groups. Check and fix")
  check_group <- unique(temp[, nonvarying_fields, with = F])
  check_group$check_uid <- c(1:nrow(check_group))
  check_group <- merge(temp, check_group, by = nonvarying_fields, all = T)
  problematic_groups <- table(check_group[, .(check_uid, group)]) %>%
    as.data.table() %>%
    .[N == 1, group] %>%
    unique()
  warning(paste0("check these groups: ", problematic_groups))
} else {
  processing[["sex"]] <- unique(temp[, nonvarying_fields, with = F])
  for (g in unique(temp$group)) {
    total_n <- sum(as.integer(temp[group == g, N]))
    total_cases <- sum(as.integer(temp[group == g, cases]))
    row_nums <- paste(temp[group == g, row_num], collapse = ";")
    Master_UIDs <- paste(temp[group == g, Master_UID], collapse = ";")
    processing[["sex"]][group == g, c("N", "cases", "sex", "row_num", "Master_UID") := list(total_n, total_cases, "both", row_nums, Master_UIDs)]
  }
}

### sex/age specificity ##########################################################
temp <- copy(processing[["sex/age"]])
if (nrow(temp[(is.na(group) == T | group == ""), ]) > 0) warning("Error: Group field was incorrectly filled")
varying_fields <- c("row_num", "mean", "lower", "upper", "cases", "N", "unique", "sex", "age_start", "age_end", "Master_UID")
nonvarying_fields <- c(
  "nid", "group", "source", "country", "data_collect_method", "point", "shapefile", "location_code",
  "latitude", "longitude", "geo_accuracy", "geo_precision", "year_survey", "MDA", "MDA start", "diagnostic"
)

if (nrow(unique(temp[, nonvarying_fields, with = F])) != length(unique(temp$group))) {
  warning("Fields that should be non-varying are varying in some groups. Check and fix")
} else {
  processing[["sex/age"]] <- unique(temp[, nonvarying_fields, with = F])
  for (g in unique(temp$group)) {
    min_age <- min(as.integer(temp[group == g, age_start]), na.rm = T)
    max_age <- max(as.integer(temp[group == g, age_end]), na.rm = T)
    total_n <- sum(as.integer(temp[group == g, N]))
    total_cases <- sum(as.integer(temp[group == g, cases]))
    row_nums <- paste(temp[group == g, row_num], collapse = ";")
    Master_UIDs <- paste(temp[group == g, Master_UID], collapse = ";")
    processing[["sex/age"]][group == g, c("N", "cases", "sex", "row_num", "age_start", "age_end", "Master_UID") := list(total_n, total_cases, "both", row_nums, min_age, max_age, Master_UIDs)]
  }
}

### diagnostic specificity ##########################################################
# NOTE: selection criteria for this is highly case-specific
temp <- copy(processing[["diagnostic"]])

# drop group 128 (improper extraction)
temp <- temp[!(group == 128), ]

### Drop membrane filtration row for group 273
temp <- temp[!((group == 273) & (case_diagnostics == "membrane filtration"))]

if (nrow(temp[(is.na(group) == T | group == ""), ]) > 0) warning("Error: Group field was incorrectly filled")
varying_fields <- c("row_num", "mean", "lower", "upper", "cases", "N", "unique", "diagnostic", "Master_UID")
nonvarying_fields <- c(
  "nid", "group", "source", "country", "data_collect_method", "point", "shapefile", "location_code",
  "latitude", "longitude", "geo_accuracy", "geo_precision", "year_survey", "MDA", "MDA start", "sex", "age_start", "age_end"
)

### Skip this check for now (i.e., run the code in the else statement) due to strange behavior for group 142
if (nrow(unique(temp[, nonvarying_fields, with = F])) != length(unique(temp$group))) {
  warning("Fields that should be non-varying are varying in some groups. Check and fix")
} else {
  processing[["diagnostic"]] <- unique(temp[, nonvarying_fields, with = F])
  for (g in unique(temp$group)) {
    # for this situation: All groups are antibody ELISA (BM-14) and test strip; both for W. Bancrofti; FTS preferred
    # Update: One group has ICT; prefer this over FTS
    if ("ict" %in% temp[group == g, diagnostic]) {
      processing[["diagnostic"]][group == g, diagnostic := "ict"]
      processing[["diagnostic"]][group == g, N := temp[group == g & diagnostic == "ict", N]]
      processing[["diagnostic"]][group == g, row_num := temp[group == g & diagnostic == "ict", row_num]]
      processing[["diagnostic"]][group == g, cases := temp[group == g & diagnostic == "ict", cases]]
      processing[["diagnostic"]][group == g, Master_UID := temp[group == g & diagnostic == "ict", Master_UID]]
    } else {
      if ("fts" %in% temp[group == g, diagnostic]) {
        processing[["diagnostic"]][group == g, diagnostic := "fts"]
        processing[["diagnostic"]][group == g, N := temp[group == g & diagnostic == "fts", N]]
        processing[["diagnostic"]][group == g, row_num := temp[group == g & diagnostic == "fts", row_num]]
        processing[["diagnostic"]][group == g, cases := temp[group == g & diagnostic == "fts", cases]]
        processing[["diagnostic"]][group == g, Master_UID := temp[group == g & diagnostic == "fts", Master_UID]]
      } else {
        if ("mf" %in% temp[group == g, diagnostic]) {
          processing[["diagnostic"]][group == g, diagnostic := "mf"]
          processing[["diagnostic"]][group == g, N := temp[group == g & diagnostic == "mf", N]]
          processing[["diagnostic"]][group == g, row_num := temp[group == g & diagnostic == "mf", row_num]]
          processing[["diagnostic"]][group == g, cases := temp[group == g & diagnostic == "mf", cases]]
          processing[["diagnostic"]][group == g, Master_UID := temp[group == g & diagnostic == "mf", Master_UID]]
        }
      }
    }
  }
}

### age/diagnostic specificity ##########################################################
# NOTE: selection criteria for this is highly case-specific
temp <- copy(processing[["age/diagnostic"]])
if (nrow(temp[(is.na(group) == T | group == ""), ]) > 0) warning("Error: Group field was incorrectly filled")
varying_fields <- c("row_num", "mean", "cases", "N", "unique", "age_start", "age_end", "Master_UID")
nonvarying_fields <- c(
  "nid", "group", "source", "country", "data_collect_method", "point", "shapefile", "location_code",
  "latitude", "longitude", "geo_accuracy", "geo_precision", "year_survey", "MDA", "MDA start", "sex", "diagnostic"
)

processing[["age/diagnostic"]] <- unique(temp[, nonvarying_fields, with = F])
for (g in unique(temp$group)) {
  # for this situation: All groups are mf and br (differnt species); keep both
  for (d in unique(temp[group == group, diagnostic])) {
    min_age <- min(as.integer(temp[group == g & diagnostic == d, age_start]), na.rm = T)
    max_age <- max(as.integer(temp[group == g & diagnostic == d, age_end]), na.rm = T)
    total_n <- sum(as.integer(temp[group == g & diagnostic == d, N]))
    total_cases <- sum(as.integer(temp[group == g & diagnostic == d, cases]))
    row_nums <- paste(temp[group == g & diagnostic == d, row_num], collapse = ";")
    Master_UIDs <- paste(temp[group == g, Master_UID], collapse = ";")
    processing[["age/diagnostic"]][group == g & diagnostic == d, c("N", "cases", "row_num", "age_start", "age_end", "diagnostic", "Master_UID") := list(total_n, total_cases, row_nums, min_age, max_age, d, Master_UIDs)]
  }
}

### age/sex/diagnostic specificity ###########################################################################
# NOTE: selection criteria for this is highly case-specific
temp <- copy(processing[["age/sex/diagnostic"]])
if (nrow(temp[(is.na(group) == T | group == ""), ]) > 0) warning("Error: Group field was incorrectly filled")
varying_fields <- c("row_num", "mean", "lower", "upper", "cases", "N", "unique", "age_start", "age_end", "sex", "Master_UID")
nonvarying_fields <- c(
  "nid", "group", "source", "country", "data_collect_method", "point", "shapefile", "location_code",
  "latitude", "longitude", "geo_accuracy", "geo_precision", "year_survey", "MDA", "MDA start", "diagnostic"
)

processing[["age/sex/diagnostic"]] <- unique(temp[, nonvarying_fields, with = F])
for (g in unique(temp$group)) {
  if ("ict" %in% temp[group == g, diagnostic]) {
    min_age <- min(as.integer(temp[group == g & diagnostic == "ict", age_start]), na.rm = T)
    max_age <- max(as.integer(temp[group == g & diagnostic == "ict", age_end]), na.rm = T)
    total_n <- sum(as.integer(temp[group == g & diagnostic == "ict", N]))
    total_cases <- sum(as.integer(temp[group == g & diagnostic == "ict", cases]))
    row_nums <- paste(temp[group == g & diagnostic == "ict", row_num], collapse = ";")
    Master_UIDs <- paste(temp[group == g, Master_UID], collapse = ";")
    processing[["age/sex/diagnostic"]][group == g, c("N", "cases", "row_num", "age_start", "age_end", "diagnostic", "sex", "Master_UID") := list(total_n, total_cases, row_nums, min_age, max_age, "ict", "both", Master_UIDs)]
  }
  if (!("ict" %in% temp[group == g, diagnostic]) & ("mf" %in% temp[group == g, diagnostic])) {
    min_age <- min(as.integer(temp[group == g & diagnostic == "mf", age_start]), na.rm = T)
    max_age <- max(as.integer(temp[group == g & diagnostic == "mf", age_end]), na.rm = T)
    total_n <- sum(as.integer(temp[group == g & diagnostic == "mf", N]))
    total_cases <- sum(as.integer(temp[group == g & diagnostic == "mf", cases]))
    row_nums <- paste(temp[group == g & diagnostic == "mf", row_num], collapse = ";")
    Master_UIDs <- paste(temp[group == g, Master_UID], collapse = ";")
    processing[["age/sex/diagnostic"]][group == g, c("N", "cases", "row_num", "age_start", "age_end", "diagnostic", "sex", "Master_UID") := list(total_n, total_cases, row_nums, min_age, max_age, "ict", "both", Master_UIDs)]
  }
}

### location & no specificity ##################################################################################
# for this extraction: location specificity is equivalent to no specificity... filled in incorrectly
processing[["none"]] <- rbind(processing[["none"]], processing[["location"]])

### Subset to consistent fields with other list items and rbindlist back together ###############################
processing[["none"]] <- processing[["none"]][, names(processing[["sex"]]), with = F]

if (length(unique(unlist(lapply(processing[names(processing)[!(names(processing) %in% c("incorrect", "location"))]], FUN = function(x) length(names(x)))))) == 1) {
  processed <- rbindlist(processing[names(processing)[!(names(processing) %in% c("incorrect", "location"))]], use.names = T, fill = T)
} else {
  warning("# of fields between different specificity subsets do not line up. Please review and fix.")
}

### process to include standard fields for combining with potential_full object in clean_espen.R prior to crosswalk ##############################
processed <- rbindlist(processing, use.names = T, fill = T)
processed$age_start <- as.integer(processed$age_start)
processed$age_end <- as.integer(processed$age_end)
age_cats <- unique(processed[, .(age_start, age_end)])
age_cats[age_start >= 15 & age_end >= 15, age_group := "Adults"]
age_cats[age_start < 15 & age_end <= 15, age_group := "Children"]
age_cats[is.na(age_group), age_group := "Adults/Children"]
processed <- merge(processed, age_cats, by = c("age_start", "age_end"), all.x = T)

processed[, c("cluster_id", "lf_prev", "original_year", "resamp_test_case", "weight", "sampling") := list(1, cases / N, year_survey, NA, 1, 1)]
processed[, row_num := paste0("kd_lit_", row_num)]
setnames(processed, c("cases", "year_survey"), c("had_lf_poly", "year"))

# report out changes from data input
unique_obs <- lapply(location_link$country, FUN = function(x) (length(unique(extraction[country == x & group != "" & specificity != "location", group], na.rm = T)) + nrow(extraction[country == x & (is.na(group) == T | group == "" | specificity == "location"), ]))) %>% unlist()
summary <- data.table(country = location_link$country, input_obs = unique_obs)
summary$processed_obs <- lapply(location_link$country, FUN = function(x) nrow(processed[country == x])) %>% unlist()

if (nrow(summary[processed_obs < input_obs]) > 0) {
  message("CAUTION: some data might have been dropped during processing for the following countries. Please check.")
  print(summary[processed_obs < input_obs])
}

# NOTE: in summary object, if processed_obs > input_obs, it's probably because there's specificity in group by diagnostic in which case we currently retain all (because most of these multi-diagnostic extractions are BR + ict/mf)
multi_diag <- unique(processed[group != "" & is.na(group) == F, group])[which(unlist(lapply(unique(processed[group != "" & is.na(group) == F, group]), FUN = function(x) length(unique(processed[group == x, diagnostic])))) != 1)]
if (length(multi_diag) > 0) {
  message("The following groups have multiple diagnostics remaining after processing")
  for (i in multi_diag) {
    print(i)
    print(unique(processed[group == i, diagnostic]))
  }
}

# Fix groups still containing multiple diagnostics
processed <- processed[!(group == 102 & diagnostic == "mf")]
processed <- processed[!(group == 106 & diagnostic == "mf")]
processed <- processed[!(group == 107 & diagnostic == "mf")]
processed <- processed[!(group == 312 & diagnostic == "mf")]
processed <- processed[!(group == 313 & diagnostic == "mf")]
processed <- processed[!(group == 302 & diagnostic == "Wb123 multiplex bead assay using recombinant antigens of W. bancrofti")]
processed <- processed[!(group == 314 & diagnostic == "mf")]

write.csv(processed, file = paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)

full_output <- rbind(full_output, processed, use.names = TRUE, fill = TRUE)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Add new TAS data ########################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
new_TAS <- fread(<<<< FILEPATH REDACTED >>>>)
new_TAS$data_collect_method <- "TAS"
new_TAS$source <- "WHO"
setnames(new_TAS, c("lat", "long", "examined", "positive", "poly_reference", "poly_id"), c("latitude", "longitude", "N", "had_lf_poly", "shapefile", "location_code"))
new_TAS$N <- as.integer(gsub(",", "", new_TAS$N))
new_TAS$lf_prev <- new_TAS$had_lf_poly / new_TAS$N
new_TAS[, c("original_year", "age_group", "weight", "resamp_test_case", "cluster_id", "sampling", "nid") := list(year, "", 1, NA, 1, NA, NA)]

new_TAS[diagnostic == "FTS (Ag)", diagnostic := "FTS"]

cols <- grep("V", names(new_TAS), value = TRUE, invert = TRUE)
new_TAS <- new_TAS[, ..cols]

full_output <- rbind(full_output, new_TAS, use.names = TRUE, fill = TRUE)

full_output[, lf_prev := had_lf_poly / N]

# add new extractions
asia_extract <- fread(<<<< FILEPATH REDACTED >>>>)
asia_extract[, c("year", "cluster_id") := list(original_year, 1)]
asia_extract <- asia_extract[point == 0]
asia_extract$latitude <- as.numeric(asia_extract$latitude)
asia_extract$longitude <- as.numeric(asia_extract$longitude)
full_output <- rbind(full_output, asia_extract, use.names = TRUE, fill = TRUE)

india_review <- fread(<<<< FILEPATH REDACTED >>>>)
india_review[, c("year", "cluster_id", "original_year", "weight", "resamp_test_case", "sampling", "source", "age_group")
             := list(year_start, NA, year_start, 1, NA, NA, "sys_rev", "Adults/Children")]
setnames(india_review, c("NID", "N_pos"), c("nid", "had_lf_poly"))
full_output <- rbind(full_output, india_review, use.names = TRUE, fill = TRUE)

# Add new data from taskforce
taskforce_data <- fread(<<<< FILEPATH REDACTED >>>>)
taskforce_data[, resamp_test_case := NA]
taskforce_data$nid <- NULL
taskforce_data[country == "PHL", nid := "412563"]
taskforce_data[country == "TZA", nid := "412585"]
taskforce_data[country == "HTI", nid := "412586"]
full_output <- rbind(full_output, taskforce_data, use.names = TRUE, fill = TRUE)

### Add new TAS data
TAS_2018 <- fread(<<<< FILEPATH REDACTED >>>>)
lookup_table <- fread(<<<< FILEPATH REDACTED >>>>)
TAS_2018 <- merge(TAS_2018, lookup_table[, c("iso3", "iso2")], by.x = "Country", by.y = "iso2", all.x = TRUE)
TAS_2018[!is.na(pop_test_strip) | (pop_test_strip == 0 & sample_test_strip > 0), c("diagnostic") := list("FTS")]
TAS_2018[is.na(invalid_test_strip), invalid_test_strip := 0]
TAS_2018[!is.na(sample_test_strip) & pop_test_strip == 0 & sample_test_strip > 0, pop_test_strip := sample_test_strip]
TAS_2018[, N := pop_test_strip - invalid_test_strip]
TAS_2018 <- TAS_2018[!is.na(N)]
setnames(TAS_2018, c("iso3", "np_test_strip", "Year_survey", "lat", "long"), c("country", "had_lf_poly", "year", "latitude", "longitude"))
TAS_2018[, c("nid", "data_collect_method", "source", "weight") := list(327584, "TAS", "WHO", 1)]
TAS_2018[, lf_prev := had_lf_poly / N]
TAS_2018[is.na(latitude) & is.na(longitude), point := 0]
TAS_2018[!is.na(latitude) & !is.na(longitude), point := 1]
TAS_2018[, sampling := NULL]
TAS_2018[, sampling := 1]
TAS_2018 <- TAS_2018[!is.na(had_lf_poly) & !is.na(N)]
setnames(TAS_2018, c("shapefile_ref", "GAUL_CODE"), c("shapefile", "location_code"))

full_output <- rbind(full_output, TAS_2018, use.names = TRUE, fill = TRUE)

### Subset columns
full_output_backup <- copy(full_output)
full_output <- full_output[, c("Master_UID", "nid", "data_collect_method", "source", "latitude", "longitude", "year", "N", "had_lf_poly", "lf_prev", "country", "diagnostic", "age_start", "age_end", "point", "age_group", "weight", "shapefile", "location_code", "resamp_test_case", "cluster_id", "sampling")]

# Fix Master_UID column
full_output[, Master_UID_new := do.call(paste, .SD), .SDcols = "Master_UID"]
full_output$Master_UID_new <- gsub("c\\(", "", full_output$Master_UID_new)
full_output$Master_UID_new <- gsub("\\)", "", full_output$Master_UID_new)
full_output$Master_UID_new <- gsub(",", ";", full_output$Master_UID_new)
full_output$Master_UID_new <- gsub("; ", ";", full_output$Master_UID_new)
full_output$Master_UID_new <- gsub(";", "; ", full_output$Master_UID_new)
full_output$Master_UID <- full_output$Master_UID_new
full_output$Master_UID_new <- NULL


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Small corrections #######################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Fix coord
full_output[Master_UID == 17037, latitude := 17.569676]
full_output[Master_UID == 17037, longitude := 97.912934]

## Corrections to age extractions
full_output[country == "BGD" & Master_UID %in% 70:126, c("age_start", "age_end") := list(6, 7)]
full_output[country == "NGA" & Master_UID %in% c(14397:14398, 21880:21883), c("age_start", "age_end") := list(1, 99)]
full_output[country == "NGA" & Master_UID %in% c(23930:23933), c("age_start", "age_end") := list(6, 7)]

# Standardize diagnostic categories
full_output[full_output$diagnostic %in% c("Blood smear", "mf"), diagnostic := "adjusted_mf"]
full_output[full_output$diagnostic %in% c("ict", "ICT"), diagnostic := "ICT"]
full_output[full_output$diagnostic %in% c("FTS", "fts", "test_strip"), diagnostic := "FTS"]

### Where cases are not reported but prevalence is, calculate cases
full_output[is.na(had_lf_poly) & !is.na(lf_prev), had_lf_poly := round(lf_prev * N, 0)]

# Check that point field is properly set; change where necessary
full_output[!is.na(latitude) & !is.na(longitude) & !(latitude == "") & !(longitude == "") & point == 0, point := 1]
full_output[!is.na(latitude) & !is.na(longitude) & !(latitude == "") & !(longitude == "") & is.na(point), point := 1]

# Post-hoc shapefile corrections
full_output[shapefile == "gadm28_adm0" & country == "COK", c("shapefile", "location_code") := list("gadm_36_ad0", "49")]
full_output[shapefile == "gadm28_adm0" & country == "NIU", c("shapefile", "location_code") := list("gadm_36_ad0", "166")]
full_output[shapefile == "gadm28_adm0" & country == "PLW", c("shapefile", "location_code") := list("gadm_36_ad0", "178")]
full_output[shapefile == "gadm28_adm0" & country == "WLF", c("shapefile", "location_code") := list("gadm_36_ad0", "244")]
full_output[shapefile == "lf_gadm28_adm2" & country == "MLI" & location_code %in% c(57500, 57532, 57549, 57705, 57711), shapefile := "lf_gadm28_adm3"]

# Reference TZA MOH data to ESPEN instead (exists in there too)
full_output[source == "country" & country == "TZA", source := "ESPEN"]

# Assign NIDs for large data sources (ESPEN, WHO, Taskforce)
full_output[grep("DOS", Master_UID), c("source", "nid") := list("WHO", NA)]
full_output[nid == 9999, c("source", "nid") := list("WHO", NA)]
full_output[is.na(nid) == T & source == "WHO", nid := 327584]

# ESPEN (nids by country and point/polygon)
espen_nids <- fread(<<<< FILEPATH REDACTED >>>>)
espen_nids <- espen_nids[, .(Nid, ihme_loc_id, `Secondary data type`)]
espen_nids[`Secondary data type` == "", point := 0]
espen_nids[is.na(point), point := 1]
setnames(espen_nids, c("Nid", "ihme_loc_id"), c("nid", "country"))
full_output[source == "Espen", source := "ESPEN"]
output_espen <- full_output[source == "ESPEN",]
output_espen$nid <- NULL
output_espen <- merge(output_espen, espen_nids[,.(nid, country, point)], by = c("country", "point"), all.x = T)
output_espen[stri_detect_fixed(Master_UID, "espen_ID"), nid := 431787]
full_output <- rbind(full_output[source != "ESPEN",], output_espen, use.names = T, fill = T)

# Use general GAHI NID for remaining rows without underlying NIDs
full_output[(is.na(nid) | (nid == "")) & source == "GAHI", nid := "143009"]

full_output[nid == "", nid := NA]

# Review India TAS data
full_output[source == "WHO_TAS", c("source", "nid") := list("WHO", "327584")]
full_output[source == "gahi_both", c("source", "nid") := list("GAHI", "143009")]

full_output <- full_output[!(country == "IND" & data_collect_method == "SS/Spot" & N > 980000)]

full_output[diagnostic %in% c("Filtration", "Knott"), diagnostic := "adjusted_mf"]
full_output <- full_output[diagnostic %in% c("adjusted_mf", "br", "FTS", "ICT")]
full_output <- full_output[year > 1989]
full_output_backup <- copy(full_output)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Check that data have the correct country assigned
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gaul_list <- get_adm0_codes("lf_endem_afr", shapefile_version = shapefile_version)
gaul_list <- c(gaul_list, get_adm0_codes("lf_s_asia", shapefile_version = shapefile_version))
gaul_list <- c(gaul_list, get_adm0_codes("lf_se_asia", shapefile_version = shapefile_version))
gaul_list <- c(gaul_list, get_adm0_codes("lf_non_mbg", shapefile_version = shapefile_version))
gaul_list <- c(gaul_list, get_adm0_codes("lf_hispaniola", shapefile_version = shapefile_version))

simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = F, shapefile_version = shapefile_version)
subset_shape <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]

modeling_shapefile_version <- shapefile_version

## Build administrative and population rasters
raster_list <- build_simple_raster_pop(subset_shape)
simple_raster <- raster_list[["simple_raster"]]

# Set TZZ to TZA
full_output[country == "TZZ", country := "TZA"]

full_output$ihme_loc_id <- full_output$country

nrow(full_output) # 17240

l <- identify_wrong_adm0(data = full_output, s_poly = subset_shape, rast = simple_raster)
full_output <- l[[1]]
full_output_wrong_adm0 <- l[[2]]

# Separate out those to retain and those to drop
retain <- full_output_wrong_adm0[country %in% c("ASM", "COK", "FJI", "FSM", "MDV", "MHL", "NCL", "NIU", "NRU", "PYF", "TON", "TUV", "VUT", "WLF")]
drop <- full_output_wrong_adm0[!(country %in% c("ASM", "COK", "FJI", "FSM", "MDV", "MHL", "NCL", "NIU", "NRU", "PYF", "TON", "TUV", "VUT", "WLF"))]

# Retain Burundi data erroneously referenced to Burkina Faso
drop[country == "BFA" & Master_UID %in% full_output_wrong_adm0[country == "BFA" & ihme_lc_id.x == "BDI" & source == "ESPEN" & year == 2007, Master_UID], c("country", "latitude", "longitude") := list("BDI", latitude_old, longitude_old)]
retain <- rbind(retain, drop[!is.na(latitude) & !is.na(longitude)])

# Add back retained data
full_output <- rbind(full_output, retain, fill = TRUE, use.names = TRUE)

### Fix georeferencing errors
# Extraction labeled as Democratic Republic of Congo but actually Republic of Congo
full_output[nid == 389675, country := "COG"]
# Benin TAS data with faulty georeferencing; assign to polygon of admin2
full_output[Master_UID %in% c("espen_ID_294", "espen_ID_277"), c("point", "shapefile", "location_code", "latitude", "longitude") := list(0, "lf_gadm28_adm2", 4348, NA, NA)]
full_output[Master_UID %in% c("espen_ID_220"), c("point", "shapefile", "location_code", "latitude", "longitude") := list(0, "lf_gadm28_adm2", 4328, NA, NA)]
# South Sudan point with erroneous coords and duplicate Master_UID
full_output <- full_output[!(Master_UID == "temp_7841" & latitude < 0)]
# Fix erreneous coords for a Malaysia data point
full_output[Master_UID %in% c("12712"), c("latitude", "longitude") := list(5.566159, 100.662761)]
# Drop duplicate data
full_output <- full_output[!(Master_UID %in% c("13815", "13778"))]
# Drop PHL data with erroneous coords
full_output <- full_output[!(Master_UID %in% c("15038", "15039"))]
# Fix wrong country
full_output <- full_output[Master_UID %in% c("SR_A_295", "SR_A_312"), country := "WSM"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### De-duplication #############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Systematically assess duplicates by country
## Drop countries that we're not modeling
full_output <- full_output[!(country %in% c("BWA", "CHN", "NRU"))]

# Remove truly duplicate rows
full_output <- unique(full_output)

table(full_output$country)

full_output$latitude <- as.numeric(full_output$latitude)
full_output$longitude <- as.numeric(full_output$longitude)

full_output$latitude_rounded <- round(full_output$latitude, 2)
full_output$longitude_rounded <- round(full_output$longitude, 2)
full_output$latitude_rounded_1 <- round(full_output$latitude, 1)
full_output$longitude_rounded_1 <- round(full_output$longitude, 1)

#### Now deduplicate by country
## AGO
# Only one data row; nothing to be done

## ASM
check_dups <- create_check_dups("ASM")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% c("29"))]
full_output[(Master_UID %in% c("SR_A_312", "SR_A_295")), country := "WSM"] # Fix incorrect country

## BEN
full_output <- full_output[!(country == "BEN" & stri_detect_fixed(Master_UID, "new_TAS"))] # Drop new_TAS rows for Benin (now present in ESPEN data)
check_dups <- create_check_dups("BEN")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]

## BDI
check_dups <- create_check_dups("BDI")

## BFA
check_dups <- create_check_dups("BFA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("1346", "1166", "1194", "1145", "1343", "1170", "1320", "1098", "1211"))] # Manually deduplicate ESPEN and old WHO rows with identical N, cases, and year (within 1 year) and close coords (assume changes in georeferencing and correction to year column)

## BGD
check_dups <- create_check_dups("BGD")

## BRN
check_dups <- create_check_dups("BRN")

## BRA
check_dups <- create_check_dups("BRA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output[Master_UID == "825", N := 23673]
full_output <- full_output[!(Master_UID %in% c("825"))] # Drop sys_rev polygon data that are replicated by GAHI point data
full_output <- full_output[!(nid %in% c(288875, 288880))] # Drop sys_rev data that have incorrect year but appear to be based on old data that are already represented in our data set

## CAF
check_dups <- create_check_dups("CAF")
full_output <- full_output[!(Master_UID %in% c("temp_1420"))] # Drop old:new ESPEN duplicate (with year change)

## CIV
check_dups <- create_check_dups("CIV")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("3132", "3184", "3053", "3112"))] # Drop apparent ESPEN Mapping:WHO SS/SC duplicates (the latter are indicated as being ICT), plus other apparent duplicates

## CMR
check_dups <- create_check_dups("CMR")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% c("1782"))] # retain other obv_dups as they have the same source and non-identical coords
full_output <- full_output[!(Master_UID %in% c("2454"))] # drop WHO data point without survey data type, identical sample size and cases as other data with identical coords, but different year

## COD
check_dups <- create_check_dups("COD")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]] # Don't use
check_dups <- create_check_dups_full_coords_no_year("COD")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "country", "year"), preferred_source = "GAHI")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups & data_collect_method == "Mapping" & source == "WHO")]
check_dups <- create_check_dups_full_coords_no_year("COD")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups & data_collect_method == "Mapping")]

# Correct 2015 COD mapping data
library(raster)
drc_mapping_2015 <- full_output[year == 2015 & country == "COD" & data_collect_method == "Mapping" & !is.na(longitude) & !is.na(latitude), c("longitude", "latitude", "source", "Master_UID")]
lfmda_2014 <- raster(<<<< FILEPATH REDACTED >>>>)
locations <- SpatialPoints(coords = as.matrix(drc_mapping_2015[, .(longitude,latitude)]), proj4string = CRS(proj4string(lfmda_2014)))
drc_mapping_2015$mda_rounds <- raster::extract(lfmda_2014, locations)
full_output <- full_output[!(Master_UID %in% drc_mapping_2015[mda_rounds == 1, Master_UID])]

## COG
check_dups <- create_check_dups("COG")

## COK
check_dups <- create_check_dups("COK")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
obv_dups <- dups[[1]] # Don't use
full_output <- full_output[!(Master_UID %in% c("DOSS_A_41", "2975", "DOSS_A_38", "DOSS_A_40", "2957", "DOSS_A_27", "DOSS_A_32", "DOSS_A_35"))]

## COM
check_dups <- create_check_dups("COM")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

## DOM
check_dups <- create_check_dups("DOM")

## EGY
check_dups <- create_check_dups("EGY")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
check_dups <- create_check_dups_poly("EGY")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("5210", "4930", "4833", "5074", "4931", "4985", "4832", "4831"))] # Drop manually-identified duplicates (including WHO TAS surveys that differ by one year from sys_rev TAS data; prefer the latter)

## ERI
check_dups <- create_check_dups("ERI")

## ETH
check_dups <- create_check_dups("ETH")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

## FJI
check_dups <- create_check_dups("FJI")
check_dups <- create_check_dups_poly("FJI")

## FSM
check_dups <- create_check_dups("FSM")
check_dups <- create_check_dups_poly("FSM")

## GAB
check_dups <- create_check_dups("GAB")

## GHA
check_dups <- create_check_dups("GHA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("7805", "7806", "8006", "8630", "7104", "SR_A_813", "SR_A_817", "7701", "7698", "7702", "7638", "7909", "8033", "8013", "8019", "SR_A_814", "SR_A_806"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("GHA")
full_output <- full_output[!(stri_detect_fixed(Master_UID, "SR_A") & nid == 408334),] # Remove duplicate GHA extractions

## GIN
check_dups <- create_check_dups("GIN")

## GMB
check_dups <- create_check_dups("GMB")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("GMB")

## GNB
check_dups <- create_check_dups("GNB")

## GNQ
check_dups <- create_check_dups("GNQ")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

## GUY
check_dups <- create_check_dups("GUY")

## HTI
check_dups <- create_check_dups("HTI")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% c("8809", "8796", "8779"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("HTI")

## IDN
check_dups <- create_check_dups("IDN")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
full_output <- full_output[!(Master_UID %in% c("10847", "23658"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("IDN")

## IND
check_dups <- create_check_dups("IND")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = NULL)
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% c("8938", "8983", "8984", "SR_A_727", "SR_A_732", "SR_A_731", "SR_A_736"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("IND")

## KEN
check_dups <- create_check_dups("KEN")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("11438", "SR_A_270", "SR_A_271"))] # Manually identified duplicates

## KHM
check_dups <- create_check_dups("KHM")
check_dups <- create_check_dups_poly("KHM")

## KIR
check_dups <- create_check_dups("KIR")
check_dups <- create_check_dups_poly("KIR")

## LAO
check_dups <- create_check_dups("LAO")
check_dups <- create_check_dups_poly("LAO")

## LBR
check_dups <- create_check_dups("LBR")

## LKA
check_dups <- create_check_dups("LKA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("15838", "19774", "19329"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("LKA")

## MDG
check_dups <- create_check_dups("MDG")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("11866", "11864", "11796", "11567", "11809", "11794", "11814", "11569", "11820", "11857", "11793", "11791", "11792", "11666", "11665", "11662", "11645", "11631", "11628", "11556"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("MDG")
# Deal with GAHI:ESPEN duplicate mapping
full_output[source == "GAHI" & data_collect_method == "Mapping" & country == "MDG", year := 2005]
check_dups <- create_check_dups("MDG")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
# Remove ESPEN duplicates (old vs. new API downloads)
check_dups <- create_check_dups_same_source("MDG")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

## MDV
check_dups <- create_check_dups("MDV")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"))
check_dups <- create_check_dups_poly("MDV")

## MHL
check_dups <- create_check_dups("MHL")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"))
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("13461"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("MHL")

## MLI
check_dups <- create_check_dups("MLI")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("MLI")

## MMR
check_dups <- create_check_dups("MMR")
check_dups <- create_check_dups_poly("MMR")

## MOZ
check_dups <- create_check_dups("MOZ")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("13668", "13672", "13667", "13670", "13671", "13674", "13549"))] # Manually identified duplicates

## MWI
check_dups <- create_check_dups("MWI")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("11872"))] # Manually identified duplicates
full_output <- full_output[!(Master_UID %in% c("temp_5911", "12027", "12039", "temp_5912", "12022", "temp_5892", "temp_5913", "12035", "12019", "temp_5898", "12026", "temp_5884", "12044", "temp_5886"))] # Deduplicate mapping data (differ in years across sources but identical coords, cases, and N; prefer most recent ESPEN source)
check_dups <- create_check_dups_poly("MWI")

## MYS
check_dups <- create_check_dups("MYS")
check_dups <- create_check_dups_poly("MYS")

## NCL
check_dups <- create_check_dups("NCL")
check_dups <- create_check_dups_poly("NCL")

## NER
check_dups <- create_check_dups("NER")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("NER")
full_output <- full_output[!(Master_UID %in% c("14021", "14019", "14020", "14022", "14023"))] # Manually identified duplicates

## NGA
check_dups <- create_check_dups("NGA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("NGA")
full_output <- full_output[!(Master_UID %in% c("SR_A_256", "SR_A_255", "SR_A_254", "SR_A_253", "14476", "14459", "14081"))] # Manually identified duplicates

## NIU
check_dups <- create_check_dups("NIU")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("NIU")
full_output <- full_output[!(Master_UID %in% c("DOSS_A_21", "SR_B_296; SR_B_297; SR_B_298", "DOSS_A_20", "SR_B_299", "SR_B_301", "SR_B_300"))] # Manually identified duplicates

## NPL
check_dups <- create_check_dups("NPL")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "WHO")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("NPL")
full_output <- full_output[!(Master_UID %in% c("13717", "13716", "13714", "13713", "13712", "13789", "13697", "13696"))] # Manually identified duplicates

## PHL
check_dups <- create_check_dups("PHL")
check_dups <- create_check_dups_poly("PHL")

## PLW
check_dups <- create_check_dups_poly("PLW")

## PNG
check_dups <- create_check_dups("PNG")
check_dups <- create_check_dups_poly("PNG")

## PYF
check_dups <- create_check_dups("PYF")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"))
check_dups <- create_check_dups_poly("PYF")

## RWA
check_dups <- create_check_dups("RWA")

## SEN
check_dups <- create_check_dups("SEN")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("SEN")

## SLE
check_dups <- create_check_dups("SLE")
check_dups <- create_check_dups_poly("SLE")

## SSD
check_dups <- create_check_dups("SSD")

## STP
check_dups <- create_check_dups("STP")

## TCD
check_dups <- create_check_dups("TCD")

## TGO
check_dups <- create_check_dups("TGO")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
full_output <- full_output[!(stri_detect_fixed(Master_UID, "SR_A") & country == "TGO" & data_collect_method == "TAS")]
check_dups <- create_check_dups_poly("TGO")

## THA
check_dups <- create_check_dups("THA")
check_dups <- create_check_dups_poly("THA")

## TLS
check_dups <- create_check_dups("TLS")
full_output <- full_output[!(Master_UID %in% c("17789"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("TLS")

## TON
check_dups <- create_check_dups("TON")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"))
full_output <- full_output[!(Master_UID %in% c("SR_F_55", "SR_F_53"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("TON")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"))

## TUV
check_dups <- create_check_dups("TUV")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"))

## TZA
check_dups <- create_check_dups("TZA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("16187", "16152", "16542", "16207", "16469", "16541", "16685", "temp_8103"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("TZA")

## UGA
check_dups <- create_check_dups("UGA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("18302", "18386", "18289", "18374", "18530"))] # Manually identified duplicates
check_dups <- create_check_dups_rounded_1("UGA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_poly("UGA")

## VNM
check_dups <- create_check_dups_poly("VNM")

## VUT
check_dups <- create_check_dups("VUT")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("18658", "18663", "18682", "18628", "18631", "18662", "18627", "18626", "18653", "18662", "18661", "18690", "18634", "18633", "SR_A_447"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("VUT")

## WLF
check_dups <- create_check_dups("WLF")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
full_output <- full_output[!(Master_UID %in% c("SR_B_293", "SR_B_292"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("WLF")

## WSM
check_dups <- create_check_dups("WSM")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
full_output <- full_output[!(Master_UID %in% c("15350", "15420"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("WSM")

## YEM
check_dups <- create_check_dups("YEM")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
full_output <- full_output[!(Master_UID %in% c("SR_C_488", "SR_C_485", "SR_C_472", "18775"))] # Manually identified duplicates
check_dups <- create_check_dups_poly("YEM")

## ZMB
check_dups <- create_check_dups("ZMB")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "sys_rev")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
full_output <- full_output[!(Master_UID %in% c("18891", "18893", "18902", "18901", "temp_9629", "temp_9721", "temp_9658", "temp_9722", "temp_9735", "temp_9737", "temp_9670", "temp_9736", "temp_9708", "temp_9628", "temp_9657", "temp_9707", "temp_9697", "temp_9650", "temp_9728", "temp_9631", "temp_9826", "temp_9681"))] # Manually identified duplicates
# check_dups <- create_check_dups_poly("ZMB")
check_dups <- create_check_dups_full_coords_no_year("ZMB")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]
check_dups <- create_check_dups_same_source_no_year("ZMB")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

## ZWE
check_dups <- create_check_dups_same_source_no_year("ZWE")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_lf_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "data_collect_method", "N", "had_lf_poly", "lf_prev", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

table(full_output$country)

####### Final cleanups
full_output[point == 1 & is.na(latitude) & !is.na(shapefile) & !is.na(location_code), point := 0]
full_output <- full_output[!(point == 1 & is.na(latitude) & country %in% c("YEM", "GHA", "LKA", "NGA", "IND", "TZA"))]
full_output <- full_output[!is.na(N) & !is.na(had_lf_poly)]
full_output <- full_output[N > 0]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Perform diagnostic & age adjustments ###################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Clean up age_start and age_end
# Set missing ages according to programmatic stage, when available
full_output[is.na(age_start) & (data_collect_method == "Mapping"), age_start := 15]
full_output[is.na(age_end) & (data_collect_method == "Mapping"), age_end := 94]
full_output[is.na(age_start) & (data_collect_method == "SS/Spot") & (year < 2011), age_start := 2]
full_output[is.na(age_start) & (data_collect_method == "SS/Spot") & (year >= 2011), age_start := 5]
full_output[is.na(age_end) & (data_collect_method == "SS/Spot"), age_end := 94]
full_output[is.na(age_start) & (data_collect_method == "TAS"), age_start := 6]
full_output[is.na(age_end) & (data_collect_method == "TAS"), age_end := 7]

full_output[age_end > 94, age_end := 94] # Cap population at 94 due to availability of single-year GBD estimates only to that age
full_output[age_group == "Adults/Children" & age_start == 0 & age_end == 0, age_end := 94] # Set adult/children to all ages if no age range is provided
full_output[age_group == "Adults" & age_start == 0 & age_end == 0, age_start := 15] # Set adults to 15-94 if no age range is provided
full_output[age_group == "Adults" & age_end == 0, age_end := 94] # Set adults to 15-94 if no age range is provided
full_output[age_group == "Unknown" & age_start == 0 & age_end == 0, age_end := 94] # If age group is unknown, assume all ages
full_output[age_group == "Children" & age_start == 0 & age_end == 0, age_end := 14] # If age group is children with no given end age, assume 0-14
full_output[age_group == "Adults/Children" & is.na(age_start) & is.na(age_end), age_start := 0] # If adults/children and no age or start age, assume 0-94
full_output[age_group == "Adults/Children" & is.na(age_end), age_end := 94] # If adults/children and no age or start age, assume 0-94
full_output[is.na(age_start), age_start := 0] # If no age_start provided, assume 0
full_output[is.na(age_end), age_end := 94] # If no age_end provided, assume 94
full_output[age_end == 0, age_end := 94] # If age_end is given as zero, assume 94

# Be sure that end age is >= start_age; these studies just had their age_end and age_start swapped
temp_age_data <- full_output[age_end < age_start, ]
temp_end <- temp_age_data$age_end
temp_start <- temp_age_data$age_start
full_output[age_end < age_start, "age_end"] <- temp_start
full_output[age_end < age_start, "age_start"] <- temp_end

###### Create a new unique ID
full_output$new_unique_id <- 1:nrow(full_output)

###### Subset years
full_output <- full_output[year %in% 1990:2018]

if (perform_crosswalk_new) { ###### Perform age and diagnostic crosswalks
  ### Setup
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(MASS, boot, data.table, ggplot2)
  path <- paste0(<<<< FILEPATH REDACTED >>>>)
  library(fda, lib.loc = path)
  library(faraway, lib.loc = path)
  library(subplex, lib.loc = <<<< FILEPATH REDACTED >>>>)
  library(BayesianTools, lib.loc = <<<< FILEPATH REDACTED >>>>)
  library(matrixStats, lib.loc = path)

  ### Load shared functions
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)

  ### Define needed functions
  ## Calculate the prevalence-age curve in logit space; can put something other than splines in here
  applySpline_dx <- function(vec, Basis, avec) {
    AgeSpline <- fd(vec, Basis)
    return(predict(AgeSpline, avec))
  }

  ## Apply constraints to shape of spline
  applyBasisConstraints_dx <- function(betas) {
    alphas <- betas
    alphas[1] <- betas[1]
    return(alphas)
  }

  ## Calculate prevalence
  calculatePrevalence_dx <- function(AgeSpline) {
    return((ilogit(AgeSpline)))
  }

  ## Calculate the prevalence-age curve in logit space; can put something other than splines in here
  applySpline <- function(vec, Basis, avec) {
    AgeSpline <- fd(vec, Basis)
    return(predict(AgeSpline, avec))
  }

  ## Apply constraints to shape of spline
  applyBasisConstraints <- function(betas, diag) {
    alphas <- betas
    alphas[1] <- betas[1]

    if (diag == "mf") {
      for (a in 2:2) {
        alphas[a] <- betas[1] + sum(exp(betas[2:a]))
      }
    }

    return(alphas)
  }

  ## Calculate prevalence
  calculatePrevalence <- function(AgeSpline) {
    return((ilogit(logit(0.00001) + AgeSpline)))
  }

  crosswalk_data <- full_output

  country_years <- unique(crosswalk_data[, c("country", "year")])

  age_ids <- c(49:142) # Age 1 = ID 49, Age 94 = ID 142; get_population does not return single_year_age for ages 95-99
  gbd_round_id <- 6 ### 2019

  location_hierarchy <- get_location_metadata(location_set_id = 35, gbd_round_id = gbd_round_id)
  age_metadata <- get_age_metadata(age_group_set_id = 12, gbd_round_id = gbd_round_id)

  country_years <- merge(country_years, location_hierarchy[, c("ihme_loc_id", "location_id")], by.x = "country", by.y = "ihme_loc_id", all.x = TRUE)
  country_years[country %in% c("NCL", "TUV", "PYF", "NRU", "NIU", "WLF"), location_id := 21] ### Use population distribution for Oceania, as GBD doesn't break out these countries separately

  popNumbers <- data.table(get_population(age_group_id = age_ids, location_id = unique(country_years$location_id), year_id = unique(country_years$year), sex_id = 3, single_year_age = T, gbd_round_id = gbd_round_id, decomp_step = "iterative"))
  popNumbers_under_1 <- data.table(get_population(age_group_id = 28, location_id = unique(country_years$location_id), year_id = unique(country_years$year), sex_id = 3, single_year_age = F, gbd_round_id = gbd_round_id, decomp_step = "iterative"))
  popNumbers_concat <- rbind(popNumbers, popNumbers_under_1)

  popNumbers_concat$population <- ceiling(popNumbers_concat$population)

  ## Get age IDs and combine with age group names
  ages <- get_ids("age_group")
  pop_age_merged <- merge(popNumbers_concat, ages, by = "age_group_id")

  ## Calculate total population (0-94 years of age)
  pop_agg <- aggregate(population ~ location_id + year_id, data = pop_age_merged, sum)
  colnames(pop_agg)[3] <- "total_population"
  pop_age_merged <- merge(pop_age_merged, pop_agg, by = c("location_id", "year_id"))

  ## Calculate % of pop. in each age year
  pop_age_merged$age_prop <- pop_age_merged$population / pop_age_merged$total_population

  ## Convert to wide format
  pop_age_merged_wide <- reshape(subset(pop_age_merged, select = -c(age_group_id, sex_id, run_id, population, total_population)), idvar = c("location_id", "year_id"), timevar = "age_group_name", direction = "wide")
  colnames(pop_age_merged_wide)[2] <- "year"

  ## Add location_id to crosswalk_data
  crosswalk_data <- merge(crosswalk_data, unique(country_years[, c("country", "location_id")]), by = "country", all.x = TRUE)
  crosswalk_data$location_id <- as.integer(crosswalk_data$location_id)

  ## Merge age distribution with crosswalk data
  newData <- merge(crosswalk_data, pop_age_merged_wide, by = c("location_id", "year"), all.x = TRUE)
  names(newData)[names(newData) == "age_prop.<1 year"] <- "age_prop.0" # Rename age_prop.<1 year

  # Retrieve crosswalk model objects
  load(<<<< FILEPATH REDACTED >>>>)
  load(<<<< FILEPATH REDACTED >>>>)
  load(<<<< FILEPATH REDACTED >>>>)

  age_BREAKS_mf <- c(3, 6, 15, 34, 65)

  freq_coefs <- dx_crosswalk_model$Optim$par

  if (diag %in% c("ict", "ICT")) {
    data_to_adjust <- newData[diagnostic == "adjusted_mf"]
  } else if (diag == "mf") {
    data_to_adjust <- newData[diagnostic %in% c("ICT", "FTS")]
  }

  # Setup
  SplineNum <- length(dx_crosswalk_model$BREAKS) + 2
  BREAKS <- dx_crosswalk_model$BREAKS
  Basis <- create.bspline.basis(rangeval = c(0, 94.5), breaks = BREAKS)
  avec <- dx_crosswalk_model$avec
  
  age_threshold <- 64.5
  age_threshold_lower <- 5.5

  #### Retrieve parameter values
  alphas <- applyBasisConstraints_dx(freq_coefs[1:SplineNum])

  #### Calculate prevalences
  AgeSpline <- applySpline_dx(alphas, Basis, avec)

  AgeLoc <- data.table("age" = avec, "prev" = calculatePrevalence_dx(AgeSpline))
  colnames(AgeLoc)[2] <- "prev"
  AgeLoc[age > age_threshold, prev := AgeLoc[age == age_threshold, prev]]
  AgeLoc[age < age_threshold_lower, prev := AgeLoc[age == age_threshold_lower, prev]]

  AgeLoc <- AgeLoc[age > 0]

  for (i in 1:nrow(data_to_adjust)) {
    message(paste0(i, " of ", nrow(data_to_adjust)))
    Start <- data_to_adjust[i, age_start]
    End <- data_to_adjust[i, age_end]
    prev <- data_to_adjust[i, had_lf_poly / N]

    AgeLocTemp <- AgeLoc[(age - 0.5) %in% Start:End, ]
    prev.mf.logit <- logit(prev)
    if (prev.mf.logit == "-Inf") {
      prev.mf.logit <- logit(0.001)
    } else if (prev.mf.logit == "Inf") {
      prev.mf.logit <- logit(0.999)
    }

    if (diag %in% c("ict", "ICT")) {
      AgeLocTemp$prev_mod <- inv.logit(logit(AgeLocTemp$prev) + prev.mf.logit)
    } else if (diag == "mf") {
      AgeLocTemp$prev_mod <- inv.logit(prev.mf.logit - logit(AgeLocTemp$prev))
    }

    PGroup <- data_to_adjust[i, which(colnames(data_to_adjust) == paste0("age_prop.", Start)):which(colnames(data_to_adjust) == paste0("age_prop.", End))]
    crosswalked_prev <- sum(AgeLocTemp[(age - 0.5) %in% Start:End, prev_mod] * PGroup / sum(PGroup))
    if (data_to_adjust[i, had_lf_poly] == 0) {
      crosswalked_prev <- 0
    }
    dx_adjusted <- crosswalked_prev * data_to_adjust[i, N]

    data_to_adjust[i, "had_lf_poly_dx_adjusted"] <- dx_adjusted
  }

  ###### Now perform age crosswalk
  newData <- merge(newData, data_to_adjust[, c("new_unique_id", "had_lf_poly_dx_adjusted")], by = "new_unique_id", all.x = TRUE)
  newData[is.na(had_lf_poly_dx_adjusted), had_lf_poly_dx_adjusted := had_lf_poly]

  new_data_to_adjust <- newData[(diagnostic %in% c("adjusted_mf", "ICT", "FTS")) & (!(age_start == 0 & age_end == 94))]

  # Setup
  if (diag %in% c("ict", "ICT")) {
    freq_coefs <- ICT_crosswalk_model$Optim$par
    SplineNum <- length(ICT_crosswalk_model$BREAKS) + 2
    BREAKS <- ICT_crosswalk_model$BREAKS
    Basis <- create.bspline.basis(rangeval = c(0, 94.5), breaks = BREAKS)
    avec <- ICT_crosswalk_model$avec
    age_threshold <- 64.5
  } else if (diag == "mf") {
    freq_coefs <- freq_coefs_mf
    BREAKS <- age_BREAKS_mf
    SplineNum <- length(BREAKS) + 2
    Basis <- create.bspline.basis(rangeval = c(0, 94.5), breaks = BREAKS)
    age_threshold <- 94.5
  }
  
  if (nrow(new_data_to_adjust) > 0) {
    for (i in 1:nrow(new_data_to_adjust)) {
      message(paste0(i, " of ", nrow(new_data_to_adjust)))

      current <- new_data_to_adjust[i,]

      #### Retrieve parameter values
      alphas <- applyBasisConstraints(freq_coefs[1:SplineNum], diag)

      #### Calculate prevalences
      AgeSpline <- applySpline(alphas, Basis, avec)

      AgeLoc <- data.table("age" = avec, "prev" = calculatePrevalence(AgeSpline))
      colnames(AgeLoc)[2] <- "prev"

      #### Constrain prevalence to be flat beyond a threshold age
      AgeLoc[age > age_threshold, prev := AgeLoc[age == age_threshold, prev]]

      Start <- current$age_start
      End <- current$age_end
      prev <- current[, had_lf_poly_dx_adjusted / N]

      # Calculate baseline prevalence in sampled age range
      PGroup <- current[, which(colnames(current) == paste0("age_prop.", Start)):which(colnames(current) == paste0("age_prop.", End))]
      baseline_prev_sample_range <- sum(PGroup * AgeLoc[(age - 0.5) %in% Start:End, prev]) / sum(PGroup)

      # Calculate baseline prevalence in all ages (0-94)
      PGroup_all <- current[, which(colnames(current) == paste0("age_prop.", 0)):which(colnames(current) == paste0("age_prop.", 94))]
      baseline_prev_all <- sum(PGroup_all * AgeLoc$prev) / sum(PGroup_all)

      # Calculate scaling factor and adjusted all-age prevalence
      scaling <- logit(prev) - logit(baseline_prev_sample_range)
      adjusted_prev_all_age <- inv.logit(logit(baseline_prev_all) + scaling)
      age_adjusted <- adjusted_prev_all_age * current$N

      new_data_to_adjust[i, "had_lf_poly_dx_age_adjusted"] <- age_adjusted
    }
  }

  newData <- merge(newData, new_data_to_adjust[, c("new_unique_id", "had_lf_poly_dx_age_adjusted")], by = "new_unique_id", all.x = TRUE)
  newData[is.na(had_lf_poly_dx_age_adjusted), had_lf_poly_dx_age_adjusted := had_lf_poly_dx_adjusted]

  full_output <- copy(newData[, grep("age_prop.", colnames(newData), invert=TRUE), with=FALSE])
  full_output$had_lf_poly_original <- full_output$had_lf_poly
  full_output$had_lf_poly <- as.integer(round(full_output$had_lf_poly_dx_age_adjusted, 0))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Saving outputs #########################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Final processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Miscellaneous post hoc corrections
full_output[country == "BEN" & nid == "431787" & latitude == "11.35543", point := 1]
full_output <- full_output[!(point == 0 & is.na(location_code))]

###### Drop new id col
full_output <- unique(full_output[, -c("new_unique_id", "resamp_test_case", "latitude_rounded", "latitude_rounded_1", "longitude_rounded", "longitude_rounded_1")]) # drop completely duplicate rows (excluding new_unique_id)

### Additional ad hoc data cleaning
full_output[shapefile == "lf_Tas_1323", shapefile := "lf_tas_1323"]
full_output[shapefile == "lf_Tas_1325", shapefile := "lf_tas_1325"]

full_output <- full_output[!(nid %in% c(409420))] # Deduplication
full_output <- full_output[!(Master_UID %in% c("TAS_2018_1349", "TAS_2018_1344", "TAS_2018_1345"))] # Outlier
full_output <- full_output[!(nid == "147731" & source == "sys_rev")] # Deduplication

if (save_outputs) {
  full_output_backup <- full_output
  full_output <- full_output[diagnostic %in% c("ICT", "adjusted_mf", "FTS", "br")]
  
  full_output[, lf_prev := as.numeric(round(had_lf_poly, 0) / round(N, 0))]
  write.csv(full_output, paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
}
