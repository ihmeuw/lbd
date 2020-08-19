#####################################################################
## Processing and combining various data sources of LF prevalence data
##
#####################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Step 1: Script Prep ################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Set directory and load libraries
rm(list = ls())
pacman::p_load(
  data.table,
  ggplot2,
  magrittr,
  stringr,
  plyr,
  Hmisc,
  dplyr
)

# string to add to filename for exporting
pathadd <- <<<< FILEPATH REDACTED >>>>

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Step 2: Clean geo-referencing ######################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (exists("loc_meta") == F) { # get central location database for ISO
  loc_meta <- fread(<<<< FILEPATH REDACTED >>>>) %>% as.data.table()
  loc_meta <- loc_meta[, c("short_name", "iso3")]
  setnames(loc_meta, c("short_name", "iso3"), c("location_name_short", "ihme_loc_id"))
}
locations <- copy(loc_meta)

# separate out rows with problematic shapefile references
check_for_missing_shapefiles <- function(df) {
  # Pull a list of shapefiles for the polygon data for the region in question
  shapefiles <- unique(df$shapefile) %>% as.character()

  ## Check that all shapefile entries are real shapefiles and report bad entries
  real_shapefiles <- gsub(".shp", "", grep(".shp", value = T, list.files(<<<< FILEPATH REDACTED >>>>)))

  not_real_shapefiles <- shapefiles[!(shapefiles %in% real_shapefiles) & is.na(shapefiles) == F & shapefiles != ""]

  shapefiles <- shapefiles[(shapefiles %in% real_shapefiles)]

  return(list(
    "shapefiles" = shapefiles,
    "not_real_shapefiles" = not_real_shapefiles
  ))
}

if (FALSE) {
  geo_ref <- fread(<<<< FILEPATH REDACTED >>>>)

  ## drop unnecessary columns to cut comp time
  geo_ref <- geo_ref[, 1:grep("_merge", names(geo_ref))]

  # clean up points that have polygon reference data entered (this can cause issues with coverage plot code and other code)
  geo_ref[shapefile == "", shapefile := NA]

  # keep polygon information for reference in notes field
  geo_ref[point == 1 & is.na(shapefile) == F, notes := paste0(notes, "/", shapefile, "/", location_code)]
  geo_ref[point == 1 & is.na(shapefile) == F, shapefile := NA]
  geo_ref[point == 1 & is.na(shapefile) == F, location_code := NA]

  # fix known polygon reference problems
  fix_poly_ref <- fread(<<<< FILEPATH REDACTED >>>>)
  for (i in 1:length(fix_poly_ref$old_shapefile)) {
    geo_ref[shapefile == fix_poly_ref$old_shapefile[i], shapefile := fix_poly_ref$correct_shapefile[i]]
  }

  shapefile_check <- check_for_missing_shapefiles(geo_ref)
  bad_shapefiles <- shapefile_check[["not_real_shapefiles"]]
  check_shapefiles <- geo_ref[shapefile %in% bad_shapefiles]
  check_shapefiles <- rbind(check_shapefiles, geo_ref[is.na(shapefile) == F & shapefile != "" & (is.na(location_code) == T | location_code == "")])

  geo_ref <- geo_ref[!(shapefile %in% bad_shapefiles)]
  geo_ref <- geo_ref[!(is.na(shapefile) == F & shapefile != "" & (is.na(location_code) == T | location_code == ""))]
  fix_sf_loc <- fread(<<<< FILEPATH REDACTED >>>>)
  geo_ref <- rbind(geo_ref, fix_sf_loc)
  geo_ref[shapefile == "lf_TAS_Mali_2015_BANKASS_BANDIAGARA_KORO\r\n", shapefile := "lf_TAS_Mali_2015_BANKASS_BANDIAGARA_KORO"]
  geo_ref[shapefile == "lf_TAS_Mali_2015_KOULIKORO_BANAMBA_KANGABA\r\n", shapefile := "lf_TAS_Mali_2015_KOULIKORO_BANAMBA_KANGABA"]
  write.csv(geo_ref, <<<< FILEPATH REDACTED >>>>, row.names = F)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Step 3: Loading in current version of master doc ###################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

geo_ref <- fread(<<<< FILEPATH REDACTED >>>>)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Step 4: Standardize fields & drop duplicates and rows with unresolvable missingness ################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

shapefile_check <- check_for_missing_shapefiles(geo_ref)
bad_shapefiles <- shapefile_check[["not_real_shapefiles"]]
nrow(bad_shapefiles)
check_shapefiles <- geo_ref[shapefile %in% bad_shapefiles]

# get ISO from adm0
geo_ref$adm0 <- tolower(geo_ref$adm0)
locations$location_name_short <- tolower(locations$location_name_short)
change_names <- unique(geo_ref$adm0)[(unique(geo_ref$adm0) %in% locations$location_name_short) == F] %>% sort()
if (length(change_names) > 0) {
  # set up vector of corresponding standard names
  change_names_std <- c("central african republic", "c??te d'ivoire", "micronesia", "sao tome and principe", "wallis & futuna")
  for (i in 1:length(change_names)) {
    geo_ref[adm0 == change_names[i], adm0 := change_names_std[i]]
  }
}

setnames(locations, "location_name_short", "adm0")
geo_ref <- merge(geo_ref, locations, by = "adm0", all.x = T)

# get all geo-referenced place-names
geo_ref_locs <- geo_ref[point == 1 & (is.na(latitude) == F | is.na(longitude) == F | latitude != "" | longitude != ""), .(adm0, adm1, adm2, adm3, adm4, latitude, longitude)] %>% unique()
setnames(geo_ref_locs, c("adm0", "adm1", "adm3"), c("adm0_geo_ref", "adm1_geo_ref", "adm3_geo_ref"))

# overall cleaning
geo_ref[is.na(point) == T] %>% nrow()
dropped <- geo_ref[is.na(point) == T] %>% .[, drop_reason := "Not classified point/polygon"]
geo_ref <- geo_ref[is.na(point) == F]

geo_ref$data_collect_method <- tolower(geo_ref$data_collect_method)
table(geo_ref$data_collect_method)
geo_ref[grep("check", geo_ref$data_collect_method), data_collect_method := "SS/Spot"]
geo_ref[grep("site", geo_ref$data_collect_method), data_collect_method := "SS/Spot"]
geo_ref[unique(c(grep("lqas", geo_ref$data_collect_method), grep("mapping", geo_ref$data_collect_method))), data_collect_method := "Mapping"]
geo_ref[data_collect_method %in% c("c survey", "village surveillance", "other", "cl_tas", "operational research", "surveillance", "census", "post-endemic", "special survey", "post survey", "research"), data_collect_method := "Other"]
geo_ref[data_collect_method %in% c("ss", "sentinal site", "ss/spot", "ss/sc", "sc"), data_collect_method := "SS/Spot"]
geo_ref[data_collect_method %in% c("mapping", "baseline"), data_collect_method := "Mapping"]
geo_ref[grep("tas", geo_ref$data_collect_method), data_collect_method := "TAS"]
geo_ref[data_collect_method %in% c(
  "", "5-year programme to control bancroftian filariasis",
  "analyze 2007-2008 data to also conform to the 2011 antigen guidelines for the 6 to 7 yr old age group"
), data_collect_method := NA]
table(geo_ref$data_collect_method)

# clean/impute year_survey
geo_ref[is.na(year_survey) == T & year_end == 0 & is.na(year_start) == F & year_start != 0, year_survey := year_start]
geo_ref[is.na(year_survey) == T & is.na(year_end) == F & is.na(year_start) == F & year_start != 0 & year_end != 0, year_survey := round((year_start + year_end) / 2)]
geo_ref[is.na(year_survey) == T & is.na(year_start) == F & year_start != 0, year_survey := year_start]

geo_ref[(is.na(year_survey) == T | year_survey == 0 | year_survey == "") & reference_year != "" & is.na(reference_year) == F, year_survey := as.integer(reference_year) - 2]

### number being dropped due to missing year
table(is.na(geo_ref$year_survey) | geo_ref$year_survey == 0)
dropped <- rbind(dropped, geo_ref[is.na(year_survey) == T | year_survey == 0] %>% .[, drop_reason := "year unspecified"], fill = T)
geo_ref <- geo_ref[is.na(year_survey) == F & year_survey != 0]

## drop if WHO_DUP or EPI_REF_DUP
geo_ref[(WHO_DUP %in% c(1, 2)) | (EPI_REF_DUP %in% c(1)), ] %>% nrow()
dropped <- rbind(dropped, geo_ref[(WHO_DUP %in% c(1, 2)) | (EPI_REF_DUP %in% c(1)), ] %>% .[, drop_reason := "duplicated"])
geo_ref <- geo_ref[!(WHO_DUP %in% c(1, 2)) & !(EPI_REF_DUP %in% c(1)), ]

# drop out of range coordinates
geo_ref$latitude <- as.numeric(geo_ref$latitude)
geo_ref$longitude <- as.numeric(geo_ref$longitude)
geo_ref[(longitude > 180 | longitude < -180 | latitude > 180 | latitude < -180)] %>% nrow()
dropped <- rbind(dropped, geo_ref[point == 1 & is.na(longitude) == F & is.na(latitude) == F &
  (longitude > 180 | longitude < -180 | latitude > 180 | latitude < -180), ] %>% .[, ihme_loc_id := NULL] %>%
  .[, drop_reason := "coordinates out of bounds"], fill = T)
geo_ref <- geo_ref[!(point == 1 & is.na(longitude) == F & is.na(latitude) == F & (longitude > 180 | longitude < -180 | latitude > 180 | latitude < -180))]

# standardize data source field
geo_ref[grep("gahi", geo_ref$data_source), new_source := "GAHI"]
geo_ref[data_source == "sys_rev", new_source := "sys_rev"]
geo_ref[grep("WHO", geo_ref$data_source), new_source := "WHO"]
geo_ref[grep("Dossier", geo_ref$data_source), new_source := "WHO"]
geo_ref[data_source %in% c("Ben_TZ", "Brazil_MOH"), new_source := "country"]

# geo_review field
geo_ref[, geo_review := 1]

# sort out prevalences
geo_ref$pop_mf <- as.numeric(gsub(",", "", geo_ref$pop_mf)) %>% round()
geo_ref$pop_ict <- as.numeric(gsub(",", "", geo_ref$pop_ict)) %>% round()
geo_ref$np_ict <- as.numeric(gsub(",", "", geo_ref$np_ict)) %>% round()
geo_ref$np_mf <- as.numeric(gsub(",", "", geo_ref$np_mf)) %>% round()
geo_ref$prev_mf <- as.double(gsub(",", ".", geo_ref$prev_mf))
geo_ref$prev_ict <- as.double(gsub(",", ".", geo_ref$prev_ict))
geo_ref$pop_test_strip <- as.numeric(geo_ref$pop_test_strip)

geo_ref[is.na(prev_ict) == T & pop_ict != 0 & is.na(np_ict) == F, prev_ict := np_ict / pop_ict]
geo_ref[is.na(prev_mf) == T & pop_mf != 0 & is.na(np_mf) == F, prev_mf := np_mf / pop_mf]
geo_ref[is.na(prev_test_strip) == T & pop_test_strip != 0 & is.na(np_test_strip) == F, prev_test_strip := np_test_strip / pop_test_strip]

geo_ref[is.na(np_ict) == T & is.na(prev_ict) == F & prev_ict < 1 & is.na(pop_ict) == F, np_ict := prev_ict * pop_ict]
geo_ref[is.na(np_mf) == T & is.na(prev_mf) == F & prev_mf < 1 & is.na(pop_mf) == F, np_mf := prev_mf * pop_mf]
geo_ref[is.na(np_test_strip) == T & is.na(prev_test_strip) == F & prev_test_strip < 1 & is.na(pop_test_strip) == F, np_test_strip := prev_test_strip * pop_test_strip]

geo_ref[is.na(np_ict) == T & is.na(prev_ict) == F & prev_ict >= 1 & is.na(pop_ict) == F, np_ict := (prev_ict / 100) * pop_ict]
geo_ref[is.na(np_mf) == T & is.na(prev_mf) == F & prev_mf >= 1 & is.na(pop_mf) == F, np_mf := (prev_mf / 100) * pop_mf]
geo_ref[is.na(np_test_strip) == T & is.na(prev_test_strip) == F & prev_test_strip >= 1 & is.na(pop_test_strip) == F, np_test_strip := (prev_test_strip / 100) * pop_test_strip]

## output crosswalk analysis dataset
ict_mf <- geo_ref[is.na(prev_mf) == F & is.na(prev_ict) == F & pop_mf != 0 & pop_ict != 0, ]
ict_mf[, c("prev_ict", "prev_mf") := list(np_ict / pop_ict, np_mf / pop_mf)]

## generate analysis fields using hierarchy for diagnosis preference (ICT > MF > Test Strip > BR > ELISA > PCR)
geo_ref[
  is.na(prev_ict) == F & pop_ict != 0 & is.na(pop_ict) == F & pop_ict != "" &
    ((pop_ict > 20 & pop_ict > pop_mf) | (pop_mf != 0 & is.na(pop_mf) == F & pop_mf / pop_ict < 2) | is.na(pop_mf) == T),
  c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_ict, np_ict, prev_ict, "ict")
]
geo_ref[is.na(lf_prev) == T & pop_mf != 0 & is.na(pop_mf) == F & pop_mf != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_mf, np_mf, prev_mf, "mf")]
geo_ref[is.na(lf_prev) == T & pop_test_strip != 0 & is.na(pop_test_strip) == F & pop_test_strip != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_test_strip, np_test_strip, prev_test_strip, "test_strip")]
geo_ref[is.na(lf_prev) == T & pop_br != 0 & is.na(pop_br) == F & pop_br != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_br, np_br, prev_br, "br")]
geo_ref[is.na(lf_prev) == T & pop_elisa != 0 & is.na(pop_elisa) == F & pop_elisa != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_elisa, np_elisa, prev_elisa, "elisa")]
geo_ref[is.na(lf_prev) == T & pop_pcr != 0 & is.na(pop_pcr) == F & pop_pcr != "", c("N", "lf_pos", "lf_prev", "diagnostic") := list(pop_pcr, np_pcr, prev_pcr, "pcr")]

dropped <- rbind(dropped, geo_ref[is.na(N) == T, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "no valid sample size"], fill = T)

geo_ref <- geo_ref[is.na(N) == F]
dropped <- rbind(dropped, geo_ref[N < 20, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "N < 20"], fill = T)
geo_ref <- geo_ref[N > 20]
geo_ref$lf_pos <- geo_ref$lf_pos %>% as.integer()
dropped <- rbind(dropped, geo_ref[is.na(lf_pos) == T & is.na(lf_prev) == T, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "No numerator/prevalence"], fill = T)
geo_ref <- geo_ref[!(is.na(lf_pos) == T & is.na(lf_prev) == T)]

dropped <- rbind(dropped, geo_ref[lf_pos > N, ] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Illogical sample size"], fill = T)
geo_ref <- geo_ref[N >= lf_pos, ]
geo_ref[, lf_prev := lf_pos / N]

# change one polygon to point (known issue, check this if running again)
# fixing other discovered issues
geo_ref[point == 0 & (is.na(shapefile) == F | shapefile == "") & (location_code == "" | is.na(location_code) == T), point := 1]
geo_ref[location_code == "BINAH", location_code := 1]
geo_ref[shapefile == "lf_g2015_2001_0", shapefile := "lf_g2015_2010_0"]
geo_ref[shapefile == "lf_g2015_2010_0" & adm0 %in% c("togo", "vietnam") & adm1 != "", shapefile := "lf_g2015_2011_2"]
dropped <- rbind(dropped, geo_ref[grep(123.6693, geo_ref$latitude)] %>% .[, drop_reason := "coordinates out of bounds"], fill = T)
geo_ref <- geo_ref[!(grep(123.6693, geo_ref$latitude))]

### drop non-georeferenced data ###############################################
dropped <- rbind(dropped, geo_ref[point == 0 & (is.na(shapefile) == T | shapefile == "")] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Polygon not geolocated"], fill = T)
geo_ref <- geo_ref[!(point == 0 & (is.na(shapefile) == T | shapefile == ""))]

# check for presence of geo-referencing in other datasets
if (FALSE) {
  gahi <- fread(<<<< FILEPATH REDACTED >>>>)
  gahi_locs <- gahi[is.na(Latitude) == F & is.na(Longitude) == F, .(Latitude, Longitude, ADM0, ADM1, ADM2, Site_name)] %>% unique()
  gahi_locs$ADM0 <- tolower(gahi_locs$ADM0)
  gahi_locs$ADM1 <- tolower(gahi_locs$ADM1)
  gahi_locs$ADM2 <- tolower(gahi_locs$ADM2)
  gahi_locs$Site_name <- tolower(gahi_locs$Site_name)
  setnames(gahi_locs, c("ADM0", "ADM1", "ADM2", "Site_name"), c("adm0_gahi", "adm1_gahi", "adm2", "adm4"))

  m_locs <- geo_ref[point == 1 & (is.na(latitude) == T | is.na(longitude) == T | latitude == "" | longitude == ""), .(Master_UID, adm0, adm1, adm2, adm3, adm4)]
  m_locs <- apply(m_locs, 2, tolower) %>% as.data.table()
  check_merge <- merge(m_locs, gahi_locs, by = c("adm2", "adm4")) %>%
    as.data.table() %>%
    .[adm2 != "", ]

  if (nrow(check_merge) > 0) {
    mergeable <- geo_ref[Master_UID %in% as.integer(check_merge$Master_UID), ]
    check_merge$Master_UID <- as.integer(check_merge$Master_UID)
    mergeable <- merge(mergeable, check_merge[, .(Master_UID, Latitude, Longitude)], by = "Master_UID")
    mergeable[, latitude := Latitude]
    mergeable[, longitude := Longitude]
    mergeable[, c("Latitude", "Longitude") := NULL]

    geo_ref <- geo_ref[!(Master_UID %in% as.integer(check_merge$Master_UID)), ]
    geo_ref <- rbind(geo_ref, mergeable, use.names = T)
  }
}

dropped <- rbind(dropped, geo_ref[point == 1 & (is.na(latitude) == T | is.na(longitude) == T | latitude == "" | longitude == "")] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Point not geolocated"], fill = T)
geo_ref <- geo_ref[!(point == 1 & (is.na(latitude) == T | is.na(longitude) == T | latitude == "" | longitude == ""))]

# drop outliered data

dropped <- rbind(dropped, geo_ref[outlier == 1] %>% .[, ihme_loc_id := NULL] %>% .[, drop_reason := "Outliered"], fill = T)
geo_ref <- geo_ref[is.na(outlier) == T]

# tag test polygons
geo_ref[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Kenema", "Kailahun", "Kono"), location_code := 2654]
geo_ref[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Bombali", "Koinadugu", "Port Loko", "Tonkolili", "Kambia"), location_code := 2655]
geo_ref[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Pujehun", "Bonthe", "Moyamba", "Bo"), location_code := 2656]
geo_ref[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & adm1 %in% c("Western Rural", "Freetown"), location_code := 2657]
geo_ref[adm0 == "sierra leone" & year_survey == 2005 & data_collect_method == "Mapping" & location_code == 2655, resamp_test_case := 3]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Step 5: Final subsetting and export ################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

geo_ref_final <- geo_ref[is.na(lf_prev) == F, .(
  Master_UID, new_source, year_survey, ihme_loc_id, point, latitude, longitude, shapefile, location_code, N,
  lf_pos, lf_prev, diagnostic, data_collect_method, sampling_prod, resamp_test_case, age_start, age_end, age_group
)]
setnames(geo_ref_final, c("new_source", "year_survey", "ihme_loc_id"), c("source", "original_year", "country"))
geo_ref_final[, c("weight", "cluster_id", "year") := list(1, 1, original_year)]

# currently imputating population based survey design for unknowns (since these would generally be TAS)
geo_ref_final[sampling_prod == 1 | sampling_prod == 9, sampling := 1]
geo_ref_final[sampling_prod == 0, sampling := 0]
geo_ref_final[point == 0 & is.na(sampling_prod) == T, sampling := 1]
geo_ref_final$sampling_prod <- NULL

setnames(geo_ref_final, c("lf_pos"), c("had_lf_poly"))
write.csv(geo_ref_final, paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
setnames(geo_ref_final, "had_lf_poly", "had_lf_points")
write.csv(geo_ref_final[point == 1], paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)

if (TRUE) { # deal with dropped rows
  drop_bar <- ggplot(data = dropped[drop_reason != "duplicated" & drop_reason != "Outliered"]) + geom_bar(aes(x = drop_reason, fill = drop_reason)) +
    coord_flip() + guides(fill = "none") + ylab("Count") + xlab("Reason") +
    ggtitle("Why were data dropped?") + theme_minimal(base_size = 20)

  drop_table <- table(dropped$drop_reason) %>% as.data.table()
  setnames(drop_table, "V1", "Reason")
  write.csv(drop_table, <<<< FILEPATH REDACTED >>>>)
  write.csv(dropped, <<<< FILEPATH REDACTED >>>>)
  png(filename = <<<< FILEPATH REDACTED >>>>, height = 1000, width = 1000)
  drop_bar
  dev.off()
}
