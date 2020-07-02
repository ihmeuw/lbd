
# Subset geomatched data functions for covariates (easier to use)
subset_geomatched_male_circumcision <- function(geo_data){
  
  # Drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, year, strata, psu, point, shapefile, location_code, latitude, longitude,
                              sex_id, age_year, int_year, int_month, male_circumcision, pweight, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # rename variables
  setnames(geo_data, c("survey_series"), c("source"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  geo_data[, longitude := as.numeric(longitude)]
  
  # Drop NID's related to the country specific Nigeria survey (too different from the DHS)
  # Looking into Nigeria
  drop_nigeria_nid <- c(325046, 151719, 324443)
  geo_data <- geo_data[!nid %in% drop_nigeria_nid]
  
  #Drop Timor-Leste NIDS until further Non-SSA vetting
  drop_tls_nids <- c(21274)
  geo_data <- geo_data[!nid %in% drop_tls_nids]
  
  #Drop DHS_SP survey matched to polygons(see docs)
  geo_data <- geo_data[!(nid == 21198)]
  
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
  
  # subset to ages 15-49 (and only surveys with this full range)
  geo_data <- geo_data[between(age_year, 15, 49)]
  drop <- geo_data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
  geo_data <- geo_data[!nid %in% drop]
  
  # For paper, see how many are missing MC status, location and survey weight
  pct_total_missing <- 
    geo_data %>%
    filter(is.na(male_circumcision) | is.na(pweight) | is.na(point) |
             (point == 1 & (is.na(latitude) | is.na(longitude))) | 
             (point == 0 & (is.na(shapefile) | is.na(location_code)))) %>% 
    nrow() %>% 
    divide_by(nrow(geo_data)) %>% 
    multiply_by(100) %>% round(1)
  
  # drop observations with missing circumcision response or individual weight
  pct_missing_circ <- 
    geo_data %>% 
    filter(is.na(male_circumcision) | is.na(pweight)) %>% 
    nrow() %>% 
    divide_by(nrow(geo_data)) %>% 
    multiply_by(100) %>% round(1)
  message(paste0(pct_missing_circ, "% of data are missing male circumcision indicator or weight response"))
  geo_data <- geo_data[!is.na(male_circumcision) & !is.na(pweight)]
  
  # Drop data missing information on point or polygon 
  pct_missing <- 
    geo_data %>% 
    filter(is.na(point)) %>% 
    nrow() %>% 
    divide_by(nrow(geo_data)) %>% 
    multiply_by(100) %>% round(2)
  message(paste0(pct_missing, "% of data are missing point indicator"))
  geo_data <- geo_data[!is.na(point)]
  
  # drop point data with missing latitude/longitude
  pct_missing <- 
    geo_data %>% 
    filter(point == 1, is.na(latitude) | is.na(longitude)) %>% 
    nrow() %>% 
    divide_by(nrow(geo_data %>% filter(point == 1))) %>% 
    multiply_by(100) %>% round(2)
  message(paste0(pct_missing, "% of point data are missing either latitude or longitude and are being dropped"))
  
  # drop polygon data with missing shapefile/location code
  pct_missing <- 
    geo_data %>% 
    filter(point == 0, is.na(shapefile) | is.na(location_code)) %>% 
    nrow() %>% 
    divide_by(nrow(geo_data %>% filter(point == 0))) %>% 
    multiply_by(100) %>% round(2)
  message(paste0(pct_missing, "% of polygon data is missing either shapefile or location codes and are being dropped"))
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | point == 0]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1]
  return(geo_data)
  
}

# Subset for STI symptoms
subset_geomatched_sti_symptoms <- function(geo_data) {
  # Drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, year_n,
                      strata, psu, point, shapefile, location_code, latitude, longitude, sex_id, age_year, 
                      int_year, had_intercourse, discharge_or_sore, pweight, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # Drop Timor-Leste NIDS until further Non-SSA vetting
  drop_tls_nids <- c(21274)
  geo_data <- geo_data[!nid %in% drop_tls_nids,]
  
  #Drop DHS_SP surveys matched to polygons in same years as regular DHS data(see docs)
  geo_data <- geo_data[!(nid == 21198)]
 
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
  
  # rename variables
  setnames(geo_data, c("survey_series", "year_n"), c("source", "year"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  geo_data[, longitude := as.numeric(longitude)]
  
  # only want surveys 1998 or later for now
  geo_data <- geo_data[year >= 1998]
  
  # only want surveys that asked both men and women
  svys = geo_data[,list(sex_split = mean(sex_id, na.rm = T)), by = 'nid,country']
  geo_data <- geo_data[nid %in% svys[sex_split > 1 & sex_split < 2,nid]]
  
  # only want respondents who have had intercourse
  geo_data <- geo_data[had_intercourse == 1]
  
  # subset to ages 15-49 (and only surveys with this full range)
  geo_data <- geo_data[between(age_year, 15, 49),]
  drop <- geo_data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
  geo_data <- geo_data[!nid %in% drop,]
  
  # add a new variable combining the sti_symptoms: discharge_OR_sore
  geo_data[, sti_symptoms := discharge_or_sore]
  
  # drop observations with missing sti_symptoms or pweight
  geo_data <- geo_data[!is.na(pweight) & !is.na(sti_symptoms)]
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | point == 0]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1]

}

# RSP mismatch 
#



# In union 
subset_geomatched_in_union <- function(geo_data) {
  # Drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, survey_module, year,
                      strata, psu, point, shapefile, location_code, latitude, longitude,
                      sex_id, age_year, int_year, pweight, marital_status, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # rename variables
  setnames(geo_data, c("survey_series"), c("source"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  
  geo_data[, longitude := as.numeric(longitude)]
  
  # only want surveys 1998 and later for now
  geo_data <- geo_data[year >= 1998]
  
  # drop observations with missing marital_status or pweight
  geo_data <- geo_data[!is.na(pweight) & !is.na(marital_status)]
  
  # create a new indicator in_union that's true for respondents who are currently married/living with partner, false otherwise
  geo_data[, in_union := ifelse(marital_status < 3, 1, 0)]
  
  # subset to ages 15-49 (and only surveys with this full range)
  geo_data <- geo_data[between(age_year, 15, 49)]
  drop <- geo_data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
  geo_data <- geo_data[!nid %in% drop,]
  
  # Only keep surveys that asked both men and women about marital status
  svys = geo_data[, list(sex_split = mean(sex_id, na.rm = T)), by = 'nid,country']
  geo_data <- geo_data[nid %in% svys[sex_split > 1 & sex_split < 2, nid]]
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | point == 0]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1]
  
  # drop PMA WN mods to prevent duplication (see docs)
  geo_data <- geo_data[!(grepl("*PMA2020*", survey_name) & survey_module == "WN")]
  
  # drop COD PMA 2013-14, 2014, 2015, 2015-2016, 2016, 2018 surveys because of data oddities (see docs)
  #Drop DHS_SP survey matched to polygons(see docs) (nid = 21198, 21274) 
  geo_data <- geo_data[!nid %in% c(257822, 257823, 257826, 286019, 286054, 286020, 21198, 21274, 408008)]
  
  #drop IPUMS_Census census data
  geo_data <- geo_data[!nid %in% c(367640, 106473, 151296, 151304, 43726, 56492, 35329)]
  
  # drop MLI LSMS 2014-2015, KEN WHO STEPS_NCD 2015,
  # NAM HH INCOME & EXPENDITURE 2009-2010, UGA LSMS ISA 2009-2010,
  # BEN HH Living survey 2009 for data oddities (see docs)
  geo_data <- geo_data[!nid %in% c(260407, 165492, 134371, 81004, 151768)]
  
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
  
  #drop VNM AIS 13544 to prevent duplication(see docs)
  geo_data <- geo_data[!(nid == 13544 & survey_module == "WN")]
  
  #drop DHS VNM 2002 to prevent duplication(see docs)
  geo_data <- geo_data[!(nid == 21058 & survey_module == "WN")]
  
  #drop  MOZ UNICEF MICS 27031 2008 to having only WN asked (see docs)
  geo_data <- geo_data[!(nid == 27031 & survey_module == "WN")]
  
  # drop DHS WN mods where HHM mod was extracted to prevent duplication (see docs)
  geo_data <- geo_data[!((nid == 19539 | nid == 111432) & survey_module == "WN")]
  
  # drop ZAF ISSP surveys for data oddities (see docs)
  geo_data <- geo_data[!(source == "ISSP" & country == "ZAF")]
  
  #drop VNM STEP Skills 2012 survey where it only surveys non-granular urban areas(Hanoi, Ho Chi Minh City)(See doc)
  geo_data <- geo_data[!(nid == 299045)]
  
  # Dont keep this part of survey module
  geo_data[, survey_module := NULL]
}

subset_geomatched_partner_away <- function(geo_data) {
  # Drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, year,
                      strata, psu, point, shapefile, location_code, latitude, longitude,
                      sex_id, age_year, int_year, pweight, marital_status, partner_away, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # rename variables
  setnames(geo_data, c("survey_series"), c("source"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  geo_data[, longitude := as.numeric(longitude)]
  
  # only want surveys 1998 or later for now
  geo_data <- geo_data[year >= 1998]
  
  # only want women
  geo_data <- geo_data[sex_id == 2]
  
  # only want currently married/living with partner women
  geo_data <- geo_data[marital_status < 3, ]
  
  # drop observations with missing in_union or pweight
  geo_data <- geo_data[!is.na(pweight) & !is.na(partner_away)]
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | point == 0]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1]
  
  # subset to ages 15-49 (and only surveys with this full range)
  # doing this pre marital status restriction because some surveys
  # don't have the complete age range once you restrict to currently married
  geo_data <- geo_data[between(age_year, 15, 49)]
  drop <- geo_data[, .(low = range(age_year)[1], high = range(age_year)[2]), by = nid][low != 15 | high != 49, nid]
  geo_data <- geo_data[!nid %in% drop]
  
  # drop NGA PMA 2016 and BWA AIDS INDICATOR SURVEYS 2001, 2004 because of data oddities (see docs)
  geo_data <- geo_data[!nid %in% c(286022, 22114, 22112)]
  
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
  
  #Drop DHS_SP survey matched to polygons(see docs)
  geo_data <- geo_data[!(nid == 21198)]

}


subset_geomatched_had_intercourse <- function(geo_data) {
  
  # drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, year,
                      strata, psu, point, shapefile, location_code, latitude, longitude,
                      sex_id, age_year, int_year, pweight, had_intercourse, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # rename variables
  setnames(geo_data, c("survey_series"), c("source"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  geo_data[, longitude := as.numeric(longitude)]
  
  # only want women
  geo_data <- geo_data[sex_id == 2]
  
  # only want surveys 1998 and later for now
  geo_data <- geo_data[year >= 1998]
  
  # drop observations with missing had_intercourse or pweight
  geo_data <- geo_data[!is.na(pweight) & !is.na(had_intercourse)]
  
  # subset to ages 15-24 (and only surveys with this full range)
  geo_data <- geo_data[between(age_year, 15, 24)]
  drop <- geo_data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 24, nid]
  geo_data <- geo_data[!nid %in% drop]
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | (point == 0 & shapefile != ""),]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1,]
  
  #Drop Timor-Leste NIDS until further Non-SSA vetting
  drop_tls_nids <- c(21274)
  geo_data <- geo_data[!nid %in% drop_tls_nids]
  
  #Drop DHS_SP survey matched to polygons(see docs)
  geo_data <- geo_data[!(nid == 21198)]
  
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
  
  # drop GHA 2007-08 MICS WN (nid 160576) & CMR 2006 MICS WN (nid 2063) because of data quality issues (see docs)
  geo_data <- geo_data[!(nid %in% c(160576, 2063) & sex_id == 2)]
}


subset_geomatched_condom_last_time <- function(geo_data) {
  
  # drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, year, strata, psu, 
                      point, shapefile, location_code, latitude, longitude, sex_id, age_year,
                      int_year, pweight, sex_in_last_12_months, condom_last_time, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # rename variables
  setnames(geo_data, c("survey_series"), c("source"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  geo_data[, longitude := as.numeric(longitude)]
  
  # subset to ages 15-49 (and only surveys with this full range)
  # doing pre sex in last year restriction because one survey has
  # all 15-17 year-old-years who have not had sex in last year
  geo_data <- geo_data[between(age_year, 15, 49)]
  drop <- geo_data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
  geo_data <- geo_data[!nid %in% drop]
  
  # only want those who've had sex in the last 12 months
  geo_data <- geo_data[sex_in_last_12_months == 1]
  
  # drop observations with missing condom_last_time or pweight
  geo_data <- geo_data[!is.na(pweight) & !is.na(condom_last_time)]
  
  # only want surveys 1998 or later for now
  geo_data <- geo_data[year >= 1998]
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | point == 0]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1]
  
  #Drop Timor-Leste NIDS until further Non-SSA vetting
  drop_tls_nids <- c(21274)
  geo_data <- geo_data[!nid %in% drop_tls_nids]
  
  #Drop DHS_SP survey matched to polygons(see docs)
  geo_data <- geo_data[!(nid == 21198)]
  
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
  
  # Only keep surveys that asked both men and women about condom last time
  svys = geo_data[, list(sex_split = mean(sex_id, na.rm = T)), by = 'nid,country']
  geo_data <- geo_data[nid %in% svys[sex_split > 1 & sex_split < 2, nid]]
}

subset_geomatched_multiple_partner_year_mn <- function(geo_data) {
  # drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, year,
                      strata, psu, point, shapefile, location_code, latitude, longitude,
                      sex_id, age_year, int_year, pweight, num_partners_year, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # rename variables
  setnames(geo_data, c("survey_series"), c("source"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  geo_data[, longitude := as.numeric(longitude)]
  
  # only want surveys 1998 or later for now and only want men
  geo_data <- geo_data[year >= 1998 & sex_id == 1]
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | point == 0]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1]
  
  # subset to ages 15-49 (and only surveys with this full range)
  geo_data <- geo_data[between(age_year, 15, 49),]
  drop <- geo_data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
  geo_data <- geo_data[!nid %in% drop,]
  
  # Create indicator
  geo_data[, multiple_partners_year := ifelse(num_partners_year <= 1, 0, 1)]
  
  # drop observations with missing multiple_partners_year or pweight
  geo_data <- geo_data[!is.na(pweight) & !is.na(multiple_partners_year)]
  
  # drop all DHS3 surveys (plus DHS4 GIN & DHS4 GHA which also required a spousal adjust)
  # for differing question construction (see docs)
  dhs3_nids <- c(19076, 18531, 20537, 20132, 19198, 20909, 20796, 19305, 21090, 20382,
                 20780, 20212, 20301, 19493, 19546, 20852, 18938, 19370, 19292, 20976,
                 19604, 21139, 18519, 20120, 19670, 19614)
  geo_data <- geo_data[!(nid %in% dhs3_nids),]
  
  # drop GAB DHS 2000 because of data oddities
  geo_data <- geo_data[nid != 19579]
  
  #Drop DHS_SP survey matched to polygons(see docs)
  geo_data <- geo_data[!(nid == 21198)]
  
  #Drop Timor-Leste NIDS until further Non-SSA vetting
  drop_tls_nids <- c(21274)
  geo_data <- geo_data[!nid %in% drop_tls_nids,]
  
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
}

subset_geomatched_multiple_partner_year_wn <- function(geo_data) {
  # drop unneeded variables and remove annoying stata attributes
  geo_data <- geo_data[, list(nid, country, survey_series, survey_name, year,
                              strata, psu, point, shapefile, location_code, latitude, longitude,
                              sex_id, age_year, int_year, pweight, num_partners_year, hh_id, line_id, geospatial_id)]
  for (ii in 1:ncol(geo_data)) attributes(geo_data[[ ii]]) <- NULL
  
  # rename variables
  setnames(geo_data, c("survey_series"), c("source"))
  
  # fix variable class
  geo_data[, year := as.numeric(year)]
  geo_data[, point := as.numeric(point)]
  geo_data[, latitude := as.numeric(latitude)]
  geo_data[, longitude := as.numeric(longitude)]
  
  # only want surveys 1998 or later for now and only want men
  geo_data <- geo_data[year >= 1998 & sex_id == 2]
  
  # drop point data with missing latitude/longitude
  geo_data <- geo_data[(!is.na(latitude) & !is.na(longitude)) | point == 0]
  
  # drop polygon data with missing shapefile/location code
  geo_data <- geo_data[(!is.na(shapefile) & !is.na(location_code)) | point == 1]
  
  # subset to ages 15-49 (and only surveys with this full range)
  geo_data <- geo_data[between(age_year, 15, 49),]
  drop <- geo_data[, as.list(range(age_year)), nid][V1 != 15 | V2 != 49, nid]
  geo_data <- geo_data[!nid %in% drop,]
  
  # Create indicator
  geo_data[, multiple_partners_year := ifelse(num_partners_year <= 1, 0, 1)]
  
  # drop observations with missing multiple_partners_year or pweight
  geo_data <- geo_data[!is.na(pweight) & !is.na(multiple_partners_year)]
  
  # drop all DHS3 surveys (plus DHS4 GIN & DHS4 GHA which also required a spousal adjust)
  # for differing question construction (see docs)
  dhs3_nids <- c(19076, 18531, 20537, 20132, 19198, 20909, 20796, 19305, 21090, 20382,
                 20780, 20212, 20301, 19493, 19546, 20852, 18938, 19370, 19292, 20976,
                 19604, 21139, 18519, 20120, 19670, 19614)
  geo_data <- geo_data[!(nid %in% dhs3_nids),]
  
  # drop GAB DHS 2000 because of data oddities
  geo_data <- geo_data[nid != 19579]
  
  #Drop DHS_SP survey matched to polygons(see docs)
  geo_data <- geo_data[!(nid == 21198)]
  
  #Drop DR special survey that is not representative
  geo_data <- geo_data[!(nid == 165645)]
  
  #Drop Timor-Leste NIDS until further Non-SSA vetting
  drop_tls_nids <- c(21274)
  geo_data <- geo_data[!nid %in% drop_tls_nids,]
}

