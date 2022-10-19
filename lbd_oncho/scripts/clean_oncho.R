### Clean oncho data
### Author: LBD Focal 3
### Purpose: Oncho data from systematic review + ESPEN + OCP programs. This scipt checks for data quality issues and prepares it for use in MBG modeling
###
###          This involves:
###               - Standardizing fields
###               - Checking for, and removing duplicates
###               - Applying necessary data adjustments (age and diagnostic)
#############################################################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Script Parameters ######################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

produce_crosswalk_data_set <- FALSE ### produce crosswalk data set? If false, the final modeling data set will be produced instead
calculate_age_midpoints <- FALSE ### calculate age midpoints? (for INLA crosswalk)
add_remo_espen_year <- TRUE ### update survey years for pre-2001 nodule mapping data in ESPEN data set

save_outputs <- TRUE ### save cleaned data files?
drop_outliered_data <- TRUE ### drop outliered data?
output_crosswalk_quantiles <- FALSE ### save data files with had_oncho_poly set at lower and upper 95% of crosswalked values?

pathadd <- <<<< FILEPATH REDACTED >>>>
modeling_shapefile_version <- shapefile_version <- "2020_05_21"
raster_agg_factor <- 1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Script Prep #############################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Set repo locations
## Get current user name
user <- Sys.info()[["user"]] ## Get current user name

## Set repo location and indicator group
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))
package_list <- c(package_list, "fasterize", "sf")

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

# NOTE: package load errors likely due to R version problems
#       check .libPaths() and R.version. Change sing_imagename in config to singularity image with correct internals.
mbg_setup(package_list = package_list, repos = core_repo)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Load data ##############################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data_ocp <- fread(<<<< FILEPATH REDACTED >>>>)
data_espen <- fread(<<<< FILEPATH REDACTED >>>>)
data_sr <- fread(<<<< FILEPATH REDACTED >>>>)
data_epirf_new <- fread(<<<< FILEPATH REDACTED >>>>)

###### Combine SR data sets

###### Add LBD_group and LBD_uniq to data_sr rows lacking them
decideGroup <- function(data, groupNum) {
  
  ### Points ### ----
  
  data$LBD_group <- as.integer(data$LBD_group)
  data[shape_type == "point", LBD_group := .GRP
       , by = list(lat, long, year_start, year_end, month_start, month_end, site_memo)]
  data$LBD_group <- as.double(data$LBD_group)
  data[shape_type == "point", LBD_group := LBD_group + groupNum]
  # Update group number for polygons
  groupNum <- max(data$LBD_group, na.rm = TRUE) + 1
  
  ### Polygons ### ----
  
  data$LBD_group <- as.integer(data$LBD_group)
  data[shape_type == "polygon", LBD_group := .GRP
       , by = list(poly_type, poly_reference, poly_id, year_start, year_end, month_start, month_end, site_memo)]
  data$LBD_group <- as.double(data$LBD_group)
  data[shape_type == "polygon", LBD_group := LBD_group + groupNum]
  
  return(data)
}

decideUniq <- function(data) {
  
  # Loop through groups
  for (g in unique(data$LBD_group)) {
    
    # Check if only one site_memo
    if (length(unique(data[LBD_group == g]$site_memo)) == 1) {
      data[LBD_group == g, LBD_uniq := 1]
      
      # Otherwise, do a bunch of stuff
    } else {
      
      # Initiate unique numbering
      uNum <- 2
      
      # Loop through individual sites
      for (s in unique(data[LBD_group == g, ]$site_memo)) {
        
        # Pull T/F of whether current site should be an aggregate of the others
        aggregated <- aggCheck(s, unique(data[LBD_group == g, ]$site_memo))
        
        # Not aggregated - it's an underling, give it a uNum
        if (aggregated == FALSE) {
          data[LBD_group == g & site_memo == s, LBD_uniq := uNum]
          uNum <- uNum + 1
          
          # Yes aggregated - is overarching, give a unique 1 number
        } else {
          data[LBD_group == g & site_memo == s, LBD_uniq := 1]
        }
      }
    }
  }
  
  # Return Data
  message("...Unique Finished")
  return(data)
}

decideSpec <- function(data) {
  
  # Looping through by groups
  for (g in unique(data$LBD_group)) {
    
    # Pull individual group
    cGroup <- data[LBD_group == g, ]
    
    ### Check for DX Spec ### ----
    if(length(unique(cGroup$cv_DX_type)) != 1) {
      data[LBD_group == g, LBD_spec := "DX "]
    }
    
    # Looping through uniques to test sex and age
    for (u in unique(data$LBD_uniq)) {
      
      # Pull unique combo
      cUniq <- cGroup[LBD_uniq == u]
      
      # Looping through Diagnostics as well
      for (dx in unique(cUniq$cv_DX_type)) {
        
        cDX <- cUniq[cv_DX_type == dx]
        
        ### Test for Sex Spec ### ----
        if (length(unique(cDX$sex)) != 1 ) {
          data[LBD_group == g & LBD_uniq == u & sex != "Both" & cv_DX_type == dx
               , LBD_spec := paste0(LBD_spec, "Sex ")]
        }
        
        ### Test for Age Spec ### ----
        maxA <- max(cDX$age_end)
        minA <- min(cDX$age_start)
        data[LBD_group == g & LBD_uniq == u & cv_DX_type == dx
             & (age_start != minA | age_end != maxA)
             , LBD_spec := paste0(LBD_spec, "Age ")]
      }
    }
  }
  
  # Return new specified data
  message("...Specificity Finished")
  return(data)
}

srData <- decideGroup(data_sr[is.na(LBD_group)], max(data_sr$LBD_group, na.rm = TRUE))
srData <- decideUniq(srData)
# srData <- decideSpec(srData)
data_sr <- rbind(data_sr[!is.na(LBD_group)], srData)

# Fix wrong ISO3 in ESPEN data set (use country field to set ISO3)
data_espen[country == "burundi", iso3 := "BDI"]
data_espen[country == "ethiopia", iso3 := "ETH"]
data_espen[country == "malawi", iso3 := "MWI"]
data_espen[country == "mozambique", iso3 := "MOZ"]
data_espen[country %in% c("tanzania (mainland)"), iso3 := "TZA"]
data_espen[country == "uganda", iso3 := "UGA"]
data_espen[country == "angola", iso3 := "AGO"]
data_espen[country == "cameroon", iso3 := "CMR"]
data_espen[country == "equatorial guinea", iso3 := "GNQ"]
data_espen[country == "central african republic", iso3 := "CAF"]
data_espen[country == "chad", iso3 := "TCD"]
data_espen[country == "congo (brazzaville", iso3 := "COG"]
data_espen[country == "d.r. congo", iso3 := "COD"]
data_espen[country == "gabon", iso3 := "GAB"]
data_espen[country == "sudan", iso3 := "SDN"]
data_espen[country == "south sudan", iso3 := "SSD"]
data_espen[country == "benin", iso3 := "BEN"]
data_espen[country == "nigeria", iso3 := "NGA"]
data_espen[country == "liberia", iso3 := "LBR"]
data_espen[country == "sierra leone", iso3 := "SLE"]
data_espen[country == "togo", iso3 := "TGO"]
data_espen[country == "kenya", iso3 := "KEN"]
data_espen[country == "rwanda", iso3 := "RWA"]
data_espen[country == "cote d'ivoire", iso3 := "CIV"]

#### Update survey years for pre-2001 nodule mapping data in ESPEN data set
if (add_remo_espen_year) {
  # First, use onchomda covariate to extract first year of MDA at each data point
  gaul_list <- get_adm0_codes("oncho_endem_afr", shapefile_version = modeling_shapefile_version)
  
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = FALSE, shapefile_version = modeling_shapefile_version)
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  
  ## Build administrative and population rasters
  raster_list <- build_simple_raster_pop(subset_shape)
  simple_raster <- raster_list[["simple_raster"]]
  
  onchomda_config <- data.table(covariate = "onchomda", measure = "numrounds", release = "", year_lag = 0)
  interval_mo = 12
  onchomda <- load_lagged_covariates(covariate_config = onchomda_config,
                                template = simple_raster,
                                start_year = 2001,
                                end_year = 2001,
                                raster_agg = 1)$onchomda
  
  period_map <- make_period_map(modeling_periods = 1988:2018)
  
  nodule_year_missing <- data_espen[dx == "nod" & survey_year == 0 & source == "ESPEN"] # as of 20200504, these correspond to pre-2001 nodule mapping surveys without identified survey years
  nodule_year_not_missing <- data_espen[!(dx == "nod" & survey_year == 0 & source == "ESPEN")] # save other data rows for later appending
  nodule_year_missing$a_rowid <- 1:nrow(nodule_year_missing)
  setnames(nodule_year_missing, "survey_year", "year")
  
  cs_covs_no_cs <- extract_covariates(nodule_year_missing,
                                      onchomda,
                                      id_col = "a_rowid",
                                      return_only_results = TRUE,
                                      centre_scale = FALSE,
                                      period_var = "year",
                                      period_map = period_map)
  nodule_year_missing <- merge(nodule_year_missing, cs_covs_no_cs, by = "a_rowid", all = TRUE)
  setnames(nodule_year_missing, "year", "survey_year")
  nodule_year_missing$mda_start_year <- 2002 - nodule_year_missing$onchomda
  nodule_year_missing[onchomda == 0, mda_start_year := NA]
  nodule_year_missing$year_lower_bound <- as.integer(NA)
  nodule_year_missing$year_upper_bound <- nodule_year_missing$mda_start_year - 1
  nodule_year_missing[is.na(mda_start_year), year_upper_bound := 2000]
  
  # Assume 1992 as lower bound for CAF data rows
  nodule_year_missing[iso3 == "CAF", year_lower_bound := 1992]
  
  # Assume 1993 as lower bound for CMR data rows
  nodule_year_missing[iso3 == "CMR", year_lower_bound := 1993]
  
  # Assume 1997 as lower bound for COD data rows
  nodule_year_missing[iso3 == "COD", year_lower_bound := 1997]
  
  # Assume 1996 as lower bound for COG data rows
  nodule_year_missing[iso3 == "COG", year_lower_bound := 1996]
  
  # Assume 1997 as upper and lower bound for ETH data rows
  nodule_year_missing[iso3 == "ETH", c("year_lower_bound", "year_upper_bound") := 1997]
  
  # Assume 1996 as lower bound for KEN data rows
  nodule_year_missing[iso3 == "KEN", year_lower_bound := 1996]
  
  # Assume 1994 and 1995 as lower and upper bounds for NGA data rows
  nodule_year_missing[iso3 == "NGA", c("year_lower_bound", "year_upper_bound") := list(1994, 1995)]
  
  # Assume 1999 as lower and upper bounds for RWA data without MDA data
  nodule_year_missing[iso3 == "RWA", year_lower_bound := 1999]
  nodule_year_missing[iso3 == "RWA" & is.na(mda_start_year), year_upper_bound := 1999]
  
  # Assume 1996 as lower bound for SDN data rows
  nodule_year_missing[iso3 == "SDN", year_lower_bound := 1996]
  
  # Assume 1996 as lower bound for SSD data rows and 1998 as max. upper bound
  nodule_year_missing[iso3 == "SSD", year_lower_bound := 1996]
  nodule_year_missing[iso3 == "SSD" & year_upper_bound > 1998, year_upper_bound := 1998]
  
  # Assume 1994 as lower bound for TZA data rows
  nodule_year_missing[iso3 == "TZA", year_lower_bound := 1994]
  
  # Assume 1996 and 1998 as lower and upper bounds for TCD data rows, respectively
  nodule_year_missing[iso3 == "TCD", year_lower_bound := 1996]
  nodule_year_missing[iso3 == "TCD" & year_upper_bound > 1998, year_upper_bound := 1998]
  
  # Assume 1993 and 1997 as lower and upper bounds for UGA data rows, respectively
  nodule_year_missing[iso3 == "UGA", year_lower_bound := 1993]
  nodule_year_missing[iso3 == "UGA" & year_upper_bound > 1997, year_upper_bound := 1997]
  
  # Set survey year to midpoint of lower and upper bound (rounding down)
  nodule_year_missing$new_survey_year <- floor((nodule_year_missing$year_lower_bound + nodule_year_missing$year_lower_bound) / 2)
  
  # Set age range to 20+, if not otherwise provided
  nodule_year_missing[is.na(age_start), age_start := 20]
  nodule_year_missing[is.na(age_end), age_end := 94]
  
  nodule_year_missing[, survey_year := new_survey_year]
  
  data_espen <- rbind(nodule_year_missing, nodule_year_not_missing, fill = TRUE, use.names = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Standardize fields #####################################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Combine OCP and ESPEN data

setnames(data_ocp, "iso3", "ihme_loc_id")
setnames(data_espen, "iso3", "ihme_loc_id")
setnames(data_espen, "survey_year", "year")
setnames(data_espen, "survey_month", "month")

# Set ages for REMO data to 20+ (if not provided)
data_espen[source == "REMO" & dx == "nod" & is.na(age_start) & is.na(age_end), c("age_start", "age_end") := list(20, 94)]

# Also set ages for ESPEN mapping nodule data from APOC countries to 20+ (if not provided), assuming REMO; don't make the same assumption for CIV (unable to find information on methodology for 2014 nodule mapping surveys, but these surveys also reported mf prevalence, and sample sizes are larger than are typical for REMO)
data_espen[source == "ESPEN" & dx == "nod" & is.na(age_start) & is.na(age_end) & !(ihme_loc_id == "CIV"), c("age_start", "age_end") := list(20, 94)]

data_ocp_espen <- rbindlist(list(data_ocp, data_espen), fill=TRUE)


# Standardize column names
setnames(data_ocp_espen, "examined", "sample_size")
setnames(data_ocp_espen, "positive", "cases")
setnames(data_sr, "lat", "latitude")
setnames(data_sr, "long", "longitude")
setnames(data_sr, "poly_reference", "shapefile")
setnames(data_sr, "DX", "diagnostic")
setnames(data_ocp_espen, "dx", "diagnostic")

data_sr[diagnostic == "" & case_definition %in% c("OV16-positive", "Ov16 positive"), diagnostic := "sero"]
data_sr[diagnostic == "" & case_definition %in% c("presence of onchocercal nodules", "presence of o.v nodules"), diagnostic := "nod"]
data_sr$source <- "sys_rev"

## Note: the OCP data set does not contain age information.
all_data_rbind <- rbindlist(list(data_ocp_espen, data_sr), fill = TRUE)

# Retain only a subset of columns
retain_cols <- c("row_id", "nid",  "ihme_loc_id", "source", "year", "location_name", "sample_size", "cases", "latitude", "longitude", "shapefile", "poly_id", "diagnostic", "sex", "MDA_status", "LBD_group", "LBD_uniq", "LBD_spec", "site_memo", "age_start", "age_end", "field_citation_value", "file_path", "page_num", "table_num", "location_id", "source_type", "representative_name", "urbanicity_type", "recall_type", "latlong_source", "shape_type", "site_notes", "sampling_type", "sampling_description", "sample_strategy", "measure", "unit_type", "case_name", "case_definition", "case_diagnostics", "year_start", "year_end", "mean", "group", "unique", "specificity", "cv_MDA_years", "note_SR", "issue", "loc_dup")
all_data <- all_data_rbind[, ..retain_cols]

## Process new EPIRF data
colnames(data_epirf_new)
colnames(all_data)
setnames(data_epirf_new, c("ID", "Source", "ISO3", "LocationName", "Longitude", "Latitude", "SurveyYear", "Method_2", "Age_start", "Age_end", "Examined", "Positive"), c("row_id", "source", "ihme_loc_id", "location_name", "longitude", "latitude", "year", "diagnostic", "age_start", "age_end", "sample_size", "cases"))
data_epirf_new[diagnostic %in% c("microscopie", "Microscopy"), diagnostic := "ss"]
data_epirf_new[Method_0 == "BCE", diagnostic := "ss"]
data_epirf_new[Method_0 == "REMO", diagnostic := "nod"]
data_epirf_new$source <- "EPIRF"
data_epirf_new[ihme_loc_id == "Angola", ihme_loc_id := "AGO"]

all_data <- rbindlist(list(all_data, data_epirf_new), fill = TRUE)
all_data <- all_data[, ..retain_cols]

#### Data cleaning
all_data <- all_data[!(is.na(sample_size) & is.na(cases))]
all_data[nid == 0, nid := NA]
all_data[year == 0,  year := NA]
all_data$year_original <- all_data$year
all_data$year <- as.integer(all_data$year)
all_data[is.na(year) & !is.na(year_start) & !is.na(year_end), year := as.integer(round((year_start + year_end) / 2, 0))]
all_data[grep("_", ihme_loc_id), ihme_loc_id := substr(ihme_loc_id, 1, 3)]
all_data <- all_data[is.na(loc_dup) | loc_dup != 1]
all_data[, loc_dup := NULL]
all_data[, c("year_original", "year_start", "year_end") := NULL]
all_data <- all_data[!is.na(year)]
all_data <- all_data[!is.na(sample_size) & sample_size > 0]
all_data <- all_data[!(diagnostic %in% c("fly", "skin", "eye", "otherPrev", "PCR"))]
all_data[is.na(cases) & !is.na(mean) & !is.na(sample_size) & (mean > 1), cases := (sample_size * mean) / 100]
all_data[is.na(cases) & !is.na(mean) & !is.na(sample_size), cases := round(sample_size * mean, 0)]
setnames(all_data, "shape_type", "point")
all_data[!is.na(point) & point == "point", point := "1"]
all_data[!is.na(point) & point == "polygon", point := "0"]
all_data[is.na(point) & !is.na(latitude) & !is.na(longitude), point := "1"] # Some rows still need to be georeferenced; return to this
all_data[point == "", point := NA]
all_data[is.na(point) & recall_type == "point", point := "1"]
all_data[shapefile == "", shapefile := NA]
all_data[poly_id == "-1", poly_id := NA]
setnames(all_data, c("cases"), c("had_oncho_poly"))
all_data$weight <- 1
all_data[, issue := NULL]
all_data$LBD_spec <- tolower(all_data$LBD_spec)
setnames(all_data, c("sample_size"), c("N"))
all_data[, unique := NULL]
all_data[, group := NULL]
all_data[is.na(age_start), age_start := 0]
all_data[is.na(age_end), age_end := 94]
all_data[, mean := NULL]
all_data[, c("note_SR", "case_diagnostics", "case_definition", "sampling_description", "page_num", "table_num") := NULL]
all_data[is.na(site_memo), site_memo := ""]
all_data[is.na(LBD_uniq), LBD_uniq := ""]

### Check sample sizes and cases
nrow(all_data[had_oncho_poly > N | is.na(had_oncho_poly)]) # Should be zero
all_data <- all_data[had_oncho_poly <= N]

## Drop pre-1988 data
all_data <- all_data[year >= 1988]

## Drop non-African data
table(all_data$ihme_loc_id)
all_data <- all_data[!(ihme_loc_id %in% c("BRA", "COL", "ECU", "GTM", "MEX", "VEN"))]

### Georeferencing corrections
all_data[ihme_loc_id == "NGR", ihme_loc_id := "NER"] # Update to correct iso3 for Niger
all_data[row_id == "156_OC1", c("latitude", "longitude") := list(10.79583, -0.6169444)] # Latitude and longitude were switched

### Check for incorrect adm0
gaul_list <- get_adm0_codes("oncho_endem_afr", shapefile_version = modeling_shapefile_version)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = F, shapefile_version = modeling_shapefile_version)
subset_shape <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]

## Build administrative raster
raster_list <- build_simple_raster_pop(subset_shape)
simple_raster <- raster_list[["simple_raster"]]

l <- identify_wrong_adm0(data = all_data, s_poly = subset_shape, rast = simple_raster)
all_data <- l[[1]]
all_data_wrong_adm0 <- l[[2]]

## Retain corrected rows
all_data <- rbind(all_data, all_data_wrong_adm0[!is.na(latitude) & !is.na(longitude)], fill = TRUE)

### Check for extractions with both sex specificity and sex-age, sex-diagnostic, or sex-age-diagnostic 
all_data$sex <- tolower(all_data$sex)
all_data[sex == "", sex := "both"]

### Add unique LBD_group values to OCP data lacking them
max_LBD_group <- max(all_data[source == "OCP" & !is.na(LBD_group), LBD_group])
uniq_locs <- as.data.table(unique(all_data[source == "OCP" & is.na(LBD_group), list(latitude, longitude)]))
uniq_locs$LBD_group <- c((max_LBD_group + 1):(max_LBD_group + nrow(uniq_locs)))
for (i in 1:nrow(uniq_locs)) {
  all_data[source == "OCP" & latitude == uniq_locs[i, latitude] & longitude == uniq_locs[i, longitude] & is.na(LBD_group), LBD_group := uniq_locs[i, LBD_group]]
}

### Sex collapse (some data checks are needed; return to this)
no_sex_spec <- all_data[!(sex %in% c("male", "female"))]
sex_spec <- all_data[sex %in% c("male", "female")]

agg_sex <- merge(aggregate(N ~ LBD_group + source + age_start + age_end + diagnostic + site_memo + LBD_uniq + year, sex_spec, sum), aggregate(had_oncho_poly ~ LBD_group + source + age_start + age_end + diagnostic + site_memo + LBD_uniq + year, sex_spec, sum), all = TRUE)
agg_sex <- merge(agg_sex, aggregate(row_id ~ LBD_group + source + age_start + age_end + diagnostic + site_memo + LBD_uniq + year, sex_spec, paste0, collapse='_'), all = TRUE)
sex_collapsed <- unique(sex_spec[, -c("sex", "N", "had_oncho_poly", "row_id")])

merged <- as.data.table(merge(agg_sex, sex_collapsed, all = TRUE))
merged$sex <- "both"
merged[, LBD_spec := gsub("sex", "", LBD_spec)]
merged[, LBD_spec := gsub("  ", " ", LBD_spec)]
merged[LBD_spec == "dx ", LBD_spec := "dx"]
merged[LBD_spec == "", LBD_spec := NA]
merged <- setDT(merged)[, head(.SD, 1), by = row_id]

sex_collapsed_all <- rbindlist(list(no_sex_spec, merged), fill = TRUE)
sex_collapsed_all <- unique(sex_collapsed_all[,])
sex_collapsed_all$Master_UID <- sex_collapsed_all$row_id

#### Initial collapse of data rows by source, location, age, sex, year, and diagnostic, for specific nids
sex_collapsed_all$country <- sex_collapsed_all$ihme_loc_id
#### De-duplication (note: this does not include "duplication" based on different diagnostic tests or age groupings)
cyl <- unique(sex_collapsed_all[, .(latitude, longitude, shapefile, poly_id, year, ihme_loc_id)])
cyl$cyl.id <- c(1:nrow(cyl))
sex_collapsed_all <- merge(sex_collapsed_all, cyl, all.x = T)
cyl.tab <- sex_collapsed_all$cyl.id %>%
  table() %>%
  as.data.table()
setnames(cyl.tab, ".", "cyl.id")
oncho_check_dups <- cyl.tab[N > 1]
oncho_check_dups <- sex_collapsed_all[cyl.id %in% oncho_check_dups$cyl.id]

sex_collapsed_all$had_oncho_poly <- round(sex_collapsed_all$had_oncho_poly, 0)
sex_collapsed_all$N <- round(sex_collapsed_all$N, 0)

oncho_check_dups$country <- oncho_check_dups$ihme_loc_id
oncho_split_dups <- split_dups(data = oncho_check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = c("N", "had_oncho_poly"),
                               get_unique = c("diagnostic", "age_start", "age_end", "LBD_uniq"), preferred_source = "ESPEN")
oncho_obvi_drop <- oncho_split_dups[[1]]
oncho_nonobvi <- oncho_split_dups[[2]]
oncho_nonobvi_summary <- oncho_split_dups[[3]]

oncho <- sex_collapsed_all[!(Master_UID %in% oncho_obvi_drop)]

#### Collapse data rows by source, location, age, sex, year, and diagnostic
### Check for extractions with both dx specificity and dx-age or age specificity
oncho[LBD_spec == " age", LBD_spec := "age"]

oncho$age_start <- as.integer(round(oncho$age_start, 0))
oncho$age_end <- as.integer(round(oncho$age_end, 0))

oncho[age_end == 0, age_end := 94]
oncho[LBD_spec == "", LBD_spec := NA]

dx <- cbind(unique(oncho[LBD_spec == "dx", c("LBD_group", "LBD_uniq", "source", "diagnostic")]), "LBD_spec1" = "dx")
dx_age <- cbind(unique(oncho[LBD_spec == "dx age", c("LBD_group", "LBD_uniq", "source", "diagnostic")]), "LBD_spec2" = "dx age")
age <- cbind(unique(oncho[LBD_spec == "age", c("LBD_group", "LBD_uniq", "source", "diagnostic")]), "LBD_spec3" = "age")

LBD_spec_dx_age <- merge(dx, dx_age, all = TRUE)
LBD_spec_dx_age <- merge(LBD_spec_dx_age, age, all = TRUE)

LBD_spec_dx_age <- LBD_spec_dx_age[!(!is.na(LBD_spec1) & is.na(LBD_spec2) & is.na(LBD_spec3))]
LBD_spec_dx_age <- LBD_spec_dx_age[!(is.na(LBD_spec1) & !is.na(LBD_spec2) & is.na(LBD_spec3))]
LBD_spec_dx_age <- LBD_spec_dx_age[!(is.na(LBD_spec1) & is.na(LBD_spec2) & !is.na(LBD_spec3))]

drop_dx <- LBD_spec_dx_age[!is.na(LBD_spec1) & !is.na(LBD_spec2) & is.na(LBD_spec3)]

for (i in 1:nrow(drop_dx)) {
  drop_dx[i, "dx_age_min"] <- min(oncho[LBD_group == drop_dx[i, LBD_group] & LBD_uniq == drop_dx[i, LBD_uniq] & source == drop_dx[i, source] & diagnostic == drop_dx[i, diagnostic] & LBD_spec == "dx", age_start])
  drop_dx[i, "dx_age_max"] <- max(oncho[LBD_group == drop_dx[i, LBD_group] & LBD_uniq == drop_dx[i, LBD_uniq] & source == drop_dx[i, source] & diagnostic == drop_dx[i, diagnostic] & LBD_spec == "dx", age_end])
  drop_dx[i, "dx_age_age_min"] <- min(oncho[LBD_group == drop_dx[i, LBD_group] & LBD_uniq == drop_dx[i, LBD_uniq] & source == drop_dx[i, source] & diagnostic == drop_dx[i, diagnostic] & LBD_spec == "dx age", age_start])
  drop_dx[i, "dx_age_age_max"] <- max(oncho[LBD_group == drop_dx[i, LBD_group] & LBD_uniq == drop_dx[i, LBD_uniq] & source == drop_dx[i, source] & diagnostic == drop_dx[i, diagnostic] & LBD_spec == "dx age", age_end])
}

drop_dx_check <- drop_dx[(dx_age_min != dx_age_age_min) | (dx_age_max != dx_age_age_max)]
drop_dx <- drop_dx[(dx_age_min == dx_age_age_min) & (dx_age_max == dx_age_age_max)]

for (i in 1:nrow(drop_dx)) {
  oncho <- oncho[!((LBD_group == drop_dx[i, LBD_group]) & (LBD_uniq == drop_dx[i, LBD_uniq]) & (source == drop_dx[i, source]) & (diagnostic == drop_dx[i, diagnostic]) & (LBD_spec == "dx"))]
}

# For nid 332703, calculate prevalence for individuals 10+ by subtracting cases and sample sizes for children (<10) from all-age counts
temp_all_age <- as.data.table(oncho[nid == 332703 & age_end == 99])
temp_children <- as.data.table(oncho[nid == 332703 & age_end == 9])

for (i in 1:nrow(temp_all_age)) {
  if (temp_all_age[i, site_memo] != "Oke (village)|Edo State (state)|Nigeria (country)") {
    temp_all_age[i, "had_oncho_poly"] <- temp_all_age[i, had_oncho_poly] - temp_children[site_memo == temp_all_age[i, site_memo] & diagnostic == temp_all_age[i, diagnostic], had_oncho_poly]
    temp_all_age[i, "N"] <- temp_all_age[i, N] - temp_children[site_memo == temp_all_age[i, site_memo] & diagnostic == temp_all_age[i, diagnostic], N]
  }
  temp_all_age[i, age_start := 10]
  temp_all_age[i, LBD_spec := "dx age"]
}

oncho <- rbindlist(list(oncho[is.na(nid) | (nid != 332703)], temp_all_age, temp_children))

### Drop dx rows for nids 332709, 332730
oncho <- oncho[is.na(nid) | !(nid == 332709 & LBD_spec == "dx")]
oncho <- oncho[is.na(nid) | !(nid == 332730 & LBD_spec == "dx")]

# For nid 327919, calculate prevalence for individuals 20+ by subtracting cases and sample sizes for individuals <20 from all-age counts
temp_all_age <- as.data.table(oncho[nid == 327919 & age_start == 5 & age_end == 99])
temp_children <- as.data.table(oncho[nid == 327919 & age_start == 5 & age_end == 19])

for (i in 1:nrow(temp_all_age)) {
  temp_all_age[i, "had_oncho_poly"] <- temp_all_age[i, had_oncho_poly] - temp_children[site_memo == temp_all_age[i, site_memo] & diagnostic == temp_all_age[i, diagnostic], had_oncho_poly]
  temp_all_age[i, "N"] <- temp_all_age[i, N] - temp_children[site_memo == temp_all_age[i, site_memo] & diagnostic == temp_all_age[i, diagnostic], N]
  temp_all_age[i, age_start := 20]
  temp_all_age[i, LBD_spec := "dx age"]
}

oncho <- rbindlist(list(oncho[!(row_id %in% c(temp_all_age$row_id, temp_children$row_id))], temp_all_age, temp_children))

## Ad hoc deduplication
oncho <- oncho[!(Master_UID %in% c("14340_SR"))]

## Remove duplicate rows that represent different skin snip diagnostics from the same individuals (DX1 vs. DX10; retain the latter (more sensitive))
oncho <- oncho[!(Master_UID %in% c("4842_SR", "4844_SR", "4846_SR", "4848_SR", "4850_SR"))]

#### Additional deduplication and cleaning
# Run additional automated deduplication using site_memo, for ESPEN data
oncho[, cyl.id := NULL]
cyl <- unique(oncho[source == "ESPEN", .(site_memo, year, diagnostic)])
cyl$cyl.id <- c(1:nrow(cyl))
oncho_dedup <- merge(oncho[source == "ESPEN"], cyl, all.x = T)
cyl.tab <- oncho_dedup$cyl.id %>%
  table() %>%
  as.data.table()
setnames(cyl.tab, ".", "cyl.id")
oncho_check_dups <- cyl.tab[N > 1]
oncho_check_dups <- oncho_dedup[cyl.id %in% oncho_check_dups$cyl.id]
oncho_check_dups$Master_UID <- oncho_check_dups$row_id

oncho_split_dups <- split_dups(data = oncho_check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = c("N", "had_oncho_poly"),
                               get_unique = c("diagnostic"), preferred_source = "ESPEN")
oncho_obvi_drop <- oncho_split_dups[[1]]
oncho_nonobvi <- oncho_split_dups[[2]]
oncho_nonobvi_summary <- oncho_split_dups[[3]]

oncho <- oncho[!(row_id %in% oncho_obvi_drop)]

# Run additional automated deduplication using location_name, for ESPEN and OCP data
oncho[, cyl.id := NULL]
cyl <- unique(oncho[source %in% c("ESPEN", "OCP"), .(location_name, year, diagnostic)])
cyl$cyl.id <- c(1:nrow(cyl))
oncho_dedup <- merge(oncho[source %in% c("ESPEN", "OCP")], cyl, all.x = T)
cyl.tab <- oncho_dedup$cyl.id %>%
  table() %>%
  as.data.table()
setnames(cyl.tab, ".", "cyl.id")
oncho_check_dups <- cyl.tab[N > 1]
oncho_check_dups <- oncho_dedup[cyl.id %in% oncho_check_dups$cyl.id]
oncho_check_dups$Master_UID <- oncho_check_dups$row_id

oncho_split_dups <- split_dups(data = oncho_check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = c("N", "had_oncho_poly"),
                               get_unique = c("diagnostic"), preferred_source = "ESPEN")
oncho_obvi_drop <- oncho_split_dups[[1]]
oncho_nonobvi <- oncho_split_dups[[2]]
oncho_nonobvi_summary <- oncho_split_dups[[3]]

oncho <- oncho[!(row_id %in% oncho_obvi_drop)]

oncho$location_name <- tolower(oncho$location_name)

# Remove OCP/ESPEN duplicates
oncho <- oncho[!(row_id %in% c("275_OC2_275_OC2", "276_OC2_276_OC2", "281_OC2_281_OC2", "286_OC2_286_OC2", "278_OC2_278_OC2", "384_OC2_384_OC2", "290_OC2_290_OC2", "288_OC2_288_OC2", "287_OC2_287_OC2", "291_OC2_291_OC2", "171_OC1", "289_OC2_289_OC2", "156_OC2_156_OC2", "285_OC2_285_OC2"))]

# Remove additional duplicates
oncho <- oncho[!(row_id %in% c("14580_SR_14581_SR", "52_OC2_52_OC2", "254_OC2_254_OC2", "14583_SR_14584_SR", "32_OC1"))]

# Run additional automated deduplication to remove duplicates between REMO and ESPEN sources for CAF (and probably other countries)
oncho[, cyl.id := NULL]
cyl <- unique(oncho[source %in% c("ESPEN", "REMO"), .(location_name, year, diagnostic)])
cyl$cyl.id <- c(1:nrow(cyl))
oncho_dedup <- merge(oncho[source %in% c("ESPEN", "REMO")], cyl, all.x = T)
cyl.tab <- oncho_dedup$cyl.id %>%
  table() %>%
  as.data.table()
setnames(cyl.tab, ".", "cyl.id")
oncho_check_dups <- cyl.tab[N > 1]
oncho_check_dups <- oncho_dedup[cyl.id %in% oncho_check_dups$cyl.id]
oncho_check_dups$Master_UID <- oncho_check_dups$row_id

oncho_split_dups <- split_dups(data = oncho_check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = c("N", "had_oncho_poly"),
                               get_unique = c("diagnostic", "age_start", "age_end"), preferred_source = "ESPEN")
oncho_obvi_drop <- oncho_split_dups[[1]]
oncho_nonobvi <- oncho_split_dups[[2]]
oncho_nonobvi_summary <- oncho_split_dups[[3]]

oncho <- oncho[!(row_id %in% oncho_obvi_drop)]

# Remove additional duplicates
oncho <- oncho[!(row_id %in% c("157_OC1", "14586_SR_14587_SR", "14589_SR_14590_SR", "10386_SR_10387_SR", "7257_SR_7263_SR", "11339_SR", "5332_SR", "5243_SR", "11270_SR", "11336_SR", "3331_REM", "3358_REM", "2536_REM", "3337_REM", "3350_REM", "3344_REM", "2621_REM", "3383_REM", "3390_REM", "3367_REM", "10735_SR", "2605_REM", "2789_REM"))]

oncho <- oncho[!(row_id %in% c("15894_SR", "15895_SR", "2497_REM", "12701_REM", "12616_REM", "12267_REM", "12387_REM", "12376_REM", "12383_REM", "11882_REM", "11998_REM", "12011_REM", "11909_REM", "7726_SR", "7721_SR", "10190_SR", "10187_SR", "10195_SR", "7505_SR_7517_SR", "7498_SR_7510_SR", "7620_SR_7621_SR", "7499_SR_7511_SR", "7435_SR_7442_SR", "7434_SR_7441_SR", "7497_SR_7509_SR", "7433_SR_7440_SR", "7436_SR_7443_SR", "7454_SR_7461_SR", "7500_SR_7512_SR", "7501_SR_7513_SR", "7437_SR_7444_SR"))]

# Run additional automated deduplication
oncho[, cyl.id := NULL]
cyl <- unique(oncho[, .(longitude, latitude, year, N, had_oncho_poly)])
cyl$cyl.id <- c(1:nrow(cyl))
oncho_dedup <- merge(oncho, cyl, all.x = T)
cyl.tab <- oncho_dedup$cyl.id %>%
  table() %>%
  as.data.table()
setnames(cyl.tab, ".", "cyl.id")
oncho_check_dups <- cyl.tab[N > 1]
oncho_check_dups <- oncho_dedup[cyl.id %in% oncho_check_dups$cyl.id]
oncho_check_dups$Master_UID <- oncho_check_dups$row_id

oncho_split_dups <- split_dups(data = oncho_check_dups, group_id = "cyl.id", check_cols = c("diagnostic", "age_start", "age_end"), calc_var_cols = c("N", "had_oncho_poly"),
                               get_unique = c(), preferred_source = "ESPEN")
oncho_obvi_drop <- oncho_split_dups[[1]]
oncho_nonobvi <- oncho_split_dups[[2]]
oncho_nonobvi_summary <- oncho_split_dups[[3]]

oncho <- oncho[!(row_id %in% oncho_obvi_drop & ihme_loc_id != "TGO")]

# Remove additional duplicates
oncho <- oncho[!(row_id %in% c("156_OC1", "377_OC2_377_OC2", "14769_SR_14783_SR", "14763_SR_14777_SR", "14757_SR_14771_SR", "6953_SR", "14767_SR_14781_SR", "14759_SR_14773_SR", "6950_SR", "5697_REM", "5696_REM", "5718_REM", "12254_SR_12266_SR", "12218_SR_12230_SR", "12249_SR_12261_SR", "12213_SR_12225_SR", "13478_SR_13484_SR", "13476_SR_13482_SR", "12350_SR", "12210_SR_12222_SR"))]
oncho <- oncho[!(row_id %in% c("57_OC2_57_OC2", "53_OC2_53_OC2", "8119_SR", "8114_SR", "8112_SR", "8534_SR_8543_SR", "8120_SR", "8115_SR", "8122_SR", "9848_SR", "8123_SR", "15306_SR_15311_SR", "334_OC2_334_OC2", "125_OC2_125_OC2", "13522_SR_13523_SR", "13506_SR_13507_SR", "7345_SR", "7357_SR", "7353_SR", "7343_SR", "7344_SR", "7358_SR", "7359_SR"))]

### Deal with OCP/ESPEN rows having identical locations, N and cases but different years
# Run additional automated deduplication to remove duplicates between REMO and ESPEN sources for CAF (and probably other countries)
oncho[, cyl.id := NULL]
cyl <- unique(oncho[source %in% c("ESPEN", "OCP") & ihme_loc_id %in% c("BEN", "GHA", "GIN"), .(location_name, N, had_oncho_poly, age_start, age_end, diagnostic)])
cyl$cyl.id <- c(1:nrow(cyl))
oncho_dedup <- merge(oncho[source %in% c("ESPEN", "OCP") & ihme_loc_id %in% c("BEN", "GHA", "GIN")], cyl, all.x = T)
cyl.tab <- oncho_dedup$cyl.id %>%
  table() %>%
  as.data.table()
setnames(cyl.tab, ".", "cyl.id")
oncho_check_dups <- cyl.tab[N > 1]
oncho_check_dups <- oncho_dedup[cyl.id %in% oncho_check_dups$cyl.id]
oncho_check_dups$Master_UID <- oncho_check_dups$row_id

oncho_split_dups <- split_dups(data = oncho_check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = c("N", "had_oncho_poly"),
                               get_unique = c("diagnostic", "age_start", "age_end"), preferred_source = "ESPEN")
oncho_obvi_drop <- oncho_split_dups[[1]]
oncho_nonobvi <- oncho_split_dups[[2]]
oncho_nonobvi_summary <- oncho_split_dups[[3]]

oncho <- oncho[!(row_id %in% oncho_obvi_drop)]

oncho <- oncho[!(row_id %in% c("274_OC2_274_OC2", "369_OC2_369_OC2", "349_OC2_349_OC2", "534_OC2_534_OC2", "693_OC2_693_OC2", "538_OC2_538_OC2", "692_OC2_692_OC2", "516_OC2_516_OC2", "699_OC2_699_OC2", "518_OC2_518_OC2", "515_OC2_515_OC2", "539_OC2_539_OC2", "697_OC2_697_OC2", "537_OC2_537_OC2", "696_OC2_696_OC2", "540_OC2_540_OC2", "695_OC2_695_OC2"))]

oncho <- oncho[!(row_id %in% c("183_OC1", "117_OC1"))]

### Additional deduplication
## Drop serological data
oncho <- oncho[diagnostic != "sero"]

# Assess countries individually
table(oncho$ihme_loc_id)

full_output <- copy(oncho)
full_output$latitude_rounded <- round(full_output$latitude, 2)
full_output$longitude_rounded <- round(full_output$longitude, 2)
full_output$latitude_rounded_1 <- round(full_output$latitude, 1)
full_output$longitude_rounded_1 <- round(full_output$longitude, 1)

### Compare REMO and ESPEN to find dups that differ in sample size and cases but have identical coords
check_dups <- create_check_dups_no_country(full_output[source %in% c("ESPEN", "REMO")])
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "N", "had_oncho_poly", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

### Drop REMO data when there are corresponding ESPEN data (same year and coordinates; ESPEN only reports males, while the REMO [APOC] data have both sexes; prefer the male data for consistency)
check_dups <- create_check_dups_no_country(full_output[source %in% c("ESPEN", "REMO")])
full_output <- full_output[!(Master_UID %in% check_dups[source == "REMO", Master_UID])]

## AGO
check_dups <- create_check_dups("AGO")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "N", "had_oncho_poly", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
obv_dups <- dups[[1]]
full_output <- full_output[!(Master_UID %in% obv_dups)]

check_dups <- create_check_dups("AGO")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "N", "had_oncho_poly", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
full_output <- full_output[!(Master_UID %in% check_dups[source == "EPIRF", Master_UID])] # Drop EPIRF

## BDI
check_dups <- create_check_dups("BDI")

## BEN
check_dups <- create_check_dups("BEN")

## BFA
check_dups <- create_check_dups("BFA")

## CAF
check_dups <- create_check_dups("CAF")

## CIV
check_dups <- create_check_dups("CIV")

## CMR
check_dups <- create_check_dups("CMR")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "N", "had_oncho_poly", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")
full_output <- full_output[!(Master_UID %in% c("4584_SR", "4585_SR", "4586_SR", "4583_SR", "4591_SR"))] # prefer ESPEN rows

## COD
check_dups <- create_check_dups("COD")
check_dups <- create_check_dups_no_year_raw_coords_location_name("COD")

## COG
check_dups <- create_check_dups("COG")

## ETH
check_dups <- create_check_dups("ETH")
full_output <- full_output[!(nid %in% c(327885))] # prefer ESPEN rows

## GAB
check_dups <- create_check_dups("GAB")

## GHA
check_dups <- create_check_dups("GHA")

## GIN
check_dups <- create_check_dups("GIN")
check_dups <- create_check_dups_no_year("GIN")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "N", "had_oncho_poly", "age_start", "age_end", "country"), preferred_source = "ESPEN")

## GNQ
check_dups <- create_check_dups("GNQ")

## KEN
check_dups <- create_check_dups("KEN")

## LBR
check_dups <- create_check_dups("LBR")

## MLI
check_dups <- create_check_dups("MLI")

## MOZ
check_dups <- create_check_dups("MOZ")

## MWI
check_dups <- create_check_dups("MWI")

## NGA
check_dups <- create_check_dups("NGA")
full_output <- full_output[!(Master_UID %in% c("5743_SR", "5748_SR", "5747_SR", "5729_SR"))] # prefer ESPEN rows
full_output <- full_output[!(Master_UID %in% c("9360_SR"))] # sex deduplication
full_output <- full_output[!(Master_UID %in% c("28_SR"))] # diagnostic deduplication

## RWA
check_dups <- create_check_dups("RWA")

## SDN
check_dups <- create_check_dups("SDN")

## SEN
check_dups <- create_check_dups("SEN")

## SLE
check_dups <- create_check_dups("SLE")

## SSD
check_dups <- create_check_dups("SSD")

## TCD
check_dups <- create_check_dups("TCD")

## TGO
check_dups <- create_check_dups("TGO")

## TZA
check_dups <- create_check_dups("TZA")

## UGA
check_dups <- create_check_dups("UGA")
dups <- split_dups_obvi_only(data = check_dups, group_id = "cyl.id", check_cols = c("N", "had_oncho_poly"), calc_var_cols = NULL, get_unique = c("diagnostic", "N", "had_oncho_poly", "age_start", "age_end", "country", "year"), preferred_source = "ESPEN")

oncho <- copy(full_output)

## Create temporary nids for ESPEN and REMO data
oncho[is.na(nid) & source == "ESPEN", nid := -1]
oncho[is.na(nid) & source == "REMO", nid := -2]

oncho$counter <- 1

oncho <- oncho[!(Master_UID %in% c("8864_SR_8867_SR", "7262_SR_7268_SR", "7259_SR_7265_SR", "7261_SR_7267_SR", "7258_SR_7264_SR", "7260_SR_7266_SR", "9844_SR", "9135_SR_9162_SR", "9136_SR_9163_SR", "9138_SR_9165_SR", "9139_SR_9166_SR", "9140_SR_9167_SR", "8897_SR_8924_SR", "9142_SR_9169_SR", "9143_SR_9170_SR", "9144_SR_9171_SR", "9145_SR_9172_SR", "9146_SR_9173_SR", "9147_SR_9174_SR", "9148_SR_9175_SR", "9153_SR_9180_SR", "9157_SR_9184_SR", "9159_SR_9186_SR", "9150_SR_9177_SR", "9151_SR_9178_SR", "9156_SR_9183_SR", "6496_SR_6503_SR", "9367_SR_9368_SR", "5050_SR"))]
oncho[Master_UID %in% c("5274_SR", "5340_SR"), LBD_uniq := 2] # most likely a data entry error; sample size is different from other LBD_uniq == 1
oncho[nid == 324729 & case_name == "presence of nodules", diagnostic := "nod"] # correct inaccurate diagnostics

### Prep for crosswalk
LBD_diag_agg <- as.data.table(aggregate(diagnostic ~ LBD_group + LBD_uniq + source + age_start + age_end, oncho, unique))
LBD_diag_agg <- LBD_diag_agg[!(diagnostic %in% c("ss", "nod"))]

oncho$use_in_xwalk <- FALSE

for (i in 1:nrow(LBD_diag_agg)) {
  oncho[(LBD_group == LBD_diag_agg[i, LBD_group]) & (LBD_uniq == LBD_diag_agg[i, LBD_uniq]) & (source == LBD_diag_agg[i, source]) & (age_start == LBD_diag_agg[i, age_start]) & (age_end == LBD_diag_agg[i, age_end]), weight := weight / 2]
  oncho[(LBD_group == LBD_diag_agg[i, LBD_group]) & (LBD_uniq == LBD_diag_agg[i, LBD_uniq]) & (source == LBD_diag_agg[i, source]) & (age_start == LBD_diag_agg[i, age_start]) & (age_end == LBD_diag_agg[i, age_end]), use_in_xwalk := TRUE]
}

### Now assign studies with multiple age groups to crosswalk data set
oncho$age_range <- paste0(oncho$age_start, "_", oncho$age_end)
age_group_agg <- as.data.table(aggregate(age_range ~ LBD_group + LBD_uniq + source + diagnostic, oncho, unique))
age_group_agg <- age_group_agg[age_range != "0_94"]
for (i in 1:nrow(age_group_agg)) {
  age_group_agg[i, "age_range_length"] <- length(unique(age_group_agg[i, age_range][[1]]))
}

age_group_agg <- age_group_agg[age_range_length > 1]
for (i in 1:nrow(age_group_agg)) {
  oncho[(LBD_group == age_group_agg[i, LBD_group]) & (LBD_uniq == age_group_agg[i, LBD_uniq]) & (source == age_group_agg[i, source]) & (diagnostic == age_group_agg[i, diagnostic]), use_in_xwalk := TRUE]
}

oncho[is.na(LBD_group) & source == "ESPEN", c("LBD_group", "LBD_uniq") := list((max(oncho[source == "ESPEN", LBD_group], na.rm = TRUE) + 1):(max(oncho[source == "ESPEN", LBD_group], na.rm = TRUE) + nrow(oncho[is.na(LBD_group) & source == "ESPEN"])), 1)]
oncho[is.na(LBD_group) & source == "EPIRF", c("LBD_group", "LBD_uniq") := list(1:nrow(oncho[is.na(LBD_group) & source == "EPIRF"]), 1)]
oncho[, "cohort_id" := paste0(source, "_", LBD_group, "_", LBD_uniq)]

retain_cols <- c("nid", "Master_UID", "ihme_loc_id", "source", "year", "had_oncho_poly", "N", "latitude", "longitude", "shapefile", "poly_id", "diagnostic", "age_start", "age_end", "weight", "point", "use_in_xwalk", "cohort_id")
oncho <- oncho[, ..retain_cols]

setnames(oncho, "ihme_loc_id", "country")
setnames(oncho, "poly_id", "location_code")
oncho$sampling <- 1
oncho$cluster_id <- 1
oncho$data_collect_method <- NA
oncho$original_year <- NA
oncho[, shapefile := gsub(<<<< FILEPATH REDACTED >>>>)]
oncho[, shapefile := gsub(<<<< FILEPATH REDACTED >>>>)]

# Save serological data (will append to crosswalk data set later)
sero_data <- oncho[diagnostic == "sero"]
oncho <- oncho[diagnostic != "sero"]

#### Generate crosswalk training data set, or apply crosswalk and generate final modeling data set
if (produce_crosswalk_data_set) {
  ## Check for any remaining deduplication needs
  # First check ss
  Data <- oncho[oncho$diagnostic == "ss" & use_in_xwalk == TRUE]
  Data$uid_crosswalk <- 1:nrow(Data)
  
  ### Subset to studies reporting > 1 age group for diagnostic
  age_groups_by_nid <- as.data.table(aggregate(uid_crosswalk ~ cohort_id, data = Data, length))
  Data <- Data[cohort_id %in% age_groups_by_nid[uid_crosswalk > 1, cohort_id]]
  Data$prev <- Data$had_oncho_poly / Data$N
  
  ## Drop all cohorts with 0 total cases
  total_cases_by_cohort <- as.data.table(aggregate(had_oncho_poly ~ cohort_id, Data, sum))
  Data <- Data[cohort_id %in% total_cases_by_cohort[had_oncho_poly != 0, cohort_id]]
  
  Data <- Data[!(Master_UID %in% c("7618_SR", "7619_SR", "5478_SR", "8405_SR", "7303_SR_7304_SR", "9862_SR_9868_SR", "9360_SR", "15909_SR_15910_SR"))]
  Data[Master_UID == "5483_SR", c("Master_UID", "had_oncho_poly", "N") := list("5483_SR_5478_SR", 94 + 44, 117 + 63)]
  Data[Master_UID == "5756_SR", c("had_oncho_poly", "N", "age_start", "age_end") := list(379 - 350, 1024 - 815, 1, 9)]
  Data[Master_UID == "5758_SR", c("had_oncho_poly", "N", "age_start", "age_end") := list(54 - 54, 1525 - 1249, 1, 9)]
  
  oncho <- oncho[!(Master_UID %in% c("7618_SR", "7619_SR", "5478_SR", "8405_SR", "7303_SR_7304_SR", "9862_SR_9868_SR", "9360_SR", "15909_SR_15910_SR"))]
  oncho[Master_UID == "5483_SR", c("Master_UID", "had_oncho_poly", "N") := list("5483_SR_5478_SR", 94 + 44, 117 + 63)]
  oncho[Master_UID == "5756_SR", c("had_oncho_poly", "N", "age_start", "age_end") := list(379 - 350, 1024 - 815, 1, 9)]
  oncho[Master_UID == "5758_SR", c("had_oncho_poly", "N", "age_start", "age_end") := list(54 - 54, 1525 - 1249, 1, 9)]
  
  # Then check nod
  Data <- oncho[oncho$diagnostic == "nod" & use_in_xwalk == TRUE]
  Data$uid_crosswalk <- 1:nrow(Data)
  
  ### Subset to studies reporting > 1 age group for diagnostic
  age_groups_by_nid <- as.data.table(aggregate(uid_crosswalk ~ cohort_id, data = Data, length))
  Data <- Data[cohort_id %in% age_groups_by_nid[uid_crosswalk > 1, cohort_id]]
  Data$prev <- Data$had_oncho_poly / Data$N
  
  ## Drop all cohorts with 0 total cases
  total_cases_by_cohort <- as.data.table(aggregate(had_oncho_poly ~ cohort_id, Data, sum))
  Data <- Data[cohort_id %in% total_cases_by_cohort[had_oncho_poly != 0, cohort_id]]
  
  Data <- Data[!(nid == "332893" & age_start == 15 & age_end == 99)]
  oncho <- oncho[!(nid == "332893" & age_start == 15 & age_end == 99)]
  
  # Generate final crosswalk data set
  crosswalk_data <- oncho[oncho$diagnostic %in% c("ss", "nod") & use_in_xwalk == TRUE & source == "sys_rev"]
  
  # Create unique IDs and rename some columns for consistency with crosswalk scripts
  crosswalk_data$uid_crosswalk <- 1:nrow(crosswalk_data)
  crosswalk_data$a_rowid <- 1:nrow(crosswalk_data)
  setnames(crosswalk_data, "cohort_id", "cohort")
  
  ## Subset to studies reporting > 1 age group for either diagnostic (ss or nod)
  age_groups_by_nid_ss <- as.data.table(aggregate(uid_crosswalk ~ cohort, data = crosswalk_data[diagnostic == "ss"], length))
  age_groups_by_nid_nod <- as.data.table(aggregate(uid_crosswalk ~ cohort, data = crosswalk_data[diagnostic == "nod"], length))
  crosswalk_data <- crosswalk_data[(cohort %in% age_groups_by_nid_ss[uid_crosswalk > 1, cohort]) | (cohort %in% age_groups_by_nid_nod[uid_crosswalk > 1, cohort])]
  
  #### Load shared db functions
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  
  country_years <- unique(crosswalk_data[, c("country", "year")])
  
  age_ids <- c(49:142) # Age 1 = ID 49, Age 94 = ID 142; get_population does not return single_year_age for ages 95-99
  gbd_round_id <- 6 ### 2019
  
  location_hierarchy <- get_location_metadata(location_set_id = 35, gbd_round_id = gbd_round_id)
  age_metadata <- get_age_metadata(age_group_set_id = 12, gbd_round_id = gbd_round_id)
  
  country_years <- merge(country_years, location_hierarchy[, c("ihme_loc_id", "location_id")], by.x = "country", by.y = "ihme_loc_id", all.x = TRUE)
  
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
  
  ## Add location_id to oncho
  crosswalk_data <- merge(crosswalk_data, unique(country_years[, c("country", "location_id")]), by = "country", all.x = TRUE)
  crosswalk_data$location_id <- as.integer(crosswalk_data$location_id)
  
  ## Merge age distribution with crosswalk data
  newData <- merge(crosswalk_data, pop_age_merged_wide, by = c("location_id", "year"), all.x = TRUE)
  names(newData)[names(newData) == "age_prop.<1 year"] <- "age_prop.0" # Rename age_prop.<1 year
  
  newData[age_end > 94, age_end := 94]
  newData[is.na(age_start), age_start := 0] # If no age_start provided, assume 0
  newData[is.na(age_end), age_end := 94] # If no age_end provided, assume 94
  
  crosswalk_data <- copy(newData)
  setnames(crosswalk_data, "had_oncho_poly", "cases")
  setnames(crosswalk_data, "N", "sample_size")
  
  ## Final clean-ups
  setnames(crosswalk_data, "age_prop.<1 year", "age_prop.0", skip_absent = TRUE)
  crosswalk_data$age_start <- as.integer(crosswalk_data$age_start)
  crosswalk_data$age_end <- as.integer(crosswalk_data$age_end)
  crosswalk_data[age_start > 94, age_start := 94]
  crosswalk_data[age_end > 94, age_end := 94]
  crosswalk_data <- crosswalk_data[!is.na(sample_size)]
  
  ## Drop all cohorts with 0 total cases
  total_cases_by_cohort <- as.data.table(aggregate(cases ~ cohort, crosswalk_data, sum))
  crosswalk_data <- crosswalk_data[cohort %in% total_cases_by_cohort[cases != 0, cohort]]
  
  ## Save crosswalk training data set to file
  write.csv(crosswalk_data, file=<<<< FILEPATH REDACTED >>>>, row.names=F)
} else { # Apply crosswalks and generate final data set for MBG modeling
  ## Apply data fixes identified via crosswalk prep, above
  oncho <- oncho[!(Master_UID %in% c("7618_SR", "7619_SR", "5478_SR", "8405_SR", "7303_SR_7304_SR", "9862_SR_9868_SR", "9360_SR", "15909_SR_15910_SR"))]
  oncho[Master_UID == "5483_SR", c("Master_UID", "had_oncho_poly", "N") := list("5483_SR_5478_SR", 94 + 44, 117 + 63)]
  oncho[Master_UID == "5756_SR", c("had_oncho_poly", "N", "age_start", "age_end") := list(379 - 350, 1024 - 815, 1, 9)]
  oncho[Master_UID == "5758_SR", c("had_oncho_poly", "N", "age_start", "age_end") := list(54 - 54, 1525 - 1249, 1, 9)]
  
  oncho <- oncho[!(nid == "332893" & age_start == 15 & age_end == 99)]
  
  ## Check for cohorts reporting both ss and nod, and drop the latter
  dx_by_group <- as.data.table(aggregate(diagnostic ~ cohort_id, data = oncho[order(diagnostic), ], unique))
  multiple_dx <- oncho[cohort_id %in% dx_by_group[diagnostic %in% c("c(\"nod\", \"ss\")"), cohort_id]]
  
  oncho <- rbind(oncho[!(cohort_id %in% multiple_dx$cohort_id)], oncho[(cohort_id %in% multiple_dx$cohort_id) & diagnostic == "ss"])
  
  ## Now apply crosswalks
  ### Setup
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(MASS, boot, data.table, ggplot2)
  path <- <<<< FILEPATH REDACTED >>>>
  library(fda, lib.loc = path) ### NOTE: Change path as necessary
  library(faraway, lib.loc = path) ### NOTE: Change path as necessary
  library(subplex, lib.loc = <<<< FILEPATH REDACTED >>>>)
  library(BayesianTools, lib.loc = <<<< FILEPATH REDACTED >>>>)
  library(matrixStats, lib.loc = path)
  
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  source(<<<< FILEPATH REDACTED >>>>)
  
  ### Apply diagnostic crosswalk
  ## Define needed functions
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
  
  calculatePrevalence_dx <- function(AgeSpline, scale, diag_fe) {
    return((ilogit(logit(0.00001) + AgeSpline + scale + diag_fe))) # Specify near-zero intercept
  }
  
  ## Calculate the prevalence-age curve in logit space; can put something other than splines in here
  applySpline <- function(vec, Basis, avec) {
    AgeSpline <- fd(vec, Basis)
    return(predict(AgeSpline, avec))
  }
  
  ### Apply constraints to shape of spline
  applyBasisConstraints <- function(betas) {
    alphas <- copy(betas)
    return(alphas)
  }
  
  #### Calculate prevalence
  calculatePrevalence <- function(AgeSpline, scale, diag_fe) {
    return((ilogit(logit(0.00001) + AgeSpline + scale + diag_fe))) # Specify near-zero intercept
  }
  
  crosswalk_data <- copy(oncho)
  
  country_years <- unique(crosswalk_data[, c("country", "year")])
  
  age_ids <- c(49:142) # Age 1 = ID 49, Age 94 = ID 142; get_population is not returning single_year_age for ages 95-99
  gbd_round_id <- 6 ### 2019
  
  location_hierarchy <- get_location_metadata(location_set_id = 35, gbd_round_id = gbd_round_id, decomp_step = "iterative")
  
  country_years <- merge(country_years, location_hierarchy[, c("ihme_loc_id", "location_id")], by.x = "country", by.y = "ihme_loc_id", all.x = TRUE)
  
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
  
  ###### Create a new unique ID
  crosswalk_data$new_unique_id <- 1:nrow(crosswalk_data)
  
  ## Merge age distribution with crosswalk data
  newData <- merge(crosswalk_data, pop_age_merged_wide, by = c("location_id", "year"), all.x = TRUE)
  names(newData)[names(newData) == "age_prop.<1 year"] <- "age_prop.0" # Rename age_prop.<1 year
  
  # Retrieve crosswalk model objects
  load(<<<< FILEPATH REDACTED >>>>)
  
  freq_coefs <- oncho_crosswalk_model$Optim$par
  
  BREAKS_ss <- oncho_crosswalk_model$BREAKS_ss
  BREAKS_nod <- oncho_crosswalk_model$BREAKS_nod
  SplineNum <- length(BREAKS_ss) + 2
  
  Basis_ss <- create.bspline.basis(rangeval = c(0, 94.5), breaks = BREAKS_ss)
  Basis_nod <- create.bspline.basis(rangeval = c(0, 94.5), breaks = BREAKS_nod)
  
  avec <- oncho_crosswalk_model$avec
  
  diag <- "ss"
  
  data_to_adjust <- copy(newData)
  
  # Setup
  age_threshold <- 99
  age_threshold_lower <- 0
  
  ss_fe <- freq_coefs[length(freq_coefs)]
  
  data_to_adjust[age_end > 94, age_end := 94]
  
  #### Calculate negative log likelihood
  NegLL_dx <- function(vec, site) {return(-LL_dx(vec, site))}
  
  #### Calculate likelihood
  LL_dx <- function(vec, site) {
    alphas <- applyBasisConstraints_dx(freq_coefs[(SplineNum + 1):(SplineNum * 2)])
    AgeSpline <- applySpline(alphas, Basis_nod, avec)
    
    Scales <- as.numeric(vec)
    cohort_number <- site
    
    #### Calculate likelihood for each study-site
    ll <- 0
    AgeLoc <- data.table("age" = avec, "prev" = calculatePrevalence_dx(AgeSpline, Scales, 0))
    colnames(AgeLoc)[2] <- "prev"
    
    #### Constrain prevalence to be flat beyond a threshold age
    AgeLoc[age > age_threshold, prev := AgeLoc[age == age_threshold, prev]]
    
    locs <- which(USiteTemp == cohort_number)
    
    #### For each location, calculate the likelihood of the data given the prevalance-by-age curve formed from the parameters
    for (j in 1:length(locs)) {
      Start <- data_to_adjust$age_start[locs[j]]
      End <- data_to_adjust$age_end[locs[j]]
      Test <- data_to_adjust$N[locs[j]]
      Pos <- round(data_to_adjust$had_oncho_poly[locs[j]], 0)
      PGroup <- as.vector(as.matrix(data_to_adjust[locs[j], which(colnames(data_to_adjust) == paste0("age_prop.", Start)):which(colnames(data_to_adjust) == paste0("age_prop.", End))]))
      P <- sum(AgeLoc[(age - 0.5) %in% Start:End, prev] * PGroup / sum(PGroup))
      
      ll <- ll + dbinom(Pos, Test, P, log = TRUE)
    }
    
    return(ll)
  }
  
  #### Calculate negative log likelihood
  NegLL <- function(vec, site, diag_fe) {return(-LL(vec, site, diag_fe))}
  
  #### Calculate likelihood
  LL <- function(vec, site, diag_fe) {
    alphas <- applyBasisConstraints(freq_coefs[1:SplineNum])
    AgeSpline <- applySpline(alphas, Basis_ss, avec)
    
    Scales <- as.numeric(vec)
    cohort_number <- site
    
    #### Calculate likelihood for each study-site
    ll <- 0
    AgeLoc <- data.table("age" = avec, "prev" = calculatePrevalence(AgeSpline, Scales, diag_fe))
    colnames(AgeLoc)[2] <- "prev"
    
    #### Constrain prevalence to be flat beyond a threshold age
    AgeLoc[age > age_threshold, prev := AgeLoc[age == age_threshold, prev]]
    
    locs <- which(USiteTemp == cohort_number)
    
    #### For each location, calculate the likelihood of the data given the prevalance-by-age curve formed from the parameters
    for (j in 1:length(locs)) {
      Start <- data_to_adjust$age_start[locs[j]]
      End <- data_to_adjust$age_end[locs[j]]
      Test <- data_to_adjust$N[locs[j]]
      Pos <- round(data_to_adjust$had_oncho_poly[locs[j]], 0)
      PGroup <- as.vector(as.matrix(data_to_adjust[locs[j], which(colnames(data_to_adjust) == paste0("age_prop.", Start)):which(colnames(data_to_adjust) == paste0("age_prop.", End))]))
      P <- sum(AgeLoc[(age - 0.5) %in% Start:End, prev] * PGroup / sum(PGroup))
      
      ll <- ll + dbinom(Pos, Test, P, log = TRUE)
    }
    
    return(ll)
  }
  
  #### Grab unique study-site levels
  USiteTemp <- data_to_adjust$cohort_id
  USite <- as.data.table(unique(USiteTemp))
  setnames(USite, "V1", "cohort")
  
  ## Calculating scaling factor for each cohort
  for (i in 1:nrow(USite)) {
    message(paste0("Calculating scaling factor for cohort ", USite[i, cohort], ", ", i, " of ", nrow(USite)))
    if ((sum((data_to_adjust[cohort_id == USite[i, cohort], had_oncho_poly])) > 0) & data_to_adjust[cohort_id == USite[i, cohort], had_oncho_poly] != data_to_adjust[cohort_id == USite[i, cohort], N]) {
      if (data_to_adjust[cohort_id == USite[i, cohort], diagnostic][1] == "nod") {
        message("nod")
        Optim <- optim(0, NegLL_dx, control = list(maxit = 10000), method = "L-BFGS-B", site = USite[i, cohort])
        USite[i, "scale"] <- Optim$par
        message(Optim$par)
      } else {
        message("ss")
        Optim <- optim(0, NegLL, control = list(maxit = 10000), method = "L-BFGS-B", site = USite[i, cohort], diag_fe = ss_fe)
        USite[i, "scale"] <- Optim$par
        message(Optim$par)
      }
    } else {
      USite[i, "scale"] <- NA
    }
  }
  
  ## Apply scaling factor to calculate all-age prevalence
  if (nrow(data_to_adjust) > 0) {
    for (i in 1:nrow(USite)) {
      message(paste0(i, " of ", nrow(USite)))
      current <- data_to_adjust[cohort_id == USite[i, cohort],]
      
      #### Retrieve parameter values
      alphas <- applyBasisConstraints(freq_coefs[1:SplineNum])
      
      #### Calculate age-prevalence (logit) curve
      AgeSpline <- applySpline(alphas, Basis_ss, avec)
      
      #### Retrieve scaling factor
      scaling <- USite[i, scale]
      if (!is.na(scaling)) {
        #### Calculate scale-adjusted age curve
        AgeLoc_scaled <- data.table("age" = avec, "prev" = calculatePrevalence(AgeSpline, scaling, ss_fe))
        colnames(AgeLoc_scaled)[2] <- "prev"
        
        #### Constrain prevalence to be flat beyond a threshold age
        AgeLoc_scaled[age > age_threshold, prev := AgeLoc_scaled[age == age_threshold, prev]]
        
        #### Calculate final all-age prevalence
        Start <- 0
        End <- 94
        PGroup <- as.vector(as.matrix(current[, which(colnames(current) == paste0("age_prop.", Start)):which(colnames(current) == paste0("age_prop.", End))]))
        P <- sum(AgeLoc_scaled[(age - 0.5) %in% Start:End, prev] * PGroup / sum(PGroup))
        
        USite[i, "adjusted_prev"] <- P
        rm(scaling)
        rm(P)
      } else {
        USite[i, "adjusted_prev"] <- NA
      }
    }
  }
  
  #### Merge adjusted prevalence to newData
  newData <- merge(newData, USite[, c("cohort", "adjusted_prev")], by.x = "cohort_id", by.y = "cohort", all.x = TRUE)
  newData[is.na(adjusted_prev), adjusted_prev := had_oncho_poly / N]
  
  #### Now collapse rows with multiple age bins
  newData$counter <- 1
  rows_by_cohort <- aggregate(counter ~ cohort_id, data = newData, sum)
  
  unique_cohorts <- newData[J(unique(cohort_id)), mult="first"]
  summed_sample_size <- merge(merge(merge(newData[, sum(N), by = cohort_id], newData[, min(age_start), by = cohort_id]), newData[, max(age_end), by = cohort_id]), aggregate(Master_UID ~ cohort_id, data = newData, paste, collapse = "_"))
  summed_cases <- merge(merge(merge(newData[, sum(had_oncho_poly), by = cohort_id], newData[, min(age_start), by = cohort_id]), newData[, max(age_end), by = cohort_id]), aggregate(Master_UID ~ cohort_id, data = newData, paste, collapse = "_"))
  colnames(summed_cases) <- c("cohort_id", "had_oncho_poly_collapsed", "age_start_collapsed", "age_end_collapsed", "Master_UID_collapsed")
  colnames(summed_sample_size) <- c("cohort_id", "N_collapsed", "age_start_collapsed", "age_end_collapsed", "Master_UID_collapsed")
  summed <- merge(summed_cases, summed_sample_size, by = c("cohort_id", "age_start_collapsed", "age_end_collapsed", "Master_UID_collapsed"))
  unique_cohorts <- merge(unique_cohorts, summed, by = "cohort_id", all.x = TRUE)
  unique_cohorts$had_oncho_poly_dx_age_adjusted <- unique_cohorts$adjusted_prev * unique_cohorts$N_collapsed
  
  ### Set up final columns
  unique_cohorts$N <- unique_cohorts$N_collapsed
  unique_cohorts$had_oncho_poly <- round(unique_cohorts$had_oncho_poly_dx_age_adjusted, 0)
  unique_cohorts$age_start <- unique_cohorts$age_start_collapsed
  unique_cohorts$age_end <- unique_cohorts$age_end_collapsed
  unique_cohorts$Master_UID <- unique_cohorts$Master_UID_collapsed
  unique_cohorts$weight <- 1
  
  oncho <- copy(unique_cohorts)
  oncho <- oncho[, grep("age_prop.", colnames(oncho), invert=TRUE), with=FALSE]
  oncho$age_midpoint <- NA
  oncho$crosswalk_weight <- NA
}

oncho[point == 0 & !is.na(latitude) & !is.na(longitude), point := 1]
oncho[Master_UID == "16031_SR", shapefile := "lf_gadm28_adm2"]

## Ad hoc deduplication
oncho <- oncho[!(nid %in% c(203861))]
oncho[source == "REMO", source := "APOC"]

write.csv(oncho, <<<< FILEPATH REDACTED >>>>, row.names = F)

oncho <- rbindlist(list(oncho, sero_data), use.names = TRUE, fill = TRUE)
write.csv(oncho, <<<< FILEPATH REDACTED >>>>, row.names = F)

### Output crosswalk training data set
xwalk <- oncho[use_in_xwalk == TRUE]

write.csv(xwalk, <<<< FILEPATH REDACTED >>>>, row.names = F)
