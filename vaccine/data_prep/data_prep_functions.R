# HEADER ------------------------------------------------------------------
# data_prep_functions.R
# Purpose: Functions for use in prepping vaccine coverage data for 
#          MBG models
#**************************************************************************

# Geopositioning function - geoposition survey data using a codebook that matches
# each PSU to either a polygon within a shapefile or a lat/long. 

geoposition_data <- function(survey_type, df, country_specific = F, error_dir) {
  
  # INPUTS: 
  #          survey_type: "MACRO_DHS", etc
  #          df: data frame (should only contain data of survey_type)
  #          country_specific: if the data source is country-specific (T) or general (F)
  #          error_dir: where to output any error messages
  #
  # OUTPUTS: list with two items:
  #          [[1]] = DT of geopositioned data
  #          [[2]] = DT summarizing failed matches (surveys without match)
  #
  #          also writes a series of files in `error_dir`:
  #               dups_[svy_name].csv contains duplicated nid/psu combos for [svy_name]
  
  # 1. Create geographies database -----------------

  if (country_specific == T) {
   geo_file <- "<<<< FILEPATH REDACTED >>>>/COUNTRY_SPECIFIC.csv"
   geographies_raw <- read.csv(geo_file, stringsAsFactors = F) %>% as.data.table
   geographies_raw$svy_type <- geographies_raw$file_path %>% as.character %>%
                                                                gsub("\\\\", "/", .) %>%
                                                                gsub("<<<< FILEPATH REDACTED >>>>|/[0-9].*", "", .)

   # Subset to only survey of interest
   geographies_raw <- geographies_raw[svy_type == survey_type]       

  } else if (country_specific == F) {
  
    geo_file <- paste0("<<<< FILEPATH REDACTED >>>>/", survey_type, ".csv")
    geographies_raw <- read.csv(geo_file, stringsAsFactors = F) %>% as.data.table
  
  }

  # rename cluster_number as psu
  names(geographies_raw)[names(geographies_raw) == "cluster_number"] <- "psu"

   # Ensure numeric
  df$nid <- as.numeric(df$nid)
  df$psu <- as.numeric(df$psu)
  geographies_raw$nid <- as.numeric(geographies_raw$nid)
  geographies_raw$psu <- as.numeric(geographies_raw$psu)

  # Drop if NID or PSU is missing - misspecified in codebooks
  geographies_raw <- geographies_raw[(!is.na(nid) & !is.na(psu)),]
  geographies_raw <- unique(geographies_raw)

  # custom fixes for nids, etc.
  if (survey_type == "UNICEF_MICS") {
    # recode the SSD portion of (UNICEF/MICS 2000) from 12243 to 12232
    geographies_raw[iso3 == "SSD" & nid == 12243, nid := 12232]
  }
  
  # remove unwanted fields
  keepvars <- c("nid", "psu", "uncertain_point", "buffer", 
                "location_name", "location_code",
                "lat", "long", "file_path",
                "admin_level", "shapefile")
  
  geographies <- subset(geographies_raw, select = keepvars)
 
  # 2. Match data with geographies ------------------

  # Remove any duplicate NID/PSU combinations
  n_pre <- nrow(geographies)
  dup_rows <- geographies[duplicated(geographies[,c("nid", "psu")]) | duplicated(geographies[,c("nid", "psu")], fromLast = T)]
  geographies <- geographies[!duplicated(geographies[, c("nid", "psu")])]
  n_post <- nrow(geographies)

  if (n_pre > n_post) {
    message(paste0("Warning: dropping ", n_pre - n_post, " rows from geography codebook with duplicate nid/psu combos."))
    message("By default using the first row, but check assumptions")
    message(paste0("See ", error_dir, "dups_",survey_type,".csv"))
    write.csv(dup_rows, file=paste0(error_dir, "/dups_",survey_type,".csv"))
  }

  # Set up custom merges for specific data sets
  nid_match_admin1 <- c(12243, # UNICEF_MICS SDN 2000
                        8932,  # UNICEF_MICS MMR 2000
                        1404   # UNICEF_MICS BWA 2000
                       )

  # Split out custom match data
  df_match_admin1 <- df[nid %in% nid_match_admin1]
  df <- df[!(nid %in% nid_match_admin1)]

  # Do the main merge
  matched_data <- merge(df, geographies, by = c("nid", "psu")) 

  # Do the custom merge

  if (nrow(df_match_admin1) > 0) {

    # Set up custom geographies database
    geo_admin1 <- geographies
    geo_admin1[,psu := NULL]
    geo_admin1 <- geographies[nid %in% nid_match_admin1,]

    # Harmonize case & trim white space    
    geo_admin1$matchvar <- tolower(geo_admin1$location_name) %>% trimws
    df_match_admin1$matchvar <- tolower(df_match_admin1$admin_1) %>% trimws

    # Custom merge on new matchvar
    matched_admin1 <- merge(df_match_admin1, geo_admin1, by = c("nid", "matchvar"))

    # Drop matchvar
    matched_admin1[, matchvar := NULL]

    # Rbind with the main merge results
    matched_data <- rbind(matched_data, matched_admin1)
  }

  message(paste0("Dropping ", nrow(df) + nrow(df_match_admin1) - nrow(matched_data),
                  " rows during merge"))
  
  # Replace blank shapefiles, lat, long with NA for consistency
  matched_data[shapefile == "", ]$shapefile <- NA
  matched_data[long == "", ]$long <- NA
  matched_data[lat == "", ]$lat <- NA
  
  # 3. Generate additional variables, fix names ------
  
  # Drop pointpoly latitude/longitude
  # Use geography codebook lat/longs instead
  
  dropvars <- c("latitude", "longitude")
  matched_data <- subset(matched_data, select = !names(matched_data) %in% dropvars)
  
  # point = 1 if lat/long, 0 if polygon only, NA if not matched
  matched_data[!is.na(lat) & !is.na(long), point := 1]
  matched_data[!is.na(shapefile) & is.na(lat) & is.na(long), point := 0]
  
  # 4. Drop unmatched rows ----------------------------
    
  # Drop unmatched rows drop!

  message(paste0("Dropping ", nrow(matched_data[is.na(point)]),
                 " rows in geography database but not geolocated"))
  matched_data <- matched_data[!is.na(point)]
  
  # Drop any extraneous shapefiles
  
  matched_data[point == 1, shapefile := NA]
  
  message("Geopositioned ", nrow(matched_data), 
          " out of ", nrow(df), " rows (", 
          round(nrow(matched_data)/nrow(df) * 100, 2), 
          " %)")

  # 5. Find unmatched surveys
    # Create list of concatenated nid/psu
    nid_match <- unique(geographies_raw$nid)
    nidpsu_match <- unique(paste0(geographies_raw$nid, "_", geographies_raw$psu))
    umd <- subset(df, select = c("nid", "psu", "survey_name", "ihme_loc_id", "year_start"))

    # First, figure out which clusters don't have an entry
    umd[, nidpsu := paste0(nid, "_", psu)]
    umd[!(nidpsu %in% nidpsu_match), no_nidpsu := 1]
    
    # Next, figure out which don't have any nid entry
    umd[!(nid %in% nid_match), no_nid := 1]

    # Drop the matched ones
    umd <- umd[!is.na(no_nidpsu) | !is.na(no_nid), ] %>% unique

    # Figure out what's affected
    umd[no_nidpsu == 1 & no_nid == 1, affected := "Entire survey"]
    umd[no_nidpsu == 1 & is.na(no_nid), affected := "Certain PSUs only"]
    umd <- unique(umd)

    # Create psu column
    psu_table <- umd[, .(psus = paste(psu, collapse = ",")), by = c("nid")]
    umd <- merge(umd, psu_table)
    umd[affected == "Entire survey", psus := NA]

    keepvars <- c("nid", "survey_name", "ihme_loc_id", "year_start", "affected", "psus")
    umd <- subset(umd, select = keepvars) %>% unique

    # replace for those where entire survey missing
    umd[affected == "Entire survey", psus := NA]

  # 6. Return results
   
  return(list(matched_data, umd))
  
}

# Add a vaccine to the list of vaccines to process in the data preparation process
add_vaccine <- function(prefix, title, doses, age_min = 12, age_max = 59, vaccine_list = NULL) {
   
    # Note ages in months
    cond_doses <- 1:(doses - 1)
    cond_vaccines <- paste0(prefix, cond_doses, "_cond")
    all_doses <- 0:doses
    all_vaccines <- paste0(prefix, all_doses, "_cov")
    all_titles <- paste0(title, all_doses)
    age_range <- c(age_min, age_max)
    
    return_list <- list(prefix = prefix,
                        title = title,
                        cond_doses = cond_doses, 
                        cond_vaccines = cond_vaccines, 
                        all_doses = all_doses, 
                        all_vaccines = all_vaccines, 
                        all_titles = all_titles, 
                        age_range = age_range)
    
    if(is.null(vaccine_list)) vaccine_list <- list()
    vaccine_list[[prefix]] <- return_list
   
    return(vaccine_list)
  }

# Convenience function to print a header to the console
print_header <- function() {
  message("\n#######################################################################")
  message("#######################################################################")
}