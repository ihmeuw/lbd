#----HEADER-------------------------------------------------------------------------------------------------------------
# Date:    November 2018
# Purpose: Tabulate unit record data extracted with UbCov; saves dataset by NID
#               LBD - Collapse by lat/long or admin unit
# Inputs:  nid  ->  NID from dataset to tabulate
#***********************************************************************************************************************

source(paste0(vaccines_repo, "FILEPATH/process.r"))
#----START FUNCTION-----------------------------------------------------------------------------------------------------
username        <- 'USERNAME'
### load functions
core_repo      <- paste0("FILEPATH"); setwd(core_repo)
commondir      <- sprintf("FILEPATH")

# source(paste0(core_repo, "mbg_central/graph_data_coverage.R"))

### load extra packages
package_list   <- c(t(read.csv(sprintf('%s/package_list.csv', commondir), header=FALSE)))
source('FILEPATH/setup.R')
source("FILEPATH/geocodebook_functions.R")
invisible(mbg_setup(package_list=package_list, repos=core_repo))

### specific function loads (avoid namespace conflicts)
str_match <- stringr::str_match


#----HELPER FUNCTIONS------------------------------------------------------------------------------------------------

# Given pair of antigens, return dataset with ratio of antigens created at cluster-level.
# Ratios not generated where denominator (ratio_antigen_2) is 0.
add_ratios <- function(data, ratio_antigen_1, ratio_antigen_2) {
  vaccines_in_dataset    <- unique(data$me_name)
  ratio_antigen_1_exists <- ratio_antigen_1 %in% vaccines_in_dataset
  ratio_antigen_2_exists <- ratio_antigen_2 %in% vaccines_in_dataset
  
  if (ratio_antigen_1_exists & ratio_antigen_2_exists) {
    
    ### Reshape long to wide in order to add ratios, then wide back to long
    # value   me_name   (Long to Wide)   DPT3   HepB3  HepB3_DPT3_Ratio  (Wide to Long)   value   me_name   
    #   2      DPT3          --->          2      1          .5              --->          2      DPT        
    #   1      HepB3         --->                                            --->          1      HepB       
    #                                                                                     .5      HepB3_DPT3_Ratio 

    ### Reshape long to wide
    # In reshape, LHS of formula is determined by all columns excluding me_name and value
    lhs_names <- names(data)[!names(data) %in% c("me_name", "value")]
    lhs <- paste(lhs_names, collapse = " + ")
    rhs <- "me_name"
    lhs_rhs <- paste(c(lhs, rhs), collapse = " ~ ")
    
    data <- data.table(data.table::dcast(data, eval(parse(text = lhs_rhs))))
    
    ### Add ratios
    # Set ratio column name
    ratio_column_name <- gsub("vacc_", "", paste0(ratio_antigen_1, "_", ratio_antigen_2, "_ratio"))
    data[, (ratio_column_name) := numeric() ]
    # Set ratio ignoring rows where denominator is
    data[!is.na(get(ratio_antigen_2)) & get(ratio_antigen_2) != 0 & !is.na(get(ratio_antigen_1)), (ratio_column_name) := get(ratio_antigen_1) / get(ratio_antigen_2)]
    
    ### Reshape wide back to long
    vax_cols     <- names(data)[grep("vacc|ratio", names(data))]
    id_variables <- names(data)[!names(data) %in% vax_cols]
    data <- melt.data.table(data, id.vars=id_variables, measure.vars = vax_cols, variable.name="me_name", variable.factor = F, value.name = "value")
    
    # Drop NAs resulting from reshaping process  
    data <- data[!is.na(value), ]
  }
  return(data)
}


### read function
read_add_name_col <- function(file) {
  rn  <- gsub(".csv", "", file, ignore.case=TRUE)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df  <- fread(file)
  df$survey_series <- svy
  return(df)
}

#----TABULATE--------------------------------------------------------------------------------------------------------

tabulate_lbd <- function(nid, crop_to_region="all", collapse_method="kish") {
  
  ### Prepare for geomatching
  ## Get all geo codebooks and package them together
  files <- list.files("FILEPATH", pattern=".csv$", ignore.case=TRUE, full.names=TRUE)
  files <- grep("Copy|together|linkage|IPUMS", files, value=TRUE, invert=TRUE) #IPUMS is handled separately
  ### Transition to geo_codebook database
  geo <- get_geocodebooks(nids = nid)
  setnames(geo, "iso3", "ihme_loc_id")
  geo$ihme_loc_id <- substr(geo$ihme_loc_id, 0, 3)
  
  if(record_data_drop){
    filter_table <- load_filter_table()
  } else {
    filter_table <- NULL
  }
  
  
  ##################################################################################
  ########################### GEOGRAPHY MATCHING ###################################
  ##################################################################################
  
  ### get survey data
  message("Prep for tabulation")
  dataset <- prep_for_tabulation(nid, team="lbd", vaccines.=c(vaccines, "rotac", "dpt3_timeliness_ratio"), filter_table) #In process.r, checks for key variables before going on to tabulation 
  
  # Data from prep for tabulation is returned as a list. The first element is T/F indicating whether the data is usable,
  # the 2nd element is the data, and the optional 3rd element is a "filter table" for recording data drops at each filtering point
  if (dataset[[1]] != FALSE){ 
    data <- dataset[[2]] %>% as.data.table 
    
    if(record_data_drop){
      filter_table <- dataset[[3]]
    }
    
    
    for(n in c("pweight", "psu", "strata")){   #If we are missing these variables, add them as NAs to make it through reshapes/tabulations. Only drop if these are missing for polygons later
      if(!n %in% names(data)){
        data[ ,(n):=NA]
      }
    }
    if(max(data[ ,year_id]) >= 1999){  #year id = birth year. Born in 1999 = vaccine coverage cohort is 2000
      
      ### Remove any characters after ISO3 code
      data$ihme_loc_id <- as.character(data$ihme_loc_id)
      data$ihme_loc_id[data$ihme_loc_id != "KOSOVO"] <- substr(data$ihme_loc_id, 0, 3)[data$ihme_loc_id != "KOSOVO"]
      data$ihme_loc_id[data$ihme_loc_id %in% c("KOSOVO", "RKS")] <- "XKO"  ###This is Kosovo's iso3 in the codebooks
      
      data$ihme_loc_id[data$nid == 7688] <- "LBN" #fix for Palestinians living in Lebanon survey, can drop this if extraction is fixed
      
      # data[, nid := as.character(nid)]
      data[, geospatial_id := as.character(geospatial_id)]
      geo[, geospatial_id := as.character(geospatial_id)]
      
      #merge with geo codebooks. this could fail if geospatial IDs have been changed, or if the iso code is wrong (subnational iso3s, wrong country, mistyped iso3)
      
      data <- merge(geo, data, by=c("nid", "geospatial_id", "ihme_loc_id"), all.x=FALSE, all.y=TRUE, allow.cartesian=TRUE)  
      
      ### rename
      rename_table <- data.table(rbind(c("ihme_loc_id", "country"),
                                       c("nid", "svy_id"),
                                       c("geospatial_id", "geo_id"), 
                                       c("lat", "latitude"),
                                       c("long", "longitude"), 
                                       c("loc_name1", "admin1_name"),
                                       c("loc_code1", "loc_code_admin1"),
                                       c("admin_level1", "admin_level1"),
                                       c("point", "point")))
      names(rename_table) <- c("old", "new")
      invisible(lapply(1:nrow(rename_table), function(x) check_and_rename_col(data, rename_table$old[x], rename_table$new[x])))
      
      ###SPECIFIC RECODES
      
      fix <- copy(data)
      fix[, num_clusters:= uniqueN(geo_id), by = c("svy_id", "country", "year_id", "survey_name", "file_path")]
      
      fix <- data.table(subset(fix, is.na(shapefile) & (is.na(latitude) | is.na(longitude)) ))   
      fix_collapse <- fix[, uniqueN(geo_id), by = c("svy_id", "country", "year_id", "survey_name", "file_path", "num_clusters")]
      names(fix_collapse)[names(fix_collapse) == "V1"] = "num_unmatched_clusters"
      
      fix_total <- fix[, unique(geo_id), by=.(svy_id, country, year_id, survey_name, file_path)]
      names(fix_total)[names(fix_total) == "V1"] = "geo_id"
      
      ##################################################################################
      ########################### END GEOGRAPHY MATCHING ###############################
      ##################################################################################
      
      # recode shapefile
      data[shapefile == "", shapefile := NA]

      # Placeholder value to sum over
      data[, N := 1]
      
      # Ensure column formats correct
      if (class(data$latitude) != "numeric") data[, latitude := as.numeric(as.character(latitude))]
      if (class(data$longitude) != "numeric") data[, longitude := as.numeric(as.character(longitude))]
      
      ### column for missingness
      missing_shapefiles <- list()
      
      ### COLLAPSE POLYS AND POINTS ###########################################
      
      # Get shapefiles into character format (in case factor)
      data$shapefile <- as.character(data$shapefile)
      
      # Recode shapefile for points
      
      # Figure out what shapefiles are missing and exclude a priori
      # (Otherwise will break resample_polygons)
      shapefile_list <- unique(data$shapefile[!is.na(data$shapefile)])
      shapefile_dir  <- "FILEPATH"
      shapefiles_available <- list.files(shapefile_dir, ".shp$") %>% gsub(".shp", "", .)
      
      no_shapefile <- shapefile_list[!(shapefile_list %in% shapefiles_available)]
      
      if (length(no_shapefile) > 0) {
        
        message("The following shapefiles are missing from the shapefile directory. Associated data will be dropped:")
        message(paste(paste0("  ", no_shapefile, ".shp"), collapse = " \n"))
        data <- data[!(shapefile %in% no_shapefile)]
        missing_shapefiles[[as.character(nid)]] <- no_shapefile
        
      }
      
      ### geo missingness log file
      geo_missingness <- nrow(data[is.na(latitude) & is.na(shapefile)]) / nrow(data)
      vet_log <- fread(file.path(extraction_root, "log/details", paste0(nid, ".csv")))
      vet_log[, missingness_geomatch := geo_missingness]
      write.csv(vet_log, file.path(extraction_root, "log/details", paste0(nid, ".csv")), row.names=FALSE)
      
      # Fix missing points
      # Some nids were being passed through processing wrapper that still had is.na(point)
      data[is.na(point) & !is.na(latitude) & !is.na(longitude), point := 1]
      data[is.na(point) & !is.na(shapefile) & !is.na(pweight), point := 0]
      
      # **DATA DROP**
      data <- data[!(is.na(latitude) & is.na(shapefile))] # geo_drop
      data <- data[!(is.na(latitude) & is.na(pweight))]   # pweight_drop if polygon data does not have pweight
      data <- data[!(is.na(year_id))]                     # age_drop
      data <- data[!(is.na(value))]                       # missing vaccination_data
      
      # First, check if there are any polygons
      if (nrow(data[point==0]) == 0) {
        
        #If so, skip ahead
        message("No polygon data found - moving to next cycle")
        df_pointpoly <- copy(data)
        keep_vars <- c("svy_id", "survey_name", "country", "point", "svy_year", 
                       "location_code", "shapefile", "year_id", 
                       "psu", "latitude", "longitude", "outlier", "me_name", "value", "N")
        df_pointpoly <- subset(df_pointpoly, select = names(df_pointpoly) %in% keep_vars)
        
        df_pointpoly[ ,N_obs := 1]
        
        # Collapse the point data to clusters
        df_pointpoly <- df_pointpoly[, lapply(.SD, sum), 
                                     by = .(svy_id, survey_name, country, point, 
                                            location_code, shapefile, year_id, 
                                            psu, latitude, longitude, me_name),
                                     .SDcols = c("value", "N", "N_obs")]
        
      } else {
        
        # Collapse the polygon data
        data[, cluster_id := NULL]
        data[, strata := as.numeric(strata)]
        
        DHS_subset <- copy(data)
        DHS_subset <- DHS_subset[survey_name=="MACRO_DHS" & !is.na(shapefile)]
        
        df_point <- data[point == 1] %>% copy
        df_poly  <- data[point == 0] %>% copy
        
        df_poly$pweight <- as.numeric(df_poly$pweight)
        
        
        # Collapse polygons -------------------------------------------------------------
        
        if (collapse_method == "survey") {
          
          
          # need to work this in
          
        } else if (collapse_method == "kish") {
          
          # Collapse using the Kish approximation to generate an effective sample size (N)
          # Keeps the original sample size in N_obs
          
          df_poly <-   df_poly[, c(list(N_obs = sum(N),
                                        SSW = sum(pweight),
                                        N = sum(pweight) ^ 2 / sum(pweight ^ 2)),
                                   lapply(.SD, FUN = weighted.mean, w = pweight)), 
                               by = list(svy_id, survey_name, country, point, 
                                         location_code, shapefile, year_id, me_name),
                               .SDcols = "value"]
          
          # Convert these back to counts
          df_poly[, ("value") := lapply(.SD, '*', N), .SDcols = "value"]                            
          
        }
        
        # Collapse the point data to clusters
        if(nrow(df_point) > 0){
          df_point[ ,N_obs := 1]
          df_point <- df_point[, lapply(.SD, sum), 
                               by = .(svy_id, survey_name, country, point, 
                                      location_code, shapefile, year_id, pweight, #need to figure out SSW
                                      strata, psu, latitude, longitude, me_name),
                               .SDcols = c("value", "N", "N_obs")]
        }
        # Harmonize names & combine
        #drop_vars <- c("strata", "pweight")
        #df_point <- subset(df_point, select = !(names(df_point) %in% drop_vars))
        
        df_pointpoly <- rbind(df_point, df_poly, fill=TRUE)
        
      }
      
      # Add a weights column - this will be resampled in resample_polygons
      df_pointpoly$weight <- 1
      
      
      
      # Add Ratios -------------------------------------------------------------
      
      df_pointpoly[, value := as.double(value)]
      
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_mcv2",  "vacc_mcv1")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_pcv3",  "vacc_dpt3")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_rotac", "vacc_dpt3")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_hib3",  "vacc_dpt3")
      df_pointpoly <- add_ratios(df_pointpoly, "vacc_hepb3", "vacc_dpt3")
      
      
      
      # End Ratios -------------------------------------------------------------
      
      
      ### Final renaming stuff
      #setnames(df_pointpoly, "svy_year", "original_year")
      
      df_pointpoly[, country := tstrsplit(country, "_")[1]]
      
      write.csv(df_pointpoly, paste0(extraction_root, "/tabulated/lbd/", nid, ".csv"), row.names=FALSE)
      
      if(record_data_drop & !is.null(filter_table)){
        filter_table[, survey_in_year_range := TRUE]
        write.csv(filter_table, file = paste0(filter_table_path, nid, ".csv"), row.names = FALSE)
      }
      
      message("Finished LBD tabulation")
      
    } else {
      
      df_pointpoly <- NULL
      message(paste0(nid, ": All data is pre-2000, skipping tabulation"))
      if (paste0(nid, ".txt") %in% list.files(file.path(extraction_root, "log/datestamp", date))) {
        cat("\nDataPre2000\n", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
      }
    }
    
  } else {
    df_pointpoly <- NULL
    message(paste0(nid, ": Tabulation prep failed"))
    
    if(record_data_drop){
      
      message("\n **Survey missing crucial indicator. Ending** \n")
      filter_table[, c("missing_indicator", "record_filter_successful") := list(TRUE)]
      filter_table[, setdiff(names(filter_table), c("survey_id", "missing_indicator", 
                                                    "survey_in_year_range", "outlier", "record_filter_successful")) := NA_integer_]
      
      record_data_drop <- FALSE
      write.csv(filter_table, file = paste0(filter_table_path, nid, ".csv"), row.names = FALSE)
    }
  }
}


#***********************************************************************************************************************
