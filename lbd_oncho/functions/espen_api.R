###
### ESPEN API DATA / MAPS CALLS
###    https://admin.espen.afro.who.int/docs/api
###
### Written to automate the process of pulling data from the ESPEN Portal through the ESPEN API
###
### Packages needed: httr, jsonlite, lubridate, sf, geojsonsf, geojson, data.table
###################################################################################################

###################################################################################################
### API CALLS #####################################################################################

### ESPEN API DATA CALL ###########################################################################
#
### FUNCTION ###
# Interacts with ESPEN API to pull specific data for the user
# Currently saves data to two locations with the current date :
#   <<<< FILEPATH REDACTED >>>>
#   <<<< FILEPATH REDACTED >>>>
#
### PARAMETERS ###
#'@disease - which disease to pull data on { lf | oncho | loa | sch | sth | trachoma }
#'@level - which geographic size level to pull data on. { iu | sitelevel }

espenAPIData <- function(disease, level) {
  
  ### Base inputs ### ----
  options(stringsAsFactors = F)
  diseases <- c("lf", "oncho", "loa", "sch", "sth", "trachoma")
  levels <- c("sitelevel", "iu")
  urlBase <- "https://admin.espen.afro.who.int/api/data?"
  num <- 1
  iso <- pullISO()
  dirs <- pullDir(disease)
  outDir <- dirs[[1]]
  incDir <- dirs[[2]]
  
  ### Test user inputs ### ----
  if (disease %in% diseases == F) { stop(paste0(disease, " is not in the disease list"))}
  if (level %in% levels == F) { stop(paste0(disease, " is not in the level list"))}
  
  ### Create API URL ### ----
  url <- paste0(urlBase, "disease=", disease, "&level=", level)
  
  ### Loop through iso3s pulling each data table ### ----
  for (i in 1:nrow(iso)) {
    
    i <- iso[i]
    message(paste0("Pulling Data for : ", i$ADMIN0ISO3))
    
    # Update URL for current iso2
    currURL <- paste0(url, "&iso2=", i$ADMIN0ISO2)
    
    # Pull data and translate it to a data table
    initPull <- GET(currURL)
    if (initPull$status_code != 200) { stop(paste0("Status Code Does Not Equal 200: ", i$ADMIN0ISO3, " ", currURL))}
    charPull <- rawToChar(initPull$content)
    isoData <- as.data.table(fromJSON(charPull))
    
    # Make sure something is in the file
    if (nrow(isoData) > 0 ) {
      
      ### Save country data to incoming folder ### ----
      if (file.exists(paste0(incDir, "/", i$ADMIN0ISO3)) == F) {
        dir.create(paste0(incDir, "/", i$ADMIN0ISO3))
      }
      fwrite(isoData, paste0(incDir, "/", i$ADMIN0ISO3, "/DATA_", i$ADMIN0ISO3, "_", disease, "_", level, ".csv"))
      
      ### Processing the end ### ----
      if (num == 1) {
        data <- isoData
        num <- num + 1
      } else {
        data <- rbindlist(list(data, isoData), fill = T)
      }
    } else {
      message(paste0("...No ", level, " data for ", i$ADMIN0ISO3))
    }
  }
  
  
  message("SAVING ALL")
  ### Save Africa Data Incoming ### ----
  
  if (file.exists(paste0(incDir, "/ALL")) == F) {
    dir.create(<<<< FILEPATH REDACTED >>>>)
  }
  fwrite(data, <<<< FILEPATH REDACTED >>>>)
  
  ### Any extra processing or formatting? ### ----
  
  ### Need add data check for differences
  
  ### Saving Full Data ### ----
  if (outDir != "None") {
    fwrite(data, paste0(outDir, "/ESPEN_Data_", disease, "_", level, "_", Sys.Date(), ".csv"))
    message(paste0("Data Saved Here : ", outDir, "/ESPEN_Data_", disease, "_", level, "_", Sys.Date(), ".csv"))
  }
  # options(stringsAsFactors = T)
}


### ESPEN API MAPS CALL ###########################################################################
#
### FUNCTION ###
# Interacts with ESPEN API to pull specific geojson data for the user. Usually running only
# disease and type = endemicity. Have not clarified a use for other types
# Currently saves data to two locations with the current date :
#   <<<< FILEPATH REDACTED >>>>
#   <<<< FILEPATH REDACTED >>>>
#
# ONLY WORKS CURRENTLY TO SAVE ENDEMICITY TYPE
#
### PARAMETERS ###
#'@disease - which disease to pull data on { lf | oncho | loa | sch | sth | trachoma }
#'@level - which geographic size level to pull data on. You'll usually want iu { iu | sitelevel }
#'@type - check the ESPEN API website for which type can be pulled with specific disease
#'@year - which year of classifications would you like to pull for - haven't checked what years
#' are available
#'@subtype - check the ESPEN API website for which subtype can be pulled with specific disease
#' and type combination

espenAPIMap <- function(disease, level, type, year = 2017, subtype = "nosubtype") {
  
  ### Base inputs ### ----
  options(stringsAsFactors = F)
  diseases <- c("lf", "oncho", "loa", "sch", "sth", "trachoma")
  levels <- c("sitelevel", "iu")
  urlBase <- "https://admin.espen.afro.who.int/api/maps?"
  
  ### LF Possible Inputs
  lfSLtype <- c("tas", "sentinel_sites", "mapping_surveys")
  lfIUtype <- c("endemicity", "mda_pc_rounds")
  lfIUMDAsubtype <- c("projections", "therapeutic", "geographic")
  ### Oncho Possible Inputs
  onSLtype <- c("impact_assessment", "mapping_surveys")
  onSLsubtype <- c("skin_biopsy", "anti_ov16_test", "nodule_palpation")
  onIUtype <- c("endemicity", "mda_pc_coverage", "mda_pc_rounds")
  onIUMDAsubtype <- c("geographic", "therapeutic")
  ### Loa Possible Inputs
  loaSLtype <- c("mapping_surveys")
  loaSLsubtype <- c("ewh_questionnaire", "blood_smear")
  loaIUtype <- c("endemicity")
  ### Schisto Possible Inputs
  schSLtype <- c("mapping_surveys")
  schSLsubtype <- c("all_species", "s_haematobium", "s_mansoni")
  schIUtype <- c("endemicity",  "mda_pc_coverage", "mda_pc_rounds")
  schIUMDAsubtype <- c("geographic_sac", "geographic_total", "therapeutic_sac", "therapeutic_total")
  ### STH Possible Inputs
  sthSLtype <- c("mapping_surveys")
  sthSLsubtype <- c("all_species", "ascaris", "hookworms", "trichuris")
  sthIUtype <- c("endemicity",  "mda_pc_coverage", "mda_pc_rounds")
  sthIUMDAsubtype <- c("geographic_sac", "geographic_total", "therapeutic_sac", "therapeutic_total")
  ### Trachoma Possible Inputs
  # currently nothing
  
  ### Test user inputs ### ----
  if (disease %in% diseases == F) { stop(paste0(disease, " is not in the disease list"))}
  if (level %in% levels == F) { stop(paste0(disease, " is not in the level list"))}
  ### The rest of these checkes aren't super necessary as this function currently is only
  ### really used for the disease + endemicity pulls
  ### Oncho Input Checks
  if (disease == "oncho" & level == "sitelevel") {
    if (type %in% onSLtype == F) { stop(paste0("The type ", type , " is not in the Oncho Site Level Type list: ", paste(onSLtype, collapse= ", ")))}
    if (subtype %in% onSLsubtype == F) { stop(paste0("The subtype", subtype, " is not in the Oncho Site Level Subtype list: ", paste(onSLsubtype, collapse= ", ")))}
  }
  if (disease == "oncho" & level == "iu") {
    if (type %in% onIUtype == F) { stop(paste0("The type ", type, " is not in the Oncho IU Type list: ", paste(onIUtype, collapse= ", ")))}
    if (grepl("mda", type)) {
      if (subtype %in% onIUMDAsubtype == F) { stop(paste0("The subtype ", subtype, " is not in the Oncho IU MDA/PC Coverage and Rounds Subtype list: ", paste(onIUMDAsubtype, collapse= ", ")))}
    }
  }
  ### LF Input Checks
  if (disease == "lf" & level == "sitelevel") {
    if (type %in% lfSLtype == F) { stop(paste0("The type ", type , " is not in the LF Site Level Type list: ", paste(lfSLtype, collapse= ", ")))}
  }
  if (disease == "lf" & level == "iu") {
    if (type %in% lfIUtype == F) { stop(paste0("The type ", type, " is not in the LF IU Type list: ", paste(lfIUtype, collapse= ", ")))}
    if (grepl("mda", type)) {
      if (subtype %in% lfIUMDAsubtype == F) { stop(paste0("The subtype ", subtype, " is not in the LF IU MDA/PC Coverage and Rounds Subtype list: ", paste(lfIUMDAsubtype, collapse= ", ")))}
    }
  }
  ### Loa Input Checks
  if (disease == "loa" & level == "sitelevel") {
    if (type %in% loaSLtype == F) { stop(paste0("The type ", type , " is not in the Loa Site Level Type list: ", paste(loaSLtype, collapse= ", ")))}
    if (subtype %in% loaSLsubtype == F) { stop(paste0("The subtype ", subtype , " is not in the Loa Site Level Subtype list: ", paste(loaSLsubtype, collapse= ", ")))}
  }
  if (disease == "loa" & level == "iu") {
    if (type %in% loaIUtype == F) { stop(paste0("The type ", type, " is not in the Loa IU Type list: ", paste(loaIUtype, collapse= ", ")))}
  }
  ### Schisto Input Checks
  if (disease == "sch" & level == "sitelevel") {
    if (type %in% schSLtype == F) { stop(paste0("The type ", type , " is not in the Schisto Site Level Type list: ", paste(schSLtype, collapse= ", ")))}
    if (subtype %in% schSLsubtype == F) { stop(paste0("The subtype ", subtype , " is not in the Schisto Site Level Subtype list: ", paste(schSLsubtype, collapse= ", ")))}
  }
  if (disease == "sch" & level == "iu") {
    if (type %in% schIUtype == F) { stop(paste0("The type ", type, " is not in the Schisto IU Type list: ", paste(schIUtype, collapse= ", ")))}
    if (grepl("mda", type)) {
      if (subtype %in% schIUMDAsubtype == F) { stop(paste0("The subtype ", subtype, " is not in the Schisto IU MDA/PC Coverage and Rounds Subtype list: ", paste(schIUMDAsubtype, collapse= ", ")))}
    }
  }
  ### STH Input Checks
  if (disease == "sth" & level == "sitelevel") {
    if (type %in% sthSLtype == F) { stop(paste0("The type ", type , " is not in the STH Site Level Type list: ", paste(schSLtype, collapse= ", ")))}
    if (subtype %in% sthSLsubtype == F) { stop(paste0("The subtype ", subtype , " is not in the STH Site Level Subtype list: ", paste(schSLsubtype, collapse= ", ")))}
  }
  if (disease == "sth" & level == "iu") {
    if (type %in% sthIUtype == F) { stop(paste0("The type ", type, " is not in the STH IU Type list: ", paste(schIUtype, collapse= ", ")))}
    if (grepl("mda", type)) {
      if (subtype %in% sthIUMDAsubtype == F) { stop(paste0("The subtype ", subtype, " is not in the STH IU MDA/PC Coverage and Rounds Subtype list: ", paste(schIUMDAsubtype, collapse= ", ")))}
    }
  }
  
  ### Other small things needed to run ### ----
  num <- 1 # easy tracker
  iso <- pullISO()
  
  ### Create Directories for Saving ### ----
  dirs <- pullDir(disease)
  outDir <- dirs[[1]]
  incDir <- dirs[[2]]
  
  ### Create Values Vectors ### ----
  
  # {"categories":[{"name":"Non-endemic","value":0},{"name":"Endemic (MDA not started)","value":1}
  # ,{"name":"Unknown","value":2}],
  if (disease == "oncho") {
    classF <- c("Non-endemic","Endemic (MDA not started)","Endemic (under MDA)","Unknown (under LF MDA)"
                ,"Endemic (post-MDA surveillance)","Unknown / consider Oncho Elimination mapping")
  } else if (disease == "lf") {
    classF <- c("Non-endemic", "Endemic (MDA not started)", "Endemic (post-MDA surveillance)", "Unknown")
  } else if (disease == "sth") {
    classF <- c("< 1%", "1 - 19.9%", "20 - 49.9%", ">= 50%")
  } else if (disease == "sch") {
    classF <- c("< 1%", "1 - 9.9%", "10 - 49.9%", ">= 50%")
  } else if (disease == "loa") {
    classF <- c("Non-endemic","Endemic (MDA not started)","Unknown")
  }
  
  ### Create API URL ### ----
  url <- paste0(urlBase, "disease=", disease, "&level=", level, "&type=", type)
  if (subtype != "nosubtype") {
    url <- paste0(url, "&subtype=", subtype)
  }
  
  ### Loop through iso3s pulling each one and combining all together ### ----
  for (n in 1:nrow(iso)) {
    i <- iso[n]
    message(paste0("Pulling Data for : ", i$ADMIN0ISO3, " : ", n, "/", nrow(iso)))
    
    # Update URL for current iso2
    currURL <- paste0(url, "&iso2=", i$ADMIN0ISO2, "&start_year=", year, "&end_year=", year)
    
    ### Pull data and translate it to a data table ### ----
    initPull <- GET(currURL)
    if (initPull$status_code != 200) { stop(paste0("Status Code Does Not Equal 200: ", i$ADMIN0ISO3, " - ", initPull$status_code, " - ", currURL))}
    charPull <- rawToChar(initPull$content)
    # Separate the character string by where the geojson actually starts
    starting <- regexpr("\"geojson\":", charPull)
    startN <- starting[1] + 10
    gjsonPull <- to_geojson(substr(charPull, startN, nchar(charPull)-1))
    sfPull <- unique(geojson_sf(gjsonPull))
    
    ### Check if has any data - if no value column then it has no info
    if (nrow(sfPull) > 0 & "value" %in% colnames(sfPull)) {
      
      pullVal <- data.table(value = c(sfPull$value)+1, Endem_MDA = "")
      pullVal[is.na(value), value := 1]
      
      for (v in 1:length(classF)) {
        pullVal[value == v, Endem_MDA := classF[v]]
      }
      
      ### Remove columns st_write doesn't like
      sfPull$Shape_Area <- NULL
      sfPull$Shape_Leng <- NULL
      sfPull$Shape_Le_1 <- NULL
      sfPull$orig_val <- sfPull$value
      sfPull$value <- NULL
      sfPull$Endem_MDA <- pullVal$Endem_MDA
      
      ### Save Country Level Map Info ### ----
      
      ### Create directory if not there yet
      if (file.exists(<<<< FILEPATH REDACTED >>>>) == F) { dir.create(<<<< FILEPATH REDACTED >>>>) }
      
      ### Saving shapefile and combining
      if (type == "endemicity") {
        st_write(sfPull, <<<< FILEPATH REDACTED >>>>)
        
        if (num == 1) {
          data <- sfPull
          num <- num + 1
        } else {
          data <- do.call(rbind, list(data, sfPull))
        }
        if (nrow(data) == 0) {
          stop(message("Data Is Not Combining"))
        }
      } else { ### Should not enter here
        stop(message("...Entering a non Endemcity call - this section is not made"))
      }
    }
  }
  
  ### Save Africa Data Incoming ### ----
  message("SAVING ALL")
  
  ### Create directory
  if (file.exists(<<<< FILEPATH REDACTED >>>>) == F) {
    dir.create(<<<< FILEPATH REDACTED >>>>)
  }
  
  ### Saving to <<<< FILEPATH REDACTED >>>> the full shapefile
  if (type == "endemicity") {
    st_write(data, <<<< FILEPATH REDACTED >>>>)
  }
  
  ### Saving to <<<< FILEPATH REDACTED >>>> the full shapefile but only for Schisto, Oncho, and LF
  if (outDir != "None") {
    if (type == "endemicity") {
      st_write(data, <<<< FILEPATH REDACTED >>>>)
      message(paste0("Data Saved Here : ", o<<<< FILEPATH REDACTED >>>>)
    }
  }
}

###################################################################################################
### HELPER FUNCTIONS ##############################################################################

### Pull ISO ### ----
# Pulls all ISO that ESPEN possibly tracks
pullISO <- function() {
  if (Sys.info()[1] == "Windows") rootJ <- <<<< FILEPATH REDACTED >>>> else rootJ <- <<<< FILEPATH REDACTED >>>>
  allData <- as.data.table(fread(<<<< FILEPATH REDACTED >>>>))
  allData[ADMIN0ISO3 == "NAM", ADMIN0ISO2 := "NA"]
  return(allData)
}

### Pull Directories ###---
# Pulls/creates save directories for a specific disease
pullDir <- function(disease, date = "today") {
  
  if (Sys.info()[1] == "Windows") root <- <<<< FILEPATH REDACTED >>>> else root <- <<<< FILEPATH REDACTED >>>>
  if (date == "today") { date <- Sys.Date() }
  
  ### Incoming Data Directory
  incDir <- <<<< FILEPATH REDACTED >>>>
  
  ### Output of full combined data
  outDir <- <<<< FILEPATH REDACTED >>>>
  
  if (disease == "oncho") {
    outDir <- <<<< FILEPATH REDACTED >>>>
    incDir <- <<<< FILEPATH REDACTED >>>>
  } else if (disease == "lf") {
    outDir <- <<<< FILEPATH REDACTED >>>>
    incDir <- <<<< FILEPATH REDACTED >>>>
  } else if (disease == "sch") {
    outDir <- <<<< FILEPATH REDACTED >>>>
    incDir <- <<<< FILEPATH REDACTED >>>>
  } else if (disease == "sth") {
    outDir <- "None"
    incDir <- <<<< FILEPATH REDACTED >>>>
  } else if (disease == "loa") {
    outDir <- "None"
    incDir <- <<<< FILEPATH REDACTED >>>>
  } else if (disease == "trachoma") {
    outDir <- "None"
    incDir <- <<<< FILEPATH REDACTED >>>>
  }
  
  if(file.exists(incDir) == F) {
    # make file
    dir.create(incDir)
  }
  
  return(list(outDir, incDir))
}

###################################################################################################
### COMPARE FUNCTION ##############################################################################

### Function ###
### Unable to really test currently

comparePulls <- function(disease, level, new_date, old_date
                         , type = "notype", subtype = "nosubtype") {
  
  ### Create Directories ### ----
  dirs <- pullDir(disease, new_date)
  incDirN <- dirs[[2]]
  dirs <- pullDir(disease, old_date)
  incDirO <- dirs[[2]]
  
  ### Pull Data ### ----
  
  if (type == "endemicity") {
    # Pulling shapefile
    stop(message("Not built to compare shapefile currently"))
  } else if (level == "iu") {
    oldData <- fread(<<<< FILEPATH REDACTED >>>>)
    newData <- fread(<<<< FILEPATH REDACTED >>>>)
    
    if ((sort(colnames(oldData)) == sort(colnames(newData))) == F) {
      stop(paste0("Column names do not match up, need do manual comparison check"))
    }
    
    if (nrows(oldData) > nrows(newData)) {
      message(paste0("There are ", nrows(oldData)-nrows(newData), " less rows in the new data set"))
    }
    if (nrows(oldData) < nrows(newData)) {
      message(paste0("There are ", nrows(newData)-nrows(oldData), " more rows in the new data set"))
    }
    if (nrows(oldData) == nrows(newData)) {
      message("Same number of rows in both datasets")
    }
    
    ### Check for Differences ### ----
    
    comData <- unique(rbindlist(list(newData, oldData)))
    
    diffID <- comData[, .(nRows = .N), by = list(IU_ID)]
    diffID <- diffID[nRows > 1]
    
    diffID <- unique(diffID$IU_ID)
    message(paste0("There are ", length(diffID), " IUs that have changed"))
    
    ### Pull ID for new rows or lost rows ### ----
    
    newData$changed <- "new"
    oldData$changed <- "old"
    
    comData <- rbindlist(list(newData, oldData))
    
    updID <- comData[, .(nRows = .N), by = list(IU_ID)]
    updID <- updID[nRows == 1]$IU_ID
    
    ### Pull all rows that are new/old/changed ### ----
    
    diffData <- comData[IU_ID %in% updID | IU_ID %in% diffID, ]
    
  } else if (level == "sitelevel") {
    
    oldData <- fread(<<<< FILEPATH REDACTED >>>>)
    newData <- fread(<<<< FILEPATH REDACTED >>>>)
    
    if ((sort(colnames(oldData)) == sort(colnames(newData))) == F) {
      stop(paste0("Column names do not match up, need do manual comparison check"))
    }
    
    if (nrows(oldData) > nrows(newData)) {
      message(paste0("There are ", nrows(oldData)-nrows(newData), " less rows in the new data set"))
    }
    if (nrows(oldData) < nrows(newData)) {
      message(paste0("There are ", nrows(newData)-nrows(oldData), " more rows in the new data set"))
    }
    if (nrows(oldData) == nrows(newData)) {
      message("...Same number of rows in both datasets")
    }
    
    ### Check for Differences ### ----
    
    comData <- unique(rbindlist(list(newData, oldData)))
    
    diffLL <- comData[, .(nRows = .N), by = list(Latitude, Longitude, Diagnostic)]
    diffLL <- diffLL[nRows > 1]
    
    diffLL <- unique(diffLL[nRows == 1, .(Latitude, Longitude, Diagnostic)])
    message(paste0("There are ", nrow(diffLL), " Sites that have changed"))
    
    ### Pull ID for new rows or lost rows ### ----
    
    newData$changed <- "new"
    oldData$changed <- "old"
    
    comData <- rbindlist(list(newData, oldData))
    
    updLL <- comData[, .(nRows = .N), by = list(Latitude, Longitude, Diagnostic)]
    updLL <- unique(updLL[nRows == 1, .(Latitude, Longitude, Diagnostic)])
    
    ### Pull all rows that are new/old/changed ### ----
    
    pullLL <- rbindlist(list(updLL, diffLL))
    
    for (n in 1:nrow(pullLL)) {
      currR <- pullLL[n]
      pullR <- comData[Latitude == currR$Latitude & Longitude == currR$Longitude
                       & Diagnostic == currR$Diagnostic, ]
      if (n == 1) {
        diffData <- pullR
      } else {
        diffData <- rbindlist(list(diffData, pullR))
      }
    }
  }
  
  
  # Actually want to do all of them as csv because easier
  
  if (type == "endemicity") {
    # Save shapefile
  } else {
    # Saving csv
    fwrite(diffData, <<<< FILEPATH REDACTED >>>>)
  }
}
