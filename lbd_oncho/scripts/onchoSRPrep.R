#
# SCRIPT FOR ONCHO SR DATA
#
# Purpose: All Oncho Systematic Review Data combined/cleaned/hard code fixes of data
###################################################################################################

### Necessary Packages ## ----
library(data.table)
library(readxl)
library(dplyr)


### Example Call ### ----
# You should just run this whole thing like a script 

###################################################################################################
### FUNCTIONS #####################################################################################

### PULL DATA FUNCTION ### ------------------------------------------------------------------------
# Pulls data and returns the data table of it
# Also removes unncessary rows as based on needs
# Only enters this function if needs to be repulled from extraction sheets

pullData <- function(root) {
  # Establishing file paths
  folderRoot <- <<<< FILEPATH REDACTED >>>>
  part1 <- <<<< FILEPATH REDACTED >>>>
  part2 <- <<<< FILEPATH REDACTED >>>>
  part3 <- <<<< FILEPATH REDACTED >>>>
  part4 <- <<<< FILEPATH REDACTED >>>>
  part5 <- <<<< FILEPATH REDACTED >>>>
  hilleleTier1 <- <<<< FILEPATH REDACTED >>>>
  hilleleTier2 <- <<<< FILEPATH REDACTED >>>>
  hilleleTier3 <- <<<< FILEPATH REDACTED >>>>
  hilleleTier4 <- <<<< FILEPATH REDACTED >>>>
  jmbensonTier1 <- <<<< FILEPATH REDACTED >>>>
  jmbensonTier2 <- <<<< FILEPATH REDACTED >>>>
  jmbensonTier3 <- <<<< FILEPATH REDACTED >>>>
  jmbensonTier4 <- <<<< FILEPATH REDACTED >>>>
  jbhallTier1 <- <<<< FILEPATH REDACTED >>>>
  jbhallTier3 <- <<<< FILEPATH REDACTED >>>>
  jbhallTier4 <- <<<< FILEPATH REDACTED >>>>
  cooperProject <- <<<< FILEPATH REDACTED >>>>
  jmbensonNonEnglish <- <<<< FILEPATH REDACTED >>>>
  upd_2019hillele <- <<<< FILEPATH REDACTED >>>>
  upd_2019donkers <- <<<< FILEPATH REDACTED >>>>
  upd_2019donkers_2 <- <<<< FILEPATH REDACTED >>>>

  # Splitting into groups to help trouble shooting / somehow increases runtime
  partsFiles <- c(part1, part2, part3, part4, part5)
  hillFiles <- c(hilleleTier1, hilleleTier2, hilleleTier3, hilleleTier4)
  jmFiles <- c(jmbensonTier1, jmbensonTier2, jmbensonTier3, jmbensonTier4, jmbensonNonEnglish)
  jbFiles <- c(jbhallTier1, jbhallTier3, jbhallTier4, cooperProject)
  upd_2019Files <- c(upd_2019hillele, upd_2019donkers, upd_2019donkers_2)
  
  # READ OUR DATA -----------------------------------------------------------------------------------
  
  # FUNCTION: Read data using fread
  readData <- function(fn){
    dt_temp <- as.data.table(read_excel(fn))
    # Add ext_sheet row to mark where each row came from
    currFile <- unlist(strsplit(fn, "/"))
    currFile <- currFile[length(currFile)]
    dt_temp[!is.na(nid), ext_sheet := currFile]
    message(paste0("done: ", fn))
    return(dt_temp)
  }
  
  # Utilize above function to combine all parts
  parts <- lapply(partsFiles, readData)
  hillele <- lapply(hillFiles, readData)
  jmbenson <- lapply(jmFiles, readData)
  jbhall <- lapply(jbFiles, readData)
  upd_2019 <- lapply(upd_2019Files, readData)
  mylist <- c(parts, hillele, jmbenson, jbhall, upd_2019) # list all together
  mydata <- rbindlist(mylist, fill = TRUE) # Mesh all together
  
  # Light cleaning - Remove no NIDs and rows that start w descriptions
  mydata <- mydata[!is.na(nid), ]
  mydata <- mydata[!grepl("F", nid), ]
  
  message("...pullData() Complete")
  
  return(mydata)
}


### CLEANING DATA ### -----------------------------------------------------------------------------
# Takes in data frame and cleans it up, Tracking issue rows to later be removed
# Marks rows that have been changed/added for future comparison needs

cleanData <- function(fullData) {
  
  ### BASIC CLEANING ### Data types, legitimate row removals --------------------------------------
  
  # Remove rows that aren't actual extraction rows and add column to mark issue, if any
  fullData <- fullData[bundle_name == "ntd_oncho" , issue := ""]
  # Change types of the rows to match what it should be
  fullData$lat <- as.numeric(as.character(fullData$lat))
  fullData$long <- as.numeric(as.character(fullData$long))
  fullData$approx_lat <- as.numeric(as.character(fullData$approx_lat))
  fullData$approx_long <- as.numeric(as.character(fullData$approx_long))
  fullData$age_start <- as.numeric(as.character(fullData$age_start))
  fullData$age_end <- as.numeric(as.character(fullData$age_end))
  fullData$shape_type <- as.character(fullData$shape_type)
  
  # Some adjustments
  fullData$poly_id <- as.numeric(as.character(fullData$poly_id))
  fullData[is.na(poly_id), poly_id := -1]
  
  fullData$poly_reference <- as.character(fullData$poly_reference)
  fullData[is.na(poly_reference) | poly_reference == "", poly_reference := "NA"]
  
  # Fix some month columns
  # Month_start
  fullData[month_start == "June", month_start := 6]
  fullData[month_start == "December", month_start := 12]
  fullData[month_start == "May", month_start := 5]
  fullData$month_start <- as.integer(as.character(fullData$month_start))
  fullData[is.na(month_start), month_start := 0]
  
  # Month_end
  fullData[month_end == "July", month_end := 7]
  fullData$month_end <- as.integer(as.character(fullData$month_end))
  fullData[is.na(month_end), month_end := 0]
  
  # Change NAs because it doesn't like me comparing NAs later
  fullData[is.na(year_start) | year_start == "" | year_start == "NA", issue := paste0(issue, "No year start; ")]
  fullData[is.na(year_end) | year_end == "" | year_end == "NA", issue := paste0(issue, "No year end; ")]
  # Year to numberic
  fullData$year_start <- as.numeric(as.character(fullData$year_start))
  fullData$year_end <- as.numeric(as.character(fullData$year_end))
  
  ### GEOREFERENCING CLEANING ### Standardizing shape_type, poly_type, etc. -----------------------
  
  # Change "Polygon Buffers" - should be marked as points, none of these have a poly_reference
  fullData[shape_type == "polygon buffer", poly_type := "buffer"]
  fullData[shape_type == "polygon buffer", shape_type := "point"]
  
  # Change rows of shape_type point, poly_type to buffer or point
  fullData[shape_type == "point" & buffer_radius != "" & !is.na(buffer_radius), poly_type := "buffer"]
  fullData[shape_type == "point" & (poly_type != "buffer" | is.na(poly_type)), poly_type := "point"]
  
  # Check poly_reference == buffer - these are points with buffers
  fullData[poly_reference == "buffer", shape_type := "point"]
  fullData[poly_reference == "buffer", poly_type := "buffer"]
  fullData[poly_reference == "buffer", poly_reference := NA]
  
  # Check poly_type == buffer - know they are all points now
  fullData[poly_type == "buffer" & shape_type != "point", shape_type := "point"]
  
  # Replace "" with NAs
  fullData[shape_type == "" | shape_type == "NA", shape_type := NA]
  fullData[poly_type == "" | poly_type == "NA", poly_type := NA]
  fullData[poly_reference == "" | poly_reference == "NA", poly_reference := NA]
  fullData[lat == "" | lat == "NA", lat := NA]
  fullData[long == "" | long == "NA", long := NA]
  fullData[buffer_radius == "" | buffer_radius == "NA", buffer_radius := NA]
  
  # Push all approx_lat, longs, buffers into exact lat long columns
  fullData[, latlong_change := 0]
  fullData[shape_type == "point" & (is.na(lat) | is.na(long)), latlong_change := 1]
  fullData[latlong_change == 1 & is.na(lat), lat := approx_lat]
  fullData[latlong_change == 1 & is.na(long), long := approx_long]
  fullData[, buffer_change := 0]
  fullData[shape_type == "point" & poly_type == "buffer" & is.na(buffer_radius), buffer_change := 1]
  fullData[buffer_change == 1, buffer_radius := approx_buffer]
  
  # Track rows where shape_type or poly_type have NAs
  fullData[is.na(shape_type) | shape_type == "" | shape_type == "NA", issue := paste0(issue, "Missing shape_type; ") ]
  fullData[is.na(poly_type ) | poly_type == ""  | poly_type == "NA" , issue := paste0(issue, "Missing poly_type; ") ]
  
  # Mark invalid lat/longs
  fullData[shape_type == "point" & (lat < -90 | lat > 90), issue := paste0(issue, "Invalid Latitude; ") ]
  fullData[shape_type == "point" & (long < -180 | long > 180), issue := paste0(issue, "Invalid Longitude; ") ]
  
  # Track rows where shape_type = polygon, but poly_reference is NA
  fullData[shape_type == "polygon" & (is.na(poly_reference) | poly_reference == "" | poly_reference == "NA")
           , issue := paste0(issue, "Marked as polygon in shape_type but has no poly_reference; ") ]
  
  # Normalize poly_reference
  fullData[!is.na(poly_reference), poly_reference := gsub("\\\\", "/", poly_reference)]
  fullData[!is.na(poly_reference) & !endsWith(poly_reference, ".shp"), poly_reference := paste0(poly_reference, ".shp")]
  
  ### AGE CLEANING ### Adding / tracking ages by group for comparison needs -----------------------
  
  # Add row to track when changed
  fullData[ , age_change := 0]
  # Mark which rows need to be changed
  fullData[is.na(age_start) | is.na(age_end) | age_start == "" | age_end == "", age_change := 1]
  fullData[age_change == 1 & age_group_category == "adults", age_start := 20 ]
  fullData[age_change == 1 & age_group_category == "adults", age_end := 99 ]
  fullData[age_change == 1 & grepl("children", age_group_category, fixed = TRUE), age_start := 5]
  fullData[age_change == 1 & grepl("children", age_group_category, fixed = TRUE), age_end := 15]
  fullData[age_change == 1 & age_group_category == "newborns", age_start := 0 ]
  fullData[age_change == 1 & age_group_category == "newborns", age_end := 1 ]
  fullData[age_change == 1 & age_group_category == "post-puberty", age_start := 12 ]
  fullData[age_change == 1 & age_group_category == "post-puberty", age_end := 18 ]
  fullData[age_change == 1 & age_group_category == "SAC", age_start := 4 ]
  fullData[age_change == 1 & age_group_category == "SAC", age_end := 14 ]
  fullData[age_change == 1 & (unit_type == "Fly" | unit_type == "parous fly" | unit_type == "Parous Fly"), age_start := 0]
  fullData[age_change == 1 & (unit_type == "Fly" | unit_type == "parous fly" | unit_type == "Parous Fly"), age_end := 99]
  
  ## Hardcoded checks ##
  # nid 327974 is unpublished REMO data
  fullData[nid == 327974 & unit_type == "Person", age_change := 1]
  fullData[nid == 327974 & age_change == 1, age_start := 20]
  fullData[nid == 327974 & age_change == 1, age_end := 99]
  # nid 286983 is non descript ministry of health data
  fullData[nid == 286983 & unit_type == "Person" & is.na(age_start), age_change := 1]
  fullData[nid == 286983 & age_change == 1, age_start := 0]
  fullData[nid == 286983 & age_change == 1, age_end := 99]
  
  # Track issues with rows that have no age start or end marked
  fullData[is.na(age_start) | age_start == "" | age_start == "NA", issue := paste0(issue, "No age_start; ") ]
  fullData[is.na(age_end) | age_end == "" | age_end == "NA", issue := paste0(issue, "No age_end; " ) ]
  
  ### ORGANIZE BY DX ### --------------------------------------------------------------------------
  
  # DX groupings
  ss <- c(1, 2, 9, 10)
  nod <- c(30, 92, 95, 39)
  sero <- c(100, 102, 416, 426, 427, 429)
  fly <- c(200, 201, 202, 417)
  eye <- c(136, 29, 131, 130, 150, 13, 141, 143, 139, 140, 25, 411, 434, 151, 145, 146, 300, 26, 301, 15, 409, 148, 147, 153, 432, 20, 133, 21, 24, 28, 415, 19, 27, 17, 22, 23, 149, 152, 16, 403, 405, 406, 135, 436, 437, 138, 142, 438, 439, 410, 435, 14, 404, 144, 18, 436, 407, 137, 441)
  skin <- c(423, 40, 419, 420, 52, 431, 47, 55, 60, 63, 34, 35, 64, 65, 46, 400, 430, 53, 424, 402, 44, 49, 61, 62, 401, 418, 41, 66, 42, 36, 433, 70, 71, 68, 75, 74, 94, 58, 413, 67, 440, 442)
  otherPrev <- c(91, 93, 94, 425, 428, 422)
  
  # Create new column for DX type and sort through them by the groupings above
  fullData[, DX := "NA"]
  fullData[cv_DX_type %in% ss, DX := "ss"]
  fullData[cv_DX_type %in% nod, DX := "nod"]
  fullData[cv_DX_type %in% sero, DX := "sero"]
  fullData[cv_DX_type %in% fly, DX := "fly"]
  fullData[cv_DX_type %in% eye, DX := "eye"]
  fullData[cv_DX_type %in% skin, DX := "skin"]
  fullData[cv_DX_type %in% otherPrev, DX := "otherPrev"]
  fullData[DX == "NA", DX := NA]
  
  # Track any rows that didn't get marked with a DX
  fullData[is.na(DX), issue := paste0(issue, "No DX match; ")]
  
  # Track the fly data rows
  fullData[DX == "fly", issue := paste0(issue, "Fly data; ")]
  
  message("...cleanData() Complete")
  
  return(fullData)
}

### REMOVE UNNECESSARY ROWS FUNCTION ### ----------------------------------------------------------
# Marks rows of data based on NID / Site_memo combos that have been manually checked and
# determined to be duplicate data by location in previous checks, tied together through a csv
# Leaving in these marks for the untouched so GBD side knows what is duplicated
#
# PARAMETERS:
#'@fullData : data table, full in tact 
#'@removingFile : a file with nid/site_memo combos previously found to be duplicates of other rows

markLocDup <- function(fullData, removingFile) {
  
  # Give column to mark duplicates
  fullData[!is.na(nid), loc_dup := 0]
  
  # Remove specific NID overarching / underling issues
  removeThese <- as.data.table(read_excel(removingFile))
  
  # Loop through and pull the nid and site_memo combo by row and then mark them
  for (i in 1:nrow(removeThese)) {
    currRow <- removeThese[i]
    fullData[nid == currRow$nid & site_memo == currRow$remove_site_memo
             , loc_dup := 1]
    fullData[nid == currRow$nid & site_memo == currRow$remove_site_memo
             , issue := paste0(issue, "Location Duplicate Removal; ") ]
  }
  
  message("...markLocDup() Complete")
  
  return(fullData)
}

### CALCULATE NUMBERS ### -----------------------------------------------------------------------
# Takes in the data table and makes sure columns mean, cases, and sample_size are calculated
# Really only cares that we have cases and sample_size numbers calculated

calculateNumbers <- function(fullData) {
  
  # Normalize the variables
  fullData$mean <- as.numeric(as.character(fullData$mean))
  fullData$cases <- as.numeric(as.character(fullData$cases))
  fullData$sample_size <- as.numeric(as.character(fullData$sample_size))
  
  # Mark rows that are impossible to calculate
  fullData[is.na(mean) & is.na(cases) & is.na(sample_size), issue := paste0(issue, "Missing mean, cases, and sample_size; ")]
  fullData[sample_size == 0, issue := paste0(issue, "sample_size is 0, impossible; ")]
  
  # Missing sample_size, but not enough info to calculate it
  fullData[issue == "" & is.na(sample_size) & (is.na(mean) | is.na(cases))
           , issue := paste0("No sample_size but not enough info to calculate; ")]
  # Missing cases, but not enough info to calculate it
  fullData[issue == "" & is.na(cases) & (is.na(mean) | is.na(sample_size))
           , issue := paste0("No cases but not enough info to calculate; ")]
  
  # Calculate columns - using ceiling to round up to next whole number. Don't want to underestimate
  fullData[issue == "" & is.na(sample_size)
           , calculated := 1]
  fullData[issue == "" & is.na(sample_size)
           , sample_size := ceiling( cases / mean )]
  fullData[issue == "" & is.na(cases)
           , calculated := 1]
  fullData[issue == "" & is.na(cases)
           , cases := ceiling( mean * sample_size )]
  
  
  # Message and return the data.table
  message("...calculateNumbers() Complete")
  return(fullData)
}

### HARDCODE ISSUES MARKING ### -----------------------------------------------------------------
# Used for some hardcoded issues marking like ESPEN data, outliers, etc

hardcodeIssues <- function(fullData) {
  
  # NID: 327881
  # Issue: Already included in ESPEN data
  fullData[nid == 327881, ESPEN_data := 1]
  fullData[nid == 327881, issue := paste0(issue, "ESPEN Data - already marked; ")]
  
  # NID: 334473
  # Issue: Too large of a space, want to find smaller units
  fullData[nid == 334473, outlier := 1]
  fullData[nid == 334473, issue := paste0(issue, "Double check georeferencing - need to find smaller unit; ")]
  
  # NID: 338594
  # Double extracted, but did not want to remove from the extraction sheet.
  # Redid georeferencing
  fullData[nid == 338594 & extractor == "c8louie", issue := paste0(issue, "Already extracted from before in another extraction; ")]
  
  # NID : 332842
  # No way to georeference
  fullData[nid == 332842 & site_memo == "Kwaka | Gurara | Niger State | Nigeria"
           , issue := paste0(issue, "Unable to georeference; ")]
  
}

### DECIDE GROUP ### ------------------------------------------------------------------------------
# Needs a properly formatted shape_type

decideGroup <- function(data, groupNum) {
  
  ### Points ### ----
  
  data$LBD_group <- as.integer(data$LBD_group)
  data[shape_type == "point", LBD_group := .GRP
       , by = list(lat, long, year_start, year_end, month_start, month_end)]
  data$LBD_group <- as.double(data$LBD_group)
  data[shape_type == "point", LBD_group := LBD_group + groupNum]
  # Update group number for polygons
  groupNum <- max(data$LBD_group) + 1
  
  ### Polygons ### ----
  
  data$LBD_group <- as.integer(data$LBD_group)
  data[shape_type == "polygon", LBD_group := .GRP
       , by = list(poly_type, poly_reference, poly_id, year_start, year_end, month_start, month_end)]
  data$LBD_group <- as.double(data$LBD_group)
  data[shape_type == "polygon", LBD_group := LBD_group + groupNum]
  
  return(data)
}

### DECIDE UNIQUE ### -----------------------------------------------------------------------------

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

### CHECK AGGREGATED OR NOT ### -------------------------------------------------------------------

aggCheck <- function(currSite, unqSites) {
  
  # Set up small table to check all site_memos
  unqSites <- unique(data.table(
    sites = unqSites,
    under = FALSE
  ))
  
  # Pull out first part of the site_memo
  currSite <- strsplit(as.character(currSite), "\\|")[[1]][1]
  
  # Set any row that does not start with the currSite but does contain the currSite to true
  unqSites[!startsWith(sites, currSite) & grepl(currSite, sites, fixed = TRUE)
           , under := TRUE]
  
  # If any are marked as true, then the currSite memo is considered an overarching site_memo
  if (nrow(unqSites[under == TRUE, ]) > 0 ) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

### DECIDE SPECIFICITY ### ------------------------------------------------------------------------

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

### FIX NUMBERING ISSUES ### --------------------------------------------------------------------
# Hardcode group or unique number fixing
#

adjustNumbers <- function(fullData) {
  
  ### Hardcode fixes, checks matches NID / site_memo / year_start / year_end to ensure
  ### looking at same row, without trying to handle changing group numbers
  
  # Issue: Mistakes Garmadi as being an overarching
  # Solution: hardcode to change the unique number
  fullData[nid == 328063 & site_memo == "Garmadi (village)|Kauru (local government area)|Kaduna State (state)|Nigeria (country)"
           & year_start == 1987 & year_end == 1987
           , LBD_uniq := 12]
  
  # Issue: Mistakes Garmadi as being an overarching, diff years from above
  # Solution: hardcode to change the unique number
  fullData[nid == 328063 & site_memo == "Garmadi (village)|Kauru (local government area)|Kaduna State (state)|Nigeria (country)"
           & year_start == 2008 & year_end == 2008
           , LBD_uniq := 12]
  
  # Issue: Mistakes El Recuerdo as being an overarching
  # Solution: hardcode the unique number
  fullData[nid == 287122 & site_memo == "El Recuerdo (community)|Chicacao (municipality)|Suchitepequez (province)|Central Endemic Zone (foci)|Guatemala (country)"
           & year_start == 2010 & year_end == 2010
           , LBD_uniq := 30]
  
  # Issue: Mistaken numbers as overarching
  # Solution: hardcode the unique number
  fullData[nid == 332786 & site_memo == "1 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 25]
  fullData[nid == 332786 & site_memo == "2 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 26]
  fullData[nid == 332786 & site_memo == "3 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 27]
  fullData[nid == 332786 & site_memo == "4 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 28]
  fullData[nid == 332786 & site_memo == "5 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 29]
  fullData[nid == 332786 & site_memo == "6 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 30]
  fullData[nid == 332786 & site_memo == "7 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 31]
  fullData[nid == 332786 & site_memo == "8 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 32]
  fullData[nid == 332786 & site_memo == "9 (village number)|heavily forested (area)|Dja et Lobo (division)|Cameroon (country)"
           & year_start == 1992 & year_end == 1992
           , LBD_uniq := 33]
  
  # Issue: Marks Endemic Villages as ovearching
  # Solution: Change the unique number
  fullData[nid == 334422 & site_memo == "Endemic Villages (Village)|Nebbi (District)|Uganda (Country)"
           & year_start == 1994 & year_end == 1994
           , LBD_uniq := 3]
  
  # Issue:
  # Solution:
  #fullData[nid == 0 & site_memo == ""
  #         & year_start == 0 & year_end == 0
  #         , ]
  
  
  message("...Hardcode Adjustments Finished")
  return(fullData)
}

###################################################################################################
### RUNNING #######################################################################################

# For other script 
if (Sys.info()[1] == "Windows") { root <- <<<< FILEPATH REDACTED >>>> } else { root <- <<<< FILEPATH REDACTED >>>>}
rePull <- F

prepOnchoSR <- function(root, rePull = F) {
  
  ### No more repull 
  rePull <- F
  
  ### Pull and/or combine all oncho extraction data ### ----
  
  untDir <- <<<< FILEPATH REDACTED >>>>
  
  # If user specifies to pull extraction sheets info, otherwise pull previous untouched version
  if (rePull) {
    message("Begin Pulling data from extraction sheets")
    srData <- pullData(root)
    fwrite(srData, <<<< FILEPATH REDACTED >>>>)
    fwrite(srData, <<<< FILEPATH REDACTED >>>>)
    message(paste0("Initial data pull has been saved at : \n...", <<<< FILEPATH REDACTED >>>>))
  } else {
    srData <- as.data.table(read.csv(<<<< FILEPATH REDACTED >>>>))
    srData$site_memo <- as.character(srData$site_memo)
    message("Data pulled from previously compiled extraction sheets")
  }
  
  ### Update MDA Status ### ----
  
  ### Pull Katie Made fixes, slight tidying 
  mdaAdd <- as.data.table(read_excel(<<<< FILEPATH REDACTED >>>>))
  mdaAdd <- mdaAdd[, .(count = .N), by = list(nid, site_memo, year_start, MDA_adjust)]
  message(paste0("NID / Site_memo undefined MDA status remaining : ", nrow(mdaAdd[MDA_adjust > 2, ])))
  mdaAdd[MDA_adjust == 0, new_MDA := "No MDA"]
  mdaAdd[MDA_adjust == 1 | MDA_adjust == 2, new_MDA := "MDA has started"]
  mdaAdd <- mdaAdd[MDA_adjust < 3, ]
  
  ### Loop through MDA updates and fix or add classifications for MDA in the main data set
  for (i in 1:nrow(mdaAdd)) {
    curr <- mdaAdd[i]
    srData[nid == curr$nid & site_memo == curr$site_memo & year_start == curr$year_start
           & (cv_MDA_status == "NA" | is.na(cv_MDA_status)), cv_MDA_status := curr$new_MDA]
  }
  
  ### Clean and tidy data ### -----
  
  ### Issue marking prep
  srData[ , ESPEN_data := 0]
  srData[ , outlier := 0]
  srData[ , calculated := 0]
  
  ### Issue marking
  srData <- cleanData(srData)              # marks issues with data formatting
  srData <- markLocDup(srData              # marks in issue box the location duplicates
                       , <<<< FILEPATH REDACTED >>>>)
  srData <- calculateNumbers(srData)       # marks rows that can't calculate Sample Size, mean, or cases
  srData <- hardcodeIssues(srData)         # does some hardcode fixes / markings
  
  ### Remove marked issue rows
  issues  <- srData[issue != "", ]
  srData <- srData[issue == "", ]
  srData[, source := "sys_rev"]
  
  ### Numbering ### ----
  
  ### Column Prep
  srData[ , LBD_group := 0]
  srData[ , LBD_uniq := 0]
  srData[ , LBD_spec := ""]
  
  ### Number the rows
  print("Begin Numbering")
  srData <- decideGroup(srData, 500)
  srData <- decideUniq(srData)
  srData <- decideSpec(srData)
  srData <- adjustNumbers(srData)
  print("Finish Numbering")
  
  ### Finishing up ### ----
  
  srData[month_start == 0, month_start := NA]
  srData[month_end == 0, month_end := NA]
  
  ### Save cleaned/numbered SR data 
  fwrite(srData, <<<< FILEPATH REDACTED >>>>)
  
  ### Save issues next to cleaned/numbered data ### ---- 
  
  issues[, confirmed := "no"]
  issues[nid == 332212 & issue == "No sample_size but not enough info to calculate; "
         , confirmed := "hillele"]
  issues[grepl("Fly data", issue) | grepl("Location Duplicate", issue)
         | grepl("Already extracted", issue) | grepl("ESPEN", issue), confirmed := "hillele"]
  issues[nid == 286983 & issue == "No sample_size but not enough info to calculate; "
         , confirmed := "hillele"]
  issues[nid == 332842 & issue == "Unable to georeference; "
         , confirmed := "hillele"]
  issues[nid == 327885 & issue == "Missing shape_type; Missing poly_type; "
         , confirmed := "hillele"]
  issues[nid == 327988 & issue == "No sample_size but not enough info to calculate; "
         , confirmed := "hillele"]
  
  message(paste0("Issue rows without being confirmed : ", nrow(issues[confirmed == "no", ])))
  
  fwrite(issues, <<<< FILEPATH REDACTED >>>>)

}

prepOnchoSR(root, rePull)


