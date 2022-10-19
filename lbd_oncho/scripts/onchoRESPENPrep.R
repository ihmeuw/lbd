#
# SCRIPT FOR ONCHO REMO ESPEN DATA PREP
#
# Purpose: All Oncho REMO and ESPEN Data combined/cleaned/hard code fixes of data
###################################################################################################

### Necessary Packages ## ----
library(data.table)
library(readxl)
library(stringr)

### Example Call ### ----
# You should just run this whole thing like a script 

###################################################################################################
### FUNCTIONS #####################################################################################

### PULL ESPEN'S LATEST DATA PULL DATE ### --------------------------------------------------------

espenPullDate <- function(root) {
  
  # Pull all possible files to be used 
  fp <- <<<< FILEPATH REDACTED >>>>
  files <- list.files(fp, "sitelevel")
  
  # Sort just the dates in the files 
  dates <- str_replace(files, "ESPEN_Data_oncho_sitelevel_", "")
  dates <- str_replace(dates, ".csv", "")
  dates <- sort(dates, decreasing = T)
  
  # Return newest date
  return(dates[1])
  
}

### CLEAN REMO DATA ### ---------------------------------------------------------------------------

cleanREMO <- function(allREMO, respenISSUES) {
  
  ### Hardcode iso3 name updates to match ours ### ----
  allREMO[COUNTRY == "D.R. CONGO", COUNTRY := "Democratic Republic of the Congo"]
  allREMO[COUNTRY == "CAMEROUN", COUNTRY := "Cameroon"]
  allREMO[COUNTRY == "CAR", COUNTRY := "Central African Republic"]
  
  ### Pull not usable rows ### ----
  allREMO[is.na(LAT) | is.na(LON), issue := "Missing lat and/or long"]
  respenISSUES <- rbindlist(list(respenISSUES, allREMO[!is.na(issue) ]), fill = TRUE)
  allREMO <- allREMO[is.na(issue), ]
  
  ### Add iso3 codes ### ----
  iso3Codes <- pullISO(root)
  for (curr in unique(allREMO$COUNTRY)) {
    currISO <- unique(iso3Codes[tolower(location_name) == tolower(curr), ]$iso_code)
    allREMO[COUNTRY == curr, iso3 := currISO]
  }
  
  # Fill in TOT_EXAM / TOT_NOD so we don't have empty values
  allREMO[is.na(TOT_NOD), TOT_NOD := Final_NOD]
  allREMO[is.na(TOT_EXAM), TOT_EXAM := Final_EXAM]
  
  allREMO[is.na(YEAR_REMO), YEAR_REMO := 0]
  
  # Create new data table of only relevant columns
  remoData <- data.table(
    row_id = allREMO$row_id,
    nid = 0, 
    country = tolower(allREMO$COUNTRY),
    iso3 = allREMO$iso3,
    source = "REMO",
    survey_year = allREMO$YEAR_REMO,
    survey_month = 0, 
    location_name = tolower(allREMO$VILLAGE),
    site_memo = tolower(allREMO$CountryVillage),
    examined = allREMO$TOT_EXAM,
    positive = allREMO$TOT_NOD,
    latitude = allREMO$LAT,
    longitude = allREMO$LON,
    dx_code = 30,
    dx = "nod",
    sex = "Both", 
    MDA_status = 0
  )
  
  return(list(remoData, respenISSUES))
}

### PULL ISO3 CODES ### ---------------------------------------------------------------------------

pullISO <- function(root) {
  iso3Codes <- fread(<<<< FILEPATH REDACTED >>>>)
  setnames(iso3Codes, old = "country", new = "iso_code")
  iso3Codes$region_name <- NULL
  return(iso3Codes)
}

### CLEAN ESPEN DATA ### --------------------------------------------------------------------------

cleanESPEN <- function(allESPEN, respenISSUES) {
  
  # Update id to include letters for ESPEN (_ESP) 
  allESPEN[, row_id := paste0(ID, "_ESP")]
  
  # Determine MDA status from Assessment Type
  # Mapping = before MDA; otherwise it is after MDA
  allESPEN[, MDA_status := 0]
  allESPEN[Method_0 != "Mapping", MDA_status := 1]
  
  # Determine diagnostic shorts
  allESPEN[Method_2 == "Nodule palpation", dx := "nod"]
  allESPEN[Method_2 %in% c("Immunology (Ov16)", "RDT-Ov16"), dx := "sero"]
  allESPEN[Method_2 == "Parasitology (Skin biopsy)", dx := "ss"]
  allESPEN[Method_2 == "PCR", dx := "Fly"]
  # Adding dx code
  allESPEN[Method_2 == "Nodule palpation", dx_code := 30]
  allESPEN[Method_2 == "Immunology (Ov16)", dx_code := 100]
  allESPEN[Method_2 == "RDT-Ov16", dx_code := 102]
  allESPEN[Method_2 == "Parasitology (Skin biopsy)", dx_code := 1]

  # Format a site_memo to be used in comparison
  allESPEN[, site_memo := paste0(Country, "_", LocationName)]
  
  # Change months to numbers 
  months <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
  allESPEN[tolower(SurveyMonth) %in% tolower(months), survey_month := match(tolower(SurveyMonth), tolower(months))]
  allESPEN[SurveyMonth == "Novembre", survey_month := 11]
  allESPEN[is.na(survey_month), survey_month := 0]

  allESPEN[is.na(SurveyYear), SurveyYear := 0]
  
  # Mark and remove issues
  # Mark issues 
  allESPEN[ dx == "Fly", issue := "Fly data"]
  allESPEN[is.na(Latitude) | is.na(Longitude), issue := "Missing lat and/or long"]
  respenISSUES <- rbindlist(list(respenISSUES, allESPEN[!is.na(issue) ]), fill = TRUE)
  allESPEN <- allESPEN[is.na(issue), ]
  
  # Create new table properly formatted for merging
  espenData <- data.table(
    row_id = allESPEN$row_id,
    nid = 0, 
    country = tolower(allESPEN$Country),
    iso3 = allESPEN$ISO3,
    source = "ESPEN",
    survey_year = allESPEN$SurveyYear,
    survey_month = allESPEN$survey_month,
    location_name = tolower(allESPEN$LocationName),
    site_memo = tolower(allESPEN$site_memo),
    examined = allESPEN$Examined,
    positive = allESPEN$Positive,
    latitude = allESPEN$Latitude,
    longitude = allESPEN$Longitude,
    age_start = allESPEN$Age_start,
    age_end = allESPEN$Age_end,
    dx_code = allESPEN$dx_code, 
    dx = allESPEN$dx,
    sex = "Both", 
    MDA_status = allESPEN$MDA_status
  )
  
  return(list(espenData, respenISSUES))
}

### COMPARE MORE THAN TWO ROWS ### ----------------------------------------------------------------

compMore <- function(moreR, allRESPEN) {
  
  # Remove easy one rows (ESPEN > REMO) 
  cleanM <- allRESPEN[LBD_group %in% moreR & source == "ESPEN", ]
  grpsW_ESPEN <- unique(cleanM$LBD_group) # for use later 
  cleanM <- cleanM[LBD_group %in% cleanM[, .(count = .N), by = LBD_group][count == 1]$LBD_group, ]
  cleanM[, LBD_uniq := 1]
  cleanM[, LBD_spec := ""]

  mRows <- allRESPEN[LBD_group %in% unique(moreR) & !(LBD_group %in% unique(cleanM$LBD_group)), ]
  mGrps <- unique(mRows$LBD_group)

  for (cGrp in mGrps) {
    cRows <- mRows[LBD_group == cGrp & source == "ESPEN", ]

    cLocs <- unique(cRows$location_name)

    for (cLoc in cLocs) {
      
      cLRows <- cRows[location_name == cLoc, ]
      if (nrow(cLRows) == 1 & length(cLocs) == 1) {
        cleanM <- rbindlist(list(cleanM, cLRows[, LBD_uniq := 1]), fill = T)
      } else if (nrow(cLRows) == 1 & length(cLocs) > 1) {
        cleanM <- rbindlist(list(cleanM, cLRows[, LBD_uniq := match(cLoc, cLocs) + 1]), fill = T)
      } else {
        
        cLRows[, LBD_uniq := match(cLoc, cLocs) + 1]
        cDXs <- unique(cLRows$dx)

        for (cDX in cDXs) {
          
          cDRows <- cLRows[dx == cDX, ]
          if (nrow(cDRows) == 1) {
            cleanM <- rbindlist(list(cleanM, cDRows[, LBD_spec := "DX "]), fill = T)
          } else {
            stop(paste0(cGrp, " ", cLoc, " ", cDX))
          }
        }
      }
    }
  }
  
  return(cleanM)
}


### COMPARE TWO ROWS ### --------------------------------------------------------------------------

compTwo <- function(twoR, allRESPEN) {
  
  # Remove easy one rows (ESPEN > REMO) 
  cleanTwo <- allRESPEN[LBD_group %in% twoR & source == "ESPEN", ]
  cleanTwo <- cleanTwo[LBD_group %in% cleanTwo[, .(count = .N), by = LBD_group][count == 1]$LBD_group, ]
  cleanTwo[, LBD_uniq := 1]
  cleanTwo[, LBD_spec := ""]

  twoRows <- allRESPEN[LBD_group %in% unique(twoR) & !(LBD_group %in% unique(cleanTwo$LBD_group)), ]
  
  ### Finding groups with two location_names, marking them as unique
  twoLocs <- twoRows[, .(locations = length(unique(location_name))), by = list(LBD_group)][locations == 2]$LBD_group
  
  # giving unique numberings
  for (cGrp in twoLocs) {
    currLocations <- unique(twoRows[LBD_group == cGrp]$location_name)
    twoRows[LBD_group == cGrp & location_name == currLocations[1], LBD_uniq := 2]
    twoRows[LBD_group == cGrp & location_name == currLocations[2], LBD_uniq := 3]
  }
  
  ### Finding groups with two location_names, marking them as unique
  twoDX <- twoRows[, .(dxs = length(unique(dx))), by = list(LBD_group)][dxs == 2]$LBD_group
  twoRows[, LBD_spec := ""]
  twoRows[LBD_group %in% twoDX, LBD_spec := "DX "]
  
  if (nrow(twoRows[is.na(LBD_uniq) & is.na(LBD_spec), ]) > 0) { stop("Problem in two rows section")}
  
  twoRows[is.na(LBD_uniq), LBD_uniq := 1]
  
  return(rbindlist(list(cleanTwo, twoRows), fill = T)) 

}

###################################################################################################
### RUNNING #######################################################################################

# For other script 
root <- "<REDACTED>"

prepOnchoRESPEN <- function(root) {
  
  ### Pull all data ### ----
  
  ### REMO
  remoF <- <<<< FILEPATH REDACTED >>>>
  allREMO <- as.data.table(read_excel(remoF))
  
  ### ESPEN
  espenDate <- espenPullDate(root)
  espenF <- <<<< FILEPATH REDACTED >>>>
  allESPEN <- as.data.table(fread(espenF))
  
  ### Prep Data for Merge ### ----

  # REMO
  remoData <- cleanREMO(allREMO, data.table())
  respenISSUES <- remoData[[2]]
  remoData <- remoData[[1]]
  # ESPEN 
  espenData <- cleanESPEN(allESPEN, respenISSUES)
  respenISSUES <- espenData[[2]]
  espenData <- espenData[[1]]
  
  ### Merge and light cleaning ### ----
  
  allRESPEN <- rbindlist(list(espenData, remoData), fill = TRUE)
  
  # Standardize Congo
  congoNames <- c("democratic republic of the congo", "congo (kinshasa)")
  allRESPEN[country %in% congoNames, country := "d.r. congo"]
  allRESPEN[country == "d.r. congo", site_memo := paste0(country, "_", location_name)]
  
  
  ### Adding group / counting duplicate rows ### ----
  
  # Adding Group 
  allRESPEN[, LBD_group := .GRP, by = list(latitude, longitude, survey_year)]
  allRESPEN$LBD_group <- as.double(allRESPEN$LBD_group)
  
  # Pulling counts of rows in a single group and separating
  groupCounts <- allRESPEN[, .(count = .N), by = list(LBD_group)]
  
  ### Separate and sort by possible duplicate rows ### ----
  
  singleR <- unique(groupCounts[count == 1]$LBD_group)
  twoR <- unique(groupCounts[count == 2]$LBD_group)
  moreR <- unique(groupCounts[count > 2]$LBD_group)
  
  # Single rows == no duplicates 
  cleanData <- allRESPEN[LBD_group %in% singleR, ]
  cleanData[, LBD_uniq := 1]
  cleanData[, LBD_spec := ""]
  
  # Two Rowss
  twoRows <- compTwo(twoR, allRESPEN)
  
  # Three or more rows 
  moreRows <- compMore(moreR, allRESPEN)
  
  ### Combine all together and save ### ----
  
  cleanData <- rbindlist(list(cleanData, twoRows, moreRows), fill = T)
  
  # Save in 
  fwrite(cleanData, <<<< FILEPATH REDACTED >>>>)
  message(paste0("REMO ESPEN Deduplicated and cleaned data saved here: ", <<<< FILEPATH REDACTED >>>>))
  
  return(cleanData)
}

prepOnchoRESPEN(root)

