#
# SCRIPT FOR ONCHO OCP DATA PREP
#
# Purpose: All Oncho OCP Data combined/cleaned/hard code fixes of data
# NOTE: Some OCP data lives on L drive, but here we are
###################################################################################################

### Necessary Packages ## ----
library(data.table)
library(readxl)

### Example Call ### ----
# You should just run this whole thing like a script

###################################################################################################
### FUNCTIONS #####################################################################################


### PREP OCP DATA ### ------------------------------------------------------------

prepOCPnotL <- function(root) {

  ocp <- as.data.table(read_excel(<<<< FILEPATH REDACTED >>>>))
  
  # Update column names
  colnames(ocp) <- as.matrix(ocp[2, ])[1, ]

  # Still has column row in the data itself
  ocp <- ocp[Basin != "Basin", ]

  setnames(ocp, old=c("COUNTRY NAME")
           , new=c("countryName") )

  # Determine ISO3
  ocp[countryName == "TOGO", iso3 := "TGO"]
  ocp[countryName == "BENIN", iso3 := "BEN"]
  ocp[countryName == "COTE D'IVOIRE" | countryName == "COTE IVOIRE", iso3 := "CIV"]
  ocp[countryName == "GHANA", iso3 := "GHA"]
  ocp[countryName == "BURKINA FASO", iso3 := "BFA"]
  ocp[countryName == "SIERRA LEONE", iso3 := "SLE"]
  ocp[countryName == "GUINEE CONAKRY", iso3 := "GIN"]

  # Remove
  ocp <- ocp[!is.na(`VILAGE LATITUDE`) & !is.na(`VILLAGE LONGITUDE`), ]

  ### Grouping ### ----

  ocp[, LBD_group := .GRP, by = list(`VILAGE LATITUDE`, `VILLAGE LONGITUDE`, YEARS, MONTHS)]
  groupCounts <- ocp[, .(count = .N), by = LBD_group]

  newOCP <- ocp[LBD_group %in% groupCounts[count == 1]$LBD_group, ]
  ocp <- ocp[LBD_group %in% groupCounts[count > 1]$LBD_group, ]

  newOCP <- rbindlist(list(newOCP
                           , ocp[LBD_group == 47][1]
                           , ocp[LBD_group == 48][1]
                           , ocp[LBD_group == 103][1]), fill = T)

  newOCP[, LBD_uniq := 1]
  newOCP[, LBD_spec := ""]


  ocpMergePrep <- data.table(
    row_id = newOCP$row_id,
    nid = 364663,
    country = newOCP$countryName,
    iso3 = newOCP$iso3,
    source = "OCP",
    year = newOCP$YEARS,
    month = newOCP$MONTHS,
    location_name = newOCP$`VILLAGE NAME`,
    examined = newOCP$Exa,
    positive = newOCP$Pos,
    latitude = newOCP$`VILAGE LATITUDE`,
    longitude = newOCP$`VILLAGE LONGITUDE`,
    shapefile = NA,
    poly_id = NA,
    dx_code = 9,
    dx = "ss",
    sex = "Both",
    MDA_status = 0,
    LBD_group = newOCP$LBD_group,
    LBD_uniq = newOCP$LBD_uniq,
    LBD_spec = newOCP$LBD_spec
  )

  return(ocpMergePrep)

}

### PREP OTHER OCP DATA ### --------------------------------------------------------

prepOCPonL <- function(root) {

  ### Initial pull and group numbering ### ----

  # Temp loc to pull from
  ocp <- as.data.table(fread(<<<< FILEPATH REDACTED >>>>))
  
  # Determine ISO3
  ocp[PAYNOM == "TOGO", iso3 := "TGO"]
  ocp[PAYNOM == "BENIN", iso3 := "BEN"]
  ocp[PAYNOM == "COTE D'IVOIRE" | PAYNOM == "COTE IVOIRE", iso3 := "CIV"]
  ocp[PAYNOM == "GHANA", iso3 := "GHA"]
  ocp[PAYNOM == "BURKINA FASO", iso3 := "BFA"]
  ocp[PAYNOM == "SIERRA LEONE", iso3 := "SLE"]
  ocp[PAYNOM == "GUINEE CONAKRY", iso3 := "GIN"]
  ocp[PAYNOM == "NIGER", iso3 := "NGR"]
  ocp[PAYNOM == "MALI", iso3 := "MLI"]
  ocp[PAYNOM == "SENEGAL", iso3 := "SEN"]

  # group numbering
  ocp[, LBD_group := .GRP, by = list(VILATI, VILONG, PHANNE, PHAMOIS)]

  ### Remove not single rows ### ----

  groupCounts <- ocp[, .(count = .N), by = list(LBD_group)]
  mOCP <- ocp[LBD_group %in% groupCounts[count > 1]$LBD_group, ]
  ocp <- ocp[!(LBD_group %in% groupCounts[count > 1]$LBD_group), ]

  # Hardcode back in one row each
  ocp <- rbindlist(list(ocp, mOCP[LBD_group == 7][1]))
  ocp <- rbindlist(list(ocp, mOCP[LBD_group == 436][1]))

  ### Loop through each row and separating out for m/f and nod, acute lesion, blindness ### ----
  newOCP <- data.table()

  for (i in 1:nrow(ocp)) {

    curr <- ocp[i]

    dxList <- c("ss", "ss", "eye", "eye", "eye", "eye")
    dxCList <- c(9, 9, 26, 26, 25, 25)
    sexList <- c("Male", "Female")
    posList <- c(curr$Posm, curr$Posf, curr$Acum, curr$Acuf, curr$Blindf, curr$Blindm)
    exaList <- c(curr$Exam, curr$Exaf)

    newOCPMerge <- data.table(
      row_id = curr$row_id,
      nid = 391071,
      country = curr$PAYNOM,
      iso3 = curr$iso3,
      source = "OCP",
      year = curr$PHANNE,
      month = curr$PHAMOIS,
      location_name = curr$VILNOM,
      examined = exaList,
      positive = posList,
      latitude = curr$VILATI,
      longitude = curr$VILONG,
      shapefile = NA,
      poly_id = NA,
      dx_code = dxCList,
      dx = dxList,
      sex = sexList,
      MDA_status = 0, # pre-control implies no MDA
      LBD_group = curr$LBD_group,
      LBD_uniq = 1,
      LBD_spec = "DX Sex "
    )

    #
    if (nrow(newOCP) == 0) {
      newOCP <- newOCPMerge
    } else {
      newOCP <- rbindlist(list(newOCP, newOCPMerge))
    }
  }

  return(newOCP)

}


###################################################################################################
### RUNNING #######################################################################################

# For other script
if (Sys.info()[1] == "Windows") { root <- <<<< FILEPATH REDACTED >>>> } else { root <- <<<< FILEPATH REDACTED >>>>}

prepOnchoOCP <- function(root) {

  ocp1 <- prepOCPnotL(root)
  ocp2 <- prepOCPonL(root)

  # Adjust group numbers
  ocp2[, LBD_group := LBD_group + max(ocp1$LBD_group)]

  fullOCP <- rbindlist(list(ocp1, ocp2))

  # Saving
  fwrite(fullOCP, <<<< FILEPATH REDACTED >>>>)
  message(paste0("OCP cleaned data saved here: ", <<<< FILEPATH REDACTED >>>>))

  return(fullOCP)

}


prepOnchoOCP(root)
