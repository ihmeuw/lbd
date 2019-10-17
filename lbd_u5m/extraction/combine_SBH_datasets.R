## ##############################################################
## Author: Michael Collison
## Inputs: Extracted, non-geomatched DTAs from Ubcov for SBH
## Outputs: Geomatched SBH dataset
## Description: This script combines extracted .DTA files 
##    and geomatches them using the geography codebook
## Dependencies: data.table, haven, plyr
## ##############################################################

################################################################
######################### Setup ################################
################################################################
source("<<<< FILEPATH REDACTED >>>>")

package_list <- c('data.table', 'dplyr', "haven")
load_R_packages(package_list)

## set the working directory
setwd("<<<< FILEPATH REDACTED >>>>")

update <- T

################################################################
#################### Combine Datasets ##########################
################################################################

## Ensure that all DTA or dta files are selected
file.list <- list.files(pattern = '*dta$')
file.list2 <- list.files(pattern = '*DTA$')
file.list <- c(file.list, file.list2)
rm(file.list2)

## Loop through and bind everything together
for(i in 1:length(file.list)) {
  file.name <- file.list[i]
  mydata <- data.table(read_dta(file.name))
  message(paste("processing file ", i, " of ", length(file.list) ))
  if(file.name == file.list[1]){
    combined <- mydata 
  } else{
    combined <- rbind(combined, mydata, fill=TRUE)
  }
}

######################################################################
#################### Helper Functions ################################
######################################################################

#FOR GEOGRAPHY CODEBOOKS. READS THEM IN AND ADDS A COLUMN WITH THEIR CORRESPONDING SURVEY_SERIES
read_add_name_col <- function(file){
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file, integer64="character", data.table=T)
  df[, survey_series := svy]
  df <- lapply(df, as.character)
  return(df)
}

#####################################################################
#################### Merge Geographies ##############################
#####################################################################

#Read in geographies
message("get all geography codebooks")
files <- list.files( "<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T)
message("lapply geographies into list")
geogs <- lapply(files, read_add_name_col)
message("rbind geography codebooks together")
geo <- rbindlist(geogs, fill=T, use.names=T)

#dedupe the geography codebook by geospatial_id and nid
setkey(geo, nid, geospatial_id)
geo <- unique(geo, use.key=T)

#match nid and geospatial_id data types between combined and geo files
message("make types between merging datasets match")
if (class(combined$nid) == "numeric"){
  geo[, nid := as.numeric(nid)]
} else if (class(combined$nid) == "character"){
  geo[, nid := as.character(nid)]
} else{
  message("update code to accomodate combined nid as")
  message(class(combined$nid))
}

if (class(combined$cluster_number) == "numeric"){
  geo[, geospatial_id := as.numeric(geospatial_id)]
} else if (class(combined$cluster_number) == "character"){
  geo[, geospatial_id := as.character(geospatial_id)]
} else{
  message("update code to accomodate combined geospatial_id as")
  message(class(combined$cluster_number))
}

#subset geography dataset to needed columns
geo_keep <- c("nid", "iso3", "geospatial_id", "point", "lat", "long", "shapefile", "location_code", "location_name", "admin_level", "survey_series")
geo_k <- geo[, geo_keep, with=F]

#merge extracted dataset and geography dataset together
sbh <- merge(geo_k, combined, by.x=c("nid", "iso3", "geospatial_id"), by.y=c("nid", "country", "cluster_number"), all.x=F, all.y=T)

#####################################################################
########################## Cleaning #################################
#####################################################################

#rename some columns for consistency with old datasets
setnames(sbh, c("iso3", "geospatial_id", "lat", "long"),
         c("country", "cluster_number", "latnum", "longnum"))

#columns to keep
sbh <- sbh[,c("nid",
              "country",
              "year",
              "source",
              "survey",
              "strata",
              "cluster_number",
              "hh_id",
              "line_id",
              "weight",
              "age_group_of_woman",
              "age_year",
              "ceb",
              "ced",
              "interview_date_cmc",
              "ceb_female",
              "ceb_male",
              "ced_female",
              "ced_male",
              "point",
              "latnum",
              "longnum",
              "location_name",
              "location_code",
              "admin_level",
              "shapefile"), 
              with=FALSE]

if(update){
  old_sbh <- data.table(readRDS("<<<< FILEPATH REDACTED >>>>"))
  nids <- unique(sbh$nid)
  old_sbh <- old_sbh[!(nid %in% nids),]
  sbh <- rbind(old_sbh, sbh, fill=T)
}

## Save it as an RDS file
saveRDS(sbh, file = "<<<< FILEPATH REDACTED >>>>")
