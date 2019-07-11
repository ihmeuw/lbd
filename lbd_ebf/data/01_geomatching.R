##########################################################################################
# POST UBCOV EXTRACTION DATA CLEANING FOR GEOSPATIAL DATA EXTRACTIONS & GEOGRAPHY MATCHING
##########################################################################################

## SETUP----------------------------------------------------------------------------------
rm(list = ls())

#Define values
topic <- "breastfeeding"
nid_vec <- c() # list the NIDs of surveys (if any) to merge using admin_1 column
cluster <- TRUE # cluster true/false
cores <- 5

#setup
folder_in <- "<<<< FILEPATH REDACTED >>>>"
folder_out <- "<<<< FILEPATH REDACTED >>>>"
options(warn = -1)
module_date <- "<<<< FILEPATH REDACTED >>>>"
module_date <- gsub("-", "_", module_date)

#Load packages 
if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}
p_load(haven, stringr, plyr, data.table, magrittr, parallel, doParallel, feather)

read_add_name_col <- function(file){
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file, integer64="character", data.table=T)
  df[, survey_series := svy]
  df <- lapply(df, as.character)
  return(df)
}

all_to_char_df <- function(df){
  df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
  return(df)
}

extractions <- list.files(folder_in, full.names = T, pattern = ".dta", recursive = F)
extractions <- grep("IPUMS_CENSUS", extractions, invert = T, value = T) # IPUMS is handled separately

if(cluster == TRUE) {
  message("cluster is set to TRUE")
  message("Make cluster")
  cl <- makeCluster(cores)
  message("Register cluster")
  registerDoParallel(cl)
  message("Start foreach")
  
  top <- foreach(i = 1:length(extractions), .packages = c("haven")) %dopar% {
    dta <- read_dta(extractions[i], encoding = "latin1")
    return(dta)
  }
  message("Foreach finished")
  message("Closing cluster")
  stopCluster(cl)
} else if(cluster == FALSE) {
  message("cluster is set to FALSE")
  top <- foreach(i = 1:length(extractions)) %do% {
    message(paste0("Reading in: ", extractions[i]))
    dta <- read_dta(extractions[i], encoding = "latin1")
    return(dta)
  }
}

message("rbindlist all extractions together")
topics <- rbindlist(top, fill=T, use.names=T)

message("Retrieve geo codebook filepaths")
files <- list.files("<<<< FILEPATH REDACTED >>>>", pattern = ".csv$", ignore.case = T, full.names = T)
files <- grep("IPUMS|special", files, value = T, invert = T) 

message("Read geo codebooks into list")
geogs <- lapply(files, read_add_name_col)

message("Bind geo codebooks together")
geo <- rbindlist(geogs, fill = T, use.names = T)
rm(geogs)

setkey(geo, nid, geospatial_id)
geo <- unique(geo, use.key=T)

topics$geospatial_id <- as.character(topics$geospatial_id)
geo$geospatial_id <- as.character(geo$geospatial_id)

geo_keep <- c("nid", "iso3", "geospatial_id", "point", "lat", "long", "shapefile", "location_code", "survey_series")
geo_k <- geo[, geo_keep, with=F]

topics[, geospatial_id := gsub("[^[:alnum:] | _ ]", "", geospatial_id)]
geo_k[, geospatial_id := gsub("[^[:alnum:] | _ ]", "", geospatial_id)]

message("Merge ubCov outputs & geo codebooks together")
all <- merge(geo_k, topics, by.x=c("nid", "iso3", "geospatial_id"), by.y=c("nid", "ihme_loc_id", "geospatial_id"), all.x=F, all.y=T)

message("Saving as .Rdata")
save(all, file=paste0(folder_out, "/", module_date, ".Rdata"))

