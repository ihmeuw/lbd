###################################
# POST UBCOV EXTRACTION DATA CLEANING FOR GEOSPATIAL DATA EXTRACTIONS & GEOGRAPHY MATCHING
###################################

## Setup----------------------------------------------------------------------
rm(list=ls())


topic <- "rsp"
folder_in <- "<<<< FILEPATH REDACTED >>>>"
folder_out <- "<<<< FILEPATH REDACTED >>>>"

if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}
p_load(haven, stringr, plyr, data.table, magrittr, parallel, doParallel, rgdal)

options(warn=-1)
module_date <- format(Sys.Date(), "%Y_%m_%d")

read_add_name_col <- function(file){
  #FOR GEOGRAPHY CODEBOOKS. READS THEM IN AND ADDS A COLUMN WITH THEIR CORRESPONDING SURVEY_SERIES
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file, integer64="character", data.table=T)
  df[, survey_series := svy]
  df <- lapply(df, as.character)
  return(df)
}


## Read data and bind together----------------------------------------------------------------------
extractions <- list.files(folder_in, full.names=T, pattern="*.dta")

#bind all extracted surveys together
top <- data.frame()
for (i in 1:length(extractions)){
  dta <- read_dta(extractions[i])
  top <- rbind.fill(top, dta)
  message(paste("processing file ", i, " of ", length(extractions) ))
  try(nrow(top$currently_married))
}

topics <- top


## Read in and prepare geocodebooks----------------------------------------------------------------------
message("get all geography codebooks")
files <- list.files("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T)
message("lapply geographies into list")
geogs <- lapply(files, read_add_name_col)
message("rbind geography codebooks together")
geo <- rbindlist(geogs, fill=T, use.names=T)

#dedupe the geography codebook by geospatial_id and nid
setkey(geo, nid, iso3, geospatial_id)
geo <- unique(geo, use.key=T)
#make everything the same data type to merge
message("make types between merging datasets match")
if (class(topics$nid) == "numeric"){
  geo[, nid := as.numeric(nid)]
} else if (class(topics$nid) == "character"){
  geo[, nid := as.character(nid)]
} else{
  message("update code to accomodate topics nid as")
  message(class(topics$nid))
}

if (class(topics$geospatial_id) == "numeric"){
  geo[, geospatial_id := as.numeric(geospatial_id)]
} else if (class(topics$geospatial_id) == "character"){
  geo[, geospatial_id := as.character(geospatial_id)]
} else{
  message("update code to accomodate topics geospatial_id as")
  message(class(topics$geospatial_id))
}
geo_keep <- c("nid", "iso3", "geospatial_id", "point", "lat", "long", "shapefile", "location_code", "location_name", "admin_level", "survey_series")
geo_k <- geo[, geo_keep, with=F]


## Merge extracted data & geocodebooks----------------------------------------------------------------------
all <- merge(geo_k, topics, by.x=c("nid", "iso3", "geospatial_id"), by.y=c("nid", "ihme_loc_id", "geospatial_id"), all.x=F, all.y=T)

missing_shapefile <- subset(all, is.na(point))
missing_shapefile$shapefile <- NA
all <- subset(all, !(is.na(point)))
all <- rbind(all, missing_shapefile, fill=T)

setnames(all, c("iso3", "lat", "long"), c("country", "latitude", "longitude"))
all$ihme_loc_id <- all$country
all$country <- substr(all$country, 1, 3)


## Save RDS file of matched extractions----------------------------------------------------------------------
saveRDS(all, file=paste0(folder_out, "/", topic, "_", module_date, ".RDS"))


