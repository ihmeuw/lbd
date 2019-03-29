## ########################################################################################################################################
##
## HIV_TEST POST UBCOV EXTRACTION DATA CLEANING & GEOGRAPHY MATCHING
## 
## Purpose: Prepare microdatasets extracted from UbCov for collapse.
##            Bind all extratced survey microdata, 
##            merge with point and polygon information in the Geography codebooks,
##            and format and subset to prepare data for collpase.
##          
## ########################################################################################################################################


## SET_UP ENVIRONMENT ---------------------------------------------------------------------------------------------------------------------

rm(list=ls())

topic <- "hiv"
folder_in <- "<<<< FILEPATH REDACTED >>>>"
folder_out <- "<<<< FILEPATH REDACTED >>>>"


if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}
p_load(haven, stringr, plyr, data.table, magrittr, parallel, doParallel)

options(warn=-1)
module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)

## DEFINE FUNCTIONS ------------------------------------------------------------------------------------------------------------------------

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

all_to_char_df <- function(df){
  df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
  return(df)
}

say_file_and_read_dta <- function(file, keep){
  #FOR READING IN UBCOV EXTRACTIONS
  message(file)
  dta <- try(read_dta(file))
  if ("data.frame" %in% class(dta)){
    return(dta)
  }
}

## BIND ALL GEOGRAPHY CODEBOOKS ------------------------------------------------------------------------------------------------------------

##Get all geo codebooks and package them together

message("get all geography codebooks")
files <- list.files("<<<< FILEPATH REDACTED >>>>", pattern=".csv$", ignore.case = T, full.names = T)
files <- grep("Copy|together|linkage|IPUMS|special", files, value = T, invert = T) #IPUMS is handled separately
message("lapply geographies into list")
geogs <- lapply(files, read_add_name_col)
message("rbind geography codebooks together")
geo <- rbindlist(geogs, fill=T, use.names=T)
rm(geogs)

# dedupe the geography codebook by geospatial_id and nid
setkey(geo, nid, geospatial_id, iso3)
geo <- unique(geo, use.key=T)

## BIND EXTRACTED UBCOV FILES ---------------------------------------------------------------------------------------------------------------

extractions <- list.files(folder_in, full.names=T)

#bind all extracted surveys together
top <- data.frame()
for (i in 1:length(extractions)){
  dta <- read_dta(extractions[i])
  top <- rbind.fill(top, dta)
  message(paste("processing file ", i, " of ", length(extractions) ))
}

topics <- top
  

## MERGE GEO-CODEBOOKS AND SURVEY DATA ------------------------------------------------------------------------------------------------------

# Prepare for merge by aligning data type of `topics` and `geo`
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
geo_keep <- c("nid", "iso3", "geospatial_id", "point", "lat", "long",
              "shapefile", "location_code", "location_name", "admin_level", "survey_series")
geo_k <- geo[, geo_keep, with=F]

# Merge extracted dataset and geography dataset together
all <- merge(geo_k, topics, by.x=c("nid", "iso3", "geospatial_id"), 
                            by.y=c("nid", "ihme_loc_id", "geospatial_id"), all.x=F, all.y=T)

missing_geography <- subset(all, is.na(point))
missing_geography$shapefile <- NA
all <- subset(all, !(is.na(point)))
all <- rbind(all, missing_geography, fill=T)
rm(missing_shapefile)

## SUBSETTING -------------------------------------------------------------------------------------------------------------------------------

# Drop lines with rejected/missing blood samples
all <- subset(all, all$rejected_missing == 0 | is.na(all$rejected_missing))

# Drop people known to not be tested
all <- subset(all, !(all$tested_hiv == 0 & is.na(all$hiv_test)))

# Drop rows where missing hiv result is accounted for (via reports, etc.)
accounted_nids <- c(12630, 20145, 313076, 19557, 19627, 324443, 325046)
all <- subset(all, !(all$nid %in% accounted_nids & is.na(all$hiv_test)))

# Change int_year for ZMB 2001 and ZAF 2004 to better reflect survey
all$int_year[all$nid == 21097] <- rep(2002,length(all$int_year))
all$int_year[all$nid == 313074] <- rep(2005,length(all$int_year))

# Make adjustments for ACDIS (nid: 11780) (have 'year' reflect round and remove 2001 and 2002 years [see cohort documentation])
all <- subset(all, !(nid == 11780 & int_year %in% c(2001, 2002))) 
all$year[all$nid == 11780] <- all$int_year[all$nid == 11780]
all$end_year[all$nid == 11780] <- all$int_year[all$nid == 11780]
all$year[all$nid == 11780 & all$int_year == 2004] <- 2003
all$end_year[all$nid == 11780 & all$int_year ==2003] <- 2004
# Assign weights as 1
all$pweight[all$nid == 11780] <- 1
all$hiv_weight[all$nid == 11780] <- 1

# Rename columns to prep for collapse code
colnames(all)[which(names(all) == "lat")] <- "latitude"
colnames(all)[which(names(all) == "long")] <- "longitude"
colnames(all)[which(names(all) == "iso3")] <- "country"
colnames(all)[which(names(all) == "geospatial_id")] <- "cluster_id"

## SAVE DATA ---------------------------------------------------------------------------------------------------------------------------------
# Save .Rdata file of matched extractions
saveRDS(all, file=paste0(folder_out, "/", module_date, ".RDS"))

