## ########################################################################################################################################
##
## HIV_INCIDENCE POST UBCOV EXTRACTION DATA CLEANING & GEOGRAPHY MATCHING
##          
## ########################################################################################################################################


## SET_UP ENVIRONMENT ---------------------------------------------------------------------------------------------------------------------

rm(list=ls())


topic <- "hiv_incidence"
folder_in <- paste0("<<<< FILEPATH REDACTED >>>>")
folder_out <- paste0("<<<< FILEPATH REDACTED >>>>")
folder_out_gbd <- paste0("<<<< FILEPATH REDACTED >>>>")

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
rm(top, dta)
save(topics, file=paste0(folder_out, "/topics_no_geogs_", module_date, ".Rdata"))



# MERGE GEO-CODEBOOKS AND SURVEY DATA -----------------------------------------------------------------------------------------------------
# Read in and prepare geocodebooks
# Enter your sys username!
source("<<<< FILEPATH REDACTED >>>>")
geo_k <- get_geocodebooks(nids = unique(topics$nid))
geo_k$iso3 <- str_trim(geo_k$iso3)
geo_k$nid <- as.numeric(geo_k$nid)


# Merge extracted dataset and geography dataset together
all <- merge(geo_k, topics, by.x=c("nid", "iso3", "geospatial_id"),
             by.y=c("nid", "ihme_loc_id", "geospatial_id"), all.x=F, all.y=T)

missing_geography <- subset(all, is.na(point))
missing_geography$shapefile <- NA
all <- subset(all, !(is.na(point)))
all <- rbind(all, missing_geography, fill=T)
rm(missing_shapefile)

setnames(all, c("iso3", "lat", "long"), c("country", "latitude", "longitude"))
all$ihme_loc_id <- all$country
all$country <- substr(all$country, 1, 3)


## SUBSETTING & Cleaning -------------------------------------------------------------------------------------------------------------------------------

# Drop lines with rejected/missing blood samples
all <- subset(all, all$rejected_missing == 0 | is.na(all$rejected_missing))

# Drop people known to not be tested - where hiv_test is NA (via microdata)
all <- subset(all, !(is.na(all$hiv_test)))

# Rename columns to prep for collapse code
colnames(all)[which(names(all) == "lat")] <- "latitude"
colnames(all)[which(names(all) == "long")] <- "longitude"
colnames(all)[which(names(all) == "iso3")] <- "country"
colnames(all)[which(names(all) == "geospatial_id")] <- "cluster_id"



## SAVE DATA ---------------------------------------------------------------------------------------------------------------------------------
# Save .Rdata file of matched extractions
saveRDS(all, file=paste0(folder_out, "/", module_date, ".RDS"))

# Subset out and save csv of unmatched geographies
gnid <- unique(topics$nid)
fix <- subset(all, !(all$nid %in% gnid))
fix <- as.data.frame(fix)
fix_collapse <- unique(fix[c("nid", "country", "year_end", "year_start", "survey_name")])
fix_collapse$country <- strsplit(fix_collapse$country, 1, 3)
write.csv(fix_collapse, file=paste0("<<<< FILEPATH REDACTED >>>>"))