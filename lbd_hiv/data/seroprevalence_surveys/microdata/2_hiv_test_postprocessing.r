## ########################################################################################################################################
##
## HIV_TEST POST UBCOV EXTRACTION DATA CLEANING & GEOGRAPHY MATCHING
## ########################################################################################################################################


## SET_UP ENVIRONMENT ---------------------------------------------------------------------------------------------------------------------

rm(list=ls())

topic <- "hiv"
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

# bind all extracted surveys together
top <- data.frame()
for (i in 1:length(extractions)){
  dta <- read_dta(extractions[i])
  top <- rbind.fill(top, dta)
  message(paste("processing file ", i, " of ", length(extractions) ))
}



topics <- top
rm(top, dta)
save(topics, file=paste0("<<<< FILEPATH REDACTED >>>>"))



# MERGE GEO-CODEBOOKS AND SURVEY DATA -----------------------------------------------------------------------------------------------------

# Read in and prepare geocodebooks----------------------------------------------------------------------
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



## SUBSETTING -------------------------------------------------------------------------------------------------------------------------------

# Dropping all NGA_NARHS data for NAIIS reports when available
all <- subset(all, !(all$survey_name == "NGA_NARHS"))

#Dropping DHS Special Stratified survey that masks DOM 2013 MACRO_DHS
all <- subset(all, !(all$nid == 165645))


# Drop lines with rejected/missing blood samples
all <- subset(all, all$rejected_missing == 0 | is.na(all$rejected_missing))

# Drop people known to not be tested - including where hiv_test and tested_hiv NA (via microdata)
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

# Assign weights as 1 (study attempted tio account for everyone in study region)
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

# Subset out and save csv of unmatched geographies
gnid <- unique(topics$nid)
fix <- subset(all, !(all$nid %in% gnid))
fix <- as.data.frame(fix)
fix_collapse <- unique(fix[c("nid", "country", "year", "survey_name")])
fix_collapse$country <- strsplit(fix_collapse$country, 1, 3)
write.csv(fix_collapse, file=paste0("<<<< FILEPATH REDACTED >>>>"))

# For GBD use, remove `Limited Use` data
gbd_all <- subset(all, !(nid %in% c(13084, 11780, 425256, 287630, 287629)))
save(gbd_all, file=paste0(folder_out_gbd, "/", module_date,".Rdata"))
saveRDS(gbd_all, file=paste0(folder_out_gbd, "/", module_date,".RDS"))

