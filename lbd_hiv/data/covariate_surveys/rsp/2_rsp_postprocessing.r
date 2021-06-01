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
topics <- rbindlist(lapply(extractions, read_dta), fill = T)

#Special Case: fix iso3 code for 408008(cleaning before geomatching)
topics$ihme_loc_id[topics$nid == "408008" & topics$ihme_loc_id == "COD_KIN"] <- "COD"
topics$ihme_loc_id[topics$nid == "408008" & topics$ihme_loc_id == "COD_KON"] <- "COD"
topics$ihme_loc_id[topics$nid == "408008" & topics$ihme_loc_id == "COD_KONGO"] <- "COD"

            
## Read in and prepare geocodebooks----------------------------------------------------------------------
# Enter your sys username!
source("<<<< FILEPATH REDACTED >>>>")
geo_k <- get_geocodebooks(nids = unique(topics$nid))
geo_k$iso3 <- str_trim(geo_k$iso3)

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
all$nid  <- as.numeric(all$nid)

## Save RDS file of matched extractions----------------------------------------------------------------------
saveRDS(all, file=paste0(folder_out, "/", topic, "_", module_date, ".RDS"))


## Save csv of unmatched surveys----------------------------------------------------------------------
gnid <- unique(geo$nid)
fix <- subset(all, !(all$nid %in% gnid))
fix <- as.data.frame(fix)
fix_collapse <- unique(fix[c("nid", "country", "year", "survey_name")])
fix_collapse$country <- strsplit(fix_collapse$country, 1, 3)
if (nrow(fix_collapse) > 0) {
  write.csv(fix_collapse, file=paste0("<<<< FILEPATH REDACTED >>>>"), row.names = F)
}

