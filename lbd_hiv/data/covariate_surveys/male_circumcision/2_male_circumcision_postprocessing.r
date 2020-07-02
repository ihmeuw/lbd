###################################
# POST UBCOV EXTRACTION DATA CLEANING FOR GEOSPATIAL DATA EXTRACTIONS & GEOGRAPHY MATCHING
# PIONEERED BY ANNIE BROWNE
# UPDATED & OVERHAULED BY MANNY GARCIA
# EMAIL ABROWNE@WELL.OX.AC.UK
# EMAIL GMANNY@UW.EDU

# INSTRUCTIONS: 
# UBCOV OUTPUTS MUST BE SAVED IN LIMITED USE DIRECTORY
###################################

## Setup----------------------------------------------------------------------


rm(list=ls())

topic <- "male_circumcision"
folder_in <- paste0("<<<< FILEPATH REDACTED >>>>")
folder_out <- paste0("<<<< FILEPATH REDACTED >>>>")

if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}
p_load(haven, stringr, plyr, data.table, magrittr, parallel, doParallel, rgdal)

options(warn=-1)
module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)

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
extractions <- list.files(folder_in, full.names=T)

#bind all extracted surveys together
top <- data.frame()
for (i in 1:length(extractions)){
  dta <- read_dta(extractions[i])
  top <- rbind.fill(top, dta)
  message(paste("processing file ", i, " of ", length(extractions) ))
  try(nrow(top$male_circumcision))
}

topics <- top

##set names that changed in new extraction process
setnames(topics, "int_year", "year")
setnames(topics, "year_end", "end_year")


## Topic specific modifications----------------------------------------------------------------------

# NID 27987 has additional clusters appended that don't correspond to interviews
topics <- topics[!(topics$nid == 27987 & topics$psu > 105), ]

#Assign sex_id to 411301 MN mod 
#topics$sex_id[topics$nid == "411301" & is.na(topics$sex_id)] <- "1"
# Some of the svys ask circumcision questions of both men and women so this is dropping 
# responses from women so it's just males
topics <- subset(topics, sex_id == 1)
topics <- data.table(topics)

## Read in and prepare geocodebooks----------------------------------------------------------------------
source("/lbd_core/data_central/geocodebook_functions.R")
geo_k <- get_geocodebooks(nids = unique(topics$nid))


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


all[all$nid == '411301']


## Save RDS file of matched extractions----------------------------------------------------------------------
saveRDS(all, file=paste0(folder_out, "/", topic, "_", module_date, ".RDS"))


## Save csv of unmatched surveys----------------------------------------------------------------------
# gnid <- unique(geo$nid)
# fix <- subset(all, !(all$nid %in% gnid))
# fix <- as.data.frame(fix)
# fix_collapse <- unique(fix[c("nid", "country", "year", "survey_name")])
# fix_collapse$country <- strsplit(fix_collapse$country, 1, 3)
# if (nrow(fix_collapse) > 0) {
#   write.csv(fix_collapse, file=paste0(folder_out, "/", module_date, "_", topic, "_", "geography_matching_to_do", ".csv"), row.names = F)
# }


