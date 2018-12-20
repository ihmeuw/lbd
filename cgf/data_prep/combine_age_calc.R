### Master Data Prep to combine all surveys for geospatial data prep

# initiate timing
start_time <- Sys.time()

root <- ifelse(Sys.info()[1]=="Windows", <<<< FILEPATH REDACTED >>>>>, <<<< FILEPATH REDACTED >>>>>)
## Load libraries and miscellaneous MBG project functions.

repo <- '<<<<< FILEPATH REDACTED >>>>>'
setwd(repo)

# set directories
data_dir <- <<<< FILEPATH REDACTED >>>>>
mbg_dir  <- <<<< FILEPATH REDACTED >>>>>
save_dir <- <<<< FILEPATH REDACTED >>>>>
root_dir <- <<<< FILEPATH REDACTED >>>>>


package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                      paste0(root, <<<< FILEPATH REDACTED >>>>>),
                      paste0(root, <<<< FILEPATH REDACTED >>>>>))

.libPaths(package_lib)                  

## load functions
source('mbg_central/mbg_functions.R')                   
source('mbg_central/prep_functions.R')                  
source('mbg_central/covariate_functions.R')             
source('mbg_central/misc_functions.R')                  
source('mbg_central/post_estimation_functions.R')
source('mbg_central/gbd_functions.R')
source('mbg_central/shiny_functions.R')
source('mbg_central/holdout_functions.R')
source('mbg_central/polygon_functions.R')
source('mbg_central/collapse_functions.R')
source('mbg_central/seegMBG_transform_functions.R')     

package_list <- c('survey', 'pbapply', 'readstata13', 'foreign',
                  'rgeos', 'data.table','raster','rgdal','INLA',
                  'seegSDM','seegMBG','plyr','dplyr', 'foreach',
                  'doParallel')

for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

## Load iso3 to country name map
iso3_to_country <- fread(paste0(root_dir,"ref/iso3_to_country.csv"))

##### 1. ################################################################################################################
## Read all data 

## ######################################################################################################
# load, prep and format DHS/MICS
dhs_mics <- read.csv(paste0(root,<<<< FILEPATH REDACTED >>>>>), stringsAsFactors = FALSE)
names(dhs_mics)[names(dhs_mics)=="cluster_number"] <- "psu"
names(dhs_mics)[names(dhs_mics)=="survey_series"] <- "source"
names(dhs_mics)[names(dhs_mics)=="iso3"] <- "country"
names(dhs_mics)[names(dhs_mics)=="lat"] <- "latitude"
names(dhs_mics)[names(dhs_mics)=="long"] <- "longitude"
names(dhs_mics)[names(dhs_mics)=="sample_weight"] <- "pweight"
names(dhs_mics)[names(dhs_mics)=="sex_id"] <- "sex" 
names(dhs_mics)[names(dhs_mics)=="year_end"] <- "end_year" 
names(dhs_mics)[names(dhs_mics)=="height"] <- "child_height" 
names(dhs_mics)[names(dhs_mics)=="weight"] <- "child_weight" 
dhs_mics$master.ind <- 1:nrow(dhs_mics)

## then divide days by 7.0192 to get weeks (365/52=7.0192)
dhs_mics$age_wks <- NA
dhs_mics$age_wks <- as.numeric(dhs_mics$age_wks)
dhs_mics$age_wks <- (dhs_mics$age_day)/7.0192
dhs_mics$age_wks[is.na(dhs_mics$age_wks)] <- (dhs_mics$age_month[is.na(dhs_mics$age_wks)]*4.333)
dhs_mics$age_wks <- as.integer(round_any(dhs_mics$age_wks, 1)) 

dhs_mics$orig.psu <- dhs_mics$psu
dhs_mics$cluster_number <- dhs_mics$psu 
dhs_mics$psu <- dhs_mics$geospatial_id

na.geo <- aggregate(is.na(geospatial_id) ~ nid, dhs_mics, mean)
na.psu <- aggregate(is.na(orig.psu) ~ nid, dhs_mics, mean)
colnames(na.geo)[2] <- "nas.g"
colnames(na.psu)[2] <- "nas.p"
nas <- merge(na.geo, na.psu)
use.psu <- nas$nid[which(nas$nas.g > nas$nas.p)]
psu.r <- which(dhs_mics$nid %in% use.psu)
dhs_mics[psu.r, 'psu'] <- dhs_mics[psu.r, 'orig.psu']

na.pweight  <- aggregate(is.na(pweight) ~ nid, dhs_mics, mean)
colnames(na.pweight)[2] <- "nas.p"
na.hhweight <- aggregate(is.na(hhweight) ~ nid, dhs_mics, mean)
colnames(na.hhweight)[2] <- "nas.h"
nas <- merge(na.hhweight, na.pweight)
use.hh <- nas$nid[which(nas$nas.h < nas$nas.p)]
hh.r <- which(dhs_mics$nid %in% use.hh)
dhs_mics[hh.r, 'pweight'] <- dhs_mics[hh.r, 'hhweight']

## set variables to keep and subset
setDF(dhs_mics)
vars_to_keep <- c("nid", "psu", "source", "country", "start_year",
                  "end_year", "pweight", "strata", "sex", "age_wks",
                  "birth_weight", "child_weight", "child_height",
                  "point", "latitude", "longitude", "uncertain_point",
                  "buffer", "location_name", "location_code",
                  "admin_level", "shapefile", "int_month", "int_year",
                  "hh_id", "birth_weight_card", "geospatial_id", "cluster_number", "master.ind")
dhs_mics <- dhs_mics[vars_to_keep]

## format and save for main polygon resampling
all_data <- dhs_mics 

all_data$age_mo <- all_data$age_wks/4.333
all_data$age_mo <- as.integer(round_any(all_data$age_mo, 1))
all_data$age_wks <- as.integer(round_any(all_data$age_wks, 1))

all_data$age_cat_1 <- ifelse(all_data$age_wks <= 104, "0-2", "2-5")
all_data$age_cat_2 <- "0-5"

all_data <- subset(all_data, !is.na(sex))

all_data <- subset(all_data, all_data$age_wks >= 0)
all_data <- subset(all_data, all_data$age_wks <= 260)

## ensure that we don't have subnational codes in country ISOs: e.g. KEN_351235 -> KEN
all_data$country <- substr(all_data$country, 1, 3) ## take the first 3 letters only

write.csv(all_data, paste0(root, <<<< FILEPATH REDACTED >>>>>))
