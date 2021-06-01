###################################
# POST UBCOV EXTRACTION DATA CLEANING FOR GEOSPATIAL DATA EXTRACTIONS & GEOGRAPHY MATCHING

# INSTRUCTIONS:

# IF ANY SURVEYS HAVE BEEN ADDED SINCE THE LAST RUN:
# If the survey does not ask the heard about sti question, add its nid to missing_heard_nids
# If the survey asks heard about sti but NOT had sti, add it to heard_not_had_nids
##########################################

## Setup----------------------------------------------------------------------
rm(list=ls())

topic <- "sti_symptoms"
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


say_file_and_read_dta <- function(file){
  #FOR UBCOV EXTRACTIONS. READ THEM IN AND DROP COLUMNS THAT DON'T COME FROM YOUR EXTRACTIONS
  message(file)
  dta <- try(read_dta(file))
  # drop all the rsp columns that come from the topic merge in extraction except had_intercourse
  bye <- c("partner_away", "debut_age", "condom_last_time", "condom_every_time_3_partner",
           "client_sex_worker_3_partner", "client_sex_worker_2_partner", "num_partners_year",
           "num_partners_lifetime", "paid_for_sex_ever_m", "paid_for_sex_year_m", "condom_paid_sex",
           "client_sex_worker", "last_sex_years", "paid_for_sex_year_w", "age_at_first_union",
           "age_at_first_union_year", "marital_status")
  dta <- dta[, !(names(dta) %in% bye)]
  if ("data.frame" %in% class(dta)){
    return(dta)
  }
}



## Read data and bind together----------------------------------------------------------------------
extractions <- list.files(folder_in, full.names=T)


# Bind all extracted surveys together
top <- data.frame()
for (i in 1:length(extractions)){
  dta <- say_file_and_read_dta(extractions[i])
  top <- rbind.fill(top, dta)
  message(paste("processing file ", i, " of ", length(extractions) ))
  try(nrow(top$currently_married))
}

topics <- top


## Topic specific modifications----------------------------------------------------------------------

# OVERWRITE HAD STI TO 0 IF HEARD STI IS 0 (don't do this for surveys that ask heard but not had)
heard_not_had_nids <- c(22114, 22116, 325046, 324443)
topics$had_sti <- ifelse(topics$heard_about_stis == 0 & !(topics$nid %in% heard_not_had_nids), 0, topics$had_sti)

# If discharge and sore/ulcer asked as one question, keep. Otherwise build from individual discharge and sore questions:
# NA if either discharge or sore is missing. If neither missing, OR them together
topics$discharge_or_sore <- as.numeric(ifelse(is.na(topics$discharge_or_sore),
                                              ifelse(is.na(topics$genital_discharge) | is.na(topics$genital_sore), NA, (topics$genital_discharge | topics$genital_sore)),
                                              topics$discharge_or_sore))


# DROP ALL HEARD_ABOUT_STI IS N/A (except for svys that didn't ask), HAD INTERCOURSE IS NOT 1
missing_heard_nids <- c(27924, 27952, 27987, 59339, 80790, 133304, 134753, 153643, 228102, 287630, 287629, 450419)
topics <- subset(topics, !is.na(heard_about_stis | nid %in% missing_heard_nids))


## Read in and prepare geocodebooks----------------------------------------------------------------------
source("<<<< FILEPATH REDACTED >>>>")
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


## Save RDS file of matched extractions----------------------------------------------------------------------
saveRDS(all, file=paste0(folder_out, "/", topic, "_", module_date, ".RDS"))


## Save csv of unmatched surveys----------------------------------------------------------------------
gnid <- unique(geo$nid)
fix <- subset(all, !(all$nid %in% gnid))
fix <- as.data.frame(fix)
fix_collapse <- unique(fix[c("nid", "country", "year_n", "survey_series")])
fix_collapse$country <- strsplit(fix_collapse$country, 1, 3)
if (nrow(fix_collapse) > 0) {
  write.csv("<<<< FILEPATH REDACTED >>>>")
}

