#####################################################################
# POST UBCOV EXTRACTION DATA CLEANING FOR GEOSPATIAL DATA EXTRACTIONS 
# & GEOGRAPHY MATCHING
#####################################################################
############################## SETUP ################################
#####################################################################
rm(list=ls())

#number of cores
cores <- 5

#Setup
folder_in <- '<<<< FILEPATH REDACTED >>>>'
folder_out <- '<<<< FILEPATH REDACTED >>>>'

today <- gsub("-", "_", Sys.Date())

#Load packages
package_list <- c('haven', 'stringr', 'data.table', 'dplyr', 'magrittr', 
                  'parallel', 'doParallel')
lapply(package_list, library, character.only = T)



#####################################################################
######################## DEFINE FUNCTIONS ###########################
#####################################################################
#Read in geo codebook, add column with corresponding survey series
read_add_name_col <- function(file){
  message(file)
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file)
  df <- as.data.table(df)
  df[, survey_series := svy]
  df <- lapply(df, as.character, stringsAsFactors = FALSE)
  return(df)
}

#####################################################################
######################## BIND UBCOV EXTRACTS ########################
#####################################################################
#Generate list of extraction filepaths
extractions <- list.files(folder_in, full.names=T, pattern = ".dta", 
                          ignore.case=T, recursive = F)

#append all ubcov extracts together
message("Make cluster")
cl <- makeCluster(cores)
clusterEvalQ(cl, .libPaths('<<<< FILEPATH REDACTED >>>>'))
message("Register cluster")
registerDoParallel(cl)
message("Start foreach")
#Read in each .dta file in parallel - returns a list of data frames
top <- foreach(i=1:length(extractions), .packages = c('haven')) %dopar% {
  dta <- read_dta(extractions[i], encoding = 'latin1')
  return(dta)
}
message("Foreach finished")
message("Closing cluster")
stopCluster(cl)



message("rbindlist all extractions together")
topics <- rbindlist(top, fill=T, use.names=T)
rm(top)

#####################################################################
######################## PULL IN GEO CODEBOOKS ######################
#####################################################################

message("Retrieve geo codebook filepaths")
files <- list.files('<<<< FILEPATH REDACTED >>>>', pattern=".csv$", 
                    ignore.case = T, full.names = T)
message("Read geo codebooks into list")
geogs <- lapply(files, read_add_name_col)

message("Append geo codebooks together")
geo <- rbindlist(geogs, fill=T, use.names=T)
geo[is.na(admin_level), admin_level := "NA"]
geo <- geo[admin_level != "0", ] #drop anything matched to admin0
rm(geogs)

#Dedupe the geography codebook by geospatial_id, iso3, and nid
geo <- distinct(geo, nid, iso3, geospatial_id, .keep_all=T)

#coerce lat/longs to numeric
geo <- geo[, lat := as.numeric(lat)]
geo <- geo[, long := as.numeric(long)]

#####################################################################
######################## PREP DATA FOR MERGE ########################
#####################################################################

#Drop unnecessary geo codebook columns
geo_keep <- c("nid", "iso3", "geospatial_id", "point", "lat", "long", 
              "shapefile", "location_code", "survey_series")
geo_k <- geo[, geo_keep, with=F]

#####################################################################
############################### MERGE ###############################
#####################################################################

message("Merge ubCov outputs & geo codebooks together")
names(topics)[names(topics) == 'ihme_loc_id'] <- 'iso3'
geo_k$geospatial_id <- as.character(geo_k$geospatial_id)
topics$geospatial_id <- as.character(topics$geospatial_id)
all <- merge(geo_k, topics, by = c("nid", "iso3", "geospatial_id"), 
             all.x=F, all.y=T)

#####################################################################
######################### MAKE year_experiment COLUMN ###############
#####################################################################

message("Adding year_experiment column")
all[, start_year := year_start]
all[, year_dummy := start_year]
all[, year_experiment := year_dummy]

all[, year_experiment := round(mean(x=year_dummy, na.rm=T)), by=.(nid, iso3)]

all[(!is.na(int_year) & int_year <= year_start+5 & int_year >= year_start), 
    year_experiment := round(mean(int_year, na.rm=T)), by=c("nid", "iso3")]

#make point-level clusters annually representative of themselves
all[!is.na(lat) & !is.na(long) & !is.na(int_year) & int_year <= year_start+5 
    & int_year >= year_start, year_experiment := round(mean(x=int_year, na.rm=T)), 
    by=.(nid, iso3, lat, long)]

bad_year_nids <- unique(all[is.na(year_experiment), nid])
for (bad_nid in bad_year_nids){
  message(bad_nid)
  only_year <- unique(all[nid==bad_nid, year_experiment])
  only_year <- only_year[!is.na(only_year)]
  all[nid == bad_nid, year_experiment := only_year]
}
message("if a table longer than 0 rows appears here diagnose issues with year_experiment")
unique(all[is.na(year_experiment), .(nid, iso3)])
message("end of table")

rm(geo_k)
rm(geo)
rm(topics)


#####################################################################
######################### TOPIC-SPECIFIC CODE #######################
#####################################################################

message("WaSH-specific Fixes")
source('wash_specific_post_extract.R')

