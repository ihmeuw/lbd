#####################################################################
############################## SETUP ################################
#####################################################################
rm(list=ls())

#Define values
topic <- "lri"
priority <- '1' #designate '1', '2a', '2b', or '3' to output missing geogs from a given priority level
nid_vec <- c() #list the NIDs of surveys (if any) to merge using admin_1 column
cluster <- TRUE #cluster true/false
cores <- 5

#Setup
root <- '<<<< FILEPATH REDACTED >>>>' #filepath root
folder_in <- '<<<< FILEPATH REDACTED >>>>' #where your extractions are stored
folder_out <- '<<<< FILEPATH REDACTED >>>>' #where you want to save the big csv of all your extractions together
options(warn=-1)
module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)

#load packages and core functions
package_list <- c('haven', 'stringr', 'plyr', 'data.table', 'magrittr', 'parallel', 'doParallel', 'feather')
source("'<<<< FILEPATH REDACTED >>>>'")
mbg_setup(package_list = package_list, repos="")


#####################################################################
######################## DEFINE FUNCTIONS ###########################
#####################################################################
#Read in geo codebook, add column with corresponding survey series
read_add_name_col <- function(file){
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- strsplit(rn, "/") %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file, integer64="character", data.table=T)
  df[, survey_series := svy]
  df <- lapply(df, as.character)
  return(df)
}

#Convert all values of a data frame to characters
all_to_char_df <- function(df){
  df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
  return(df)
}


#####################################################################
######################## BIND UBCOV EXTRACTS ########################
#####################################################################
#Generate list of extraction filepaths
extractions <- list.files(folder_in, full.names = T, pattern = ".dta", recursive = F)

#rbindlist all ubcov extracts together
if(cluster == TRUE) {
  message("Make cluster")
  cl <- makeCluster(cores)
  clusterEvalQ(cl, .libPaths('')) #path to library
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
} else if(cluster == FALSE) {
  top <- foreach(i=1:length(extractions)) %do% {
    message(paste0("Reading in: ", extractions[i]))
    dta <- read_dta(extractions[i])
    return(dta)
  }
}

message("rbindlist all extractions together")
topics <- rbindlist(top, fill=T, use.names=T)
rm(top)
topics <- data.table(topics)
#####################################################################
######################## PULL IN GEO CODEBOOKS ######################
#####################################################################
#Read table with country priority stages
iso <- fread('', stringsAsFactors=F)

#Get all geog codebooks and package them together
message("Retrieve geo codebook filepaths")
files <- list.files('', pattern=".csv$", ignore.case = T, full.names = T) #grab filepath to all geography information codebooks

message("Read geo codebooks into list")
geogs <- lapply(files, read_add_name_col)

message("Bind geo codebooks together")
geo <- rbindlist(geogs, fill=T, use.names=T)
rm(geogs)

#Dedupe the geography codebook by geospatial_id and nid
setkey(geo, nid, geospatial_id)
geo <- unique(geo, use.key=T)


#####################################################################
######################## PREP DATA FOR MERGE ########################
#####################################################################
#Reconcile ubCov & geo codebook data types
message("Force NID data types to match")
if (class(topics$nid) == "numeric"){
  geo[, nid := as.numeric(nid)]
} else if (class(topics$nid) == "character"){
  geo[, nid := as.character(nid)]
} else{
  message("update code to accomodate topics nid as")
  message(class(topics$nid))
}
topics$geospatial_id <- as.character(topics$geospatial_id)
geo$geospatial_id <- as.character(geo$geospatial_id)

#Drop unnecessary geo codebook columns
geo_keep <- c("nid", "iso3", "geospatial_id", "point", "lat", "long", "shapefile", "location_code", "survey_series")
geo_k <- geo[, geo_keep, with=F]

#Drop geospatial_id special characters
#redone in data table for speed
topics[, geospatial_id := gsub("[^[:alnum:] | _ ]", "", geospatial_id)]

geo_k[, geospatial_id := gsub("[^[:alnum:] | _ ]", "", geospatial_id)]

#####################################################################
############################### MERGE ###############################
#####################################################################

message("Merge ubCov outputs & geo codebooks together")
all <- merge(geo_k, topics, by.x=c("nid", "iso3", "geospatial_id"), by.y=c("nid", "ihme_loc_id", "geospatial_id"), all.x=F, all.y=T)


##################### OPTIONAL - SPECIFY SURVEYS TO MERGE ON admin_1 NAME RATHER THAN GEOSPATIAL_ID #######################

if(length(nid_vec) > 0) {
  message("SPECIFIC RECODES")
  match_on_admin1_string <- function(nid_vec){
    nid_vec <- as.character(nid_vec)
    for (nid_i in nid_vec){
      if(nid_i %in% all$nid){
        message(paste("fixing", nid_i))
        nid_geo <- subset(geo, nid == nid_i)
        nid_dta <- subset(topics, nid == nid_i)
        nid_geo[, "svy_area1"] <- as.character(nid_geo[, "svy_area1"])
        nid_dta[, "admin_1"] <- str_trim(nid_dta[, "admin_1"])
        nid_all <- merge(nid_geo, nid_dta, by.x = c("nid", "svy_area1"), by.y=c("nid", "admin_1"), all.x=F, all.y = T)
        all <<- all[all$nid != nid_i, ] #clear from all
        all <<- rbind.fill(all, nid_all) #rbind to all
      }
    }
  }
  match_on_admin1_string(nid_vec)
}

##########################################################################################################################

#####################################################################
######################### CLEAN UP & SAVE ###########################
#####################################################################

#Create & export a list of all surveys that have not yet been matched & added to the geo codebooks
message("Exporting a list of surveys that need to be geo matched")
gnid <- unique(geo$nid)
fix <- subset(all, !(all$nid %in% gnid))
fix <- as.data.frame(fix)
fix_collapse <- unique(fix[c("nid", "iso3", "year_start", "survey_name")])
write.csv(fix_collapse, paste0(folder_out, "/geographies_to_match_", module_date, ".csv"), row.names=F)

#Same as above, only for countries of a designated priority level
if(!is.null(priority)) {
  message(paste0("writing csv of unmatched extractions from priority ", priority, " countries"))
  stage = iso[Stage == priority,]$iso3
  fix_collapse$iso3 <- strsplit(fix_collapse$iso3, 1, 3)
  fix_collapse$iso3 = as.character(fix_collapse$iso3)
  fix_stage <- subset(fix_collapse, iso3 %in% stage)
  write.csv(fix_stage, file=paste0(folder_out, "/priority", priority, "_geography_matching_to_do_", module_date, ".csv"), row.names = F)
}

#Diagnostic for completeness of the geography match for all surveys
message("Generate diagnostic table of geo matching completeness")
all = all[, n_row := .N, by=nid]
all = all[, n_matched := length(na.omit(point)), by=nid]
thin = all[, c("nid", "iso3", "year_start", "survey_series", "n_row", "n_matched")]
short = unique(thin)
short = short[, pct_matched := (n_matched/n_row)*100]
#Compare the number of unique geographies in the geo cb to the number in the extract
short$n_in_merge = NA
short$n_in_geocb = NA
for(i in 1:nrow(short)) {
  n = short$nid[i]
  short$n_in_merge[i] = length(unique(all[nid == n,]$geospatial_id))
  short$n_in_geocb[i] = length(unique(geo_k[nid == n,]$geospatial_id))
}
short = short[, pct_present := (n_in_merge/n_in_geocb)*100]
short = short[pct_present == Inf, pct_present := 0]
#saving geomatching merge diagnostics
write.csv(short, '', row.names = F)


#Save
message("Saving as output")
write_feather(all, '<<<< FILEPATH REDACTED >>>>')