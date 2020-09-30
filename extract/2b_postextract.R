#####################################################################
#Redacted
#####################################################################

#####################################################################
############################## SETUP ################################
#####################################################################
rm(list=ls())

#Define values
#Redacted
redownload <- F #update the codebook from google drive
cores <- 15
year_cutoff <- 2000 #only modelling >2000

#Setup
#Redacted

package_lib    <- file.path(h, '_code/_lib/pkg')
  .libPaths(package_lib)

#Load packages
pacman::p_load(haven, stringr, data.table, dplyr, magrittr, feather, fst, parallel, doParallel, googledrive, readxl)

#timestamp
today <- Sys.Date() %>% gsub("-", "_", .)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

## if running interactively, set arguments
if (interactive) {
  
  ## set arguments
  this_stage <- '1'
  core_repo <- "#Redacted"
  indicator_group <- 'cooking'
  
} else {
  
  ## otherwise, grab arguments from qsub
  ## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
  this_stage <- as.character(commandArgs()[4])
  core_repo <- as.character(commandArgs()[5])
  indicator_group <- as.character(commandArgs()[6])
  
} 

share.model.dir  <- file.path('#Redacted')
mod.dt <- file.path(share.model.dir, 'cooking_fuel_dirty.csv') %>% fread

#read info about stages, subset to current stage
stages <- file.path(j, '#Redacted/stage_master_list.csv') %>% 
  fread %>% 
  .[Stage==this_stage]
these_countries <- unique(stages$iso3)
#####################################################################
######################## DEFINE FUNCTIONS ###########################
#####################################################################
# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
package_list <- paste0(core_repo, '#Redacted/package_list.csv') %>% fread %>% t %>% c
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

source("#Redacted/geocodebook_functions.R")
#####################################################################
######################## BIND UBCOV EXTRACTS ########################
#####################################################################
#Redacted
codebook <- read_xlsx('hap.xlsx', sheet='codebook') %>% as.data.table

#subset codebook based on outliers, or files that are too large to handle here
#census is handled separately
codebook <- codebook[!(survey_name %like% 'IPUMS')]

#subset the codebook to ONLY the files that our extractors have worked on
codebook <- codebook[assigned %in% extractor_ids]

#if requested, subset to working indicators
if (indicator_group=='cooking') codebook <- codebook[!is.na(cooking_fuel)]

#subset the codebook to working stage and year
codebook <- codebook[ihme_loc_id %in% these_countries & year_end>=year_cutoff]
these_tails <- unique(codebook$nid) %>% 
  paste0(., '.csv', collapse='|') 

#Generate list of extraction filepaths, subset to region & to our extractions
extractions <- list.files(folder_in, full.names=T, pattern = ".csv", ignore.case=T, recursive = F) %>% 
  .[basename(.) %>% substr(., 1, 3) %in% these_countries] %>% #check that string starts with iso3 for stage
  .[str_ends(., these_tails)] #check that string ends in NID.csv with our extracted NIDs for stage

#subset to a country if testing
if (exists('test_country')) extractions <- extractions %>% .[. %like% paste0(test_country, '_')]

#read in all extractions
message('freading and appending all extractions')
topics <- mclapply(extractions, 
                   fread, 
                   integer64='character',
                   mc.cores = cores) %>% 
  rbindlist(fill=T, use.names=T)
#####################################################################
######################## VERIFY EXTRACTIONS##########################
#####################################################################
#Double check your compiled extractions against the codebook to make sure 
#that all variables are being extracted as expected.
#Return CSV with errors for each NID

#this step is only valid if the data is not being subset by test_country
if (!exists('test_country')) { 
  
  #also return a list of all the NIDs that are present in the codebook but not in the extracted topics
  #subset to make sure they are not stage3 or <2000
  message('writing csv of broken extractions')
  codebook.nids <- codebook[!(year_end < 2000 | 
                                ihme_loc_id %in% stages[Stage %in% c('2x', '3'), iso3] | 
                                is.na(cooking_fuel) |
                                is.na(hh_size)), nid] %>% unique
  broken_extractions <- codebook.nids[!(codebook.nids %in% unique(topics$nid))]
  write.csv(broken_extractions, paste0(folder_out, "/stage_",  this_stage, "_broken_extractions.csv"), na="", row.names=F)
  
  #make a vector of the expected variables
  var.list <- c('cooking_fuel', 'cooking_location', 'cooking_type', 'cooking_type_chimney',
                'heating_fuel', 'heating_type', 'lighting_fuel', 
                'electricity',
                'housing_roof', 'housing_wall', 'housing_floor',
                'housing_roof_num', 'housing_wall_num', 'housing_floor_num')
  
  #summarize the codebook and extracted topics inversely
  #for the topics, identify any columns in our list that are entirely missing
  topics.miss <- topics[, lapply(.SD, function(x) is.na(x) %>% all), .SDcols=var.list, by='nid']
  #for the codebooks, identify any columns in our list that are present
  #note that we subset by those NIDs that didnt fail to extract, because obviously those are missing
  codebook.prez <- codebook[!(nid %in% broken_extractions), lapply(.SD, function(x) !is.na(x)), .SDcols=var.list, by='nid']
  
  #bind together these two tables and then summarize using the product of each NID
  #this will have the effect of returning true for any NIDs that were present in the codebook but missing in the topics
  verification <- rbindlist(list(topics.miss,codebook.prez))[, lapply(.SD, prod, na.rm = TRUE), .SDcols=var.list, by='nid']
  verification[, failures := rowSums(.SD), .SDcols=var.list]
  #also merge on the notes column from the codebook, to help ID cases where there was missingess in the raw data
  verification <- merge(verification, codebook[, .(nid, notes)], by='nid', all.x=T)
  write.csv(verification[failures>0], paste0(folder_out, "/stage_",  this_stage, "_problem_extractions.csv"), na="", row.names=F)

}

#####################################################################
######################## PULL IN GEO CODEBOOKS ######################
#####################################################################
#use new geodata database 
geo_k <- get_geocodebooks(nids = unique(topics$nid))

#####################################################################
############################### MERGE ###############################
#####################################################################
message("Merge ubCov outputs & geo codebooks together")
geo_k[, geospatial_id := as.character(geospatial_id)]
topics[, geospatial_id := as.character(geospatial_id)]
all <- merge(geo_k, topics, by.x=c("nid", "iso3", "geospatial_id"), by.y=c("nid", "ihme_loc_id", "geospatial_id"), 
             all.x=F, all.y=T)
all[iso3 == "KOSOVO", iso3 := "SRB"] #GBD rolls Kosovo data into Serbia
message(nrow(topics)-nrow(all), ' rows lost in geocodebook merge!')
message(uniqueN(topics$nid)-uniqueN(all$nid), ' NIDs lost in geocodebook merge!')

#####################################################################
############################### MERGE DIAGNOSTIC ####################
#####################################################################
geo_nids <- as.data.table(geo_k$nid)
geo_nids <- unique(geo_nids)
topic_nids <- unique(topics$nid)
merged_correctly <- all[(!is.na(shapefile) & !is.na(location_code)) | (!is.na(lat) & !is.na(long)),]
merged_nids <- unique(merged_correctly$nid)

missing_nids <- topic_nids[(topic_nids %in% geo_nids) & !(topic_nids %in% merged_nids)]
if (length(missing_nids) > 0){
  message(paste("Writing csv of the", length(missing_nids), "surveys that are not properly merging"))
  merge_issues <- all[nid %in% missing_nids, .(nid, iso3, survey_name)] %>% distinct
  merge_issues <- merge(merge_issues, stages, by="iso3", all.x=T)
  if (!exists('test_country')) write.csv(merge_issues, paste0(folder_out, "/merge_issues.csv"), na="", row.names=F)
} else{
  message("All nids merged correctly. You are so thorough.")
  #Once R can handle unicode please add the clap emoji to this message.
}

#####################################################################
######################### MAKE year_experiment COLUMN ###############
#####################################################################
message("Adding year_experiment column")

all[, year_experiment := round(mean(year_start, na.rm=T)), by=.(nid, iso3)]

all[(!is.na(int_year) & int_year <= year_start+5 & int_year >= year_start), year_experiment := round(mean(int_year, na.rm=T)), by=c("nid", "iso3")]
print(nrow(all))
#make point-level clusters annually representative of themselves
all[!is.na(lat) & !is.na(long) & !is.na(int_year) & int_year <= year_start+5 & int_year >= year_start, year_experiment := round(mean(x=int_year, na.rm=T)), by=.(nid, iso3, lat, long)]
print(nrow(all))
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

#####################################################################
######################### TOPIC-SPECIFIC CODE #######################
#####################################################################
message("HAP-specific Fixes")
#replace missing hhweight with pweights
all[is.na(hhweight) & !is.na(pweight), hhweight := pweight]
print(nrow(all))
#drop useless vars
drop <- c("line_id", "sex_id", "age_year", "age_month", "pweight", "latitude", "longitude")
all <- all[, (drop):= NULL]
print(nrow(all))

#cleanup women module
message("dropping duplicate women in single household (so household sizes aren't duplicated)")
wn <- all[survey_module == "WN", ]
wn_key <- c("psu", "hh_id")
wn <- distinct(wn, psu, hh_id, .keep_all=T) 
all <- list(all[survey_module != "WN", ], wn) %>% 
  rbindlist(fill=T, use.names=T)

#create indicator for hh_size missingness
all[, missingHHsize := sum(is.na(hh_size)), by=nid]
all[, obs := .N, by=nid]
all[, pct_miss_hh_size := 100 * missingHHsize / obs]
all[, is_hh := pct_miss_hh_size > 0]
print(nrow(all))
#subset cases where all hh_sizes are present. Make sure each Row is a HH
has_hh_size_no_id <- all[!is.na(hh_size) & is.na(hh_id), ]
has_hh_size_id <- all[!is.na(hh_size) & !is.na(hh_id), ]
has_hh_size_id[, uq_id := paste(nid, psu, geospatial_id, hh_id, year_start, lat, long, shapefile, location_code, sep="_")] #includes space-time
has_hh_size_id[, prev_uq_id := paste(nid, psu, hh_id, sep="_")]
length(unique(has_hh_size_id$uq_id)) - length(unique(has_hh_size_id$prev_uq_id)) %>% 
message(paste("\n There are", ., "more unique households from including spacetime than excluding."))
hhhs <- distinct(has_hh_size_id, uq_id, .keep_all=T)
print(nrow(all))

#subset cases where all hh_sizes are missing and each row is not a HH. Set hh_size to 1
missing_hh_size <- all[is.na(hh_size) & survey_module != 'HH', ]
missing_hh_size[, hh_size := 1]

missing_hh_size_hh <- all[is.na(hh_size) & survey_module == 'HH', ]

all <- list(hhhs, has_hh_size_no_id, missing_hh_size, missing_hh_size_hh) %>% 
  rbindlist(fill=T, use.names=T)

#check to see if any NIDs were lost in the process by comparing to model dt
if (exists('test_country')) { 
  
  mod.dt[ihme_loc_id %like% test_country & survey_series!='IPUMS_CENSUS', unique(nid)] %>% 
    .[!(. %in% unique(all$nid))]
  
} else mod.dt[survey_series!='IPUMS_CENSUS' & ihme_loc_id %in% these_countries, unique(nid)] %>% .[!(. %in% unique(all$nid))]

message("Saving data")
if (!exists('test_country')) write.fst(all, path=paste0(folder_out, "/stage_",  this_stage, "_", today, ".fst"))

#####################################################################
######################### CLEAN UP & SAVE ###########################
#####################################################################

#Create & export a list of all surveys that have not yet been matched & added to the geo codebooks
message("Exporting a list of surveys that need to be geo matched")

if (!exists('test_country')) {
  
  all[!(nid %in% unique(geo_k$nid))] %>% 
    .[, .(nid, iso3, year_start, survey_name)] %>% 
    distinct %>% 
    merge(., stages, by='iso3', all.x=T) %>% 
    write.csv(., 
              file.path("#Redacted", topic,
                        paste0("stage_",  this_stage, "_new_geographies_to_match.csv")), 
              row.names=F, na="")
  print(nrow(all))

}