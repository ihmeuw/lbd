#####################################################################
## POST UBCOV EXTRACTION DATA CLEANING FOR GEOSPATIAL DATA EXTRACTIONS & GEOGRAPHY MATCHING
#####################################################################

## Setup, load & clean data -------------------------------------------------------------------------
rm(list=ls())

#Define values
topic <- 'ort'
nid_vec <- c() #list the NIDs of surveys (if any) to merge using admin_1 column
cluster <- TRUE #cluster true/false
sing <- T # running on RStudio singularity image? T/F
repo <- '<<<< FILEPATH REDACTED >>>>'
cores <- 20
module_date <- Sys.Date()
module_date <- gsub('-', '_', module_date)

# directories
l <- '<<<< FILEPATH REDACTED >>>>'
folder_in <- paste0(l, '<<<< FILEPATH REDACTED >>>>') #where your extractions are stored
folder_out <- paste0(l, '<<<< FILEPATH REDACTED >>>>') #where you want to save the big csv of all your extractions together

# load libraries
package_list <- readLines('<<<< FILEPATH REDACTED >>>>/package_list.csv')

# load functions
message('Loading in required R packages and MBG functions')
source(paste0(repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = repo)
library('feather')
library('stringr')


## Define functions -------------------------------------------------------------------------

# read in geo codebook, add column with corresponding survey series
read_add_name_col <- function(file){
  rn <- gsub('.csv', '', file, ignore.case=T)
  spl <- strsplit(rn, '/') %>% unlist()
  svy <- spl[length(spl)]
  df <- fread(file, integer64='character', data.table=T)
  df[, survey_series := svy]
  df <- lapply(df, as.character)
  return(df)
}

# convert all values of a data frame to characters
all_to_char_df <- function(df){
  df <- data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
  return(df)
}


## Bind UbCov extracts -------------------------------------------------------------------------

# generate list of extraction filepaths
extractions <- list.files(folder_in, full.names = T, pattern = '.dta', recursive = F)

# rbindlist all ubcov extracts together
if(cluster == TRUE) {
  message('Make cluster')
  cl <- makeCluster(cores)
  message('Register cluster')
  registerDoParallel(cl)
  message('Start foreach')
    #Read in each .dta file in parallel - returns a list of data frames
    top <- foreach(i=1:length(extractions), .packages = c('haven')) %dopar% {
      dta <- read_dta(extractions[i], encoding = 'latin1')
      return(dta)
    }
  message('Foreach finished')
  message('Closing cluster')
  stopCluster(cl)
} else if(cluster == FALSE) {
  top <- foreach(i=1:length(extractions)) %do% {
    message(paste0('Reading in: ', extractions[i]))
    dta <- read_dta(extractions[i])
    return(dta)
  }
}

message('rbindlist all extractions together')
topics <- rbindlist(top, fill=T, use.names=T)
rm(top)

##Save raw data file, if desired
save(topics, file=paste0(folder_out, '/ort_no_geogs_', module_date, '.Rdata')) 

# Check number of rows in data table
message(paste0('There are ', nrow(topics), ' rows in the rbound data files.'))


## Pull in geo-codebooks ------------------------------------------------------------------------

# read table with country priority stages
iso <- fread('<<<< FILEPATH REDACTED >>>>/stage_master_list.csv', stringsAsFactors=F) # <- this is out of date
iso <- iso[, c('location_name', 'loc_id', 'sdi_2017', 'Stage', 'iso3')]

# get all geog codebooks and package them together
message('Retrieve geo codebook filepaths')
files <- list.files('<<<< FILEPATH REDACTED >>>>', pattern='.csv$', ignore.case = T, full.names = T)

message('Read geo codebooks into list')
geogs <- lapply(files, read_add_name_col)

message('Bind geo codebooks together')
geo <- rbindlist(geogs, fill=T, use.names=T)
rm(geogs)

# deduplicate the geography codebook by geospatial_id and nid
setkey(geo, nid, geospatial_id)
geo <- unique(geo, use.key=T)


## Prep data for merge -------------------------------------------------------------------------

# reconcile ubCov & geo codebook data types
message('Force NID data types to match')
if (class(topics$nid) == 'numeric'){
  geo[, nid := as.numeric(nid)]
} else if (class(topics$nid) == 'character'){
  geo[, nid := as.character(nid)]
} else{
  message('update code to accomodate topics nid as')
  message(class(topics$nid))
}
topics$geospatial_id <- as.character(topics$geospatial_id)
geo$geospatial_id <- as.character(geo$geospatial_id)

# drop unnecessary geo codebook columns
geo_keep <- c('nid', 'iso3', 'geospatial_id', 'point', 'lat', 'long', 'shapefile', 'location_code', 'survey_series', 'user')
geo_k <- geo[, geo_keep, with=F]

# drop geospatial_id special characters
topics[, geospatial_id := gsub('[^[:alnum:] | _ ]', '', geospatial_id)]
geo_k[, geospatial_id := gsub('[^[:alnum:] | _ ]', '', geospatial_id)]


## Merge -------------------------------------------------------------------------

# record number of observations prior to merge
topics[, observations_pre_merge :=  .N, by = c('nid', 'ihme_loc_id', 'year_start')]

# merge
message('Merge ubCov outputs & geo codebooks together')
all <- merge(geo_k, topics, by.x=c('nid', 'iso3', 'geospatial_id'), by.y=c('nid', 'ihme_loc_id', 'geospatial_id'), all.x=F, all.y=T)

# Check number of rows in data table
message(paste0('There are ', nrow(all), ' rows in the merged data file.'))

# record number of observations by shapefile and after merge
all[, observations_post_merge :=  .N, by = c('nid', 'iso3', 'year_start')]
all[, observations_pct_change := (observations_post_merge - observations_pre_merge)/observations_pre_merge]
all[, observations_by_shapefile := .N, by = c('nid', 'iso3', 'year_start', 'point', 'shapefile')]
all[, pct_by_shapefile := observations_by_shapefile/observations_pre_merge]


## Identify surveys that have been geomatched twice or to more than one point/shapefile ------------------------------------------------------

# save record of duplicated studies
dup_summary <- unique(all[observations_pct_change > 0.05, c('nid', 'iso3', 'year_start', 'shapefile', 'point', 'observations_pct_change', 'pct_by_shapefile', 'observations_pre_merge', 'observations_post_merge', 'observations_by_shapefile', 'user')])
write.csv(dup_summary, paste0(folder_out, '/twice_matched_nids_', module_date, '.csv'))

# clean up data file
all[, c('observations_pct_change', 'pct_by_shapefile', 'observations_pre_merge', 'observations_post_merge', 'observations_by_shapefile', 'user') := NULL]


## Clean up and save -------------------------------------------------------------------------

# create & export a list of all surveys that have not yet been matched & added to the geo codebooks
message('Exporting a list of surveys that need to be geo matched')
gnid <- unique(geo$nid)
fix <- subset(all, !(all$nid %in% gnid))
fix <- as.data.frame(fix)
fix_collapse <- unique(fix[c('nid', 'iso3', 'year_start', 'survey_name')])
write.csv(fix_collapse, paste0(folder_out, '/geographies_to_match_', module_date, '.csv'), row.names=F)

# diagnostic for completeness of the geography match for all surveys
message('Generate diagnostic table of geo matching completeness')
all = all[, n_row := .N, by=nid]
all = all[, n_matched := length(na.omit(point)), by=nid]
thin = all[, c('nid', 'iso3', 'year_start', 'survey_series', 'n_row', 'n_matched')]
short = unique(thin)
short = short[, pct_matched := (n_matched/n_row)*100]

# compare the number of unique geographies in the geo cb to the number in the extract
short$n_in_merge = NA
short$n_in_geocb = NA
for(i in 1:nrow(short)) {
  n = short$nid[i]
  short$n_in_merge[i] = length(unique(all[nid == n,]$geospatial_id))
  short$n_in_geocb[i] = length(unique(geo_k[nid == n,]$geospatial_id))
}
short = short[, pct_present := (n_in_merge/n_in_geocb)*100]
short = short[pct_present == Inf, pct_present := 0]
write.csv(short, paste0(folder_out, '/merge_diagnostic_', module_date, '.csv'), row.names = F)

# save
message('Saving as .Rdata')
save(all, file=paste0(folder_out, '/', module_date, '.Rdata'))
message('Saving as .csv')
write.csv(all, file=paste0(folder_out, '/', module_date, '.csv'))
message('Saving as feather')
write_feather(all, paste0(folder_out, '/', module_date, '.feather'))


