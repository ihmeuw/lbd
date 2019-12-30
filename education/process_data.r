################ process_data.R #######################################################
### Revised 5/2018 by Kate Wilson (kwilson7@uw.edu) 
### Postprocessing for LBD Education.
### Reads in and binds together all single year data, all level 7 binned data, and all
### level 8 binned data.
### As of 5/2018, drops all WHO_WHS surveys for data quality issues. 
### Should be run after parallel_collapse_binned in order for binned data to be updated.
######## Usually launched by submit_collapse before collapsing all indicators #########
########################################################################################


##################################### PREP ################################################

## Set repo locations and indicator group
core_repo <- paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_core/')
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", '/edu/')
cores <- 15

## Load libraries and  MBG project functions.
package_list <- readLines("<<<< FILEPATH REDACTED >>>>")

if (Sys.info()[1] == "Windows") stop("Don't run from Windows or bad things happen!")


# Source package imports function
source(paste0(core_repo, "/mbg_central/setup.R"))
load_R_packages(package_list)


######################## PULL IN GEO CODEBOOKS ############################################

#Read in geo codebook, add column with corresponding survey series
read_add_name_col <- function(file){
  #FOR GEOGRAPHY CODEBOOKS. READS THEM IN AND ADDS A COLUMN WITH THEIR CORRESPONDING SURVEY_SERIES
  message(file)
  rn <- gsub(".csv", "", file, ignore.case=T)
  spl <- unlist(strsplit(rn, "/"))
  svy <- spl[length(spl)]
  df <- fread(file, integer64="character", data.table=T)
  df <- as.data.table(df)
  df[, survey_series := svy]
  df <- lapply(df, as.character, stringsAsFactors = FALSE)
  return(df)
}

#Get all geog codebooks and package them together
message("Retrieve geo codebook filepaths")
files <- list.files(paste0("<<<< FILEPATH REDACTED >>>>"), pattern=".csv$", ignore.case = T, full.names = T)
message("Read geo codebooks into list")
geogs <- mclapply(files, read_add_name_col,mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores), mc.preschedule = F) # parallelizing this
message("Append geo codebooks together")
geo <- rbindlist(geogs, fill=T, use.names=T)

geo[is.na(admin_level), admin_level := "NA"] #set NA values for admin_level to "NA" as a string to keep the following line from dropping them because of bad R logic
geo <- geo[admin_level != "0", ] #drop anything matched to admin0
geo <- geo[survey_module != "GSEC1"] #drop this geomatch which is creating a m:m mismatch on different keys issue
rm(geogs)

#dedupe the geography codebook by geospatial_id, nid, and iso3
setkey(geo, nid, iso3, geospatial_id)
geo <- unique(geo, use.key=T)

# change types so that NA checks later work
geo[, location_code := as.numeric(location_code)]
geo[, lat := as.numeric(lat)]
geo[, long := as.numeric(long)]
geo[, nid := as.numeric(nid)]


##### Read in SINGLE YEAR DATA ###################################

source(paste0(indic_repo, '/education/edu_collapse_functions.R'))
list.of.packages <- c("data.table","gtools","ggplot2","haven","parallel","SDMTools", "readstata13")
lapply(list.of.packages, require, character.only = TRUE)
single_year_dt <- collapse_single_year(cores = cores,
                                       in_dir = "<<<< FILEPATH REDACTED >>>>")
single_year_dt[, nid := as.numeric(nid)]
single_year_dt <- single_year_dt[year >1997]

#print nids of surveys that didn't merge
all_sy <- merge(single_year_dt, geo, by = c('geospatial_id', 'nid', 'iso3'), all.x=T)

# drop missing data
all_sy <- all_sy[!is.na(edu_min), ]

## Some diagnostics on failed geography merges.
## It seems the creation of level_7 or level_8 went awry in some. For example, some IPUMS seem to require the merging of two raw variables
## to create the geospatial_id in the geography codebook but only have one raw variable as "level_7".
failed_nids <- length(unique(single_year_dt[!(nid %in% all_sy[, nid]), nid]))
message(paste0(failed_nids, '/', length(unique(single_year_dt[, nid])), ' of extracted NIDs failed to merge.'))
if (failed_nids > 0) {
  for(t in unique(single_year_dt$survey_name)) {
    message(paste0('--- Failed ', t, ' nids: ', paste(unique(single_year_dt[survey_name==t & !(nid %in% all_sy[, nid]), nid]), collapse=" ")))
  }
}

############### SAVE SINGLE YEAR DATA to be read in by various collapse codes ###################
# drop WHS
all_sy <- all_sy[survey_name != "WHS" & survey_series != "WHO_WHS",]

all_sy <- all_sy[year > 1996]


saveRDS(all_sy, paste0("<<<< FILEPATH REDACTED >>>>", format(Sys.Date(), "%Y_%m_%d"), ".rds"))


####################### Read in BINNED DATA #############################################

## Collapse binned point data.
binned_point_dir <- "<<<< FILEPATH REDACTED >>>>"
binned_point_files <- list.files(binned_point_dir)

collapse_binned_points <- function(x) {
  cluster_data <- fread(paste0(binned_point_dir, '/', x))
  print(paste0("Working on file: ", x))
  if(!"geospatial_id" %in% colnames(cluster_data)) cluster_data[, geospatial_id := ihme_loc_id]
  if(!"iso3" %in% colnames(cluster_data)) cluster_data[, iso3 := level_3]
  # Keep only necessary columns
  cols <- c("age", "sex", "year", "nid", "survey_name", "sample_size", "bin", "edu_yrs", "prop", "rake", "total", "total_ss", "proportion", "orig_proportion_se", "prop_se", "proportion_se", "mean", "sd", "mean_se", "geospatial_id", "iso3", "kish_N")
  subset_data <- cluster_data[, ..cols]
  return(subset_data)
}
binned_level_8 <- rbindlist(mclapply(binned_point_files, collapse_binned_points, mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores), mc.preschedule = F))

## Collapse binned polygon data.
binned_poly_dir <- "<<<< FILEPATH REDACTED >>>>"
binned_poly_files <- list.files(binned_poly_dir)
collapse_binned_polys <- function(x) {
  cluster_data <- fread(paste0(binned_poly_dir, '/', x))
  if(!"geospatial_id" %in% colnames(cluster_data)) cluster_data[, geospatial_id := ihme_loc_id]
  if(!"iso3" %in% colnames(cluster_data)) cluster_data[, iso3 := level_3]
  # Keep only necessary columns
  cols <- c("age", "sex", "year", "nid", "survey_name", "sample_size", "bin", "edu_yrs", "prop", "rake", "total", "total_ss", "proportion", "orig_proportion_se", "prop_se", "proportion_se", "mean", "sd", "mean_se", "geospatial_id", "iso3", "kish_N")
  subset_data <- cluster_data[, ..cols]
  return(subset_data)
  
}
binned_level_7 <- rbindlist(mclapply(binned_poly_files, collapse_binned_polys,mc.cores=ifelse(Sys.info()[1]=="Windows", 1, cores), mc.preschedule = F))

## Combine crosswalked extracts, merge geography codebooks, and collapse to points or polygons.
all_binned <- rbind(binned_level_7, binned_level_8)
all_binned[, nid := as.numeric(nid)]
all_binned[, year := as.numeric(year)]
all_binned <- all_binned[year > 1997]
all_by <- merge(all_binned, geo, by = c('geospatial_id', 'nid', 'iso3'))

## Some diagnostics on failed geography merges.
failed_nids <- length(unique(all_binned[!(nid %in% all_by[, nid]), nid]))
message(paste0(failed_nids, '/', length(unique(all_binned[, nid])), ' of extracted NIDs failed to merge.'))
if (failed_nids > 0) {
  for(t in unique(all_binned$survey_name)) {
    message(paste0('--- Failed ', t, ' nids: ', paste(unique(all_binned[survey_name==t & !(nid %in% all_by[, nid]), nid]), collapse=" ")))
  }
}

############### SAVE BINNED DATA to be read in by various collapse codes ###################
# drop WHS
all_by <- all_by[survey_name != "WHS" & survey_series != "WHO_WHS",]

all_by <- all_by[year > 1996]


saveRDS(all_by, paste0("<<<< FILEPATH REDACTED >>>>", format(Sys.Date(), "%Y_%m_%d"), ".rds"))


