# HEADER ------------------------------------------------------------------
# data_statistics.R
# Purpose: Produce descriptive statistics for data flow diagram - including
# where data is lost along the processing pathway.
#**************************************************************************

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

# Specific function loads (avoid namespace conflicts)
str_match <- stringr::str_match

### II: DEFINE VARIABLES, ETC. #####################################################

# Variables ------------------------------------------------------------------------
run_date     <- NULL # Specify here as "YYYY_MM_DD" if running on specific date
                     # Otherwise will use today if re-reading data
                     # Or the most recent existing if not
cores        <- 30                    

# Vaccines & vaccine titles --------------------------------------------------------

# Set up parameters for DPT
vaccine_list <- add_vaccine(prefix = "dpt", 
                            title = "DPT",
                            doses = 3,
                            age_min = 12,
                            age_max = 59)

vaccines <- names(vaccine_list)

graph_vaccines <- c("dpt3_cov")
final_model_vaccines <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout", "dpt1_3_rel_dropout")
graph_vaccine_titles <- c("DPT3")
log_dir <- "<<<< FILEPATH REDACTED >>>>"
dir.create(log_dir, showWarnings = F, recursive = T)


# Toggles --------------------------------------------------------------------------

drop_pre_2000   <- T     # Drop surveys and individuals before 2000

read_data       <- T     # Read in data and combine
merge_data      <- T     # Merge all the files together to vaccine_clean_list_file (below)
make_graphs     <- T     # Make data coverage graphs
do_pointpoly    <- T     # Run point/poly processing
make_csvs       <- T     # Generate CSVs

# Directories & files --------------------------------------------------------------

# directory for data prep code
script_dir <- paste0(repo, "vaccine/data_prep/")

# filepath holding extracted survey data in .csv format
in_dir <- "<<<< FILEPATH REDACTED >>>>"

# output directory
out_dir <- "<<<< FILEPATH REDACTED >>>>"

# Prep a run date-------------------------------------------------------------------
if (is.null(run_date) & read_data == F) {
  # make list of dates already  generated
  dates <- list.files(in_dir)
  dates <- sort(dates[grep(".*_.*_.*", dates)], decreasing = T)
  # pull most recent date folder & use that
  run_date <- dates[1]
} else if (is.null(run_date) & read_data == T) {
  run_date <- Sys.Date() #date-stamp 
  run_date <- gsub("-", "_", run_date)
}

in_dir <- paste0(in_dir, run_date, "/")
dir.create(in_dir, showWarnings = F)

# File for temporary, cleaned data set of processed vaccines - in same data as input data
# Will have the following added to it:  "..._[vaccine_prefix].rds"
vaccine_cleaned_file_prefix <- paste0(in_dir, "combined_vaccine_clean_")

# File for temporary, *resampled* data set of processed vaccines - in same data as input data
# Will have the following added to it:  "..._[vaccine_prefix].rds"
vaccine_resampled_file_prefix <- paste0(in_dir, "combined_vaccine_resampled_")

#############################################################################
## Data processing and cleaning (from 02_merge_vaccine.R)
#############################################################################

# Load all vaccine data
message("Loading vaccine data... \n")
vaccine_data <- readRDS(paste0(in_dir, run_date, '.Rda')) %>% 
  as.data.table

head(vaccine_data[!(is.na(shapefile)) & !(is.na(psu)),])

intro_dt <- readRDS(intro_file) %>% 
  as.data.table %>%
  subset(., select = c("ihme_loc_id", "cv_intro", "me_name")) %>%
  unique

#### Load in and apply exclusions saved in .csv file
message(paste0("Loading exclusion data... \n", "<<<< FILEPATH REDACTED >>>>/dpt_vaccine_exclusions.csv"))
exclusions <- fread('<<<< FILEPATH REDACTED >>>>/dpt_vaccine_exclusions.csv')
message("Applying Exclusions... \n")
exclusions <- exclusions[GEO_exclude==1,]

message("Excluding the following NIDs  \n")
for (nid in exclusions[,nid]){
  message(nid)
}

vaccine_data<- vaccine_data[!(nid %in% exclusions[,nid]),]

########################################################################################
### Exclude pre-2000 surveys
########################################################################################
vaccine_data <- vaccine_data[year_start >= 2000,]

########################################################################################
### Exclude non-Africa data
########################################################################################

gaul_list <- get_gaul_codes('africa') 
gaul_list <- gaul_list[gaul_list != 269] # Remove Yemen
gaul_list <- c(gaul_list, 40764) # Add SDN

# create a matching variable that excludes everything after "_" (subnationals)
vaccine_data[, match_country := ihme_loc_id]
vaccine_data[grep("_", vaccine_data$ihme_loc_id), match_country :=  str_match(ihme_loc_id, "([^_]+)_.*")[,2]]

t_lookup <- unique(vaccine_data$match_country) %>% as.character %>% as.data.table
names(t_lookup) <- "match_country" 
t_lookup[, gaul := gaul_convert(match_country)]

# Merge lookup table
vaccine_data <- merge(vaccine_data, t_lookup, all.x = T, by = "match_country")
vaccine_data[, match_country := NULL]

# Drop if no gaul code found
if (length(unique(vaccine_data[is.na(gaul)]$ihme_loc_id)) > 0) {
  message(paste0("No gaul codes found for ", 
                 paste(unique(vaccine_data[is.na(gaul)]$ihme_loc_id), collapse = ", "),
                 " - dropping..."))
}

vaccine_data <- vaccine_data[!is.na(gaul)]

# Drop if not in gaul list & remove gaul variable
keep_countries <- unique(vaccine_data[gaul %in% gaul_list]$ihme_loc_id)
drop_countries <- unique(vaccine_data[!(gaul %in% gaul_list)]$ihme_loc_id)
message(paste0("Dropping the following countries (not in gaul_list): ", paste(drop_countries, collapse = ", ")))
vaccine_data <- vaccine_data[gaul %in% gaul_list]
vaccine_data[, gaul := NULL]

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- nrow(vaccine_data)
tmp <- subset(vaccine_data, select = c("nid", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(vaccine_data$nid))

data_flow_list[["total_n_after_exclusion"]] <- data.table(title = "Total N after exclusions",
                                                          kept_n = total_n,
                                                          kept_psu = total_psu,
                                                          kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

### II: DATA CLEANING ##############################################################
# General data cleaning (applies to all vaccines)

message("cleaning data... \n")
# combine age categories
if ("child_age_month" %in% names(vaccine_data)) {
  vaccine_data[!is.na(child_age_month) & is.na(age_month), age_month := child_age_month]
  vaccine_data[, child_age_month := NULL]
}

if ("child_age_year" %in% names(vaccine_data)) {
  vaccine_data[!is.na(child_age_year) & is.na(age_year), age_year := child_age_year]
  vaccine_data[, child_age_year := NULL]
}

if ("child_sex_id" %in% names(vaccine_data)) {
  vaccine_data[, child_sex_id := NULL]
  vaccine_data[!is.na(child_sex_id) & is.na(sex_id), sex_id := child_sex_id]
}

# drop some variables
drop_vars <- c("smaller_site_unit", "line_id", "hh_id", 
               "admin_2_id", "admin_2_mapped", "admin_2",
               "admin_3_id", "admin_3_mapped", "admin_3",
               "admin_4_id", "admin_4_mapped", "admin_4",
               "admin_1_id", "admin_1_mapped", "file_path", 
               "nid.x", "iso3")

vaccine_data <- vaccine_data[,!(names(vaccine_data) %in% drop_vars), with = F]

# variable for year midpoint
# now switched to year_start
vaccine_data[, year_mid := (year_start + year_end) / 2]

#### calculate pweight from hhweight if pweight is missing
vaccine_data[is.na(pweight) & is.na(latitude) & !(is.na(hhweight)), pweight := hhweight]

# drop if not successfully geopositioned
vaccine_data[shapefile == "", shapefile := NA]
nrow_all <- nrow(vaccine_data)
geopos_pre <- vaccine_data[,.N, by=.(nid)]
names(geopos_pre) <- c("nid","nrow_pre_geoposition")

vaccine_data <- vaccine_data[!(is.na(lat) & is.na(shapefile))]
nrow_post_drop <- nrow(vaccine_data)
geopos_post <- vaccine_data[,.N, by=.(nid)]
names(geopos_post) <- c("nid","nrow_post_geoposition")

geopos_drop_log <- merge(geopos_pre, geopos_post, by="nid", all.x=T)
geopos_drop_log[is.na(nrow_post_geoposition), nrow_post_geoposition:=0]
message(paste0("Dropping ", nrow_all-nrow_post_drop,
               " rows not geopositioned (",
               round(((nrow_all-nrow_post_drop)/nrow_all)*100,2),
               "%)"))

geopos_drop_log[, percent_lost := ((nrow_pre_geoposition - nrow_post_geoposition)/nrow_pre_geoposition)*100]

message(paste0("Saving log file for rows dropped due to lack of geopositioning at: ", log_dir, "geopositioning_drop.csv"))
write.csv(geopos_drop_log, paste0(log_dir, "geopositioning_drop.csv"))

if(!"psu" %in% names(vaccine_data)){
  vaccine_data$psu <- vaccine_data$psu.y
}

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- nrow(vaccine_data)
tmp <- subset(vaccine_data, select = c("nid", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(vaccine_data$nid))

data_flow_list[["total_n_after_geopositioning"]] <- data.table(title = "After dropping non-geopositioned",
                                                               kept_n = total_n,
                                                               kept_psu = total_psu,
                                                               kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

vax <-"dpt"

### LOAD VACCINE-SPECIFIC PARAMETERS -------------------------------------------
vax_df <- copy(vaccine_data)

vax_prefix <- vaccine_list[[vax]][["prefix"]]
vax_title  <- vaccine_list[[vax]][["title"]]
vax_doses  <- vaccine_list[[vax]][["all_doses"]]
age_range  <- vaccine_list[[vax]][["age_range"]]
age_min <- as.numeric(min(age_range))
age_max <- as.numeric(max(age_range))

message("\n################################################")
message(paste0("Working on ", vax_title, "...\n"))

### ELIGIBILITY AND AGE COHORTS ------------------------------------------------

# Process age eligibility 
message("Processing age eligibility ...")

# Age information present?
vax_df[!is.na(age_month), age_info := 1]
vax_df[!is.na(age_year), age_info := 1]
vax_df[!is.na(age_categorical), age_info := 1]

total_ages <- vax_df[,.N,by="nid"]
names(total_ages) <- c("nid","total_rows")
no_age <- vax_df[!(is.na(age_info)),.N,by="nid"]
names(no_age) <- c("nid","no_age_rows")

no_age <- merge(no_age, total_ages, by="nid", all.y=T)
no_age[is.na(no_age_rows),no_age_rows:=0]
no_age[,percent_missing_ages := (total_rows - no_age_rows)/total_rows*100]

# Drop all without ages
message(paste0("Dropping ", nrow(vax_df[is.na(age_info)]), 
               " out of ", nrow(vax_df), " rows (",
               round(nrow(vax_df[is.na(age_info)])/nrow(vax_df) * 100, 1),
               "%) without age information."))

vax_df <- vax_df[!(is.na(age_info))]

message(paste0("Saving log file for rows dropped due to lack of age information at: ", log_dir, "age_data_drop.csv"))
write.csv(no_age, paste0(log_dir, "age_data_drop.csv"))

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- nrow(vax_df)
tmp <- subset(vax_df, select = c("nid", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(vax_df$nid))

data_flow_list[["total_n_after_age_info"]] <- data.table(title = "After dropping those without age info",
                                                         kept_n = total_n,
                                                         kept_psu = total_psu,
                                                         kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

# Drop those outside of age eligibility ranges
n_before <- nrow(vax_df)
age_range_pre <- vax_df[,.N, by=.(nid)]
names(age_range_pre)<-c("nid","pre_age_range_drop")
if (age_min == 12 & age_max == 59) {
  vax_df <- vax_df[(!is.na(age_month) & age_month >= 12) | (is.na(age_month) & age_year >= 1)]
  vax_df <- vax_df[(!is.na(age_month) & age_month <= 59) | (is.na(age_month) & age_year <= 4)]

  # Also drop a small number of observations with age_month < 60 but age_year > 4 
  vax_df <- vax_df[!(age_year > 4)]

} else {
  stop("Only the 12-59 month age range is currently supported.")
}

message(paste0("Dropping ", n_before - nrow(vax_df),
               " out of ", n_before, " rows (", 
               round((n_before - nrow(vax_df))/n_before * 100, 2), 
               "%) under ", age_min, " or over ", age_max, " months."))
age_range_post <- vax_df[,.N, by=.(nid)]
names(age_range_post)<-c("nid","post_age_range_drop")
age_range_post <- merge(age_range_post, age_range_pre, by="nid", all.y=T)
age_range_post[is.na(post_age_range_drop),post_age_range_drop := 0]
age_range_post[, percent_dropped_age_range := (pre_age_range_drop - post_age_range_drop)/pre_age_range_drop*100]

message(paste0("Saving log file for rows dropped from being out of the age range: ", log_dir, "age_range_drop.csv"))
write.csv(age_range_post, paste0(log_dir, "age_range_drop.csv"))

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- nrow(vax_df)
tmp <- subset(vax_df, select = c("nid", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(vax_df$nid))

data_flow_list[["total_n_after_age_eligibility"]] <- data.table(title = "After dropping those outside of age eligibility",
                                                                kept_n = total_n,
                                                                kept_psu = total_psu,
                                                                kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

### Generate birth cohorts ------------------------------------------------

message("\nGenerating birth cohorts ...")

if (age_min == 12 & age_max == 59) {
  vax_df[!is.na(age_month), birth_cohort := floor(age_month / 12)]
  vax_df[is.na(age_month), birth_cohort := as.numeric(age_year)]

  # Note: add one to year because estimating 12-23 month as target age range
  vax_df[, year_start_cohort := year_start - birth_cohort + 1] 

} else {
  stop("Only the 12-59 month age range is currently supported.")
}

### Exclude children in a cohort before 2000 -------------------------------------------------
if (drop_pre_2000 == T){
  # drop all children with a year_start_cohort before 2000
  message(paste0("\nDropping ", vax_df[year_start_cohort<2000,.N], " rows out of ",vax_df[,.N]," that start before 2000"))
  vax_df <- vax_df[year_start_cohort >= 2000]
}

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- nrow(vax_df)
tmp <- subset(vax_df, select = c("nid", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(vax_df$nid))

data_flow_list[["total_n_after_drop_pre2000"]] <- data.table(title = "After dropping pre-2000 cohorts",
                                                             kept_n = total_n,
                                                             kept_psu = total_psu,
                                                             kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

### Generate a categorical dummy variable for doses of vaccine given
dose_col <- paste0(vax_prefix, "_dose")

# Convert to numeric prior to vacc_ever check
vax_df[[dose_col]] <- as.numeric(vax_df[[dose_col]])
vax_df[[paste0(vax_prefix, "_ever")]] <- as.numeric(vax_df[[paste0(vax_prefix, "_ever")]])

# Ensure numeric
vax_df[[dose_col]] <- as.numeric(vax_df[[dose_col]])

for (i in 0:(max(vax_doses)-1)) {
  vax_df[get(dose_col) == i & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", i) := 1]
  vax_df[get(dose_col) != i & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", i) := 0]
}

vax_df[get(dose_col) >= max(vax_doses) & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", max(vax_doses)) := 1]
vax_df[get(dose_col) <  max(vax_doses) & !is.na(get(dose_col)), paste0(vax_prefix, "_dose_", max(vax_doses)) := 0]

# List for missing shapefiles
missing_shapefiles <- list()

### Collapse 
  
message(paste0("Collapsing ", vax_title))

# Only keep those with vaccine information
vax_df <- vax_df[!is.na(get(dose_col))]

# Keep only useful variables
keep_vars <- c("survey_name", "ihme_loc_id", "year_start", "nid", "geospatial_id", "point",
               "lat", "long", "location_code", "shapefile", "year_start_cohort", "shapefile1","loc_name1","loc_code1",
               "admin_level", "pweight", "strata","psu",paste0(vax_prefix, "_dose_", vax_doses))
vax_df <- subset(vax_df, select=keep_vars)

# Rename columns
rename_table <- data.table(rbind(c("survey_name", "source"), 
                                 c("ihme_loc_id", "country"), 
                                 c("year_start", "svy_year"), 
                                 c("nid", "svy_id"), 
                                 c("geospatial_id", "psu"), 
                                 c("lat", "latitude"),
                                 c("long", "longitude"), 
                                 c("year_start_cohort", "year"),
                                 c("shapefile1","admin1_shapefile"),
                                 c("loc_name1","admin1_name"),
                                 c("loc_code1","loc_code_admin1"),
                                 c("admin_level","admin_level"),
                                 c("point","point"),
                                 c("pweight","pweight"), 
                                 c("strata", "strata"),
                                 c("psu","psu")))

names(rename_table) <- c("old", "new")
setnames(vax_df, rename_table$old, rename_table$new)

# recode shapefile
vax_df[shapefile == "", shapefile := NA]
vax_df[admin1_shapefile == "", admin1_shapefile := NA]
vax_df[admin1_name == "", admin1_name := NA]
vax_df[loc_code_admin1 == "", loc_code_admin1 := NA]
vax_df[admin_level == "", admin_level := NA]
vax_df[,admin1_shapefile := as.character(admin1_shapefile)]

# Placeholder value to sum over
vax_df[,N := 1]

# Ensure column formats correct
if (class(vax_df$latitude) != "numeric") vax_df[, latitude := as.numeric(as.character(latitude))]
if (class(vax_df$longitude) != "numeric") vax_df[, longitude := as.numeric(as.character(longitude))]

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- nrow(vax_df)
tmp <- subset(vax_df, select = c("svy_id", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(vax_df$svy_id))

data_flow_list[["total_n_after_dropping_no_vax"]] <- data.table(title = "After dropping those without vaccine info",
                                                                kept_n = total_n,
                                                                kept_psu = total_psu,
                                                                kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

### COLLAPSE POLYGONS AND POINTS ###########################################

# Get shapefiles into character format (in case factor)
vax_df$shapefile <- as.character(vax_df$shapefile)

# Recode shapefile for points

# Figure out what shapefiles are missing and exclude a priori

shapefile_list <- unique(vax_df$shapefile[!is.na(vax_df$shapefile)])
shapefile_dir <- "<<<< FILEPATH REDACTED >>>>"
shapefiles_available <- list.files(shapefile_dir, ".shp$") %>% 
                          gsub(".shp", "", .)

no_shapefile <- shapefile_list[!(shapefile_list %in% shapefiles_available)]

if (length(no_shapefile > 0)) {
  message("The following shapefiles are missing from the shapefile directory. Associated data will be dropped:")
  message(paste(paste0(no_shapefile, ".shp"), collapse = " "))
  vax_df <- vax_df[!(shapefile %in% no_shapefile)]
  missing_shapefiles[[vax]] <- no_shapefile
}

# First, check if there are any polygons

if (nrow(vax_df[point == 0]) == 0) {

  #If so, skip ahead
  message("No polygon data found - moving to next cycle")
  df_pointpoly <- vax_df
  keep_vars <- c("svy_id", "source", "country", "point", "svy_year", 
                 "location_code", "shapefile", "year", 
                 "psu", "latitude", "longitude", paste0(vax_prefix, "_dose_", vax_doses), "N")
  df_pointpoly <- subset(df_pointpoly, select = names(df_pointpoly %in% keep_vars))

} else {

  # Collapse the polygon data
  vax_df[ ,cluster_id := NULL]
  vax_df[, strata := as.numeric(strata)]

  df_point <- vax_df[point == 1]
  df_poly <- vax_df[point == 0]

  pre_pweight_drop <- df_poly[,.N, .(svy_id)]
  names(pre_pweight_drop) <- c("nid","pre_pweight_drop")

  # Drop all rows without pweight
  warning(paste0('Dropping ', nrow(df_poly[is.na(pweight) & is.na(longitude),]), 
                 ' rows of ',nrow(df_poly),  ' due to missing pweights'))
  
  df_poly <- df_poly[!is.na(longitude) | is.na(longitude) & !is.na(pweight),]
  post_pweight_drop <- df_poly[,.N, .(svy_id)]
  names(post_pweight_drop) <- c("nid","post_pweight_drop")
  post_pweight_drop <- merge(pre_pweight_drop, post_pweight_drop,by="nid", all.x=T)
  post_pweight_drop[is.na(post_pweight_drop), post_pweight_drop:=0]

  post_pweight_drop[,percent_pweight_missing := (pre_pweight_drop - post_pweight_drop)/pre_pweight_drop*100]

  message(paste0("Saving log file for rows dropped from having no pweight: ", log_dir, "pweight_drop.csv"))
  write.csv(post_pweight_drop, paste0(log_dir, "pweight_drop.csv"))

  # create "lonely" column to indicate whether a polygon_year only has 1 individual
  by_vars <- c("svy_id", "source", "country", "point", 
               "svy_year", "location_code", "shapefile", "year")
  df_poly[is.na(latitude) & !is.na(shapefile), lonely := lapply(.SD, length), .SDcols="N", by=by_vars]
  
  warning(paste0('Dropping ', nrow(df_poly[lonely==1,]), 
                 " out of ", uniqueN(df_poly[!is.na(lonely),], by=c("location_code", "shapefile", "year")),
                 " polygon-years that contain just 1 individual"))
      
  ## For the purposes of these descriptive statistics, don't collapse but rather treat as if point data

  sum_cols <- c("N", paste0(vax_prefix, "_dose_", vax_doses))
  df_poly <- df_poly[, lapply(.SD, sum), 
                         by = .(svy_id, source, country, point, svy_year, 
                                location_code, shapefile, year, pweight, 
                                strata, psu, latitude, longitude, admin_level),
                         .SDcols = sum_cols]                   

  # Collapse point data
  sum_cols <- c("N", paste0(vax_prefix, "_dose_", vax_doses))
  df_point <- df_point[, lapply(.SD, sum), 
                         by = .(svy_id, source, country, point, svy_year, 
                                location_code, shapefile, year, pweight, 
                                strata, psu, latitude, longitude, admin_level),
                         .SDcols = sum_cols]

  # Harmonize names & combine
  drop_vars <- c("strata", "pweight")
  df_point <- subset(df_point, select = !(names(df_point) %in% drop_vars))
  df_poly  <- subset(df_poly, select = !(names(df_point) %in% drop_vars))

  df_pointpoly <- rbind(df_point, df_poly, fill = T)

}

# Add a weights column
df_pointpoly$weight <- 1

### Final renaming 
setnames(df_pointpoly, "svy_year", "original_year")

# For subnational locations, take just the country code
df_pointpoly[grepl("_", country),country:=stringr::str_match(country, "(.*)_")[,2]]

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- sum(df_pointpoly$N)
tmp <- subset(df_pointpoly, select = c("svy_id", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(df_pointpoly$svy_id))

data_flow_list[["total_n_after_pointpoly_process"]] <- data.table(title = "After pointpoly process",
                                                                  kept_n = total_n,
                                                                  kept_psu = total_psu,
                                                                  kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

#############################################################################
## Subsetting to Africa as per 04_pointpoly.R
#############################################################################
df_input <- copy(df_pointpoly)

## Crop to Africa -----------------------------------------
gaul_list <- get_gaul_codes('africa') 
gaul_list <- gaul_list[gaul_list != 269] # Remove Yemen
gaul_list <- c(gaul_list, 40764) # Add SDN

# create a matching variable that excludes everything after "_" (subnationals)
df_input[, match_country := country]
df_input[grep("_", df_input$country), match_country :=  str_match(country, "([^_]+)_.*")[,2]]

t_lookup <- unique(df_input$match_country) %>% as.character %>% as.data.table
names(t_lookup) <- "match_country" 
t_lookup[, gaul := gaul_convert(match_country)]

# Merge lookup table
df_input <- merge(df_input, t_lookup, all.x = T, by = "match_country")
df_input[, match_country := NULL]

# Drop if no gaul code found
if (length(unique(df_input[is.na(gaul)]$country)) > 0) {
  message(paste0("No gaul codes found for ", 
                 paste(unique(df_input[is.na(gaul)]$country), collapse = ", "),
                 " - dropping..."))
}
df_input <- df_input[!is.na(gaul)]

# Drop if not in gaul list & remove gaul variable
keep_countries <- unique(df_input[gaul %in% gaul_list]$country)
drop_countries <- unique(df_input[!(gaul %in% gaul_list)]$country)
message(paste0("Dropping the following countries (not in gaul_list): ", paste(drop_countries, collapse = ", ")))
df_input <- df_input[gaul %in% gaul_list]
df_input[, gaul := NULL]

###### PULL DATA NUMBERS HERE ##############################################################

total_n <- sum(df_input$N)
tmp <- subset(df_input, select = c("svy_id", "psu"))
total_psu <- nrow(unique(tmp))
total_nid <- length(unique(df_input$svy_id))

data_flow_list[["total_n_final"]] <- data.table(title = "Final (should match previous row)",
                                                                  kept_n = total_n,
                                                                  kept_psu = total_psu,
                                                                  kept_nid = total_nid)

rm(tmp, total_n, total_psu, total_nid)

###### END DATA NUMBERS PULLING ############################################################

#############################################################################
## Tabulate data here
#############################################################################

# Create separate columns for point & polys

df_input[point == 0, N_poly := N]
df_input[point == 1, N_poly := 0]
df_input[point == 0, N_point := 0]
df_input[point == 1, N_point := N]

# Generate unique svyid/cluster identifier
df_input[, unique_cluster_id := paste0(svy_id, "_", psu)]

df_input[point == 0, unique_cluster_id_poly := unique_cluster_id]
df_input[point == 1, unique_cluster_id_point := unique_cluster_id]
df_input[point == 1, unique_cluster_id_poly := NA]
df_input[point == 0, unique_cluster_id_point := NA]

### Collapse and count individuals and clusters
df_svy <- df_input[, .(N_clusters = uniqueN(unique_cluster_id, na.rm = T), 
                       N_clusters_poly = uniqueN(unique_cluster_id_poly, na.rm = T),
                       N_clusters_point = uniqueN(unique_cluster_id_point, na.rm = T),
                       N_individuals = sum(N),
                       N_individuals_poly = sum(N_poly),
                       N_individuals_point = sum(N_point)), 
                  by=c("svy_id","original_year", "source", "country")]

### Collapse and count individuals and clusters - by country
df_country <- df_input[, .(N_clusters = uniqueN(unique_cluster_id, na.rm = T), 
                           N_clusters_poly = uniqueN(unique_cluster_id_poly, na.rm = T),
                           N_clusters_point = uniqueN(unique_cluster_id_point, na.rm = T),
                           N_individuals = sum(N),
                           N_individuals_poly = sum(N_poly),
                           N_individuals_point = sum(N_point)), 
                    by=c("country")]

### Collapse and count individuals and clusters - by source type
df_source <- df_input[, .(N_clusters = uniqueN(unique_cluster_id, na.rm = T), 
                          N_clusters_poly = uniqueN(unique_cluster_id_poly, na.rm = T),
                          N_clusters_point = uniqueN(unique_cluster_id_point, na.rm = T),
                          N_individuals = sum(N),
                          N_individuals_poly = sum(N_poly),
                          N_individuals_point = sum(N_point)), 
                    by=c("source")]
            
##############################
# Code to generate survey counts for main text
##############################

##number of individuals
N_individuals       <- df_svy[,sum(N_individuals)]
N_individuals_poly  <- df_svy[, sum(N_individuals_poly)]
N_individuals_point <- df_svy[, sum(N_individuals_point)]

#number of countries
N_countries <- df_svy[, uniqueN(country)]

#unique sources
N_sources <- df_svy[,uniqueN(source)]

#unique surveys
N_surveys <- df_svy[,uniqueN(svy_id)]

#countries with data
countries_with_data <- unique(df_svy)

##number of unique psus 
N_clusters <- df_svy[, sum(N_clusters)]
N_clusters_poly <- df_svy[, sum(N_clusters_poly)]
N_clusters_point <- df_svy[, sum(N_clusters_point)]

##spatial resolution of polygon data
df_input[point == 1, admin_level := "point"]
df_input[, .(N_clusters = uniqueN(psu, na.rm = T),
             N_individuals = sum(N)), 
        by = c("admin_level", "point")]

###############################
# Write outputs
###############################
out_dir <- paste0(in_dir, "number_plugging/")
dir.create(out_dir, showWarnings = F)

# Write data flow table
data_flow <- rbindlist(data_flow_list)
fwrite(data_flow, file = paste0(out_dir, "data_flow.csv"))

# Write table of Ns
N_table <- as.data.table(rbind(c("Individuals", N_individuals, N_individuals_point, N_individuals_poly),
                               c("Clusters", N_clusters, N_clusters_point, N_clusters_poly),
                               c("Countries", N_countries, NA, NA),
                               c("Sources", N_sources, NA, NA), 
                               c("Surveys", N_surveys, NA, NA)))

names(N_table) <- c("Quantity", "Total", "Point", "Polygon")
fwrite(N_table, file = paste0(out_dir, "n_table.csv"))

# Write table of countries with data  
fwrite(df_country, file = paste0(out_dir, "all_countries.csv"))                             

# Write table of sources (and include %s)
tot_indv <- copy(N_individuals)
tot_clusters <- copy(N_clusters)         

df_source[, perc_individuals := (N_individuals / tot_indv)*100]
df_source[, perc_clusters := (N_clusters / tot_clusters)*100] 

fwrite(df_source, file = paste0(out_dir, "sources.csv"))

# Finally merge on the survey key to create a nicely formatted table
svy_key <- fread("<<<< FILEPATH REDACTED >>>>/survey_key.csv")
df_svy[, svy_id := as.numeric(svy_id)]                  
setnames(svy_key, "country", "country_long")  
df_svy_formatted <- merge(df_svy, svy_key, all.x = T, by.x = "svy_id", by.y = "nid")

cols <- c("country_long", "N_clusters_point", "N_clusters_poly", "N_individuals", "series", "years", "svy_id", "citation")
newcols <- c("Country", "GPS-located clusters", "Areally-located clusters", "Number of children included", "Series", "Year(s)", "NID for record in GHDx", "Citation")
df_svy_formatted <- subset(df_svy_formatted, select = cols)
setnames(df_svy_formatted, cols, newcols)
setorderv(df_svy_formatted, c("Country", "Year(s)"))                        

fwrite(df_svy_formatted, file = paste0(out_dir, "formatted_svy_table.csv"))

df_svy_missing <-  subset(df_svy_formatted, is.na(Citation))
df_svy_missing <- merge(df_svy_missing,
                        subset(df_svy, select = c("country", "original_year", "source", "svy_id")),
                        all.x = T, all.y = F,
                        by.x = "NID for record in GHDx", by.y = "svy_id")

fwrite(df_svy_missing, file = paste0(out_dir, "surveys_missing_descriptions.csv"))                      