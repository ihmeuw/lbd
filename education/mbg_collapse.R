# Collapse education extracts to MBG input dataset for given indicator passed from submission script.
# Nick Graetz --> Lauren Woyczynski
# Notes:
#   - Right now, this script requires running submit_collapse_binned.R first and waiting for all jobs to finish.

###########################################################################################
########### Set up. ####################################################################### 
###########################################################################################

if (Sys.info()[1] == "Windows") stop("Don't run from Windows or bad things happen!")

## I'm running this script on the RStudio singularity, please adjust cores based on environment.
## Could change this to be an argument passed from submit_mbg_collapse.R
rm(list=ls())
cores <- 15
indicator_group <- 'education'

## Set repo locations and indicator group
core_repo <- paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_core/')
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", '/edu/')
commondir      <- sprintf("<<<< FILEPATH REDACTED >>>>")

## Load libraries and  MBG project functions.
package_list <- c(readLines("<<<< FILEPATH REDACTED >>>>"), 'ggrepel')

# Source package imports function
source(paste0(core_repo, "/mbg_central/setup.R"))
mbg_setup(package_list = package_list, repos = c(core_repo))

## Setting for compiling
indicator_name <- as.character(commandArgs()[4])
print(as.character(commandArgs()[4]))


## if running interactively make sure to source this file to read in and save the most recent data
# source(paste0(indic_repo, "/education/process_data.r"))

## Get sex and age groups for collapsing from indicator_name.
target_sex <- ifelse(grepl('female', indicator_name), 2, 1)
age_min <- ifelse(grepl('15_49', indicator_name), 15, 20)
age_max <- ifelse(grepl('15_49', indicator_name), 49, 24)

# function to find most recent dated file in folder. used to read in processed pre-collapse data
most_recent_date <- function(dir, date_format = "%Y_%m_%d", out_format = "%Y_%m_%d", file_pattern = NULL) {
  date_pattern <- gsub("%y|%m|%d", "[[:digit:]]{2}", gsub("%Y", "[[:digit:]]{4}", date_format))
  dates <- dir(dir, pattern = date_pattern)
  if (!is.null(file_pattern)) dates <- grep(file_pattern, dates, value = T)
  dates <- gsub(paste0("(.*)(", date_pattern, ")(.*)"), "\\2", dates)
  dates <- as.Date(dates, date_format)
  format(max(dates), out_format)
}

current_date <- format(Sys.Date(), "%Y_%m_%d") # used in determining the folder in /share/geospatial/mbg/education/mapped_datasets to save things into
mapped_datasets_dir <- paste0("<<<< FILEPATH REDACTED >>>>", current_date, "/")
if (!dir.exists(mapped_datasets_dir)) dir.create(mapped_datasets_dir, showWarnings = TRUE)

###########################################################################################
###### 1. Collapse single-year data to MBG geographies. ###################################
###########################################################################################

all_sy_date <- most_recent_date(dir="<<<< FILEPATH REDACTED >>>>", file_pattern="all_sy_*")
all_sy <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>", all_sy_date, ".rds"))
setnames(all_sy, 'edu_min', 'edu_yrs')

## Subset to target age/sex group.
all_sy <- all_sy[age >= age_min & age <= age_max & sex == target_sex, ]

## Drop any surveys that don't contain complete age range
# subset to ages age_start-age_max (and only surveys with this full range)
orig_nids <- unique(all_sy$nid)
drop_sy <- all_sy[, as.list(range(age)), nid][V1 != age_min | V2 != age_max,]

all_sy <- all_sy[!nid %in% drop_sy$nid,]
missing_nid <- setdiff(orig_nids, all_sy$nid)
if (length(missing_nid)>0) {
  print(paste("NIDs have been dropped for incomplete age-range:", paste(missing_nid, collapse=", ")))
}


############ Create edu_yrs for proportional indicators and subset data accordingly
## If this is a proportion indicator rather than a mean indicator, convert edu_yrs to a 0/1 and take the weighted.mean.
## edu_no_primary_prop = (N_no_primary) / (N_total)
## edu_primary_prop = (N_primary) / (N_total - N_no_primary)
## edu_secondary_prop = residual from MBG models (implicitly modeled via edu_primary_prop, because this is the proportion with primary
##                      of those who are above no primary)
if(grepl('prop',indicator_name)) {
  
  if(grepl('edu_zero_prop',indicator_name)) {
    all_sy[, edu_yrs := ifelse(edu_yrs == 0, 1, 0)]
    
  }
  if(grepl('edu_no_primary_prop',indicator_name)) {
    all_sy <- all_sy[edu_yrs >= 1, ] ## Remove those who have zero education from denominator.
    all_sy[, edu_yrs := ifelse(edu_yrs < 6, 1, 0)]
    
  }
  if(grepl('edu_primary_prop',indicator_name)) {
    all_sy <- all_sy[edu_yrs >= 6, ] ## Remove those who have less than primary from denominator.
    all_sy[, edu_yrs := ifelse(edu_yrs < 12, 1, 0)] ## OF THOSE WHO HAVE AT LEAST PRIMARY, what proportion have less than secondary?
  }
}

## Collapse all single-year data to either point (lat/long) or polygon (location_code/shapefile).

## Collapse point data. Here use "num_persons" to weight, the actual sample size.
point_sy_micro <- all_sy[!is.na(lat) & !is.na(long) & lat!="" & long!="", ]
# Drop any observations missing education or weight
point_sy_micro <- point_sy_micro[!is.na(edu_yrs) & !is.na(num_persons),]
point_sy <- point_sy_micro[,.(N=sum(num_persons),
                              kish_N = sum(kish_N, na.rm=T),
                              sum_sample_weight=sum(num_persons,na.rm=T),
                              mean=weighted.mean(x=edu_yrs,w=num_persons)),
                           by=.(year,iso3,survey_name,nid,lat,long,survey_series)]


## Collapse polygon data. Here use "count" to weight, the survey-weighted sample size.
poly_sy_micro <- all_sy[!is.na(location_code) & !is.na(shapefile) & location_code!="" & shapefile!="" & (is.na(lat) | lat == ""), ]
## Drop any observations missing education or weight
poly_sy_micro <- poly_sy_micro[!is.na(edu_yrs) & !is.na(count),]
poly_sy <- poly_sy_micro[,.(N=sum(num_persons),
                            kish_N = sum(kish_N, na.rm=T),
                            sum_sample_weight=sum(count,na.rm=T),
                            mean=weighted.mean(x=edu_yrs,w=count)),
                         by=.(year,iso3,survey_name,nid,location_code,shapefile,survey_series)]

## Make national aggregates for diagnostics
natl_sy <- all_sy[,.(N=sum(num_persons),
                     kish_N = sum(kish_N, na.rm=T),
                     mean=weighted.mean(x=edu_yrs,w=count)),
                  by=.(year,iso3,survey_name,nid,survey_series)]

## Make DHS adm1 aggregates for comparison
all_dhs_sy <- all_sy[survey_series == "MACRO_DHS" & !is.na(lat) & !is.na(shapefile1) & shapefile1 != ""] #Keep points that are matched to both points and an admin1 poly from DHS
all_dhs_sy <- all_dhs_sy[,.(N=sum(num_persons),
                            kish_N = sum(kish_N, na.rm=T),
                            mean=weighted.mean(x=edu_yrs,w=count)),
                         by=.(year,iso3,survey_name,nid,survey_series, loc_code1, shapefile1)]

rm(all_sy)

###########################################################################################
###### 2. Collapse binned-year data to MBG geographies. ###################################
###########################################################################################

all_by_date <- most_recent_date("<<<< FILEPATH REDACTED >>>>", file_pattern="all_by_")
all_by <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>", all_by_date, ".rds"))

## Subset to target age/sex group.
## For CHN surveillance surveys, capping at 35 (which corresponds to age group 35-44) for 15-49, 15 (corresponds to 15-24) for 20-24
## For ZAF census, capping at 45 (which corresponds to age group 45-49) for 15-49, 20 (corresponds to 20-24) for 20-24
all_by <- all_by[age >= age_min & age <= ifelse(age_max == 49, 
                                                ifelse(nid %in% c(120994, 111916, 103981), 35, 
                                                       ifelse(nid == 12144, 45, age_max)), 
                                                ifelse(nid %in% c(120994, 111916, 103981), 15, 
                                                       ifelse(nid == 12144, 20, age_max))) 
                 & sex == target_sex, ]


## Drop any surveys that don't contain complete age range
# subset to ages age_start-age_max (and only surveys with this full range)
orig_nids <- unique(all_by$nid)
drop_by <- all_by[, as.list(range(age)), nid][V1 != age_min | V2 != age_max,]

all_by <- all_by[!nid %in% drop_by$nid | nid %in% c(124070, 124067, 120994, 111916, 103981, 283812, 283815, 283816, 12144),] # keep the tabulated data even though the age ranges are a little weird
missing_nid <- setdiff(orig_nids, all_by$nid)
if (length(missing_nid)>0) {
  print(paste("NIDs have been dropped for incomplete age-range:", paste(missing_nid, collapse=", ")))
}

############ Create edu_yrs for proportional indicators and subset data accordingly
if(grepl('prop',indicator_name)) {
  
  if(grepl('edu_zero_prop',indicator_name)) {
    all_by[, edu_yrs := ifelse(edu_yrs == 0, 1, 0)]
  }
  if(grepl('edu_no_primary_prop',indicator_name)) {
    all_by <- all_by[edu_yrs >= 1, ] ## Remove those who have zero education from denominator.
    all_by[, edu_yrs := ifelse(edu_yrs < 6, 1, 0)]
    
  }
  if(grepl('edu_primary_prop',indicator_name)) {
    all_by <- all_by[edu_yrs >= 6, ] ## Remove those who have less than primary from denominator.
    all_by[, edu_yrs := ifelse(edu_yrs < 12, 1, 0)] ## OF THOSE WHO HAVE AT LEAST PRIMARY, what proportion have less than secondary?
  }
} 

## Drop any observations missing education or weight
all_by <- all_by[!is.na(edu_yrs) & !is.na(proportion),]

## Collapse point data. Here use "proportion" to weight.
point_by_micro <- all_by[!is.na(lat) & !is.na(long) & lat!="" & long!="", ]
point_by <- point_by_micro[,.(N=sum(sample_size),
                              kish_N = sum(kish_N, na.rm=T),
                              sum_sample_weight=sum(proportion,na.rm=T),
                              mean=weighted.mean(x=edu_yrs,w=proportion)),
                           by=.(year,iso3,survey_name,nid,lat,long,survey_series)]

## Collapse polygon data. Here use "proportion" to weight.
poly_by <- all_by[!is.na(location_code) & !is.na(shapefile) & location_code!="" & shapefile!="" & (is.na(lat) | lat == ""), ]
poly_by <- poly_by[,.(N=sum(sample_size),
                      kish_N = sum(kish_N, na.rm=T),
                      sum_sample_weight=sum(proportion,na.rm=T),
                      mean=weighted.mean(x=edu_yrs,w=proportion)),
                   by=.(year,iso3,survey_name,nid,location_code,shapefile,survey_series)]

## Make national aggregates for diagnostics
natl_by <- all_by[,.(N=sum(sample_size),
                     kish_N = sum(kish_N, na.rm=T),
                     mean=weighted.mean(x=edu_yrs,w=proportion)),
                  by=.(year,iso3,survey_name,nid,survey_series)]

## Make DHS adm1 aggregates for comparison
all_dhs_by <- all_by[survey_series == "MACRO_DHS" & !is.na(lat) & !is.na(shapefile1) & shapefile1 != ""] #Keep points that are matched to both points and an admin1 poly from DHS
all_dhs_by <- all_dhs_by[,.(N=sum(sample_size),
                            kish_N = sum(kish_N, na.rm=T),
                            mean=weighted.mean(x=edu_yrs,w=proportion)),
                         by=.(year,iso3,survey_name,nid,survey_series, loc_code1, shapefile1)]

rm(all_by)

###########################################################################################
###### 3. Combine and run data coverage plots. ############################################
###########################################################################################

## Combine our four types of collapsed/mapped data: single-year point/poly, crosswalked (binned-year) point/poly.
all_data <- rbind(point_sy, poly_sy, point_by, poly_by, fill = T)
all_data <- all_data[, lat := as.numeric(lat)]
all_data <- all_data[, long := as.numeric(long)]
all_data <- all_data[, year := as.numeric(as.character(year))]

#### drop surveys that are over 50% missingness (only from ubCov, DHS, and WHS as of 4/26/2018)
to_be_dropped <- fread(paste0(indic_repo, "education/drop_for_missingness.csv"))
row <- ifelse(target_sex == 1, "drop.mn", "drop.wn")
nids_to_drop <- to_be_dropped[get(row) == 1, nid]
all_data <- all_data[!(nid %in% nids_to_drop)]
print(paste("NIDs have been dropped for missingness above 50% threshold:", paste(nids_to_drop, collapse=", ")))

## Make data aggregate plots over time by country for diagnostics.
all_natl <- all_data[, .(mean=weighted.mean(mean, sum_sample_weight)),by=.(year,iso3,survey_series,nid)]
all_natl <- all_natl[!is.na(survey_series) & !is.na(mean)]
all_natl <- all_natl[!(nid %in% nids_to_drop)]
all_natl[, survey_series := ifelse(survey_series == "JHSPH_PERFORMANCE_MONITORING_ACCOUNTABILITY_SURVEY_PMA2020", "PMA2020", survey_series)]
natl_polys <- rbind(poly_sy, poly_by)
natl_points <- rbind(point_sy, point_by)
natl_polys[, type := 'polygons']
natl_points[, type := 'points']
natl_geo_data <- rbind(natl_polys, natl_points, fill=TRUE)
natl_geo_data <- natl_geo_data[!(nid %in% nids_to_drop)]
pdf(paste0(mapped_datasets_dir, '/national_diagnostics_combined_', indicator_name, '.pdf'), height=6, width=10)
for(this_iso3 in sort(unique(all_natl[, iso3]))) {
  filler <- data.table(year = 1998:2017,
                       iso3 = rep(this_iso3, 19))
  this_all_natl <- all_natl[iso3==this_iso3, ]
  this_all_natl <- rbind(this_all_natl, filler, fill=TRUE)
  natl_labels <- all_natl[iso3==this_iso3 & !is.na(mean),c("year","mean","nid")]
  natl_labels <- natl_labels[!duplicated(natl_labels,by=c("nid","year"))]
  agg_gg <- ggplot() +
    geom_boxplot(data = natl_geo_data[iso3==this_iso3,],
                 aes(x = as.factor(as.numeric(year)),
                     y = mean,
                     color = type,
                     weight = N)) +
    geom_point(data = this_all_natl[iso3==this_iso3,],
               aes(x = as.factor(as.numeric(year)),
                   y = mean,
                   fill = survey_series),
               size = 10,
               shape = 21) +
    geom_text_repel(data = natl_labels,
                    aes(x = as.factor(as.numeric(year)),
                        y=mean,
                        label=nid, angle=90)) + 
    ylim(c(0,ifelse(grepl('prop',indicator_name), 1, 12))) + 
    labs(y=indicator_name, x='Year') + 
    ggtitle(this_iso3) + 
    theme_minimal()
  print(agg_gg)
}
dev.off()

rm(all_natl, natl_by, natl_sy, natl_polys, natl_points, point_sy, poly_sy, point_by, poly_by)

###########################################################################################
###### 4. Resample polygons to population-weighted pseudo-points. #########################
###########################################################################################

## Dropping shapefiles flagged by resample_polygons() as missing from central library.
## Please delete this drop once these shapefiles are added and work properly.
stages <- fread("<<<< FILEPATH REDACTED >>>>")
gauls <- stages[stages$Stage == "1" | stages$Stage == "2a" | stages$Stage == "2b",iso3]
# reformat to match resample polygons expected layout
all_data$location_code <- as.numeric(as.character(all_data$location_code))
all_data$shapefile <- as.character(all_data$shapefile)
setnames(all_data, 'lat', 'latitude')
setnames(all_data, 'long', 'longitude')
all_data[, cluster_id := 1:nrow(all_data)]

resampled_poly_data <- resample_polygons(data = all_data,
                                         cores = cores,
                                         indic = 'mean',
                                         pull_poly_method = "fast",
                                         gaul_list = get_adm0_codes(gauls, shapefile_version = 'current'),
                                         seed = 98121)

###########################################################################################
###### 5. Format and save as MBG input dataset. ###########################################
###########################################################################################

## Save a copy of full mapped dataset for diagnostics, dashboards, etc.
write.csv(all_data, paste0(mapped_datasets_dir, indicator_name, '.csv'), row.names=FALSE)

final_data <- copy(resampled_poly_data)
final_data <- final_data[!is.na(mean),]
setnames(final_data, 'mean', indicator_name)
setnames(final_data, 'iso3', 'country')
setnames(final_data, 'survey_series', 'source')

## If proportion, the indicator_name variable needs to be the count to be modelled.
if(grepl('prop',indicator_name)) final_data[, (indicator_name) := get(indicator_name) * N]

## Save to central MBG input folder on /share.
write.csv(final_data, file = paste0("<<<< FILEPATH REDACTED >>>>", indicator_name, ".csv"), row.names = FALSE)

## map resampled points to admin data
source(paste0(indic_repo, "/education/point_to_admin.r"))
admin_data <-  collapse_to_admin(final_data)
write.csv(admin_data, file = paste0(mapped_datasets_dir, indicator_name, "_resampled_points_with_admins.csv"), row.names = FALSE)
rm(final_data, admin_data)

## Format a version of the dataset for data coverage plot function.
coverage_data <- copy(all_data)
rm(all_data)
coverage_data$country <- substr(coverage_data$iso3, 1, 3)
coverage_data[, cluster_id := 1:nrow(coverage_data)]
setnames(coverage_data, 'survey_series', 'source')
coverage_data <- coverage_data[, nid := gsub('.DTA','',nid)]
coverage_data <- coverage_data[, nid := as.numeric(nid)]
coverage_data <- coverage_data[, survey := paste0(source,'_',country,'_',year)]
coverage_data <- coverage_data[source != "PMA", ]
coverage_data <- coverage_data[!is.na(mean),]

stages <- read.csv("<<<< FILEPATH REDACTED >>>>")
stages <- stages[stages$Stage == "1" | stages$Stage == "2a" | stages$Stage == "2b",]
gaul_to_loc <- read.csv("<<<< FILEPATH REDACTED >>>>")

stages <- merge(stages, gaul_to_loc, by ='loc_id')

###########################################################################################
### 6.  Data coverage plots ###############################################################
###########################################################################################
# only making coverage plots for these 4 indicators by default
if (indicator_name %in% c('edu_mean_15_49_female', 'edu_mean_15_49_male','edu_zero_prop_15_49_female','edu_no_primary_prop_15_49_female')) {
  ## Data coverage plot.
  library(stringr)
  source("<<<< FILEPATH REDACTED >>>>")
  source(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
  regions <- c('latin_america', 'middle_east', 'south_asia', 'se_asia', 'africa')
  for(r in regions) {
    if(!grepl('prop',indicator_name)) {
      coverage_maps <- graph_data_coverage_values(df = coverage_data,
                                                  var = 'mean',
                                                  title = 'Education',
                                                  legend_title = "Mean Years",
                                                  year_min = 1998,
                                                  year_max = 2017,
                                                  year_var = 'year',
                                                  region = r,
                                                  
                                                  cores = cores,
                                                  indicator = indicator_name,
                                                  core_repo= core_repo,
                                                  
                                                  color_scheme = "classic",
                                                  high_is_bad = FALSE,
                                                  cap = 12,
                                                  cap_type = "absolute",
                                                  stage_3_gray = TRUE,
                                                  simplify_polys = TRUE,
                                                  tolerance = 0.03,
                                                  
                                                  extra_file_tag = '',
                                                  out_dir = NULL,
                                                  save_on_share = FALSE,
                                                  log_dir = NULL,
                                                  
                                                  fast_shapefiles = TRUE,
                                                  new_data_plots = FALSE,
                                                  since_date = NULL,
                                                  annual_period_maps= FALSE,
                                                  save_period_maps = FALSE,
                                                  prep_shiny = FALSE,
                                                  return_maps = FALSE,
                                                  debug = FALSE,
                                                  
                                                  color_scheme_scatter = "brewer",
                                                  legend_min = NA,
                                                  legend_max = NA,
                                                  endemic_gauls = NULL,
                                                  base_font_size = 18,
                                                  map_point_size = 0.8,
                                                  poly_line_width = 0.2)
      
    }
    
    if(grepl('prop',indicator_name)) {
      coverage_maps <- graph_data_coverage_values(df = coverage_data,
                                                  var = 'mean',
                                                  title = 'Education',
                                                  legend_title = "Proportion",
                                                  year_min = 1998,
                                                  year_max = 2017,
                                                  year_var = 'year',
                                                  region = r,
                                                  
                                                  cores = cores,
                                                  indicator = indicator_name,
                                                  core_repo= core_repo,
                                                  
                                                  color_scheme = "classic",
                                                  high_is_bad = TRUE,
                                                  cap = 1,
                                                  cap_type = "absolute",
                                                  stage_3_gray = TRUE,
                                                  simplify_polys = TRUE,
                                                  tolerance = 0.03,
                                                  
                                                  extra_file_tag = '',
                                                  out_dir = NULL,
                                                  save_on_share = FALSE,
                                                  log_dir = NULL,
                                                  
                                                  fast_shapefiles = TRUE,
                                                  new_data_plots = FALSE,
                                                  since_date = NULL,
                                                  annual_period_maps= FALSE,
                                                  save_period_maps = FALSE,
                                                  prep_shiny = FALSE,
                                                  return_maps = FALSE,
                                                  debug = FALSE,
                                                  
                                                  color_scheme_scatter = "brewer",
                                                  legend_min = NA,
                                                  legend_max = NA,
                                                  endemic_gauls = NULL,
                                                  base_font_size = 18,
                                                  map_point_size = 0.8,
                                                  poly_line_width = 0.2,
                                                  
                                                  remove_rank = TRUE)
    }
    
  }
}



###########################################################################################
###### END ################################################################################
###########################################################################################