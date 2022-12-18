####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   circumcision in men, ages 15-49
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- paste0('<<<< FILEPATH REDACTED >>>>')
indic_repo <- paste0('<<<< FILEPATH REDACTED >>>>')
cores <- 10
shapefile_version = 'current'
modeling_shapefile_version <- shapefile_version

# set arguments for this indicator/topic
indicator <- "male_circumcision"
topic <- "male_circumcision"
shapefile_version <- "current"

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, '<<<< FILEPATH REDACTED >>>>'))
mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))

# find most recent versions
source(paste0(indic_repo, "mbg/functions/collapse_functions.r"))
geomatched_version <- most_recent_date(paste0('<<<< FILEPATH REDACTED >>>>', topic, "/"))
prev_version <- most_recent_date("<<<< FILEPATH REDACTED >>>>/archive/hiv/", file_pattern = indicator)
prev_version <- "2020_02_26"
collapse_version <- format(Sys.time(), "%Y_%m_%d")

## Load and subset data ----------------------------------------------------------------------------
# load geomatched data
data <- readRDS(paste0('<<<< FILEPATH REDACTED >>>>', topic, "/",
                       geomatched_version, "/", topic, "_", geomatched_version, ".RDS"))
data <- data.table(data)
data_for_crosswalk <- copy(data) # used for cross-walking function

# Subset geomatched data
data <- subset_geomatched_male_circumcision(data)

# store NIDs and sample size to compare at the end
nid_list <- unique(data$nid)

#run this start_n
start_n <- data[, list(start = .N), by = 'country,nid']

# save microdata for indicator after all processing but before collapse
saveRDS(data, file = paste0('<<<< FILEPATH REDACTED >>>>', topic, "/",
                            "pre_collapse/", indicator, "_pre_collapse_", collapse_version, ".rds"))


## Collapse all data (polygon and point) -----------------------------------------------------------
data <- data[, .(int_year = floor(median(int_year, na.rm = T)),
                 male_circumcision = weighted.mean(male_circumcision, pweight),
                 N = sum(pweight) ^ 2 / sum(pweight ^ 2),
                 sum_of_sample_weights = sum(pweight),
                 N_obs = .N),
             by = .(nid, country, source, year, point, shapefile, location_code, latitude, longitude)]

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")

## Add survey report data ---------------------------------------------------------------------------
report <- fread(paste0(indic_repo, "<<<< FILEPATH REDACTED >>>>reports/circumcision_extraction.csv")) 

# Get male circumcision in the correct format, divide by known since we drop unknown in micro data
report <- 
  report %>% 
  filter(!is.na(location_code)) %>% # Drop those missing location code, those should only be admin0
  mutate(mc_unknown = ifelse(is.na(mc_unknown), 0, mc_unknown)) %>% 
  mutate(male_circumcision = male_circumcision / (100 - mc_unknown)) %>% 
  data.table()

# Remove survey reports that don't overlap 15-49
report <- report[(start_age <= 15 & end_age >= 15) | (start_age >= 15 & start_age <= 49), ]

# Remove if sample size is 0
report <- report[N != 0]

# Add survey data that is not in standard format
# identify all age ranges other than 15-49, subset report to covering those age_ranges
age_list <- unique(report[,c("start_age", "end_age")])

#Only add age lists that cover at least 20 years, this is just a rule of thumb since crosswalking for less than that does not seem to work
age_list <-
  age_list %>%
  filter(start_age != 15 | end_age != 49)

crosswalk_report_list <- apply(age_list, 1, function(x) return(report[start_age == x[1] & end_age == x[2]]))

# Prepare microdata for crosswalk
source(paste0(indic_repo, "/mbg/functions/crosswalk_functions.r"))

# only keep data from Africa since reports are in SSA (check to make sure this doesn't change)
africa_countries <- 
  data.table(country = unique(data_for_crosswalk$country)) %>% 
  mutate(ad0_code = get_adm0_codes(country)) %>% 
  filter(ad0_code %in% get_adm0_codes("africa")) %>% 
  pull(country)

data_for_crosswalk <- data_for_crosswalk[country %in% africa_countries, ]
data_for_crosswalk <- crosswalk_data_prepare(data_for_crosswalk, indicator, "pweight", shapefile_version = shapefile_version)

# Apply crosswalk
crosswalked <-
  rbindlist(lapply(copy(crosswalk_report_list), function(x) {
    crosswalk_reports(data_for_crosswalk, x, indicator, "pweight", indicator)
  }))

setnames(crosswalked, c("prev_old", "prev_new"), c("male_circumcision_old", "male_circumcision_new"))
rm(data_for_crosswalk)

# save cross-walk data for future reference
write.csv(crosswalked, file = paste0('<<<< FILEPATH REDACTED >>>>/male_circumcision/crosswalk_output/',
                                     collapse_version, '_crosswalked_data.csv'), row.names = F)

# quick plot to show N vs N_eff for cross-walked report data
pdf('<<<< FILEPATH REDACTED >>>>crosswalk_output/sample_size.pdf')
crosswalked %>%
  ggplot(aes(N, N_eff, fill = male_circumcision_old)) +
  geom_point(size = 2, shape = 21, color = "grey") +
  geom_abline(slope = 1, intercept = 0) +
  coord_equal(ylim = c(0, max(crosswalked$N)), xlim = c(0, max(crosswalked$N))) +
  scale_fill_distiller("Previous circumcision\nprevalence", palette = "YlOrRd",
                       direction = 1, limits = c(0, 1.01)) +
  labs(x = "Previous sample size", y = "Rescaled sample size to reflect uncertainty") +
  theme_bw()
dev.off()

# recombine cross-walked report data with data for 15-49
crosswalked[, c("male_circumcision", "N") := list(male_circumcision_new, N_eff)]
report <- rbind(report[start_age == 15 & end_age == 49,], crosswalked, fill = T)

# drop national estimates
report <- report[location_type != "National" | admin_level == 0,]

# # drop countries outside Africa
# report <- report[country %in% loc,]

# drop any NIDs we have microdata for
report <- report[!nid %in% data$nid,]

# approximate the effective sample size using the design effects from the polygon microdata
# Note: this is *very* rough; the purpose is to make sure the report data doesn't have an unfair
# advantage compared to microdata; the calculation below excludes surveys where there is no apparent
# design effect -- this is not super plausible, and we don't want these to skew the distribution.
de <- data[point == 0, list(de = sum(N)/sum(N_obs)), by='nid,shapefile,location_code'][de < 0.99999, median(de)]
report[, N_obs := N]
report[, N := N_obs * de]

# update the survey list
nid_list <- c(nid_list, unique(report$nid))

# combine with collapsed microdata
report <- report[, list(nid, country, survey_series, int_year, point, shapefile, location_code, latitude, longitude, male_circumcision, N, N_obs)]
setnames(report, c("int_year", "survey_series"), c("year", "source"))

# Mark report data, and make sum of sample weights equal to N to approximate later in visualization plots
report[, sum_of_sample_weights := N]
report[, report := 1]
data <- rbind(data, report, fill=T)

# now that all data is compiled, assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]


## Make data coverage plots ------------------------------------------------------------------------
# load the data coverage plot functions
#setup
  source(paste0(core_repo, 'mbg_central/graph_data_coverage.R'))
  coverage <- copy(data)


# # # make the data coverage maps
#   #regions <- c('africa', 'latin_america', 'south_asia', 'se_asia')
#     for(reg in regions) {
      coverage_maps <- graph_data_coverage_values(df = coverage,
                                           var = 'male_circumcision',
                                           title = 'Male Circumcision',
                                           legend_title = "Prevalence",
                                           year_min = 2000,
                                           year_max = 2017,
                                           year_var = 'year',
                                           region = 'africa',

                                           cores = 1,
                                           indicator = 'male_circumcision',
                                           core_repo= core_repo,

                                           color_scheme = "RdPu",
                                           high_is_bad = TRUE,
                                           cap = 1,
                                           cap_type = "none",
                                           stage_3_gray = TRUE,
                                           simplify_polys = TRUE,
                                           tolerance = 0.01,

                                           extra_file_tag = '',
                                           save_on_share = FALSE,
                                           log_dir = NULL,
                                           out_dir = "'<<<< FILEPATH REDACTED >>>>'",

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



rm(coverage)

## Convert to counts and resample the polygon data -------------------------------------------------
data[, male_circumcision := male_circumcision * N]

#For some reason we need this?
modeling_shapefile_version <- shapefile_version

#Separate point and polygon data
pt_data   <- data[!is.na(latitude) & !is.na(longitude)]
poly_data <- data[is.na(latitude) & is.na(longitude)]

gaul_list = get_adm0_codes(unique(poly_data$country), shapefile_version = 'current')
resampled_poly_data <- resample_polygons(data = poly_data, 
                                         indic = "male_circumcision", 
                                         ignore_warnings = FALSE,
                                         cores = cores,
                                         pull_poly_method = "fast",
                                         seed = 98121,
                                         shapefile_version = shapefile_version,
                                         gaul_list = gaul_list)

#Add columns to original point data to align with resampled polygon data
pt_data[, weight := 1]
pt_data[, pseudocluster:=FALSE]

#Recombine point and resampled polygon data
data <- rbind(pt_data, resampled_poly_data)

# use resampling weights to down-weight N and male_circumcision, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, male_circumcision := male_circumcision * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, male_circumcision, weight, sum_of_sample_weights, report, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# This saves a text file to check for differences
compare_most_recent(data, indicator, collapse_version, prev_version)

## save input data for model runs & log dates-------------------------------------------
save_data(data, indicator, geomatched_version, collapse_version)
