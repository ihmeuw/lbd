####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   multiple sexual partners in the last year among women, ages 15-49
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- paste0('<<<< FILEPATH REDACTED >>>>', 'lbd_core')
indic_repo <- paste0('<<<< FILEPATH REDACTED >>>>', 'lbd_hiv')
cores <- 10
shapefile_version <- 'current'
modeling_shapefile_version<-shapefile_version

# set arguments for this indicator/topic
indicator <- "multiple_partners_year_WN"
topic <- "rsp"

# load libraries & functions
source(paste0(core_repo,'/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))

# find most recent versions
source(paste0(indic_repo, "mbg/functions/collapse_functions.r"))
geomatched_version <- most_recent_date(paste0('<<<< FILEPATH REDACTED >>>>', topic, "/"))
prev_version <- most_recent_date("<<<< FILEPATH REDACTED >>>>", file_pattern = indicator)
prev_version <- "2020_02_26"
collapse_version <- format(Sys.time(), "%Y_%m_%d")


## Load and subset data ----------------------------------------------------------------------------
# load geomatched data
data <- readRDS(paste0('<<<< FILEPATH REDACTED >>>>', topic, "/",
                       geomatched_version, "/", topic, "_", geomatched_version, ".RDS"))
data <- data.table(data)
data <- subset_geomatched_multiple_partner_year_wn(data)

# store NIDs and sample size to compare at the end
nid_list <- unique(data$nid)
start_n <- data[, list(start = .N), by='country,nid']

# save microdata for indicator after all processing but before collapse
saveRDS(data, file = paste0('<<<< FILEPATH REDACTED >>>>', topic, "/",
                            "pre_collapse/", indicator, "_pre_collapse_", collapse_version, ".rds"))


## Collapse all data (polygon and point) -----------------------------------------------------------
# collapse data by survey and location
data <- data[, list(int_year = floor(median(int_year, na.rm=T)),
                    multiple_partners_year = weighted.mean(multiple_partners_year, pweight),
                    N = sum(pweight)^2/sum(pweight^2),
                    sum_of_sample_weights = sum(pweight),
                    N_obs = .N),
             by='nid,country,source,year,point,shapefile,location_code,latitude,longitude']

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")

# assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]


## Make data coverage plots ------------------------------------------------------------------------
# load the data coverage plot functions
source(paste0(core_repo, '/mbg_central/setup.R'))

coverage <- copy(data)

# make the data coverage maps
coverage_maps <- graph_data_coverage_values(df = coverage,
                                            var = 'multiple_partners_year',
                                            title = 'Multiple sexual partners (women)',
                                            legend_title = "Prevalence",
                                            year_min = 1997,
                                            year_max = 2017,
                                            year_var = 'year',
                                            region = 'africa',

                                            cores = 1,
                                            indicator = 'multiple_partners_year',
                                            core_repo= core_repo,

                                            color_scheme = "RdPu",
                                            high_is_bad = TRUE,
                                            cap = 0.15,
                                            cap_type = "absolute",
                                            stage_3_gray = TRUE,
                                            simplify_polys = TRUE,
                                            tolerance = 0.01,

                                            extra_file_tag = '',
                                            out_dir = "<<<< FILEPATH REDACTED >>>>",
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
                                            legend_max = NA,
                                            legend_min = NA,
                                            endemic_gauls = NULL,
                                            base_font_size = 18,
                                            map_point_size = 0.8,
                                            poly_line_width = 0.2,

                                            remove_rank = TRUE)
rm(coverage)


## --------------------------------------------------------------------------------------------------

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, '<<<< FILEPATH REDACTED >>>>'))
mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))

# Convert to counts and resample the polygon data -------------------------------------------------
data[, multiple_partners_year := multiple_partners_year * N]

#Separate point and polygon data
pt_data   <- data[!is.na(latitude) & !is.na(longitude)]
poly_data <- data[is.na(latitude) & is.na(longitude)]

# run polygon resampling
gaul_list = get_adm0_codes(unique(poly_data$country), shapefile_version = 'current')
resampled_poly_data <- resample_polygons(data = poly_data, indic = "multiple_partners_year",
                                         cores = cores, pull_poly_method = "fast", seed = 98121,
                                         shapefile_version = shapefile_version,
                                         gaul_list = gaul_list)


#Add columns to original point data to align with resampled polygon data
pt_data[, weight:=1]
pt_data[, pseudocluster:=FALSE]

#Recombine point and resampled polygon data
data <- rbind(pt_data, resampled_poly_data)

# use resampling weights to down-weight N and multiple_partners_year, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, multiple_partners_year := multiple_partners_year * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, multiple_partners_year, weight, sum_of_sample_weights, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# report differences compared to the most recent version
compare_most_recent(data, indicator, collapse_version, prev_version)


## save input data for model runs & log dates-------------------------------------------
save_data(data, indicator, geomatched_version, collapse_version)
