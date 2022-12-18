####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, resample
##              polygons, and create data coverage plots
## Indicator:   married or living as married among men and women, ages 15-49
####################################################################################################


## Setup -------------------------------------------------------------------------------------------
rm(list=ls())

# set arguments for this run and user.
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")
cores <- 10
shapefile_version <- 'current'
modeling_shapefile_version<-shapefile_version

# set arguments for this indicator/topic
indicator <- "in_union"
topic <- "rsp"

# load libraries & functions
source(paste0(core_repo,'/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
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

# Subset appropriately, now a function 
data <- subset_geomatched_in_union(data)

# store NIDs and sample size to compare at the end
nid_list <- unique(data$nid)
start_n <- data[, list(start = .N), by='country,nid']

# save microdata for indicator after all processing but before collapse
saveRDS(data, file = paste0('<<<< FILEPATH REDACTED >>>>', topic, "/",
                            "pre_collapse/", indicator, "_pre_collapse_", collapse_version, ".rds"))

## Collapse all data (polygon and point) -----------------------------------------------------------

data <- data[, .(int_year = floor(median(int_year, na.rm = T)),
                    in_union = weighted.mean(in_union, pweight),
                    N = sum(pweight) ^ 2 / sum(pweight ^ 2),
                    N_obs = .N,
                    sum_of_sample_weights = sum(pweight)),
             by=.(nid, country, source, year, point, shapefile, location_code, latitude, longitude)]

# replace year (the year the survey started in) with int_year (the median interview year for each location)
data[, year := NULL]
setnames(data, "int_year", "year")

# assign cluster_id (uniquely identifies a location-source pair)
data[, cluster_id := 1:.N]


# Make data coverage plots ------------------------------------------------------------------------
#load the data coverage plot functions
source(paste0(core_repo,'/mbg_central/setup.R'))

coverage <- copy(data)

# make the data coverage maps
coverage_maps <- graph_data_coverage_values(df = coverage,
                                            var = 'in_union',
                                            title = 'Currently married or living with partner',
                                            legend_title = "Prevalence",
                                            year_min = 1997,
                                            year_max = 2017,
                                            year_var = 'year',
                                            region = 'africa',

                                            cores = 1,
                                            indicator = 'in_union',
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

# Convert to counts and resample the polygon data -------------------------------------------------
data[, in_union := in_union * N]

#Separate point and polygon data
pt_data   <- data[!is.na(latitude) & !is.na(longitude)]
poly_data <- data[is.na(latitude) & is.na(longitude)]

# Drop these two problematic location_codes in shapefile (Audrey included in notes)
poly_data <- poly_data[shapefile != "NGA_IPUMS_adm2_2009" | !location_code %in% c(1099, 14099)]

# run polygon resampling
gaul_list = get_adm0_codes(unique(poly_data$country), shapefile_version = 'current')

resampled_poly_data <- 
  resample_polygons(data = poly_data,
                    indic = "in_union",
                    ignore_warnings = FALSE,
                    cores = cores,
                    pull_poly_method = "mclapply",
                    seed = 98121,
                    shapefile_version = shapefile_version,
                    gaul_list = gaul_list)


#Add columns to original point data to align with resampled polygon data
pt_data[, weight := 1]
pt_data[, pseudocluster := FALSE]

#Recombine point and resampled polygon data
data <- rbind(pt_data, resampled_poly_data)

# use resampling weights to down-weight N and in_union, then reset to 1. This is as an alternative
# to using the resampling weights as analytic weights in INLA.
data[, N := N * weight]
data[, in_union := in_union * weight]
data[, sum_of_sample_weights := sum_of_sample_weights * weight]
data[, weight := 1]

## Format, check, and save model input data --------------------------------------------------------
# subset and rename variables
data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, in_union, weight, sum_of_sample_weights, point, cluster_id)]

# check that no surveys were dropped
missing_nid <- setdiff(nid_list, data$nid)
if (length(missing_nid)>0) stop(paste("NIDs have been dropped:", paste(missing_nid, collapse=", ")))

# check that the sample size has not changed and the effective sample size is roughly in range
check_sample_sizes(data, start_n)

# # report differences compared to the most recent version
compare_most_recent(data, paste0(indicator, "_nonconservative"), prev_version)
# 
# ## save input data for model runs & log dates-------------------------------------------
save_data(data, paste0(indicator, "_nonconservative"), geomatched_version, collapse_version)

# drop nids that are borderline outliers for conservative version
data <- data[!nid %in% c(1912, 22950, 7375, 224223, 151797, 1981, 2209)]

# report differences compared to the most recent version
compare_most_recent(data, paste0(indicator, "_conservative"), collapse_version, prev_version)

save_data(data, paste0(indicator, "_conservative"), geomatched_version, collapse_version)

