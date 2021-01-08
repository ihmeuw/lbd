### This script is used to prepare all model inputs for HIV small area estimation model
## This script is used for all data prep functions 
## Includes completeness functions

rm(list = ls())

# Load libraries
library(dplyr)
library(data.table)
library(sf)
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)
library(scales)
library(rgdal)
library(parallel)

# Source functions from core repo
core_repo <- paste0("<<<< FILEPATH REDACTED >>>>")
source(paste0(core_repo, "data_central/vr_functions.R"))
source(paste0(core_repo, "data_central/sae_functions.R"))
source(paste0(core_repo, "mbg_central/fractional_raking_functions.R"))
source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/post_estimation_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/prep_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/shapefile_functions.R')
source(paste0(core_repo, "mbg_central/misc_functions.R"))

# Source functions from HIV repo 
function_files = list.files(paste0("<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/functions/"), full.names = T)
sapply(function_files, source)

# Remap some functions
rename <- dplyr::rename
summarize <- dplyr::summarize
select <- dplyr::select

# Source central function
source("<<<< FILEPATH REDACTED >>>>/get_age_metadata.R")

## Pick country to run over 
country <- "gtm"
rerun_completeness <- T

## All other parameters are fixed for this analysis
admin <- "adm2"
sexes <- c(1, 2)
ages <- seq(0, 80, 5)
shapefile_field <- "uid"
cores <- 15


## Begin script --------------------------------------------
if (length(get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)$location_id) > 1) {
  rake_level = "admin1"
} else {
  rake_level <- "national"
}

# Pull years from formatted events data (had an error in the metadata before)
years <- 
  vr_pull_cod(paste0(country,"_", admin), "formatted") %>% 
  pull(year_id) %>% 
  sort() %>% 
  unique() 

# Grab shapefile path and load shapefile
shapefile_path <- vr_pull_shp(paste0(country, "_adm2"), "stable", "full")$shapefile
shape <- st_read(shapefile_path, quiet = T) 

# Mark date of data prep and create directory to save data inputs
data_prep_date <- format(Sys.time(), '%Y_%m_%d')
out_dir   <- paste0("<<<< FILEPATH REDACTED >>>>")
figure_dir <- paste0(out_dir, "/figures/")
dir.create(out_dir, recursive = T)
dir.create(figure_dir)

## Prep population data ------------------------------------------------------------------------------
# Pop over admin 2  --------------
area_shape_path <-  
  vr_pull_shp(source = paste0(country, "_adm2"),
              type = "annual",
              year = "last",
              admin_level = NULL,
              link = F)$shapefile 

pop_list_area <- 
  vr_get_population(ages = ages,
                    country = country,
                    admin = admin,
                    sexes = sexes,
                    years = years,
                    core_repo = core_repo,
                    raked = T,
                    shapefile_path = area_shape_path,
                    cores = cores,
                    worldpop_release = "2020_03_20",
                    shapefile_field = "GAUL_CODE",
                    link_table = "custom")

pop_admin2    <- pop_list_area[["pop_raked"]]
setnames(pop_admin2, "area", "admin2")
pop_rf <- pop_list_area[["rf"]]

# Compare pop over uid and pop over admin2 ----------------
# Map admin2 to uids
uid_area_map <- 
  vr_pull_loc_meta(paste0(country, "_adm2")) %>% 
  filter(level == 2, year_end == max(year_end)) %>% 
  dplyr::select(uid, admin2 = location_id) %>% 
  arrange(uid, admin2) %>% data.table()

pop <- 
  merge(pop_admin2, uid_area_map)[, .(pop = sum(pop)), by = c("uid", "year", "age", "sex")] %>% 
  rename(area = uid) %>% 
  arrange(area, year, age, sex) %>% 
  data.table()

# Make sure not missing any admin units data 
are_equal(sort(as.numeric(as.character(shape[[shapefile_field]]))), sort(unique(pop$area)))

# Check if any population is zero, see if its worth dropping 
if (nrow(filter(pop, pop == 0))) {
  admin_nopop <- filter(pop, pop == 0) %>% pull(area) %>% unique()
  message(paste0("Area ", admin_nopop, " is missing worldpop information, consider dropping from analysis"))
} 

# Save population and population by admin2 
save(pop, file = paste0(out_dir, "pop.rdata"))
save(pop_admin2, file = paste0(out_dir, "pop_admin2.rdata"))

# Save population raking factor plot
vr_plot_raking_factors(pop_rf, paste0(figure_dir, "raking_factor.pdf"))

# Plot population change over time
total_pop <- pop[, .(pop = sum(pop)), .(area, year)]
model_results_compare_sae(country = country,
                          years = years,
                          save_file = paste0(figure_dir, "pop_over_time.pdf"),
                          covariate_data = total_pop)

## Create geographic aggregation file ------------------------------------------------------------------
nat_filepath <- paste0(out_dir, "national_estimates_age_sex.rdata")
admin1_filepath <- paste0(out_dir, "admin1_estimates_age_sex.rdata")

# Make national weights from populatoin data, have location id be the number for the weights
nat_loc_id <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = F)$location_id
weights <- copy(pop)
weights[, "national" := nat_loc_id]
save(weights, file = nat_filepath)
  
# Make admin1 aggregation file
weights <- 
  vr_pull_loc_meta(paste0(country, "_", admin)) %>% 
  filter(level == 2, year_end == max(year_end)) %>% # Gets most recent map to 
  dplyr::select(area = uid, admin2 = location_id, admin1 = location_parent_id) %>% 
  left_join(pop_admin2, by = "admin2") %>% # Add admin2 populations
  group_by(area, year, age, sex, admin1) %>% # Summarize across areas for case when area is two admin2's that either are either
  summarize(pop = sum(pop)) %>%              # in the same admin1 or not
  ungroup() %>% 
  data.table()

save(weights, file = admin1_filepath)

# Prep covariates -------------------------------------------------------------------------
settings_folder <- "<<<< FILEPATH REDACTED >>>>"

# Get covariate config
cov_config <-
  fread(paste0(settings_folder, "cov_list.csv")) %>%
  filter(include == T) %>%
  data.table()

# Load fractionally aggregated covariates using stable shapefile
covar <- 
  frac_agg_covs(cov_config = cov_config,
                years = years,
                shapefile_path = shapefile_path,
                shapefile_field = shapefile_field,
                core_repo = core_repo,
                cores = cores,
                worldpop_age_sex_release = "2020_03_20")

# Transform worldpop into pop_density
shapefile_area <- 
  shape %>% 
  mutate(area_km2 = st_area(geometry) %>% units::set_units("km^2") %>% as.numeric) %>% 
  dplyr::select(shapefile_field, area_km2) %>%
  st_set_geometry(NULL) %>%
  mutate_at(vars(shapefile_field), funs(as.numeric(as.character(.)))) %>% 
  data.table()

covar <- 
  covar %>% 
  left_join(shapefile_area, by = shapefile_field) %>% 
  mutate(pop_density = worldpop / area_km2) %>% 
  dplyr::select(area = shapefile_field, year, access, nightlights = dmspntl, urbanicity = ghslurbanicity, pop_density) %>% 
  data.table()

save(covar, file = paste0(out_dir, "covariates_sae.rdata"))

#Plot aggregated covariates
vr_plot_aggregated_covs(covar, shapefile_path, shapefile_field, paste0(figure_dir, "/aggregated_covs.pdf"))

# Plot covariates changes by year, similar to population
mclapply(setdiff(names(covar), c("area", "year")), mc.cores = cores, function(cov) {
  covariate_data <- dplyr::select(covar, area, year, cov)
  model_results_compare_sae(country = country,
                            years = years,
                            save_file = paste0(figure_dir, cov, "_over_time.pdf"),
                            covariate_data = covariate_data)
})

# Prep mortality data --------------------------------------------------------------------------------
events <- vr_prep_hiv_mortality(country, admin, years)
save(events, file = paste0(out_dir, "mort_events.rdata"))

# Check for missing deaths
# Make sure there are no age-year-sex combinations with no deaths, this might be a red flag
events[, .(events = sum(events)), .(sex, year, age)] %>% filter(events == 0) 

# Make sure there are no UID missing data, if there are inspect which ones
events[, .(events = sum(events)), .(area)] %>% filter(events == 0)

# Plot mortality data 
vr_plot_mortality_data(events, save_file = paste0(figure_dir, "mortality_data_visualize.pdf"))

# Save shapefile -----------------------------------------------------------------
# Need to load as readOGR here for some reason 
shapefile <- readOGR(shapefile_path)[shapefile_field]
save(shapefile, file = paste0(out_dir, "shape.rdata"))

# Prepare adjacency matrix ------------------------------------------------------
neighbors <- spdep::poly2nb(shapefile)
adjmat <- as(spdep::nb2mat(neighbors, style="B", zero.policy = T), "dgTMatrix")
colnames(adjmat) <- rownames(adjmat) # So they are symmetric
save(adjmat, file = paste0(out_dir, "adj_matrix.rdata"))

# Make plot of adjacency matrix --------------------
coords <- st_coordinates(st_centroid(st_geometry(shape)))
pdf(paste0(figure_dir, "adjacency_matrix.pdf"), height = 10, width = 10)
plot(st_geometry(shape), border = "grey")
plot(neighbors, coords, add = TRUE, col = "blue")
dev.off()

# Prepare adjusted age weights --------------------------------------------------
age_wts <- 
  get_age_metadata(age_group_set_id = 12, gbd_round_id = 5) %>% 
  dplyr::select(wt  = age_group_weight_value, age = age_group_years_start) %>%
  mutate(age = case_when(age < 5 ~ 0,
                         between(age, 5, 75) ~ age,
                         age >= 80 ~ 80)) %>% 
  group_by(age) %>% 
  dplyr::summarize(wt = sum(wt)) %>% 
  ungroup() %>% 
  write_csv(paste0(out_dir, "age_wts.csv"))


# Save settings folder -----------------------------------------------------
# Add sexes and age grouping (Should this be read in in settings?)

settings <- read_csv("<<<< FILEPATH REDACTED >>>>", col_names = F) %>% data.table()

settings[X1 == "area_var",     X2 := "'area'"]
settings[X1 == "years",        X2 := paste0(min(years), ":", max(years))]
settings[X1 == "pop_file",     X2 := paste0(out_dir, "pop.rdata")]
settings[X1 == "covar_file",   X2 := paste0(out_dir, "covariates_sae.rdata")]
settings[X1 == "covar_file",   X2 := paste0(out_dir, "covariates_sae.rdata")]
settings[X1 == "events_file",  X2 := paste0(out_dir, "mort_events.rdata")] 
settings[X1 == "covars",       X2 :=  "c('access', 'urbanicity', 'nightlights', 'pop_density')"]
settings[X1 == "age_std_file", X2 :=  paste0(out_dir, "age_wts.csv")]
settings[X1 == "shape_file",   X2 :=  paste0(out_dir, "shape.rdata")]
settings[X1 == "adjmat_file",  X2 :=  paste0(out_dir, "adj_matrix.rdata")]
settings[X1 == "geoagg_files", X2 :=  paste0('c(national = \"', nat_filepath, '\", admin1 = \"', admin1_filepath, '\")')]
settings[X1 == "raked", X2 :=  paste0('c(country = \"', country, '\", cause_id = \"', 298,
                                      '\", gbd_round_id = \"5\", rake_level = \"', rake_level, '\")')]
fwrite(settings, paste0(out_dir, "/settings.csv"))

# Prep completeness data ----------------------------------------------------------------------------
# Run completeness by country, this takes a long time so only run if needed
# Decide on adult or child completeness based on country
if (country == "gtm") {
  calc_child_comp <- F
  calc_adult_comp <- T
} else if (country %in% c("col", "mex")) {
  calc_child_comp <- T
  calc_adult_comp <- F
} else if (country %in% c("ecu", "bra")) {
  calc_child_comp <- T
  calc_adult_comp <- T
} else {
  calc_child_comp <- F
  calc_adult_comp <- F
}

if (calc_child_comp | calc_adult_comp) {
  if (rerun_completeness) mem <- 500 else mem <- 30
  if (!country %in% c("bra", "mex")) {
    qsub <- paste('qsub -e <<<< FILEPATH REDACTED >>>> -l',
                paste0('m_mem_free=', mem, 'G -l fthread=5 -N'),
                paste0('vr_completeness_prep_', country),
                '-P proj_geo_nodes -q geospatial.q -l h_rt=08:00:00:00 -l archive=TRUE',
                '-v sing_image=<<<< FILEPATH REDACTED >>>> -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0('<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/vr/1_data_processing/vr_completeness_prep.r'),
                country, rerun_completeness, calc_child_comp, calc_adult_comp)
    system(qsub)
    settings <- rbind(settings, data.table(X1 = "admin1_completeness", X2 = "FALSE"))
  } else {
    qsub <- paste('qsub -e <<<< FILEPATH REDACTED >>>> -o <<<< FILEPATH REDACTED >>>> -l',
                  paste0('m_mem_free=', mem, 'G -l fthread=5 -N'),
                  paste0('vr_completeness_prep_', country),
                  '-P proj_geo_nodes -q geospatial.q -l h_rt=08:00:00:00 -l archive=TRUE',
                  '-v sing_image=<<<< FILEPATH REDACTED >>>> -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                  paste0('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                  paste0('<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/vr/1_data_processing/vr_mex_bra_completeness.r'),
                  country, rerun_completeness, calc_child_comp, calc_adult_comp)
    system(qsub)
    settings <- rbind(settings, data.table(X1 = "admin1_completeness", X2 = "TRUE"))
  }
  # Additional settings for completeness
  settings <- rbind(settings, data.table(X1 = "completeness_prior_file", X2 = paste0(out_dir, "/comp_prior_dist.RDS")))
  settings <- rbind(settings, data.table(X1 = "calc_child_comp", X2 = calc_child_comp))
  settings <- rbind(settings, data.table(X1 = "calc_adult_comp", X2 = calc_adult_comp))
  
  fwrite(settings, paste0(out_dir, "/settings.csv"))
}

# Run final data check on settings (will only work after completeness prep is run, may take awhile) ----------------------------------------------------
source("<<<< FILEPATH REDACTED >>>>/lbd_core/sae_central/settings.r")
setwd("<<<< FILEPATH REDACTED >>>>/lbd_core/sae_central/")
get_settings(out_dir)
check_settings(out_dir)

# Make dated log
date_log <- data.table(input_date = data_prep_date)
if (file.exists(paste0("<<<< FILEPATH REDACTED >>>>"))) {
  write_csv(date_log, path = paste0("<<<< FILEPATH REDACTED >>>>"))
} else {
  write_csv(date_log, path = paste0("<<<< FILEPATH REDACTED >>>>"), append = TRUE)
}

message("Done with model prep for ", country)

