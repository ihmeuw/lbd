####################################################################################################
## Description: Rake specific cause to GBD estimates for a specific cause by year, age and sex
##
## Passed args: main_dir [character] -- home directory for settings and final output
##              sex [numeric] -- sex of estimates to be raked
##              year [numeric] -- year of estimates to be raked
##              
## Requires:    output from aggregating to rake_level ([temp_dir]/est_[rake_level]_[year]_[sex].rdata)
##              raked argument with country iso3, GBD cause_id, GBD round id, and rake level that has been aggregated
##
## Outputs:     raking factors by year, sex, and age ([temp_dir]/rf.csv)
##              draws of raked estimates at area_var level ('[temp_dir]/draws_[area_var]_[year]_[sex]_raked.rdata')
##              mean, upper and lower of raked estimates at area_var level ('[temp_dir]/est_[area_var]_[year]_[sex]_raked.rdata')
##
## Run from within 'sae_central' directory!
####################################################################################################

library(data.table)
library(dplyr)
library(parallel)

rm(list=ls())

## Get settings and functions ----------------------------------------------------------------------
main_dir <- commandArgs()[4]
sex      <- as.numeric(commandArgs()[5])
year     <- as.numeric(commandArgs()[6])

source("settings.r")
get_settings(main_dir)

source("<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/post_estimation_functions.R")
source('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/misc_functions.R')
source("<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/prep_functions.R")
source("<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/shapefile_functions.R")
source("<<<< FILEPATH REDACTED >>>>/get_outputs.R")
source("<<<< FILEPATH REDACTED >>>>/get_age_metadata.R")
source("<<<< FILEPATH REDACTED >>>>/get_population.R")
source("post_estimation/calc_all_ages.r")
source("post_estimation/collapse.r")

# Get GBD output for country and cause -----------------
country      <- raked[["country"]]
cause_id     <- as.numeric(raked[["cause_id"]])
gbd_round_id <- as.numeric(raked[["gbd_round_id"]])
rake_level   <- raked[["rake_level"]]

# Match location ID to country 
loc_ids <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)$location_id

# Link between GBD age group and SAE ages
age_link <- 
  get_age_metadata(age_group_set_id = 12, gbd_round_id = 6) %>% 
  dplyr::filter(between(age_group_years_start, 5, 75)) %>% 
  dplyr::select(age_group_id, age = age_group_years_start)

# Grab GBD ids for under 5 and over 80 year olds
under5 <- data.table(age_group_id = 1, age = 0)
over80 <- data.table(age_group_id = 21, age = 80)
age_link <- rbind(under5, age_link, over80)

# Get GBD outputs and convert to proper age link
gbd_output <- 
  get_outputs(topic = "cause",
              version = "best",
              # decomp_step = "step5",
              gbd_round_id = 5,
              year_id = year,
              location_id = loc_ids,
              sex_id = sex,
              age_group_id = age_link$age_group_id,
              cause_id = cause_id,
              metric_id = 3,
              measure = 1) %>% 
  left_join(age_link, by = "age_group_id") %>%
  dplyr::select(age, year = year_id, sex = sex_id, loc_id = location_id, gbd = val)

# Grab unraked small area estimates
sae_est <-
  get(load(paste0(temp_dir, "est_", rake_level, "_", year, "_", sex, ".rdata"))) %>% 
  filter(age < 98) %>%
  dplyr::select(year, sex, age, loc_id = area, mean)

# Calculate raking factors 
rf <- 
  sae_est %>% 
  left_join(gbd_output, by = c("year", "sex", "age", "loc_id")) %>% 
  dplyr::mutate(rf = gbd / mean) %>% 
  dplyr::select(year, sex, age, loc_id, rf) %>% 
  data.table()

fwrite(rf, file = paste0(temp_dir, "rf_", year, "_", sex, ".csv"))

# Load draws without crude and age-standardized rates
load(paste0(temp_dir, "draws_area_", year, "_", sex, ".rdata"))
draws <- draws[age < 98,]

# Load area to rake level link (most important for subnationals) and merge
load(geoagg_files[[rake_level]])
area_to_rake_level <- unique(weights[, .(area, loc_id = get(rake_level))])
draws <- merge(draws, area_to_rake_level, by = "area")

# Now merge raking factors on draws and calculate raked draws 
draws <- merge(draws, rf, by = c("year", "sex", "age", "loc_id"))
draws[, value := value * rf]
draws[, rf := NULL]
draws[, loc_id := NULL]

# Load standard weight and pop for age-standardizing weights
std_wt <- fread(age_std_file)[, list(age, wt = wt/sum(wt))]
load(pop_file)

# Add crude and age-standardized rates to draws using pop
draws <- merge(draws, pop, by = c("area", "year", "sex", "age"))
draws <- calc_all_ages(draws, std_wt)
    
# Collapse to make estimates
est <- collapse_draws(draws, c("level", "area", "year", "sex", "age"))

# Save draws and estimates 
save(draws, file = paste0(temp_dir, "draws_area_", year, "_", sex, "_raked.rdata"))
save(est, file = paste0(temp_dir, "est_area_", year, "_", sex, "_raked.rdata"))
message(paste0("Done saving raked draws for sex ", sex, " in year ", year))


