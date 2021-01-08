####################################################################################################
## Description: Calculate draws for aggregated geographies based on a cross-walk file from
##              'area_var' to 'level'.
##
## Passed args: main_dir [character] -- home directory for settings and final output
##              year [integer] -- year to process
##              sex [integer] -- sex to process
##              level [character] -- the geographic level to aggregate draws to
##              rake [logical] -- are estimates raked? If raked, output file is tagged with `_raked.rdata`
##
## Requires:    crosswalks from 'area_var' to 'level' (geoagg_files[level])
##              draws ('[temp_dir]/draws_[area_var]_[year]_[sex].rdata')
##
## Outputs:     draws at the 'level' level, formatted identically to the draws input
##                ('[temp_dir]/draws_[level]_[year]_[sex].rdata')
##
## Run from within 'sae_central' directory!!
####################################################################################################
library(data.table)

rm(list=ls())

## Get settings and functions ----------------------------------------------------------------------
main_dir <- commandArgs()[4]
year     <- as.integer(commandArgs()[5])
sex      <- as.integer(commandArgs()[6])
level    <- commandArgs()[7]
rake     <- as.logical(commandArgs()[8])

source("settings.r")
get_settings(main_dir)

source("post_estimation/calc_all_ages.r")
source("post_estimation/collapse.r")

## Load and format the crosswalk files for 'level' -------------------------------------------------
# load and subset the crosswalk file
load(geoagg_files[level])

weights <- weights[sex == get("sex", .GlobalEnv) & year == get("year", .GlobalEnv),
                   list(area = get(area_var), area2 = get(level), age, pop)]

# calculate aggregation weights
weights[pop < 1e-5, pop := 0] # change effective 0s to 0s to avoid some annoying edge cases where one area gets all the weight even if it also has effectively no population
weights[, wt := pop/sum(pop), by='area2,age']

# when the total population in an area2-age is 0, the weights are undefined; in this case, derive
# weights from the area-specific population totals across all ages
weights[, total := sum(pop), by='area2,area']
weights[is.na(wt), wt := total/sum(total), by='area2,age']
weights[, total := NULL]

# load the age standard
std_wt <- fread(age_std_file)[, list(age, wt = wt/sum(wt))]

## Aggregate draws to 'level' and collapse to generate point estimates -----------------------------
# load draws
load(paste0(temp_dir, "/draws_", area_var, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))

# drop all ages, we recalculate these after aggregating by geography
draws <- draws[age < 98,]

# aggregate draws to 'level'
draws[, level := get("level", .GlobalEnv)]
draws <- merge(draws, weights, by=c("area", "age"), allow.cartesian=T)
draws <- draws[, list(value = sum(value*wt), pop = sum(pop)), by='level,area2,year,sex,age,sim']
setnames(draws, "area2", "area")

# add crude and age-standardized rates
draws <- calc_all_ages(draws, std_wt)

# collapse to get estimates
est <- collapse_draws(draws, c("level", "area", "year", "sex", "age"))

# save draws and estimates
save(draws, file = paste0(temp_dir, "/draws_", level, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))
save(est, file = paste0(temp_dir, "/est_", level, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))

# Save pi*mort to compare to data, run only when unraked to not repeat
if (grepl("_c", model) & !rake) {
  load(paste0(temp_dir, "/draws_pi_mort_", area_var, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))
  
  # drop all ages, we recalculate these after aggregating by geography
  draws <- draws[age < 98,]
  
  # aggregate draws to 'level'
  draws[, level := get("level", .GlobalEnv)]
  draws <- merge(draws, weights, by=c("area", "age"), allow.cartesian=T)
  draws <- draws[, list(value = sum(value*wt), pop = sum(pop)), by='level,area2,year,sex,age,sim']
  setnames(draws, "area2", "area")
  rm(weights)
  
  # add crude and age-standardized rates
  draws <- calc_all_ages(draws, std_wt)
  
  # collapse to get estimates
  est <- collapse_draws(draws, c("level", "area", "year", "sex", "age"))
  
  # save draws and estimates
  save(draws, file = paste0(temp_dir, "/draws_pi_mort_", level, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))
  save(est, file = paste0(temp_dir, "/est_pi_mort_", level, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))
}
