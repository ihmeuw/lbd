####################################################################################################
## Description: Calculate draws for both sexes combined by aggregating the draws for males and
##              females using population-weighted averages.
##
## Passed args: main_dir [character] -- home directory for settings and final output
##              year [integer] -- year to calculate both sexes combined draws for
##              level [character] -- the geographic level to aggregate draws for
##              rake [logical] -- are estimates raked? If raked, output file is tagged with `_raked.rdata`
##
## Requires:    populations (pop_file)
##              draws ('[temp_dir]/draws_[level]_[year]_[sex].rdata')
##
## Outputs:     draws for both sexes combined ('[temp_dir]/draws_[level]_[year]_3.rdata')
##              estimates for both sexes combined ('[temp_dir]/est_[level]_[year]_3.rdata')
##
## Run from within 'sae_central' directory!!
####################################################################################################

library(data.table)

rm(list=ls())

## Get settings and functions ----------------------------------------------------------------------
main_dir <- commandArgs()[4]
year     <- as.integer(commandArgs()[5])
level    <- as.character(commandArgs()[6])
rake     <- as.logical(commandArgs()[7])

source("settings.r")
get_settings(main_dir)

source("post_estimation/calc_all_ages.r")
source("post_estimation/collapse.r")

## Load and format populations ---------------------------------------------------------------------
# load the populations file (pop_file for area_var, geo_agg_file for anything else)
if (level == area_var) {
  load(pop_file)
  weights <- pop[year == get("year", .GlobalEnv), c(level, "sex", "age", "pop"), with=F]
  rm(pop)
} else {
  load(geoagg_files[level])
  weights <- weights[year == get("year", .GlobalEnv), list(pop = sum(pop)), by=c(level, "sex", "age")]
}
setnames(weights, level, "area")

# calculate aggregation weights
weights[pop < 1e-5, pop := 0] # change effective 0s to 0s to avoid some annoying edge cases where one sex gets all the weight even if it also has effectively no population
weights[, wt := pop/sum(pop), by='area,age']

# when the total population in an area-age is 0, the weights are undefined; in this case, derive
# weights from the sex-specific population totals across all areas
weights[, total := sum(pop), by='sex,age']
weights[is.na(wt), wt := total/sum(total), by='area,age']
weights[, total := NULL]

# load the age standard
std_wt <- fread(age_std_file)[, list(age, wt = wt/sum(wt))]

## Aggregate draws to both sexes combined and collapse to generate point estimates -----------------
# load draws
draws <- rbindlist(lapply(sexes, function(sex) {
  load(paste0(temp_dir, "/draws_", level, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))
  draws
}))

# drop all-ages rates, we recalculate these after aggregating by sex
draws <- draws[age < 98,]

# if necessary, add in 0s for sexes we don't model
if (length(sexes) == 1) {
  temp <- copy(draws)
  temp[, value := 0]
  if (sexes == 1) temp[, sex := 2] else temp[, sex := 1]
  draws <- rbind(draws, temp)
  setkeyv(draws, c("level", "area", "year", "sex", "age", "sim"))
  rm(temp); gc()
}

# aggregate draws to both sexes combined
draws <- merge(draws, weights, by=c("area", "sex", "age"))
draws <- draws[, list(value = sum(value*wt), pop = sum(pop)), by='level,area,year,age,sim']
draws[, sex := 3]

# add crude and age-standardized rates
draws <- calc_all_ages(draws, std_wt)

# collapse to get estimates
est <- collapse_draws(draws, c("level", "area", "year", "sex", "age"))

# save draws and estimates
save(draws, file = paste0(temp_dir, "/draws_", level, "_", year, "_3", if (rake) "_raked" else "", ".rdata"))
save(est, file = paste0(temp_dir, "/est_", level, "_", year, "_3", if (rake) "_raked" else "", ".rdata"))

if (grepl("_c", model) & !rake) {
  ## Aggregate draws to both sexes combined and collapse to generate point estimates -----------------
  # load draws
  draws <- rbindlist(lapply(sexes, function(sex) {
    load(paste0(temp_dir, "/draws_pi_mort_", level, "_", year, "_", sex, if (rake) "_raked" else "", ".rdata"))
    draws
  }))
  
  # drop all-ages rates, we recalculate these after aggregating by sex
  draws <- draws[age < 98,]
  
  # if necessary, add in 0s for sexes we don't model
  if (length(sexes) == 1) {
    temp <- copy(draws)
    temp[, value := 0]
    if (sexes == 1) temp[, sex := 2] else temp[, sex := 1]
    draws <- rbind(draws, temp)
    setkeyv(draws, c("level", "area", "year", "sex", "age", "sim"))
    rm(temp); gc()
  }
  
  # aggregate draws to both sexes combined
  draws <- merge(draws, weights, by=c("area", "sex", "age"))
  draws <- draws[, list(value = sum(value*wt), pop = sum(pop)), by='level,area,year,age,sim']
  draws[, sex := 3]
  rm(weights)
  
  # add crude and age-standardized rates
  draws <- calc_all_ages(draws, std_wt)
  
  # collapse to get estimates
  est <- collapse_draws(draws, c("level", "area", "year", "sex", "age"))
  
  # save draws and estimates
  save(draws, file = paste0(temp_dir, "/draws_pi_mort_", level, "_", year, "_3", if (rake) "_raked" else "", ".rdata"))
  save(est, file = paste0(temp_dir, "/est_pi_mort_", level, "_", year, "_3", if (rake) "_raked" else "", ".rdata"))
}