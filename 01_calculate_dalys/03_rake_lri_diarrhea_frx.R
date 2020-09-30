## Setup ----------------------------------------------------------------
interactive <- FALSE 

if (interactive) {
  
  ## Set repo location, indicator group, and some arguments
  user <- '<<< USERNAME REDACTED >>>'
  core_repo <- '<<< FILEPATH REDACTED >>>'
  indicator_group <- 'lri'
  indicator <- 'had_lri'
  holdout <- 0
  age <- 0
  run_date <- '2019_07_18_14_28_10'
  measure <- 'daly'
  reg <- 'name'
  modeling_shapefile_version <- '2019_02_27'
  raking_shapefile_version <- '2019_02_27'
  year_list <- c(2000:2017)
  pop_release <- '2019_08_29'
  
} else {
  
  ## Set repo location, indicator group, and some arguments
  user                       <- commandArgs()[4]
  core_repo                  <- commandArgs()[5]
  indicator_group            <- commandArgs()[6]
  indicator                  <- commandArgs()[7]
  reg                        <- commandArgs()[8]
  run_date                   <- commandArgs()[9]
  measure                    <- commandArgs()[10]
  holdout                    <- as.numeric(commandArgs()[11])
  modeling_shapefile_version <- commandArgs()[12]
  raking_shapefile_version   <- commandArgs()[13]
  year_list                  <- as.numeric(unlist(strsplit(commandArgs()[14], '~')))
  dalys_run_date             <- commandArgs()[15]
  pop_release                <- commandArgs()[16]
  age                        <- 0
  
  
}

## Load MBG packages
package_list <- c(t(read.csv('<<< FILEPATH REDACTED >>>',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)
library(mgcv)

## source any functions with custom edits from lbd_custom folder
source('<<< FILEPATH REDACTED >>>/lbd_core_custom/fractional_raking_functions.R')
source('<<< FILEPATH REDACTED >>>/get_outputs.R')

#set config args manually rather than reading from config
rake_countries <- TRUE
rake_subnational <- TRUE
if (indicator == 'had_diarrhea') {
  countries_not_to_rake = NULL
} else {
  countries_not_to_rake <- 'ESH+GUF'
}
countries_not_to_subnat_rake <- NULL
metric_space <- 'rates'
pop_measure <- 'a0004t'
interval_mo <- '12'

## Set filepath and pathaddin
sharedir <- '<<< FILEPATH REDACTED >>>'
outputdir <- '<<< FILEPATH REDACTED >>>'
pathaddin <- '<<< FILEPATH REDACTED >>>'
dir.create('<<< FILEPATH REDACTED >>>')

# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(reg)
message(pop_measure)
message(holdout)

## Define raking parameters ----------------------------------------------------------------

#turn on subnational raking
subnational_raking <- TRUE
countries_not_to_subnat_rake <- NULL
if (subnational_raking == FALSE) rake_subnational <- FALSE else rake_subnational <- TRUE

# Determine if a crosswalk is needed
if (modeling_shapefile_version == raking_shapefile_version) crosswalk <- F else crosswalk <- T

# Force linear raking to avoid issues with logit raking
rake_method <- 'linear'

# Force countries not to subnat rake
countries_not_to_subnat_rake <- c('NGA','PAK','PHL')

# Print raking info
print(paste0('Metric Space                       : ', metric_space))
print(paste0('Subnational raking                 : ', rake_subnational))
print(paste0('Countries not to rake at all       : ', countries_not_to_rake))
print(paste0('Countries not to rake subnationally: ', countries_not_to_subnat_rake))

## Load cell pred and populations -----------------------------------------------------------------

# Prep simple and new simple rasters
message('Loading simple and prepped rasters')

#diqrrhea name fix
if (reg == 'dia_name') {
  reg_load <- 'dia_name-ESH'
} else {
  reg_load <- reg
}

raster_outputs <- prep_shapes_for_raking(reg = reg_load,
                                         modeling_shapefile_version = modeling_shapefile_version,
                                         raking_shapefile_version = raking_shapefile_version,
                                         field = 'loc_id')

## Extract needed items from list
simple_raster <- raster_outputs[['simple_raster']]
new_simple_raster <- raster_outputs[['new_simple_raster']]

simple_polygon <- raster_outputs[['simple_polygon']]
new_simple_polygon <- raster_outputs[['new_simple_polygon']]

pixel_id <- raster_outputs[['pixel_id']]

# Load cell draws
message('Loading cell pred')
load(paste0(outputdir, indicator, '_cell_draws_eb_bin0_', reg, '_', holdout, '.RData'))

# Load populations
message('Loading populations')
gbd_pops <- prep_gbd_pops_for_fraxrake(pop_measure = pop_measure, reg = reg, year_list = year_list, gbd_round_id = 5)

## Load GBD estimates -------------------------------------------------------------------------------------

#Set measure_id for get_outputs
if (measure == 'incidence'){
  measure_id <- 6
} else if (measure == 'prevalence'){
  measure_id <- 5
} else if (measure == 'mortality'){
  measure_id <- 1
} else if (measure == 'daly'){
  measure_id <- 2
}else if (measure == 'yll'){
  measure_id <- 4
}else if (measure == 'yld'){
  measure_id <- 3
}

if (indicator == 'has_lri') cause_id <- 322
if (indicator == 'had_diarrhea') cause_id <- 302

# Subset locations to National and Subnational estimates
locs <- read.csv('<<< FILEPATH REDACTED >>>')
country <- subset(locs, level >=3)$location_id

gbd <- get_outputs(topic= 'cause',
                   cause_id=cause_id,
                   location_id=country,
                   year_id=year_list,
                   age_group_id=1,
                   gbd_round_id=5,
                   metric_id=3,
                   measure_id = measure_id,
                   sex_id = 3,
                   version='latest')%>%
  
  as.data.frame() %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)), newdata = data.frame(year_id = year_list)))))

gbd$val <- ifelse(is.na(gbd$val), gbd$inter_val, gbd$val)
gbd <- gbd[complete.cases(gbd),]

setnames(gbd, c('location_id', 'year_id', 'val'), c('name', 'year', 'mean'))

## Rake and aggregate -------------------------------------------------------------------------------------

print('Using the rates raking and aggregation functions:')

## First, create all the fractional rake factors
fractional_rake_rates(cell_pred = cell_pred,
                      simple_raster = simple_raster,
                      simple_polygon = simple_polygon,
                      pixel_id = pixel_id,
                      shapefile_version = raking_shapefile_version,
                      reg = reg,
                      pop_measure = pop_measure,
                      year_list = year_list,
                      interval_mo = interval_mo,
                      rake_subnational = rake_subnational,
                      age_group = age_group,
                      sex_id = sex_id,
                      sharedir = sharedir,
                      run_date = run_date,
                      indicator = indicator,
                      gbd = gbd,
                      rake_method = rake_method,
                      gbd_pops = gbd_pops,
                      countries_not_to_rake = countries_not_to_rake,
                      countries_not_to_subnat_rake = countries_not_to_subnat_rake,
                      custom_output_folder = '<<< FILEPATH REDACTED >>>',
                      hold = holdout,
                      meas = measure,
                      etiology = 'none')

## Now, create the raked cell pred files!
outputs <- fractional_agg_rates(cell_pred = cell_pred,
                                simple_raster = simple_raster,
                                simple_polygon = simple_polygon,
                                pixel_id = pixel_id,
                                shapefile_version = raking_shapefile_version,
                                reg = reg,
                                pop_measure = pop_measure,
                                year_list = year_list,
                                interval_mo = interval_mo,
                                rake_subnational = rake_subnational,
                                sharedir = sharedir,
                                run_date = run_date,
                                indicator = indicator,
                                main_dir = outputdir,
                                rake_method = rake_method,
                                age = age,
                                holdout = holdout,
                                countries_not_to_subnat_rake = countries_not_to_subnat_rake,
                                custom_output_folder = '<<< FILEPATH REDACTED >>>',
                                return_objects = TRUE,
                                meas = measure,
                                etiology = 'none')


## Save mean, upper, lower rasters -------------------------------------------------------------------------------------

## Save mean
ras  <- insertRaster(new_simple_raster, matrix(rowMeans(outputs[['raked_cell_pred']]), ncol = 18))
writeRaster(
  ras,
  file = ('<<< FILEPATH REDACTED >>>'),
  overwrite = TRUE
)
rm(ras)

## Save lower
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, lower), ncol = 18))
writeRaster(
  ras,
  file = ('<<< FILEPATH REDACTED >>>'),
  overwrite = TRUE
)
rm(ras)

## Save upper
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, upper), ncol = 18))
writeRaster(
  ras,
  file = ('<<< FILEPATH REDACTED >>>'),
  overwrite = TRUE
)
rm(ras)
