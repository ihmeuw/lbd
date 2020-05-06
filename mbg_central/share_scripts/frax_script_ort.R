###############################################################################
## Parallel script for post-estimation and aggregation
###############################################################################


## Setup ---------------------------------------------------------------------------------------------------

interactive <- FALSE 

if (interactive) {
  
  ## Set repo location, indicator group, and some arguments
  user <- Sys.info()['user']
  core_repo <- '<<<< FILEPATH REDACTED >>>>'
  indicator_group <- 'ort'
  indicator <- 'had_diarrhea'
  config_par <- 'config_dia_best'
  cov_par <- 'covs_ort_dia_sssa'
  holdout <- 0
  age <- 0
  run_date <- '2019_09_17_14_12_53'
  measure <- 'deaths'
  reg <- 'dia_sssa'
  
} else {
  
  ## Set repo location, indicator group, and some arguments
  user            <- commandArgs()[4]
  core_repo       <- commandArgs()[5]
  indicator_group <- commandArgs()[6]
  indicator       <- commandArgs()[7]
  config_par      <- commandArgs()[8]
  cov_par         <- commandArgs()[9]
  reg             <- commandArgs()[10]
  run_date        <- commandArgs()[13]
  measure         <- commandArgs()[14]
  holdout         <- as.numeric(commandArgs()[15])
  age             <- 0

}

## Load MBG packages
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## Read config file and save all parameters in memory
config_filepath <- '<<<< FILEPATH REDACTED >>>>'
cov_filepath <- '<<<< FILEPATH REDACTED >>>>'
config <- set_up_config(repo            = core_repo,
                        indicator_group = '',
                        indicator       = '',
                        config_name     = paste0(config_filepath, config_par),
                        covs_name       = paste0(cov_filepath, cov_par))

# Get the necessary variables out from the config object into global env
rake_countries <- eval(parse(text = config[V1 == 'rake_countries', V2]))
rake_subnational <- eval(parse(text = config[V1 == 'subnational_raking', V2]))
modeling_shapefile_version <- config[V1 == 'modeling_shapefile_version', V2]
raking_shapefile_version <- config[V1 == 'raking_shapefile_version', V2]
countries_not_to_rake <- config[V1 == 'countries_not_to_rake', V2]
countries_not_to_subnat_rake <- config[V1 == 'countries_not_to_subnat_rake', V2]
year_list <- eval(parse(text = config[V1 == 'year_list', V2]))
metric_space <- config[V1 == 'metric_space', V2]

## Set filepath and pathaddin
sharedir <- '<<<< FILEPATH REDACTED >>>>'
outputdir <- '<<<< FILEPATH REDACTED >>>>'
pathaddin <- paste0('_bin0_', reg, '_', holdout)

# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(reg)
message(pop_measure)
message(holdout)
message(measure)
message(modeling_shapefile_version)
message(raking_shapefile_version)


## Define raking parameters ---------------------------------------------------------------------------

# Make sure we're subnationally raking
subnational_raking <- TRUE
countries_not_to_subnat_rake <- 'NGA+PAK+IRN+PHL'

# Make sure we're not taking South Africa incidence or prevalence
if (measure == 'incidence' | measure == 'prevalence') {
  countries_not_to_rake <- 'ZAF'
} else {
  countries_not_to_rake <- NULL
}

# Assume subnational raking unless otherwise speciified
if (subnational_raking == FALSE) rake_subnational <- FALSE else rake_subnational <- TRUE

# Determine if a crosswalk is needed
if (modeling_shapefile_version == raking_shapefile_version) crosswalk <- F else crosswalk <- T

# Force linear raking to avoid issues with logit raking
rake_method <- 'linear'

# Set population version
pop_release <- '2019_08_29'

# Print raking info
print(paste0('Metric Space                       : ', metric_space))
print(paste0('Subnational raking                 : ', rake_subnational))
print(paste0('Countries not to rake at all       : ', countries_not_to_rake))
print(paste0('Countries not to rake subnationally: ', countries_not_to_subnat_rake))


## Load GBD estimates -------------------------------------------------------------------------------------

rake_to_path <- '<<<< FILEPATH REDACTED >>>>'
gbd <- as.data.table(read.csv(paste0(rake_to_path, 'gbd_',  measure, '.csv'), stringsAsFactors = FALSE))
setnames(gbd, c('location_id', 'year_id', 'val'), c('name', 'year', 'mean'))


## Load cell pred and populations -----------------------------------------------------------------

# Get the simple and new_simple rasters prepped up for us
message('Loading simple and prepped rasters')
raster_outputs <- prep_shapes_for_raking(
  reg = reg,
  modeling_shapefile_version = modeling_shapefile_version,
  raking_shapefile_version = raking_shapefile_version,
  field = 'loc_id'
)

## Take out the objects from the list that actually matters to us:
simple_raster <- raster_outputs[['simple_raster']]
new_simple_raster <- raster_outputs[['new_simple_raster']]

simple_polygon <- raster_outputs[['simple_polygon']]
new_simple_polygon <- raster_outputs[['new_simple_polygon']]

pixel_id <- raster_outputs[['pixel_id']]

# Load cell draws
message('Loading cell pred')
load(paste0(outputdir, indicator, '_cell_draws_eb_bin0_', reg, '_', holdout, '.RData'))

# Convert to incicence using GBD approximation before raking
if (measure == 'incidence') cell_pred <- cell_pred/(4.2/365)

# Load populations
message('Loading populations')
gbd_pops <- prep_gbd_pops_for_fraxrake(pop_measure = pop_measure, reg = reg, year_list = year_list, gbd_round_id = 5)


## Rake and aggregate -----------------------------------------------------------------

print('Using the rates raking and aggregation functions:')

## First, create all the fractional rake factors
fractional_rake_rates(
  cell_pred = cell_pred,
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
  hold = holdout,
  meas = measure
)

## Now, create the raked cell pred files!
outputs <- fractional_agg_rates(
  cell_pred = cell_pred,
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
  return_objects = TRUE,
  meas = measure
)

## Save mean
ras  <- insertRaster(new_simple_raster, matrix(rowMeans(outputs[['raked_cell_pred']]), ncol = 18))
writeRaster(
  ras,
  file = (paste0(outputdir, '/', indicator, '_prediction_', measure, '_eb',pathaddin)),
  overwrite = TRUE
)
rm(ras)

## Save lower
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, lower), ncol = 18))
writeRaster(
  ras,
  file = paste0(outputdir, '/', indicator,'_lower_', measure, '_eb', pathaddin),
  overwrite = TRUE
)
rm(ras)

## Save upper
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, upper), ncol = 18))
writeRaster(
  ras,
  file = paste0(outputdir, '/', indicator,'_upper_', measure, '_eb', pathaddin),
  overwrite = TRUE
)
rm(ras)


## Finish up -----------------------------------------------------------------

## Write a file to mark done
write(NULL, file = paste0(outputdir, '/fin_', pathaddin))

## All done
message(paste0('Done with post-estimation and aggregation for ', reg))