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
  indicator <- 'ors'
  config_par <- 'ort_best'
  cov_par <- 'covs_ors/covs_ort_dia_sssa'
  holdout <- 0
  age <- 0
  run_date <- '2019_11_07_13_50_11'
  measure <- 'prevalence'
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
config_filepath <- 'ors/3_modeling/'
config <- set_up_config(repo            = core_repo,
                        indicator_group = '',
                        indicator       = '',
                        config_name     = paste0(config_filepath, 'config_', config_par),
                        covs_name       = paste0(config_filepath, cov_par))

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

# Create function to pull isos by region
get_region_isos <- function(modeling_region) {
  
  # define regions
  region_list <- list(
    'dia_afr_horn' = 'dji+eri+eth+sdn+som+ssd+yem',
    'dia_cssa-cod' = 'ago+caf+cog+gab+gnq+stp',
    'dia_wssa-nga' = 'ben+bfa+civ+cmr+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+ner+sen+sle+tcd+tgo',
    'dia_name' = 'dza+egy+lby+mar+tun',
    'dia_sssa' = 'bwa+nam+zaf',
    'dia_mcaca' = 'blz+cri+cub+dma+dom+grd+gtm+hnd+hti+jam+lca+mex+nic+pan+slv+vct',
    'dia_s_america-bra' = 'bol+col+ecu+guy+per+pry+sur+tto+ven',
    'dia_central_asia' = 'kgz+tjk+tkm+uzb',
    'dia_se_asia' = 'khm+lao+mmr+mys+tha+vnm',
    'dia_malay' = 'idn+phl+png+tls',
    'dia_south_asia-ind-pak' = 'bgd+btn+lka+npl',
    'dia_mid_east' = 'afg+irn+irq+jor+pse+syr',
    'dia_essa-zwe-ken' = 'bdi+com+lso+mdg+moz+mwi+rwa+swz+syc+tza+uga+zmb',
    'PAK' = 'pak',
    'KEN' = 'ken',
    'NGA' = 'nga',
    'COD' = 'cod',
    'IND' = 'ind',
    'ZWE' = 'zwe',
    'MNG' = 'mng',
    'dia_cssa' = 'ago+caf+cod+cog+gab+gnq+stp',
    'dia_wssa' = 'ben+bfa+civ+cmr+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+ner+nga+sen+sle+tcd+tgo',
    'dia_s_america' = 'bol+bra+col+ecu+guy+per+pry+sur+tto+ven',
    'dia_chn_mng' = 'chn+mng',
    'dia_south_asia' = 'bgd+btn+ind+lka+npl+pak',
    'dia_south_asia-ind' = 'bgd+btn+lka+npl+pak',
    'dia_essa' = 'bdi+com+ken+lso+mdg+moz+mwi+rwa+swz+syc+tza+uga+zmb+zwe'
  )
  # return region list
  return(toupper(region_list[[modeling_region]]))
}

# apply function
countries_not_to_rake <- get_region_isos(reg)
countries_not_to_subnat_rake <- get_region_isos(reg)
subnational_raking <- FALSE

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

## Delete raked files since we don't rake ORS estimates to GBD
to_del <- list.files(outputdir, pattern = paste0('_raked_', measure, '_eb_bin0_', reg, '_', holdout))
to_del <- c(to_del, list.files(outputdir, pattern = paste0('_raked_', measure, '_admin_draws_eb_bin0_', reg, '_', holdout)))
to_del <- c(to_del, list.files(outputdir, pattern = paste0('_raked_', measure, '_c_admin_draws_eb_bin0_', reg, '_', holdout)))
if (length(to_del) > 0) lapply(paste0(outputdir, to_del), file.remove)


## Finish up -----------------------------------------------------------------

## Write a file to mark done
write(NULL, file = paste0(outputdir, '/fin_', pathaddin))

## All done
message(paste0('Done with post-estimation and aggregation for ', reg))