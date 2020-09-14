######################################################################
# Fractional raking script for LRI
#####################################################################

# (1) Setup -----------------------------------------------------------------------------------------------------------------

# Pull in commandArgs
user            <- commandArgs()[4]
core_repo       <- commandArgs()[5]
indicator_group <- commandArgs()[6]
indicator       <- commandArgs()[7]
config_par      <- commandArgs()[8]
cov_par         <- commandArgs()[9]
reg             <- commandArgs()[10]
run_date        <- commandArgs()[11]
measure         <- commandArgs()[12]
holdout         <- as.numeric(commandArgs()[13])
etiology        <- commandArgs()[14]
age             <- 0

print(commandArgs())

# MBG packages
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)
library(mgcv)

# source any functions with custom edits from lbd_custom folder
source('<<<< FILEPATH REDACTED >>>>')
source('<<<< FILEPATH REDACTED >>>>')
source('<<<< FILEPATH REDACTED >>>>')

# Read config file and save all parameters in memory
config <- set_up_config(repo = '<<<< FILEPATH REDACTED >>>>',
                        core_repo = core_repo,
                        indicator_group,
                        indicator,
                        config_name = paste0('config_', config_par),
                        covs_name = paste0('covs_', cov_par),
                        run_tests = FALSE)

# Get the necessary variables out from the config object into global env
rake_countries <- eval(parse(text = config[V1 == 'rake_countries', V2]))
rake_subnational <- eval(parse(text = config[V1 == 'subnational_raking', V2]))
modeling_shapefile_version <- config[V1 == 'modeling_shapefile_version', V2]
raking_shapefile_version <- config[V1 == 'raking_shapefile_version', V2]
countries_not_to_rake <- config[V1 == 'countries_not_to_rake', V2]
countries_not_to_subnat_rake <- config[V1 == 'countries_not_to_subnat_rake', V2]
year_list <- eval(parse(text = config[V1 == 'year_list', V2]))
metric_space <- config[V1 == 'metric_space', V2]

# Set filepath and pathaddin
sharedir <- '<<<< FILEPATH REDACTED >>>>'
outputdir <- '<<<< FILEPATH REDACTED >>>>'
pathaddin <- '<<<< FILEPATH REDACTED >>>>'

# Print select settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(reg)
message(pop_measure)
message(holdout)
message(paste0('pop_release: ', pop_release))
if (!(etiology == 'none')) message('etiology = ', etiology) else message('not using specific etiology split')

# (2) Define raking parameters ----------------------------------------------------------------

# turn on subnational raking
subnational_raking <- TRUE
rake_subnational <- TRUE
countries_not_to_subnat_rake <- c('NGA','PAK','PHL')

# countries not to rake
countries_not_to_rake <- NULL

# Determine if a crosswalk is needed
if (modeling_shapefile_version == raking_shapefile_version) crosswalk <- F else crosswalk <- T

# Force linear raking to avoid issues with logit raking
rake_method <- 'linear'

# Print raking info
print(paste0('Metric Space                       : ', metric_space))
print(paste0('Subnational raking                 : ', rake_subnational))
print(paste0('Countries not to rake at all       : ', countries_not_to_rake))
print(paste0('Countries not to rake subnationally: ', countries_not_to_subnat_rake))

# (3) Load cell pred and populations -----------------------------------------------------------------

# Prep simple and new simple rasters
message('Loading simple and prepped rasters')
raster_outputs <- prep_shapes_for_raking(reg = reg,
                                         modeling_shapefile_version = modeling_shapefile_version,
                                         raking_shapefile_version = raking_shapefile_version,
                                         field = 'loc_id',
                                         pop_release = pop_release,
                                         start_year = min(year_list),
                                         end_year = max(year_list))

# Extract needed items from list
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

# (4) Load GBD estimates -------------------------------------------------------------------------------------

# Set measure_id for get_outputs
if (measure == 'incidence'){
  measure_id <- 6
} else if (measure == 'prevalence'){
  measure_id <- 5
} else if (measure == 'mortality'){
  measure_id <- 1
} else if (measure == 'daly'){
  measure_id <- 2
}

# Subset locations to national and subnational estimates
locs <- read.csv('<<<< FILEPATH REDACTED >>>>')
country <- subset(locs, level >=3)$location_id

gbd <- get_outputs(topic= 'cause',
                   cause_id=322,
                   location_id=country,
                   year_id=year_list,
                   age_group_id=1,
                   gbd_round_id=6,
                   metric_id=3,
                   measure_id = measure_id,
                   sex_id = 3,
                   version = 'latest',
                   decomp_step = 'step5') %>%
  
  as.data.frame() %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)), newdata = data.frame(year_id = year_list)))))

gbd$val <- ifelse(is.na(gbd$val), gbd$inter_val, gbd$val)
gbd <- gbd[complete.cases(gbd),]

setnames(gbd, c('location_id', 'year_id', 'val'), c('name', 'year', 'mean'))

# if using etiology, pull GBD PAFs and multiply total burder by this fraction:

if (!(etiology == 'none')){
  message(paste0('using etiology-specific GBD estimates for cause = ', etiology))
  
  # set et_cause_id
  if (etiology == 'lri_flu') rei_id <- 187
  if (etiology == 'lri_hib') rei_id <- 189
  if (etiology == 'lri_pneumo')rei_id <- 188
  if (etiology == 'lri_rsv') rei_id <- 190
  
  # set YLD (incidence) or YLL (mortality) measure IDs
  if (measure == 'mortality') et_measure_id <- 4 #YLL
  if (measure == 'incidence') et_measure_id <- 3 #YLD
  
  # pull GBD YLD or YLL PAFs
  gbd_paf <- get_outputs(topic= 'rei',
                         rei_id= rei_id,
                         location_id = country,
                         cause_id=322,
                         year_id=year_list,
                         age_group_id=1,
                         gbd_round_id=6,
                         metric_id=2,
                         measure_id = et_measure_id,
                         sex_id = 3,
                         version='latest',
                         decomp_step = 'step5') %>%
    
    as.data.frame() %>%
    group_by(location_id) %>%
    mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)), newdata = data.frame(year_id = year_list)))))
  
  gbd_paf$val <- ifelse(is.na(gbd_paf$val), gbd_paf$inter_val, gbd_paf$val)
  gbd_paf <- gbd_paf[complete.cases(gbd_paf),]
  
  setnames(gbd_paf, c('location_id', 'year_id', 'val'), c('name', 'year', 'paf'))
  
  gbd_paf <- as.data.table(gbd_paf) %>%
    dplyr::select(year, name, paf)
  
  # multiply GBD values by GBD PAFs
  gbd <- merge(gbd, gbd_paf, by = c('year','name')) %>%
    as.data.table
  gbd[,mean := paf*mean]
}

# (5) Rake and aggregate -------------------------------------------------------------------------------------

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
                      hold = holdout,
                      meas = measure,
                      etiology = etiology)

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
                                return_objects = TRUE,
                                meas = measure,
                                etiology = etiology)


# (6) Save mean, upper, lower rasters -------------------------------------------------------------------------------------

# Save mean
ras  <- insertRaster(new_simple_raster, matrix(rowMeans(outputs[['raked_cell_pred']]), ncol = length(year_list)))
writeRaster(
  ras,
  file = (paste0(outputdir, '/', indicator, if(!(etiology == 'none')) '_', if(!(etiology == 'none')) etiology, '_prediction_', measure, '_eb',pathaddin)),
  overwrite = TRUE
)
rm(ras)

# Save lower
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, lower), ncol = length(year_list)))
writeRaster(
  ras,
  file = paste0(outputdir, '/', indicator,  if(!(etiology == 'none')) '_', if(!(etiology == 'none')) etiology, '_lower_', measure, '_eb', pathaddin),
  overwrite = TRUE
)
rm(ras)

# Save upper
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, upper), ncol = length(year_list)))
writeRaster(
  ras,
  file = paste0(outputdir, '/', indicator, if(!(etiology == 'none')) '_', if(!(etiology == 'none')) etiology, '_upper_', measure, '_eb', pathaddin),
  overwrite = TRUE
)
rm(ras)
