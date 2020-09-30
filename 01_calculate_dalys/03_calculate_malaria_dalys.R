#####################################################################
# (1) Sum malaria aaked YLDs and YLLs to obtain unraked DALYs
# (2) Fractional raked and aggregate malaria DALYs
#####################################################################

# (1) Setup ---------------------------------------------------------

#read in arguments from qsub
run_date <- commandArgs()[4] 
shapefile_version <- commandArgs()[5]
year_list <- unlist(strsplit(commandArgs()[6], '~'))
pop_release <- commandArgs()[7]

#fractional raking arguments
rake_countries <- TRUE
rake_subnational <- TRUE
modeling_shapefile_version <- shapefile_version
raking_shapefile_version <- shapefile_version
countries_not_to_rake <- 'ESH'
countries_not_to_subnat_rake <- 'NGA'
year_list <- as.numeric(year_list)
metric_space <- 'rates'
indicator <- 'had_malaria'
indicator_group <- 'malaria'
reg <- 'africa-ESH'
holdout <- 0
pop_measure <- 'a0004t'
interval_mo <- 12
sex_id <- 3
age <- 0
etiology <- 'none'

#load in functions and packages
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/setup.R')
commondir      <- '<<< FILEPATH REDACTED >>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list = package_list, repos = '<<< FILEPATH REDACTED >>>')

library(mgcv)
source('<<< FILEPATH REDACTED >>>/lbd_core_custom/fractional_raking_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/shapefile_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/prep_functions.R')
source('<<< FILEPATH REDACTED >>>/get_outputs.R')

# Set filepath and pathaddin
sharedir <- '<<< FILEPATH REDACTED >>>'
outputdir <- '<<< FILEPATH REDACTED >>>'
pathaddin <- '<<< FILEPATH REDACTED >>>'

#or
load('<<< FILEPATH REDACTED >>>')
ylls <- raked_cell_pred
load('<<< FILEPATH REDACTED >>>')
ylds <- raked_cell_pred

dalys <- ylls + ylds
save(dalys, file = '<<< FILEPATH REDACTED >>>')
rm(ylls, ylds)

# (3) Pull relevant outputs from GBD ---------------------------------

#set cause_id
cause_id <- 345
measure_id <- 2 #DALYs

#subset location_ids to national and subnational estimates
locs <- read.csv('<<< FILEPATH REDACTED >>>')
country <- subset(locs, level >=3)$location_id

#pull ylds
gbd <- get_outputs(topic= 'cause',
                   cause_id=cause_id,
                   location_id=country, 
                   year_id=year_list,
                   age_group_id=1,
                   gbd_round_id=5,
                   metric_id=3,
                   measure_id = measure_id,
                   version='latest') %>%
  as.data.frame() %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017)))))

gbd <- gbd[complete.cases(gbd),]
gbd$val <- ifelse(is.na(gbd$val), gbd$inter_val, gbd$val)

# Prep GBD table 
gbd <- gbd %>% ungroup() %>%
  dplyr::select(location_id, year_id, val) %>%
  dplyr::rename(name = location_id, year = year_id, mean = val)

# (4) Rake  ----------------------------------------------------

# print message
print(paste0('Raking malaria to DALYs'))

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
if (subnational_raking == FALSE) rake_subnational <- FALSE else rake_subnational <- TRUE

# Determine if a crosswalk is needed
if (modeling_shapefile_version == raking_shapefile_version) crosswalk <- F else crosswalk <- T

# Force linear raking to avoid issues with logit raking
rake_method <- 'linear'

# Print raking info
print(paste0('Metric Space                       : ', metric_space))
print(paste0('Subnational raking                 : ', rake_subnational))
print(paste0('Countries not to rake at all       : ', countries_not_to_rake))
print(paste0('Countries not to rake subnationally: ', countries_not_to_subnat_rake))

## Load cell pred and populations -----------------------------------------------------------------

# Prep simple and new simple rasters
message('Loading simple and prepped rasters')
raster_outputs <- prep_shapes_for_raking(reg = reg,
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
cell_pred <- dalys
rm(dalys)
gc()

# Load populations
message('Loading populations')
gbd_pops <- prep_gbd_pops_for_fraxrake(pop_measure = pop_measure, reg = reg, year_list = year_list, gbd_round_id = 5)

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
                      hold = holdout,
                      meas = 'daly',
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
                                return_objects = TRUE,
                                meas = 'daly',
                                etiology = 'none')


## Save mean, upper, lower rasters -------------------------------------------------------------------------------------

## Save mean
ras  <- insertRaster(new_simple_raster, matrix(rowMeans(outputs[['raked_cell_pred']]), ncol = 18))
writeRaster(
  ras,
  file = (paste0(outputdir, '/', indicator, if(!(etiology == 'none')) '_', if(!(etiology == 'none')) etiology, '_prediction_daly_eb',pathaddin)),
  overwrite = TRUE
)
rm(ras)

## Save lower
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, lower), ncol = 18))
writeRaster(
  ras,
  file = paste0(outputdir, '/', indicator,  if(!(etiology == 'none')) '_', if(!(etiology == 'none')) etiology, '_lower_daly_eb', pathaddin),
  overwrite = TRUE
)
rm(ras)

## Save upper
ras  <- insertRaster(new_simple_raster, matrix(apply(outputs[['raked_cell_pred']], 1, upper), ncol = 18))
writeRaster(
  ras,
  file = paste0(outputdir, '/', indicator, if(!(etiology == 'none')) '_', if(!(etiology == 'none')) etiology, '_upper_daly_eb', pathaddin),
  overwrite = TRUE
)
rm(ras)
