##############################################################################
## Script to split cell pred draw-level MBG outputs by age and/or sex groups
##############################################################################


## Setup -------------------------------------------------------------------------

# Clear environment
rm(list = ls())

# Get MBG arguments
repo            <- commandArgs()[4]
indicator_group <- commandArgs()[5]
indicator       <- commandArgs()[6]
raked           <- commandArgs()[7]
measure         <- commandArgs()[8]
run_date        <- commandArgs()[9]
region          <- commandArgs()[10]
modeling_shapefile_version <- commandArgs()[11]
split_by        <- commandArgs()[12]

# Set GBD arguments
cause_id        <- 322
year_id         <- 2000:2017
metric_id       <- c(1,3)
gbd_round_id    <- 5

if (split_by == 'age'){
  age_group_id    <- c(2,3,4,5) 
  sex_id          <- 3
} else if(split_by == 'sex'){
  age_group_id    <- 1 
  sex_id          <- c(1,2)
}

# Load MBG packages and functions
package_list <- c(t(read.csv('<<< FILEPATH REDACTED >>>', header=FALSE)))
source(paste0(repo, 'lbd_core/mbg_central/setup.R'))
mbg_setup(package_list = package_list, paste0(repo, 'lbd_core'))

# Load custom functions
source(paste0(repo, 'lri/y4_age_sex_splitting/functions/format_draw_objects.R'))
source(paste0(repo, 'lri/y4_age_sex_splitting/functions/format_gbd_results.R'))

# Define share directory
share_dir <- '<<< FILEPATH REDACTED >>>'

# Create output directory
outdir <- paste0(share_dir, 'age_sex_splits/')
dir.create(outdir, showWarnings = FALSE)

# Get GBD measure_id based on measure
measure_id <- get_gbd_measure_id(measure)

# Report arguments
message('\nStarting age and/or sex splitting at the cell pred level for:')
message('indicator_group = ', indicator_group)
message('indicator = ', indicator)
message('raked = ', raked)
message('measure = ', measure)
message('run_date = ', run_date)
message('modeling_shapefile_version = ', modeling_shapefile_version)
message('region = ', paste0(region, '\n'))


## Get GBD scaling factors ------------------------------------------------------------------

# Get scaling factors
message('Getting GBD scaling factors')
results <- get_gbd_group_sf(region = region,
                            modeling_shapefile_version = modeling_shapefile_version,
                            year_id = year_id,
                            cause_id = cause_id,
                            age_group_id = age_group_id,
                            sex_id = sex_id,
                            measure_id = measure_id,
                            gbd_round_id = gbd_round_id)


## Pull cell draws -----------------------------------------------------------------

# format cell draw object
message('Formatting cell pred object')
dt <- format_cell_pred(ind_gp = indicator_group,
                       ind = indicator,
                       var_names = list('a'), #by default uses ind name, but pass 'a' to make generalizable
                       rd = run_date,
                       reg = region,
                       measure = 'count', #this is for population raster
                       pop_measure = 'a0004t',
                       year_start = year_id[1],
                       year_end = year_id[length(year_id)],
                       rk = raked,
                       rk_measure = measure,
                       shapefile_version = modeling_shapefile_version,
                       coastal_fix = T)

# remove unneeded column to free up space
dt[, 'ADM1_CODE' := NULL]

# save long cell draw object so we don't have to store it in memory
message('Saving temporary formatted cell pred file')
saveRDS(dt, file = paste0(outdir, indicator, '_', 
                          ifelse(raked, measure, 'unraked'), 
                          '_cell_pred_long_', region, '.rds'))
rm(dt)


## Split cell draws ------------------------------------------------------------

# loop over groups
for (g in unique(results$group)) {
  message('Performing splits for group: ', g)
  
  # get outputs to work with
  dt <- as.data.table(readRDS(paste0(outdir, indicator, '_', 
                                     ifelse(raked, measure, 'unraked'),  
                                     '_cell_pred_long_', region, '.rds')))
  results_by_group <- copy(results)
  
  # subset by group and free up space
  results_by_group <- results_by_group[group == g]
  results_by_group[, group := NULL]
  
  # add scaling factors onto cell pred and free up space
  dt <- merge(dt, results_by_group, by = c('ADM0_CODE', 'year'), allow.cartesian = T)
  dt[, ADM0_CODE := NULL]
  
  # multiply to get split estimate and free up space
  message('~>multiplying by scaling factor and area fraction')
  dt[, split := lapply(.SD, function(x, y) {x * y}, y = scaling_factor),
     .SDcols = 'a', by = c('year', 'draw', 'pixel_id', 'ADM2_CODE')]
  dt[, c('a', 'scaling_factor') := NULL]
  
  # get weighted mean by pixel id and clean up
  message('~~>multiplying by area fraction to get a single value by pixel')
  dt[, split := lapply(.SD, function(x, y) {x * y}, y = area_fraction),
     .SDcols = 'split', by = c('year', 'draw', 'pixel_id')]
  dt[, c('ADM0_CODE', 'ADM1_CODE', 'ADM2_CODE', 'pop', 'area_fraction') := NULL]
  dt <- unique(dt)
  
  # clean up and reshape wide
  message('~~~>reshaping wide')
  dt <- dcast(dt, ... ~ draw, value.var = 'split')
  
  # check for NAs
  message('TESTING: Percent of NA rows per column is: ', mean(is.na(dt[, V1])), '%')
  
  # save matrix object
  outpath <- paste0(outdir, indicator, '_', 
                    ifelse(raked, measure, 'unraked'),
                    '_cell_draws_', g, '_', region, '.rds')
  message('Finished making splits across cell draws for group: ', g, '. Now saving at \n', outpath)
  saveRDS(dt, outpath)
  rm(dt)
  
} # End loop over groups

# delete long cell draw object as it's no longer needed
message('Deleting temporary formatted cell pred file')
unlink(paste0(outdir, indicator, '_cell_pred_long_', region, '.rds'))

# done!
message('Finished all age and/or sex group splits for ', region)