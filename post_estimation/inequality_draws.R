# -----------------------------------------------------------------------------------------------------------
# Run inequality analysis at the draw level
#------------------------------------------------------------------------------------------------------------


# Functions and packages -------------------------------------------------------------------------

# clear memory
rm(list=ls())

#load external packages
library('pcaPP')
library('ccaPP')
library('mgsub')

# load MBG packages
core_repo <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# custom functions
source(paste0(core_repo, '/custom_functions/format_draw_objects.R'))
source(paste0(core_repo, '/custom_functions/analysis_formulas.R'))


# Arguments-------------------------------------------------------------------------

# indicate whether running interactively
interactive <- FALSE

# if running interactively, set arguments
if (interactive) {
  warning('interactive is set to TRUE - if you did not mean to run code interactively then kill the job and set interactive to FALSE')
  
  ## set arguments
  indicator_group          <- 'ort'
  indicator                <- 'had_diarrhea'
  run_date                 <- '2019_09_17_14_12_53'
  shapefile                <- '2019_09_10'
  region                   <- 'dia_sssa'
  meas                     <- 'prevalence'
  ineq_type                <- 'rel'
  modeling_shapefile_version <- shapefile
  
} else {
  
  # otherwise, grab arguments from qsub
  indicator_group          <- commandArgs()[[4]]; message(indicator_group)
  indicator                <- commandArgs()[[5]]; message(indicator)
  run_date                 <- commandArgs()[[6]]; message(run_date)
  shapefile                <- commandArgs()[[7]]; message(shapefile)
  region                   <- commandArgs()[[8]]; message(region)
  meas                     <- commandArgs()[[9]]; message(meas)
  ineq_type                <- commandArgs()[[10]]; message(ineq_type)
  modeling_shapefile_version <- shapefile
  
} 

# print out session info so we have it on record
sessionInfo()


# Setup -------------------------------------------------------------------------

PID <- Sys.getpid()
tic('Entire script') # Start master timer

## Set seed for reproducibility
message('Setting seed 98118 for reproducibility')
set.seed(98118)


# Make inequality objects ------------------------------------------------------------

# define share directory
share_dir <- '<<<< FILEPATH REDACTED >>>>'

message('beginning inequality calculations for: ', region)

tic('Make table')
dt <- format_cell_pred(ind_gp = indicator_group,
                       ind = indicator,
                       var_names = list('a'), #by default use ind name, but will pass a/b to simplify cor calc
                       rd = run_date,
                       reg = region,
                       measure = 'count', #this is for population raster
                       pop_measure = 'a0004t',
                       year_start = 2000,
                       year_end = 2017,
                       rk = TRUE,
                       rk_measure = meas,
                       shapefile_version = shapefile,
                       coastal_fix = T)

# reset key (to take calculation over year for each pixel/draw)
setkey(dt, pixel_id, draw)
toc(log = TRUE)

# get population weight and country mean
tic('Getting weighted pop and country mean')
dt[, pop := pop * area_fraction]
dt[, b := weighted.mean(a, pop, na.rm = TRUE), by = c('ADM0_CODE', 'year', 'draw')]
toc(log = TRUE)

# calculate inequality by cell ID
tic('Long inequality calculation')
if (ineq_type == 'rel') dt <- dt[, dev := rel_dev(a, b), by = key(dt)]
if (ineq_type == 'abs') dt <- dt[, dev := abs_dev(a, b), by = key(dt)]
dt[, c('a', 'b') := NULL]
toc(log = TRUE)


# Aggregate inequality objects to admins and save ------------------------------------------------------------

# create output directory
dir.create('<<<< FILEPATH REDACTED >>>>')

# loop over admins and metrics to aggregate
for (i in 0:2) {
  tic(paste0('Aggregating results to admin ', i, ' for metric ', ineq_type))
  # aggregate
  dt[, agg := weighted.mean(dev, pop, na.rm = TRUE), by = c(paste0('ADM', i, '_CODE'), 'year', 'draw')]
  # save matrix
  ad <- unique(dt[, c(paste0('ADM', i, '_CODE'), 'year', 'draw', 'agg'), with = FALSE])
  ad <- dcast(ad, ... ~ draw, value.var = 'agg')
  saveRDS(ad, paste0('<<<< FILEPATH REDACTED >>>>', region, '_', indicator,
                     ifelse(indicator == 'had_diarrhea', paste0('_', meas, '_'), '_'), 
                     ineq_type, '_dev_', 'adm', i, '_draw_matrix.RDs'))
  toc(log = TRUE)
}

# cleanup
dt[, agg := NULL]
rm(ad)


# Save cell pred objects ------------------------------------------------------------

# reshape wide
tic(paste0('Reshaping ', ineq_type, ' cell draws back to wide'))
dt <- dcast(dt, ... ~ draw, value.var = 'dev')
toc(log = TRUE)

# finish up and save
tic('Saving cell pred')
message(paste0('TESTING: Percent of NA rows per column is: ', mean(is.na(dt[, V1]))))

out_path <- paste0('<<<< FILEPATH REDACTED >>>>', region, '_', indicator,
                   ifelse(indicator == 'had_diarrhea', paste0('_raked_', meas, '_'), '_'),
                   ineq_type, '_dev_cell_draws.RDs')

message('-- finished making inequality calculations across draws. now saving as \n', out_path)

saveRDS(dt, out_path)
toc(log = TRUE)

