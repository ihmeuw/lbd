# -----------------------------------------------------------------------------------------------------------
# Perform ORS correlation analyses at the draw level
#------------------------------------------------------------------------------------------------------------


# Functions and packages -------------------------------------------------------------------------

# clear memory
rm(list=ls())

#load external packages
library('pcaPP', lib.loc = '<<<< FILEPATH REDACTED >>>>')
library('ccaPP', lib.loc = '<<<< FILEPATH REDACTED >>>>')
library('mgsub', lib.loc = '<<<< FILEPATH REDACTED >>>>')

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
  indicator_group1         <- 'ort'
  indicator_group2         <- 'ort'
  indicator1               <- 'ors'
  indicator2               <- 'rhf'
  run_date1                <- '2019_11_07_13_50_11'
  run_date2                <- '2019_11_07_13_50_11'
  shapefile                <- '2019_029_01'
  region                   <- 'dis_sssa'
  modeling_shapefile_version <- shapefile
  
} else {
  
  # otherwise, grab arguments from qsub
  indicator_group1         <- commandArgs()[[4]]
  indicator_group2         <- commandArgs()[[5]]
  indicator1               <- commandArgs()[[6]]
  indicator2               <- commandArgs()[[7]]
  run_date1                <- commandArgs()[[8]]
  run_date2                <- commandArgs()[[9]]
  shapefile                <- commandArgs()[[10]]
  region                   <- commandArgs()[[11]]
  modeling_shapefile_version <- shapefile

} 

# Set population version
pop_release <- '2019_08_29'

# print out session info so we have it on record
sessionInfo()


# Setup -------------------------------------------------------------------------

PID <- Sys.getpid()
tic("Entire script") # Start master timer

# Set seed for reproducibility
message('Setting seed 98118 for reproducibility')
set.seed(98118)

# Pull indicators and run dates
indicator_groups <- list(indicator_group1, indicator_group2)
indicators <- list(indicator1, indicator2)
run_dates <- list(run_date1, run_date2)


# Make correlation objects ------------------------------------------------------------

# define share directory
share_dir <- '<<<< FILEPATH REDACTED >>>>'
  
message('beginning correlation calculations for: ', region)

tic('Make table')
dt <- format_cell_pred(ind_gp = indicator_groups,
                       ind = indicators,
                       var_names = list('a', 'b'), #by default use ind name, but will pass a/b to simplify cor calc
                       rd = run_dates,
                       reg = region,
                       measure = 'count', #this is for population raster
                       pop_measure = 'a0004t',
                       year_start = 2000,
                       year_end = 2017,
                       rk = FALSE,
                       rk_measure = meas,
                       shapefile_version = shapefile,
                       coastal_fix = T)

# reset key (to take correlation over year for each pixel/draw)
setkey(dt, pixel_id, draw)
toc(log = TRUE)

# calculate spearmans correlation over years by cell ID
# note corSpearman from ccaPP, it is vectorized so should be faster and more suitable for DT
tic('Long calculation')
dt <- dt[, cor := corSpearman(a, b), by=key(dt)]
dt[, `:=`(a=NULL, b=NULL)] #remove indicators to save memory, no longer needed
toc(log = TRUE)

# get population weighted
tic('Getting weighted pop and country mean')
dt[, pop := pop * area_fraction]
toc(log = TRUE)


# Aggregate inequality objects to admins and save ------------------------------------------------------------

# create output directory
dir.create('<<<< FILEPATH REDACTED >>>>')
dir.create('<<<< FILEPATH REDACTED >>>>')

# loop over admins and metrics to aggregate
for (i in 0:2) {
  tic(paste0('Aggregating results to admin ', i))
  # aggregate
  dt[, agg := weighted.mean(cor, pop, na.rm = TRUE), by = c(paste0('ADM', i, '_CODE'), 'year', 'draw')]
  # save matrix
  ad <- unique(dt[, c(paste0('ADM', i, '_CODE'), 'year', 'draw', 'agg'), with = FALSE])
  ad <- dcast(ad, ... ~ draw, value.var = 'agg')
  saveRDS(ad, paste0('<<<< FILEPATH REDACTED >>>>', region, '_', indicators[[1]], '_', indicators[[2]],
                     '_corr_', 'adm', i, '_draw_matrix.RDs'))
  toc(log = TRUE)
}

# cleanup
dt[, agg := NULL]
rm(ad)


# Save cell pred objects ------------------------------------------------------------
# reshape wide
tic('Reshaping back to wide')
dt <- dcast(dt, ... ~ draw, value.var= 'cor')
toc(log = TRUE)

# finish up and save
tic('Saving cell pred')
message(sprintf('TESTING: Percent of NA rows per column is: %f%%', mean(is.na(dt[, V1]))))

out.path <- sprintf('<<<< FILEPATH REDACTED >>>>/%s_%s_%s_corr_cell_draws.rds',
                    share_dir, region, indicators[[1]], indicators[[2]])

dir.create('<<<< FILEPATH REDACTED >>>>')

message('-- finished making correlation across draws. now saving as \n', out.path)

saveRDS(dt, out.path)
toc(log = TRUE)
