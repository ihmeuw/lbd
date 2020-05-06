# -----------------------------------------------------------------------------------------------------------
# Calculate childhood diarrheal deaths averted at the draw level for ORS
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


# Arguments ------------------------------------------------------------------------------------------------------------
  
## otherwise, grab arguments from qsub
indicator_groups         <- list('ort', 'ort')
indicators               <- list('ors', 'had_diarrhea')
run_dates                <- list('2019_10_08_Lancet_RTR1', '2019_09_17_14_12_53')
measures                 <- list('', '_deaths')
suffixes                 <- list('_eb_bin0_0', '_eb_bin0_0')
rks                      <- list(F, T)
shapefile                <- '2019_09_10'
ors_rr                   <- 0.31

# print out session info so we have it on record
sessionInfo()


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PID <- Sys.getpid()
tic("Entire script") # Start master timer

## Set seed for reproducibility
message('Setting seed 98118 for reproducibility')
set.seed(98118)

# Setup -------------------------------------------------------------------------

# create outputdir
outputdir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(outputdir)

# Get formatted tables ------------------------------------------------------------

share_dir <- '<<<< FILEPATH REDACTED >>>>'

tic('Make table')
dt <- format_admin_results(ind_gp = indicator_groups,
                           ind = indicators,
                           rd = run_dates,
                           measure = measures,
                           suffix = suffixes,
                           rk = rks)

#reset key 
setkey(dt, code, year, agg_level, draw)
toc(log = TRUE)


# Deaths averted ------------------------------------------------------------

# read in population from diarrhea files
load(paste0('<<<< FILEPATH REDACTED >>>>/', indicators[[2]], ifelse(rks[[2]], '_raked', '_unraked'), measures[[2]], '_admin_draws', suffixes[[2]], '.RData'))
rm(admin_0, admin_1)

# subset data table
dt <- dt[agg_level == 'ADM2' & (year == 2000 | year == 2017)]

# merge on population counts
setnames(admin_2, 'ADM2_CODE', 'code')
dt <- merge(dt, admin_2[, c('year', 'code', 'pop'), with = FALSE], by = c('year', 'code'))
rm(admin_2)

# get diarrhea deaths
dt[, diarrhea_deaths := had_diarrhea*pop]

# convert to RR of death without ORS
no_ors_rr <- round(1/ors_rr, 2)

# multiple by factor for sensitivity analysis
no_ors_rr <- no_ors_rr*sensitivity_factor

# calculate deaths_averted
dt[, no_ors := 1 - ors, by = c('code', 'draw', 'year')]
setkey(dt, code, draw)
dt[, ors_deaths_averted := get_deaths_averted(last(no_ors), first(no_ors), last(diarrhea_deaths), no_ors_rr), by = c('code', 'draw')]

# only keep year 2017
dt <- dt[year == 2017]

# reshape wide
tic('Reshaping back to wide')
dt[, c('ors', 'had_diarrhea', 'diarrhea_deaths', 'no_ors') := NULL]
dt <- dcast(dt, ... ~ draw, value.var = c('ors_deaths_averted'))
toc(log = TRUE)

# save draws
tic('Saving draw level results')
message(sprintf('TESTING: Percent of NA rows per column is: %f%%', mean(is.na(dt[, V1]))))

out.path <- '<<<< FILEPATH REDACTED >>>>'

dir.create('<<<< FILEPATH REDACTED >>>>')

message('-- finished making deaths averted across draws. now saving as \n', out.path)

saveRDS(dt, out.path)
toc(log = TRUE)

# save summary file
tic('Saving summary of draw level results')

summary <- data.table(ADM2_CODE = dt$code,
                      year = dt$year,
                      pop = dt$pop,
                      mean = rowMeans(dt[, grep('V1', names(dt))[1]:ncol(dt)], na.rm = T),
                      lower = apply(dt[, grep('V1', names(dt))[1]:ncol(dt)], 1, lower),
                      upper = apply(dt[, grep('V1', names(dt))[1]:ncol(dt)], 1, upper))

summary <- merge(summary, sp_hierarchy_list, by = 'ADM2_CODE')

out.path <- '<<<< FILEPATH REDACTED >>>>'

write.csv(summary, out.path)
toc(log = TRUE)

# save summary file with rate averted
tic('Saving summary of draw level results for rate')

summary[, mean := mean/pop*1000]
summary[, upper := upper/pop*1000]
summary[, lower := lower/pop*1000]

out.path <- '<<<< FILEPATH REDACTED >>>>'

write.csv(summary, out.path)
toc(log = TRUE)

# ------------------------------------------------------------------------------------------------------------