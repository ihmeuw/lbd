##############################################################################
## Prepare additional diarrhea estimates for figures and viz tool
##
## Includes:
##  AROC
##  Mean
##  Admin 0, Admin 1, Admin 2
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'ort'

# Load functions
library(data.table)
library(raster)

# set indicator
indicators <- 'had_diarrhea'

# set measure arguments
indicator_run_date <- '2019_09_17_14_12_53'
measures <- c('_prevalence', '_incidence', '_deaths')
rake <- '_raked'


## Start saving -------------------------------------------------------------------------

# loop over measures
for (measure in measures) {
  
  # set directories
  in_dir <- '<<<< FILEPATH REDACTED >>>>'
  out_dir <- '<<<< FILEPATH REDACTED >>>>'
  out_dir2 <- '<<<< FILEPATH REDACTED >>>>'
  
  # create directories
  dir.create(out_dir)
  dir.create(out_dir2)
  
  
  ## Save AROC -------------------------------------------------------------------------
  
  # loop over admins
  for (i in 0:2) {
    
    # read in aggregates and save for viz team
    aroc <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '_admin_', i,
                         paste0('_raked', measure), '_aroc_summary.csv'))
    write.csv(aroc[, V1 := NULL], paste0(out_dir, indicator, measure, '_aroc_ad', i, '.csv'))
    
    # save for figures
    for (m in c('mean', 'upper', 'lower')) {
      aroc2 <- aroc[, grep(paste0('ADM', i, '_CODE|year|', m), names(aroc)), with = FALSE]
      setnames(aroc2, m, 'value')
      write.csv(aroc2, paste0(out_dir2, indicator, ifelse(measure == '_deaths', '_mortality', measure), 
                              '_aroc_', m, rake, '_ad', i, '.csv'))
    }
    
  }
  
}