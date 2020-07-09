##############################################################################
## Prepare additional ORT estimates for figures and viz tool
##
## Includes:
##  AROC and correlation analyses
##  Mean
##  Admin 0, Admin 1
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
indicators <- c('ors', 'rhf', 'ors_or_rhf')

# set measure, raked, and run_date 
indicator_run_date <- '2019_11_07_13_50_11'
measures <- ''
rake <- '_unraked'


## Set indicator information -------------------------------------------------------------------------

for (indicator in indicators) {
  message(indicator)
  
  # loop over measures
  for (measure in measures) {
    
    # set directories
    in_dir <- paste0('<<<< FILEPATH REDACTED >>>>')
    out_dir <- paste0(in_dir, 'cleaned_results/')
    out_dir2 <- paste0('<<<< FILEPATH REDACTED >>>>')
    
    # create directories
    dir.create(paste0('<<<< FILEPATH REDACTED >>>>'), showWarnings = FALSE)
    dir.create(out_dir, showWarnings = FALSE)
    
    
    ## Save AROC -------------------------------------------------------------------------

    # loop over admins
    for (i in 0:2) {

      # read in aggregates and save for viz team
      aroc <- fread(paste0('<<<< FILEPATH REDACTED >>>>', indicator, '_admin_', i,
                           '_unraked_prevalence_aroc_summary.csv'))
      write.csv(aroc[, V1 := NULL], paste0(out_dir, indicator, measure, '_aroc_ad', i, '.csv'))

      # save for figures
      for (m in c('mean', 'upper', 'lower')) {
        aroc2 <- aroc[, grep(paste0('ADM', i, '_CODE|year|', m), names(aroc)), with = FALSE]
        setnames(aroc2, m, 'value')
        write.csv(aroc2, paste0(out_dir2, indicator, measure, '_aroc_', m, rake, '_ad', i, '.csv'))
      }
    }


    ## Save correlation -------------------------------------------------------------------------

    if (indicator == 'ors') {

      # loop over admins
      for (i in 0:2) {

        # read in diarrhea vs. ors aggregates and save
        ad <- fread(paste0('<<<< FILEPATH REDACTED >>>>/dia_ors_corr_adm', i, '_summary.csv'))
        write.csv(ad[, V1 := NULL], paste0(out_dir, indicator, measure, '_dia_ors_correlation_ad', i, '.csv'))
        ad <- unique(ad[, grep(paste0('ADM', i, '_CODE|mean'), names(ad)), with = FALSE])
        setnames(ad, 'mean', 'value')
        ad[, year := 2017]
        write.csv(ad, paste0(out_dir2, 'dia_ors_correlation_mean_unraked_ad', i, '.csv'))

        # read in diarrhea vs. rhf aggregates and save
        ad <- fread(paste0('<<<< FILEPATH REDACTED >>>>/dia_rhf_corr_adm', i, '_summary.csv'))
        write.csv(ad[, V1 := NULL], paste0(out_dir, indicator, measure, '_dia_rhf_correlation_ad', i, '.csv'))
        ad <- unique(ad[, grep(paste0('ADM', i, '_CODE|mean'), names(ad)), with = FALSE])
        setnames(ad, 'mean', 'value')
        ad[, year := 2017]
        write.csv(ad, paste0(out_dir2, 'dia_rhf_correlation_mean_unraked_ad', i, '.csv'))

        # read in ors vs. rhf aggregates and save
        ad <- fread(paste0('<<<< FILEPATH REDACTED >>>>/ors_rhf_corr_adm', i, '_summary.csv'))
        write.csv(ad[, V1 := NULL], paste0(out_dir, indicator, measure, '_ors_rhf_correlation_ad', i, '.csv'))
        ad <- unique(ad[, grep(paste0('ADM', i, '_CODE|mean'), names(ad)), with = FALSE])
        setnames(ad, 'mean', 'value')
        ad[, year := 2017]
        write.csv(ad, paste0(out_dir2, 'ors_rhf_correlation_mean_unraked_ad', i, '.csv'))
      }

    }


    ## Save deaths averted -------------------------------------------------------------------------

    if (indicator == 'ors') {
      
      get_filename <- function(x) {
        names <- list('ors_rate_averted_lo' = 'ors_rate_averted_admin_2_summary_0.5', 
                      'ors_rate_averted' = 'ors_rate_averted_admin_2_summary_1', 
                      'ors_rate_averted_hi' = 'ors_rate_averted_admin_2_summary_2',
                      'ors_deaths_averted_lo' = 'ors_deaths_averted_admin_2_summary_0.5', 
                      'ors_deaths_averted' = 'ors_deaths_averted_admin_2_summary_1', 
                      'ors_deaths_averted_hi' = 'ors_deaths_averted_admin_2_summary_2',
                      'deaths_attributable_lo' = 'ors_deaths_attributable_admin_2_summary_0.5', 
                      'deaths_attributable' = 'ors_deaths_attributable_admin_2_summary_1', 
                      'deaths_attributable_hi' = 'ors_deaths_attributable_admin_2_summary_2')
        return(names[[x]])
      }

      for (m in c('ors_rate_averted_lo', 'ors_rate_averted', 'ors_rate_averted_hi',
                  'ors_deaths_averted_lo', 'ors_deaths_averted', 'ors_deaths_averted_hi',
                  'deaths_attributable_lo', 'deaths_attributable', 'deaths_attributable_hi')) {
        
        ad <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', get_filename(m), '.csv'))
        ad <- unique(ad[, c('ADM2_CODE', 'mean'), with = FALSE])
        setnames(ad, 'mean', 'value')
        ad[, year := 2017]
        write.csv(ad, paste0(out_dir2, m, '_mean_unraked_ad2.csv'))
        write.csv(ad, paste0(out_dir, m, '_mean_ad2.csv'))
      }

    }
    
  }
  
}