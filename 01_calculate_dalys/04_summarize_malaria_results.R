#######################################################################################
# (1) Summarize malaria raked YLDs, YLLs, DALYs, incidence, mortality to obtain csvs
# (2) Sum malaria raked YLDs and YLLs admin draws to obtain DALYs admin draws
#######################################################################################

# (1) Setup --------------------------------------------------------------------------
rm(list = ls())

#user arguments
run_date <- '2019_10_28' #malaria
measures <- c('incidence', 'mortality', 'yll','yld','daly')

indicator_group <- 'malaria'
indicator <- 'had_malaria'

#load in functions and packages
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/setup.R')
commondir      <- '<<< FILEPATH REDACTED >>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
mbg_setup(package_list = package_list, repos = '<<< FILEPATH REDACTED >>>')

# Load custom lbd core functions
lapply(paste0('<<< FILEPATH REDACTED >>>/lbd_core_custom/',
              list.files('<<< FILEPATH REDACTED >>>/lbd_core_custom/', pattern = 'functions')), source)


use_inla_country_fes <- FALSE

for (measure in measures){
  combine_aggregation(rd       = run_date,
                      indic    = 'had_malaria',
                      ig       = 'malaria',
                      ages     = 0,
                      regions  = 'africa-ESH',
                      holdouts = 0,
                      raked    = T,
                      measure  = measure,
                      delete_region_files = F,
                      metrics = c('rates','counts'))
  
  summarize_admins(ad_levels = c(0,1,2), raked = T, measure = measure, metrics = c('rates','counts'))
}
