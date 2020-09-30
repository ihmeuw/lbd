#########################################################################
# Convert incidence to mortality, YLL, YLD, and DAlYs 
# for malaria counterfactual scenarios
#########################################################################

# (1) Setup --------------------------------------------------------------
rm(list = ls())
library(dplyr)

#user inputs
run_date <- '2019_10_28' #tridalys
scenarios <- c('act','all','irs','itn')
admin_levels <- c(0:2)

#set directories
in_dir <- '<<< FILEPATH REDACTED >>>'

#read in adm0 incidence rates and take means
mortality_inputs <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(Year %in% c(2000:2017)) %>%
  mutate(mean_count = rowMeans(dplyr::select(mortality_inputs, starts_with('draws_')), na.rm = TRUE)) %>%
  mutate(mean_rate = mean_count / Pop) %>%
  dplyr::select(IHME_location_id, Year, metric, mean_rate, mean_count, Pop, MAP_Country_Name)

scenario <- 'act'
admin <- 1

for (scenario in scenarios){
  for (admin in admin_levels){
  
    # (2) Convert incidence to mortality ------------------------------------
    adx_inc <- fread(paste0(in_dir, 'had_malaria_', scenario, '_admin_', admin, '_raked_incidence_summary.csv'))
    ad0_mort <- filter(mortality_inputs, Year %in% c(2000:2017))
    
  }
}



