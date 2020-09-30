########################################################################################
# Aggregate LRI and diarrhea DALYs/YLDs to admin 0, 1, 2 
########################################################################################

# (1) Setup ------------------------------------------------------------------------
rm(list = ls())

#Arguments from qsub call
indicator_group <- 'lri' #lri or ort
measures <- c('yll','yld','daly') #yll, yld, daly
run_date <- '2019_10_28' #tridaly
shapefile_version <- '2019_10_28'

lri_reg_list <- c('dia_wssa','dia_sssa','dia_name-ESH','dia_essa','dia_cssa','dia_afr_horn')
dia_reg_list <- c('dia_wssa','dia_sssa','dia_name','dia_essa','dia_cssa','dia_afr_horn')

use_inla_country_fes <- F

#set arguments automatically by indicator
if (indicator_group == 'ort'){
  indicator <- 'had_diarrhea'
  regions <- dia_reg_list
}
if (indicator_group == 'lri'){
  indicator <- 'has_lri'
  regions<- lri_reg_list
}

#functions
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/misc_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core_custom/misc_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/qsub_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/post_estimation_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core_custom/post_estimation_functions.R')

library(dplyr)
library(data.table)

# (2) Summarize aggregations  ---------------------------------------
for (measure in measures){
  message('Combining aggregated results')
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = regions,
                      holdouts = 0,
                      raked    = T,
                      delete_region_files = F,
                      measure = measure,
                      metrics = c('rates','counts'),
                      dir_to_search = '<<< FILEPATH REDACTED >>>')
  
  summarize_admins_tridalys(ad_levels = c(0,1,2), 
                            raked = T, 
                            measure = measure, 
                            metrics = c('rates','counts'),
                            input_dir =  '<<< FILEPATH REDACTED >>>',
                            output_dir = '<<< FILEPATH REDACTED >>>')
}



