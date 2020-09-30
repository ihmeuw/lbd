#####################################################################
# Launcher for raking LRI/diarrhea prevalence to GBD DALYs
# uses 03_rake_lri_diarrhea_frx_launcher.R
#####################################################################

# Setup --------------------------------------------------------------
indicator_group            <- 'lri' #lri or ort
dia_run_date               <- '2019_09_17_14_12_53'
lri_run_date               <- '2019_09_16_16_39_00'

pop_release                <- '2019_08_29' 

lri_reg_list               <- c('dia_wssa','dia_cssa','dia_sssa','dia_name-ESH','dia_essa', 'dia_afr_horn')
dia_reg_list               <- c('dia_wssa','dia_sssa','dia_name','dia_essa','dia_cssa','dia_afr_horn')
user                       <- '<<< USERNAME REDACTED >>>'
core_repo                  <- '<<< FILEPATH REDACTED >>>'
measure                    <- 'yld'
holdout                    <- 0
modeling_shapefile_version <- '2019_09_10'
raking_shapefile_version   <- '2019_09_10'
year_list                  <- c(2000:2017)
dalys_run_date             <- '2019_10_28'
age                        <- 0

dia_reg_list <- 'dia_name'

#set arguments automatically by indicator
if (indicator_group == 'ort'){
  indicator <- 'had_diarrhea'
  region_list <- dia_reg_list
  run_date <- dia_run_date
}
if (indicator_group == 'lri'){
  indicator <- 'has_lri'
  region_list <- lri_reg_list
  run_date <- lri_run_date
}
# Make qsub string and submit -----------------------------------------
for (reg in region_list){
  mem <- '200G'
  rt <- '10:00:00'
  queue <- 'geospatial.q'
  proj <- 'proj_geo_nodes'
  name <- paste0('rake_',reg, '_', indicator, '_', measure)
  shell <- '<<< FILEPATH REDACTED >>>'
  code <- '<<< FILEPATH REDACTED >>>/01_calculate_dalys/03_rake_lri_diarrhea_frx.R'
  
  args <- paste(user,
                core_repo,
                indicator_group,
                indicator,
                reg,
                run_date,
                measure,
                holdout,
                modeling_shapefile_version,
                raking_shapefile_version,
                paste(year_list, collapse = '~'),
                dalys_run_date,
                pop_release)
  
  qsub <- paste0('qsub -l m_mem_free=', mem,
                 ' -l fthread=1 -l h_rt=', rt, 
                 ' -v sing_image=default -q ', queue, 
                 ' -P ', proj, 
                 ' -N ', name, 
                 ' ', shell, 
                 ' ', code, 
                 ' ', args)
  
  system(qsub)
}
