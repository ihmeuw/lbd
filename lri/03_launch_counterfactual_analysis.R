#####################################################################
# Launcher for LRI counterfactual analysis
# uses specific child scripts for each risk group
# currently set up to calculate:
      #(1) CGF PAFs (wasting, stunting, underweight)
      #(2) TAP PAFs
#####################################################################

# Setup -------------------------------------------------------------
rm(list = ls())

# user arguments
risks_to_calc <- c('tap') #currently set up for cgf, tap

lri_run_date <- '2019_10_23_16_13_17'
stunting_run_date <- '2019_09_30_11_45_54'
wasting_run_date <- '2019_09_30_11_45_54'
underweight_run_date <- '2019_09_30_11_45_54'
tap_run_date <- '2019-11-20'

# Make qsub string and submit for each risk -----------------------------------------
mem <- '20G'
rt <- '02:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
shell <- '<<<< FILEPATH REDACTED >>>>'


for (risk in risks_to_calc){
  name <- paste0(risk, '_counterfact')
  
  if (risk == 'cgf'){
    code <-'<<<< FILEPATH REDACTED >>>>'
    args <- paste(lri_run_date,
                  stunting_run_date,
                  wasting_run_date,
                  underweight_run_date)
    
  } else if (risk == 'tap'){
    code <- '<<<< FILEPATH REDACTED >>>>'
    args <- paste(lri_run_date,
                  tap_run_date)
  }
  
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
