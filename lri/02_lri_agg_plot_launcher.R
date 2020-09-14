#####################################################################
# Launcher for LRI aggregation and diagnostic plots
# uses launch_lri_agg_plot.R
#####################################################################

# Setup -------------------------------------------------------------

## clear environment
rm(list=ls())

all_stg1 <- c('cssa','wssa','essa','name','sssa')

all_stg2 <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name', 'dia_sssa',
             'dia_mcaca', 'dia_s_america', 'dia_central_asia', 'dia_chn_mng',
             'dia_se_asia', 'dia_malay', 'dia_south_asia', 'dia_mid_east', 'dia_essa')

regions <- all_stg1
run_date <- '2019_07_23_13_11_06'
holdout <- 0
age <- 0

# Make qsub string and submit -----------------------------------------

mem <- '20G'
rt <- '03:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
name <- 'lri_agg_launch'
shell <- '<<<< FILEPATH REDACTED >>>>'
code <- '<<<< FILEPATH REDACTED >>>>'

args <- paste(paste(regions, collapse = '~'),
              holdout,
              run_date,
              age)

qsub <- paste0('qsub -l m_mem_free=', mem,
               ' -l fthread=1 -l h_rt=', rt, 
               ' -v sing_image=default -q ', queue, 
               ' -P ', proj, 
               ' -N ', name, 
               ' ', shell, 
               ' ', code, 
               ' ', args)

system(qsub)
