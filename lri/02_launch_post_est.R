#####################################################################
# Launcher for LRI post-estimation
# uses lri_post_estimation.R
# Katie Welgan
#####################################################################

# Setup -------------------------------------------------------------
rm(list = ls())

# user arguments
config_par <- '87'
cov_par <- '27'
run_date <- '2020_06_24_09_58_34'
year_tag <- '2000_2019' #for raster file names, should correspond to modeling year range
two_xgboost <- FALSE #adds a correction for line plots to make sure xgboost and xgboost2 read as separate stackers

user <- Sys.info()['user']
core_repo <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'lri'
indicator <- 'has_lri'
measures <- c('prevalence','incidence','mortality')
force_plot <- FALSE

# Make qsub string and submit -----------------------------------------
mem <- '100G'
rt <- '10:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
name <- 'lri_post_est'
shell <- '<<<< FILEPATH REDACTED >>>>'
code <- '<<<< FILEPATH REDACTED >>>>'

args <- paste(user,
              core_repo,
              indicator_group,
              indicator,
              config_par,
              cov_par,
              run_date,
              paste(measures, collapse = '~'),
              force_plot,
              year_tag,
              two_xgboost)

qsub <- paste0('qsub -l m_mem_free=', mem,
               ' -l fthread=1 -l h_rt=', rt, 
               ' -v sing_image=default -q ', queue, 
               ' -P ', proj, 
               ' -N ', name, 
               ' ', shell, 
               ' ', code, 
               ' ', args)

system(qsub)

