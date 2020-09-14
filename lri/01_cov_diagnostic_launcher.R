#####################################################################
# Launcher for LRI covariate time trend plots
# uses cov_diagnostic.R
#####################################################################

# Setup -------------------------------------------------------------

## clear environment
rm(list=ls())

## set arguments
exp_num <- 4 #experiment number from LRI model tracker
test <- 0
holdout <- 0
age <- 0

## libraries
library(googlesheets)
library(data.table)
library(dplyr)

## read in list of pixels of interest from LRI model tracker
ttt <- readRDS('<<<< FILEPATH REDACTED >>>>')
gs_auth(token = ttt)
lri_tracker <- gs_title('LRI model tracker')
pixel_info <- gs_read(lri_tracker, ws = 'plot covariate time trends') %>%
  filter(experiment_number == exp_num)

## pull run date and stackers to plot
run_date <- unique(pixel_info$model_run_date)
stackers <- unique(pixel_info$stackers) %>%
  paste(sep = '',collapse = ',') %>%
  strsplit(',') %>%
  unlist() %>%
  unique()


## for each row in pixel info, make a qsub string to launch cov_diagnostic.R
for (n in 1:nrow(pixel_info)){
  info <- pixel_info[n,]
  region <- info$region
  country <- info$country
  lat <- info$lat
  long <- info$long
  pixel_id <- info$pixel_id
  
  mem <- '15G'
  rt <- '05:00:00'
  queue <- 'geospatial.q'
  proj <- 'proj_geo_nodes'
  name <- paste0(country, '_cov_diag')
  shell <- '<<<< FILEPATH REDACTED >>>>'
  code <- '<<<< FILEPATH REDACTED >>>>'
  
  args <- paste(run_date,
                region,
                country,
                lat,
                long,
                paste(stackers, collapse = '~'),
                test,
                holdout,
                age,
                pixel_id)
  
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
