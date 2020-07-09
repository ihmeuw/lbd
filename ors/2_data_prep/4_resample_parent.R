####################################################################################################
## Resample polygons parent script for oral rehydration therapy (ORT) indicators
####################################################################################################


## Setup -------------------------------------------------------------------------
rm(list = ls())
library('data.table')

# get input version from most recently collapsed data
files <- file.info(list.files(paste0('<<<< FILEPATH REDACTED >>>>'), pattern = '*.RData', full.names=TRUE))
files <- rownames(files)
files <- gsub('.RData', '', files)
files <- gsub('<<<< FILEPATH REDACTED >>>>/alldata_', '', files)
input_version <- max(files)

# list all indicators
indics <- c('ors', 'ors_or_rhf', 'rhf')

# indicate country iso3 to resample (NULL if not country-specific)
country_code <- NULL

# set arguments
user <- 'kewiens'
run_date <- Sys.Date()

# set cluster arguments
use_geos_nodes  <- TRUE
proj_arg        <- 'proj_geo_nodes_dia'

# load data
load(paste0('<<<< FILEPATH REDACTED >>>>/countdata_', input_version, '.RData'))


## Run child script to resample polygons -------------------------------------------------------------------------

for (i in indics) {
  # indicate indicator
  indicator <- i
  # subset by indicator and polygon
  polydat <- countdata[[i]]
  polydat <- polydat[point == 0]
  # subset by country (if running country-specific)
  if (!is.null(country_code)) {polydat <- polydat[country == country_code]}
  # remove rows where we're missing shapefiles
  polydat <- polydat[shapefile != '']
  # submit qsubs and run child script
  for (shp in unique(polydat$shapefile)) { 
    jname <- paste(indicator, shp, sep = '_')
    mymem <- '15G'
    sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output ', 
                      '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                      '-l fthread=1 -l h_rt=00:02:00:00 -v sing_image=default -N ', jname, ' -l archive=TRUE')
    script <- '<<<< FILEPATH REDACTED >>>>/child.R'
    r_shell <- '<<<< FILEPATH REDACTED >>>>/shell_sing.sh'
    args <- paste(shp, indicator, run_date, input_version)
    system(paste(sys.sub, r_shell, script, args))
  }
}
