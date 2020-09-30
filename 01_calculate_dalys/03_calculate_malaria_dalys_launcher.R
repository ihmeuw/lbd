#####################################################################
# Launcher for summing malaria YlDs and YLLS, raking to DALYs
# uses 03_calculate_malaria_dalys.R
#####################################################################

# Setup -------------------------------------------------------------
rm(list = ls())

# user arguments:
run_date <- '2019_10_28' #for malaria output
shapefile_version <- '2019_09_10'
year_list <- c(2000:2017)
pop_release <- '2019_08_29'

# Make qsub string and submit -----------------------------------------
mem <- '300G'
rt <- '48:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
shell <- '<<< FILEPATH REDACTED >>>'
code <- '<<< FILEPATH REDACTED >>>/01_calculate_dalys/03_calculate_malaria_dalys.R'

name <- paste0('rake_malaria_dalys')
args <- paste(run_date,
              shapefile_version,
              paste(year_list, collapse = '~'),
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

