####################################################################
# Launcher for LRI AROC
# uses make_aroc.R
#####################################################################

# Setup -------------------------------------------------------------
rm(list = ls())

# user arguments
run_date                 <- '2020_06_11_11_19_26'
modeling_shapefile_version <- '2019_09_10'
end_year <- 2019
mortality_2010_end_year <- TRUE

#other arguments
indicator_group          <- 'lri'
indicator                <- 'has_lri'


# Make qsub string and submit -----------------------------------------
mem <- '200G'
rt <- '10:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
name <- 'lri_aroc'
shell <- '<<<< FILEPATH REDACTED >>>>'
code <- '<<<< FILEPATH REDACTED >>>>'

args <- paste(indicator_group,
              indicator,
              run_date,
              modeling_shapefile_version,
              mortality_2010_end_year,
              end_year)

qsub <- paste0('qsub -l m_mem_free=', mem,
               ' -l fthread=1 -l h_rt=', rt, 
               ' -v sing_image=default -q ', queue, 
               ' -P ', proj, 
               ' -N ', name, 
               ' ', shell, 
               ' ', code, 
               ' ', args)

system(qsub)
