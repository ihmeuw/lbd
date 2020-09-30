#####################################################################
# Launcher for raking malaria incidence/mortality to YLDs/YLLs
# uses 02_rake_malaria_cell_pred.R
#####################################################################

# Setup -------------------------------------------------------------
rm(list = ls())

# user arguments:
run_date <- '2019_10_28' #for malaria output
shapefile_version <- '2019_09_10'
cell_pred_year_list <- c(2000:2016) #years in cell pred
goal_year_list <- c(2000:2017) #years to rake to
rake_to_list <- c('mortality','incidence','yll','yld') #measures to rake to
copy_2017 <- TRUE #if true, read in the 2000:2016 cell pred and copy 2017
pop_release <- '2019_08_29'

# Make qsub string and submit -----------------------------------------
mem <- '300G'
rt <- '48:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
shell <- '<<< FILEPATH REDACTED >>>'
code <- '<<< FILEPATH REDACTED >>>/01_calculate_dalys/02_rake_malaria_cell_pred.R'

for (rake_to in rake_to_list){
  name <- paste0('rake_malaria_',rake_to)
  args <- paste(run_date,
                shapefile_version,
                paste(cell_pred_year_list, collapse = '~'),
                paste(goal_year_list, collapse = '~'),
                rake_to,
                copy_2017,
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
