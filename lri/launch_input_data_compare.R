#####################################################################
# Launch input data comparison graphs
#####################################################################

# (1) User arguments  ----------------------------------------------------------

old_tag <- 'has_lri_test_2018_05_03' #name in input data folder. if custom_old_data = T, specify any name for the "old" data to propogate to the file name
new_tags <- c('has_lri_current')

custom_old_data <- F
old_data_filepath <- '<<<< FILEPATH REDACTED >>>>' #stage 1 data


# Make qsub string and submit -----------------------------------------
mem <- '10G'
rt <- '00:30:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
name <- 'input_data_plots'
shell <- '<<<< FILEPATH REDACTED >>>>'
code <- '<<<< FILEPATH REDACTED >>>>/lri/input_data_compare.R'

for (new_tag in new_tags){
  args <- paste(old_tag,
                new_tag,
                custom_old_data,
                old_data_filepath)
  
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
