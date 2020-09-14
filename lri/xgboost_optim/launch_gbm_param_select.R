rm(list = ls())

indicator = 'has_lri'
indicator_group = 'lri'
run_date = '2018_05_10_16_22_42'

region_list <- 'sssa'
region <- region_list[1]
  
for(region in region_list){
  
  qsub = build_qsub(job_name = paste0('param_',rrr),
                    output_folder = '<<< FILEPATH REDACTED >>>',
                    error_folder =  '<<< FILEPATH REDACTED >>>',
                    make_folders = T,
                    shell_path = '<<< FILEPATH REDACTED >>>',
                    script_path = '<<< FILEPATH REDACTED >>>',
                    additional_options = ' -v sing_image=default',
                    project = 'proj_geospatial_ppp',
                    slots = 20,
                    arguments_string = paste(indicator, indicator_group, run_date, rrr, 'blarg'))
  system(qsub)
  
}