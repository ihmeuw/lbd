library('data.table')

indicator = 'has_lri'
indicator_group = 'lri'
run_date = '<<<< RUN DATE REDACTED >>>>'

#load region list
source('<<<< FILEPATH REDACTED >>>>') #post proc functions
source('<<<< FILEPATH REDACTED >>>>') #misc functions     

temprun = '<<<< FILEPATH REDACTED >>>>'
region_list = fetch_from_rdata(temprun,item_name = 'region_list')

for(rrr in region_list){
  
  qsub = build_qsub(job_name = paste0('param_',rrr),
                    output_folder = '<<<< FILEPATH REDACTED >>>>',
                    error_folder =  '<<<< FILEPATH REDACTED >>>>',
                    make_folders = T,
                    shell_path = '<<<< FILEPATH REDACTED >>>>',
                    script_path = '<<<< FILEPATH REDACTED >>>>',
                    additional_options = ' -v sing_image=default',
                    project = '<<<< PROJECT NAME REDACTED >>>>',
                    slots = 20,
                    arguments_string = paste(indicator, indicator_group, run_date, rrr, 'blarg'))
  system(qsub)
  
}