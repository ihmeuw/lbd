#Launch postprocessing
library('data.table')

indicator = 'has_lri'
indicator_group = 'lri'
run_dates = '<<<< RUN DATES REDACTED >>>>'
pos_rd = '<<<< RUN DATES REDACTED >>>>'
perms =  c('_gpcovs','_gponly','_stacking','_covs',"")
run_dates = '<<<< RUN DATES REDACTED >>>>'

slots = 4
geos_nodes = T
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
shelldir = '<<<< FILEPATH REDACTED >>>>'

#load misc_functions
source('<<<< FILEPATH REDACTED >>>>')


#for each run_date; run post processing steps
for(run_date in run_dates){

  #decide steps
  shiny_diagnostics = F
  make_aggregates   = T
  make_rasters      = T
  oos_statistics    = T
  coverage_report   = F
  line_plots        = F
  combine_aggs      = T
  little_csvs       = T
  arocproj_results  = F
  arocproj_graphics = F
  oos_summaries     = T



  pp_qs = build_qsub(job_name = paste0('pos_proc_',run_date),
                     output_folder = '<<<< FILEPATH REDACTED >>>>',
                     error_folder = '<<<< FILEPATH REDACTED >>>>',
                     make_folders = T,
                     shell_path = '<<<< FILEPATH REDACTED >>>>',
                     script_path = '<<<< FILEPATH REDACTED >>>>',
                     additional_options = paste0('-v SET_MKL_THREADS=',floor(slots * .6), ' -v sing_image=default'),
                     project = '<<<< PROJECT NAME REDACTED >>>>',
                     geos_node = geos_nodes,
                     slots = slots,
                     num_tasks = 0,
                     arguments_string = paste(indicator, indicator_group, run_date,
                                              shiny_diagnostics, make_aggregates, make_rasters, oos_statistics,
                                              coverage_report, line_plots,combine_aggs, little_csvs,arocproj_results, arocproj_graphics, oos_summaries, 'blarg'))

  system(pp_qs)

}
