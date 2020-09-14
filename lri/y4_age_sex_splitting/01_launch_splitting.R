##############################################################################
## Age and/or sex splitting launcher script
## Indicators: has_lri
##
## Launches:
## age_sex_split_cell_pred.R
##		- Performs age-sex splits for cell pred draw-level results
##		- Runs in parallel over regions
##		- Requires a massive amount of memory, so includes a check to
##		  make sure that you really want to do this
## age_sex_split_raster.R
##		- Performs age-sex splits for mean, upper, and lower rasters
##    - Only performs age-sex splits for raked estimates
##		- Runs in parallel over regions
## age_sex_split_admin.R
##		- Performs age-sex splits for aggregated admin 0, 1, and 2 
##		  draw-level results
##		- Runs one script for all regions
##
## Instructions:
##		- Make sure all MBG arguments are correct in the setup block
##    - Make sure all GBD arguments are correct in the scripts themselves
##		- Only run the cell pred script if you absolutely need cell draws,
##		  otherwise the raster and admin scripts should be sufficient
##    - Currently only compatible with the measures: 
##          deaths, daly, yll, yld, deaths or mortality, prevalence
##############################################################################

## Setup -----------------------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- file.path('/homes', user, '_code/lbd/')

# set project
use_geos_nodes  <- TRUE
proj_arg <- ifelse(use_geos_nodes, 'proj_geo_nodes_dia', 'proj_geospatial_dia')

# indicate which levels to run splits for
levels          <- c('admin', 'raster')
levels          <- 'admin'
split_by        <- 'age_sex' #run both for Y5
raked           <- T
modeling_shapefile_version <- '2019_09_10'

#indicate which group we are working on (lri | ort=diarrhea)
indicator_group <- 'lri'

if(indicator_group=='lri') {
  # set MBG arguments for LRI
  indicator       <- 'has_lri'
  run_date        <- '2020_06_11_11_19_26'
  regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_sssa',
             'dia_mcaca', 'dia_s_america_n', 'dia_s_america_s', 'dia_central_asia',
             'dia_se_asia', 'dia_malay', 'dia_south_asia-IND', 'dia_mid_east', 'dia_essa', 'IND','MNG')
  #regions <- 'IND'
  measures_list   <- c('prevalence', 'incidence', 'mortality')

} else {
  # set MBG arguments for diarrhea
  indicator       <- 'had_diarrhea'
  run_date        <- '2019_09_17_14_12_54'
  regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_sssa', 'dia_name', 'dia_s_america',
               'dia_mcaca', 'dia_central_asia', 'MNG', 'IND',
               'dia_se_asia', 'dia_malay', 'dia_south_asia-ind', 'dia_mid_east', 'dia_essa')
  regions <- 'dia_s_america'
  measures_list   <- c('prevalence', 'incidence', 'deaths')
}

## Launch -----------------------------------------------------------------------------------------
# set up base qsub
sys.sub <- paste0('qsub -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                  '-v sing_image= <<< FILEPATH REDACTED >>> ',
                  '-e <<< FILEPATH REDACTED >>>',
                  '/errors -o <<< FILEPATH REDACTED >>> ',
                  '-l archive=TRUE -N')
r_shell <- '<<< FILEPATH REDACTED >>>'


## Submit qsub different splits -----------------------------------------------------
for (measure in measures_list){
 
  ## Submit qsub for raster-level splits --------------------------------------------------------
  
  if ('raster' %in% levels) {
    
    # loop over regions
    for (region in regions) {
      
      # set qsub arguments
      jname <- paste0(indicator,'_', region, '_', measure, '_raster_split_', split_by)
      mem <- '-l m_mem_free=25G'
      time <- '-l h_rt=01:00:00:00'
      thread <- '-v SET_OMP_THREADS=2 -v SET_MKL_THREADS=2 -l fthread=2'
      
      # set script arguments
      script <- '<<< FILEPATH REDACTED >>>'
      args <- paste(repo, indicator_group, indicator, measure, 
                    run_date, region, modeling_shapefile_version, split_by)
      
      # submit job
      system(paste(sys.sub, jname, mem, time, thread, r_shell, script, args))
      
    }
    
  }
  
  
  ## Submit qsub for admin-level splits -------------------------------------------------------
  
  if ('admin' %in% levels) {
    
    # set qsub arguments
    jname <- paste0(indicator, '_', measure, '_admin_split_', split_by)
    mem <- '-l m_mem_free=200G'
    time <- '-l h_rt=01:00:00:00'
    thread <- '-v SET_OMP_THREADS=4 -v SET_MKL_THREADS=4 -l fthread=4'
    
    # set script arguments
    script <- '<<< FILEPATH REDACTED >>>'
    args <- paste(repo, indicator_group, indicator, raked, measure, 
                  run_date, modeling_shapefile_version, split_by)
    
    # submit job
    system(paste(sys.sub, jname, mem, time, thread, r_shell, script, args))
    
  }
}