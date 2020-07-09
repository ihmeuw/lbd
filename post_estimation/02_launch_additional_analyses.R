##############################################################################
## Launch additional postestimation and analysis at the draw level
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'

# set project
use_geos_nodes  <- TRUE
proj_arg <- ifelse(use_geos_nodes, 'proj_geo_nodes_dia', 'proj_geospatial')

# set run- and indicator-specific arguments
indicator_group         <- 'ort'
indicator               <- 'ors'
dia_run_date            <- '2019_09_17_14_12_53'
ort_run_date            <- '2019_11_07_13_50_11'
shapefile               <- '2019_09_10'

# set measures, run dates, and regions based on indicator
run_date <- ort_run_date
meas <- 'prevalence'
regions <- c('dia_afr_horn-eth-yem', 'ETH', 'YEM', 'dia_name', 'dia_sssa', 
             'dia_mcaca', 'dia_central_asia', 'MNG', 
             'dia_se_asia', 'dia_malay', 'dia_mid_east',
             'ZWE', 'KEN', 'NGA', 'COD', 'IND', 'PAK',
             'dia_essa-zwe-ken', 'dia_cssa-cod', 'dia_south_asia-ind-pak',
             'dia_s_america_n', 'dia_s_america_s', 'dia_wssa-nga')

# set up base qsub
sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output', 
                  ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                  '-l fthread=1 -l h_rt=01:00:00:00 -v sing_image=default -l archive=TRUE -N')
r_shell <- paste0(repo, '/mbg_central/share_scripts/shell_sing.sh')

# function to get conversion factors
get_conversion_factor <- function(r) {
  region_factors <- 6
  if(r == 'MNG' | r == 'dia_central_asia' | r == 'PAK') region_factors <- 4
  if(r == 'KEN' | r == 'ZWE') region_factors <- 2
  if(r == 'dia_malay' | r == 'dia_name') region_factors <- 8
  if(r == 'dia_wssa' | r == 'dia_wssa-nga') region_factors <- 10
  if(r == 'dia_s_america') region_factors <- 12
  return(region_factors)
}


## Submit qsub for deaths averted -----------------------------------------------------

# set analysis type
analysis_types <- c('ors_deaths_averted_draws', 'ors_deaths_attributable_draws')

for (analysis_type in analysis_types) {
  
  script <- paste0(repo, '/post_estimation/', analysis_type, '.R')
  
  # loop over analyses
  for (i in c(1, 0.5, 2)) {
    
    # set memory
    mymem <- 80
    mem <- paste0('-l m_mem_free=', mymem, 'G')
    
    # set job names
    jname <- paste0(analysis_type)
    
    # set up args
    args <- paste0(ort_run_date, ' ', dia_run_date, ' ', i, ' ', shapefile)
    
    # run launch script
    system(paste(sys.sub, jname, mem, r_shell, script, args))
    
  }
  
}


## Submit qsub for AROC -------------------------------------------------------------------------

# set analysis type and script
analysis_type <- 'aroc_draws'
script <- paste0(repo, '/post_estimation/', analysis_type, '.R')

# set memory
mymem <- '80'
mem <- paste0('-l m_mem_free=', mymem, 'G')

# set job names
jname <- paste0(analysis_type, '_', indicator)

# set up args
args <- paste0(indicator_group, ' ', indicator, ' ', run_date, ' ', shapefile)

# run launch script
system(paste(sys.sub, jname, mem, r_shell, script, args))


## Submit qsub for correlation -------------------------------------------------------------------------

# set analysis type
analysis_type <- 'ort_correlation_draws'
script <- paste0(repo, '/post_estimation/', analysis_type, '.R')

# loop over regions
for (reg in regions) {
  
  # set job names
  jname <- paste0(analysis_type, '_', reg)
  
  # set memory
  mymem <- get_conversion_factor(reg)*50
  mem <- paste0('-l m_mem_free=', mymem, 'G')
  
  # set up args
  args <- paste0('ort ort ors rhf ', run_date, ' ', run_date, ' ', shapefile, ' ', reg)
  
  # run launch script
  system(paste(sys.sub, jname, mem, r_shell, script, args))
  
}
