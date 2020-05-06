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

# set run-specific arguments
indicator_group         <- 'ort'
indicator               <- 'had_diarrhea'
run_date                <- '2019_09_17_14_12_53'
shapefile               <- '2019_09_10'

# set measures
meas <- c('prevalence', 'incidence', 'deaths')

# set regions
regions <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_sssa', 'dia_name', 'dia_s_america',
             'dia_mcaca', 'dia_central_asia', 'MNG', 'IND',
             'dia_se_asia', 'dia_malay', 'dia_south_asia-ind', 'dia_mid_east', 'dia_essa')


# set up base qsub
sys.sub <- paste0('qsub -e <<<< FILEPATH REDACTED >>>>', user,'/errors -o <<<< FILEPATH REDACTED >>>>', user, '/output', 
                  ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                  '-l fthread=1 -l h_rt=01:00:00:00 -v sing_image=default -l archive=TRUE -N')
r_shell <- paste0(repo, '/mbg_central/share_scripts/shell_sing.sh')

# function to get conversio factors
get_conversion_factor <- function(r) {
  region_factors <- 6
  if(r == 'MNG' | r == 'dia_central_asia') region_factors <- 4
  if(r == 'dia_malay' | r == 'dia_name') region_factors <- 8
  if(r == 'dia_wssa') region_factors <- 10
  if(r == 'dia_s_america') region_factors <- 12
  return(region_factors)
}


## Submit qsub for deaths averted -----------------------------------------------------

# deaths averted calculations for ORS (as an example)
analysis_type <- 'ors_deaths_averted_draws'
script <- paste0(repo, '/post_estimation/', analysis_type, '.R')

# set memory
mymem <- 100
mem <- paste0('-l m_mem_free=', mymem, 'G')

# set job names
jname <- paste0(analysis_type)

# run launch script
system(paste(sys.sub, jname, mem, r_shell, script))


## Submit qsub for AROC -------------------------------------------------------------------------

# set analysis type and script
analysis_type <- 'aroc_draws'
script <- paste0(repo, '/post_estimation/', analysis_type, '.R')

# set memory
mymem <- '100'
mem <- paste0('-l m_mem_free=', mymem, 'G')

# set job names
jname <- paste0(analysis_type, '_', indicator)

# set up args
args <- paste0(indicator_group, ' ', indicator, ' ', run_date, ' ', shapefile)

# run launch script
system(paste(sys.sub, jname, mem, r_shell, script, args))


## Submit qsub for inequality -----------------------------------------------------

# set analysis type
analysis_type <- 'inequality_draws'
script <- paste0(repo, '/post_estimation/', analysis_type, '.R')

# only run for deaths for diarrhea
meas <- 'deaths'

# loop over regions
for (reg in regions) {

  # loop over inequality types
  for (i in c('rel', 'abs')) {

    # loop over measures
    for (m in meas) {

      # set memory
      mymem <- get_conversion_factor(reg)*50
      mem <- paste0('-l m_mem_free=', mymem, 'G')

      # set job names
      jname <- paste0(i, '_', analysis_type, '_', reg)

      # set up args
      args <- paste0(indicator_group, ' ', indicator, ' ', run_date, ' ', shapefile, ' ', reg, ' ', m, ' ', i)

      # run launch script
      system(paste(sys.sub, jname, mem, r_shell, script, args))

    }

  }

}
