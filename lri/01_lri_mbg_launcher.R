#####################################################################
# Launcher for LRI MBG models & raking
# uses launch_lri.R
#####################################################################

# Setup -------------------------------------------------------------

## clear environment
rm(list=ls())

## set region lists
all_stg1 <- c('cssa','wssa','essa','name','sssa')

stg_2_best <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_sssa',
               'dia_mcaca', 'dia_s_america_n', 'dia_s_america_s', 'dia_central_asia',
               'dia_se_asia', 'dia_malay', 'dia_south_asia-IND', 'dia_mid_east', 'dia_essa', 'IND','MNG')

## set arguments
Regions         <- stg_2_best
use_run_date    <- NA #NA if new run date

region_specific_cov <- TRUE
config_par      <- '82'
cov_par         <- '27'
gbm_par         <- 'NA' #'NA' if no gbm parameters

priority        <- '0' #'0' for default priority, or a number 0 > n > -1003 for low priority

make_holdouts   <- TRUE
holdout_type <- 'nids' #'nids' or 'admin'

go_to_inla      <- FALSE  #skips to INLA, reading in stackers from saved model image history
skip_to_rake    <- FALSE #skips to raking step, reading in cell preds from run date
go_to_inla_from_rundate <- use_run_date #use_run_date for same run date, or specify another
etiology        <- 'none' #options: lri_flu, lri_rsv, lri_pneumo, lri_hib or 'none' to skip etiology splitting
rake_measures   <- c('incidence', 'prevalence', 'mortality') #what measures to rake to/aggregate over (if doing etiologies, can only be mort or inc)

## several arguments generally should not change:
indicator_group <- 'lri'
indicator       <- 'has_lri'
user            <- Sys.info()['user']
core_repo       <- '<<<< FILEPATH REDACTED >>>>'
indicator_repo  <- '<<<< FILEPATH REDACTED >>>>'

# Make qsub string and submit -----------------------------------------
mem <- '1G'
rt <- '120:00:00'
queue <- 'geospatial.q'
proj <- 'proj_geo_nodes'
name <- 'lri_mbg_launch'
shell <- '<<<< FILEPATH REDACTED >>>>'
code <- '<<<< FILEPATH REDACTED >>>>'

args <- paste(user,
              use_run_date,
              indicator_group,
              indicator,
              core_repo,
              indicator_repo,
              config_par,
              cov_par,
              gbm_par,
              region_specific_cov,
              paste(Regions, collapse = '~'),
              go_to_inla,
              paste(rake_measures, collapse = '~'),
              skip_to_rake,
              etiology,
              priority,
              go_to_inla_from_rundate,
              make_holdouts,
              holdout_type)

qsub <- paste0('qsub -l m_mem_free=', mem,
               ' -l fthread=1 -l h_rt=', rt, 
               ' -v sing_image=default -q ', queue,
               ' -p ', priority,
               ' -P ', proj,  
               ' -N ', name, 
               ' ', shell, 
               ' ', code, 
               ' ', args)

system(qsub)

