
####################################################################################################
## Launch jobs for all male circumcision model
####################################################################################################
rm(list = ls())
## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0('<<<< FILEPATH REDACTED >>>>'))
  
} else {
  models <- list(
    c('Standard PC model (Model 0)',                   'default_pc'),
    c('Standard PC model, informative rho (Model 5)',  'default_pc_strict_rho'),
    c('Loose PC prior (Model 1)',                      'loose_pc'),
    c('Loose PC prior, informative rho (Model 2)',     'loose_pc_strict_rho'),
    c('INLA default prior (Model 3)',                  'inla_default'),
    c('INLA default prior, informative rho (Model 4)', 'inla_default_strict_rho'),
    c('No nugget or random effects',                   'no_nugget'),
    c('No random effects',                             'no_re'),
    c('No fixed effect on time',                       'no_time_fe'))
  
  models <- data.frame(do.call('rbind', models))
  models <- cbind("male_circumcision", models,"a1549m", "none", "FALSE")
  names(models) <- c('indicator', 'run_label', 'prior', 'population', 'data_tag', "s2_mesh")
  
}

## Define function for submitting launch scripts
qsub_launch <- function(indicator, prior, population, data_tag, run_date = NULL) {
  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  outdir <- paste0('<<<< FILEPATH REDACTED >>>>')
  dir.create(outdir, showWarnings = F)
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)
  
  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), '-l m_mem_free=15G -l fthread=5 -N', paste0('launch_', run_date),
                '-P proj_geo_nodes -q geospatial.q -l h_rt=08:00:00:00 -l archive=TRUE',
                '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0('/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0('/lbd_hiv/mbg/covariates/2_modeling/launch.r'),
                indicator, prior, population, data_tag, run_date)
  
  system(qsub)
  Sys.sleep(10)
  return(run_date)
}

## Submit jobs and save job submission info
if (!is.null(resub_date)) {
  apply(models, 1, function(x) qsub_launch(x[1], x[3], x[4], x[5], x[7]))
} else {
  models$run_date <- apply(models, 1, function(x) qsub_launch(x[1], x[3], x[4], x[5]))
  write.table(models, file = paste0('<<<< FILEPATH REDACTED >>>>'),
              sep = ",", row.names = F, append = T)
}