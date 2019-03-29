####################################################################################################
## Launch jobs to run the prior sensitivity analysis
####################################################################################################

## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0('<<<< FILEPATH REDACTED >>>>/prior_sensitivity_models_submitted_', resub_date, '.csv'))

} else {
  models <- list(
    c('Standard PC model (Model 0)',                   'config_default_pc',              'cov_list_std_and_hiv',  'new_geos'),
    c('Standard PC model, informative rho (Model 5)',  'config_default_pc_strict_rho',   'cov_list_std_and_hiv',  'new_geos'),
    c('Loose PC prior (Model 1)',                      'config_loose_pc',                'cov_list_std_and_hiv',  'new_geos'),
    c('Loose PC prior, informative rho (Model 2)',     'config_loose_pc_strict_rho',     'cov_list_std_and_hiv',  'new_geos'),
    c('INLA default prior (Model 3)',                  'config_inla_default',            'cov_list_std_and_hiv',  'new_geos'),
    c('INLA default prior, informative rho (Model 4)', 'config_inla_default_strict_rho', 'cov_list_std_and_hiv',  'new_geos'))

  models <- data.frame(do.call('rbind', models))
  names(models) <- c('run_label', 'config_file', 'covs_file', 'cluster')

}

## Define function for submitting launch scripts
qsub_launch <- function(config_file, covs_file, cluster, run_date = NULL) {

  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  outdir <- paste0('<<<< FILEPATH REDACTED >>>>')
  dir.create(outdir, showWarnings = F)
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)

  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), '-pe multi_slot 7 -N', paste0('launch_', run_date),
                '-P proj_geo_nodes_hiv -l geos_node=TRUE -q geospatial.q@@intel_e5_2660v4_1024gb',
                '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0('<<<< FILEPATH REDACTED >>>>', '/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0('<<<< FILEPATH REDACTED >>>>', '/lbd_hiv/mbg/hiv_test/2_modeling/launch.r'),
                config_file, covs_file, run_date, cluster)
  system(qsub)
  Sys.sleep(10)
  return(run_date)
}

## Submit jobs and save job submission info
if (!is.null(resub_date)) {
  apply(models, 1, function(x) qsub_launch(x[2], x[3], x[4], x[5]))

} else {
  models$run_date <- apply(models, 1, function(x) qsub_launch(x[2], x[3], x[4], x[5]))
  write.table(models, file = paste0('<<<< FILEPATH REDACTED >>>>/prior_sensitivity_models_submitted_', format(Sys.Date(), '%Y_%m_%d'), '.csv'),
              sep = ",", row.names = F, append = T)

}
