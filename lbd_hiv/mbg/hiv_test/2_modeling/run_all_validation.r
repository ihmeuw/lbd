####################################################################################################
## Launch jobs for main model and all validation/sensitivity analyses
####################################################################################################

## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0('<<<< FILEPATH REDACTED >>>>/models_submitted_', resub_date, '.csv'))

} else {
  models <- list(# standard model
                 c('Standard model',                       'config_main',               'cov_list_std_and_hiv',  'new_geos'),

                 # alternate covariate strategies
                 c('Covariate Test: GP only',              'config_gp_only',            'cov_list_std_and_hiv',  'new_geos'),
                 c('Covariate Test: Raw only, all',        'config_raw_only',           'cov_list_std_and_hiv',  'old_geos'),
                 c('Covariate Test: Stack only, all',      'config_stack_only',         'cov_list_std_and_hiv',  'old_geos'),
                 c('Covariate Test: GP + Raw, all',        'config_raw_and_gp',         'cov_list_std_and_hiv',  'new_geos'),

                 # alternate ANC data strategies
                 c('ANC Test: Survey only',                'config_no_anc',             'cov_list_std_and_hiv',  'old_geos'),
                 c('ANC Test: Survey + ANC',               'config_anc_no_correction',  'cov_list_std_and_hiv',  'old_geos'),
                 c('ANC Test: Survey + ANC (iid)',         'config_anc_iid_correction', 'cov_list_std_and_hiv',  'old_geos'),

                 # alternate polygon data strategies
                 c('Polygon Test: Point data only',        'config_no_poly',            'cov_list_std_and_hiv',  'new_geos'))

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
  write.table(models, file = paste0('<<<< FILEPATH REDACTED >>>>/models_submitted_', format(Sys.Date(), '%Y_%m_%d'), '.csv'),
              sep = ",", row.names = F, append = T)

}
