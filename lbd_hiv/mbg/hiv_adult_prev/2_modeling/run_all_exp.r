####################################################################################################
## Launch jobs for all IS and OOS jobs
####################################################################################################

## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))

} else {
  models <- list(c('Full model for EPP--lso phia adm codes fixed',
                   'epp_config',
                   'cov_list_std_and_hiv'))

  models <- data.frame(do.call('rbind', models))
  names(models) <- c('run_label', 'config_file', 'covs_file')

}

## Define function for submitting launch scripts
qsub_launch <- function(config_file, covs_file, run_date = NULL) {

  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  outdir <- paste0("<<<< FILEPATH REDACTED >>>>")
  dir.create(outdir, showWarnings = F)
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)

  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), '-l m_mem_free=150G -l fthread=5 -N', paste0('launch_', run_date),
                '-P proj_geo_nodes -q geospatial.q -l h_rt=08:00:00:00 -l archive=TRUE',
                '-v "<<<< FILEPATH REDACTED >>>>" -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0("<<<< FILEPATH REDACTED >>>>",'/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0("<<<< FILEPATH REDACTED >>>>",'/lbd_hiv/mbg/hiv_adult_prev/2_modeling/launch.r'),
                config_file, covs_file, run_date)

  system(qsub)
  Sys.sleep(10)
  return(run_date)
}

## Submit jobs and save job submission info
if (!is.null(resub_date)) {
  apply(models, 1, function(x) qsub_launch(x[2], x[3]))

} else {
  models$run_date <- apply(models, 1, function(x) qsub_launch(x[2], x[3]))
  write.table(models, file = paste0("<<<< FILEPATH REDACTED >>>>"),
              sep = ",", row.names = F, append = T)

}
