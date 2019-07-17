## clear environment
rm(list=ls()) 

## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0('<<<< FILEPATH REDACTED >>>>/models_submitted_', resub_date, '.csv'))
  
} else {
  models <- list(# standard model
    c('Standard model',                       'config_ebf_main',               'covs_list_ebf_updated',  'old_geos'),

    # alternate covariate strategies
    c('Covariate Test: GP only',              'config_ebf_gp_only',            'covs_list_ebf_updated',  'old_geos'),
    c('Covariate Test: Stack only',           'config_ebf_stack_only',         'covs_list_ebf_updated',  'prod'),
    c('Covariate Test: Raw cov only',         'config_ebf_raw_only',           'covs_list_ebf_updated',  'prod'),
    c('Covariate Test: Raw cov + GP',         'config_ebf_raw_and_gp',         'covs_list_ebf_updated',  'old_geos'),

    # alternate 24 hour recall period
    c('Recall Period Test: 24 h only',         'config_ebf_main_24h_only',     'covs_list_ebf_updated',  'new_geos')
  )
  
  models <- data.frame(do.call('rbind', models))
  names(models) <- c('run_label', 'config_file', 'covs_file', 'cluster')
}


## Define function for submitting launch scripts
qsub_launch <- function(config_file, covs_file, cluster, run_date = NULL) {
  
  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  outdir <- "<<<< FILEPATH REDACTED >>>>"
  dir.create(outdir, showWarnings = F)
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)
  
  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), 
                '-pe multi_slot 7 -N', paste0('launch_', run_date),
                '-P proj_geo_nodes -l geos_node=TRUE' , 
                '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0('<<<< FILEPATH REDACTED >>>>', '/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0('<<<< FILEPATH REDACTED >>>>', '/breastfeeding/modeling/ebf/launch_ebf.R'),
                config_file, covs_file, run_date, cluster)
  system(qsub)
  Sys.sleep(10)
  return(run_date)
}

## Submit jobs and save job submission info
if (!is.null(resub_date)) {
  apply(models, 1, function(x) qsub_launch(x[2], x[3], x[4], x[5]))
  
} else {
  models$run_date <- apply(models, 1, function(x) qsub_launch(x[2], x[3], x[4]))
  write.table(models, file = paste0("<<<< FILEPATH REDACTED >>>>", '/models_submitted_', format(Sys.Date(), '%Y_%m_%d'), '.csv'),
              sep = ",", row.names = F, append = T)
}