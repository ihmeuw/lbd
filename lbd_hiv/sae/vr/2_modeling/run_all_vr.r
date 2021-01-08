####################################################################################################
## Launch jobs for VR models 
####################################################################################################

rm(list = ls())

library(data.table)
## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0('<<<< FILEPATH REDACTED >>>>'))
  
} else {
  
  models <- rbind(
    data.table('col', c(6, "6_c")),
    data.table('ecu', c(6, "6_c")),
    data.table('gtm', c(6, "6_c")),
    data.table('bra', c(6, "6_c")),
    data.table('mex', c(6, "6_c")),
    data.table('cri', c(6)))
  names(models) <- c('country', 'model')
  
}

## Define function for submitting launch scripts
qsub_launch <- function(country, model, run_date = NULL, resub = F) {
  
  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  
  outdir <- paste0('<<<< FILEPATH REDACTED >>>>')
  dir.create(outdir, showWarnings = F)
  # Error and output are created later, so not needed
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)
  
  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), '-l m_mem_free=15G -l fthread=5 -N', paste0('launch_', run_date),
                '-P proj_geo_nodes -q geospatial.q -l h_rt=08:00:00:00 -l archive=TRUE',
                '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0('<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0('<<<< FILEPATH REDACTED >>>>/lbd_hiv/sae/vr/2_modeling/launch.r'),
                country, model, run_date, resub)
  
  system(qsub)
  Sys.sleep(10)
  return(run_date)
}

## Submit jobs and save job submission info
if (!is.null(resub_date)) {
  apply(models, 1, function(x) qsub_launch(x[1], x[2], x[3], T))
} else {
  models$run_date <- apply(models, 1, function(x) qsub_launch(x[1], x[2]))
  write.table(models, file = paste0('<<<< FILEPATH REDACTED >>>>'),
              sep = ",", row.names = F, append = T)
}

