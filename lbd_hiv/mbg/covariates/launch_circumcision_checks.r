####################################################################################################
## Launch circumcision checks
####################################################################################################

indicator<- 'circumcision'

## Define function for submitting launch scripts
qsub_launch <- function(indicator, run_date = NULL) {
  
  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  outdir <- paste0('<<<< FILEPATH REDACTED >>>>')
  dir.create(outdir, showWarnings = F)
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)
  
  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), 
                '-l m_mem_free=10G -l fthread=5 -N', 
                paste0('mc_checks_', run_date),
                '-P proj_geo_nodes -q geospatial.q -l h_rt=00:00:30:00 -l archive=TRUE',
                '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0('/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0('/lbd_hiv/mbg/covariates/1_collapse/2_', indicator, '_collapse_checks.r'))
  
  system(qsub)
  return(run_date)
}

qsub_launch(indicator)
