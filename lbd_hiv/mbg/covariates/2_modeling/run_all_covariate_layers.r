####################################################################################################
## Launch jobs for all covariate models
####################################################################################################

rm(list = ls())
# Load packages
# Set repo
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")
setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
mbg_setup(package_list = package_list, repos = c(indic_repo, core_repo))

## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))

} else {
  models <- rbind(
    # Base model
    data.table("male_circumcision",      c("CopyOfdefault_pc"),  "a1549m",  "none"),
    data.table("sti_symptoms",           c("CopyOfdefault_pc"),  "a1549t",  "none"),
    data.table("condom_last_time",       c("CopyOfdefault_pc"),  "a1549t",  "_BOTH"),
    data.table("had_intercourse",        c("CopyOfdefault_pc"),  "a1549f",   "_WN"),
    data.table("in_union",               c("CopyOfdefault_pc"),  "a1549t",  "_conservative"),
    data.table("multiple_partners_year", c("CopyOfdefault_pc"),  "a1549f",   "_WN"),
    data.table("multiple_partners_year", c("CopyOfdefault_pc"),  "a1549m",  "_MN"),
    data.table("partner_away",           c("CopyOfdefault_pc"),  "a1549f",  "none"))
  names(models) <- c('indicator', 'prior', 'population', 'data_tag')

}

## Define function for submitting launch scripts
#CHANGE BACK TO 7 FOR QSUB
qsub_launch <- function(indicator, prior, population, data_tag, run_date = NULL) {

  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  outdir <- paste0("<<<< FILEPATH REDACTED >>>>")
  dir.create(outdir, showWarnings = F)
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)

  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), '-l m_mem_free=15G -l fthread=5 -N', paste0('launch_', run_date),
                '-P proj_geo_nodes_hiv -q geospatial.q -l h_rt=08:00:00:00 -l archive=TRUE',
                '-v sing_image=<<<< FILEPATH REDACTED >>>>singularity-images/lbd/releases/lbd_full_20200128.simg -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
                paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_hiv/mbg/covariates/2_modeling/launch.r'),
                indicator, prior, population, data_tag, run_date)

  system(qsub)
  Sys.sleep(10)
  return(run_date)
}

## Submit jobs and save job submission info
if (!is.null(resub_date)) {
  apply(models, 1, function(x) qsub_launch(x[1], x[2], x[3], x[4]))

} else {
  models$run_date <- apply(models, 1, function(x) qsub_launch(x[1], x[2], x[3], x[4]))
  write.table(models, file = paste0("<<<< FILEPATH REDACTED >>>>"),
              sep = ",", row.names = F, append = T)

}

