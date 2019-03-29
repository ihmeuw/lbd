####################################################################################################
## Launch jobs for all covariate models
####################################################################################################

rm(list = ls())

# Load packages
# Set repo
core_repo  <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- "<<<< FILEPATH REDACTED >>>>/lbd_hiv/"
setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))

## Allow for re-submission with old run_dates (if resub_date is NULL, this is a new submission)
resub_date <- NULL

## List all models to run
if (!is.null(resub_date)) {
  models <- read.csv(paste0('<<<< FILEPATH REDACTED >>>>/models_submitted_', resub_date, '.csv'))

} else {
  models <- rbind(
    # Base model
    data.table("male_circumcision",      c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "a1549m",  "none",           "new_geos"),
    data.table("sti_symptoms",           c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "a1549t",  "none",           "old_geos"),
    data.table("condom_last_time",       c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "a1549t",  "_BOTH",          "new_geos"),
    data.table("had_intercourse",        c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "wocba",   "_WN",            "old_geos"),
    data.table("in_union",               c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "a1549t",  "_conservative",  "new_geos"),
    data.table("multiple_partners_year", c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "wocba",   "_WN",            "old_geos"),
    data.table("multiple_partners_year", c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "a1549m",  "_MN",            "new_geos"),
    data.table("partner_away",           c("default_pc", "default_pc_strict_rho", "loose_pc", "loose_pc_strict_rho", "inla_default", "inla_default_strict_rho"),  "wocba",  "none",            "new_geos"))
  names(models) <- c('indicator', 'prior', 'population', 'data_tag', 'cluster')

}

## Define function for submitting launch scripts
qsub_launch <- function(indicator, prior, population, data_tag, cluster, run_date = NULL) {

  if (is.null(run_date)) run_date <- format(Sys.time(), '%Y_%m_%d_%H_%M_%S')
  outdir <- paste0('<<<< FILEPATH REDACTED >>>>')
  dir.create(outdir, showWarnings = F)
  dir.create(paste0(outdir, 'errors'), showWarnings = F)
  dir.create(paste0(outdir, 'output'), showWarnings = F)

  qsub <- paste('qsub -e', paste0(outdir, 'errors'), '-o', paste0(outdir, 'output'), '-pe multi_slot 7 -N', paste0('launch_', run_date),
                '-P proj_geo_nodes_hiv -l geos_node=TRUE',
                '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
                '<<<< FILEPATH REDACTED >>>>/lbd_core/mbg_central/share_scripts/shell_sing.sh',
                '<<<< FILEPATH REDACTED >>>>/lbd_hiv/mbg/covariates/2_modeling/launch.r',
                indicator, prior, population, data_tag, cluster, run_date)
  system(qsub)
  Sys.sleep(10)
  return(run_date)
}

## Submit jobs and save job submission info
if (!is.null(resub_date)) {
  apply(models, 1, function(x) qsub_launch(x[1], x[2], x[3], x[4], x[5], x[6]))

} else {
  models$run_date <- apply(models, 1, function(x) qsub_launch(x[1], x[2], x[3], x[4], x[5]))
  write.table(models, file = paste0('<<<< FILEPATH REDACTED >>>>/models_submitted_', format(Sys.Date(), '%Y_%m_%d'), '.csv'),
              sep = ",", row.names = F, append = T)

}
