# rerake 2017_stage_1_rerun using smaller tolerance
# using rake_lri_cell_pred_experiment
rm(list = ls())

#set arguments
Regions <- 'cssa'
rake_measures <- 'prevalence'
agg_measures <- 'prevalence'
run_date <- '2017_stage_1_updates'
rakecores <- 50
makeholdouts <- FALSE
year_list <- c(2000:2017)
indicator <- 'has_lri'
indicator_group <- 'lri'
modeling_shapefile_version <- '2018_08_28'
raking_shapefile_version <- '2018_12_04'
fun_tol_try <- c(1e-10, 1e-8, 1e-7)
outputdir <- '<<<< FILEPATH REDACTED >>>>'
pop_measure <- 'a0004t'
use_geos_nodes <- FALSE
core_repo       <- '<<<< FILEPATH REDACTED >>>>'

# set variables
nodes <- ''
user <- '<<<< USERNAME REDACTED >>>>'

## drive locations
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# set repo
setwd('<<<< FILEPATH REDACTED >>>>')

#source custom functions
source('submit_aggregation_script_experiment.R')

#### RAKING
# Define nodes
if (nodes == 'geos') {
  proj <- '-P proj_geo_nodes_dia -l gn=TRUE'
  r_shell <- '<<<< FILEPATH REDACTED >>>>'
} else {
  proj <- '-P proj_geospatial'
  r_shell <- '<<<< FILEPATH REDACTED >>>>'
}

# Launch child scripts
for (fun_tol in fun_tol_try){
  for (reg in Regions) {
    for (raketo in rake_measures){
      jname <- paste0(reg, '_', raketo, '_rake')
      sys.sub <- paste0('qsub ',proj,paste0(' -e <<<< FILEPATH REDACTED >>>>',
                                            ' -o <<<< FILEPATH REDACTED >>>> '),
                        '-cwd -N ', jname, ' ', '-pe multi_slot ', rakecores, ' ')
      script <- 'rake_lri_cell_pred_experiment.R'
      args <- paste(raketo, reg, outputdir, makeholdouts, year_list, indicator, modeling_shapefile_version, raking_shapefile_version, fun_tol)
      system(paste(sys.sub, r_shell, 1, script, args))
    }
  }
}

### AGGREGATION
submit_aggregation_script_experiment (indicator       = indicator,
                                      indicator_group = indicator_group,
                                      run_date        = run_date,
                                      raked           = c(T),
                                      pop_measure     = pop_measure,
                                      overwrite       = T,
                                      ages            = 0, # Note: can take vector of ages
                                      holdouts        = 0,
                                      regions         = Regions,
                                      corerepo        = '<<<< FILEPATH REDACTED >>>>',
                                      log_dir         = paste0(sharedir, "/output/", run_date),
                                      geo_nodes       = as.logical(use_geos_nodes),
                                      singularity     = "default",
                                      slots           = 10,
                                      modeling_shapefile_version = modeling_shapefile_version,
                                      raking_shapefile_version = raking_shapefile_version,
                                      measures        = agg_measures,
                                      fun_tol_try     = fun_tol_try)
