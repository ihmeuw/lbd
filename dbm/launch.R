## ##########################################################################
## ##########################################################################
## This script is a generic launch script used for overweight and wasting indicators
## ##########################################################################
## ##########################################################################


###############################################################################
## Launch overweight and wasting prevalence in Low and Middle Income Countries (LMICs)
###############################################################################

## clear environment
rm(list=ls())

###############################################################################
## SETUP
###############################################################################


##Set indicator 
indicator_group <- 'child_growth_failure'
indicator       <- 'indicator' 


#Set repo
core_repo  <- "<<<< FILEPATH REDACTED >>>>"
indic_repo <- "<<<< FILEPATH REDACTED >>>>"
sharedir   <- "<<<< FILEPATH REDACTED >>>>"
commondir  <- "<<<< FILEPATH REDACTED >>>>"

setwd(core_repo)


# Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo) 

## Load passed arguments
config_file <- commandArgs()[4]
covs_file   <- commandArgs()[5]
run_date    <- commandArgs()[6]
cluster     <- commandArgs()[7]

## Create output folder with the run_date
outputdir <- "<<<< FILEPATH REDACTED >>>>"
dir.create(outputdir)


## Load and check config
config <- load_config(repo            = indic_repo,
                      indicator_group = "",
                      indicator       = "",
                      config_name     = paste0('modeling/', indicator, '/', config_file),
                      covs_name       = paste0('modeling/', indicator, '/', covs_file))

check_config()


if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if(length(summstats) == 1 & grepl(",", summstats)) summstats <- eval(parse(text=summstats))


## Get cluster-related arguments
stopifnot(cluster %in% c('prod', 'old_geos', 'new_geos'))
if (cluster == 'old_geos') {
  project <- 'proj_geo_nodes_cgf'
  use_geos_nodes <- T
  new_geos_nodes <- F
  slot_multiplier <- 1
}
if (cluster == 'new_geos') {
  project <- 'proj_geo_nodes_cgf'
  use_geos_nodes <- T
  new_geos_nodes <- T
  slot_multiplier <- 1.5
}
if (cluster == 'prod') {
  project <- 'proj_geospatial_cgf'
  use_geos_nodes <- F
  new_geos_nodes <- NULL
  slot_multiplier <- 2.0
}

###############################################################################
## Make Holdouts
###############################################################################

if(makeholdouts){
  
  set.seed(98121)
  
  # Load the full input data
  df <- load_input_data(indicator   = indicator,
                        simple      = NULL,
                        removeyemen = TRUE,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        yl          = year_list,
                        use_share   = as.logical(use_share))
  
  # Add in location information
  df <- merge_with_ihme_loc(df, shapefile_version = modeling_shapefile_version)
  
  # Make a list of dfs for each region, were each list element is a fold (testing dataset)
  if (holdout_strategy == "qt") {
    stratum_ho <- make_folds(data       = df,
                             n_folds    = as.numeric(n_ho_folds),
                             spat_strat = 'qt',
                             temp_strat = 'prop_comb',
                             strat_cols = 'region',
                             ts         = as.numeric(ho_ts),
                             mb         = as.numeric(ho_mb))
    
  } else if (holdout_strategy == "nids") {
    stratum_ho <- make_folds(data       = df,
                             n_folds    = as.numeric(n_ho_folds),
                             spte_strat = 'nids',
                             strat_cols = 'region',
                             ss_col     = ss_col,
                             yr_col     = yr_col)
    
  } else {
    stop("holdout strategy not recognized")
  }
}

###############################################################################
## Launch Parallel Script and Define SINGULARITY image
###############################################################################

message("###############################################################################")
message("Launch Parallel Script")
message("###############################################################################")

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(makeholdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)

# Specify number slots per region (based on observed memory use)
slot_list <- c(8, 5, 5, 9, 8, 5, 5, 8, 5, 10, 10, 10)

names(slot_list) <- c('noaf', 'cssa', 'mide', 'caca', 'wssa', 'ansa', 'sssa',
                    'trsa', 'essa', 'eaas', 'south_asia_cgf','seas_ocea_cgf')


## Loop over them, save images and submit qsubs
for(i in 1:nrow(loopvars)){
  
  # Set slots (based on observed memory use)
  slots <- slot_list[[loopvars[i,1]]]
  if (!as.logical(use_geos_nodes)) slots <- ceiling(slots * slot_multiplier)
  
  message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3], slots))
  
  # Make a qsub string
  qsub <- make_qsub_share(age           = loopvars[i,2],
                          reg           = as.character(loopvars[i,1]),
                          holdout       = loopvars[i,3],
                          test          = as.logical(test),
                          indic         = indicator,
                          saveimage     = TRUE,
                          memory        = 10,
                          cores         = slots,
                          geo_nodes     = use_geos_nodes,
                          proj          = project,
                          queue         = "geospatial.q",
                          run_time      = "04:00:00:00",
                          addl_job_name = eval(parse(text = jn)), ## from config
                          singularity   = 'default',
                          singularity_opts  = list(SET_OMP_THREADS=4, SET_MKL_THREADS=1))
  
  system(qsub)
}

## Check to make sure models are done before continuing
waitformodelstofinish(lv = cbind(as.character(loopvars[, 1]), loopvars[, 3]), sleeptime = 60)

###############################################################################
## Post-Estimation
###############################################################################

message("###############################################################################")
message("Launch Post-Estimation Script")
message("###############################################################################")

## Save strata for Shiny to use in producing aggregated fit statistics

strata <- unique(as.character(loopvars[, 1]))
dir.create(paste0(outputdir, '/fit_stats'))
save(strata, file = paste0(outputdir, '/fit_stats/strata.RData'))

# Load GBD Estimates for this indicator which will be used in raking
gbd <- fread(paste0(sharedir, "/[indicator]_gbd_2017.csv"))

## Prepare for parallel post-estimation - save file with objs to re-load in child processes
prep_postest(indicator       = indicator,
             indicator_group = indicator_group,
             run_date        = run_date,
             save_objs       = c("core_repo", "indic_repo", "gbd", "year_list", "summstats", 
                                 "rake_transform", "pop_measure", "rake_subnational",
                                 "shapefile_version", "config", "pop_release"))

# Specify number slots per region (based on observed memory use)
slot_list <- c(8, 5, 5, 9, 8, 5, 5, 8, 5, 10, 10, 10)

names(slot_list) <- c('noaf', 'cssa', 'mide', 'caca', 'wssa', 'ansa', 'sssa',
                    'trsa', 'essa', 'eaas', 'south_asia_cgf','seas_ocea_cgf')

postest_script <- "postest_nonfrax_script"


for (s in strata) {
  
  # set slots (based on observed memory use)
  slots <- slot_list[[s]]
  if (!as.logical(use_geos_nodes)) slots <- ceiling(slots * slot_multiplier)
  
  qsub <- make_qsub_postest(code          = "postest_script",
                            stratum       = s,
                            log_location  = 'sharedir',
                            memory        = 10,
                            cores         = slots,
                            proj          = project,
                            singularity   = "default",
                            queue         = "geospatial.q",
                            run_time      = "04:00:00:00",
                            proj          =  project,
                            cores         = 10,
                            subnat_raking = subnational_raking,
                            modeling_shapefile_version = modeling_shapefile_version,
                            raking_shapefile_version   = raking_shapefile_version,)
  system(qsub)
}

## Check to make sure post-est done before continuing
waitformodelstofinish(lv = cbind(strata, 0), sleeptime=60)

## Combine post est stuff across regions and save needed outputs
post_load_combine_save(summstats  = summstats)

## Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)

###############################################################################
## Aggregate to admin2, admin1, and national levels
###############################################################################

message("###############################################################################")
message("Launch Aggregation Script")
message("###############################################################################")

# Specify number slots per region (based on observed memory use)
slot_list <- c(8, 5, 5, 9, 8, 5, 5, 8, 5, 10, 10, 10)

names(slot_list) <- c('noaf', 'cssa', 'mide', 'caca', 'wssa', 'ansa', 'sssa',
                    'trsa', 'essa', 'eaas', 'south_asia_cgf','seas_ocea_cgf')

for (s in strata) {
  
  # set slots (based on observed memory use)
  slots <- slot_list[[s]]
  if (!as.logical(use_geos_nodes)) slots <- slots * 2
  
  # make and submit qsub string
  submit_aggregation_script(indicator       = indicator,
                            indicator_group = indicator_group,
                            run_date        = run_date,
                            raked           = c(TRUE, FALSE),
                            pop_measure     = pop_measure,
                            overwrite       = T,
                            ages            = 0, # Note: can take vector of ages
                            holdouts        = 0,
                            regions         = s,
                            corerepo        = core_repo,
                            log_dir         = outputdir,
                            slots           = slots,
                            geo_nodes       = use_geos_nodes,
                            use_c2_nodes    = FALSE,
                            proj            = project,
                            singularity     = 'default',
                            run_time        = "04:00:00:00",
                            queue           = "geospatial.q",
                            modeling_shapefile_version = modeling_shapefile_version,
                            raking_shapefile_version   = raking_shapefile_version)
  
}

waitforaggregation(rd       = run_date,
                   indic    = indicator,
                   ig       = indicator_group,
                   ages     = 0,
                   regions  = strata,
                   holdouts = 0,
                   raked    = c(T, F))

message("###############################################################################")
message("Combine Aggregation")
message("###############################################################################")

combine_aggregation(rd       = run_date,
                    indic    = indicator,
                    ig       = indicator_group,
                    ages     = 0,
                    regions  = strata,
                    holdouts = 0,
                    raked    = c(TRUE, FALSE))

message("###############################################################################")
message("Summarize Admins")
message("###############################################################################")

summarize_admins(summstats = summstats,
                 ad_levels = c(0, 1, 2),
                 raked     = c(TRUE, FALSE))