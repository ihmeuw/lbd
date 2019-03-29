###########################################################################################
###########################################################################################
## Launch all models
##
###########################################################################################
###########################################################################################

###########################################################################################
## SETUP
###########################################################################################

## Clear environment
rm(list=ls())

## Write message for clarity
message(commandArgs()[4])
message(commandArgs()[5])
message(commandArgs()[6])
message(commandArgs()[7])
message(commandArgs()[8])
message(commandArgs()[9])

## Set indicator
indicator_group <- 'hiv'
indicator       <- commandArgs()[4]
prior           <- commandArgs()[5]
population      <- commandArgs()[6]
data_tag        <- commandArgs()[7]
cluster         <- commandArgs()[8]
run_date        <- commandArgs()[9]

# set arguments for this run and user.
core_repo  <- "<<<< FILEPATH REDACTED >>>>/lbd_core/"
indic_repo <- "<<<< FILEPATH REDACTED >>>>/lbd_hiv/"
setwd(core_repo)

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = c(core_repo, indic_repo))

## Record git status info
outputdir <- '<<<< FILEPATH REDACTED >>>>'

#load config
config <- load_config(repo            = indic_repo,
                      indicator_group = "",
                      indicator       = "",
                      config_name     = paste0('mbg/covariates/2_modeling/config_', prior),
                      covs_name       = 'mbg/covariates/2_modeling/cov_list')

## Ensure you have defined all necessary settings in your config
check_config()

## Create a few objects from the config file loaded above
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (length(summstats) == 1 & grepl(",", summstats)) summstats <- eval(parse(text=summstats))
rake_subnational <- F

## Make sure pc prior options are correct class
if (is.character(pc_prior)) pc_prior <- as.logical(pc_prior)
if (pc_prior & is.character(sigma_upper_tail) | pc_prior & is.character(range_lower_tail)){
  range_lower_tail <- as.numeric(range_lower_tail)
  sigma_upper_tail <- as.numeric(sigma_upper_tail)
}

## Make sure correct population is specified
pop_measure <- population
if (data_tag == "none"){
  withtag <- FALSE
  datatag <- ""
} else {
  withtag <- TRUE
  datatag <- data_tag
}

## Get cluster-related arguments
stopifnot(cluster %in% c('prod', 'old_geos', 'new_geos'))
if (cluster == 'old_geos') {
  project <- 'proj_geo_nodes_hiv'
  use_geos_nodes <- T
  new_geos_nodes <- F
  slot_multiplier <- 1
}
if (cluster == 'new_geos') {
  project <- 'proj_geo_nodes_hiv'
  use_geos_nodes <- T
  new_geos_nodes <- T
  slot_multiplier <- 1.5
}
if (cluster == 'prod') {
  project <- 'proj_geospatial_hiv'
  use_geos_nodes <- F
  new_geos_nodes <- NULL
  slot_multiplier <- 2.0
}

###############################################################################
## Launch Parallel Script
###############################################################################

## Make sure bigger regions are submitted first
loopvars <- expand.grid(Regions, 0, 0)
region_order <- c("essa_sdn", "wssa", "cssa", "sssa")
loopvars <- loopvars[order(match(loopvars$Var1, region_order)),]

# Loop over them, save images and submit qsubs
for (i in 1:nrow(loopvars)) {

  message(paste(loopvars[i, 2], as.character(loopvars[i, 1]), loopvars[i, 3]))

  # set slots (based on observed memory use)
  slots <- c(cssa = 5, essa_sdn = 9, sssa = 5, wssa = 8)[as.character(loopvars[i, 1])]

  if (is.na(slots)) slots <- 7
  slots <- ceiling(slots * slot_multiplier)

  # make a qsub string
  qsub <- make_qsub_share(age            = loopvars[i, 2],
                          reg            = as.character(loopvars[i, 1]),
                          holdout        = loopvars[i, 3],
                          test           = as.logical(test),
                          indic          = indicator,
                          saveimage      = TRUE,
                          memory         = 10,
                          cores          = slots,
                          proj           = project,
                          geo_nodes      = use_geos_nodes,
                          use_c2_nodes   = TRUE,
                          new_geos_nodes = new_geos_nodes,
                          singularity    = 'default')

  system(qsub)
}


## check to make sure models are done before continuing
waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

###############################################################################
## Post-Estimation
###############################################################################

## Save strata for Shiny to use in producing aggregated fit statistics
strata <- unique(as.character(loopvars[, 1]))
dir.create(paste0(outputdir, '/fit_stats'))
save(strata, file = paste0(outputdir, '/fit_stats/strata.RData'))

## Load GBD Estimates for this indicator which will be used in raking
gbd <- NULL

#this is so we can grab these objects in temp_estimation_folder, all of our processes will need some of these objects.
prep_postest(indicator        = indicator,
             indicator_group  = indicator_group,
             run_date         = run_date,
             save_objs        = c("core_repo", "gbd", "year_list", "summstats", "pop_measure")) # can add gbd and rake_transform

## Parallelized post-estimation over region

for (s in strata) {

  # set slots (based on observed memory use)
  slots <- c(cssa = 16, essa_sdn = 23, sssa = 13, wssa = 23)[s]
  if (is.na(slots)) slots <- 23
  slots <- ceiling(slots * slot_multiplier)

  # make qsub string
  qsub <- make_qsub_postest(code           = "postest_script",
                            stratum        = s,
                            log_location   = 'sharedir',
                            memory         = 10,
                            cores          = slots,
                            proj           = project,
                            geo_nodes      = use_geos_nodes,
                            use_c2_nodes   = TRUE,
                            new_geos_nodes = new_geos_nodes,
                            singularity    = 'default')
  system(qsub)
}

## check to make sure post-est done before continuing
waitformodelstofinish(lv = cbind(strata, 0), sleeptime=60)

## Combine post est stuff across regions and save needed outputs
#This makes a combined brick of all of the regions

#make sure strata is converted back into list

post_load_combine_save(summstats = summstats,
                       raked     = 'unraked',
                       rf_table  = F,
                       run_summ  = F,
                       regions   = strata)

# Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)


# Combine csv files of input data
csvs <- list.files(outputdir, pattern = "input_data_(.*).csv", full.names = T)
csv_master <- rbindlist(lapply(csvs, fread))
csv_master[, V1 := NULL]
write.csv(csv_master, file=paste0(outputdir, '/input_data.csv'))
message('Done!')
