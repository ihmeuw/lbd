####################################################################################################
## Launch HIV prevalence (hiv_test) models for Africa
####################################################################################################

###############################################################################
## SETUP
###############################################################################

## Clear environment
rm(list=ls())

## Set indicator
indicator_group <- 'hiv'
indicator       <- 'hiv_test'

## Set repo
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

setwd(core_repo)

## Load libraries and MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = package_list, repos = core_repo)

## Load passed arguments
config_file <- commandArgs()[4]
covs_file <- commandArgs()[5]
run_date <- commandArgs()[6]
cluster <- commandArgs()[7]

## Record git status info
outputdir <- paste("<<<< FILEPATH REDACTED >>>>")
record_git_status(core_repo, indic_repo, show_diff = T, file = paste0(outputdir, "git_at_launch.txt"))

## Load and check config
config <- load_config(repo            = indic_repo,
                      indicator       = "",
                      indicator_group = "",
                      config_name     = paste0('mbg/hiv_test/2_modeling/', config_file),
                      covs_name       = paste0('mbg/hiv_test/2_modeling/', covs_file))
check_config()

if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (length(summstats) == 1 & grepl("c(", summstats, fixed =T)) summstats <- eval(parse(text=summstats))

## Make sure pc prior options are correct class
if (is.character(pc_prior)) pc_prior <- as.logical(pc_prior)
if (pc_prior & is.character(sigma_upper_tail) | pc_prior & is.character(range_lower_tail)){
  range_lower_tail <- as.numeric(range_lower_tail)
  sigma_upper_tail <- as.numeric(sigma_upper_tail)
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
## Make Holdouts
###############################################################################
if (makeholdouts) {

  set.seed(98121)

  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        simple      = NULL,
                        removeyemen = TRUE,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        yl          = year_list,
                        use_share   = as.logical(use_share))

  # add in location information
  df <- merge_with_ihme_loc(df)

  # ANC is not included in the testing data, so remove then add back after making folds
  anc <- df[type == "ANC",]
  df <- df[type != "ANC",]

  # make a list of dfs for each region, were each list element is a fold (testing dataset)
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

  # add back in ANC data
  anc[, fold := 99]
  for (r in Regions) stratum_ho[[paste0('region__', r)]] <- rbind(stratum_ho[[paste0('region__', r)]], anc[region == r,], fill = T)
}

###############################################################################
## Launch Parallel Script
###############################################################################

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if (makeholdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)
loopvars$submitted <- F

if (setequal(c('essa_sdn', 'wssa', 'cssa', 'sssa'), Regions)) {
  loopvars$Var1 <- factor(loopvars$Var1, levels = c('essa_sdn', 'wssa', 'cssa', 'sssa')) # make sure bigger models get submitted first
  loopvars <- loopvars[order(loopvars$Var1, loopvars$Var3), ]
  loopvars$Var1 <- as.character(loopvars$Var1)
}

## Loop over them, save images and submit qsubs
for (i in 1:nrow(loopvars)) {

  message(paste(loopvars[i, 2], as.character(loopvars[i, 1]), loopvars[i, 3]))

  # check that holdouts exist
  if (loopvars[i, 3] > 0) {
    count <- sum(stratum_ho[[paste0('region__', loopvars[i, 1])]]$fold %in% loopvars[i, 3])
    if (count == 0) next
  }

  # set slots (based on observed memory use)
  slots <- c(cssa = 5, essa_sdn = 9, sssa = 5, wssa = 8)[as.character(loopvars[i, 1])]
  if (is.na(slots)) slots <- 9
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
  loopvars$submitted[i] <- T
}

## Check to make sure models are done before continuing
loopvars <- loopvars[loopvars$submitted == T, ]
waitformodelstofinish(lv = cbind(as.character(loopvars[, 1]), loopvars[, 3]), sleeptime = 60)

###############################################################################
## Make summary metrics
###############################################################################

## Combine csv files
csvs <- list.files(outputdir, pattern = "input_data_(.*).csv", full.names = T)
csv_master <- rbindlist(lapply(csvs, fread))
csv_master[, V1 := NULL]
write.csv(csv_master, file=paste0(outputdir, '/input_data.csv'))

## Launch predictive validity script
if (!is.logical(use_geos_nodes)) {
  flag <- '-q all.q@@c2-nodes'
} else {
  if (new_geos_nodes) {
    flag <- '-l geos_node=TRUE -q geospatial.q@@intel_gold_6130_756gb'
  } else {
    flag <- '-l geos_node=TRUE -q geospatial.q@@intel_e5_2660v4_1024gb'
  }
}

qsub <- paste('qsub -e', paste0(outputdir, 'errors'), '-o', paste0(outputdir, 'output'),
              '-P', project, flag,
              '-pe multi_slot', ceiling(10 * slot_multiplier), '-N', paste0('pv_', run_date),
              '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
              paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
              paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_hiv/mbg/functions/predictive_validity.r'),
              indicator, indicator_group, run_date, core_repo)
system(qsub)

###############################################################################
## Post-Estimation, Aggregate, and PLHIV
###############################################################################

## Save strata for Shiny to use in producing aggregated fit statistics
strata <- unique(as.character(loopvars[, 1]))
dir.create(paste0(outputdir, '/fit_stats'))
save(strata, file = paste0(outputdir, '/fit_stats/strata.RData'))

## Load GBD Estimates for this indicator which will be used in raking
source(paste0(indic_repo, "mbg/functions/fractional_raking_functions_hiv.r"))
rake_subnational = F
age_group <- 24
sex_id <- 3
connector <- get_gbd_locs(rake_subnational = rake_subnational, reg = Regions)

source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
gbd <- get_outputs(topic = "cause",
                   version = "latest",
                   gbd_round_id = 5,
                   cause_id = 298,
                   measure_id = 5,
                   metric_id = 3,
                   age_group_id = 24,
                   location_id = connector$location_id,
                   year_id = year_list)
gbd <- gbd[, list(name = location_id, year = year_id, mean = val)]

## Prepare for parallel post-estimation - save file with objs to re-load in child processes
prep_postest(indicator       = indicator,
             indicator_group = indicator_group,
             run_date        = run_date,
             save_objs       = c("core_repo", "indic_repo", "gbd", "year_list", "summstats", "rake_transform", "pop_measure", "rake_subnational", "age_group", "sex_id"))

## Parallelized post-estimation over region
for (s in strata) {

  # set slots (based on observed memory use)
  slots <- c(cssa = 16, essa_sdn = 23, sssa = 13, wssa = 23)[s]
  if (is.na(slots)) slots <- 23
  slots <- ceiling(slots * slot_multiplier)

  # make qsub string
  qsub <- make_qsub_postest(code           = "fractional_raking_implementation_hiv",
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

## Check to make sure post-est done before continuing
waitformodelstofinish(lv = cbind(strata, 0), sleeptime=60)

## Combine post est stuff across regions and save needed outputs
source(paste0(indic_repo, "mbg/functions/accomodating_counts_functions_hiv.r"))
post_load_combine_save_counts(regions    = strata,
                              summstats  = c('mean','cirange','upper','lower'),
                              raked      = c('raked','unraked', 'raked_c'),
                              rf_table   = TRUE,
                              run_summ   = TRUE,
                              indic      = indicator,
                              ig         = indicator_group,
                              sdir       = sharedir)

## Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)

combine_aggregation_counts(rd = run_date,
                           indic = indicator,
                           ig = indicator_group,
                           ages = 0,
                           regions = strata,
                           holdouts = 0,
                           raked = c("unraked", "raked", "raked_c"),
                           dir_to_search = NULL,
                           delete_region_files = T,
                           merge_hierarchy_list = F)

summarize_admins_counts(ind = indicator,
                        ig = indicator_group,
                        summstats = c("mean", "lower", "upper", "cirange"),
                        raked    = c("raked", "raked_c", "unraked"),
                        ad_levels = c(0,1,2),
                        file_addin = NULL)

## Create change uncertainty estimates
source(paste0(indic_repo, "mbg/functions/admin_uncertainty_draws.r"))
make_admin_change_uncertainty_tables(run_date, indicator, indicator_group, strata, core_repo, outputdir)
make_pixel_difference_uncertainty_rasters(run_date, indicator, indicator_group, strata, summstats, core_repo, outputdir)
