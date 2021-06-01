####################################################################################################
## Launch HIV prevalence (hiv_adult_prev) models for Africa
####################################################################################################

###############################################################################
## SETUP
###############################################################################

## Clear environment
rm(list = ls())

## Set indicator
indicator_group <- 'hiv'
indicator       <- 'hiv_adult_prev'

## Set repo
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")

setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo,  '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
mbg_setup(package_list = package_list, repos = core_repo)

## Load passed arguments
config_name <- commandArgs()[4]
covs_name <- commandArgs()[5]
run_date <- commandArgs()[6]

## Record git status info
outputdir <- paste("<<<< FILEPATH REDACTED >>>>", indicator_group, indicator, 'output', run_date, '', sep='/')
record_git_status(core_repo, indic_repo, show_diff = T, file = paste0(outputdir, "git_at_launch.txt"))

## Load and check config
config <- set_up_config(repo            = indic_repo,
                        core_repo       = core_repo,
                        indicator       = "",
                        indicator_group = "",
                        config_file     = paste0(indic_repo, 'mbg/', indicator, '/2_modeling/', config_name,   '.csv'),
                        covs_file       = paste0(indic_repo, 'mbg/', indicator, '/2_modeling/', covs_name,   '.csv'),
                        run_tests       = TRUE)
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

##Set cluster environment for fair cluster
Sys.setenv(SGE_ENV=SGE_ENV)

## Get cluster project-related arguments
stopifnot(queue %in% c('long.q', 'geospatial.q'))
if (queue == 'geospatial.q') {
  project <- 'proj_geo_nodes_hiv'
  use_geos_nodes <- T
}
if (queue == 'long.q') {
  project <- 'proj_geospatial_hiv'
  use_geos_nodes <- F
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
  df <- merge_with_ihme_loc(df, shapefile_version = shapefile_version)

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
  loopvars$Var1 <- factor(loopvars$Var1, levels = c('essa_sdn', 'wssa', 'cssa', 'sssa')) # temporary hack to make sure bigger models get submitted first
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

  # set memory
  memory <- c(cssa = 180, essa_sdn = 260, sssa = 150, wssa = 270)[as.character(loopvars[i, 1])]
  if (is.na(memory)) memory <- 200

  # make a qsub string
  qsub <- make_qsub_share(age            = loopvars[i, 2],
                          reg            = as.character(loopvars[i, 1]),
                          holdout        = loopvars[i, 3],
                          test           = as.logical(test),
                          indic          = indicator,
                          saveimage      = TRUE,
                          memory         = memory,
                          cores          = 12,
                          proj           = project,
                          geo_nodes      = use_geos_nodes,
                          use_c2_nodes   = TRUE,
                          singularity = '"<<<< FILEPATH REDACTED >>>>"',
                          queue          = queue,
                          run_time       = run_time,
                          code           = NULL)

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

qsub <- paste('qsub -e', paste0(outputdir, 'errors'), '-o', paste0(outputdir, 'output'),
              '-P', project, paste0('-q ', queue), '-l fthread=10 -l m_mem_free=100G',
              "-l h_rt=0:05:00:00", '-N', paste0('pv_', run_date), '-l archive=TRUE',
              '-v sing_image="<<<< FILEPATH REDACTED >>>>"pkgs3.5.0_2_gcc7mklrstudioserver1.1.453.simg  -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
              paste0("<<<< FILEPATH REDACTED >>>>",'/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
              paste0("<<<< FILEPATH REDACTED >>>>",'/lbd_hiv/mbg/functions/predictive_validity.r'),
              indicator, indicator_group, run_date, core_repo, shapefile_version)
system(qsub)

###############################################################################
## Post-Estimation, Aggregate, and PLHIV
###############################################################################

## Save strata for Shiny to use in producing aggregated fit statistics
strata <- unique(as.character(loopvars[, 1]))
dir.create(paste0(outputdir, '/fit_stats'))
save(strata, file = paste0(outputdir, '/fit_stats/strata.RData'))

source("<<<< FILEPATH REDACTED >>>>/current/r/get_ids.R")
get_ids(table = "gbd_round")

## Load GBD Estimates for this indicator which will be used in raking
#source(paste0(indic_repo, "<<<< FILEPATH REDACTED >>>>"))
rake_subnational = T
age_group <- 24
sex_id <- 3
connector <- get_gbd_locs(rake_subnational = rake_subnational, reg = Regions, shapefile_version = shapefile_version)

source("<<<< FILEPATH REDACTED >>>>/current/r/get_outputs.R")
gbd <- get_outputs(topic = "cause",
                  version = "latest", # using the latest version instead of the best version on Sept 18, 2018. Revisit as other versions become available
                  gbd_round_id = 6,
                  cause_id = 298,
                  measure_id = 5,
                  metric_id = 3,
                  age_group_id = 24,
                  location_id = connector$location_id,
                  year_id = year_list)
gbd <- gbd[, list(name = location_id, year = year_id, mean = val)]

###########################################################################################################
###########################################################################################################
## Prepare for parallel post-estimation - save file with objs to re-load in child processes
prep_postest(indicator       = indicator,
             indicator_group = indicator_group,
             run_date        = run_date,
             save_objs       = c("core_repo", "indic_repo", "gbd", "year_list",
                                 "summstats", "rake_transform", "pop_measure",
                                 "rake_subnational", "age_group", "sex_id",
                                 "shapefile_version", "config", "pop_release"))

## Parallelized post-estimation over region
for (s in strata) {

  # set memory (based on observed memory use)
  memory <- c(cssa = 500, essa_sdn = 700, sssa = 450, wssa = 750)[s]
  if (is.na(memory)) memory <- 600

  # make qsub string
  qsub <- make_qsub_postest(code           = "postest_frax_script",
                            stratum        = s,
                            log_location   = 'sharedir',
                            memory         = memory,
                            cores          = 10,
                            run_time       = '01:15:00:00',
                            proj           = project,
                            geo_nodes      = use_geos_nodes,
                            use_c2_nodes   = TRUE,
                            singularity    = "<<<< FILEPATH REDACTED >>>>",
                            modeling_shapefile_version = shapefile_version,
                            raking_shapefile_version = shapefile_version)
  system(qsub)
}

## Check to make sure post-est done before continuing
waitformodelstofinish(lv = cbind(strata, 0), sleeptime = 60)

## Combine post est stuff across regions and save needed outputs
#source(paste0(indic_repo, "<<<< FILEPATH REDACTED >>>>"))
post_load_combine_save(regions    = strata,
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

combine_aggregation(rd = run_date,
                           indic = indicator,
                           ig = indicator_group,
                           ages = 0,
                           regions = strata,
                           holdouts = 0,
                           raked = c("unraked", "raked", "raked_c"),
                           dir_to_search = NULL,
                           delete_region_files = T,
                           merge_hierarchy_list = F)

summarize_admins(ind = indicator,
                        ig = indicator_group,
                        summstats = c("mean", "lower", "upper", "cirange"),
                        raked    = c("raked", "raked_c", "unraked"),
                        ad_levels = c(0,1,2),
                        file_addin = NULL)

## Create aggregated stackers
if (use_stacking_covs) {
  aggregate_stackers_admin0(indicator                 = indicator,
                            indicator_group           = indicator_group,
                            run_date                  = run_date,
                            age                       = 0,
                            holdout                   = 0,
                            regions                   = Regions,
                            year_list                 = year_list,
                            pop_measure               = pop_measure,
                            stackers_logit_transform  = stackers_in_transform_space,
                            predictor_logit_transform = TRUE,
                            shapefile_version         = shapefile_version,
                            results_file              = paste0(outputdir, "<<<< FILEPATH REDACTED >>>>"))
}

## Create change uncertainty estimates
source(paste0(indic_repo, "<<<< FILEPATH REDACTED >>>>"))
make_admin_change_uncertainty_tables(run_date, indicator, indicator_group, strata, core_repo, outputdir)
make_pixel_difference_uncertainty_rasters(run_date, indicator, indicator_group, strata, summstats, core_repo, outputdir)

###############################################################################
## Make maps and plots
###############################################################################

dir.create(paste0(outputdir, '/diagnostic_plots/'))

## Covariate importance plots
if (use_stacking_covs) {
  source(paste0(indic_repo, "mbg/functions/get_cov_weights.r"))
  get_cov_weights(indicator, indicator_group, run_date, Regions, paste0(outputdir, '/diagnostic_plots/'))
}

## Hyperparameters
source(paste0(indic_repo, "mbg/functions/plot_hyperparameters.r"))
plot_hyperparameters(indicator, indicator_group, run_date, 0, 0, paste0(outputdir, "/diagnostic_plots/hyperparameters.pdf"))

## Data-and-estimates maps
source(paste0(indic_repo, "mbg/", indicator,"/3_post_estimation/data_and_estimates_maps.r"))
try(data_and_estimates_maps_simplified(indicator, indicator_group, run_date, anc_correction != 'none'))
try(data_and_estimates_maps(indicator, indicator_group, run_date, anc_correction != 'none'))

## Aggregated data-and-estimates plots
source(paste0(indic_repo,  "mbg/", indicator, "/3_post_estimation/make_ts_plots.r"))
try(make_ts_plots(indicator, indicator_group, run_date, Regions, paste0(outputdir, '/diagnostic_plots/'), shapefile_version, counts = F))

## National data-and-estimates w/ stackers plots
source(paste0(indic_repo, "mbg/", indicator, "/3_post_estimation/admin0_data_and_estimates_plots.r"))
try(admin0_data_and_estimates_plots(run_date))

## Levels and differences maps
source(paste0(indic_repo,  "mbg/functions/map_model_results.r"))
try(map_model_results(indicator, indicator_group, run_date, type = "mean", raked = 'raked', plot_by_year = F))
try(map_model_results(indicator, indicator_group, run_date, type = "cirange", raked = 'raked', plot_by_year = F, include_diff = F))
try(map_model_results(indicator, indicator_group, run_date, type = 'mean', raked = 'raked_c', plot_by_year = F, lvl_limits = c(0,1000)))

## Raked vs unraked comparison maps
source(paste0(indic_repo,"mbg/functions/model_results_compare.r"))
try(model_results_compare(indicator, indicator_group, c(run_date, run_date), c("Unraked", "Raked"), 'mean',
                          c(F, T), floor(seq(min(year_list), max(year_list), length.out = 4)), 'upper',
                          paste0(outputdir, "/diagnostic_plots/", indicator, "_raked_vs_unraked.pdf")))

## ANC correction plots
if (anc_correction != 'none' & as.logical(use_gp)) {
  source(paste0(indic_repo, 'mbg/', indicator, '/3_post_estimation/anc_correction_plots.r'))
  try(anc_correction_plots(indicator, indicator_group, run_date))
}

## ANC spatial field
if (anc_correction == 'gp') {
  source(paste0(indic_repo, 'mbg/', indicator, '/3_post_estimation/plot_anc_spatial_field.r'))
  try(plot_anc_spatial_field(indicator, indicator_group, run_date))
}

###############################################################################
## Close out
###############################################################################

run_date_file <- paste0("<<<< FILEPATH REDACTED >>>>")
if (file.exists(run_date_file)) {
  temp <- fread(run_date_file)
  temp[run_date == get('run_date', .GlobalEnv), done := 1]
  write.table(temp, file = run_date_file, sep = ",", row.names = F, na = '')
}

message('Done!')
