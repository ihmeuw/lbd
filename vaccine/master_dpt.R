###############################################################################
###############################################################################
## MBG Master Launch Script
##
## Purpose: This is the master launch script.  Will run individual models for 
## each vaccine or conditional vaccine data set for ordinal regression using a
## continuation-ratio model, then combine these arithmatically in order to
## the desired vaccine coverage metrics. 
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific code begins here ##########################################
indicator_group    <- 'vaccine'
vaccine            <- "dpt"
indicators         <- c('dpt3_cov', 'dpt2_cond', 'dpt1_cond')         # Indicators to model with INLA, etc.
combine_indicators <- c('dpt3_cov', 'dpt2_cov', 'dpt1_cov')           # Final _cov indicators to mark end of combine script
dose_indicators    <- c('dpt2_dose', 'dpt1_dose', 'dpt0_dose')        # Dose indicators (p(exactly that # doses))
                                                                      # Doesn't include last dose (_cov) as redundant
rake_indicators    <- c('dpt3_cov', 'dpt1_3_abs_dropout', 'dpt1_3_rel_dropout', 'dpt1_cov')
postest_indicators <- c("dpt3_cov", "dpt1_cov", 'dpt1_3_abs_dropout', 'dpt1_3_rel_dropout') # Indicators for post-estimation (i.e. final outputs)
all_indicators     <- unique(c(indicators, combine_indicators, 
                               rake_indicators, dose_indicators, 
                               postest_indicators))
slots              <- 4

## Create run date in correct format
run_date <- make_time_stamp(TRUE)
log_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", user, "/logs/", run_date, "/")
dir.create(log_dir, recursive = T, showWarnings = F)

## Read config file and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = "",
                      indicator       = NULL,
                      config_name     = paste0("config_", vaccine), 
                      covs_name       = paste0("covs_", vaccine))

# distribute config to all indicator directories
distribute_config(cfg = config, indicators = all_indicators)

# parse formatting
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

indicator <- ""

###############################################################################
## Make Holdouts
###############################################################################
if(as.logical(makeholdouts)){

  if(as.logical(load_other_holdouts)){
    message(paste0("Loading holdouts from ", holdout_rundate, "..."))

    stratum_filename <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", vaccine, "3_cov/output/",
                               holdout_rundate, "/stratum.rds") 
    
    stratum_ho <- readRDS(stratum_filename)

    saveRDS(stratum_ho,
            file = paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", vaccine, "3_cov/output/",
                               run_date, "/stratum.rds"))

  } else {

    # load the full input data for biggest vaccine
    df <- load_input_data(indicator   = paste0(vaccine, "3_cov"),
                          simple      = NULL,
                          removeyemen = TRUE,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list)

    # Add in location information
    df <- merge_with_ihme_loc(df)

    # Make admin2 folds in space and time
    shape_dir <- '<<<< FILEPATH REDACTED >>>>' # admin rasters and shapes used for spatial holdouts

    indicator <- paste0(vaccine, "3_cov") # referenced internally in make_folds

    stratum_ho <- make_folds(data       = df,
                             n_folds    = as.numeric(n_ho_folds),
                             spat_strat = 'poly',
                             temp_strat = 'prop',
                             strat_cols = 'region',
                             ts         = as.numeric(ho_ts),
                             mb         = as.numeric(ho_mb),
                             lat_col    = "latitude",
                             long_col   = "longitude",
                             ss_col     = "N",
                             admin_shps   = paste0(shape_dir, 'ad2_raster.grd'),
                             admin_raster = paste0(shape_dir, 'africa_ad2.shp'),
                             mask_shape   = paste0(shape_dir, 'africa_simple.shp'),
                             mask_raster  = paste0(shape_dir, 'ad0_raster'))

  }

  # Recreate & distribute holdouts for all indicators

  invisible(lapply((unique(c(indicators, postest_indicators))), function(ind) {

    sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>/%s/%s',indicator_group, ind)

    ## Create holdouts if not already present
    if(as.logical(makeholdouts) & !(file.exists(paste0(sharedir,'/output/',run_date,"/stratum.rds")))){

      message(paste0("Recreating holdouts for ", ind, "..."))
      # load the full input data
      df <- load_input_data(indicator   = ind,
                            simple      = NULL,
                            removeyemen = TRUE,
                            years       = yearload,
                            withtag     = as.logical(withtag),
                            datatag     = datatag,
                            use_share   = as.logical(use_share),
                            yl          = year_list)

      # add in location information
      df <- merge_with_ihme_loc(df)

      # make a list of dfs for each region, with 5 qt folds identified in each
      stratum_ho <- recreate_holdouts(data = df,
                                      row_id_col = "row_id",
                                      load_from_indic = paste0(vaccine, "3_cov"),
                                      rd = run_date,
                                      ig = indicator_group)

      # Save stratum_ho object for reference
      save_stratum_ho(indic       = ind,
                      ig          = indicator_group,
                      rd          = run_date,
                      stratum_obj = stratum_ho)

    } # if(as.logical(makeholdouts) & !(file.exists(...
  })) # invisible(lapply...
} #if(as.logical(makeholdouts))

###############################################################################
## Launch parallel models
###############################################################################

launch_lv = list(indicator = indicators,
                 use_gn = F)

launch_qsub_output <- parallelize(script = "launch_dpt",
                                  log_location = paste0(log_dir, "launch/"),
                                  expand_vars = launch_lv,
                                  save_objs = c("run_date", "indicator_group", "vaccine", "core_repo"),
                                  prefix = "launch",
                                  slots = 2,
                                  memory = 4,
                                  script_dir = indic_repo,
                                  geo_nodes = F, 
                                  singularity = 'default')

monitor_jobs(launch_qsub_output)

###############################################################################
## Combine to create individual dose objects (the _cov objects)
###############################################################################
message("Master script: done with parallel; starting combination")

if (as.logical(makeholdouts) == T) holdout_vector <- 0:as.numeric(n_ho_folds)
if (as.logical(makeholdouts) == F) holdout_vector <- 0

combine_lv<- list(region = Regions,
                  holdout = holdout_vector,
                  doses = 3)

combine_qsub_output <- parallelize(script = "combine_dpt",
                                   log_location = paste0(log_dir, "combine/"),
                                   expand_vars = combine_lv,
                                   save_objs = c("indicator_group", "run_date", "vaccine"),
                                   prefix = "combine",
                                   slots = 16,
                                   memory = 100,
                                   script_dir = indic_repo,
                                   geo_nodes = TRUE, 
                                   singularity = 'default')

monitor_jobs(combine_qsub_output)

###############################################################################
## Rake indicators: dpt1_cov, dpt3_cov, absolute & relative dropout
###############################################################################

# Skip some steps if validation metrics only
if (!as.logical(validation_metrics_only)) {

  message("Master script: done with combination; starting raking")

  rake_lv <- list(region = Regions,
                  holdout = 0)

  rake_qsub_output <- parallelize(script = "rake_dpt",
                                  log_location = paste0(log_dir, "rake/"),
                                  expand_vars = rake_lv,
                                  save_objs = c("indicator_group", "run_date", "vaccine", "rake_indicators", "gbd_date"),
                                  prefix = "rake",
                                  slots = 17,
                                  memory = 100,
                                  geo_nodes = TRUE,
                                  script_dir = indic_repo,
                                  singularity = 'default')

  monitor_jobs(rake_qsub_output)
  message("Master script: done with raking; starting summarization")

} #end validation_metrics_only check

###############################################################################
## Summarize indicators
###############################################################################

if (!as.logical(validation_metrics_only)) {

  # Summarize raked and unraked indicators
  summarize_lv <- rbind(expand.grid(rake_indicators, Regions, T, stringsAsFactors = F),
                        expand.grid(all_indicators, Regions, F, stringsAsFactors = F)) %>%
                  as.data.table %>%
                  setnames(c("Var1", "Var2", "Var3"), c("indicator", "region", "raked"))  

  summarize_qsub_output <- parallelize(script = "summarize_dpt",
                                       log_location = paste0(log_dir, "summarize/"),
                                       lv_table = summarize_lv,
                                       save_objs = c("indicator_group", "run_date", "vaccine", "summstats", "gbd_date"),
                                       prefix = "summ",
                                       slots = 6,
                                       memory = 60,
                                       geo_nodes = TRUE,
                                       #proj = proj_geo_nodes_vac,
                                       script_dir = indic_repo,
                                       singularity = 'default')

  monitor_jobs(summarize_qsub_output)

} #end validation_metrics_only check

###############################################################################
# Merge regional rasters together 

if (!as.logical(validation_metrics_only)) {
  
  merge_lv <- rbind(expand.grid(all_indicators, FALSE, stringsAsFactors = F),
                  expand.grid(rake_indicators, TRUE, stringsAsFactors = F))

  merge_lv <- merge_lv %>% as.data.table %>% setnames(c("Var1", "Var2"), c("indicator", "raked"))
  merge_lv[(indicator %in% indicators) & (raked == F), modeled := T]
  merge_lv[!((indicator %in% indicators) & (raked == F)), modeled := F]

  merge_qsub_output <- parallelize(script = "merge_regions",
                                   log_location = paste0(log_dir, "merge/"),
                                   lv_table = merge_lv,
                                   save_objs = c("indicator_group", "run_date", "vaccine", "summstats"),
                                   prefix = "merge",
                                   slots = 3,
                                   memory = 6,
                                   script_dir = indic_repo,
                                   geo_nodes = T,
                                   #proj = proj_geo_nodes_vac, 
                                   singularity = 'default')

  monitor_jobs(merge_qsub_output)

} # end validation_metrics_only check

###############################################################################
## Make plots of rasters
###############################################################################
if (!as.logical(validation_metrics_only)) {

  plot_lv <- list(indicator = rake_indicators,
                  raked = c(T,F),
                  summstat = summstats,
                  slots = 10)

  plot_raster_output <- parallelize(script = "plot_rasters",
                                    log_location = paste0(log_dir, "plot/"),
                                    expand_vars = plot_lv,
                                    save_objs = c("indicator_group", "run_date", "vaccine", "year_list"),
                                    prefix = "plot",
                                    slots = 14,
                                    memory = 14,
                                    script_dir = indic_repo, 
                                    geo_nodes = F,
                                    singularity = 'default')
}

###############################################################################
## Launch post-estimation for post-estimation indicators
###############################################################################

post_predict_lv <- list(indicator = postest_indicators)

post_pred_qsub_output <- parallelize(script = "post_predict_dpt_rerake",
                                     log_location = paste0(log_dir, "post_predict/"),
                                     expand_vars = post_predict_lv,
                                     save_objs = c("run_date", "indicator_group", "vaccine", "validation_metrics_only"),
                                     prefix = "post_pred",
                                     slots = 10,
                                     memory = 64,
                                     script_dir = indic_repo,
                                     geo_nodes = T,
                                     singularity = 'default')

monitor_jobs(post_pred_qsub_output)

###############################################################################
## Combine & aggregate non-post-estimation indicators
###############################################################################

if (!as.logical(validation_metrics_only)) {

  ## Aggregation for doses (for ternary plots)
  agg_lv <- list(indicator = dose_indicators,
                 region = Regions,
                 holdout = 0,
                 age = 0,
                 raked = c(FALSE),
                 overwrite = TRUE)

  agg_qsub_output <- parallelize(script = "aggregate_results",
                                 script_dir = paste0(core_repo, "/mbg_central/share_scripts/"), 
                                 log_location = paste0(log_dir, "aggregate/"),
                                 expand_vars = agg_lv,
                                 save_objs = c("indicator_group", "run_date", "pop_measure", "year_list"),
                                 prefix = paste0("agg_", indicator),
                                 slots = 8,
                                 memory = 60,
                                 geo_nodes = TRUE, 
                                 singularity = 'default')

  monitor_jobs(agg_qsub_output)

  for (ind in dose_indicators) {
    combine_aggregation(rd = run_date, 
                        indic = ind, 
                        ig = indicator_group,
                        ages = 0, 
                        regions = Regions,
                        holdouts = 0,
                        raked = c(F))
  }

} # end validation_metrics_only check

###############################################################################
## Save outputs in standard MBG directory for mapping
###############################################################################

if (!as.logical(validation_metrics_only)) {

  # Save tifs

  tifs_to_copy <- data.table(expand.grid(summstats, c(T,F), postest_indicators, stringsAsFactors = F))
  names(tifs_to_copy) <- c("measure", "raked", "indicator") 

  for (i in 1:nrow(tifs_to_copy)) {
    copy_tif_to_map_input_dir(ind = tifs_to_copy[i, "indicator"],
                              ig = indicator_group,
                              measure = tifs_to_copy[i, "measure"],
                              rd = run_date,
                              raked = tifs_to_copy[i, "raked"],
                              yl = year_list)
  }

  # Save ads
  admin_lvs <- data.table(expand.grid(summstats, c(T,F), c(0,1,2), postest_indicators, stringsAsFactors = F))
  admin_lvs <- rbind(admin_lvs,
                     data.table(expand.grid("psup80", 
                                            c(T,F), 
                                            c(0,1,2),
                                            postest_indicators[postest_indicators != paste0(vaccine, "1_3_rel_dropout")], 
                                stringsAsFactors = F)),
                     data.table(expand.grid("diff_2000-2016",
                                            c(T),
                                            c(0,1,2),
                                            postest_indicators,
                                stringsAsFactors = F)))

  names(admin_lvs) <- c("measure", "raked", "ad_level", "indicator") 

  for (i in 1:nrow(admin_lvs)) {
    copy_admins_to_map_input_dir(ind = admin_lvs[i, "indicator"],
                                 ig = indicator_group, 
                                 measure = as.character(admin_lvs[i, "measure"]), 
                                 rd = run_date, 
                                 raked = admin_lvs[i, "raked"], 
                                 yl = year_list, 
                                 ad_level = admin_lvs[i, "ad_level"]) 
  }

} #end validation metrics only check

###############################################################################
## Launch number plugging script & other similar publication scripts
###############################################################################
if (!as.logical(validation_metrics_only)) {

  # Number plugging script
  source(paste0(indic_repo, "number_plugging.R"))

  # Plot dfferences
  source(paste0(indic_repo, "share_scripts/plot_dpt3_admin_changes.R"))

  # Make GBD comparisons
  str_match <- stringr::str_match
  gbd_lv <- data.table(expand.grid(paste0(vaccine, c("3_cov", "1_cov")),
                                   c(T,F), stringsAsFactors=F))
  names(gbd_lv) <- c("indicator", "raked")
  gbd_lv[, gbd_ind := str_match(indicator, "(.*)_.*")[,2]]
  gbd_lv[, title := toupper(gbd_ind)]

  for (i in 1:nrow(gbd_lv)) {
    plot_mbg_vs_gbd(ind = gbd_lv[i, indicator],
                    ig = indicator_group,
                    rd = run_date,
                    gbd_ind = gbd_lv[i, gbd_ind],
                    yl = year_list,
                    vax_title = gbd_lv[i, title],
                    raked = gbd_lv[i, raked])  
  }


  # Grab priors and write to CSV for SI
  message("Grabbing thetas...")
  for (indicator in indicators) {
    message(paste0("  ", indicator))
    theta.res <- matrix(ncol = 4, nrow = length(Regions))
    colnames(theta.res) <- c('Theta1 mean', 'Theta1 var.',  'Theta2 mean', 'Theta2 var.')
    rownames(theta.res) <- Regions
    for(reg in Regions){
      message(paste0("     ", reg))
      age <- holdout <- 0
      pathaddin  <-  paste0('_bin',age,'_',reg,'_',holdout)
      load(paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/',
                  indicator, '/model_image_history/', run_date, pathaddin,
                  '.RData'))
      tmp <- param2.matern.orig(mesh_s)
      theta.res[which(Regions == reg), ] <- c(tmp$theta.prior.mean[1],
                                              1 / tmp$theta.prior.prec[1, 1],
                                              tmp$theta.prior.mean[2],
                                              1 / tmp$theta.prior.prec[2, 2])
    }
    num_plugging_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", 
                               indicator, "/output/", run_date, "/number_plugging/")
    dir.create(num_plugging_dir, showWarnings = F)
    write.csv(theta.res, file = paste0(num_plugging_dir, 'theta_priors_for_SI.csv'))
  }

  message("Plotting summary metrics...")
  if (as.logical(makeholdouts)) oos_arg <- c(T,F)
  if (!as.logical(makeholdouts)) oos_arg <- c(F)
    
  # Multiple indicator summary metric plots
  plot_summary_metrics_multindic(indics = c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout"),
                                 ig = indicator_group,
                                 rd = run_date,
                                 use_oos = oos_arg,
                                 ad_levels = c("country", "ad1", "ad2"))

  # Single indicator summary metric plots
  for (i in postest_indicators) {
    plot_summary_metrics(indic = i, 
                         ig = "vaccine", 
                         rd = run_date,
                         use_oos = oos_arg,
                         by_regions = c(T,F),
                         ad_levels = c("country", "ad1", "ad2"))
  }
  
} #end validation metrics only check

###############################################################################
## Notify and finish up
###############################################################################

message("DPT master script done")

# Note that additional scripts for supplementary analyses for the manuscript 
# can be found in the `share_scripts` folder, to be launched independently 

## END OF FILE
###############################################################################