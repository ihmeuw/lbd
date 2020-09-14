#####################################################################
# Post estimation script for LRI
# (1) combine aggregation
# (2) aggregate data and stackers
# (3) make line plots and covariate weight plots
# (4) plot priors for spatial hyperparameters
# (5) merge global rasters
#####################################################################

# Setup -------------------------------------------------------------

# agruments from qsub
user            <- commandArgs()[4]
core_repo       <- commandArgs()[5]
indicator_group <- commandArgs()[6]
indicator       <- commandArgs()[7]
config_par      <- commandArgs()[8]
cov_par         <- commandArgs()[9]
run_date        <- commandArgs()[10]
measures        <- unlist(strsplit(commandArgs()[11], '~'))
force_plot      <- as.logical(commandArgs()[12])
year_tag        <- commandArgs()[13]
two_xgboost     <- as.logical(commandArgs()[14])

#set force plot
if (force_plot == T){
  message('Force plot is set to TRUE: replacing all NAs in mbg data frame with 0 and NaN coefficients with 1.')
}

# Load MBG packages
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Load custom lbd_core functions
lapply(paste0('<<<< FILEPATH REDACTED >>>>',
              list.files('<<<< FILEPATH REDACTED >>>>', pattern = 'functions')),
       source)

# Load custom post-estimation functions
lapply(paste0('<<<< FILEPATH REDACTED >>>>',
              list.files('<<<< FILEPATH REDACTED >>>>')),
       source)

# Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')

# Read config file and save all parameters in memory
config <- set_up_config(repo = '<<<< FILEPATH REDACTED >>>>',
                        core_repo = core_repo,
                        indicator_group,
                        indicator,
                        config_name = paste0('config_', config_par),
                        covs_name = paste0('covs_', cov_par),
                        run_tests = FALSE)

# Create output folder with the run_date
outputdir      <- '<<<< FILEPATH REDACTED >>>>'

# Create proper year list object
if (class(year_list) == 'character') year_list <- eval(parse(text=year_list))

# Get regions
Regions <- get_output_regions(outputdir)

if (NA %in% Regions){
  Regions <- Regions[!is.na(Regions)]
}

# Set holdout to 0: code currently set up for diagnostics on full model only
holdouts <- 0


# (1) Combine and summarize aggregated results --------------------------

# combine unraked results (rates)
if ('prevalence' %in% measures){
  message('Combining unraked aggregated results')
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = Regions,
                      holdouts = holdouts,
                      raked    = F,
                      delete_region_files = F,
                      measure = 'prevalence',
                      metrics = c('rates', 'counts'))

# summarize admins
  summarize_admins(ad_levels = c(0,1,2), raked = F, measure = 'prevalence', metrics = c('rates','counts'))
}

for (measure in measures){
  # combine raked results (rates)
  message(paste0('Combining raked aggregated results for ', measure))
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = Regions,
                      holdouts = holdouts,
                      raked    = T,
                      measure  = measure,
                      delete_region_files = F,
                      metrics = c('rates','counts'))

  # summarize admins
  summarize_admins(ad_levels = c(0,1,2), raked = T, measure = measure, metrics = c('rates','counts'))
}

# (2) Aggregate data and stackers ------------------------------------------------------

# Aggregate data to admin 0 and 1
dat <- aggregate_input_data(indicator,
                            indicator_group,
                            run_date,
                            Regions,
                            modeling_shapefile_version)

# Aggregate stackers to admin 0 and 1
stack <- aggregate_child_stackers(indicator,
                                  indicator_group,
                                  run_date,
                                  Regions,
                                  modeling_shapefile_version,
                                  pop_release = pop_release)


# Combine data and stackers with summary results ------------------------------------------

# Load unraked estimates
mbg <- list(fread('<<<< FILEPATH REDACTED >>>>'),
            fread('<<<< FILEPATH REDACTED >>>>'))
names(mbg) <- c('ad0', 'ad1')

# Load raked estimates
mbg_raked <- list(fread('<<<< FILEPATH REDACTED >>>>'),
                  fread('<<<< FILEPATH REDACTED >>>>'))
names(mbg_raked) <- c('ad0', 'ad1')


# Combine
for (a in c('ad0', 'ad1')) {

  # raked results
  rake_cols <- c('mean', 'upper', 'lower', 'cirange')
  setnames(mbg_raked[[a]], rake_cols, paste0(rake_cols, '_raked'))
  mbg[[a]] <- merge(mbg[[a]], mbg_raked[[a]],
                    by = names(mbg[[a]])[grep('ADM|region|year|pop', names(mbg[[a]]))],
                    all.x = T)

  # stackers
  mbg[[a]] <- merge(mbg[[a]], stack[[a]],
                    by = names(mbg[[a]])[grep('CODE|year', names(mbg[[a]]))],
                    all.x = T)

  # data
  mbg[[a]] <- merge(mbg[[a]], dat[[a]],
                    by = names(mbg[[a]])[grep('CODE|year', names(mbg[[a]]))],
                    all.x = T)

  # save
  write.csv(mbg[[a]], '<<<< FILEPATH REDACTED >>>>')
}

if (force_plot == TRUE){
  mbg$ad1[is.na(mbg$ad1)] <- 0
  mbg$ad0[is.na(mbg$ad0)] <- 0
}

# (3) Plot stackers and covariates ------------------------------------------------------

dir.create(paste0(outputdir, '/diagnostic_plots/'))

if (use_stacking_covs) {

  # plot covariate weights
  message('Making covariate weight plots')
  get_cov_weights(indicator,
                  indicator_group,
                  run_date,
                  Regions,
                  outputdir)

  source('<<<< FILEPATH REDACTED >>>>/stacker_admin_time_series.R')

  # plot stackers over time aggregated to admins
  message('Making time series plots for stackers by admin unit')
  plot_stackers_by_adm01(admin_data = mbg,
                         shapefile_version = raking_shapefile_version,
                         indicator,
                         indicator_group,
                         run_date,
                         Regions,
                         measure = measure,
                         raked = T,
                         credible_interval = 0.95,
                         force_plot = force_plot,
                         two_xgboost = two_xgboost)

}


# (4) Plot priors for spatial hyperparameters --------------------------------------------------

message('Plotting spatial hyperparameters prior and posteriors')
for (reg in Regions){
  plot_hyperparameters(indicator = indicator,
                       indicator_group = indicator_group,
                       run_date = run_date,
                       age = 0,
                       holdout = holdouts,
                       save_file = NULL,
                       reg = reg)
}

# (5) Merge global rasters ---------------------------------------------------------------------

# make folder to save out to
dir.create(paste0(outputdir, 'pred_derivatives/global_rasters/'))

# loop over measures
for (measure in measures) {

  # loop over stats
  for (s in c('prediction', 'upper', 'lower')) {

    # read in raster files
    r <- list.files(paste0(outputdir),
                    pattern = glob2rx(paste0('*', s, '*', measure, '*.grd*')))
    filepaths <- paste0(outputdir, r)
    rasters <- lapply(filepaths, brick)

    # merge and save raster file
    global_raster <- do.call(raster::merge, rasters)
    writeRaster(global_raster,
                paste0(outputdir, 'pred_derivatives/global_rasters/', indicator, '_', measure, '_', ifelse(s == 'prediction', 'mean', s), '_raked_', year_tag, '.tif'),
                overwrite = TRUE)
  }
}
