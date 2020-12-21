##############################################################################
## MBG diagnostics functions and plots for ORT
##############################################################################


## Setup -------------------------------------------------------------------------

## clear environment
rm(list=ls())

## Set repo location, indicator group, and some arguments
user            <- commandArgs()[4]
core_repo       <- commandArgs()[5]
indicator_group <- commandArgs()[6]
indicator       <- commandArgs()[7]
config_par      <- commandArgs()[8]
config_file     <- commandArgs()[9]
cov_par         <- commandArgs()[10]
cov_file        <- commandArgs()[11]
message(indicator)

## Load MBG packages
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## Load custom post-estimation functions
lapply(paste0(core_repo,
              'custom_functions/',
              list.files(paste0(core_repo, 'custom_functions/'))),
       source)

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')

## Read config file and save all parameters in memory
config <- set_up_config(repo            = core_repo,
                        indicator_group = '',
                        indicator       = '',
                        config_name     = paste0(config_file, config_par),
                        covs_name       = paste0(cov_file, cov_par))
                      
## Set run date(s)
run_date <- commandArgs()[12]

## Set measure
measure <- commandArgs()[13]

## Create output folder with the run_date
outputdir      <- paste('<<<< FILEPATH REDACTED >>>>')

## Create proper year list object
if (class(year_list) == 'character') year_list <- eval(parse(text=year_list))

## Get regions
Regions <- get_output_regions(outputdir)

## Set holdout to 0 because for now we'll just run the cleaning and stacker line plots on the full model
holdouts <- 0


## Combine and summarize aggregated results --------------------------

# combine unraked results
message('Combining unraked aggregated results')
if (measure == 'prevalence') {
  combine_aggregation(rd       = run_date,
                      indic    = indicator,
                      ig       = indicator_group,
                      ages     = 0,
                      regions  = Regions,
                      holdouts = holdouts,
                      raked    = F,
                      delete_region_files = F,
                      metrics = 'rates')
  
  # summarize admins
  summarize_admins(ad_levels = c(0,1,2), raked = F, measure = measure, metrics = 'rates')
}


## Aggregate data and stackers ------------------------------------------------------

if (measure == 'prevalence') {

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
                                    modeling_shapefile_version)

}


## Combine data and stackers with summary results ------------------------------------------


if (measure == 'prevalence') {

  # Load unraked estimates
  mbg <- list(fread(paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, '_admin_0_unraked_summary.csv')),
              fread(paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, '_admin_1_unraked_summary.csv')))
  names(mbg) <- c('ad0', 'ad1')

  # Combine
  for (a in c('ad0', 'ad1')) {

    # stackers
    mbg[[a]] <- merge(mbg[[a]], stack[[a]],
                      by = names(mbg[[a]])[grep('CODE|year', names(mbg[[a]]))],
                      all.x = T)

    # data
    mbg[[a]] <- merge(mbg[[a]], dat[[a]],
                      by = names(mbg[[a]])[grep('CODE|year', names(mbg[[a]]))],
                      all.x = T)

    # save
    write.csv(mbg[[a]], paste0(outputdir, '/pred_derivatives/admin_summaries/', indicator, '_mbg_data_stackers_', a, '.csv' ))
  }
}


## Plot stackers and covariates ------------------------------------------------------

dir.create(paste0(outputdir, '/diagnostic_plots/'))

if (use_stacking_covs) {

  if (measure == 'prevalence') {

    # plot sackers over time aggregated to admins
    message('Making time series plots for stackers by admin unit')
    plot_stackers_by_adm01(admin_data = mbg,
                           indicator,
                           indicator_group,
                           run_date,
                           Regions,
                           measure = measure,
                           raked = FALSE,
                           credible_interval = 0.95,
                           N_breaks = c(0, 10, 50, 100, 500, 1000, 2000, 4000))

    # plot covariate weights
    message('Making covariate weight plots')
    get_cov_weights(indicator,
                    indicator_group,
                    run_date,
                    Regions,
                    outputdir)

  }

}


# Plot priors for spatial hyperparameters --------------------------------------------------

message('Plotting spatial hyperparameters prior and posteriors')

if (measure == 'prevalence') {

  # Plot a la HIV team
  plot_hyperparameters(indicator = indicator,
                       indicator_group = indicator_group,
                       run_date = run_date,
                       age = 0,
                       holdout = holdouts,
                       save_file = NULL,
                       regs = Regions)
}


# Make model fit statistics --------------------------------------------------

# fit statistics for prevalence
if (measure == 'prevalence') {

  # Combine input data for IS plotting
  csvs <- list.files(outputdir, pattern = 'input_data_(.*).csv', full.names = T)
  csv_master <- rbindlist(lapply(csvs, fread))
  csv_master[, grep('V1', names(csv_master)) := NULL]
  csv_master <- unique(csv_master)
  write.csv(csv_master, file=paste0(outputdir, '/input_data.csv'))

  # Get in and out of sample draws
  run_in_oos <- get_is_oos_draws(ind_gp        = indicator_group,
                                 ind           = indicator,
                                 rd            = run_date,
                                 ind_fm        = 'binomial',
                                 age           = 0,
                                 nperiod       = length(year_list),
                                 yrs           = year_list,
                                 get.oos       = as.logical(makeholdouts),
                                 year_col      = 'year',
                                 write.to.file = TRUE,
                                 shapefile_version = modeling_shapefile_version)

  # Set out_dir
  dir.create(paste0(outputdir, '/summary_metrics/'), recursive = T, showWarnings = F)

  # Calculate and save PV summary statistics
  draws.df <- fread(paste0(outputdir, '/output_draws_data.csv'))
  pvtab <- rbindlist(lapply(list(c('oos'), c('oos', 'year'), c('oos', 'region')), function(results_by) {
    pv <- get_pv_table(d               = draws.df,
                       indicator_group = indicator_group,
                       rd              = run_date,
                       indicator       = indicator,
                       aggregate_on    = c('country', 'ad1', 'ad2'),
                       result_agg_over = results_by,
                       coverage_probs  = c(25, 50, 75, 90, 95),
                       plot_ci         = TRUE,
                       draws           = as.numeric(samples),
                       save_csv        = FALSE,
                       out.dir         = paste0(outputdir, '/summary_metrics/'))
    ldply(pv, .id = 'Level')
  }), fill = T)
  write.csv(pvtab, file = paste0(outputdir, '/summary_metrics/pv_metrics.csv'), row.names=F)

  # plot results
  source('<<<< FILEPATH REDACTED >>>>/plot_pv_results.R')
  plot_pv_results(indicator,
                  indicator_group,
                  run_dates = run_date,
                  save_file = paste0(outputdir, '/summary_metrics/pv_metrics.pdf'))
}