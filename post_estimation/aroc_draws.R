# ---------------------------------------------------------------------------------------------
# AROC analysis of ORS, RHF, and ORT change over time
#
# (1) Get and save AROC for all indicators
# (2) Plot AROC for all indicators
# ---------------------------------------------------------------------------------------------


# Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# indicate whether running interactively
interactive <- FALSE

# if running interactively, set arguments
if (interactive) {
  warning('interactive is set to TRUE - if you did not mean to run code interactively then kill the job and set interactive to FALSE')
  
  ## set arguments
  indicator_group          <- 'ort'
  indicator                <- 'ors'
  run_date                 <- '2019_11_07_13_50_11'
  modeling_shapefile_version <- '2019_09_10'
  
} else {
  
  # otherwise, grab arguments from qsub
  indicator_group          <- commandArgs()[4]; message(indicator_group)
  indicator                <- commandArgs()[5]; message(indicator)
  run_date                 <- commandArgs()[6]; message(run_date)
  modeling_shapefile_version <- commandArgs()[7]; message(modeling_shapefile_version)
  
} 


# set population version
pop_release <- '2019_08_29'

# load packages
library(data.table)
library('png')
library('viridis')
library(data.table)

# set repo
core_repo <- '<<<< FILEPATH REDACTED >>>>'

# indicate whether to plot maps 
plot_maps <- T

# load MBG packages
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# load custom functions
lapply(paste0(core_repo,
              'custom_functions/',
              list.files(paste0(core_repo, '/custom_functions/'))),
       source)

# read in ADM0, ADM1, and ADM2 codes
ad_codes <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '_admin_2_unraked_summary.csv'),
                  select = c('ADM0_NAME', 'ADM0_CODE', 'ADM1_NAME', 'ADM1_CODE', 'ADM2_NAME', 'ADM2_CODE'))

# create outputdir
plot_date <- Sys.Date()
outputdir <- paste0('<<<< FILEPATH REDACTED >>>>')
dir.create(outputdir)

# print out session info so we have it on record
sessionInfo()


# Get AROC -------------------------------------------------------------------------------

message(paste0('Getting AROC for ', indicator))

# set measures
measures <- c('prevalence')

for (m in measures) {
  
  # set input directory
  in_dir <- paste0('<<<< FILEPATH REDACTED >>>>')
  
  # load collapsed microdata
  clean_collapsed_micro_data <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '.csv'))
  
  # make AROC objects
  make_aroc(ind_gp = indicator_group,
            ind = indicator,
            rd = run_date,
            matrix_pred_name = NULL,
            type = c('admin'),
            measure = m,
            year_list = c(2000:2017),
            uselogit = FALSE,
            raked = FALSE,
            weighting_res = 'domain',
            weighting_type = 'empirical',
            pow = 1,
            input_data = clean_collapsed_micro_data,
            mult_emp_exp = FALSE,
            shapefile_version = modeling_shapefile_version)
  
  # read in and summarize AROC for admin 0
  ad0 <- readRDS(paste0(in_dir, indicator, '_', m, '_aroc_adm0_draw_matrix.RDs'))
  ad0_summary <- data.table(ADM0_CODE = ad0[1:nrow(ad0)],
                            mean = rowMeans(ad0[, 2:ncol(ad0)]),
                            lower = apply(ad0[, 2:ncol(ad0)], 1, lower),
                            upper = apply(ad0[, 2:ncol(ad0)], 1, upper),
                            year = 2017)
  ad0_summary <- unique(merge(ad0_summary, ad_codes[, c('ADM0_CODE', 'ADM0_NAME')], by = 'ADM0_CODE'))
  write.csv(ad0_summary, paste0(in_dir, indicator, '_admin_0_', 'unraked_', m, '_aroc_summary.csv'))
  
  # read in and summarize AROC for admin 1
  ad1 <- readRDS(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '_', m, '_aroc_adm1_draw_matrix.RDs'))
  ad1_summary <- data.table(ADM1_CODE = ad1[1:nrow(ad1)],
                            mean = rowMeans(ad1[, 2:ncol(ad1)]),
                            lower = apply(ad1, 1, lower),
                            upper = apply(ad1, 1, upper),
                            year = 2017)
  ad1_summary <- unique(merge(ad1_summary, ad_codes[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME')], by = 'ADM1_CODE'))
  write.csv(ad1_summary, paste0(in_dir, indicator, '_admin_1_', 'unraked_', m, '_aroc_summary.csv'))
  
  # read in and summarize AROC for admin 2
  ad2 <- readRDS(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '_', m, '_aroc_adm2_draw_matrix.RDs'))
  ad2_summary <- data.table(ADM2_CODE = ad2[1:nrow(ad2)],
                            mean = rowMeans(ad2[, 2:ncol(ad2)]),
                            lower = apply(ad2[, 2:ncol(ad2)], 1, lower),
                            upper = apply(ad2[, 2:ncol(ad2)], 1, upper),
                            year = 2017)
  ad2_summary <- unique(merge(ad2_summary, ad_codes, by = 'ADM2_CODE'))
  
}


# Plot AROC -------------------------------------------------------------------------------------------

if (plot_maps == T) {
  
  message(paste0('Plotting AROC for ', indicator))
  
  # combine files
  aroc_summaries <- list(ad0_summary, ad1_summary, ad2_summary)
  names(aroc_summaries) <- c('admin0', 'admin1', 'admin2')
  
  # clean up
  rm(ad0); rm(ad1); rm(ad2); rm(ad0_summary); rm(ad1_summary); rm(ad2_summary)
  
  # set arguments
  map_levels <- c('admin0', 'admin2')
  raked_map <- ifelse(indicator == 'had_diarrhea', T, F)
  limit_levels <- c(-0.1, 0.1)
  
  # map mean model results
  map_model_results(indicator = 'exploratory_analysis',
                    indicator_group = 'ort',
                    run_date = paste0('AROC_', plot_date),
                    type = 'mean',
                    raked = raked_map,
                    raked_measure = 'prevalence',
                    lvl_years = 2017,
                    lvl_colors = 'magma',
                    high_is_bad = FALSE,
                    lvl_limits = limit_levels,
                    include_diff = FALSE,
                    geo_levels = map_levels,
                    plot_by_year = TRUE,
                    plot_combined = FALSE,
                    file_type = 'png',
                    custom_data = aroc_summaries,
                    file_addin = indicator)
  
  # map uppper model results
  map_model_results(indicator = 'exploratory_analysis',
                    indicator_group = 'ort',
                    run_date = paste0('AROC_', plot_date),
                    type = 'upper',
                    raked = raked_map,
                    raked_measure = 'prevalence',
                    lvl_years = 2017,
                    lvl_colors = 'magma',
                    high_is_bad = FALSE,
                    lvl_limits = limit_levels,
                    include_diff = FALSE,
                    geo_levels = map_levels,
                    plot_by_year = TRUE,
                    plot_combined = FALSE,
                    file_type = 'png',
                    custom_data = aroc_summaries,
                    file_addin = indicator)
  
  # map lower model results
  map_model_results(indicator = 'exploratory_analysis',
                    indicator_group = 'ort',
                    run_date = paste0('AROC_', plot_date),
                    type = 'lower',
                    raked = raked_map,
                    raked_measure = 'prevalence',
                    lvl_years = 2017,
                    lvl_colors = 'magma',
                    high_is_bad = FALSE,
                    lvl_limits = limit_levels,
                    include_diff = FALSE,
                    geo_levels = map_levels,
                    plot_by_year = TRUE,
                    plot_combined = FALSE,
                    file_type = 'png',
                    custom_data = aroc_summaries,
                    file_addin = indicator)
}

