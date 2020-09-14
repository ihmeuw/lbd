# ---------------------------------------------------------------------------------------------
# AROC of LRI change over time
#
# (1) Get and save AROC for all indicators
# (2) Generate cell-level AROC rasters for each region
# (3) Merge region-level rasters into global raster
#
# ---------------------------------------------------------------------------------------------


# (1) Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# indicate whether running interactively
interactive <- FALSE

# if running interactively, set arguments
if (interactive) {
  warning('interactive is set to TRUE - if you did not mean to run code interactively then kill the job and set interactive to FALSE')
  
  ## set arguments
  indicator_group          <- 'lri'
  indicator                <- 'has_lri'
  run_date                 <- '2019_10_23_16_13_17'
  modeling_shapefile_version <- '2019_09_10'
  
} else {
  
  # otherwise, grab arguments from qsub
  indicator_group          <- commandArgs()[4]; message(indicator_group)
  indicator                <- commandArgs()[5]; message(indicator)
  run_date                 <- commandArgs()[6]; message(run_date)
  modeling_shapefile_version <- commandArgs()[7]; message(modeling_shapefile_version)
  mortality_2010_end_year <- as.logical(commandArgs()[8])
  end_year <- as.numeric(commandArgs()[9])
  
} 

# load packages
library(data.table)
library('png')
library('viridis')
library(data.table)

# set repo
core_repo <- '<<<< FILEPATH REDACTED >>>>'

# load MBG packages
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#load custom functions
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/aroc_proj_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/prep_functions.R')

# read in ADM0, ADM1, and ADM2 codes
ad_codes <- fread('<<<< FILEPATH REDACTED >>>>',
                  select = c('ADM0_NAME', 'ADM0_CODE', 'ADM1_NAME', 'ADM1_CODE', 'ADM2_NAME', 'ADM2_CODE'))

# create outputdir
outputdir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(outputdir)

# print out session info so we have it on record
sessionInfo()


# (2) Get AROC -------------------------------------------------------------------------------

message(paste0('Getting AROC for ', indicator))
 
# set measures
measures <- c('prevalence', 'incidence', 'mortality')

for (m in measures) {

  # set input directory
  in_dir <- '<<<< FILEPATH REDACTED >>>>'

  # load collapsed microdata
  clean_collapsed_micro_data <- fread('<<<< FILEPATH REDACTED >>>>')

  # make AROC objects
  make_aroc(ind_gp = indicator_group,
            ind = indicator,
            rd = run_date,
            matrix_pred_name = NULL,
            type = c('cell','admin'),
            measure = m,
            year_list = c(2000:end_year),
            uselogit = FALSE,
            raked = TRUE,
            weighting_res = 'domain',
            weighting_type = 'exponential',
            pow = 1,
            mult_emp_exp = FALSE,
            shapefile_version = modeling_shapefile_version)

  # read in and summarize AROC for admin 0
  ad0 <- readRDS(paste0(in_dir, indicator, '_', m, '_aroc_adm0_draw_matrix.RDs'))
  ad0_summary <- data.table(ADM0_CODE = ad0[1:nrow(ad0)],
                            mean = rowMeans(ad0[, 2:ncol(ad0)]),
                            lower = apply(ad0[, 2:ncol(ad0)], 1, lower),
                            upper = apply(ad0[, 2:ncol(ad0)], 1, upper),
                            year = end_year)
  ad0_summary <- unique(merge(ad0_summary, ad_codes[, c('ADM0_CODE', 'ADM0_NAME')], by = 'ADM0_CODE'))
  write.csv(ad0_summary, paste0(in_dir, indicator, '_admin_0_raked_', m, '_aroc_summary.csv'))

  # read in and summarize AROC for admin 1
  ad1 <- readRDS(paste0(in_dir, indicator, '_', m, '_aroc_adm1_draw_matrix.RDs'))
  ad1_summary <- data.table(ADM1_CODE = ad1[1:nrow(ad1)],
                            mean = rowMeans(ad1[, 2:ncol(ad1)]),
                            lower = apply(ad1, 1, lower),
                            upper = apply(ad1, 1, upper),
                            year = end_year)
  ad1_summary <- unique(merge(ad1_summary, ad_codes[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME')], by = 'ADM1_CODE'))
  write.csv(ad1_summary, paste0(in_dir, indicator, '_admin_1_raked_', m, '_aroc_summary.csv'))

  # read in and summarize AROC for admin 2
  ad2 <- readRDS(paste0(in_dir, indicator, '_', m, '_aroc_adm2_draw_matrix.RDs'))
  ad2_summary <- data.table(ADM2_CODE = ad2[1:nrow(ad2)],
                            mean = rowMeans(ad2[, 2:ncol(ad2)]),
                            lower = apply(ad2[, 2:ncol(ad2)], 1, lower),
                            upper = apply(ad2[, 2:ncol(ad2)], 1, upper),
                            year = end_year)
  ad2_summary <- unique(merge(ad2_summary, ad_codes, by = 'ADM2_CODE'))
  write.csv(ad2_summary, paste0(in_dir, indicator, '_admin_2_raked_', m, '_aroc_summary.csv'))

}

#create 2010-end_year mortality aroc for LRI stg2 figure
if (mortality_2010_end_year == TRUE){
  
  # make AROC objects
  make_aroc(ind_gp = indicator_group,
            ind = indicator,
            rd = run_date,
            matrix_pred_name = NULL,
            type = c('admin'),
            measure = 'mortality',
            year_list = c(2010:end_year),
            uselogit = FALSE,
            raked = TRUE,
            weighting_res = 'domain',
            weighting_type = 'exponential',
            pow = 1,
            mult_emp_exp = FALSE,
            shapefile_version = modeling_shapefile_version,
            extra_file_tag = paste0('_2010_', end_year),
            alt_baseline_year = TRUE,
            baseline_year = '2010')
  
  # read in and summarize AROC for admin 0
  ad0 <- readRDS(paste0(in_dir, indicator, '_', m, '_aroc_adm0_draw_matrix_2010_', end_year, '.RDs'))
  ad0_summary <- data.table(ADM0_CODE = ad0[1:nrow(ad0)],
                            mean = rowMeans(ad0[, 2:ncol(ad0)]),
                            lower = apply(ad0[, 2:ncol(ad0)], 1, lower),
                            upper = apply(ad0[, 2:ncol(ad0)], 1, upper),
                            year = paste0('2010_', end_year))
  ad0_summary <- unique(merge(ad0_summary, ad_codes[, c('ADM0_CODE', 'ADM0_NAME')], by = 'ADM0_CODE'))
  write.csv(ad0_summary, paste0(in_dir, indicator, '_admin_0_raked_', m, '_aroc_summary_2010_', end_year, '.csv'))
  
  # read in and summarize AROC for admin 1
  ad1 <- readRDS(paste0(in_dir, indicator, '_', m, '_aroc_adm1_draw_matrix_2010_', end_year, '.RDs'))
  ad1_summary <- data.table(ADM1_CODE = ad1[1:nrow(ad1)],
                            mean = rowMeans(ad1[, 2:ncol(ad1)]),
                            lower = apply(ad1, 1, lower),
                            upper = apply(ad1, 1, upper),
                            year = paste0('2010_', end_year))
  ad1_summary <- unique(merge(ad1_summary, ad_codes[, c('ADM0_CODE', 'ADM0_NAME', 'ADM1_CODE', 'ADM1_NAME')], by = 'ADM1_CODE'))
  write.csv(ad1_summary, paste0(in_dir, indicator, '_admin_1_raked_', m, '_aroc_summary_2010_', end_year, '.csv'))
  
  # read in and summarize AROC for admin 2
  ad2 <- readRDS(paste0(in_dir, indicator, '_', m, '_aroc_adm2_draw_matrix_2010_', end_year, '.RDs'))
  ad2_summary <- data.table(ADM2_CODE = ad2[1:nrow(ad2)],
                            mean = rowMeans(ad2[, 2:ncol(ad2)]),
                            lower = apply(ad2[, 2:ncol(ad2)], 1, lower),
                            upper = apply(ad2[, 2:ncol(ad2)], 1, upper),
                            year = paste0('2010_', end_year))
  ad2_summary <- unique(merge(ad2_summary, ad_codes, by = 'ADM2_CODE'))
  write.csv(ad2_summary, paste0(in_dir, indicator, '_admin_2_raked_', m, '_aroc_summary_2010_', end_year, '.csv'))
}

# (3) Create regional-level rasters -------------------------------------------------------------------------

regions <- get_output_regions(in_dir = '<<<< FILEPATH REDACTED >>>>')

for (region in regions){

  #simple raster
  raster_outputs <- prep_shapes_for_raking(
    reg = region,
    modeling_shapefile_version = modeling_shapefile_version,
    raking_shapefile_version = modeling_shapefile_version,
    field = 'loc_id'
  )
  simple_raster <- raster_outputs[['simple_raster']]

  for (measure in measures){
    draw_results <- readRDS(paste0(outputdir, 'has_lri_', measure,
                              '_aroc_cell_draw_matrix_', region, '.RDs'))
    results <- rowMeans(draw_results)
    raster <- insertRaster(simple_raster, as.data.frame(results))

    writeRaster(raster,
                file = paste0(outputdir, 'has_lri_', measure,
                                     '_aroc_summary_', region),
                overwrite = TRUE)
  }
}


# (7) Merge regional <- global rasters for cell-level goals ############################################################

for (measure in measures){
  # read in raster files
  r <- list.files(outputdir,
                  pattern = paste0('has_lri_', measure,
                                   '_aroc_summary_+\\w+\\.grd'))
  filepaths <- paste0(outputdir, r)
  rasters <- lapply(filepaths, brick)

  # merge and save raster file
  global_raster <- do.call(raster::merge, rasters)

  writeRaster(global_raster,
              file = paste0(outputdir, 'has_lri_', measure,
                            '_aroc_summary.tif'),
              overwrite = TRUE)
}