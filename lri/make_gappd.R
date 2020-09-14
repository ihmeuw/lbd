###################################################################################
# compare LRI stage 2 results to GAPPD goals
###################################################################################

# for reference, the GAPPD goals for LRI are listed below:
# mortality below 3/1000 by 2025
# incidence reduced by 75% relative to 2010 by 2025

# calculate the folowing at pixel and a0,1,2 level
# incidence goal by 2017, 2025
# mortality goal by 2010, 2017, 2025

# (1) Setup #######################################################################
rm(list = ls())

#user inputs
run_date <- '2019_10_23_16_13_17'
year_list <- c(2000:2017)
shapefile_version <- '2019_09_10'
modeling_shapefile_version <- shapefile_version

#packages and lbd core functions
core_repo       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#custom functions
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/aroc_proj_functions.R')

#directory
outdir <- '<<<< FILEPATH REDACTED >>>>'

# (2) Define goals ################################################################
g1 = add_goal(target_year = 2025,
              target = 0.75, #75% reduction since 2010
              baseline_year = 2010,
              target_type = "less",
              abs_rel = "relative",
              pred_type = 'admin')
g1$measure ='incidence'

g2 = add_goal(target_year = 2025,
              target = 0.75, #75% reduction since 2010
              baseline_year = 2010,
              target_type = "less",
              abs_rel = "relative",
              pred_type = 'cell')
g2$measure ='incidence'

g3 <- add_goal(target_year = 2025,
               target = 3/1000, #1 per 1000
               baseline_year = 2010,
               target_type = "less",
               abs_rel = "absolute",
               pred_type = 'admin')
g3$measure ='mortality'

g4 <- add_goal(target_year = 2025,
               target = 3/1000, #1 per 1000
               baseline_year = 2010,
               target_type = "less",
               abs_rel = "absolute",
               pred_type = 'cell')
g4$measure ='mortality'

# (3) Make projections to 2025 ############################################################
for (measure in c('incidence', 'mortality')){
  make_proj(ind_gp = 'lri',
            ind = 'has_lri',
            rd = run_date,
            type = c("cell", "admin"),
            proj_years = 2025,
            measure = measure,
            raked = TRUE,
            matrix_pred_name = NULL,
            year_list = year_list,
            uselogit = FALSE,
            shapefile_version = shapefile_version)
}


#(4) Compare to goal #############################################################
compare_to_target(ind_gp = 'lri',
                  ind = 'has_lri',
                  rd = run_date,
                  goal_obj = g1,
                  measure = 'incidence',
                  year_list = year_list,
                  uselogit = FALSE,
                  raked = TRUE,
                  shapefile_version = shapefile_version)

compare_to_target(ind_gp = 'lri',
                  ind = 'has_lri',
                  rd = run_date,
                  goal_obj = g2,
                  measure = 'incidence',
                  year_list = year_list,
                  uselogit = FALSE,
                  raked = TRUE,
                  shapefile_version = shapefile_version)

compare_to_target(ind_gp = 'lri',
                  ind = 'has_lri',
                  rd = run_date,
                  goal_obj = g3,
                  measure = 'mortality',
                  year_list = year_list,
                  uselogit = FALSE,
                  raked = TRUE,
                  shapefile_version = shapefile_version)

compare_to_target(ind_gp = 'lri',
                  ind = 'has_lri',
                  rd = run_date,
                  goal_obj = g4,
                  measure = 'mortality',
                  year_list = year_list,
                  uselogit = FALSE,
                  raked = TRUE,
                  shapefile_version = shapefile_version)


# (5) Make 2010 and 2017 probabilities ############################################################

# Central functions for GAPPD currently do not accomodate non-projections as cell pred input: calculate these manually

# admin level: consider the draws <- what percent of draws met the goals?

for (aa in c(0:2)){

  # mort below 0.003 in 2010, 2017
  load('<<<< FILEPATH REDACTED >>>>')
  admin_draws <- get(paste0('admin_', aa)) %>%
    as.data.frame()

  draw_cols <- grep(names(admin_draws), pattern = 'V', value = TRUE)

  absolute_goal_draws <- ifelse(admin_draws[,draw_cols] <= 0.003, 1, 0) %>%
    cbind(admin_draws[,c(1:2)])

  absolute_goal_draws$absolute_goal_prob <- rowMeans(absolute_goal_draws[,draw_cols])

  #subset to 2010 and 2017 goal years and save out
  for (goal_year in c(2010,2017)){
    absolute_year <- filter(absolute_goal_draws, year == goal_year) %>%
      dplyr::select(-draw_cols)

    write.csv(absolute_year, paste0(outdir, 'has_lri_mortality_', goal_year, '_absolute_less_0.003_adm_', aa, '_target_probs.csv'))
    }

  #inc reduced by 75% in 2017 relative to 2010
  load('<<<< FILEPATH REDACTED >>>>')
  admin_draws <- get(paste0('admin_', aa)) %>%
    melt(id.vars = c('year', paste0('ADM', aa, '_CODE')), measure.vars = draw_cols) %>%
    filter(year %in% c(2010,2017)) %>%
    dcast(paste0('ADM', aa, '_CODE + variable ~ year'), value.var = 'value') %>%
    plyr::rename(c('2010' = 'inc_2010', '2017' = 'inc_2017')) %>%
    as.data.table()

  # a 75% redution means inc_2017/inc_2010 <= 0.25
  relative_goal_draws <- admin_draws[,ratio := inc_2017/inc_2010]
  relative_goal_draws <- relative_goal_draws[,goal_met := ifelse(ratio <= 0.25, 1, 0)] %>%
    group_by_at(paste0('ADM', aa, '_CODE')) %>%
    dplyr::summarize(relative_goal_prob = mean(goal_met))

  relative_goal_draws$year = 2017

  write.csv(relative_goal_draws, paste0(outdir, 'has_lri_incidence_', goal_year, '_vs_2010_relative_less_0.75_adm_', aa, '_target_probs.csv'))
}

# pixel level: evaluate cell preds
regions <- get_output_regions(in_dir = '<<<< FILEPATH REDACTED >>>>')

# mort below 0.003 in 2010, 2017
goal_years <- c(2010, 2017)

for (region in regions){
  cell_pred <- get_cell_pred_for_aroc(ind_gp = 'lri',
                                      ind = 'has_lri',
                                      rd = run_date,
                                      reg = region,
                                      measure = 'mortality',
                                      rk = TRUE,
                                      shapefile_version = shapefile_version)

  for (goal_year in goal_years){

    ## grab goal year preds
    year_idx = which(year_list == goal_year)
    year_draws <- cell_pred[which(cell_pred[, 1] == year_idx), ] ## first col is year
    year_draws <- year_draws[, -(1:2)]

    absolute_goal_draws <- ifelse(year_draws <= 0.003, 1, 0)
    absolute_goal_draws[is.na(absolute_goal_draws)] <- 0
    absolute_goal_prob <- rowMeans(absolute_goal_draws)

    # save
    saveRDS(object = absolute_goal_prob,
            file = paste0(outdir, 'has_lri_mortality_', goal_year, '_absolute_less_0.003_cell_target_probs_', region, '.RDs'))

  }
}

# inc reduction of 75% relative to 2010 by 2017
baseline_year <- 2010
goal_year <- 2017

for (region in regions){
  cell_pred <- get_cell_pred_for_aroc(ind_gp = 'lri',
                                      ind = 'has_lri',
                                      rd = run_date,
                                      reg = region,
                                      measure = 'incidence',
                                      rk = TRUE,
                                      shapefile_version = shapefile_version)

  #grab baseline year preds
  baseline_year_idx = which(year_list == baseline_year)
  baseline_year_draws <- cell_pred[which(cell_pred[, 1] == baseline_year_idx), ] ## first col is year
  baseline_year_draws <- baseline_year_draws[, -(1:2)]

  #grab goal year preds
  goal_year_idx = which(year_list == goal_year)
  goal_year_draws <- cell_pred[which(cell_pred[, 1] == goal_year_idx), ] ## first col is year
  goal_year_draws <- goal_year_draws[, -(1:2)]

  #compare goal to baseline year: want goal/baseline <= 0.25 to correspond to 75% reduction
  relative_proj_draws <- ifelse(goal_year_draws / baseline_year_draws <= 0.25, 1, 0)
  relative_proj_draws[is.na(relative_proj_draws)] <- 0
  relative_goal_prob <- rowMeans(relative_proj_draws)

  #save
  saveRDS(object = relative_goal_prob,
          file = paste0(outdir, 'has_lri_incidence_', goal_year, '_vs_2010_relative_less_0.75_cell_target_probs_', region, '.RDs'))

}

# (6) Write regional rasters for cell-level goals ############################################################

for (region in regions){

  #simple raster
  raster_outputs <- prep_shapes_for_raking(
    reg = region,
    modeling_shapefile_version = shapefile_version,
    raking_shapefile_version = shapefile_version,
    field = 'loc_id'
  )
  simple_raster <- raster_outputs[['simple_raster']]

  for (measure in c('incidence', 'mortality')){
    if (measure == 'incidence') goal_years <- c(2017, 2025)
    if (measure == 'mortality') goal_years <- c(2010, 2017, 2025)

    for (goal_year in goal_years){
      results <- readRDS('<<<< FILEPATH REDACTED >>>>')

      raster <- insertRaster(simple_raster, as.data.frame(results))
      writeRaster(raster,
                  file = '<<<< FILEPATH REDACTED >>>>',
                  overwrite = TRUE)
    }
  }
}

# (7) Merge regional <- global rasters for cell-level goals ############################################################

for (measure in c('incidence', 'mortality')){
  if (measure == 'incidence') goal_years <- c(2017, 2025)
  if (measure == 'mortality') goal_years <- c(2010, 2017, 2025)
  
  for (goal_year in goal_years){
    
    # read in raster files
    r <- list.files('<<<< FILEPATH REDACTED >>>>',
                    pattern = paste0('(has_lri_', measure,
                              '_', goal_year,
                              ifelse(measure == 'incidence', '_vs_2010_relative_less_0.75', '_absolute_less_0.003'), '_cell_target_probs_+)\\w+\\.grd'))
    filepaths <- paste0(outdir, r)
    rasters <- lapply(filepaths, brick)
    
    # merge and save raster file
    global_raster <- do.call(raster::merge, rasters)
    
    writeRaster(global_raster,
                file = paste0(outdir, 'has_lri_', measure,
                              '_', goal_year,
                              ifelse(measure == 'incidence', '_vs_2010_relative_less_0.75', '_absolute_less_0.003'), '_cell_target_probs.tif'),
                overwrite = TRUE)
  }
}

