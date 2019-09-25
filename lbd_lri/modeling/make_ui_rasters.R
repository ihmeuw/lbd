#makes mean, lower, upper rasters from cell draws

#load in required functions
source('<<<< FILEPATH REDACTED >>>>') #raking functions
source('<<<< FILEPATH REDACTED >>>>') #prep functions
commondir      <- '<<<< FILEPATH REDACTED >>>>'
core_repo = '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>', header = FALSE)))
source('<<<< FILEPATH REDACTED >>>>') #setup
mbg_setup(package_list = package_list, repos = core_repo)

#set parameters
region_list <- c('wssa','essa','cssa','sssa','name')
shapefile_version <- 'current'
year_list <- c(2000:2017)
indicator <- 'lriincidence'
measure <- 'incidence'
indicator_group <- 'lri'
run_date <- '<<<< RUN DATE REDACTED >>>>'


for (reg in region_list){
  #load raked cell_pred
  message(paste0("starting summaries for ", reg))
  
  cell_pred <- readRDS('<<<< FILEPATH REDACTED >>>>')
  simple_raster <- cell_pred$new_simple_raster
  cell_pred <- cell_pred$raked_cell_pred
  
  message(paste0("making mean raster for ", reg))
  
  #make mean raster
  mean_raster <- make_cell_pred_summary(draw_level_cell_pred = cell_pred,
                                        mask                 = simple_raster,
                                        return_as_raster     = TRUE,
                                        summary_stat         = 'mean')
  gc()
  message(paste0("making lower raster for ", reg))
  
  #make lower raster
  lower_raster <- make_cell_pred_summary(draw_level_cell_pred = cell_pred,
                                         mask                 = simple_raster,
                                         return_as_raster     = TRUE,
                                         summary_stat         = 'lower')
  gc()
  message(paste0("making upper raster for ", reg))
  
  #make upper raster
  upper_raster <- make_cell_pred_summary(draw_level_cell_pred = cell_pred,
                                         mask                 = simple_raster,
                                         return_as_raster     = TRUE,
                                         summary_stat         = 'upper')
  gc()
  
  assign(sprintf('%s_lower_raked_raster',reg),lower_raster)
  assign(sprintf('%s_upper_raked_raster',reg),upper_raster)
  assign(sprintf('%s_mean_raked_raster',reg),mean_raster)
  rm(lower_raster)
  rm(upper_raster)
  rm(mean_raster)
  rm(cell_pred)
  rm(simple_raster)
  gc()
}

#merge regional rasters
m_mean = do.call(raster::merge,list(name_mean_raked_raster,
                                    essa_mean_raked_raster,
                                    sssa_mean_raked_raster,
                                    wssa_mean_raked_raster,
                                    cssa_mean_raked_raster))
gc()

m_lower = do.call(raster::merge,list(name_lower_raked_raster,
                                     essa_lower_raked_raster,
                                     sssa_lower_raked_raster,
                                     wssa_lower_raked_raster,
                                     cssa_lower_raked_raster))
gc()
m_upper = do.call(raster::merge,list(name_upper_raked_raster,
                                     essa_upper_raked_raster,
                                     sssa_upper_raked_raster,
                                     wssa_upper_raked_raster,
                                     cssa_upper_raked_raster))
gc()
#save merged lower and upper rasters
save_post_est(m_lower,'raster','lower_raked_raster')
save_post_est(m_upper,'raster','upper_raked_raster')
save_post_est(m_mean, 'raster','mean_raked_raster')