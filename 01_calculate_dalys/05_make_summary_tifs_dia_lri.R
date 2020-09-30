#####################################################################
# Combine region-level tifs into single results tifs
#####################################################################
#(1) Setup ---------------------------------------------------------
rm(list = ls())

run_date <- '2019_10_28'
measures <- c('daly','yll','yld')
indicators <- c('had_malaria')
mal_run_date <- '2019_10_28'

# set directories
tridaly_in_dir <- '<<< FILEPATH REDACTED >>>'
out_dir <- paste0(tridaly_in_dir, 'results/')
mal_in_dir <- '<<< FILEPATH REDACTED >>>'
dir.create(out_dir)

require(raster)

for (indicator in indicators){
  if (indicator == 'has_lri' | indicator == 'had_diarrhea') in_dir <- tridaly_in_dir
  if (indicator == 'had_malaria') in_dir <- mal_in_dir
  
  for (measure in  measures){
    print(paste('working on', indicator, measure))
    
    if (indicator == 'had_malaria' & measure %in% c('yll','yld')) mean_name <- 'mean' else mean_name <- 'prediction'
    
    # read in mean estimate raster files and merge
    g <- list.files(in_dir, pattern = glob2rx(paste0(indicator, '_', mean_name, '_', measure, '*.grd*')))
    filepaths <- paste0(in_dir, g)
    rasters <- lapply(filepaths, brick)
    if (indicator != 'had_malaria') global_raster <- do.call(raster::merge, rasters)
    if (indicator == 'had_malaria') global_raster <- rasters[[1]]
    
    # read in upper estimate raster files and merge
    g <- list.files(in_dir, pattern = glob2rx(paste0(indicator, '_upper_', measure, '*.grd*')))
    filepaths <- paste0(in_dir, g)
    rasters <- lapply(filepaths, brick)
    if (indicator != 'had_malaria') upper_global_raster <- do.call(raster::merge, rasters)
    if (indicator == 'had_malaria') upper_global_raster <- rasters[[1]]
    
    # read in lower estimate raster files and merge
    g <- list.files(in_dir, pattern = glob2rx(paste0(indicator, '_lower_', measure, '*.grd*')))
    filepaths <- paste0(in_dir, g)
    rasters <- lapply(filepaths, brick)
    if (indicator != 'had_malaria') lower_global_raster <- do.call(raster::merge, rasters)
    if (indicator == 'had_malaria') lower_global_raster <- rasters[[1]]
    
    # save mean, upper, and lower rasters for viz team
    writeRaster(global_raster, paste0(out_dir, indicator, '_', measure, '_mean_raked_2000_2017.tif'), overwrite = TRUE)
    writeRaster(upper_global_raster, paste0(out_dir, indicator, '_', measure, '_upper_raked_2000_2017.tif'), overwrite = TRUE)
    writeRaster(lower_global_raster, paste0(out_dir, indicator, '_', measure, '_lower_raked_2000_2017.tif'), overwrite = TRUE)
  }
}  
