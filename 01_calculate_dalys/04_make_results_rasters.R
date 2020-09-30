#####################################################################
# Combine region-level tifs into single results tifs
#####################################################################

#(1) Setup ---------------------------------------------------------
rm(list = ls())

#set user arguments
causes <- c('malaria')
mal_run_date <- '2019_10_28'

tridaly_run_date <- '2019_10_28'

tridaly_dir <- '<<< FILEPATH REDACTED >>>'

# loop over stats
for (s in c('prediction', 'upper', 'lower')) {
  
  # read in mean raster files
  r <- list.files(paste0(tridaly_dir), 
                  pattern = glob2rx(paste0('*', s, '*daly*.grd*')))
  filepaths <- paste0(share_dir, r)
  rasters <- lapply(filepaths, brick)
  
  if (s == 'prediction') s <- 'mean'
  
  # merge and save mean raster file
  global_raster <- do.call(raster::merge, rasters)
  writeRaster(global_raster, 
              paste0(save_dir, indicator, '_', s, '_daly_raster.tif'),
              overwrite = TRUE)
}