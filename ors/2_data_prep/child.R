rm(list = ls())
## drive locations
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# TBD: Remve all 'setwd()'
core_repo <- repo <-  '<<<< FILEPATH REDACTED >>>>'
setwd(repo)

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = repo)



shp <- commandArgs()[4]
indicator <- commandArgs()[5]
run_date <- commandArgs()[6]
input_version <- commandArgs()[7]
message(paste(shp,indicator,run_date,input_version))

library('rgdal')
library('raster')
library('dplyr')
library('seegMBG')
library('rgeos')
library('data.table')

# set seed
set.seed(98112)

## Load data
load(paste0('<<<< FILEPATH REDACTED >>>>/countdata_', input_version, '.RData'))
polydat <- countdata[[indicator]]
polydat <- subset(polydat, shapefile == shp)

polydat <- subset(polydat, point == 0)
subset <- polydat[which(polydat$shapefile == shp),]

shape_master <- readRDS(paste0('<<<< FILEPATH REDACTED >>>>',shp,'.rds'))
generated_pts <- list()

## Resample polygon data
for (loc in unique(subset$location_code)) {
  shape <- shape_master[shape_master$GAUL_CODE == loc,]
  subset_loc2 <- subset(subset, location_code == loc)
  
  for (q in 1:nrow(subset_loc2)) {
    
    subset_loc3 <- subset_loc2[q,]
    
    year <- subset_loc3$year
    if (year <= 2000) {
      pop_raster <- raster('<<<< FILEPATH REDACTED >>>>/WorldPop_total_global_stack.tif', band = 1)
    } else {
      if (year > 2000 & year <= 2005) {
        pop_raster <- raster('<<<< FILEPATH REDACTED >>>>/WorldPop_total_global_stack.tif', band = 2)
      } else {
        if (year > 2005 & year <= 2010) {
          pop_raster <- raster('<<<< FILEPATH REDACTED >>>>/WorldPop_total_global_stack.tif', band = 3)
        } else {
          pop_raster <- raster('<<<< FILEPATH REDACTED >>>>/WorldPop_total_global_stack.tif', band = 4)
        }
      } 
    }
    
    raster_crop <- mask(crop(x = pop_raster, y = shape), shape)
    if (length(unique(raster_crop)) < 1) {
      samp_pts <- gCentroid(shape)@coords
      samp_pts <- as.data.frame(samp_pts)
      samp_pts$weight <- 1
      
    } else {
      samp_pts <- getPoints(shape = shape, raster = raster_crop, n = 0.01, perpixel = T, prob = T)  
      samp_pts <- as.data.frame(samp_pts)
    }
    
    names(samp_pts) <- c('long', 'lat','weight')
    samp_pts$shapefile <- shp
    
    subset_loc3 <- left_join(samp_pts, subset_loc3, by = 'shapefile')
    subset_loc3$point <- 0
    
    generated_pts[[length(generated_pts) + 1]] <- subset_loc3
  }
    
}

generated_pts2 <- do.call(rbind, generated_pts)
write.csv(generated_pts2, file = paste0('<<<< FILEPATH REDACTED >>>>', indicator, '/', shp,'.csv'))

