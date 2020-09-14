rm(list = ls())
# Set library and load packages
## drive locations
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

setwd('<<<< FILEPATH REDACTED >>>>')
package_lib <- '<<<< FILEPATH REDACTED >>>>'
.libPaths(package_lib)

for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}



shp <- commandArgs()[6]
indic <- commandArgs()[7]
run_date <- commandArgs()[8]
latest_collapse <- commandArgs()[9]

polydat <- fread('<<<< FILEPATH REDACTED >>>>')
polydat <- subset(polydat, shapefile == shp)

polydat <- subset(polydat, point == 0)
subset <- copy(polydat)

shape_master <- readRDS()
generated_pts <- list()
  
for (loc in unique(subset$location_code)) {
  shape <- shape_master[shape_master$GAUL_CODE == loc,]
  subset_loc2 <- subset(subset, location_code == loc)
  
  for (q in 1:nrow(subset_loc2)) {
    
    subset_loc3 <- subset_loc2[q,]
    
    year <- subset_loc3$start_year
    if (year <= 2000) {
      pop_raster <- raster('<<<< FILEPATH REDACTED >>>>', band = 1)
    } else {
      if (year > 2000 & year <= 2005) {
        pop_raster <- raster('<<<< FILEPATH REDACTED >>>>', band = 2)
      } else {
        if (year > 2005 & year <= 2010) {
          pop_raster <- raster('<<<< FILEPATH REDACTED >>>>', band = 3)
        } else {
          pop_raster <- raster('<<<< FILEPATH REDACTED >>>>', band = 4)
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
    
    names(samp_pts) <- c("long", "lat","weight")
    samp_pts$shapefile <- shp
    
    subset_loc3 <- left_join(samp_pts, subset_loc3, by = 'shapefile')
    subset_loc3$point <- 0
    
    generated_pts[[length(generated_pts) + 1]] <- subset_loc3
  }
    
}

generated_pts2 <- do.call(rbind, generated_pts)
if (!(run_date %in% list.files())) {dir.create(paste0(run_date))}
setwd(as.character(run_date))
write.csv(generated_pts2, file = paste0(shp,'.csv'), row.names = F)

