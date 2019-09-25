rm(list = ls())
# Set library and load packages
root <- '<<<< FILEPATH REDACTED >>>>' #filepath root
## drive locations
package_list <- c() #list of packages
repo <- '<<<< FILEPATH REDACTED >>>>' #Code code repo
setwd(repo) #directory to write resampled CSVs

# Load MBG packages and functions
message("Loading in required R packages and MBG functions")
source(paste0(repo, "setup.R"))
mbg_setup(package_list = package_list, repos = repo)

#grab system arguments
shp <- commandArgs()[6]
indic <- commandArgs()[7]
run_date <- commandArgs()[8]
latest_collapse <- commandArgs()[9]

#read in ready to resample data
polydat <- '<<<< FILEPATH REDACTED >>>>'

#make sure you are only working with polygon data
polydat <- subset(polydat, point == 0)
#grab only data for the shapefile of the job
subset <- polydat[which(polydat$shapefile == shp),]

#read in the shapefile
shape_master <- readRDS('<<<< FILEPATH REDACTED >>>>')

#resample into points
generated_pts <- list()
for (loc in unique(subset$location_code)) {
  shape <- shape_master[shape_master$GAUL_CODE == loc,]
  subset_loc2 <- subset(subset, location_code == loc)

  for (q in 1:nrow(subset_loc2)) {

    subset_loc3 <- subset_loc2[q,]

    #choose population raster by year bins
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

    subset_loc3 <- left_join(samp_pts, subset_loc3, by = "shapefile")
    subset_loc3$point <- 0

    generated_pts[[length(generated_pts) + 1]] <- subset_loc3
  }

}

#write out resampled CSV
generated_pts2 <- do.call(rbind, generated_pts)
if (!(run_date %in% list.files())) {dir.create('<<<< FILEPATH REDACTED >>>>')}
setwd('<<<< FILEPATH REDACTED >>>>')
write.csv(generated_pts2, file = '<<<< FILEPATH REDACTED >>>>', row.names = F)

