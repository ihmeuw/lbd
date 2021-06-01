####################################################################################################
## Create a look up table to relate different subnationals
####################################################################################################



library(Hmisc)

library(data.table)
library(rgeos)
library(rgdal)
library(maptools)
library(raster)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(grid)
library(dplyr)


rm(list = ls())

## Settings ----------------------------------------------------------------------------------------

# Get Functions
## Set repo
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")

setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo,  '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
mbg_setup(package_list = package_list, repos = c(indic_repo, core_repo))

## Set Country and non standard shp ----------------------------------------------------------------------------

c <- "NER"
iso3 <- c("NER")
other <- readOGR("<<<< FILEPATH REDACTED >>>>")

## Load shps ----------------------------------------------------------------------------

standard <- readOGR("<<<< FILEPATH REDACTED >>>>")
code <- get_adm0_codes(iso3)
standard <- standard[which(standard@data$ADM0_CODE == code), ]

shapefile_version = 'epp_shapes'
modeling_shapefile_version = shapefile_version
## Load population for weights ---------------------------------------------------------------------

#load regional simple raster
simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(paste0(c)), buffer = 0.4, subset_only = FALSE, shapefile_version = "epp_shapes")
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
message("Building simple raster from subset_shape")
raster_list    <- build_simple_raster_pop(subset_shape, link_table = NULL)
simple_raster  <- raster_list[['simple_raster']]

# get world pop
pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',
                                                     measures = 'a1549t',
                                                     simple_polygon = simple_polygon,
                                                     start_year  = 2016,
                                                     end_year    = 2016,
                                                     interval_mo = 12,
                                                     agebin = 1)

# extend and crop pop raster to ensure it matches the simple raster #not convinced this is totally needed
pop <- pop_raster_annual[[1]]
pop  <- extend(pop, simple_raster, values = NA)
pop  <- crop(pop, extent(simple_raster))
pop  <- setExtent(pop, simple_raster)
pop  <- raster::mask(pop, simple_raster)

other <- spTransform(other, crs(pop))

## check to ensure the pop raster matches the simple raster in extent and resolution
if (extent(pop) != extent(simple_raster)) {
  stop("population raster extent does not match simple raster")
}
if (any(res(pop) != res(simple_raster))) {
  stop("population raster resolution does not match simple raster")
}
pop_df <- data.table(rasterToPoints(pop))
setnames(pop_df, c("long", 'lat', 'pop'))
setkey(pop_df, long, lat)
pop_df <- as.data.frame(pop_df)
pop_df$long <- round(pop_df$long, 6)
pop_df$lat <- round(pop_df$lat, 6)

## Load geographies ----------------------------------------------------------------------------------

raster_other <- raster::rasterize(other, simple_raster, field = "OBJECTID_1")
raster_other <- data.table(rasterToPoints(raster_other))
setnames(raster_other, c("long", 'lat', "other_id"))
setkey(raster_other, long, lat)
raster_other <- as.data.frame(raster_other)
raster_other$long <- round(raster_other$long, 6)
raster_other$lat <- round(raster_other$lat, 6)


raster_standard <- raster::rasterize(standard, simple_raster, field = "geo_id")
raster_standard <- data.table(rasterToPoints(raster_standard))
setnames(raster_standard, c("long", 'lat', "ADM2_CODE"))
setkey(raster_standard, long, lat)
raster_standard <- as.data.frame(raster_standard)
raster_standard$long <- round(raster_standard$long, 6)
raster_standard$lat <- round(raster_standard$lat, 6)



full <- merge(raster_standard, raster_other, by = c("lat", "long"))
full <- merge(full, pop_df, by = c("lat", "long"))

full <- as.data.table(full)
full_c <- full[,lapply(c("pop"), function(x) sum(get(x))), by = c("ADM2_CODE", "other_id")]


look_up <- list()
for (id in c(1:length(standard@data$ADM2_CODE))) {
  c <- standard@data$geo_id[[id]]
  sub <- full_c[which(full_c$ADM2_CODE == c), ]
  r <- sub[which(sub$V1 == max(sub$V1)), ]
  c <- as.character(c)
  look_up[[c]] <- r
}

look_up_final <- bind_rows(look_up, .id = "id")

write.csv(look_up_final, file = paste0(indic_repo,"<<<< FILEPATH REDACTED >>>>"))


