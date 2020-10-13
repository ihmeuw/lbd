
## define function for producing Onchocerciasis bias grid for background and polygon sampling

## IU_shapefile: shapefile of implementation units, a SpatialPolygonsDataFrame
## custom:       shapefile OR raster/GeoTIFF of custom polygons to be excluded from background and polygon sampling
## base:         shapefile OR raster/GeoTIFF defining extent of sampling - inverse of NA region
## field:        IU_shapefile@data field name for criteria of exclusion from background and polygon sampling
## exclude:      values within field to exclude from background and polygon sampling
## ref:          reference raster for rasterizing shapefiles where necessary

## example inputs:

# ref <- raster('FILEPATH.tif')
# IU_shapefile <- shapefile('FILEPATH.shp')
# base <- shapefile("FILEPATH.shp")
# field <- 'Endem_MDA'
# exclude <- c("Endemic (under MDA)", "Endemic (MDA not started)", "Endemic, Post-MDA Surveillance")

make_bias_grid <- function(IU_shapefile, custom=NULL, base, field, exclude, ref){
  require(raster)
  require(rgdal)
  # check for map projection and change to CRS of base if needed
  crs_base <- proj4string(base)
  crs_IU <- proj4string(IU_shapefile)
  if (crs_IU!=crs_base) IU_shapefile <- spTransform(IU_shapefile, CRS(crs_base))
  # add raster value to include base in background and polygon sampling
  base@data <- cbind(base@data, rast_val=rep(0, nrow(base)))
  # exclude indicated IUs
  IU_exclude <- c()
  for (i in 1:length(exclude)){
    IU_exclude <- append(IU_exclude, which(IU_shapefile@data[,field] == exclude[i]))
  }
  # make binary-valued raster of IUs for 1=exclude, 0=include
  IU_shapefile$rast_val <- c(rep(0, nrow(IU_shapefile)))
  IU_shapefile$rast_val[IU_exclude] <- 1
  IU_raster <- rasterize(x=IU_shapefile, y=ref, field=IU_shapefile$rast_val)
  # include base region and rasterize
  IU_shapefile_base<-bind(base, IU_shapefile)
  IU_raster <- rasterize(x=IU_shapefile_base, y=ref, field=IU_shapefile_base$rast_val)
  # if given a custom shapefile/SPDF, make a raster of 1's; combine rasters to form bias_grid
  if (!is.null(custom)){
    if (class(custom)=='SpatialPolygonsDataFrame'){
      custom@data <- cbind(custom@data, rast_val=rep(1, nrow(custom)))
      custom <- rasterize(x=custom, y=ref, field=custom@data$rast_val)
    }
    bias_grid <- IU_raster + custom
  } else {
    bias_grid <- IU_raster
  }
  # force overlapping areas to 1
  bias_grid[bias_grid > 1] <- 1
  # invert bias grid to bias away from excluded regions 
  #  -- bg_sampling() bias parameter is distribution => make "square" distribution
  bias_grid <- 1 - bias_grid

  return(bias_grid)
}
