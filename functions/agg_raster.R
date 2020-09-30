# Convert raster to numeric vector
raster_numeric <- function(raster, 
						   na.rm = TRUE) {

	# Load required libraries quietly
	libs <- c('sp', 'raster')
	lapply(libs, library, character.only = TRUE, quietly = TRUE)

	# Convert raster to a numeric vector
	output <- as.numeric(as.matrix(raster))
	
	if (na.rm) {
		# Subset to non NA values
		output <- output[!is.na(output)]
	}
}

# Aggregate raster by a shape
agg_raster <- function(raster,
					   shapefile,
					   field,
					   id,
					   weights = NULL) {

	# Load required libraries quietly
	libs <- c('sp', 'raster', 'sf', 'fasterize')
	lapply(libs, library, character.only = TRUE, quietly = TRUE)

	# Subset and convert shapefile to sf object
	shapefile <- shapefile[shapefile[[field]] == id,]
	poly_sf <- st_as_sf(shapefile)
	
	#Create weights raster to use in aggregation
	if(is.null(weights)) {
		wts <- raster
		wts[!is.na(wts)] <- 1
	} else {
		wts <- weights
	}
	
	# Rasterize shapefile
	shp_raster <- fasterize(poly_sf, wts, field = field)  
	shp_raster[!is.na(shp_raster)] <- 1
	# Calculate vector of weighted average 
	raster <- 	std_rasters(raster, shp_raster)
	numerator <- sum(raster_numeric(shp_raster * raster * wts), 
				   na.rm = TRUE)
	denominator <- sum(raster_numeric(shp_raster * wts), 
					 na.rm = TRUE)
	if (is.null(weights)) {
		denominator <- 1
	}
	output <- (numerator/denominator)

	return(output)
}