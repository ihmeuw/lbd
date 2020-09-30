# Take two rasters and standardize their extent and non-NA cells
std_rasters <- function(grid, standard) {

	# Load required libraries quietly
	libs <- c('sp', 'raster')
	lapply(libs, library, character.only = TRUE, quietly = TRUE)

	# Check if grids is in the correct class, otherwise error out
	ok_class <- c('RasterBrick', 'RasterLayer')
	if (!(class(grid) %in% ok_class)) {
		stop("grids must be a RasterBrick or a RasterLayer.")
	}

	# Check if standard is a rasterlayer, otherwise error out
	if (class(standard) != 'RasterLayer') {
		stop("standard must be a RasterLayer.")
	}

	# Standardize raster
    output  <- mask(
    	  	 	 setExtent(
    			   crop(grid, extent(standard)),
    			   standard),
    			 standard)

    return(output)
}
