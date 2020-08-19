# Set up the environment
rm(list = ls())
libs <- c('sp', 'raster', 'rgeos', 'rgdal', 'data.table', 'ggplot2',
	'viridis', 'dplyr')
lapply(libs, library, character.only = TRUE)


# User inputs
## model indicator
indicator <- 'w_piped'
outdir <- "<<<< FILEPATH REDACTED >>>>"

# Custom functions
grab_input <- function(indi) {
	data_dir <- "<<<< FILEPATH REDACTED >>>>"
	return(fread(paste0(data_dir, indi, '.csv')))
}


# read in and format shapefile
myshp <- readRDS("<<<< FILEPATH REDACTED >>>>")
myshp$ADM0_CODE <- as.numeric(as.character(myshp$ADM0_CODE))
myshp$ADM0_NAME <- as.character(myshp$ADM0_NAME)
mysf <- st_as_sf(myshp)
mysf <- rename(mysf, location_name = ADM0_NAME) %>%
			select(location_name)

# Read and format lookup table for iso3 to country names
lookup_file <- fread("<<<< FILEPATH REDACTED >>>>")
lookup_file <- select(lookup_file, location_name, iso3) %>%
	rename(country = iso3)

# Read and format input data
mydat <- grab_input('w_piped')
clean_dat <- mydat %>%
	group_by(country) %>%
	summarize(n = n_distinct(nid))
clean_dat <- left_join(clean_dat, lookup_file, by = 'country')

# Merge input data onto shapefile and plot
mysf <- left_join(mysf, clean_dat, by = 'location_name')

gg1 <- ggplot() +
	geom_sf(data = mysf, aes(fill = n), lwd = 0.005) +
	scale_fill_viridis(name = 'Number of NIDs') +
	ggtitle(indicator) +
	theme_bw() +
	coord_sf()

pdf(paste0(outdir, indicator, '_', Sys.Date(), '.pdf'),
	width = 11, height = 8)
print(gg1)
dev.off()
