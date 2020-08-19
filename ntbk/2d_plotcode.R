# Load libraries
libs <- c('ggplot2', 'sf', 'raster', 'dplyr')


# Load shapefile and convert to sf object
# a2 shapefile
my_shp <- shapefile()
my_sf <- st_as_sf(shp)
#a0_shp
my_shp_a0 <- shapefile()
my_sf_a0 <- st_as_sf(my_shp_a0)

# Load data
my_data <- read.csv()

# map uncertainty and mean to a single vector
my_data$col_var <- ifelse(my_data$mean <= 0.25, 1,
	ifelse(my_data$mean <= 0.5, 2,
		ifelse(my_data$mean <= 0.75, 3, 4)))

my_data$ui_col <- ifelse(my_data$ui_width <= 0.25, 1,
	ifelse(my_data$ui_width <= 0.5, 2,
		ifelse(my_data$ui_width <= 0.75, 3, 4)))

my_data <- mutate(my_data, col_var = paste0(col_var, ui_col))

# merge data onto shapefile
my_sf <- left_join(my_sf, my_data, by = 'ADM2_CODE')

# plot
cols <- c('11' = '#89B5DE', '12'= '#B1CAE6' , '13'='#CCD7EA', '14'='#F5F0F4', '21'='#CCD7EA',
          '22'='#ACAFCD', '23'='#CABED1', '24'= '#EFD3DA', '31'= '#828EBC', '32'= '#9E92B6',
          '33'= '#C4A3BA', '34'= '#EAA9B2', '41'= '#7978A8', '42'= '#9F81A5', '43'='#C0829B',

gg_y_end <- ggplot() +
  # Map all countries with gray as a canvas
 geom_sf(data = my_sf_a0, color = 'black', fill = 'white', lwd = 0.05) +
  # Map data; "mean" is the variable you want to map
  geom_sf(data = my_data, aes(fill = col_var), lwd = 0) +
  # Map country borders
  geom_sf(data = my_sf_a0, color = 'black', fill = NA, lwd = 0.05) +
  # Use viridis color scale to map the variable
  scale_fill_manual(
    values = cols,
    na.value = 'grey') +
  theme_classic() +
  # # Makes the map pretty
  theme(legend.position = 'none',
        plot.margin=unit(c(0, 0, 0, 0), "in")) +
  coord_sf(datum=NA)
