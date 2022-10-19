###### Produce polygon resampling figure for oncho paper
###### For a selected country, show panels: (a) Locations of polygon data, colored by prevalence; (b) WorldPop population raster; (c) locations of resampling points; (d) final resampled points, colored by prevalence

user <- Sys.info()[["user"]]

## Set repo locations
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

path <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(matrixStats)
library(sf)

path <- <<<< FILEPATH REDACTED >>>>
library(ggspatial, lib.loc=path)
library(cowplot)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

##########
### Select country to display: Cameroon, 2015
# Get prevalence data
df <- fread(<<<< FILEPATH REDACTED >>>>)
df <- df[country == "CMR" & year == 2015]

# Produce country shapefile
modeling_shapefile_version <- "2020_02_20"
stage_master_list <- fread(<<<< FILEPATH REDACTED >>>>)
subset_shape <- load_simple_polygon(gaul_list = stage_master_list[iso3 == "CMR", gadm_geoid], buffer = 1, tolerance = 0.4, use_premade = FALSE, shapefile_version = modeling_shapefile_version)[[1]]

## Build administrative and population rasters
raster_list <- build_simple_raster_pop(subset_shape)
pop_raster <- raster_list[["pop_raster"]]
simple_raster <- raster_list[["simple_raster"]]

## Mask aggregated administrative raster by aggregated population raster
simple_raster <- raster::mask(simple_raster, pop_raster[[1]])

# Retrieve WorldPop data
worldpop_current <- data.table(covariate = "worldpop", measure = "total", release = "2020_03_20", year_lag = 0)

interval_mo <- 12
raster_agg_factor <- 1

pop <- load_lagged_covariates(covariate_config = worldpop_current,
                              template = simple_raster,
                              start_year = 2015,
                              end_year = 2015)$worldpop

## Load custom shapefiles
shapefiles <- unique(df[!is.na(shapefile), shapefile])
shapes <- list()
ref_data <-  vector(mode = "list", length = nrow(unique(df[!is.na(shapefile), list(shapefile, location_code)])))
count <- 0
for (i in shapefiles) {
  ids <- unique(df[shapefile == i, location_code])
  current <- readRDS(<<<< FILEPATH REDACTED >>>>)
  for (j in ids) {
    count <- count + 1
    shapes <- c(shapes, current[current@data$GAUL_CODE == j,])
    ref_data[[count]] <- data.table("shape" = i, "code" = j)
  }
}

setnames(df, "weight", "Weight")
subset_shape@data$id <- rownames(subset_shape@data)
subset_shape.df <- fortify(subset_shape, region = "id")
subset_shape.df <- left_join(subset_shape.df, subset_shape@data, by="id")
subset_shape.df$Prevalence <- as.numeric(NA)

### Panel a: Cameroon map, with inset
gg <- ggplot() + theme_classic() + coord_equal() + theme(plot.margin = unit(c(0, 60, 6, 0), "points")) + geom_polygon(data=subset_shape.df, aes(long, lat, group = group), fill = NA, color = "#666666") + labs(x = "Longitude", y = "Latitude")
gg <- gg + geom_rect(aes(xmin = 9, xmax = 12, ymin = 3.5, ymax = 6), color = "red", fill = NA)
gg
panel_a <- copy(gg)

### Panel b: Plot custom polygons
gg <- ggplot() + theme_classic() + coord_equal() + scale_fill_gradient(low = "olivedrab2", high = "mediumpurple1", na.value = "white") + geom_polygon(data=subset_shape.df, aes(long, lat, group = group, fill = Prevalence), color = "#666666") + coord_cartesian(xlim = c(9, 12), ylim = c(3.5, 6)) + labs(x = "Longitude", y = "Latitude")

for (i in 1:length(shapes)) {
  current <- shapes[[i]]
  current@data$id <- rownames(current@data)
  current.df <- fortify(current, region = "id")
  current.df <- left_join(current.df, current@data, by="id")
  current.df$Prevalence <- df[shapefile == ref_data[[i]]$shape & location_code == ref_data[[i]]$code, had_oncho_w_resamp / N][1]
  # print(prev)
  gg <- gg + geom_polygon(data=current.df, aes(long, lat, group = group, fill = Prevalence), color = "#666666")
}

# Add point data

df[, "Prevalence" := had_oncho_w_resamp / N]
gg <- gg + geom_point(data=df[is.na(shapefile)], aes(x = longitude, y = latitude, fill = Prevalence), pch = 21, size = 2, color = "#666666")
gg

panel_b <- copy(gg)

#### Panel c: Plot population raster
pop_spdf <- as(pop, "SpatialPixelsDataFrame")
pop_df <- as.data.frame(pop_spdf)
colnames(pop_df) <- c("Population", "x", "y")

gg <- ggplot() + theme_classic() + coord_equal() + geom_tile(data=pop_df, aes(x=x, y=y, fill=Population), alpha=0.8) + theme(legend.position="right") + coord_cartesian(xlim = c(9, 12), ylim = c(3.5, 6)) + scale_fill_gradient(low = "white", high = "black", na.value = "white", trans = "log10") + geom_polygon(data=subset_shape.df, aes(long, lat, group = group), color = "#666666", fill = NA) + labs(x = "Longitude", y = "Latitude")

for (i in 1:length(shapes)) {
  current <- shapes[[i]]
  current@data$id <- rownames(current@data)
  current.df <- fortify(current, region = "id")
  current.df <- left_join(current.df, current@data, by="id")
  current.df$Prevalence <- df[shapefile == ref_data[[i]]$shape & location_code == ref_data[[i]]$code, had_oncho_w_resamp / N][1]
  # print(prev)
  gg <- gg + geom_polygon(data=current.df, aes(long, lat, group = group), color = "#000000", fill = NA)
}

gg

panel_c <- copy(gg)

#### Panel d: Plot resampled data
# resampled_data <- df[!is.na(shapefile)]
gg <- ggplot() + theme_classic() + coord_equal() + scale_fill_gradient(low = "olivedrab2", high = "mediumpurple1", na.value = "white") + geom_polygon(data=subset_shape.df, aes(long, lat, group = group, fill = Prevalence), color = "#666666") + coord_cartesian(xlim = c(9, 12), ylim = c(3.5, 6))
gg <- gg + geom_point(data=df, aes(x = longitude, y = latitude, fill = Prevalence, alpha = Weight), pch = 21, size = 2, color = "#666666") + labs(x = "Longitude", y = "Latitude")
gg

panel_d <- copy(gg)

plot_grid(panel_a, panel_b, panel_c, panel_d, labels = c('a', 'b', 'c', 'd'), label_size = 24, scale = 0.9)
