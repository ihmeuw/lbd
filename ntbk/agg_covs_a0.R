# user inputs
user            <- "<<<< USERNAME REDACTED >>>>"
reg <- 'tls'
modeling_shapefile_version <- '2019_09_10'
pop_release <- '2020_03_20'

# set up environment
repo            <- "<<<< FILEPATH REDACTED >>>>"
root <- "<<<< FILEPATH REDACTED >>>>"
commondir <- "<<<< FILEPATH REDACTED >>>>"
package_list <- c(t(read.csv(sprintf("%s/package_list.csv", commondir), header = FALSE)))

setwd(repo)
core_repo <- repo
for (p in package_list) {
  try(library(p, character.only = T))
}
library(seegSDM)
library(seegMBG)
library(mgcv)
library(sf)
library(raster)

# Load MBG packages and functions
message("Loading in required R packages and MBG functions")
source(paste0(repo, "/mbg_central/setup.R"))
mbg_setup(package_list = package_list, repos = repo)

## Load simple polygon template to model over
gaul_list <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
simple_polygon_list <- load_simple_polygon(
gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
shapefile_version = modeling_shapefile_version
)
subset_shape <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
## Load list of raster inputs (pop and simple)
raster_list <- build_simple_raster_pop(subset_shape, pop_release = pop_release)
simple_raster <- raster_list[["simple_raster"]]
pop_raster <- raster_list[["pop_raster"]]

cov_layers <- load_and_crop_covariates_annual(
  covs = c('access2', 'aridity', 'elevation', 'ghslurbanicity',
  	'growingseason', 'irrigation'),
  measures = rep('mean', 6),
  simple_polygon = simple_raster,
  start_year = 2000,
  end_year = 2017,
  interval_mo = 12
)

for (ii in 1:length(cov_layers)) {
	layers <- nlayerr(cov_layers[[ii]])
	for (jj in 1:layers) {
		sum(cov_layers[[ii]][[jj]] * pop_raster[[jj]])/sum(pop_raster[[jj]])
	}
}
