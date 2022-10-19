########## Examination of LF or onchoprevalence in areas considered by ESPEN to be medium-high risk for loiasis
##### Basic setup
user <- Sys.info()[["user"]] ## Get current user name
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))
package_list <- c(package_list, "sf")
path <- paste0(<<<< FILEPATH REDACTED >>>>)

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(sf)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

#### Loa setup
loa_shape <- readOGR(<<<< FILEPATH REDACTED >>>>)
loa_iu <- fread(<<<< FILEPATH REDACTED >>>>)

plot(loa_shape)
loa_iu <- loa_iu[Endemicity %in% c("Hyper-endemic", "Meso-endemic")]
table(loa_iu$Country)

#### Load spatial templates
load(<<<< FILEPATH REDACTED >>>>)

#### Get worldpop
pop_raster <- load_worldpop_covariate(template_raster = simple_raster,
                                      covariate = "worldpop",
                                      pop_measure = "total",
                                      pop_release = "2020_03_20",
                                      start_year = 2018,
                                      end_year = 2018,
                                      interval = 12)$worldpop


#### Calculate population of loa-endemic localities
# Rasterize loa shapefile
loa_shape$Composite_IU_ID <- paste0("A", loa_shape$IU_ID)
loa_shape <- merge(loa_shape, loa_iu, by = "Composite_IU_ID", all = TRUE)
loa_shape <- loa_shape[loa_shape$Endemicity %in% c("Hyper-endemic", "Meso-endemic"),]
loa_shape_raster <- rasterize(loa_shape, pop_raster)
loa_pop <- raster::mask(pop_raster, loa_shape_raster)
loa_pop_total <- sum(values(loa_pop), na.rm = TRUE)

#### Calculate cases
ind_prev_2018 <- raster(<<<< FILEPATH REDACTED >>>>)
ind_prev_2018_in_loa <- raster::mask(ind_prev_2018, loa_pop)
ind_cases_in_loa <- loa_pop * ind_prev_2018_in_loa
cases_total <- sum(values(ind_cases_in_loa), na.rm = TRUE)
