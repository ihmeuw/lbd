########## Examination of LF prevalence in areas considered by ESPEN to be medium-high risk for loiasis
##### Basic setup
user <- Sys.info()[["user"]] ## Get current user name
core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- sprintf(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))
package_list <- c(package_list, "sf")
path <- paste0(<<<< FILEPATH REDACTED >>>>)

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(sf)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

#### Loa setup
loa_shape <- readOGR(<<<< FILEPATH REDACTED >>>>)
loa_iu <- fread(<<<< FILEPATH REDACTED >>>>)

plot(loa_shape)
loa_iu <- loa_iu[Endemicity %in% c("Hyper-endemic", "Meso-endemic")]
table(loa_iu$Country)

#### Load spatial templates
load(<<<< FILEPATH REDACTED >>>>)

#### Calculate population of loa-endemic localities
# Rasterize loa shapefile
loa_shape <- loa_shape[loa_shape$csv_2015 %in% c("Hyper-endemic", "Meso-endemic"),]
loa_shape_raster <- rasterize(loa_shape, pop_raster)
loa_pop <- raster::mask(pop_raster, loa_shape_raster)
loa_pop <- raster::subset(loa_pop, 19)
loa_pop_total <- sum(values(loa_pop), na.rm = TRUE)

#### Calculate cases
lf_prev_2018 <- raster(<<<< FILEPATH REDACTED >>>>)
lf_prev_2018_loa <- raster::mask(lf_prev_2018, loa_pop)
lf_cases_loa <- loa_pop * lf_prev_2018_loa
loa_cases_total <- sum(values(lf_cases_loa), na.rm = TRUE)

####### Identify data from loa-endemic (hyper or meso) locations
df <- fread(<<<< FILEPATH REDACTED >>>>)
spdf <- SpatialPointsDataFrame(coords = df[, c("longitude", "latitude")], data = df, proj4string = crs(subset_shape))
df_loa <- over(spdf, loa_shape)
df <- cbind(df, df_loa)
df <- df[!is.na(csv_2015)]
df <- df[diagnostic %in% c("ICT", "FTS")]
unique_Master_UID <- unique(df$Master_UID)
length(unique_Master_UID)
table(df$country)
aggregate(Master_UID ~ diagnostic, df, length)
unique_diag <- unique(df[, c("Master_UID", "diagnostic")])
table(unique_diag$diagnostic)
df <- df[!(country %in% c("BFA", "MDG", "UGA"))] # drop countries that are not loiasis-endemic
unique_diag <- unique(df[, c("Master_UID", "diagnostic")])
table(unique_diag$diagnostic)
