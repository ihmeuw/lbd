core_repo <- "<<<< FILEPATH REDACTED >>>>"
ig_repo <- "<<<< FILEPATH REDACTED >>>>"
run_date <- "<<<< FILEPATH REDACTED >>>>"
shapefile_version <- "<<<< FILEPATH REDACTED >>>>"

output_path <- "<<<< FILEPATH REDACTED >>>>"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                              SETUP
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(paste0(core_repo, '/mbg_central/setup.R'))
pl <- c("rgeos",      "data.table", "raster",     "rgdal",      "INLA",
        "seegSDM",    "seegMBG",    "dismo",      "gbm",        "foreign",
        "parallel",   "doParallel", "grid",       "gridExtra",  "pacman",
        "gtools",     "glmnet",     "ggplot2",    "RMySQL",     "plyr",
        "tictoc",     "dplyr",      "magrittr",   "tidyr",      "sp",
        "sf",         "matrixStats", "fasterize")
mbg_setup(package_list = pl, repos = core_repo)

admin_shp <- sf::st_read(get_admin_shapefile(admin_level = 2, version = shapefile_version), quiet=T)
admin_shp_data <- as.data.table(admin_shp)

#prep data csvs
for(measure in c("Q", "D")){
  for(age in c("under5", "infant", "neonatal")){
    for(lvl in c(0,1,2)){
      adm <- fread("<<<< FILEPATH REDACTED >>>>", drop = "V1")
      adm[, sex_id := 3][, sex := "Both"]
      
      if(age == "under5") {
        adm[, age_group_id := 1][, age_group_name := "Under 5"]
      } else if(age == "infant"){
        adm[, age_group_id := 28][, age_group_name := "Infant"]
      } else {
        adm[, age_group_id := 42][, age_group_name := "Neonatal"]
      }
      
      if(measure == "Q"){
        adm[, measure_id := 27][, measure := "Probability of death"]
      } else {
        adm[, measure_id := 1][, measure := "Deaths"]
      }
      
      if(lvl == 0) {
        loc_merge <- "ADM0_CODE"
        loc <- unique(admin_shp_data[, .(ADM0_CODE, ADM0_NAME)])
        com <- merge(loc, adm, by = loc_merge)
        com <- com[,.(ADM0_NAME, ADM0_CODE, year, age_group_id, age_group_name, sex_id, sex, measure_id, measure, mean, lower, upper)]
      } else if (lvl == 1) {
        loc_merge <- "ADM1_CODE"
        loc <- unique(admin_shp_data[, .(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME)])
        com <- merge(loc, adm, by = loc_merge)
        com <- com[,.(ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, year, age_group_id, age_group_name, sex_id, sex, measure_id, measure, mean, lower, upper)]
      } else {
        loc_merge <- "ADM2_CODE"
        loc <- unique(admin_shp_data[, .(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME)])
        com <- merge(loc, adm, by = loc_merge)
        com <- com[,.(ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE, year, age_group_id, age_group_name, sex_id, sex, measure_id, measure, mean, lower, upper)]
        #removing india from Adm2 estimates
        com <- com[ADM0_CODE != 105,]
      }
      
      write.csv(com, sprintf("%s/IHME_LMICS_U5M_2000_2017_%s_%s_ADM%i_Y2019M09D17.CSV", output_path, measure, toupper(age), lvl), row.names = F)
    }
  }
}

#load in masks
mask_folder <- "<<<< FILEPATH REDACTED >>>>"
lakes <- raster("<<<< FILEPATH REDACTED >>>>")
pop   <- raster("<<<< FILEPATH REDACTED >>>>")
extent_example <- brick("<<<< FILEPATH REDACTED >>>>")

e2 <- new.env()
load(paste0("<<<< FILEPATH REDACTED >>>>"), e2)
simple_raster <- get('simple_raster', e2)
rm(e2)

simple_raster[simple_raster != 105] <- NA
simple_raster[simple_raster == 105] <- 1


#prep masks to fit rasters
lakes <- crop(lakes, extent(extent_example))
lakes <- setExtent(lakes, extent(extent_example))
pop <- crop(pop, extent(extent_example))
pop <- setExtent(pop, extent(extent_example))
#each lake in this raster is assigned a number, changing them all to 1 so raster::mask works
lakes[lakes > 0] <- 1

#removing india from pixel level estimates
ind <- extend(simple_raster, extent(extent_example))
ind <- crop(ind, extent(extent_example))
ind <- setExtent(ind, extent(extent_example))
              

year_list <- 2000:2017

#make raster bricks 
for(measure in c("Q", "D")){
  for(age in c("under5", "infant", "neonatal")){
    for(stat in c("mean", "upper", "lower")) {
      adm <- brick(paste0("<<<< FILEPATH REDACTED >>>>")
      
      adm <- mask(adm, lakes, maskvalue = 1)
      adm <- mask(adm, pop, maskvalue = 1)
      adm <- mask(adm, ind, maskvalue = 1)
      
      save_filename <- sprintf("<<<< FILEPATH REDACTED >>>>")
      writeRaster(adm, "<<<< FILEPATH REDACTED >>>>", format = "GTiff", overwrite = T)
    }
  }
}

#make data inclusion sheet - this requires some work outside of R
extras <- read.csv("<<<< FILEPATH REDACTED >>>>")
inclusion <- read.csv("<<<< FILEPATH REDACTED >>>>")
ghdx_meta <- ghdx_construct_pub_table(
  nids         = inclusion$NID,
  core_repo    = core_repo,
  resolve_urls = T
)
ghdx_meta[, citation := NULL]

data <- merge(ghdx_meta, extras, by="nid")
data <- merge(data, inclusion, by.x = "nid", by.y = "NID")
write.csv(data, paste0("<<<< FILEPATH REDACTED >>>>"))

