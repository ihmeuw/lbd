indicator       = as.character(commandArgs()[3])
indicator_group = as.character(commandArgs()[4])
run_date        = as.character(commandArgs()[5])
repo            = as.character(commandArgs()[6])
baseline_year   = as.numeric(commandArgs()[7])
year_to_map     = as.numeric(commandArgs()[8])

##################################################################################
## Setup.
##################################################################################
message('Setting up...')

## Define main directories.
results_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/')
val_dir <- paste0(results_dir, 'validation_report')
dir.create(val_dir, showWarnings = FALSE)

## Load libraries and miscellaneous MBG project functions.
setwd(repo)
root <- "<<<< FILEPATH REDACTED >>>>>"
package_lib <- paste0(root,'<<<< FILEPATH REDACTED >>>>>') 
.libPaths(package_lib)                                  

source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
source('mbg_central/post_estimation_functions.R')
source('mbg_central/gbd_functions.R')
source('mbg_central/shiny_functions.R')
source('mbg_central/holdout_functions.R')
source('education/pop_splitting_functions.R')
source('mbg_central/seegMBG_transform_functions.R')
source('mbg_central/validation_report_functions.R')
package_list <- c('fields', 'gridGraphics' ,'grid', 'gridExtra',
                  'gstat', 'magrittr', 'ggplot2', 'doParallel',
                  'SDMTools', 'foreign', 'rgeos', 'data.table',
                  'raster','rgdal','INLA','seegSDM','seegMBG','plyr',
                  'dplyr','leaflet')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

for(admin_level in c(0,1,2)) {
  if(admin_level %in% c(0,1)) shapes <- shapefile(paste0("<<<<< FILEPATH REDACTED >>>>>/GAUL_admin/admin", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, ".shp"))
  if(admin_level == 2) shapes <- shapefile(paste0("<<<<< FILEPATH REDACTED >>>>>/GAUL_admin/admin", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, "_modified.shp"))
  ## Fix several african disputed territories (assign them to a country so we can include them in our estimates and they dont drop out)
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==61013]=133
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40760]=40765
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40762]=145
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==102  ]=6
  assign(paste0('admin_shape_', admin_level), shapes)
}

results_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/table_', baseline_year)
all_files <- list.files(results_dir, pattern = measure, full.names = TRUE)
all_admins <- rbindlist(lapply(all_files, fread))
all_admins <- all_admins[year == year_to_map, ]
all_admins <- all_admins[is.na(mean), mean := goal_threshold]
all_admins <- all_admins[is.na(upper), upper := goal_threshold]
all_admins <- all_admins[is.na(lower), lower := goal_threshold]
setnames(all_admins, 'admin2', 'ADM2_CODE')
setnames(all_admins, 'admin1', 'ADM1_CODE')
setnames(all_admins, 'admin0', 'ADM0_CODE')
gaul_to_loc <- fread('<<<< FILEPATH REDACTED >>>>>/gaul_to_loc_id.csv')
all_gauls <- unique(all_admins[, ADM0_CODE])

redblue <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac')
library(ggplot2)
library(ggthemes)

## Save copy of all_admins (by level) for Lucas
for(adm in c(0,1,2)) {
  if(adm == 0) {
    this_admin_data <- all_admins[is.na(ADM1_CODE) & is.na(ADM2_CODE), ]
  }
  if(adm == 1) {
    this_admin_data <- all_admins[!is.na(ADM1_CODE) & is.na(ADM2_CODE), ]
  }
  if(adm == 2) {
    this_admin_data <- all_admins[!is.na(ADM1_CODE) & !is.na(ADM2_CODE), ]
  }
  write.csv(this_admin_data, paste0(results_dir, '/pixels/adm', adm, '_', indicator, '_', measure, '_sdgprob', year_to_map, '_summary_table.csv'))
}

make_africa_raster <- function(gaul) {

  # Convert raster to SpatialPointsDataFrame
  pixel_probs <- fread(paste0(results_dir, '/pixels/mortality', gaul, '_probs.csv'))
  simple_raster <- raster(paste0(results_dir, '/simple/', gaul, '.tif'))
  probs_raster <- insertRaster(simple_raster, matrix(pixel_probs[, p_goal], ncol = 1))
  return(probs_raster)

}
africa_raster <- lapply(unique(all_admins[, ADM0_CODE]), make_africa_raster)
final_africa_raster = do.call(raster::merge, africa_raster)
writeRaster(final_africa_raster,
            file = paste0(results_dir, '/pixels/', indicator, '_', measure, '_sdgprob', year_to_map),
            format='GTiff',
            overwrite = TRUE)
