##############################################################################
## Prepare core diarrhea estimates for figures and viz tool
##
## Includes:
##  Counts and proportions/rates
##  Mean, upper, lower
##  Raster, Admin 0, Admin 1, Admin 2
##############################################################################


## Setup -------------------------------------------------------------------------

# clear environment
rm(list = ls())

# set general arguments
user            <- Sys.info()['user']
repo            <- '<<<< FILEPATH REDACTED >>>>'
indicator_group <- 'ort'

# Load functions
library(data.table)
library(raster)

# set indicator
indicators <- 'had_diarrhea'

# set indicator arguments
indicator_run_date <- '2019_09_17_14_12_53'
measures <- c('_prevalence', '_deaths', '_incidence')
rake <- '_raked'


## Start saving -------------------------------------------------------------------------
  
# loop over measures
for (measure in measures) {
  message(measure)
  
  # set directories
  in_dir <- '<<<< FILEPATH REDACTED >>>>'
  out_dir <- '<<<< FILEPATH REDACTED >>>>'
  out_dir2 <- '<<<< FILEPATH REDACTED >>>>'
  
  # create directories
  dir.create(out_dir)
  dir.create(out_dir2)
  
  
  ## Save rasters -------------------------------------------------------------------------
  
  # read in mean estimate raster files and merge
  g <- list.files(in_dir, pattern = glob2rx(paste0('*prediction', measure, '*.grd*')))
  filepaths <- paste0(in_dir, g)
  rasters <- lapply(filepaths, brick)
  global_raster <- do.call(raster::merge, rasters)
  
  # read in upper estimate raster files and merge
  g <- list.files(in_dir, pattern = glob2rx(paste0('*upper', measure, '*.grd*')))
  filepaths <- paste0(in_dir, g)
  rasters <- lapply(filepaths, brick)
  upper_global_raster <- do.call(raster::merge, rasters)
  
  # read in lower estimate raster files and merge
  g <- list.files(in_dir, pattern = glob2rx(paste0('*lower', measure, '*.grd*')))
  filepaths <- paste0(in_dir, g)
  rasters <- lapply(filepaths, brick)
  lower_global_raster <- do.call(raster::merge, rasters)
  
  # save mean, upper, and lower rasters for viz team
  writeRaster(global_raster, paste0(out_dir, indicator, '_mean', measure, '_estimate_raster.tif'), overwrite = TRUE)
  writeRaster(upper_global_raster, paste0(out_dir, indicator, '_upper', measure, '_estimate_raster.tif'), overwrite = TRUE)
  writeRaster(lower_global_raster, paste0(out_dir, indicator, '_lower', measure, '_estimate_raster.tif'), overwrite = TRUE)
  
  # save mean, upper, and lower rasters for figures
  writeRaster(global_raster, paste0(out_dir2, indicator, ifelse(measure == '_deaths', '_mortality', measure), '_mean', rake, '_2000_2017.tif'), overwrite = TRUE)
  writeRaster(upper_global_raster, paste0(out_dir2, indicator, ifelse(measure == '_deaths', '_mortality', measure), '_upper', rake, '_2000_2017.tif'), overwrite = TRUE)
  writeRaster(lower_global_raster, paste0(out_dir2, indicator, ifelse(measure == '_deaths', '_mortality', measure), '_lower', rake, '_2000_2017.tif'), overwrite = TRUE)
  
  # read in population raster files and crop
  g <- list.files('<<<< FILEPATH REDACTED >>>>')
  g <- g[1:18]
  filepaths <- paste0('<<<< FILEPATH REDACTED >>>>', g)
  pop_rasters <- lapply(filepaths, raster)
  pop_raster <- stack(pop_rasters)
  pop_raster <- crop(pop_raster, extent(global_raster))
  
  # number of children with diarrhea
  counts_raster <- round(global_raster*pop_raster)
  upper_counts_raster <- round(upper_global_raster*pop_raster)
  lower_counts_raster <- round(lower_global_raster*pop_raster)
  
  # save counts, upper, and lower rasters for viz team
  writeRaster(counts_raster, paste0(out_dir, indicator, '_mean', measure, '_counts_raster.tif'), overwrite = TRUE)
  writeRaster(upper_counts_raster, paste0(out_dir, indicator, '_upper', measure, '_counts_raster.tif'), overwrite = TRUE)
  writeRaster(lower_counts_raster, paste0(out_dir, indicator, '_lower', measure, '_counts_raster.tif'), overwrite = TRUE)
  
}


## Save aggregates -------------------------------------------------------------------------

# loop over admins
for (i in 0:2) {
  
  # read in aggregates and save
  ad <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '_admin_', i, rake, measure, '_summary.csv'))
  ad[, V1 := NULL]
  write.csv(ad, paste0(out_dir, indicator, measure, '_estimate_ad', i, '.csv'))
    
  # format for figures and save
  if (i != 0) {
    for (m in c('mean', 'upper', 'lower')) {
      ad2 <- ad[, grep(paste0('ADM', i, '_CODE|year|', m), names(ad)), with = FALSE]
      setnames(ad2, m, 'value')
      write.csv(ad2, paste0(out_dir2, indicator,
                            ifelse(measure == '_deaths', '_mortality', measure),
                            '_', m, rake, '_ad', i, '.csv'))
    }
  }
  
  # get counts of children with diarrhea
  ad[, mean := round(mean*pop)]
  ad[, upper := round(upper*pop)]
  ad[, lower := round(lower*pop)]
  ad[, cirange := NULL]
  
  # save
  write.csv(ad, paste0(out_dir, indicator, measure, '_counts_ad', i, '.csv'))
  
  # format for figures and save
  if (i != 0) {
    ad2 <- ad[, grep(paste0('ADM', i, '_CODE|year|mean'), names(ad)), with = FALSE]
    setnames(ad2, 'mean', 'value')
    write.csv(ad2, paste0(out_dir2, indicator,
                          ifelse(measure == '_deaths', '_mortality', measure),
                          '_mean_counts', rake, '_ad', i, '.csv'))
  }
  
}
  
  
## Save absolute and relative difference -------------------------------------------------------------------------

# read in data and merge
ad0 <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '_admin_0', rake, measure, '_summary.csv'))
ad0[, V1 := NULL]
ad2 <- fread(paste0('<<<< FILEPATH REDACTED >>>>/', indicator, '_admin_2', rake, measure, '_summary.csv'))
ad2[, V1 := NULL]
ad <- merge(ad0, ad2, by = c('ADM0_CODE', 'ADM0_NAME', 'region', 'year'))

# calculate relative difference and save
ad_rel <- copy(ad)
ad_rel[, value := (mean.y - mean.x)/mean.x]
write.csv(ad_rel[, grep(paste0('ADM|year|value'), names(ad_rel)), with = FALSE],
          paste0(out_dir, indicator, measure, '_rel_dev_mean_ad2.csv'))
write.csv(ad_rel[, grep(paste0('ADM2_CODE|year|value'), names(ad_rel)), with = FALSE],
          paste0(out_dir2, indicator, ifelse(measure == '_deaths', '_mortality', measure), '_rel_dev_mean', rake, '_ad2.csv'))

# calculate absolute difference and save
ad_abs <- copy(ad)
ad_abs[, value := mean.y - mean.x]
write.csv(ad_abs[, grep(paste0('ADM|year|value'), names(ad_abs)), with = FALSE],
          paste0(out_dir, indicator, measure, '_abs_dev_mean_ad2.csv'))
write.csv(ad_abs[, grep(paste0('ADM2_CODE|year|value'), names(ad_abs)), with = FALSE],
          paste0(out_dir2, indicator, ifelse(measure == '_deaths', '_mortality', measure), '_abs_dev_mean', rake, '_ad2.csv'))

