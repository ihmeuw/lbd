#code to make a usable file out of oxford tifs
rm(list=ls())

#set variables
measure <- 'incidence'
shapefile_version <- '2019_09_10'
run_date <- '2019_10_28'
extraction_date <- '2019_08_14'
raked <- FALSE
pop_release <- '2019_08_29'

#other variables
modeling_shapefile_version <- shapefile_version

#load libraries
library(seegMBG)
library(seegSDM)
library(doParallel)
library(data.table)
library(rgdal)

#load lbd core functions
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/prep_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/shapefile_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/covariate_functions.R')

#create output directory
dir.create('<<< FILEPATH REDACTED >>>')

#get list of relevant files
#incidence files do not have a sex split, while mortality files do
#mortality files exist in rate/count, unraked/rake form
if (measure == 'mortality'){
  folder_in <- '<<< FILEPATH REDACTED >>>'
  if (raked) {
    files <- list.files(folder_in, pattern="raked.death.rate.infant", full.names=T)
  } else {
    files <- list.files(folder_in, pattern="[0-9].death.rate.infant", full.names=T)
  }
} else if (measure == 'incidence'){
  folder_in <- '<<< FILEPATH REDACTED >>>'
  files <- list.files(folder_in, pattern="infants.full.tif$", full.names=T)
}


#check to make sure all files present: expect 100 per sex/year
for (year in c(2000:2016)){
  test <-Filter(function(x) grepl(year, x), files)
  if (measure == 'mortality'){
    if (length(test) != 200){
      message('STOP! Year ', year, ' is missing ', (200 - length(test)), ' files!')
    } else {
      message(paste0('all tif files are present for year ', year))
    }
  } else if (measure == 'incidence'){
    if (length(test) != 100){
      message('STOP! Year ', year, ' is missing ', (100 - length(test)), ' rows!')
    } else {
      message(paste0('all tif files are present for year ', year))
    }
  }
}
rm(test)

gaul_list           <- get_adm0_codes('africa-ESH', shapefile_version = shapefile_version)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 0.4)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

raster_list    <- build_simple_raster_pop(subset_shape,
                                          pop_release = pop_release)
simple_raster <- raster_list[['simple_raster']]


tif_fix <- function(malaria){
  malaria  <- extend(malaria, simple_raster)
  malaria  <- crop(malaria, simple_raster)
  malaria[is.na(malaria)] = 0
  malaria  <- mask(malaria, simple_raster)
  malaria <- as.numeric(as.matrix(malaria[cellIdx(malaria)]))
  return(malaria)
}


if (measure == 'mortality'){
  #make csv pixel years long and draws wide
  preproc_data <- data.table()
  for (year in c(2000:2016)) {
    message(paste0('Starting year: ', year))
    files_by_year <- Filter(function(x) grepl(year, x), files)
    cl <- makeCluster(40)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    registerDoParallel(cl)
    #Read in each .dta file in parallel - returns a list of data frames
    top <- foreach(i=1:(length(files_by_year)/2), .packages=c('raster', 'seegMBG', 'seegSDM')) %dopar% {
      tif <- rowMeans(cbind(as.data.frame(tif_fix(raster(files_by_year[2*i - 1]))),as.data.frame(tif_fix(raster(files_by_year[2*i]))))) #deal with the sex split
    }
    stopCluster(cl)


    message("cbindlist all draws together")
    topics <- as.data.table(do.call(cbind, top))
    rm(top)


    preproc_data <- rbind(preproc_data, topics)
  }

  preproc_data <- as.matrix(preproc_data)
  if (raked){
    saveRDS(preproc_data, '<<< FILEPATH REDACTED >>>')
  } else {
    saveRDS(preproc_data, '<<< FILEPATH REDACTED >>>')
  }
}else if (measure == 'incidence'){
  message(measure)
  
  #make csv pixel years long and draws wide
  preproc_data <- data.table()
  for (year in c(2000:2016)) {
    message(paste0('Starting year: ', year))
    files_by_year <- Filter(function(x) grepl(year, x), files)
    cl <- makeCluster(40)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    registerDoParallel(cl)
    
    #Read in each .dta file in parallel - returns a list of data frames
    top <- foreach(i=1:(length(files_by_year)), .packages=c('raster', 'seegMBG', 'seegSDM')) %dopar% {
      tif <- tif_fix(raster(files_by_year[i]))
    }
    stopCluster(cl)


    message("cbindlist all draws together")
    topics <- as.data.table(do.call(cbind, top))
    rm(top)


    preproc_data <- rbind(preproc_data, topics)
  }

  preproc_data <- as.matrix(preproc_data)
  saveRDS(preproc_data, '<<< FILEPATH REDACTED >>>')
}

