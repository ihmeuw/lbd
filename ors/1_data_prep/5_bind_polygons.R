####################################################################################################
## Clean collapsed and resample data for for oral rehydration therapy (ORT) modeling
####################################################################################################

## Setup -------------------------------------------------------------------------
rm(list=ls())

# list all indicators
indics <- c('any_ors', 'rhf_only', 'no_ort')

# get input version from most recently collapsed data
files <- file.info(list.files(paste0('<<<< FILEPATH REDACTED >>>>'), pattern = '*.RData', full.names=TRUE))
files <- rownames(files)
files <- gsub('.RData', '', files)
files <- gsub('<<<< FILEPATH REDACTED >>>>/alldata_', '', files)
input_version <- max(files)

# load packages
library(data.table)

# load collapsed data and subset by points
load(paste0('<<<< FILEPATH REDACTED >>>>/countdata_', input_version, '.RData'))

for (i in indics) {

  ## Load and bind resampled polygon data ---------------------------------------------------------------------------
  message(i)
  
  # bind polygon-resampled shapefiles back together
  filenames <- list.files(paste0('<<<< FILEPATH REDACTED >>>>', i), pattern = '*.csv', full.names = TRUE)
  shps <- lapply(filenames, fread)
  all_poly_data <- rbindlist(shps, use.names = TRUE)
  
  # clean up resampled polygon data
  all_poly_data <- all_poly_data[, c('latitude', 'longitude', 'V1') := NULL]
  setnames(all_poly_data, c('lat', 'long'), c('latitude', 'longitude'))
  
  # load point data
  data <- countdata[[i]]
  point_data <- data[point == 1]
  
  # bind with polygon data
  point_data[, weight := 1] # set weight for points to 1
  data <- rbind(all_poly_data, point_data, use.names = TRUE, fill = FALSE)
  
  
  ## Format data & save --------------------------------------------------------------------------------
  
  # format and create cluster ID
  data <- data[, list(nid, country, source, year, latitude, longitude, N, N_obs, mean, weight, sum_of_sample_weights, point)]
  data <- data[, cluster_id := .GRP, keyby = 'nid,country,source,year,latitude,longitude,N,mean,weight,sum_of_sample_weights,point']
  
  # create indicator variable with total count
  data <- data[, (i) := mean]
  
  # convert to proportion (rate) and check that it's between 0 and 1
  data <- data[, rate := mean/N]
  data <- data[!is.na(rate)] # need this check here until we fix the polygon matching issues
  if (min(data$rate < 0) | max(data$rate > 1)) stop('Proportion is not between 0 and 1')
  
  # save model run input data
  all_collapsed <- copy(data)
  write.csv(all_collapsed, file = paste0('<<<< FILEPATH REDACTED >>>>/', i, '.csv'), row.names=F)
  saveRDS(all_collapsed, file = paste0('<<<< FILEPATH REDACTED >>>>', i, '/', i, '_', Sys.Date(), '.rds'))
}

