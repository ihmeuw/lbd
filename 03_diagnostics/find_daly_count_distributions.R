library(data.table)
library(dplyr)
library(raster)

rm(list =ls())
admin_level <- 0

dia <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(measure == 'dalys') %>%
  select(mean_count, paste0('ADM', admin_level, '_CODE'), year) %>%
  rename(dia_count = mean_count)
mal <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(measure == 'dalys') %>%
  select(mean_count, paste0('ADM', admin_level, '_CODE'), year) %>%
  rename(mal_count = mean_count)
lri <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(measure == 'dalys') %>%
  select(mean_count, paste0('ADM', admin_level, '_CODE'), year) %>%
  rename(lri_count = mean_count)

total <- merge(dia, lri, by = c(paste0('ADM', admin_level, '_CODE'), 'year')) %>%
  merge(mal, by = c(paste0('ADM', admin_level, '_CODE'), 'year')) %>%
  as.data.table()

summary(total$dia_count)
summary(total$lri_count)
summary(total$mal_count)

total <- total[,sum := dia_count + lri_count + mal_count]

summary(total$sum)


#get max by pixel
dia_raster_00 <- raster('<<< FILEPATH REDACTED >>>')
dia_raster_17 <- raster('<<< FILEPATH REDACTED >>>')

mal_raster_00 <- raster('<<< FILEPATH REDACTED >>>')
mal_raster_17 <- raster('<<< FILEPATH REDACTED >>>')

lri_raster_00 <- raster('<<< FILEPATH REDACTED >>>')
lri_raster_17 <- raster('<<< FILEPATH REDACTED >>>')

total_raster_00 <- raster('<<< FILEPATH REDACTED >>>')
total_raster_17 <- raster('<<< FILEPATH REDACTED >>>')

max(values(dia_raster_00), na.rm = TRUE)
max(values(dia_raster_17), na.rm = TRUE)

max(values(mal_raster_00), na.rm = TRUE)
max(values(mal_raster_17), na.rm = TRUE)

max(values(lri_raster_00), na.rm = TRUE)
max(values(lri_raster_17), na.rm = TRUE)

max(values(total_raster_00), na.rm = TRUE)
max(values(total_raster_17), na.rm = TRUE)
