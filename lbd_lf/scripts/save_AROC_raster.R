library(raster)
library(seegMBG)

run_date <- <<<< FILEPATH REDACTED >>>>
region <- <<<< FILEPATH REDACTED >>>>

cells <- readRDS(paste0(<<<< FILEPATH REDACTED >>>>))
load(paste0(<<<< FILEPATH REDACTED >>>>))
mean_raster  <- insertRaster(simple_raster, matrix(rowMeans(cells), ncol=1))
writeRaster(mean_raster, file=paste0(<<<< FILEPATH REDACTED >>>>), format='GTiff', overwrite=TRUE)
plot(mean_raster)

### Make mask
mean_prev_2000 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
mean_prev_2000[mean_prev_2000 < 0.01] <- 0
mean_prev_2000[mean_prev_2000 > 0] <- 1
writeRaster(mean_prev_2000, file=paste0(<<<< FILEPATH REDACTED >>>>), format="GTiff", overwrite=TRUE)
