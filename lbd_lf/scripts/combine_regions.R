##############################################################################################################
# Combine rasters and CSVs for neighboring regions (to create composite outputs)
##############################################################################################################

user <- Sys.info()[["user"]]

## Set repo locations
core_repo <- paste0(<<<< FILEPATH REDACTED >>>>)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

path <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- paste0(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(matrixStats)
library(sf)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

run_date1 <- <<<< FILEPATH REDACTED >>>>
run_date2 <- <<<< FILEPATH REDACTED >>>>
region1 <- "lf_s_asia"
region2 <- "lf_se_asia"

years <- c(2000, 2005, 2010, 2018)

### Create composite folder in MBG maps folder
dir.create(paste0(<<<< FILEPATH REDACTED >>>>), showWarnings = TRUE)
dir.create(paste0(<<<< FILEPATH REDACTED >>>>), showWarnings = TRUE)

out_path <- paste0(<<<< FILEPATH REDACTED >>>>)

### Combine prevalence rasters
# measures <- c("mean", "upper", "lower", "range")
measures <- c("mean")
for (j in measures) {
  print(j)
  for (i in 1:length(years)) {
    print(years[i])
    raster1 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
    raster2 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
    merged <- merge(raster1, raster2)
    writeRaster(merged, filename=paste0(<<<< FILEPATH REDACTED >>>>), format='GTiff', overwrite=TRUE)
  }
}


### Combine admin summary CSVs
admin <- c(0, 1, 2)
for (i in 1:length(admin)) {
  csv1 <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  csv2 <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
  setnames(csv1, "value", "mean", skip_absent = TRUE)
  setnames(csv2, "value", "mean", skip_absent = TRUE)
  appended <- rbind(csv1, csv2)
  write.csv(appended, paste0(<<<< FILEPATH REDACTED >>>>))
}


### Combine mean_prev_0.01_binary_2000.tif rasters
binary_raster1 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
binary_raster2 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
binary_merged <- merge(binary_raster1, binary_raster2)
writeRaster(binary_merged, filename=paste0(<<<< FILEPATH REDACTED >>>>), format='GTiff', overwrite=TRUE)


### Combine pred derivatives
dir.create(paste0(<<<< FILEPATH REDACTED >>>>), showWarnings = TRUE)
dir.create(paste0(<<<< FILEPATH REDACTED >>>>), showWarnings = TRUE)
dir.create(paste0(<<<< FILEPATH REDACTED >>>>), showWarnings = TRUE)

aroc_raster1 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
aroc_raster2 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
aroc_merged <- merge(aroc_raster1, aroc_raster2)
writeRaster(aroc_merged, filename=paste0(<<<< FILEPATH REDACTED >>>>), format='GTiff', overwrite=TRUE)

targets <- c(0.01, 0.02)
for (t in targets) {
  print(t)
  for (i in 1:length(admin)) {
    print(admin[i])
    target_raster1 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
    target_raster2 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
    target_merged <- merge(target_raster1, target_raster2)
    writeRaster(target_merged, filename=paste0(<<<< FILEPATH REDACTED >>>>), format='GTiff', overwrite=TRUE)
    
    target_csv1 <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
    target_csv2 <- fread(paste0(<<<< FILEPATH REDACTED >>>>))
    setnames(target_csv1, "value", "mean", skip_absent = TRUE)
    setnames(target_csv2, "value", "mean", skip_absent = TRUE)
    target_appended <- rbind(target_csv1, target_csv2)
    write.csv(target_appended, paste0(<<<< FILEPATH REDACTED >>>>))
  }
  target_raster1 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
  target_raster2 <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
  target_merged <- merge(target_raster1, target_raster2)
  writeRaster(target_merged, filename=paste0(<<<< FILEPATH REDACTED >>>>), format='GTiff', overwrite=TRUE)
}
