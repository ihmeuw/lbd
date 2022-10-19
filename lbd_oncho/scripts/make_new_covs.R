
user <- Sys.info()[["user"]]

## Set repo locations
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

path <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

library(data.table)
library(fasterize)
library(matrixStats)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

files <- list.files(<<<< FILEPATH REDACTED >>>>, full.names = T)
files <- grep(".tif$", files, value = T)
lf_br <- lapply(files, raster)
lf_br <- raster::brick(lf_br)

files <- list.files(<<<< FILEPATH REDACTED >>>>, full.names = T)
files <- grep(".tif$", files, value = T)
oncho_br <- lapply(files, raster)
oncho_br <- raster::brick(oncho_br)

### Make combined LF & Oncho covariate #####################################
#sum each year LF + Oncho covariate, if result > 0, then assign boolean 1.
#NOTE: Not the same num of layers, line up years by name
lf_years <- lapply(names(lf_br), FUN = function(n) unlist(strsplit(n, "_"))[4]) %>% unlist %>% as.integer
oncho_years <- lapply(names(oncho_br), FUN = function(n) unlist(strsplit(n, "_"))[4]) %>% unlist %>% as.integer

lf_mda_names <- paste0("lfmda_binary_1y_", lf_years, "_00_00")

all_mda_years <- min(c(lf_years, oncho_years)):max(c(lf_years, oncho_years))
all_mda_names <- paste0("allmda_binary_1y_", all_mda_years, "_00_00")

#map values to replace in combined raster
reclass <- data.table(from = c(0, 1), to = c(0, 33), becomes = c(0, 1))
reclass <- matrix(as.numeric(unlist(reclass)), nrow = nrow(reclass))

### Fix missing pixels (NAs that should be zeroes)
# Use first layer of oncho_br as template, since it has zeroes in all land pixels
template <- copy(oncho_br[[1]])
lf_br_temp <- copy(lf_br)

for (i in 1:raster::nlayers(lf_br_temp)) {
  lf_br_temp[[i]] <- raster::mask(sum(lf_br_temp[[i]], template, na.rm = TRUE), template)
}

lf_br <- copy(lf_br_temp)

output_dir <- <<<< FILEPATH REDACTED >>>>
for(y in 1:length(lf_years)){
  writeRaster(
    lf_br[[y]],
    file = paste0(output_dir, "/", lf_mda_names[y]),
    format = "GTiff",
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999
  )
}

all_mda <- lapply(all_mda_years, FUN = function(y){
  print(y)
  if (!(y %in% lf_years)) {
    print("no lf")
    b <- oncho_br[[which(oncho_years == y)]]
  } else if (y == 2019 & (!(y %in% oncho_years))) {
    print("no oncho")
    b <- lf_br[[which(lf_years == y)]] + oncho_br[[which(oncho_years == (y - 1))]] # special case for 2019
    b <- raster::reclassify(b, rcl = reclass, include.lowest = T)
  } else {
    print("both present")
    b <- lf_br[[which(lf_years == y)]] + oncho_br[[which(oncho_years == y)]]
    b <- raster::reclassify(b, rcl = reclass, include.lowest = T)
  }
  return(b)
})

output_dir <- <<<< FILEPATH REDACTED >>>>
for(y in 1:length(all_mda_years)){
  writeRaster(
    all_mda[[y]],
    file = paste0(output_dir, "/", all_mda_names[y]),
    format = "GTiff",
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999
  )
}

allmda_numrounds <- vector("list", length(all_mda_years))
allmda_numrounds[[1]] <- all_mda[[1]]
for (y in 2:length(all_mda_years)) {
  print(y)
  allmda_numrounds[[y]] <- all_mda[[y]] + allmda_numrounds[[y - 1]]
}

output_dir <- <<<< FILEPATH REDACTED >>>>
for (y in 1:length(all_mda_years)) {
  raster::writeRaster(
    allmda_numrounds[[y]],
    file = <<<< FILEPATH REDACTED >>>>,
    format = "GTiff",
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999
  )
}

########## Reprocess LF numrounds rasters to replace NAs with 0 values
files <- list.files(<<<< FILEPATH REDACTED >>>>, full.names = T)
files <- grep(".tif$", files, value = T)
lf_br <- lapply(files, raster)
lf_br <- raster::brick(lf_br)

lf_br_original <- copy(lf_br)

for (i in 1:raster::nlayers(lf_br)) {
  print(i)
  lf_br[[i]][lf_br[[i]] == 128] <- NA
  lf_br[[i]] <- raster::mask(sum(lf_br[[i]], template, na.rm = TRUE), template)
}

lf_mda_years <- 2000:2019

output_dir <- <<<< FILEPATH REDACTED >>>>
for (y in 1:length(lf_mda_years)) {
  raster::writeRaster(
    lf_br[[y]],
    file = <<<< FILEPATH REDACTED >>>>,
    format = "GTiff",
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999
  )
}
