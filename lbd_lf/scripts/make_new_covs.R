########## Produce new covariates (MDA, etc.) ##########

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

library(data.table)
library(fasterize)
library(matrixStats)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

files <- list.files(paste0(<<<< FILEPATH REDACTED >>>>), full.names = T)
files <- grep(".tif$", files, value = T)
br <- lapply(files, raster)
br <- brick(br)

files <- list.files(paste0(<<<< FILEPATH REDACTED >>>>), full.names = T)
files <- grep(".tif$", files, value = T)
lf_br <- lapply(files, raster)
lf_br <- brick(lf_br)

files <- list.files(paste0(<<<< FILEPATH REDACTED >>>>), full.names = T)
files <- grep(".tif$", files, value = T)
oncho_br <- lapply(files, raster)
oncho_br <- brick(oncho_br)

### Make combined LF & Oncho covariate #####################################
#sum each year LF + Oncho covariate, if result > 0, then assign boolean 1.
#NOTE: Not the same num of layers, line up years by name
lf_years <- lapply(names(lf_br), FUN = function(n) unlist(strsplit(n, "_"))[4]) %>% unlist %>% as.integer
oncho_years <- lapply(names(oncho_br), FUN = function(n) unlist(strsplit(n, "_"))[4]) %>% unlist %>% as.integer

all_mda_years <- min(c(lf_years, oncho_years)):max(c(lf_years, oncho_years))
all_mda_names <- paste0("allmda_binary_1y_", all_mda_years, "_00_00")

#map values to replace in combined raster
reclass <- data.table(from = c(0, 1), to = c(0, 16), becomes = c(0, 1))
reclass <- matrix(as.numeric(unlist(reclass)), nrow = nrow(reclass))

all_mda <- lapply(all_mda_years, FUN = function(y){
  if(!(y %in% lf_years)){
    b <- oncho_br[[which(oncho_years == y)]]
  }else{
    b <- lf_br[[which(lf_years == y)]] + oncho_br[[which(oncho_years == y)]]
    b <- reclassify(b, rcl = reclass, include.lowest = T)
  }
  return(b)
})

output_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
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
for(y in 2:length(all_mda_years)){
  allmda_numrounds[[y]] <- do.call("sum", all_mda[1:y])
}

output_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
for(y in 1:length(all_mda_years)){
  writeRaster(
    allmda_numrounds[[y]],
    file = paste0(<<<< FILEPATH REDACTED >>>>),
    format = "GTiff",
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999
  )
}

### Make modified mdavector.numrounds covariate that respects endemicity #####################################
ad2_shape <- fast_load_shapefile(shape = "g2015_2014_2")

ad2_shape@data$id <- c(1:nrow(ad2_shape@data))

check_mda <- function(adm2, shape = ad2_shape, ref = br[[nlayers(br)]]) {
  poly <- shape[shape@data$id == adm2, ]
  review <- crop_set_mask(ref, poly)
  rev_values <- unique(values(review))
  rev_values <- rev_values[is.na(rev_values) == F]
  if (length(rev_values) == 1) {
    if (rev_values == 0) {
      return(0)
    }
    if (rev_values > 0) {
      return(1)
    }
  } else {
    return(1)
  }
}

adm2_l <- ad2_shape@data$id
results <- lapply(adm2_l, check_mda) %>% unlist()
ad2_shape@data$mda <- results

df <- fread(<<<< FILEPATH REDACTED >>>>)
df <- df[point == 1]
pts <- df[, .(latitude, longitude)]
pts$latitude <- as.numeric(pts$latitude)
pts$longitude <- as.numeric(pts$longitude)
coordinates(pts) <- ~ longitude + latitude
proj4string(pts) <- proj4string(ad2_shape)
adm_data <- sp::over(pts, ad2_shape) %>%
  as.data.table() %>%
  .[, .(GAUL_CODE, ADM2_NAME, id)]
data <- cbind(df, adm_data)

check_data <- function(adm2, shape = ad2_shape, ref = data) {
  subset <- ref[id == adm2, ]
  if (nrow(subset[as.numeric(lf_prev) > 0.01, ]) > 0) {
    return(1)
  } else {
    return(0)
  }
}

ad2_shape@data$prev_01 <- lapply(adm2_l, check_data) %>% unlist()
ad2_data <- ad2_shape@data %>% as.data.table()
ad2_data[mda == 0 & prev_01 == 0, nonintervention := 1]
ad2_data[is.na(nonintervention), nonintervention := 0]
ad2_shape@data$nonintervention <- ad2_data$nonintervention

ad2_raster <- rasterize(ad2_shape, br, "nonintervention")
numrounds.new <- raster::mask(br, ad2_raster, maskvalue = 1, updatevalue = -1)

output_dir <- <<<< FILEPATH REDACTED >>>>
years <- c(2000:2016)
for (i in 1:nlayers(numrounds.new)) {
  # mda_binary[[i]] <- stackApply(br_pre_post[[1:i]], indices = rep(1, i), fun = "max")
  #
  writeRaster(
    numrounds.new[[i]],
    file = paste0(<<<< FILEPATH REDACTED >>>>),
    format = "GTiff",
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999
  )
}

### Make lagged covariates #######################################################################

itn_files <- list.files(paste0(<<<< FILEPATH REDACTED >>>>), full.names = T)
itn_files <- grep(".tif$", itn_files, value = T)
br <- lapply(itn_files, raster)

output_dir <- <<<< FILEPATH REDACTED >>>>

years <- c(2000:2016)

br_new <- append(br, br[[1]], 1)
br_new <- brick(br_new)

for (i in 1:nlayers(br_new)) {
  writeRaster(
    br_new[[i]],
    file = paste0(<<<< FILEPATH REDACTED >>>>),
    format = "GTiff",
    overwrite = TRUE,
    datatype = "FLT4S",
    NAflag = -9999
  )
}
