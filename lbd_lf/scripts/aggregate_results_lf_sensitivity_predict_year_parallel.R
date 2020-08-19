# Make Admin 1 and 2 level summary draw-level objects

############# Setup

source(<<<< FILEPATH REDACTED >>>>)
load_from_parallelize()
message(paste0("Using ", core_repo))

user <- Sys.info()[["user"]]

## Set repo locations
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

path <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- paste0(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(sf)
library(matrixStats)

setwd(indic_repo)

## Focal 3 specific workflow: Pull in custom scripts.
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

message(paste0("threshold_index: ", threshold_index))

message(paste0("year_predict: ", year_predict))

thresholds <- list(c(0, 10000000, "_prev_0_pop_inf"), c(0, 750000, "_prev_0_pop_750k"), c(0, 250000, "_prev_0_pop_250k"), c(0, 100000, "_prev_0_pop_100k"), c(0, 50000, "_prev_0_pop_50k"), c(0, 25000, "_prev_0_pop_25k"), c(0, 10000, "_prev_0_pop_10k"), c(0, 5000, "_prev_0_pop_5k"))

message(thresholds[[threshold_index]])
pop_threshold <- as.integer(thresholds[[threshold_index]][[2]])
prev_threshold <- as.numeric(thresholds[[threshold_index]][[1]])
file_suffix <- thresholds[[threshold_index]][[3]]

message(paste0("pop_threshold: ", pop_threshold))
message(paste0("prev_threshold: ", prev_threshold))
message(paste0("file_suffix: ", file_suffix))

config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator, 
                                config_file = paste0(<<<< FILEPATH REDACTED >>>>), run_tests = FALSE)

message(paste0("\nRun date: ", run_date, "\n"))

## make a pathaddin that get used widely
pathaddin <- paste0(<<<< FILEPATH REDACTED >>>>)

reg <- region # convenience

share_dir <- paste0(<<<< FILEPATH REDACTED >>>>)

if (class(year_list) == "character") year_list <- eval(parse(text = year_list))
if (!exists("predict_years")) {
  predict_years <- year_list
}
if (class(predict_years) == "character") predict_years <- eval(parse(text = predict_years))

predict_years_full <- copy(predict_years)
predict_years <- year_predict

message(paste0("predict_years: ", predict_years))

if (exists("use_median_aggregation")) {
  use_median_aggregation <- as.logical(use_median_aggregation)
} else {
  use_median_aggregation <- FALSE
}

# Determining whether this is a re-run, stopping if overwrite is off and the file already exists.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (overwrite == T) {
  if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>)))) {
    message("This file already exists, this must be a re-run.")
  }
} else {
  if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>)))) {
    stop("Looks like you have already prepped that admin draw file; skipping this!")
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Looping through Regions, aggregating data to Admin Units
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start <- proc.time()

message(paste0("Pulling numbers for indicator: ", indicator, " region: ", reg, " rundate: ", run_date, " using the population measure ", pop_measure, " for weighting."))
pathaddin <- paste0(<<<< FILEPATH REDACTED >>>>)

# Loading Data
# ~~~~~~~~~~~~~~
message("Loading in the draw-level data for this region; this may take a while (but probably not longer than 10 minutes).")
# Load in cell-level results based on that same model run
if (raked) {
  # Different file conventions in use; try either rds or rdata until standardized
  filename_rds <- paste0(<<<< FILEPATH REDACTED >>>>)
  filename_rdata <- paste0(<<<< FILEPATH REDACTED >>>>)
  if (file.exists(filename_rds)) cell_pred <- readRDS(filename_rds)
  if (file.exists(filename_rdata)) {
    load(filename_rdata, verbose = T)
    cell_pred <- raked_cell_pred
    rm(raked_cell_pred)
  }
} else {
  message(paste0(<<<< FILEPATH REDACTED >>>>))
  load(paste0(<<<< FILEPATH REDACTED >>>>))
}

# Check if load was successful; stop if not
if (!exists("cell_pred")) stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")

# Getting the simple polygon and simple raster objects for this region alone
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Getting the spatial objects associated with this region.")
if (!(file.exists(paste0(<<<< FILEPATH REDACTED >>>>)))) {
  gaul_list <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = use_premade, shapefile_version = modeling_shapefile_version)
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  message("Building simple raster from subset_shape")
  
  simple_raster <- raster_list[["simple_raster"]]
  pop_raster <- raster_list[["pop_raster"]]
} else {
  load(paste0(<<<< FILEPATH REDACTED >>>>))
}

message("All done loading spatial template files (subset_shape, simple_polygon, simple_raster, pop_raster, etc.)")

## Determining a list of the valid pixel indices based on the simple raster template
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pixel_id <- seegSDM:::notMissingIdx(simple_raster)
pixel_spatial <- data.table(pixel_id = pixel_id)
message("Pixel ID that links cell_pred observations to raster locations created from simple raster, checking for sameness in dimensions compared to cell_pred.")
message(paste0("length(pixel_id): ", length(pixel_id)))
message(paste0("nrow(cell_pred): ", nrow(cell_pred)))
message(paste0("length(predict_years): ", length(predict_years)))
message(paste0("nrow(cell_pred) / length(predict_years): ", nrow(cell_pred) / length(predict_years)))
if (!(length(pixel_id) == nrow(cell_pred) / length(predict_years))) {
  stop("Excuse me, but the number of valid pixels in your simple raster is not the same number of valid pixels in your cell_pred object. Look into this!")
} else {
  message(paste0(
    "Check passed: There are ", length(pixel_id), " pixels in your simple_raster object, the same number of pixels in the cell_pred object for each year."
  ))
}

message("Generating a raster of the pixel_id values")
pixel_id_raster <- copy(simple_raster)
for (pix in pixel_id){ # For each of the valid pixels, replace the value with the pixel_id, not whatever country ID is used in the simple_polygon/raster
  pixel_id_raster@data@values[pix] <- pix
}

## Loading and Assigning Admin Units to Pixels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Getting the spatial location (admin unit) of each of the pixel locations.")

region_adm0_list <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version) # Getting the adm0 GAUL codes, we can use this to make sure we don't accidentally include countries from buffers that shouldn't be in this region

## Defining a function that will get the raster versions of each Admin level:
GetAdmin <- function(admin_level, simple_raster, region_adm0_list, shapefile_version, raster_agg = raster_agg_factor){
  message(paste0("Loading admin level ", admin_level))
  
  # load admin shape file
  admin_shp <- rgdal::readOGR(dsn=get_admin_shapefile(admin_level, version = shapefile_version))
  
  # ensure that the rasterize variable is a numeric
  admin_shp@data[[paste0('ADM', admin_level, '_CODE')]] <- as.numeric(as.character(admin_shp@data[[paste0('ADM', admin_level, '_CODE')]]))
  
  # if it doesn't exist, get areas of polygons. 
  if(is.null(admin_shp$Shape_Area)){
    admin_shp$Shape_Area <- raster::area(admin_shp) / 1e6 ## TODO
  }
  
  message("Rasterizing with the custom function...")
  
  admin_rast <- rasterize_check_coverage(admin_shp, simple_raster, paste0("ADM",admin_level,"_CODE"))
  
  message("Converted to raster based on simple_raster template. Cropping and masking:")
  
  admin_rast  <- crop(admin_rast, extent(simple_raster))
  admin_rast  <- setExtent(admin_rast, simple_raster)
  admin_rast <- mask(admin_rast, simple_raster)
  
  message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
  admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% region_adm0_list,]
  admin_centroids<-SpatialPointsDataFrame(gCentroid(admin_shp, byid=TRUE), admin_shp@data, match.ID=FALSE)
  
  message("Compiling and returning results.")
  admin <- list()
  admin[["spdf"]] <- admin_shp
  admin[["centroids"]] <- admin_centroids
  admin[["rast"]] <- admin_rast
  admin[["attributes"]] <- copy(data.table(admin_shp@data))
  
  return(admin)
}

message("Rasterizing shapefiles; this may take a while.")
message(paste0("Raster_agg_factor set to ", raster_agg_factor))
admin_levels <- list() # Empty list of levels that will be filled with admin levels
for (lvl in c(0, 1, 2)) {
  fieldname <- paste0("ADM",lvl,"_CODE")
  admin_info <- GetAdmin(admin_level = lvl, simple_raster, region_adm0_list, shapefile_version = modeling_shapefile_version)
  pixel_spatial[[fieldname]] <- extract(admin_info[["rast"]], pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
  admin_levels[[as.character(lvl)]] <- admin_info # Add the admin info to the list
  if (sum(is.na(pixel_spatial[[fieldname]])) > 0){ # Check to see if any of the pixels don't have a location assigned
    message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level ",lvl))
  }
}


## Loading child model predictions (code modified from Dan)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Loading child stacker model predictions")

## Get model information from RData objects
fetch_from_rdata <- function(file_location, item_name, use_grep = F) {
  load(file_location)
  
  if (use_grep) {
    ret_obj <- lapply(item_name, function(x) mget(grep(x, ls(), value = T)))
  } else {
    ret_obj <- lapply(item_name, function(x) get(x))
  }
  
  if (length(item_name) == 1) {
    ret_obj <- ret_obj[[1]]
  }
  return(ret_obj)
}

## Add stacking children to this process
covs <- fetch_from_rdata(<<<< FILEPATH REDACTED >>>>)
fes <- fetch_from_rdata(<<<< FILEPATH REDACTED >>>>)
submodels <- trimws(strsplit(fes, "+", fixed = T)[[1]])
covs <- covs[submodels]

## Make sure spatial extent is the same
covnames <- names(covs)

## Ensure the dimensions are the same
for (ccc in covs) {
  stopifnot(dim(ccc)[1:2] == dim(simple_raster)[1:2])
}

covs_full <- copy(covs)

years <- year_predict

adm0_list <- vector("list", length(years))
adm1_list <- vector("list", length(years))
adm2_list <- vector("list", length(years))

for (i in 1:length(years)) {
  current_years <- years[[i]]
  message(paste("Working on year list:", paste(current_years, collapse = " ")))
  
  start_year <- min(current_years)
  end_year <- max(current_years)
  index_start <- which(predict_years == start_year)
  index_end <- which(predict_years == end_year)
  cell_length <- nrow(cell_pred) / length(predict_years)
  cell_pred_current <- cell_pred[(cell_length * (index_start - 1) + 1):(cell_length * index_end), 1:ncol(cell_pred)]
  
  ## These are where draw-level, admin-level results will be stored
  message("Starting empty lists for results")
  regions <- list()
  admin_0 <- list()
  admin_1 <- list()
  admin_2 <- list()
  sp_hierarchy_list <- list()
  
  covs <- copy(covs_full)
  
  ## Pulling Population
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Pull annual population brick using new covariates function
  message(paste0("Pulling ", pop_measure, " population raster."))
  
  worldpop_current <- fixed_effects_config[covariate == "worldpop",]
  worldpop_current$year_lag <- 0
  
  pop <- load_lagged_covariates(covariate_config = worldpop_current,
                                template = simple_raster,
                                start_year = start_year,
                                end_year = end_year,
                                raster_agg = as.integer(raster_agg_factor))$worldpop
  
  message("Ensuring that the spatial extents of the population matches the simple raster")
  pop <- raster::crop(pop, extent(simple_raster))
  pop <- setExtent(pop, simple_raster)
  if (!(all.equal(dim(pop)[1:2], dim(simple_raster)[1:2]))) pop <- raster::resample(pop, simple_raster, method="ngb")
  pop <- raster::mask(pop, simple_raster)
  
  if (exists("rake_gbd_pop")) {
    message(paste0("rake_gbd_pop: ", rake_gbd_pop))
    if (rake_gbd_pop) { # Prep for rake to GBD national populations
      message("Raking to GBD population...")
      if (Sys.info()[1] == "Windows") root <- <<<< FILEPATH REDACTED >>>> else root <- <<<< FILEPATH REDACTED >>>>
      source(<<<< FILEPATH REDACTED >>>>)
      source(<<<< FILEPATH REDACTED >>>>)
      
      locs <- get_location_metadata(gbd_round_id = 6, location_set_id = 35)
      loc_ids <- locs[, location_id]
      gbd <- get_population(decomp_step = "iterative", gbd_round_id = 6, year_id = start_year:end_year, location_id = loc_ids, run_id=189)
      gbd <- merge(as.data.table(gbd), as.data.table(locs)[, .(location_id, ihme_loc_id)])
      
      connector <- get_gbd_locs(
        rake_subnational = T,
        reg = reg,
        shapefile_version = modeling_shapefile_version
      ) %>% as.data.table()
      
      pop_dt <- data.table(extract(pop, pixel_id)) # getting the pixels from the population raster that correspond to the valid pixels in the pop-raster
      message("Reshaping population from wide to long, and generating an id variable that corresponds to the simple raster.")
      pop_dt[, pixel_id := pixel_id]
      pop_dt <- melt(pop_dt, id.vars = "pixel_id") # Melting the dataframe such that it is long () and should match the number of rows in cell_pred
      pop_dt[, year := (min(current_years) - 1) + as.numeric(gsub("worldpop.", "", variable))] # Converting "worldpop.1" variables to actual years.
      pop_dt <- pop_dt[, list(pixel_id, year, pop = value)] # Subsetting columns
      pop_dt[is.na(pop), pop := 0] # Setting values where pop is NA to 0
      message("Checking whether the population has the same dimensions as cell_pred")
      if (!(nrow(pop_dt) == nrow(cell_pred_current))) {
        stop("Excuse me, but the number of rows in your population table is *not* the same as the number of rows in your cell predictions. Look into this!")
      } else {
        message("Check passed: The population file has the same number of rows as the cell predictions after transformation to long format.")
      }
      
      ## Subnational raking
      lbd_gbd <- fread(<<<< FILEPATH REDACTED >>>>)
      
      pop_dt[, ADM1_CODE := admin_levels$`1`$rast[pixel_id]]
      pop_dt[, ADM0_CODE := simple_raster[pixel_id]]
      
      pop_dt[ADM1_CODE %in% connector$ADM1_CODE, "ADM_CODE" := ADM1_CODE]
      pop_dt <- merge(pop_dt, lbd_gbd, by.x = "ADM1_CODE", by.y = "ADM_CODE", all.x = TRUE)
      
      pop_dt[is.na(ADM_CODE), "ADM_CODE" := ADM0_CODE]
      pop_dt <- merge(pop_dt, lbd_gbd, by.x = "ADM_CODE", by.y = "ADM_CODE", all.x = TRUE)
      
      pop_dt$iso3_1 <- substr(pop_dt$ihme_lc_id.x, 0, 3)
      pop_dt$iso3_2 <- substr(pop_dt$ihme_lc_id.y, 0, 3)
      
      pop_dt <- merge(pop_dt, connector, by = "ADM_CODE", all.x = T)
      
      pixel_ids_to_drop <- pop_dt[(iso3_1 != iso3_2) | (is.na(location_id)), pixel_id]
      pop_dt <- pop_dt[!(pixel_id %in% pixel_ids_to_drop)]
      
      rake_geo_pop <- pop_dt[, lapply(c("pop"), function(x) sum(get(x), na.rm = T)), by = c("year", "location_id")]
      rake_geo_pop$world_pop_unraked <- rake_geo_pop$V1
      loc_ids <- unique(connector$location_id)
      
      scalars <- merge(rake_geo_pop, gbd, by.x = c("location_id", "year"), by.y = c("location_id", "year_id"))
      scalars[, pop_scalar := population / world_pop_unraked]
    }
  }
  
  message("Done generating population rasterbrick; extracting relevant pixels based on simple_raster.")
  pop_dt <- data.table(raster::extract(pop, pixel_id)) # getting the pixels from the population raster that correspond to the valid pixels in the pop-raster
  message("Reshaping population from wide to long, and generating an id variable that corresponds to the simple raster.")
  pop_dt[, pixel_id := pixel_id]
  
  if (exists("pixel_ids_to_drop")) {
    pop_dt[pixel_id %in% pixel_ids_to_drop, worldpop.1 := 0]
  }
  
  pop_dt <- melt(pop_dt, id.vars = "pixel_id") # Melting the dataframe such that it is long () and should match the number of rows in cell_pred
  pop_dt[, year := (min(current_years) - 1) + as.numeric(gsub("worldpop.", "", variable))] # Converting "worldpop.1" variables to actual years.
  pop_dt <- pop_dt[, list(pixel_id, year, pop = value)] # Subsetting columns
  message("Checking whether the population has the same dimensions as cell_pred")
  if (!(nrow(pop_dt) == nrow(cell_pred_current))) {
    stop("Excuse me, but the number of rows in your population table is *not* the same as the number of rows in your cell predictions. Look into this!")
  } else {
    message("Check passed: The population file has the same number of rows as the cell predictions after transformation to long format.")
  }
  
  if (exists("rake_gbd_pop")) {
    if (rake_gbd_pop) { # apply population raking
      message("Applying population raking...")
      
      ## Subnational raking
      lbd_gbd <- fread(<<<< FILEPATH REDACTED >>>>)
      
      pop_dt[, ADM1_CODE := admin_levels$`1`$rast[pixel_id]]
      pop_dt[, ADM0_CODE := simple_raster[pixel_id]]
      
      pop_dt[ADM1_CODE %in% connector$ADM1_CODE, "ADM_CODE" := ADM1_CODE]
      pop_dt <- merge(pop_dt, lbd_gbd, by.x = "ADM1_CODE", by.y = "ADM_CODE", all.x = TRUE)
      
      pop_dt[is.na(ADM_CODE), "ADM_CODE" := ADM0_CODE]
      pop_dt <- merge(pop_dt, lbd_gbd, by.x = "ADM_CODE", by.y = "ADM_CODE", all.x = TRUE)
      
      pop_dt <- merge(pop_dt, connector, by = "ADM_CODE", all.x = T)
      for (s in 1:nrow(scalars)) {
        pop_dt[location_id == scalars[s, location_id] & year == scalars[s, year], pop := pop * scalars[s, pop_scalar]]
      }
    }
  }
  
  pop <- copy(pop_dt)
  pop <- pop[order(year, pixel_id)]
  
  if (as.logical(use_stacking_covs)) {
    ## Subset covariate rasters to prediction years
    index_start_full <- which(predict_years_full == start_year)
    index_end_full <- which(predict_years_full == end_year)
    
    for (k in 1:length(covs)) {
      covs[[k]] <- subset(covs[[k]], index_start_full:index_end_full)
    }
    
    ## Convert to datables, reshape and stuff
    brick_to_dt <- function(bbb) {
      dt <- setDT(as.data.frame(bbb))
      dt[, pxid := .I] # probably uncessary
      
      # #Drop rows now in cellIdx
      dt <- dt[pixel_id, ]
      
      dt <- melt(dt, id.vars = "pxid", variable.factor = F)
      dt <- dt[, .(value)]
      return(dt)
    }
    
    covdt <- lapply(covs, brick_to_dt)
    covdt <- do.call(what = cbind, covdt)
    
    setnames(covdt, names(covs))
    covdt[, pixel_id := rep(pixel_id, nrow(covdt)/length(pixel_id))] # woo recycling
    
    ## Add year to covdt
    yyy <- as.vector(unlist(lapply(min(current_years):max(current_years), function(x) rep.int(x, times = length(pixel_id)))))
    covdt[, year := yyy]
    covdt <- covdt[year %in% current_years]
    
    ## Free up a bit of space
    rm(covs)
  }
  
  ## Merging together cell_pred, stackers, spatial information, and population information
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  pop <- merge(pop, pixel_spatial, by = "pixel_id", all.x = T) # Merging on the spatial information to the population information.
  
  cell_pred_current <- data.table(cell_pred_current) # Converting cell_pred to a data.table
  draw_colnames <- names(cell_pred_current) # names of draw columns
  
  cell_pred_current <- cbind(cell_pred_current, pop) # Adding on the population and spatial information to cell_pred.
  
  if (as.logical(use_stacking_covs)) {
    cell_pred_current <- cbind(cell_pred_current, covdt) # Adding stackers to cell_pred.
    stopifnot(any(!(cell_pred_current[, pixel_id] != rep(pixel_id, nrow(covdt)/length(pixel_id))))) # make sure everything combined properly
    draw_colnames <- c(covnames, draw_colnames) # indicate which columns to aggregate on
  }
  
  # Loading and Assigning Pixels to Admin Units (that are missing)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # We now need to figure out what admin units don't show up when we pull the locations of individual raster pixels of our results.
  # To do this, we need to know what Admin 1 and 2 units exist in the shapefiles we want to use, but don't show up in the results.
  
  missing_admins <- list() # This will contain the cell_preds for otherwise missing admin units, by level.
  for (lvl in c(0, 1, 2)) {
    fieldname<-paste0("ADM",lvl,"_CODE")
    
    # First, we discover what GAUL codes are missing from the pixel_spatial table compared to the shapefile:
    shpfile_attributes<-admin_levels[[as.character(lvl)]][["attributes"]] # Get the attribute table from the admin level's shapefile
    shpfile_attributes = subset(shpfile_attributes,ADM0_CODE != 40762)
    attribute_codes<-shpfile_attributes[[paste0("ADM",lvl,"_CODE")]] # Get list of codes for that level from the admin level's attribute table
    admin_cols<-names(shpfile_attributes)[names(shpfile_attributes)%in%c("ADM0_CODE","ADM1_CODE","ADM2_CODE")] # Any columns in ADM0, ADM1, and ADM2_CODE that should be included as merges.
    pixel_codes<-pixel_spatial[[fieldname]] # Get list of codes based on what's in the pixel_spatial object
    missing_codes<-attribute_codes[!(attribute_codes %in% pixel_spatial[[fieldname]])] # Get list of missing codes
    missing_table<-shpfile_attributes[shpfile_attributes[[fieldname]]%in%missing_codes,admin_cols,with=F]
    
    if(length(missing_codes)==0){
      message(paste0("No missing codes at level ",lvl))
    }else{
      message(paste0("Missing codes found for level ",lvl,":"))
      ## Strategy 1: Assign a pixel location based on centroid location
      ## Develop a raster of the pixel-IDs:
      message("  Discovering centroid location")
      points<-admin_levels[[as.character(lvl)]][["centroids"]] # SpatialPointsDataFrame of that level
      missing_points<-points[(points@data[[fieldname]] %in% missing_codes),] # getting only the missing points
      missing_centroid_locs<-data.table(raster::extract(pixel_id_raster,missing_points)) # Extracting the missing location values as centroid points...
      names(missing_centroid_locs)<-"point_method"
      missing_admins_centroids<-data.table(missing_codes,missing_centroid_locs)
      
      ## Strategy 2: Assign a pixel location based on polygon extraction
      ## Develop a raster of the pixel-IDs:
      message("  Discovering first raster pixel touching polygon")
      polys<-admin_levels[[as.character(lvl)]][["spdf"]] # SpatialPointsDataFrame of that level
      missing_polys<-polys[(polys@data[[fieldname]] %in% missing_codes),] # getting only the missing points
      missing_poly_locs<-data.table(raster::extract(x=pixel_id_raster,y=missing_polys,small=T,fun=function(x,...)first(x))) # Extracting the missing location values as polygons, pulling the first raster cell that it touches...
      names(missing_poly_locs)<-"poly_method"
      missing_admins_polys<-data.table(missing_codes,missing_poly_locs) # Extracting the
      
      ## Merging strategies together: centroids and polygons; adding to list.
      missing_locs<-merge(missing_admins_polys,missing_admins_centroids,by="missing_codes")
      setnames(missing_locs,"missing_codes",fieldname)
      missing_locs<-merge(missing_locs,missing_table,by=fieldname) # Add in admin 0, 1 levels if relevant
      missing_locs[,pixel_id:=point_method] # If centroid produced something, go with centroid
      missing_locs[is.na(point_method),pixel_id:=poly_method] # Otherwise, go with the polygon method.
      
      # For those still missing, assign them to nearest non-NA pixel using centroid
      if(NA %in% missing_locs$pixel_id){
        message('  After centroids and polygon methods, there are still some missing admins.')
        message('  Now sampling to find nearest non-NA pixel and using that')
        
        for(rr in which(is.na(missing_locs$pixel_id))){
          message(sprintf('  -- finding nearest pixel to gaul_code: %i',
                          as.numeric(missing_locs[rr, sprintf('ADM%i_CODE', lvl), with = F])))
          
          ## get centroid of shape
          mp <- polys[polys[[sprintf('ADM%i_CODE', lvl)]] == missing_locs[[sprintf('ADM%i_CODE', lvl)]][rr],]
          cent <- getSpPPolygonsLabptSlots(mp)
          
          ## loop through with an increasing radius and see if nearby points are non-NA
          found <- 0
          radius <- .005
          while(found != 1 & radius < 1.5){ ## stop for max radius or match found
            ## sample 1000 nearby locs
            near <- matrix(runif(2000, -radius, radius), ncol = 2)
            near[, 1] <- near[, 1] + cent[1, 1]
            near[, 2] <- near[, 2] + cent[1, 2]
            
            ## extract raster pixels
            near <- data.table(raster::extract(pixel_id_raster, near), near)
            colnames(near) <- c('pixel_id', 'x', 'y')
            
            if(mean(is.na(near[,pixel_id])) < 1){ ## then we've found a non-NA neighbor
              found <- 1 ## end while loop
              message(sprintf('  ---- found neighbor using radius: %f', radius))
              
              ## find closest neighbor
              near <- na.omit(near) ## those with non-NA pixels
              dist <- sqrt((near[, x] - cent[1, 1]) ^ 2 +
                             (near[, y] - cent[1, 2]) ^ 2)
              min.ind <- which.min(dist) ## in case of tie, this returns 1st
              
              ## take the pixel id of the nearest sampled neighbor
              missing_locs[rr, pixel_id := near[min.ind, pixel_id] ]
            }
            
            ## increase radius in case we didn't catch anything
            radius <- radius + .005
            
          } ## end while loop
        } ## for each admin with NA pixel_id
      } ## if any NA pixel ids after centroid and poly methods
      
      ## Check to see if NAs still exist in missing locations, make a warning about missing areas.
      if(NA %in% missing_locs$pixel_id){
        message( "The following admin units appear to be getting lost forever:")
        print(missing_locs[is.na(pixel_id),c("pixel_id",admin_cols),with=F])
      }
      
      ## Merging on locations with pixel IDs
      missing_locs <- missing_locs[!is.na(pixel_id),c("pixel_id",admin_cols),with=F]
      missing_locs <- merge(missing_locs,cell_pred_current[,c("pixel_id","pop","year",draw_colnames),with=F],by="pixel_id",allow.cartesian=T)
      missing_admins[[as.character(lvl)]] <- missing_locs
    }
  } # For each level...
  
  ## Apply population density and prevalence threshold masks for sensitivity analyses
  message("Applying population density mask...")
  cell_pred_current[pop >= pop_threshold, colnames(as.data.table(cell_pred)) := 0]
  
  ## Apply standard LBD mask
  message("Applying standard LBD mask...")
  mask_master <- raster(<<<< FILEPATH REDACTED >>>>)
  mask_master <- projectRaster(mask_master, simple_raster, method = "ngb")
  mask_master <- raster::crop(mask_master, extent(simple_raster))
  mask_master <- setExtent(mask_master, simple_raster)
  
  pixel_id_mask <- seegSDM:::notMissingIdx(mask_master)
  cell_pred_current[pixel_id %in% pixel_id_mask, colnames(as.data.table(cell_pred)) := 0,]
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Collapsing down draw information using weighted means, based on population.
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  message("Generating data.tables with draws and stackers for each admin level")
  
  # Making a region-level draws object (we can use this for Africa/Global-level measures
  # without being concerned that we have duplicated pixels-- if we used Adm0, 1, 2, etc, we
  # would be double-counting pixels that we sampled for admin units that were smaller/did not
  # align with a pixel. By creating a region-level one ahead of time, we avoid that problem--
  # it means we should use "region" to aggreagate to continent/world/anything higher and nested.)
  
  if (region %in% c("lf_s_asia", "lf_se_asia", "lf_endem_afr", "lf_hispaniola", "oncho_endem_afr")) {
    cell_pred_current[, ADM0_CODE.final := ADM0_CODE.y]
    cell_pred_current[, ADM1_CODE.final := ADM1_CODE.x]
  }
  
  cell_pred_current$ADM0_CODE <- cell_pred_current$ADM0_CODE.final
  cell_pred_current$ADM1_CODE <- cell_pred_current$ADM1_CODE.final
  
  ## Adding in the missing admin information for admin 0,  1,  2...
  if ("0" %in% names(missing_admins)) {
    adm0 <- rbind(cell_pred_current, missing_admins[["0"]], fill = T)
  } else {
    adm0 <- copy(cell_pred_current)
  }
  
  ### Old behavior
  adm0 <- adm0[ADM0_CODE %in% region_adm0_list, lapply(.SD, weighted.mean, w = pop, na.rm = T), by = list(ADM0_CODE, year), .SDcols = draw_colnames]
  pop_adm0 <- cell_pred_current[ADM0_CODE %in% region_adm0_list, lapply(c("pop"), function(x) sum(get(x), na.rm = T)), by = list(year, ADM0_CODE), .SDcols = "pop"]
  setnames(pop_adm0, "V1", "pop")
  
  if ("1" %in% names(missing_admins)) {
    adm1 <- rbind(cell_pred_current, missing_admins[["1"]], fill = T)
  } else {
    adm1 <- copy(cell_pred_current)
  }
  
  adm1 <- adm1[ADM0_CODE %in% region_adm0_list, lapply(.SD, weighted.mean, w = pop, na.rm = T), by = list(ADM1_CODE, year), .SDcols = draw_colnames]
  pop_adm1 <- cell_pred_current[ADM0_CODE %in% region_adm0_list, lapply(c("pop"), function(x) sum(get(x), na.rm = T)), by = list(year, ADM1_CODE), .SDcols = "pop"]
  setnames(pop_adm1, "V1", "pop")
  
  if ("2" %in% names(missing_admins)) {
    adm2 <- rbind(cell_pred_current, missing_admins[["2"]], fill = T)
  } else {
    adm2 <- copy(cell_pred_current)
  }
  
  adm2 <- adm2[ADM0_CODE %in% region_adm0_list, lapply(.SD, weighted.mean, w = pop, na.rm = T), by = list(ADM2_CODE, year), .SDcols = draw_colnames]
  pop_adm2 <- cell_pred_current[ADM0_CODE %in% region_adm0_list, lapply(c("pop"), function(x) sum(get(x), na.rm = T)), by = list(year, ADM2_CODE), .SDcols = "pop"]
  setnames(pop_adm2, "V1", "pop")
  
  ## Generate total populations for each level, make sure to account for missing admins
  if ("0" %in% names(missing_admins)) {
    missing_adm0_pops <- subset(missing_admins[["0"]], !(ADM0_CODE %in% pop_adm0$ADM0_CODE), select = c("ADM0_CODE", "year", "pop"))
    pop_adm0 <- rbind(pop_adm0, missing_adm0_pops, use.names = T)
  }
  
  if ("1" %in% names(missing_admins)) {
    missing_adm1_pops <- subset(missing_admins[["1"]], !(ADM1_CODE %in% pop_adm1$ADM1_CODE), select = c("ADM1_CODE", "year", "pop"))
    pop_adm1 <- rbind(pop_adm1, missing_adm1_pops, use.names = T)
  }
  
  if ("2" %in% names(missing_admins)) {
    missing_adm2_pops <- subset(missing_admins[["2"]], !(ADM2_CODE %in% pop_adm2$ADM2_CODE), select = c("ADM2_CODE", "year", "pop"))
    pop_adm2 <- rbind(pop_adm2, missing_adm2_pops, use.names = T)
  }
  
  ## Adding results to the lists of results
  admin_0[[reg]] <- merge(adm0, pop_adm0, by = c("year", "ADM0_CODE"))
  admin_1[[reg]] <- merge(adm1, pop_adm1, by = c("year", "ADM1_CODE"))
  admin_2[[reg]] <- merge(adm2, pop_adm2, by = c("year", "ADM2_CODE"))
  
  ## Defining the hierarchy of what lives within what:
  sp_hierarchy <- unique(admin_levels[["2"]][["attributes"]][, list(ADM0_CODE, ADM1_CODE, ADM2_CODE, ADM0_NAME, ADM1_NAME, ADM2_NAME)])
  sp_hierarchy[, region := reg]
  sp_hierarchy_list[[reg]] <- sp_hierarchy
  
  ## Collapsing and Saving Results
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  message("Collapsing lists of data.tables into 1 data.table for analysis:")
  admin_0 <- rbindlist(admin_0)
  admin_1 <- rbindlist(admin_1)
  admin_2 <- rbindlist(admin_2)
  sp_hierarchy_list <- rbindlist(sp_hierarchy_list)
  
  adm0_list[[i]] <- copy(admin_0)
  adm1_list[[i]] <- copy(admin_1)
  adm2_list[[i]] <- copy(admin_2)
  
  rm(admin_0, admin_1, admin_2, adm0, adm1, adm2, pop, sp_hierarchy, cell_pred_current)
}

admin_0 <- rbindlist(adm0_list)
admin_1 <- rbindlist(adm1_list)
admin_2 <- rbindlist(adm2_list)

dir.create(<<<< FILEPATH REDACTED >>>>), showWarnings = FALSE)
save(admin_0, admin_1, admin_2, sp_hierarchy_list,
     file = paste0(<<<< FILEPATH REDACTED >>>>))

elapsed <- start - proc.time()
print(elapsed)
