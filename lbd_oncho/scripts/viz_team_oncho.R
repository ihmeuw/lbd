####################################################################
### Inputs for Viz Team (oncho-speific)
####################################################################

rm(list = ls())

## Get current user name
user <- Sys.info()[["user"]] ## Get current user name

## Set repo location and indicator group
core_repo <- <<<< FILEPATH REDACTED >>>>
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

library(matrixStats)

## ##############################################################
## set the output directory for all objects & other parameters ##
## ##############################################################

out_dir <- <<<< FILEPATH REDACTED >>>>

viz_date <- "viz_date"
dir.create(paste0(out_dir, viz_date))

diag <- "micro"

region <- "oncho_endem_afr"
indicator <- "had_oncho_w_resamp"
indicator_group <- "oncho"
run_date <- "run_date"
predict_years <- 2000:2018

raster_agg_factor <- 1

modeling_shapefile_version <- "2020_05_21"

## ###################################################################################################
## Run functions ##
## ###################################################################################################
message("Creating count and prev tables")

cases <- make_case_adm_tables(ind = indicator, ig = indicator_group, rd = run_date, pop_measure = "total")

ig <- "oncho"
ind <- "had_oncho_w_resamp"
sharedir <- <<<< FILEPATH REDACTED >>>>

years <- c(2000:2018)
year_list <- min(years):max(years)
interval_mo <- 12

method_agg <- "pixel"
pathadd <- "_pixel"

# saving outputs

all_mods <- list(cases)
summary_list <- c("lower_prev", "mean_prev", "upper_prev", "lower_cases", "mean_cases", "upper_cases")

message(paste0("Saving count and prev tables in: ", <<<< FILEPATH REDACTED >>>>))
dir.create(<<<< FILEPATH REDACTED >>>>)
for (summ in summary_list) {
  summ_allmods <- lapply(all_mods, FUN = function(mod) mod[[summ]]) %>% rbindlist()

  if (length(grep("cases", summ)) > 0) {
    summary_type <- "count_"
  } else {
    summary_type <- "prev_"
  }

  for (i in predict_years) {
    if (length(grep("lower", summ)) > 0) file_name <- paste0(summary_type, "lower_", i)
    if (length(grep("mean", summ)) > 0) file_name <- paste0(summary_type, "mean_", i)
    if (length(grep("upper", summ)) > 0) file_name <- paste0(summary_type, "upper_", i)

    output_file <- copy(summ_allmods[year == i]) %>% as.data.table()
    setnames(output_file, "gaul_code", "id")
    write.csv(output_file[, .(id, value)], <<<< FILEPATH REDACTED >>>>, row.names = F)
  }
}

message("Creating individual year rasters for count and prev")

dir.create(<<<< FILEPATH REDACTED >>>>)
dir.create(<<<< FILEPATH REDACTED >>>>)
lapply(run_date, make_case_rasters_parallel, diag = diag, ind = indicator, ig = indicator_group, pop_measure = "total", output_dir = <<<< FILEPATH REDACTED >>>>, indic_repo = indic_repo)

########## Remove aggregate results from adm1 and adm2 that are at least 50% masked by standard LBD population mask
admin1_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 1, version = modeling_shapefile_version))
admin2_shp <- rgdal::readOGR(dsn = get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))

msk <- raster(<<<< FILEPATH REDACTED >>>>)
lks <- raster(<<<< FILEPATH REDACTED >>>>)
lks[!is.na(lks)] <- 1
lks <- extend(lks, msk)
lks <- raster::resample(lks, msk)

combined_mask <- msk

freq_original <- data.table("id"=integer(), "raw_count"=integer(), "masked_count"=integer())

for (i in 1:length(run_date)) {
  load(<<<< FILEPATH REDACTED >>>>)
  
  cropped_mask <- crop(combined_mask, simple_polygon)
  cropped_mask <- resample(cropped_mask, simple_raster)
  
  # Loading and Assigning Admin Units to Pixels
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  message("Getting the spatial location (admin unit) of each of the pixel locations.")
  
  region_adm0_list <- get_adm0_codes(region, shapefile_version = modeling_shapefile_version) # Getting the adm0 GAUL codes, we can use this to make sure we don't accidentally include countries from buffers that shouldn't be in this region
  
  # Defining a function that will get the raster versions of each Admin level:
  GetAdmin<-function(admin_level,simple_raster, region_adm0_list, shapefile_version){
    message(paste0("Loading admin level ",admin_level))
    
    # load admin shape file
    admin_shp <- rgdal::readOGR(dsn=get_admin_shapefile(admin_level, version = shapefile_version))
    
    # ensure that the rasterize variable is a numeric
    admin_shp@data[[paste0('ADM', admin_level, '_CODE')]] <- as.numeric(as.character(admin_shp@data[[paste0('ADM', admin_level, '_CODE')]]))
    
    # if it doesn't exist, get areas of polygons. 
    if(is.null(admin_shp$Shape_Area)){
      admin_shp$Shape_Area <- area(admin_shp) / 1e6 ## TODO
    }
    
    message("Rasterizing with the custom function...")
    # we order by area so small places don't get buried under big places (e.g. Lesotho and S. Africa)
    link_table <- modeling_shapefile_version
    admin_rast<-rasterize_check_coverage(admin_shp[order(admin_shp$Shape_Area),],simple_raster,paste0("ADM",admin_level,"_CODE"), link_table = link_table, fun="first")
    
    message("Converted to raster based on simple_raster template. Cropping and masking:")
    admin_rast  <- crop(admin_rast, extent(simple_raster))
    admin_rast  <- setExtent(admin_rast, simple_raster)
    admin_rast  <- mask(admin_rast, simple_raster)
    
    message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
    admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% region_adm0_list,]
    admin_centroids<-SpatialPointsDataFrame(gCentroid(admin_shp, byid=TRUE), admin_shp@data, match.ID=FALSE)
    
    message("Compiling and returning results.")
    admin<-list()
    admin[["spdf"]]<-admin_shp
    admin[["centroids"]]<-admin_centroids
    admin[["rast"]]<-admin_rast
    admin[["attributes"]]<-copy(data.table(admin_shp@data))
    
    return(admin)
  }
  
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  pixel_spatial <- data.table(pixel_id = pixel_id)
  
  message("Rasterizing shapefiles; this may take a while.")
  admin_levels<-list() # Empty list of levels that will be filled with admin levels
  for(lvl in c(0,1,2)){
    fieldname<-paste0("ADM",lvl,"_CODE")
    admin_info<-GetAdmin(admin_level=lvl,simple_raster,region_adm0_list, shapefile_version = modeling_shapefile_version)
    pixel_spatial[[fieldname]]<-extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
    admin_levels[[as.character(lvl)]]<-admin_info # Add the admin info to the list
    if(sum(is.na(pixel_spatial[[fieldname]]))>0){ # Check to see if any of the pixels don't have a location assigned
      message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level ",lvl))
    }
  }
  
  for (lvl in c(0, 1, 2)) {
    freqs <- as.data.table(freq(admin_levels[lvl + 1][[1]]$rast))
    colnames(freqs) <- c("id", "raw_count")
    
    masked <- mask(admin_levels[lvl + 1][[1]]$rast, cropped_mask, inverse = TRUE)
    freqs2 <- as.data.table(freq(masked))
    colnames(freqs2) <- c("id", "masked_count")
    
    freqs_merged <- merge(freqs, freqs2, by="id", all=TRUE)
    
    freq_original <- rbind(freq_original, freqs_merged)
  }
}

freq_original <- freq_original[!is.na(id)]
freq_original[, prop_unmasked := (masked_count / raw_count)] 

admin0_to_drop <- c("61")
admin1_to_drop <- c("21217", "3217", "6217", "1097", "1095", "1093", "5061", "3061", "1061")
admin2_to_drop <- c("8026067", "10026067", "7026067", "13010067", "3024067", "7024067", "4014067", "4010193", "1021217", "1003217", "1006217", "3004146", "1004146", "1001153", "3011153", "3012153", "2005061", "1005061", "4005061", "1003061", "1001061")
all_admins_to_drop <- c(admin0_to_drop, admin1_to_drop, admin2_to_drop)

keep <- c("2043118")

drop <- freq_original[((prop_unmasked <= 0.5 & id > 256) | (is.na(masked_count)) | (id %in% all_admins_to_drop)) & !(id %in% keep),]

for (file in list.files(<<<< FILEPATH REDACTED >>>>, recursive = FALSE)) {
  message(file)
  if (grepl("count", file) | grepl("prev", file) | grepl("postprob", file)) {
    current <- fread(<<<< FILEPATH REDACTED >>>>)
    current <- current[!(id %in% drop$id)]
    write.csv(current, <<<< FILEPATH REDACTED >>>>, row.names = F)
  }
}

###################################################
### Post hoc: Check CSVs for duplicate admins

library(data.table)
diags <- c("micro")
for (diag in diags) {
  for (file in list.files(<<<< FILEPATH REDACTED >>>>), recursive = FALSE)) {
    message(file)
    if (grepl("count", file) | grepl("prev", file) | grepl("postprob", file)) {
      current <- fread(<<<< FILEPATH REDACTED >>>>)
      print(nrow(current))
      current <- unique(current)
      print(nrow(current))
      write.csv(current, <<<< FILEPATH REDACTED >>>>, row.names = F)
    }
  }
}
