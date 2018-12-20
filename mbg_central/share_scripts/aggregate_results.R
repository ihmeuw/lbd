# Make Admin 1 and 2 level summary draw-level objects

# Setup
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

# Getting parameters from arguments passed in the qsub call
indicator<-commandArgs()[4]; message(indicator)
indicator_group<-commandArgs()[5]; message(indicator_group)
run_date<-commandArgs()[6]; message(run_date)
raked<-commandArgs()[7]; message(raked)
pop_measure<-commandArgs()[8]; message(pop_measure)
overwrite<-commandArgs()[9]; message(overwrite)
age<-as.numeric(commandArgs()[10]); message(age)
holdout<-as.numeric(commandArgs()[11]); message(holdout)
region<- as.character(commandArgs()[12]); message(region)

reg <- region #convenience

## Load libraries and miscellaneous MBG project functions.
user            <- Sys.info()['user']
repo            <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/mbg/',user)

remote          <- 'origin'
branch          <- 'develop'
#system(sprintf('cd %s\ngit pull %s %s', repo, remote, branch))

share_dir <- paste0("<<<< FILEPATH REDACTED >>>>>",indicator_group,"/",indicator,"/output/",run_date,"/")

setwd(repo)

root <- ifelse(Sys.info()[1]=="Windows", "<<<< FILEPATH REDACTED >>>>>", "<<<< FILEPATH REDACTED >>>>>")
package_lib    <- ifelse(grepl("geos", Sys.info()[4]),
                           paste0(root,'<<<< FILEPATH REDACTED >>>>>geos_packages'),
                           paste0(root,'<<<< FILEPATH REDACTED >>>>>/packages'))#    dependent on the machine where the user runs the code.
.libPaths(package_lib)                                  # Ensures packages look for dependencies here when called with library().

source('mbg_central/mbg_functions.R')                   # Functions to run MBG model.
source('mbg_central/prep_functions.R')                  # Functions to setup MBG run
source('mbg_central/covariate_functions.R')             # Functions to prep and transform 5*5 covariates
source('mbg_central/misc_functions.R')                  # Other computational MBG-related functions.
source('mbg_central/post_estimation_functions.R')
source('mbg_central/gbd_functions.R')
source('mbg_central/shiny_functions.R')
source('mbg_central/holdout_functions.R')
source('mbg_central/categorical_variable_functions.R')
source('mbg_central/validation_functions.R')
source('mbg_central/seegMBG_transform_functions.R')     

package_list <- c('foreign', 'rgeos', 'data.table','raster','rgdal',
                  'INLA','seegSDM','seegMBG','plyr','dplyr')
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

# Determinining whether this is a re-run, stopping if overwrite is off and the file already exists.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(overwrite==T){
  if(file.exists(paste0("<<<< FILEPATH REDACTED >>>>>",indicator_group,"/",indicator,"/output/",run_date,"/",indicator,"_admin_draws",ifelse(raked,"_raked.Rdata",".RData")))){
    message("This file already exists, this must be a re-run.")
  }
}else{
  if(file.exists(paste0("<<<< FILEPATH REDACTED >>>>>",indicator_group,"/",indicator,"/output/",run_date,"/",indicator,"_admin_draws",ifelse(raked,"_raked.Rdata",".RData")))){
    stop("Looks like you have already prepped that admin draw file; skipping this!")
  }
}

# These are where draw-level, admin-level results will be stored
message("Starting empty lists for results")
regions<-list()
admin_0<-list()
admin_1<-list()
admin_2<-list()
sp_hierarchy_list<-list()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Looping through Regions, aggregating data to Admin Units
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start <-proc.time()

message(paste0("Pulling numbers for indicator: ",indicator," region: ",reg," rundate: ",run_date," using the population measure ",pop_measure," for weighting."))
pathaddin<-paste0('_bin',age,'_',reg,'_',holdout)

# Loading Data
#~~~~~~~~~~~~~~
message("Loading in the draw-level data for this region; this may take a while (but probably not longer than 10 minutes).")
# Load in cell-level results based on that same model run
if(raked){
  # Different file conventions in use; try either rds or rdata until standardized
  filename_rds <- paste0(share_dir,reg,"_raked_cell_pred.RDs")
  filename_rdata <- paste0(share_dir, indicator, "_", ifelse(raked, "raked", "unraked"), "_", "cell_draws_eb_bin0_", reg, "_", holdout, ".RData")
  if (file.exists(filename_rds)) cell_pred<-readRDS(filename_rds)
  if (file.exists(filename_rdata)) {
      load(filename_rdata, verbose = T)
      cell_pred <- raked_cell_pred
      rm(raked_cell_pred)
    } 
} else {
  load(paste0("<<<<< FILEPATH REDACTED >>>>>",indicator_group,"/",indicator,"/output/",run_date,"/",indicator,"_cell_draws_eb",pathaddin,".RData"))
}

# Check if load was successful; stop if not
if (!exists("cell_pred")) stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")

# Getting the simple polygon and simple raster objects for this region alone
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Getting the spatial objects associated with this region.")
simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = FALSE)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
message("Building simple raster from subset_shape")
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]
rm(raster_list,simple_polygon_list,pop_raster);gc()
message("All done loading spatial template files (subset_shape,simple_polygon,simple_raster,pop_raster, etc)")

# Determining a list of the valid pixel indices based on the simple raster template
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pixel_id <- seegSDM:::notMissingIdx(simple_raster)
pixel_spatial<-data.table(pixel_id=pixel_id)
message("Pixel ID that links cell_pred observations to raster locations created from simple raster, checking for sameness in dimensions compared to cell_pred.")
if(!(length(pixel_id)==nrow(cell_pred)/17)){
  stop("Excuse me, but the number of valid pixels in your simple raster is not the same number of valid pixels in your cell_pred object. Look into this!")
}else{
  message(paste0("Check passed: There are ",length(pixel_id),
                 " pixels in your simple_raster object, the same number of pixels in the cell_pred object for each year."))
}

# Pulling Population
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Pull 2000-2015 annual population brick using new covariates function
message(paste0("Pulling ",pop_measure," population raster."))
pop <- load_and_crop_covariates_annual(covs = 'worldpop',
                                       measures = pop_measure, # Defined above
                                       simple_polygon = simple_polygon,
                                       start_year  = 2000,
                                       end_year    = 2016,
                                       interval_mo = 12,
                                       agebin=1)
message("Ensuring that the spatial extents of the population matches the simple raster")
pop<-crop_set_mask(pop[[1]],simple_raster)
message("Done generating population rasterbrick; extracting relevant pixels based on simple_raster.")
pop <- data.table(extract(pop, pixel_id)) # getting the pixels from the population raster that correspond to the valid pixels in the pop-raster
message("Reshaping population from wide to long, and generating an id variable that corresponds to the simple raster.")
pop[,pixel_id:=pixel_id]
pop<-melt(pop,id.vars="pixel_id") # Melting the dataframe such that it is long () and should match the number of rows in cell_pred
pop[,year:=1999+as.numeric(substr(as.character(variable),start=10,stop=15))] # Converting "worldpop.1" variables to actual years.
pop<-pop[,list(pixel_id,year,pop=value)] # Subsetting columns
pop[is.na(pop),pop:=0] # Setting values where pop is NA to 0
message("Checking whether the population has the same dimensions as cell_pred")
if(!(nrow(pop)==nrow(cell_pred))){
  stop("Excuse me, but the number of rows in your population table is *not* the same as the number of rows in your cell predictions. Look into this!")
}else{
  message("Check passed: The population file has the same number of rows as the cell predictions after transformation to long format.")
}


# Loading and Assigning Admin Units to Pixels
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Getting the spatial location (admin unit) of each of the pixel locations.")

region_adm0_list<-get_gaul_codes(reg) # Getting the adm0 GAUL codes, we can use this to make sure we don't accidentally include countries from buffers that shouldn't be in this region

# Defining a function that will get the raster versions of each Admin level:
GetAdmin<-function(admin_level,simple_raster, region_adm0_list){
  message(paste0("Loading admin level ",admin_level))

    if (admin_level %in% c(0,1)) {
      admin_shp<-readOGR(dsn=paste0("<<<< FILEPATH REDACTED >>>>>AUL_admin/admin", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, ".shp"),
                         layer=paste0("g2015_2014_", admin_level))
    } else if (admin_level == 2) {
      # Use most up to date admin2 code
      admin_shp <- readOGR(dsn="<<<< FILEPATH REDACTED >>>>>/g2015_2014_2_modified.shp",
                           layer="g2015_2014_2_modified")
    }

    if(admin_level==2){ # Assigning disputed areas to countries for admin 2 units
      admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==61013]=133
      admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==40760]=40765
      admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==40762]=145
      admin_shp[[paste0('ADM', admin_level,'_CODE')]][admin_shp[[paste0('ADM', admin_level,'_CODE')]]==102  ]=6
    }

  message("Rasterizing...")
    admin_rast<-rasterize(admin_shp[order(admin_shp$Shape_Area),],simple_raster,paste0("ADM",admin_level,"_CODE"), fun="first")
  message("Converted to raster based on simple_raster template. Cropping and masking:")
    admin_rast  <- crop(admin_rast, extent(simple_raster))
    admin_rast  <- setExtent(admin_rast, simple_raster)
    admin_rast  <- mask(admin_rast, simple_raster)
  message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
    admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% region_adm0_list,]
    admin_centroids<-SpatialPointsDataFrame(gCentroid(admin_shp, byid=TRUE),
                                            admin_shp@data, match.ID=FALSE)
    message("Loaded in.")
  message("Compiling and returning results.")
    admin<-list()
    admin[["spdf"]]<-admin_shp
    admin[["centroids"]]<-admin_centroids
    admin[["rast"]]<-admin_rast
    admin[["attributes"]]<-copy(data.table(admin_shp@data))
  return(admin)
}

message("Rasterizing shapefiles; this may take a while.")
admin_levels<-list() # Emtpy list of levels that will be filled with admin levels
  for(lvl in c(0,1,2)){
    fieldname<-paste0("ADM",lvl,"_CODE")
    admin_info<-GetAdmin(admin_level=lvl,simple_raster,region_adm0_list)
    pixel_spatial[[fieldname]]<-extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
    admin_levels[[as.character(lvl)]]<-admin_info # Add the admin info to the list
    if(sum(is.na(pixel_spatial[[fieldname]]))>0){ # Check to see if any of the pixels don't have a location assigned
      message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level ",lvl))
      }
  }

# Merging together cell_pred and spatial information, population information.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pop<-merge(pop,pixel_spatial,by="pixel_id",all.x=T) # Merging on the spatial information to the population information.
pop<-pop[order(year,pixel_id)] # Re-ordering the pop object by year such that pixels ascent, and years ascend (same as cell_pred)
cell_pred<-data.table(cell_pred) # Converting cell_pred to a data.table
draw_colnames<-names(cell_pred)
cell_pred<-cbind(cell_pred,pop) # Adding on the population and spatial information to cell_pred.


# Loading and Assigning Pixels to Admin Units (that are missing)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# We now need to figure out what admin units don't show up when we pull the locations of individual raster pixels of our results.
# To do this, we need to know what Admin 1 and 2 units exist in the shapefiles we want to use, but don't show up in the results.

message("Generating a raster of the pixel_id values")
  pixel_id_raster<-copy(simple_raster)
  for (pix in pixel_id){ # For each of the valid pixels, replace the value with the pixel_id, not whatever country ID is used in the simple_polygon/raster
    pixel_id_raster@data@values[pix]<-pix
  }

missing_admins<-list() # This will contain the cell_preds for otherwise missing admin units, by level.
for(lvl in c(0,1,2)){
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
    # Strategy 1: Assign a pixel location based on centroid location
    # Develop a raster of the pixel-IDs:
    message("  Discovering centroid location")
      points<-admin_levels[[as.character(lvl)]][["centroids"]] # SpatialPointsDataFrame of that level
      missing_points<-points[(points@data[[fieldname]] %in% missing_codes),] # getting only the missing points
      missing_centroid_locs<-data.table(raster::extract(pixel_id_raster,missing_points)) # Extracting the missing location values as centroid points...
      names(missing_centroid_locs)<-"point_method"
      missing_admins_centroids<-data.table(missing_codes,missing_centroid_locs)

    # Strategy 2: Assign a pixel location based on polygon extraction
    # Develop a raster of the pixel-IDs:
    message("  Discovering first raster pixel touching polygon")
      polys<-admin_levels[[as.character(lvl)]][["spdf"]] # SpatialPointsDataFrame of that level
      missing_polys<-polys[(polys@data[[fieldname]] %in% missing_codes),] # getting only the missing points
      missing_poly_locs<-data.table(raster::extract(x=pixel_id_raster,y=missing_polys,small=T,fun=function(x,...)first(x))) # Extracting the missing location values as polygons, pulling the first raster cell that it touches...
      names(missing_poly_locs)<-"poly_method"
      missing_admins_polys<-data.table(missing_codes,missing_poly_locs) # Extracting the

    # Merging strategies together: centroids and polygons; adding to list.
      missing_locs<-merge(missing_admins_polys,missing_admins_centroids,by="missing_codes")
      setnames(missing_locs,"missing_codes",fieldname)
      missing_locs<-merge(missing_locs,missing_table,by=fieldname) # Add in admin 0, 1 levels if relevant
      missing_locs[,pixel_id:=point_method] # If centroid produced something, go with centroid
      missing_locs[is.na(point_method),pixel_id:=poly_method] # Otherwise, go with the polygon method.

    # Check to see if NAs still exist in missing locations, make a warning about missing areas.
      if(NA %in% missing_locs$pixel_id){
        message( "The following admin units appear to be getting lost forever:")
        print(missing_locs[is.na(pixel_id),c("pixel_id",admin_cols),with=F])
      }

    # Merging on locations with pixel IDs
      missing_locs<-missing_locs[!is.na(pixel_id),c("pixel_id",admin_cols),with=F]
      missing_locs<-merge(missing_locs,cell_pred[,c("pixel_id","pop","year",draw_colnames),with=F],by="pixel_id",allow.cartesian=T)
      missing_admins[[as.character(lvl)]]<-missing_locs
    }

} # For each level...

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Collapsing down draw information using weighted means, based on population.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message("Generating data.tables with draws for each admin level")

# Making a region-level draws object (we can use this for Africa/Global-level measures
# without being concerned that we have duplicated pixels-- if we used Adm0, 1, 2, etc, we
# would be double-counting pixels that we sampled for admin units that were smaller/did not
# align with a pixel. By creating a region-level one ahead of time, we avoid that problem--
# it means we should use "region" to aggreagate to continent/world/anything higher and nested.)
region<-cell_pred[ADM0_CODE %in% region_adm0_list,lapply(.SD,weighted.mean,w=pop, na.rm=T),by=list(year), .SDcols=draw_colnames]
region[,region_name:=reg]

# Adding in the missing admin information for admin 0, 1, 2...
if("0" %in% names(missing_admins)){adm0<-rbind(cell_pred,missing_admins[["0"]], fill=T)}else{adm0<-copy(cell_pred)}
  adm0 <- adm0[ADM0_CODE %in% region_adm0_list,lapply(.SD,weighted.mean,w=pop, na.rm=T),by=list(ADM0_CODE,year), .SDcols=draw_colnames]
if("1" %in% names(missing_admins)){adm1<-rbind(cell_pred,missing_admins[["1"]], fill=T)}else{adm1<-copy(cell_pred)}
  adm1 <- adm1[ADM0_CODE %in% region_adm0_list,lapply(.SD,weighted.mean,w=pop, na.rm=T),by=list(ADM1_CODE,year), .SDcols=draw_colnames]
if("2" %in% names(missing_admins)){adm2<-rbind(cell_pred,missing_admins[["s"]], fill=T)}else{adm2<-copy(cell_pred)}
  adm2 <- adm2[ADM0_CODE %in% region_adm0_list,lapply(.SD,weighted.mean,w=pop, na.rm=T),by=list(ADM2_CODE,year), .SDcols=draw_colnames]

# Generate total populations for each level
pop_region<-pop[ADM0_CODE %in% region_adm0_list,list(pop=sum(pop)),by=list(year)]
pop_region[,region_name:=reg]
pop_adm0<-pop[ADM0_CODE %in% region_adm0_list,list(pop=sum(pop)),by=list(ADM0_CODE,year)]
pop_adm1<-pop[ADM0_CODE %in% region_adm0_list,list(pop=sum(pop)),by=list(ADM1_CODE,year)]
pop_adm2<-pop[ADM0_CODE %in% region_adm0_list,list(pop=sum(pop)),by=list(ADM2_CODE,year)]

# Adding results to the lists of results
regions[[reg]]<-merge(region,pop_region,by=c("year","region_name"))
admin_0[[reg]]<-merge(adm0,pop_adm0,by=c("year","ADM0_CODE"))
admin_1[[reg]]<-merge(adm1,pop_adm1,by=c("year","ADM1_CODE"))
admin_2[[reg]]<-merge(adm2,pop_adm2,by=c("year","ADM2_CODE"))

# Defining the hierarchy of what lives within what:
sp_hierarchy<-unique(admin_levels[["2"]][["attributes"]][,list(ADM0_CODE,ADM1_CODE,ADM2_CODE,ADM0_NAME,ADM1_NAME,ADM2_NAME)])
sp_hierarchy[,region:=reg]
sp_hierarchy_list[[reg]]<-sp_hierarchy

elapsed <- start - proc.time()
print(elapsed)

# Collapsing and Saving Results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
message("Collapsing lists of data.tables into 1 data.table for analysis:")
admin_0<-rbindlist(admin_0)
admin_1<-rbindlist(admin_1)
admin_2<-rbindlist(admin_2)
sp_hierarchy_list<-rbindlist(sp_hierarchy_list)

save(admin_0,admin_1,admin_2,sp_hierarchy_list, 
     file=paste0(share_dir, indicator, "_", ifelse(raked, "raked", "unraked"), "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))

# Finish up
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Write file to signal that this region is done
write(NULL, file = paste0(share_dir, "/fin_agg_", reg, "_", holdout, "_", age, "_", raked))

# Checking to make sure it went through
if(file.exists(paste0(share_dir, indicator, "_", ifelse(raked, "raked", "unraked"), "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"))){
  message("DONE!")
}else{
  stop("WHOAH! File did NOT save.")
}
