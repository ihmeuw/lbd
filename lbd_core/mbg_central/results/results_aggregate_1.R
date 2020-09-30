#########################################################################
# Purpose: Generate admin-unit level draws, using pixel-level weights
# as the
#########################################################################

# FILES REQUIRED TO RUN THIS CODE
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' This code is run when the kick-off code (results_launch_0.R) kicks off this script
#' for each indicator-date of interest. Please see that script for the .csv input that
#' determines the parameters taken from the commandArgs() used in this script.
#' ``````````````````````````````````````````````````````````````````````````````````
#' Roughly, the process undertaken is as follows:
#'
#' 1.) Determine regions based on the names of the cell draw objects saved in the outputs folder to iterate over
#' 2.) Load in the draws object (raked or unraked)
#' 3.) Get pixel-level population using indexing based on the simple_raster used as a template for this region
#' 4.) Discover the Admin 0, 1, and 2 GAUL codes for each pixel, by rasterizing admin-unit shapefiles and using
#'     the cell indexing to get appropriate GAUL values.
#' 5.) Discover which admin units were excluded from the above process, and attempt to recover them, using
#'     centroid location as a first wave, and then the first cell that touches the polygon as a second
#'     strategy. Add these onto the cell-pred objects for the relevant level.
#' 6.) Collapse the different levels (regional, with all of the pixels and populations used to create aggregates),
#'     Admin 0, 1, and 2 to the population-weighted draw. Merge on the total population, and save the objects.
#'
#'                 The draw-level objects for region, admin 0, 1, and 2 are then ready to be processed by
#'                 a separate script (02_collapse_admin.R).

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Setting up the workspace
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  rm(list=ls()) # Clearing out the memory

  # Defining Parameters
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Getting parameters from arguments passed in the qsub call
      indicator_group<-commandArgs()[4]; message(indicator_group)
      indicator<-commandArgs()[5]; message(indicator)
      run_date<-commandArgs()[6]; message(run_date)
      raked<-commandArgs()[7]; message(raked)
      pop_measure<-commandArgs()[8]; message(pop_measure)
      overwrite<-commandArgs()[9]; message(overwrite)
      age<-as.numeric(commandArgs()[10]); message(age)
      holdout<-as.numeric(commandArgs()[11]); message(holdout)
      shapefile_version<-as.character(commandArgs()[12]); message(shapefile_version)

    # Testing: Uncomment these and edit them for your example.
      # indicator_group<-"child_growth_failure"
      # indicator<-"wasting_mod_b"
      # run_date<-"2017_07_14_23_59_01"
      # raked<-T
      # pop_measure<-'a0004t'
      # age<-0
      # holdout<-0
      # overwrite<-F

    # Setting the root, loading in libraries, and functions:
      setwd('<<< FILEPATH REDACTED >>>')
      for(function_script in list.files(getwd(),pattern="*_functions.R")){message(function_script);source(function_script)};message("Central Functions Loaded.")
      load_libs(c('data.table','raster','seegSDM','seegMBG','plyr','dplyr','rgdal'))

    # Determinining whether this is a re-run, stopping if overwrite is off and the file already exists.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(overwrite==T){
        if(file.exists('<<< FILEPATH REDACTED >>>')){
          message("This file already exists, this must be a re-run.")
        }
      }else{
        if(file.exists('<<< FILEPATH REDACTED >>>')){
          stop("Looks like you have already prepped that admin draw file; skipping this!")
        }
      }

    # Getting regions from cell draw file names
      cell_draw_files<-list.files('<<< FILEPATH REDACTED >>>',pattern=paste0(indicator,"_cell_draws_eb_bin0_.*_",holdout,".RData$"),full.names = TRUE)
      get_region<-function(file_path){
        return(regmatches(file_path,
                          regexec(paste0(indicator,"_cell_draws_eb_bin0_(.*?)","_",holdout,".RData"),file_path))[[1]][2])
      }

      modeling_regions<-lapply(cell_draw_files, get_region)
      modeling_regions<-unlist(modeling_regions[!is.na(modeling_regions)])
      message("Regions detected: ")
      message(paste(modeling_regions," "))

    # These are where draw-level, admin-level results for each region will be stored
      message("Starting empty lists for results")

        regions<-list()
        admin_0<-list()
        admin_1<-list()
        admin_2<-list()
        sp_hierarchy_list<-list()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Looping through Regions, aggregating data to Admin Units
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    start<-proc.time()

  for (reg in modeling_regions){
  message(paste0("Pulling numbers for indicator: ",indicator," region: ",reg," rundate: ",run_date," using the population measure ",pop_measure," for weighting."))
    pathaddin<-paste0('_bin',age,'_',reg,'_',holdout)

  # Loading Data
  #~~~~~~~~~~~~~~
    message("Loading in the draw-level data for this region; this may take a while (but probably not longer than 10 minutes).")
      # Load in cell-level results based on that same model run
      if(raked){
        cell_pred<-readRDS('<<< FILEPATH REDACTED >>>')
      }else{
        load('<<< FILEPATH REDACTED >>>')
      }

  # Getting the simple polygon and simple raster objects for this region alone
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    message("Getting the spatial objects associated with this region.")
  simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version), buffer = 0.4, subset_only = TRUE, shapefile_version = shapefile_version)
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
      if(!(length(pixel_id)==nrow(cell_pred)/16)){
        stop("Excuse me, but the number of valid pixels in your simple raster is not the same number of valid pixels in your cell_pred object. Look into this!")
      }else{
        message(paste0("Check passed: There are ",length(pixel_id),
                       " pixels in your simple_raster object, the same number of pixels in the cell_pred object for each year."))
      }

  # Pulling Population
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Pull 2000-2015 annual population brick using new covariates function
    message(paste0("Pulling ",pop_measure," population raster."))
    pop <- load_worldpop_covariate(template_raster = simple_polygon,
                                   pop_measure = pop_measure,
                                   pop_release = pop_release,
                                   start_year = 2000,
                                   end_year = 2015,
                                   interval = 12)
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

    region_adm0_list<-get_adm0_codes(reg, shapefile_version = shapefile_version) # Getting the adm0 GAUL codes, we can use this to make sure we don't accidentally include countries from buffers that shouldn't be in this region

    # Defining a function that will get the raster versions of each Admin level:
      GetAdmin<-function(admin_level,simple_raster, region_adm0_list){
        message(paste0("Loading admin level ",admin_level))
          admin_shp<-readOGR(dsn='<<< FILEPATH REDACTED >>>',
                             layer=paste0("g2015_2014_", admin_level))

        message("Rasterizing...")
          admin_rast<-rasterize_check_coverage(admin_shp,simple_raster,paste0("ADM",admin_level,"_CODE"))
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
  } # Closing loop of regions

elapsed<-start-proc.time()
print(elapsed)

# Collapsing and Saving Results
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  message("Collapsing lists of data.tables into 1 data.table for analysis:")
  regions<-rbindlist(regions)
  admin_0<-rbindlist(admin_0)
  admin_1<-rbindlist(admin_1)
  admin_2<-rbindlist(admin_2)
  sp_hierarchy_list<-rbindlist(sp_hierarchy_list)

  save(regions,admin_0,admin_1,admin_2,sp_hierarchy_list, file='<<< FILEPATH REDACTED >>>')

  # Checking to make sure it went through
    if(file.exists('<<< FILEPATH REDACTED >>>')){
      message("DONE!")
    }else{
      stop("WHOAH! File did NOT save.")
    }
