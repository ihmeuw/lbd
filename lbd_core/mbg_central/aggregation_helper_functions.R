
#' @title Read z mapping file
#' 
#' @description Read in a file that links z values to age values.
#'   
#' @param z_map_file String that is the path to a csv file that has 2 columns. 
#'   The first contains the z values and the second contains the age mapping.
#' @param type String indicating what z represents. Currently only age is supported.
#' 
#' @return A dataframe with 2 columns: the z value and corresponding age group.
read_z_mapping <- function(z_map_file, type="age") {
  ## NOTE: only works for age aggregation
  z_map <- read.csv(z_map_file)
  colnames(z_map) <- c("z", "value")
  if(type == "age") {
    ## Run a check that it matches the worldpop labeling scheme
    z_vals <- gsub("-", "", z_map$value) # remove any "-"
    z_vals <- str_pad(z_vals, width = 4, side = "left", pad = "0")
    
    pop_dir <- <<<< FILEPATH REDACTED >>>>
    pop_files <- list.files(pop_dir)
    # Restrict to "a#..." files
    pop_files <- pop_files[grep("^a[0-9]",pop_files)]
    
    valid_z_vals <- unique(substr(pop_files, 2, 5))
    
    if(!(all(z_vals %in% valid_z_vals))) {
      stop(paste0("z_vals must match worldpop age groups.\nFor example, a valid z-value is ",
                  valid_z_vals[6], ".\nSee ", pop_dir, " for more."))
    }
    z_map$value <- z_vals
    return(z_map)
  } else {
    stop("The type of z aggregation must be set to 'age'")
  }
}

#' @title Process input data
#' 
#' @description After data has been loaded, this function will take any aggregated
#'   data, determine components (e.g., points in the polygon or age groups within
#'   the larger age aggregation) and the associated weights. Additionally, adds 
#'   the pixel_id (the pixel that the long/lat falls into).
#'   
#' @param df The data frame load_input_data().
#' @param pop_release The population measure release to use.
#' @param interval_mo The number of months between time steps.
#' @param modeling_shapefile_version String identifying version of of shapefile used in modeling.
#' @param poly_ag Logical that indicates whether the new point-polygon method 
#'   should be used.
#' @param zcol String that indicates the name of the z column if present.
#' @param zcol_ag String that indicates the name of the aggregated z column if 
#'   present.
#' @param zcol_ag_id String that indicates the name of the aggregated z column 
#'   id if present.
#' @param z_map Dataframe with mapping from z values to worldpop age groups.
#' @param z_ag_mat Matrix that contains user-supplied aggregation weights for
#'   aggregated z data (those that are specified in zcol_ag_id).
#' @param shapefile_col String indicating name of shapefile column.
#' @param location_code_col String indicating name of location code column.
#' @param fast_shapefiles Logical. Should fast shapefile loading be used?
#' @param disaggregation_factor Integer. The number of cells the raster corresponding
#'   to a polygon should be broken into. This can be useful if the polygon is very
#'   small in area.
#' @param auto_disaggregate Logical. Should a disaggregation_factor of 5 (i.e. breaking
#'   a 5km x 5km pixel to 1km x 1km) be applied if no pixel centroids fall inside the 
#'   polygon?
#' 
#' @return A data frame with the non aggregated data plus for the aggregated data, 
#'   multiple rows "disaggregated". The original columns will be present plus the
#'   following additional columns: row_id, agg_weight, first_entry, and pixel_id.
process_input_data <- function(df, 
                               pop_release, 
                               interval_mo, 
                               modeling_shapefile_version,
                               poly_ag=FALSE, 
                               zcol="z_column_default_blank", 
                               zcol_ag=NULL, 
                               zcol_ag_id=NULL, 
                               z_map=NULL,
                               z_ag_mat = NULL, 
                               shapefile_col = "shapefile",
                               location_code_col = "location_code",
                               fast_shapefiles = T,
                               disaggregation_factor = 1,
                               auto_disaggregate = T) {
  ## Inputs: df: df with i rows prior to polygon resampling. Should have
  ##             shapefile, location code, zcol (if used), zcol_ag (if used),
  ##             year, other data columns
  ## Outputs: uncollapsed_df: df that has the following columns: 
  ##            - ID (i unique values) corresponding to which row of the 
  ##              collapsed data the uncollapsed data corresponds to
  ##              NOTE: for aggregated data, there will be multiple rows in 
  ##              the uncollapsed data that correspond to the same original 
  ##              data
  ##            - long/lat: for point data this is what was originally provided
  ##              For polygon data, this is the centroids of pixels that fall
  ##              within the polygon
  ##            - z (only if original df has z): For aggregated z data, this
  ##              is broken down into the relevant zs that make up the aggregation
  ##            - t: orignal df time
  ##            - agg_weight: 1 for non z-aggregated point data. For polygon
  ##              data this is based on (total) population at the given time.
  ##              For z aggregated this is based on population at given time
  ##              and age group at location OR can be user-specified. For
  ##              polygon and z-aggregated this is based on population and age
  ##              group at the given time
  
  
  #input validations
  if(!is.numeric(disaggregation_factor) || disaggregation_factor < 1) {
    stop("disaggregation factor must be a numeric value greater than or equal to 1")
  }
  if(!shapefile_col %in% names(df) & poly_ag) {
    stop("shapefile_col: ", shapefile_col, " is not a column name in df")
  }
  if(!location_code_col %in% names(df) & poly_ag) {
    stop("location_code_col: ", location_code_col, " is not a column name in df")
  }
  
  # Check if any polygon data has latitude information
  if("latitude" %in% names(df) & poly_ag) {
    if(any(!is.na(df$latitude[df$point==0]))) {
      stop("Some of the polygon data already has longitude and latitude.")
    }
  }
  
  df <- data.table(df)
  
  # get ordering so that df can be reordered at end
  df$row_id <- 1:nrow(df)
  
  #pull pop mask for reference raster
  mask_folder <- <<<< FILEPATH REDACTED >>>>
  ref_raster <- raster(paste0(<<<< FILEPATH REDACTED >>>>))
  ref_raster <- setValues(ref_raster, rep(1, length(ref_raster)))
  
  # rename zcol so it doesn't get overwritten when we add in disaggregated zs
  if(zcol != "z_column_default_blank") {
    setnames(df, zcol, "z_tmp") 
  } 
  
  # rename zcol_ag for easy referencing
  if(!is.null(zcol_ag)) {
    setnames(df, zcol_ag, "z_ag")
  } else {
    df$z_ag <- NA
  }
  # rename zcol_ag_id for easy referencing
  if(!is.null(zcol_ag_id)) {
    setnames(df, zcol_ag_id, "z_ag_id")
  } else {
    df$z_ag_id <- NA
  }
  
  if(poly_ag) { # when using new polygon aggregation method
    setnames(df, c(shapefile_col, location_code_col), c("shapefile", "location_code"))

    #save classes for type conversion at the end
    shapefile_class <- class(df$shapefile)
    location_code_class <- class(df$location_code)
    
    #shapefile and location code need to be characters to work
    df <- cast_column(df, "shapefile", "character")
    df <- cast_column(df, "location_code", "character")
    
    ## Step 1. Process the polygon data 
    #          (includes polygon-single z and polygon-aggregated z)
    shape_table <- unique(df[, .(shapefile, location_code, point, z_ag, year)])
    shape_table[shapefile == "", shapefile := NA]
    
    # Verify there are shapefiles and location codes for all polygon data
    if(any(is.na(shape_table[point==0,.(shapefile, location_code)]))) {
      stop("There is missing shapefile and/or location code information for some of the polygon data")
    }
    
    # Subset to just polygon data
    shape_table <- shape_table[!is.na(shapefile) & !is.na(location_code)]
    setorder(shape_table, shapefile, location_code)
    
    if(nrow(shape_table) > 0) {
      message("Working on polygon data")
      
      for(shape in unique(shape_table$shapefile)) {
        message(paste0("Working on shape: ", shape))
        shape_table <- process_single_shapefile(shape, shape_table, ref_raster, 
                                                fast_shapefiles, pop_release, 
                                                interval_mo, z_map,
                                                disaggregation_factor, 
                                                auto_disaggregate)
      }
      # add the new data back onto the datatable
      df <- merge(df, shape_table, by = c("shapefile", "location_code", "point", "z_ag", "year"), all.x = T)
      
    }
  }
  
  
  ## Step 2. Age only data
  age_table <- unique(df[, .(country, latitude, longitude, z_ag, z_ag_id, year)])
  age_table <- age_table[!is.na(z_ag) & !is.na(latitude)]
  
  if(nrow(age_table) > 0 ) {
    message("Working on age aggregated data")
    
    if(!is.null(z_ag_mat)) {
      age_table_mat <- age_table[!is.na(z_ag_id), ] # only rows that will use z_ag_mat
      age_table <- age_table[is.na(z_ag_id), ] # rows that will not use z_ag_mat
      for(i in 1:length(age_table_mat$z_ag)) {
        id <- age_table_mat$z_ag_id[i]
        ages_needed <- eval(parse(text=as.character(age_table_mat$z_ag[i])))
        weight.tmp <- as.numeric(z_ag_mat[id, ages_needed])
        weight.tmp <- weight.tmp/sum(weight.tmp) # ensure sum to 1
        age_vals <- ages_needed
        age_table_mat[i, c("agg_weight","z") := list(list(weight.tmp),list(age_vals))]
      } 
      # add the new data back onto the datatable
      df[age_table_mat, on= c("country", "latitude", "longitude", "z_ag", "z_ag_id", "year"), 
         c("z","agg_weight"):=list(i.z,i.agg_weight)] 
      
    } 
    if(nrow(age_table) > 0) { # will need to use worldpop data
      age_table <- unique(age_table[, .(country, latitude, longitude, z_ag, year)])
      for(reg in unique(age_table$country)) {
        adm0_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
        simple_polygon_list <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, tolerance = 0.4,
                                                   shapefile_version = modeling_shapefile_version)
        simple_polygon      <- simple_polygon_list[[2]]
        age_gps <- unique(age_table[country == reg]$z_ag)
        for(ages in age_gps) {
          age_table <- process_single_z_gp(ages, reg, age_table, z_map,
                                           simple_polygon,
                                           pop_release, interval_mo,
                                           type="age")
        }
      }
      # add the new data back onto the datatable
      df[age_table, on= c("country", "latitude", "longitude", "z_ag", "year"), 
         c("z","agg_weight"):=list(i.z,i.agg_weight)] 
      
    }
  }
  
  
  ## reorder df
  df <- df[order(row_id),]
  
  
  cols_include <- !names(df) %in% c("latitude", "longitude", "z_tmp", "z_ag",
                                    "z_ag_id", "coordinates.long", "coordinates.lat", 
                                    "agg_weight", "z")
  ## Step 3. Create uncollapsed df
  message("Finalizing...")
  uncollapsed_df <- lapply(1:nrow(df), function(i) {
    df_row <- df[i, ]
    
    agg_weight <- df_row$agg_weight[[1]]
    if(is.null(agg_weight[1])) { # point data
      agg_weight <- 1
    }
    
    if(all(is.na(agg_weight))) { # if the population were all 0s --> weights are NaN
      message(paste0("Observation ", i, " had NA weights... dropping..."))
      to_drop <- T
    } else {
      to_drop <- F
    }
    
    # remove entries that have agg_weight of 0 since they do not contribute
    # agg weights with NAs will cause problems here but this should only happen 
    #   if all are NAs and therefore will be dropped due to code above
    idx_nonzero_aw <- agg_weight > 0 | is.na(agg_weight)
    agg_weight <- agg_weight[idx_nonzero_aw]
    
    if(!is.na(df_row$latitude)) {
      latitude <- df_row$latitude
      longitude <- df_row$longitude
    } else {
      latitude <- df_row$coordinates.lat[[1]][idx_nonzero_aw]
      longitude <- df_row$coordinates.long[[1]][idx_nonzero_aw]
    }
    if(!is.na(df_row$z_ag)) {
      z <- df_row$z[[1]][idx_nonzero_aw]
    } else {
      if(zcol != "z_column_default_blank") {
        z <- df_row$z_tmp
      } else {
        z <- NA
      }
    }
    
    first_entry <- c(1, rep(0, sum(idx_nonzero_aw)-1)) # useful for not plotting duplicates
    cbind(df_row[,..cols_include], longitude, latitude, 
          z, agg_weight, first_entry, to_drop)
  })
  
  uncollapsed_df <- do.call("rbind", uncollapsed_df)
  
  uncollapsed_df <- uncollapsed_df[!uncollapsed_df$to_drop,]
  
  uncollapsed_df$pixel_id <- raster::extract(ref_raster, 
                                             cbind(uncollapsed_df$longitude,
                                                    uncollapsed_df$latitude),
                                             cellnumbers=T)[,"cells"]
  
  if(zcol != "z_column_default_blank") {
    setnames(uncollapsed_df, "z", zcol) # return to original name
  }
  
  return(data.table(uncollapsed_df))
}

#' @title Process single z group for aggregation
#' 
#' @description Determines aggregation weights for aggregated z data using 
#'   worldpop numbers. 
#'   
#' @param ages String that indicates the z values that are being aggregated over.
#' @param reg String, name of country.
#' @param age_table Data table that has unique z aggregations, countries, years,
#'   longitude, and latitude.
#' @param z_map Dataframe with mapping from z values to worldpop age groups.
#' @param simple_polygon The simple polygon of the modeling region.
#' @param pop_release The population measure release to use.
#' @param interval_mo The number of months between time steps.
#' @param type String indicating what z represents. Currently only age is supported.
#' 
#' @return age_table that has 2 additional columns: z (list of the values that 
#'   make up z_ag) and agg_weight (list of the aggregation weights for each z).
process_single_z_gp <- function(ages, reg, age_table, z_map,
                                simple_polygon,
                                pop_release, interval_mo,
                                type="age") {
  
  # User does not specify weights
  # Need to read in from worldpop
  if(type == "age") {
    
    ages_needed <- eval(parse(text=as.character(ages)))
    
    for(year_needed in unique(age_table[z_ag == ages & country == reg]$year)) {
      worldpop_cov_rasters <- lapply(ages_needed, function(age) {
        suppressMessages(load_worldpop_covariate(simple_polygon,
                                covariate = 'worldpop',
                                pop_measure = paste0("a",z_map$value[z_map$z == age],"t"),
                                pop_release = pop_release,
                                start_year = year_needed,
                                end_year = year_needed,
                                interval = as.numeric(interval_mo))$worldpop)
        
      })
      worldpop_stack <- stack(worldpop_cov_rasters)
      longlats <- age_table[country==reg & z_ag==ages & year == year_needed, .(longitude, latitude)]
      loc_rel <- SpatialPoints(as.matrix(longlats), 
                               crs(simple_polygon))
      pop_vals <- extract(worldpop_stack, loc_rel)
      if(any(is.na(pop_vals))) {
        stop(paste0("There are missing population values for country: ", reg,
                    " in year ", year_needed))
      }
      pop_vals <- pop_vals/apply(pop_vals,1,sum)
      age_vals <- ages_needed
      # iterate through each long/lat
      for(i in 1:nrow(longlats)) {
        age_table[country==reg & z_ag==ages & year == year_needed & 
                    longitude==longlats$longitude[i] & latitude==longlats$latitude[i], 
                  c("agg_weight","z") := list(list(as.vector(pop_vals[i,])),list(age_vals))]
      }
    }
    
  } else {
    stop("The type of z aggregation must be set to 'age'")
  }
  return(age_table)
}

#' @title Process single shapefile for polygon and potentially age aggregation
#' 
#' @description Determines aggregation weights and points for polygon and potentially
#'   age aggregated data.
#'   
#' @param shape String that indicates the shapefile that should be read in.
#' @param shape_table Data table that has unique shapefiles, location codes, 
#' years, and z_ag.
#' @param ref_raster The master raster.
#' @param fast_shapefiles Dataframe with mapping from z values to worldpop age groups.
#' @param pop_release The population measure release to use.
#' @param interval_mo The number of months between time steps.
#' @param z_map Dataframe with mapping from z values to worldpop age groups.
#' @param disaggregation_factor Integer. The number of cells the raster corresponding
#'   to a polygon should be broken into. This can be useful if the polygon is very
#'   small in area.
#' @param auto_disaggregate Logical. Should a disaggregation_factor of 5 (i.e. breaking
#'   a 5km x 5km pixel to 1km x 1km) be applied if no pixel centroids fall inside the 
#'   polygon?
#' 
#' @return shape_table that has 4 additional columns: coordinates.long and
#'   coordinates.lat (lists of the longitude and latitude that make up the polygon),
#'   z (list of the values that make up z_ag), and agg_weight (list of the 
#'   aggregation weights for each point and z).
process_single_shapefile <- function(shape, shape_table, ref_raster, fast_shapefiles, 
                                     pop_release, interval_mo, z_map,
                                     disaggregation_factor, auto_disaggregate) {
  codes <- unique(shape_table[shapefile == shape]$location_code)
  # All coordinates are the same for a given shape and code
  if(fast_shapefiles){
    shp <- fast_load_shapefile(shape) 
  } else {
    shp <- rgdal::readOGR(paste0(<<<< FILEPATH REDACTED >>>>))
  }
  for(code in codes) {
    #subset the shapefile to the specific polygon
    sub_shp <- subset(shp, GAUL_CODE == code)
    sub_shp$GAUL_CODE <- as.character(sub_shp$GAUL_CODE)
    sub_raster <- crop(ref_raster, extent(sub_shp), snap="out")
    coords <- obtain_coords_given_code(sub_shp, sub_raster, disaggregation_factor, auto_disaggregate, code)
    
    
    if(nrow(coords) != 0 ) {
      coords_sp <- SpatialPoints(coords, crs(sub_shp))
      ag_ages <- unique(shape_table[shapefile == shape & location_code == code]$z_ag)
      for(ag_age in ag_ages) {
        if(is.na(ag_age)) {
          years <- shape_table[shapefile == shape & location_code == code & is.na(z_ag)]$year
          ages_needed <- NA
        } else {
          years <- shape_table[shapefile == shape & location_code == code & z_ag == ag_age]$year
          ages_needed <- eval(parse(text=as.character(ag_age)))
        }
        
        for(year_needed in years) {
          # Get population values
          worldpop_cov_rasters <- lapply(ages_needed, function(age) {
            suppressMessages(load_worldpop_covariate(sub_raster,
                                                    covariate = 'worldpop',
                                                    pop_measure = ifelse(is.na(age), "total", paste0("a",z_map$value[z_map$z == age],"t")),
                                                    pop_release = pop_release,
                                                    start_year = year_needed,
                                                    end_year = year_needed,
                                                    interval = as.numeric(interval_mo))$worldpop)
          })
          worldpop_stack <- stack(worldpop_cov_rasters)
          pop_vals <- raster::extract(worldpop_stack, coords_sp)
          pop_vals <- data.frame(pop_vals)
          
          # If all population values are NA and/or 0 for any age group or total
          # Note: the cellnumber column should not have any 0s or NAs
          if(any(apply(is.na(pop_vals) | pop_vals == 0,2,all))) { 
            stop(paste0("There are no population values at any of the pixel centroids for shape: ",
                        shape, " and location code: ", code))
          # If some population values are NA for any age group
          } else if(any(apply(is.na(pop_vals),2,any))) {
            message(paste0("There is at least one pixel centroid for shape: ",
                           shape, " and location code: ", code, " that is missing population. Setting to 0."))
            pop_vals[is.na(pop_vals)] <- 0 # assume there is no one living there
          }
          
          pop_vals <- as.matrix(pop_vals/sum(pop_vals)) # scale so that it sums to 1
          
          age_vals <- rep(ages_needed, each=nrow(pop_vals))
          pop_vals <- unlist(c(pop_vals)) # change to vector
          
          # add in population values as weights
          if(is.na(ag_age)) {
            shape_table[shapefile == shape & location_code == code & is.na(z_ag) & year==year_needed, 
                        c("coordinates.long", "coordinates.lat", "agg_weight","z") := 
                          list(list(rep(coords[,1],length(ages_needed))),
                               list(rep(coords[,2], length(ages_needed))),
                               list(pop_vals),list(age_vals))]
          } else {
            shape_table[shapefile == shape & location_code == code & z_ag == ag_age & year==year_needed, 
                        c("coordinates.long", "coordinates.lat", "agg_weight","z") := 
                          list(list(rep(coords[,1], length(ages_needed))),
                               list(rep(coords[,2], length(ages_needed))),
                               list(pop_vals),list(age_vals))]
          }
        }
      }
    } else { # no coordinates found
      shape_table[(shapefile == shape) & (location_code == code), 
                  c("coordinates.long", "coordinates.lat", "agg_weight", "z")  := 
                    list(list(NA), list(NA), list(NA), list(NA))]
      
    }
    
    
  }
  return(shape_table)
}

#' @title Obtain coordinates of pixels in a polygon
#' 
#' @description Determines which pixel centroids fall inside the polygon of interest.
#'   
#' @param sub_shp The polygon of interest
#' @param sub_raster The master raster cropped to the polygon of interest.
#' @param disaggregation_factor Integer. The number of cells the raster corresponding
#'   to a polygon should be broken into. This can be useful if the polygon is very
#'   small in area.
#' @param auto_disaggregate Logical. Should a disaggregation_factor of 5 (i.e. breaking
#'   a 5km x 5km pixel to 1km x 1km) be applied if no pixel centroids fall inside the 
#'   polygon?
#' @param code The location code for the polygon
#' 
#' @return Coordinates that fall into the polygon of interest.
obtain_coords_given_code <- function(sub_shp, sub_raster, disaggregation_factor, 
                                     auto_disaggregate, code) {
  #crop the global raster to the shapefile
  missing_codes <- c()
  
  #disaggregate raster to increase resolution
  dis_raster <- disaggregate(sub_raster, disaggregation_factor)
  
  #get centroids of raster
  centers <- coordinates(dis_raster)
  #convert centroids to spatialPoints object
  centers <- SpatialPoints(centers, crs(sub_shp))
  
  #overlay the centroid points with the polygon
  overlay <- centers[sub_shp,]
  #pull out a matrix of coordinates from the overlayed object
  coords <- overlay@coords
  
  #if no raster cell centroids within polygon are found,
  #disaggregate raster by a factor of 5 and try again
  if(length(coords) == 0 & auto_disaggregate == T) {
    #disaggregate raster to increase resolution
    dis_raster <- disaggregate(dis_raster, 5)
    
    #get centroids of raster
    centers <- coordinates(dis_raster)
    #convert centroids to spatialPoints object
    centers <- SpatialPoints(centers, crs(sub_shp))
    
    #overlay the centroid points with the polygon
    overlay <- centers[sub_shp,]
    #pull out a matrix of coordinates from the overlayed object
    coords <- overlay@coords
  }

  #if no matches, take centroid of polygon
  if(length(coords) == 0) {
    centroid <- coordinates(sub_shp)
    coords <- as.matrix(centroid)
    colnames(coords) <- c("x", "y")
    rownames(coords) <- NULL
    if(length(coords) == 0) {
      missing_codes <- c(missing_codes, code)
    }
  }
  
  # print out polygon codes that had no raster centroids
  if(length(missing_codes) > 0) {
    if(auto_disaggregate == T) {
      message("     After auto_disaggregate, no raster centroids found for these codes: ", paste(as.character(missing_codes), collapse=", "))
    } else {
      message("     No raster centroids found for these codes: ", paste(as.character(missing_codes), collapse=", "))
    }
  }
  
  return(coords)
}
