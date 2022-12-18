
#####################################################################
### Takes a dataset with polygon records in specific format.
### Returns them as points using k-means resampling
### Uses the getPoints() function from seegMBG
#####################################################################



#####################################################################
## Wrapper functions to pull polygons ###############################
#####################################################################

pull_polys_foreach <- function(cores, shapefiles, shapes, shp_path) {
  
  message('Extracting all polygons from shapefiles on disk -- in parallel with foreach().\n')
  # set up cluster foreach-style
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  # distribute lib path and packages
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  clusterCall(cl, function(x)  {
    library(raster)
    library(magrittr)
  })
  
  # run in parallel by shapefile
  polys <- foreach (i = 1:length(shapefiles)) %dopar% {
    
    # pull shapefile
    shape <- shapefiles[i]
    message(paste0("Working on shapefile: ", shape))
    shp   <- shapefile(paste0(shp_path, '/', shape, '.shp'))
    
    # subsetting will break if NAs in row index (GAUL_CODE)
    shp <- shp[!is.na(shp$GAUL_CODE),]
    
    # HOTFIX for zambia shapefile with funny gaul codes
    if(shape == "GRED_Zambia"){
      shp@data$GAUL_CODE = gsub('894.2006.','',shp@data$GAUL_CODE)
    }
    
    # get location codes as numeric, not factor
    loc_codes <- unique(shapes[shapes$shapefile == shape,]$location_code) %>%
      as.character %>% as.numeric
    
    polys_subset <- list()
    problem_shapes <- c() #ensure errors due to problematic shapes are captured and reported to user
    
    for (j in 1:length(loc_codes)) {
      code <- loc_codes[j]
      poly_name <- paste0(shape, "__", code)
      
      if (code %in% as.character(shp$GAUL_CODE)) {
        poly <- shp[as.character(shp$GAUL_CODE) == code, ]
        polys_subset[[poly_name]] <- poly
      } else{
        problem_shapes <- c(problem_shapes, poly_name)
        warning(sprintf('GAUL code: %s not found in shapefile: %s',code,shape))
      }
      
    }
    
    if(length(problem_shapes) >= 1){
      polys_subset[["problem_shapes"]] <- problem_shapes
    }else{
      polys_subset[["problem_shapes"]] <- NA
    }
    return(polys_subset)
  }
  stopCluster(cl)
  return(polys)
}

pull_polys_mclapply <- function(cores, shapefiles, shapes, shp_path) {
  
  message('Extracting all polygons from shapefiles on disk -- in parallel with mclapply().\n')
  
  getpolys <- function(i){
    # pull shapefile
    shape <- shapefiles[i]
    message(paste0("Working on shapefile: ", shape))
    shp   <- shapefile(paste0(shp_path, '/', shape, '.shp'))
    
    # HOTFIX for zambia shapefile with funny gaul codes
    if(shape == "GRED_Zambia"){
      shp@data$GAUL_CODE = gsub('894.2006.','',shp@data$GAUL_CODE)
    }
    
    # get location codes as numeric, not factor
    loc_codes <- unique(shapes[shapes$shapefile == shape,]$location_code) %>%
      as.character  %>% as.numeric
    
    polys_subset <- list()
    problem_shapes <- c() #ensure errors due to problematic shapes are captured and reported to user
    
    for (j in 1:length(loc_codes)) {
      code <- loc_codes[j]
      poly_name <- paste0(shape, "__", code)
      
      if (code %in% as.character(shp$GAUL_CODE)) {
        poly <- shp[as.character(shp$GAUL_CODE) == code, ]
        polys_subset[[poly_name]] <- poly
      } else{
        problem_shapes <- c(problem_shapes, poly_name)
        warning(sprintf('GAUL code: %s not found in shapefile: %s',code,shape))
      }
      
    }
    
    if(length(problem_shapes) >= 1){
      polys_subset[["problem_shapes"]] <- problem_shapes
    }else{
      polys_subset[["problem_shapes"]] <- NA
    }
    return(polys_subset)
  }
  # Set multithreading to serial for `mclapply()`:
  set_serial_threads()
  polys <- mclapply(1:length(shapefiles),getpolys,mc.cores=cores)
  # Return to multithreading (if any):
  set_original_threads()
  return(polys)
}

pull_polys_fast <- function(shapefiles, shapes, shp_path) {
  
  message('Extracting all polygons from shapefiles on disk -- in serial from .rds files.\n')
  
  # run in serial by shapefile with fast poly loading function
  polys <- lapply(shapefiles, function(shape) {
    
    message(paste0("Working on shapefile: ", shape))
    shp   <- fast_load_shapefile(shape)
    
    # subsetting will break if NAs in row index (GAUL_CODE)
    shp <- shp[!is.na(shp$GAUL_CODE) & (shp$GAUL_CODE != ""),]
    
    # HOTFIX for zambia shapefile with funny gaul codes
    if(shape == "GRED_Zambia"){
      shp@data$GAUL_CODE = gsub('894.2006.','',shp@data$GAUL_CODE)
    }
    
    # get location codes as numeric, not factor
    loc_codes <- unique(shapes[shapes$shapefile == shape,]$location_code) %>%
      as.character %>% as.numeric
    
    polys_subset <- list()
    problem_shapes <- c() #ensure errors due to problematic shapes are captured and reported to user
    
    for (j in 1:length(loc_codes)) {
      code <- loc_codes[j]
      poly_name <- paste0(shape, "__", code)
      
      if (code %in% as.character(shp$GAUL_CODE)) {
        poly <- shp[as.character(shp$GAUL_CODE) == code, ]
        polys_subset[[poly_name]] <- poly
      } else{
        problem_shapes <- c(problem_shapes, poly_name)
        warning(sprintf('GAUL code: %s not found in shapefile: %s',code,shape))
      }
      
    }
    
    if(length(problem_shapes) >= 1){
      polys_subset[["problem_shapes"]] <- problem_shapes
    }else{
      polys_subset[["problem_shapes"]] <- NA
    }
    
    return(polys_subset)
  })
  return(polys)
}

#####################################################################
## Wrappers for old function names ##################################
#####################################################################

# To ensure backwards compatability


#' @title Add poly centroids
#'
#' @description Take a collapsed dataframe with polygon data and find the centroids of all raster cells which fall within each polygon.
#' 
#' @details Loops through all unique shapefile-location_code pairs in the input df, doing a spatial overlay of a raster and each polygon. The resolution of the raster can be controlled manually using the `disaggregation_factor`, or automatically through `auto_disaggregate`. `auto_disaggregate` will attempt a second overlay with a raster disaggregated by a factor of 5, followed by taking the centroid of the polygon if no raster centroids are found from the overlay attempts.
#' 
#' @note The `shapefile_col`, `location_code_col`, and `point` (if it exists in `df`) are typecast to characters at the beginning of the function and cast back to their original class at the end of the function.
#'
#' @param df A collapsed dataframe with a column with shapefile names and codes matching the GAUL_CODE columns of the shapefiles in the shapefile directory
#' 
#' @param shapefile_col string default `"shapefile"``, the df column name with shapefile names
#' 
#' @param location_code_col string default `"location_code"``, the df column name with location codes
#' 
#' @param fast_shapefiles boolean default `T``, pull from the rds shapefile directory? Drastically speeds up process.
#' 
#' @param disaggregation_factor int default `1`. The factor to disaggregate rasters by. The default of 1 uses a 5x5km raster. Disaggregating by 5 would result in a 1x1km raster, etc. 
#' 
#' @param auto_disaggregate boolean default `T`, If no raster centroids are found by overlaying a 5x5km raster, tries again with a 1x1km raster. If there are still no raster centroids, takes the centroid of the polygon. If `disaggregation_factor` is set to a value other than 1, the resolution of the raster in the first overlay attempt will be based on the `disaggregation_factor`. The second overlay attempt will be a raster dissagregated again by a factor of 5.
#'
#' @return df with an additional column "coordinates". Each entry in this column is a list containing a matrix with 2 columns (x and y coordinates) for each raster centroid found for the polygon in that row of df.
#'
add_poly_centroids <- function(df,
                               shapefile_col = "shapefile",
                               location_code_col = "location_code",
                               fast_shapefiles = T,
                               disaggregation_factor = 1,
                               auto_disaggregate = T) {
  
  #input validations
  if("coordinates" %in% names(df)) {
    warning("'coordinates' column found in df, overwriting")
    df[, coordinates := NULL]
  }
  if(!is.numeric(disaggregation_factor) || disaggregation_factor < 1) {
    stop("disaggregation factor must be a numeric value greater than or equal to 1")
  }
  if(!shapefile_col %in% names(df)) {
    stop("shapefile_col: ", shapefile_col, " is not a column name in df")
  }
  if(!location_code_col %in% names(df)) {
    stop("location_code_col: ", location_code_col, " is not a column name in df")
  }
  
  df <- data.table(df)
  setnames(df, c(shapefile_col, location_code_col), c("shapefile", "location_code"))
  
  #save classes for type conversion at the end
  shapefile_class <- class(df$shapefile)
  location_code_class <- class(df$location_code)
  
  #shapefile and location code need to be characters to work
  df <- cast_column(df, "shapefile", "character")
  df <- cast_column(df, "location_code", "character")
  
  #pull pop mask for reference raster
  mask_folder <- <<<< FILEPATH REDACTED >>>>
  ref_raster <- raster(paste0(mask_folder,'global_files/global_mask_master.tif'))
  ref_raster <- setValues(ref_raster, rep(1, length(ref_raster)))
  
  #make a table of unique shapefile/location_codes in df as a merge table
  if("point" %in% names(df)) {
    point_class <- class(df$point)
    df <- cast_column(df, "point", "character")
    shape_table <- unique(df[point == "0", .(shapefile, location_code)])
  } else {
    shape_table <- unique(df[, .(shapefile, location_code)])
  }
  
  shape_table[shapefile == "", shapefile := NA]
  shape_table <- shape_table[!is.na(shapefile) & !is.na(location_code)]
  setorder(shape_table, shapefile, location_code)
  
  #for each unique shapefile
  for(shape in unique(shape_table$shapefile)){
    message(paste0("Working on shape: ", shape))
    #get all the associated codes
    codes <- shape_table[shapefile == shape]$location_code
    if(fast_shapefiles){
      shp <- fast_load_shapefile(shape) 
    } else {
      shp <- rgdal::readOGR(paste0(<<<< FILEPATH REDACTED >>>>, shape))
    }
    missing_codes <- c()
    #for each unique code for a shapefile
    for(code in codes){
      #subset the shapefile to the specific polygon
      sub_shp <- subset(shp, GAUL_CODE == code)
      sub_shp$GAUL_CODE <- as.character(sub_shp$GAUL_CODE)
      #crop the global raster to the shapefile
      sub_raster <- crop(ref_raster, extent(sub_shp), snap="out")
      
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
        
        #if disaggregating by a factor of 5 is not enough, take centroid of polygon
        if(length(coords) == 0) {
          centroid <- coordinates(sub_shp)
          coords <- as.matrix(centroid)
          colnames(coords) <- c("x", "y")
          rownames(coords) <- NULL
          if(length(coords) == 0) {
            missing_codes <- c(missing_codes, code)
          }
        }
      } else if(length(coords) == 0 & auto_disaggregate == F) {
        missing_codes <- c(missing_codes, code)
      }
      
      #add the matrix to the shape table
      shape_table[(shapefile == shape) & (location_code == code), coordinates := list(list(coords))]
    }
    
    # print out polygon codes that had no raster centroids
    if(length(missing_codes) > 0) {
      if(auto_disaggregate == T) {
        message("     After auto_disaggregate, no raster centroids found for these codes: ", paste(as.character(missing_codes), collapse=", "))
      } else {
        message("     No raster centroids found for these codes: ", paste(as.character(missing_codes), collapse=", "))
      }
    }
  }
  #add the coordinates back onto the datatable
  df <- merge(df, shape_table, by = c("shapefile", "location_code"), all.x = T)
  
  #return things to the way they were
  if("point" %in% names(df)) {
    df <- cast_column(df, "point", point_class)
  }
  df <- cast_column(df, "shapefile", shapefile_class)
  df <- cast_column(df, "location_code", location_code_class)
  setnames(df, c("shapefile", "location_code"), c(shapefile_col, location_code_col))
  
  return(df)
}


#' @title Cast column
#' @description Cast a column in a dataframe given a class
#'
#' @param df a dataframe or datatable
#' @param column the name of a column in `df`
#' @param class one of "character", "numeric", "integer", "factor" or "logical"
#' 
#' @return `df` with column `column` typecast to `class`
#'
cast_column <- function(df,
                        column,
                        class) {
  if(class == "character") {
    df[[column]] <- as.character(df[[column]])
  } else if(class == "numeric"){
    df[[column]] <- as.numeric(df[[column]])
  } else if(class == "integer") {
    df[[column]] <- as.integer(df[[column]])
  } else if(class == "factor") {
    df[[column]] <- as.factor(df[[column]])
  } else if(class == "logical") {
    df[[column]] <- as.logical(df[[column]])
  } else {
    message("For column ", column, ": ", class, " is not an available type to cast to")
  }
  return(df)
}
