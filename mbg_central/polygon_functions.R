
# Functions for resampling polygons to points. 
#####################################################################
### Takes a dataset with polygon records in specific format.
### Returns them as points using k-means resampling
### Uses the getPoints() function from seegMBG
#####################################################################


#####################################################################
## Main resample_polygons() code ####################################
#####################################################################

## resample_polygons ################################################

#' Function to perform polygon resampling
#'
#' Resample polygons using the k-means clustering algorithm developed by RB et al
#'
#' @param data a data frame or data table with the following columns: cluster_id, exposed,
#'             <indicator>, latitude, longitude, shapefile, location code. Location_code refers to
#'             the shapefile polygon that the data. Data should alreaday be aggregated up
#'             to the smallest geog unit (if this is point, then lat long should have values).
#'             If lat long don't have values then shapefile and location_id should
#' @param shp_path Path to shapefiles (standard = our shapefile directory)
#' @param ignore_warnings if TRUE, and there are polys that are not found the function will warn but not stop.
#' @param cores Number of cores to use for parallel processing bits
#' @param indic Indicator being run (e.g. `died`, `dpt3_cov`, `has_lri`.  Note that this is only used for some very specific
#'              renaming at the very end of the function - this code will actually preserve all columns passed to it.
#' @param unique_vars character values of additional names that uniquely ID rows in data, for example "age_bin" for
#'                    u5m data. Defaults to NULL
#' @param density Passed to `seegMBG::getPoints()` as the `n` argument, which determines the # of integration points.
#'                Determines how many pseudopoints are laid down for each polygon. Default = 0.001.
#'                This can be interpreted as points per pixel if `perpixel = T`, and as points per square km if `perpixel = T`
#'                and `use_1k_popraster = T` (default settings)
#' @param perpixel Passed to `seegMBG::getPoints()` as the `perpixel` argument  This determines whether the `density`
#'                 arguments represents an absolute number, or if it is the number of integration points per pixel (so
#'                 that the total # of integration points = `density` * `[total number of pixels]`).  By default this is
#'                 TRUE with a density of 0.001.
#' @param prob Passed to `seegMBG::getPoints()` as the `prob` argument. Determines whether to weight the integration points
#'             by the values of the raster (here, the population raster). TRUE by default.
#' @param use_1k_popraster Should a 1km x 1km population raster be used (disaggregates the 5x5 km population raster)?
#'                         Logical, default = TRUE.  Otherwise will use 5k pop raster
#' @param pull_poly_method How should we get the polygons for resampling?  Options are parallel processing with `mclapply()`
#'                         (`pull_poly_method = "mclapply"`), parallel processing with `foreach()` (`pull_poly_method = "foreach"`)
#'                         or fast polygon pulling using RDS files (`pull_poly_method = "fast"` -- the default)
#' @param gaul_list List of gaul codes covering your modeling area, for use in pulling the population rasters
#' @param seed Optional seed that is set just prior to resampling to gaurantee the same results when rerunning on the same data.
#'             By default this is `NULL`, in which case no seed is set.
#'
#' @return A data table with your input `data` object, now with polygons duplicated and weighed according to the resampling
#'         algorithm
#'
#' @examples
#' # Will not run unless you have vaccine data
#' df_pointpoly <- resample_polygons(data = df_input,
#'                                   cores = 20,
#'                                   indic = "dpt3_cov",
#'                                   density = 0.001)
#'

resample_polygons <- function(data,
                              shp_path = '<<<< FILEPATH REDACTED >>>>>/Shapefile directory',
                              ignore_warnings  = TRUE,
                              cores            = as.numeric(slots),
                              indic            = indicator,
                              unique_vars      = NULL,
                              density          = 0.001,
                              perpixel         = TRUE,
                              prob             = TRUE,
                              use_1k_popraster = TRUE,
                              pull_poly_method = "fast",
                              gaul_list        = get_gaul_codes('africa'),
                              seed             = NULL)
{

  ############################
  # Confirm everything is set up
  require(doParallel); require(plyr); require(snow); require(doSNOW)
  data = data.frame(data) # df not dt for this one

  # Confirm necessary columns in data
  nec_cols <- c('shapefile', 'location_code')

  if (any(!(nec_cols %in% names(data)))){
    message('Missing the following columns in your data:')
    print(nec_cols[!nec_cols %in% names(data)])
  }

  # change lat long to latiude longitude if needed
  names(data)[names(data)=="lat"] <- "latitude"
  names(data)[names(data)=="long"] <- "longitude"

  # find records with shapefiles
  data$shapefile[!is.na(data$latitude)&!is.na(data$longitude)]=NA
  noloc=rep(FALSE,nrow(data))
  noloc[data$shapefile==""&!is.na(data$shapefile)]=TRUE

  # remove any spatially unidentifiable data
  message(paste('Dropping',sum(noloc),'of',nrow(data),'rows of data due to no spatial identifiers (lat, long, shapefile info missing)\n'))
  data=data[!noloc,]

  # keep index of only polygon data
  shp_idx <- which(!is.na(data$shapefile))
  message(paste(length(shp_idx),'of',nrow(data),'rows of data are polygon assigned\n'))

  # hotfix for a weird yemen thing
  data$shapefile[data$shapefile=="matched to GADM admin 1 shapefile"&data$country=="Yemen"]="YEM_adm1_GADM"

  # identify all shapefiles from the dataset
  all_shapes <- unique(data$shapefile[shp_idx])
  message(paste(length(all_shapes),'unique shapefiles found in data.\n'))

  # check they're in the directory
  message(paste0('Checking shapefile directory (',shp_path,') for matches.. '))
  if (!all(paste0(all_shapes, '.shp') %in% list.files(shp_path))){
    message('Missing the following shapefiles:')
    print(all_shapes[!(paste0(all_shapes, '.shp') %in% list.files(shp_path))])
    if(!ignore_warnings) stop('Function breaking because not all shapefiles are a match. Please')
  }  else {
    message('All shapefiles in data match a shapefile by name in the directory.\n')
  }

  ############################
  # load and extract all polygons - now in parallel

  # get unique shapefiles/ location codes
  shapes   <- unique(data[, c('shapefile', 'location_code')])
  shapes   <- shapes[!is.na(shapes$shapefile), ]
  n_shapes <- nrow(shapes)

  # sort by shapefile name (faster if they're all clumped together)
  o <- order(shapes$shapefile)
  shapes <- shapes[o, ]

  # empty, named list to store polygons
  polys <- list()
  polys[[n_shapes]] <- NA
  names(polys) <- paste(shapes$shapefile, shapes$location_code, sep = '__')

  # null first shapefile
  shape <- ''

  # report to the user
  message(sprintf('extracting %i polygons', n_shapes))

  # vector of shapefiles
  shapefiles <- unique(shapes$shapefile)

  if (pull_poly_method == "foreach") {
    polys <- pull_polys_foreach(cores, shapefiles, shapes, shp_path)
  } else if (pull_poly_method == "mclapply") {
    polys <- pull_polys_mclapply(cores, shapefiles, shapes, shp_path)
  } else if (pull_poly_method == "fast") {
    polys <- pull_polys_fast(shapefiles, shapes, shp_path)
  } else {
    stop ("pull_poly_method must be one of \"foreach\", \"mclapply\", or \"fast\"")
  }

  polys <- unlist(polys) # get to single-level list

  # find ones that didn't work
  bad_records <- which(sapply(polys, class) == 'NULL')
  if(length(bad_records)!=0){
    warning(sprintf('%i polygons could not be found:', length(bad_records)))
    print(names(bad_records))
    if(!ignore_warnings) {
      stop('Since ignore_warnings==FALSE, the function is now breaking. Please fix your bad records or set ignore_warnings==TRUE to drop them.\n')
    } else {
      warning('Since you have set ignore_warnings=TRUE, the function will drop all bad records and continue. Note: this is NOT recommended.\n')
      data = data[!paste(data$shapefile, data$location_code, sep = '__') %in%
                    names(bad_records),]
    }
  }

  ######################################
  ##
  message('In parallel, generating integration points for each shapefile and expanding records to integrate. \n')

  # get unique sets of shapefiles and location codes (bins optional)
  d       <- data[, c(unique_vars, 'shapefile', 'location_code')]
  pars    <- unique(d)
  pars    <- pars[!is.na(pars$shapefile), ]
  n_chunk <- nrow(pars)

  # get row indices
  message('Getting indices and chunking out data by polygon in parallel.\n')
  dx      <- trim(apply(d,   1,paste,collapse='--'))
  px      <- trim(apply(pars,1,paste,collapse='--'))
  indices <- mclapply(1:n_chunk,function(x){unname(which(dx==px[x]))},mc.cores=cores)

  chunks <- mclapply(1:n_chunk,function(x){data[indices[[x]],]},mc.cores=cores)

  # grab population raster
  message('Loading Population Raster.\n')
  suppressMessages(suppressWarnings(load_simple_polygon(gaul_list = gaul_list,buffer=0.4,subset_only=TRUE)))
  raster_list<-suppressMessages(suppressWarnings(build_simple_raster_pop(subset_shape)))
  if(use_1k_popraster){
    popraster <- disaggregate(raster_list[['pop_raster']],5) # needs to be 1km for density 0.001
  } else {
    popraster = raster_list[['pop_raster']]
  }

  # get points in parallel
  message('Running getPoints() on each polygon in parallel.\n')
  getPointsWrapper <- function(x){
    poly_name <- paste(pars[x,c('shapefile', 'location_code')], collapse = '__')
    poly <- polys[[poly_name]]
    
    if (is.null(poly)) {
      # make a dummy, empty points dataframe
      points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
    } else {
      if (!is.null(seed)) set.seed(seed)
      points <- try(
        getPoints(shape    = poly,
                  raster   = popraster[[3]],
                  n        = density,
                  perpixel = perpixel,
                  prob     = prob) )
      if(inherits(points,'try-error')){
        points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
      } else {
        colnames(points) <- c('longitude', 'latitude', 'weight')
      }
    }
    return(points)
  }
  chunk_points <- mclapply(1:n_chunk,getPointsWrapper,mc.cores=cores)

  # duplicate records, and add their integration points
  message('Duplicating polygon records for each of their new integration points.\n')
  getNewChunksWrapper <-function(x){
    # data for this chunk (shapefile/age bin combo)
    chunk  <- chunks[[x]]
    points <- chunk_points[[x]]

    if (nrow(points) == 0) {
      new_chunk <- chunk[0, ]
      warning(paste('Chunk',x,'missing spatial info and will be dropped.'))
    } else {

      dupRecords <- function(j) {
        record <- chunk[j, , drop = FALSE]       # pull out the record
        # drop columns also in "points"
        record <- record[, !(colnames(record) %in% colnames(points))]
        # faster than rbind / replicate
        record_dup <- cbind(record, points, row.names = NULL)

        return(record_dup)
      }

      duped_records <- lapply(1:nrow(chunk), dupRecords)

      # append the duplicated data
      new_chunk <- rbindlist(duped_records)
    }
    return(new_chunk)
  }

  new_chunks <- mclapply(1:n_chunk,getNewChunksWrapper,mc.cores=cores)

  # remove old chunks from data
  message('Finishing up..\n')
  idx_all <- unlist(indices)
  data <- data[-idx_all, ]

  # append new chunks
  new_chunks_all <- rbindlist(new_chunks, fill = TRUE)
  new_chunks_all$pseudocluster = TRUE
  if(nrow(data) > 0) {
    data$weight=1
    data$pseudocluster= FALSE
    data <- rbind(data, new_chunks_all)
  } else {
      data <- new_chunks_all
   }

  # count rows with no lat/longs and remove them
  no_ll <- which(is.na(data$latitude) | is.na(data$longitude))
  message(paste('Dropping',length(no_ll),'rows with missing spatial information.'))
  if(length(no_ll) != 0) data <- data[-no_ll, ]

  data <- data.table(data)
  return(data)

}


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

    for (j in 1:length(loc_codes)) {
      code <- loc_codes[j]

      if (code %in% as.character(shp$GAUL_CODE)) {
        poly <- shp[as.character(shp$GAUL_CODE) == code, ]
      } else{
        warning(sprintf('GAUL code: %s not found in shapefile: %s',code,shape))
        poly <- NULL
      }

      poly_name <- paste0(shape, "__", code)
      polys_subset[[poly_name]] <- poly

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

    for (j in 1:length(loc_codes)) {
      code <- loc_codes[j]

      if (code %in% as.character(shp$GAUL_CODE)) {
        poly <- shp[as.character(shp$GAUL_CODE) == code, ]
      } else{
        warning(sprintf('GAUL code: %s not found in shapefile: %s',code,shape))
        poly <- NULL
      }

      poly_name <- paste0(shape, "__", code)
      polys_subset[[poly_name]] <- poly

    }
    return(polys_subset)
  }
  polys <- mclapply(1:length(shapefiles),getpolys,mc.cores=cores)
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

    for (j in 1:length(loc_codes)) {
      code <- loc_codes[j]

      if (code %in% shp$GAUL_CODE) {
        poly <- shp[shp$GAUL_CODE == code, ]
      } else{
        warning(sprintf('GAUL code: %s not found in shapefile: %s',code,shape))
        poly <- NULL
      }

      poly_name <- paste0(shape, "__", code)
      polys_subset[[poly_name]] <- poly

    }

    return(polys_subset)
  })
  return(polys)
}

#####################################################################
## Wrappers for old function names ##################################
#####################################################################

# To ensure backwards compatability

resample_polygons_fast  <-   function(data,
                                      shp_path = '<<<< FILEPATH REDACTED >>>>>/Shapefile directory',
                                      ignore_warnings = TRUE,
                                      cores           = as.numeric(slots),
                                      indic           = indicator,
                                      unique_vars     = NULL,
                                      density         = 0.001,
                                      perpixel        = TRUE,
                                      prob            = TRUE,
                                      use_1k_popraster= TRUE) {

  warning("resample_polygons_fast() is deprecated. Use resample_polygons() instead")
  message("Running resample_polygons() with pull_poly_method = \"fast\"")
  rp <- resample_polygons(data = data,
                          shp_path = shp_path,
                          ignore_warnings = ignore_warnings,
                          cores = cores,
                          indic = indic,
                          unique_vars = unique_vars,
                          density = density,
                          perpixel = perpixel,
                          prob = prob,
                          use_1k_popraster = use_1k_popraster,
                          pull_poly_method = "fast")
  return(rp)
}

resample_polygons_dev   <-   function(data,
                                      shp_path = '<<<< FILEPATH REDACTED >>>>>/Shapefile directory',
                                      ignore_warnings = TRUE,
                                      cores           = as.numeric(slots),
                                      indic           = indicator,
                                      unique_vars     = NULL,
                                      density         = 0.001,
                                      perpixel        = TRUE,
                                      prob            = TRUE,
                                      use_1k_popraster= TRUE) {

  warning("resample_polygons_dev() is deprecated. Use resample_polygons() instead")
  message("Running resample_polygons() with pull_poly_method = \"mclapply\"")
  rp <- resample_polygons(data = data,
                          shp_path = shp_path,
                          ignore_warnings = ignore_warnings,
                          cores = cores,
                          indic = indic,
                          unique_vars = unique_vars,
                          density = density,
                          perpixel = perpixel,
                          prob = prob,
                          use_1k_popraster = use_1k_popraster,
                          pull_poly_method = "mclapply")
  return(rp)
}

#####################################################################
## Preserving for backwards compatability, although we should wrap
## this functionality into the main function at some point if still
## being used...
#####################################################################

resample_polygons_bysamplingtype  <-   function(data,
                                 shp_path = '<<<< FILEPATH REDACTED >>>>>/Shapefile directory',
                                 ignore_warnings = TRUE,
                                 cores           = as.numeric(slots),
                                 indic           = indicator,
                                 unique_vars     = NULL,
                                 density         = 0.001,
                                 perpixel        = TRUE,
                                 prob            = TRUE,
                                 use_1k_popraster= TRUE,
                                 raster_dict     = NULL,
                                 raster_col      = NULL,
                                 regions = NULL)
{

  ## ##########################
  ## Confirm everything is set up
  require(doParallel); require(plyr); require(snow); require(doSNOW)
  data = data.frame(data) # df not dt for this one

  ## make sure things are ok if using raster_dict and raster_col
  if(!is.null(raster_dict)){

    if(is.null(raster_col))
      stop("If using raster_dict must identify raster_col")

    ## check all entries in raster col match entries in raster_dict

  }

  # Confirm necessary columns in data ------------------------------------------- may need further review ers------------
  nec_cols <- c('sampling', 'shapefile', 'location_code', 'N', 'cluster_id')

  if (any(!(nec_cols %in% names(data)))){
    message('Missing the following columns in your data:')
    print(nec_cols[!nec_cols %in% names(data)])
  } else {
    message('All necessary fields in your data! Your future is bright.')
  }

  # change lat long to latiude longitude if needed
  names(data)[names(data)=="lat"] <- "latitude"
  names(data)[names(data)=="long"] <- "longitude"


  # find records with shapefiles
  data$shapefile[!is.na(data$latitude)&!is.na(data$longitude)]=NA #this was using lat and long ers
  noloc=rep(FALSE,nrow(data))
  noloc[data$shapefile==""&!is.na(data$shapefile)]=TRUE
  noloc[is.na(data$latitude) & is.na(data$longitude) & is.na(data$shapefile)] = TRUE

  # remove any spatially unidentifiable data
  message(paste('Dropping',sum(noloc),'of',nrow(data),'rows of data due to no spatial identifiers (lat, long, shapefile info missing)\n'))
  data=data[!noloc,]

  # keep index of only polygon data
  shp_idx <- which(!is.na(data$shapefile))
  message(paste(length(shp_idx),'of',nrow(data),'rows of data are polygon assigned\n'))

  # hotfix for a weird yemen thing
  data$shapefile[data$shapefile=="matched to GADM admin 1 shapefile"&data$country=="Yemen"]="YEM_adm1_GADM"

  # identify all shapefiles from the dataset
  all_shapes <- unique(data$shapefile[shp_idx])
  message(paste(length(all_shapes),'unique shapefiles found in data.\n'))

  # check they're in the directory
  message(paste0('Checking shapefile directory (',shp_path,') for matches.. '))
  if (!all(paste0(all_shapes, '.shp') %in% list.files(shp_path))){
    message('Missing the following shapefiles:')
    print(all_shapes[!(paste0(all_shapes, '.shp') %in% list.files(shp_path))])
    if(!ignore_warnings) stop('Function breaking because not all shapefiles are a match. Please')
  }  else {
    message('All shapefiles in data match a shapefile by name in the directory.\n')
  }

  ############################
  # load and extract all polygons - now in parallel

  # get unique shapefiles/ location codes
  shapes   <- unique(data[, c('shapefile', 'location_code')])
  shapes   <- shapes[!is.na(shapes$shapefile), ]
  n_shapes <- nrow(shapes)

  # sort by shapefile name (faster if they're all clumped together)
  o <- order(shapes$shapefile)
  shapes <- shapes[o, ]

  # empty, named list to store polygons
  polys <- list()
  polys[[n_shapes]] <- NA
  names(polys) <- paste(shapes$shapefile, shapes$location_code, sep = '__')

  # null first shapefile
  shape <- ''

  # report to the user
  message('Extracting all polygons from shapefiles on disk -- in parallel.\n')
  message(sprintf('extracting %i polygons', n_shapes))

  # set up cluster foreach-style
  cl <- makeCluster(cores)
  registerDoParallel(cl)

  # distribute lib path and packages
  clusterCall(cl, function(x) .libPaths(x), .libPaths())

  clusterCall(cl, function(x)  {
    library(raster)
    library(magrittr)
  })

  # vector of shapefiles
  shapefiles <- unique(shapes$shapefile)

  # run in parallel by shapefile
  polys <- foreach (i = 1:length(shapefiles)) %dopar% {

    # pull shapefile
    shape <- shapefiles[i]
    message(paste0("Working on shapefile: ", shape))
    shp   <- shapefile(paste0(shp_path, '/', shape, '.shp'))

    # subsetting will break if NAs in row index (GAUL_CODE)
    shp <- shp[!is.na(shp$GAUL_CODE),]

    # get location codes as numeric, not factor
    loc_codes <- unique(shapes[shapes$shapefile == shape,]$location_code) %>%
      as.character %>% as.numeric

    polys_subset <- list()

    for (j in 1:length(loc_codes)) {
      code <- loc_codes[j]

      if (code %in% shp$GAUL_CODE) {
        poly <- shp[shp$GAUL_CODE == code, ]
      } else{
        warning(sprintf('GAUL code: %s not found in shapefile: %s',code,shape))
        poly <- NULL
      }

      poly_name <- paste0(shape, "__", code)
      polys_subset[[poly_name]] <- poly

    }

    return(polys_subset)
  }

  polys <- unlist(polys)

  # find ones that didn't work
  bad_records <- which(sapply(polys, class) == 'NULL')
  if(length(bad_records)!=0){
    warning(sprintf('%i polygons could not be found:', length(bad_records)))
    print(names(bad_records))
    if(!ignore_warnings) {
      stop('Since ignore_warnings==FALSE, the function is now breaking. Please fix your bad records or set ignore_warnings==TRUE to drop them.\n')
    } else {
      warning('Since you have set ignore_warnings=TRUE, the function will drop all bad records and continue. Note: this is NOT recommended.\n')
      data = data[!paste(data$shapefile, data$location_code, sep = '__') %in%
                    names(bad_records),]
    }
  }

  ######################################
  ##
  message('In parallel, generating integration points for each shapefile and expanding records to integrate. \n')

  data <- as.data.table(data)

  ## get unique sets of shapefiles and location codes (bins optional)
  if(is.null(raster_dict)){
    d       <- data[, c(unique_vars, 'shapefile', 'location_code')]
    pars    <- unique(d)
    pars    <- pars[!is.na(pars$shapefile), ]
    n_chunk <- nrow(pars)
  }else{
    d       <- data[, c(unique_vars, 'shapefile', 'location_code', raster_col), with = FALSE]
    pars    <- unique(d)
    pars    <- pars[!is.na(pars$shapefile), ]
    n_chunk <- nrow(pars)
  }

  ## get row indices
  message('Getting indices and chunking out data by polygon in parallel.\n')
  dx      <- trim(apply(d,   1,paste,collapse='--'))
  px      <- trim(apply(pars,1,paste,collapse='--'))
  indices <- mclapply(1:n_chunk,function(x){unname(which(dx==px[x]))},mc.cores=cores)

  chunks <- mclapply(1:n_chunk,function(x){data[indices[[x]],]},mc.cores=cores)

  # grab population raster
  message('Loading Population Raster.\n')
  suppressMessages(suppressWarnings(load_simple_polygon(gaul_list = (regions),buffer=0.4,subset_only=TRUE))) 
  raster_list<-suppressMessages(suppressWarnings(build_simple_raster_pop(subset_shape)))
  if(use_1k_popraster){
    popraster <- disaggregate(raster_list[['pop_raster']],5) 
  } else {
    popraster = raster_list[['pop_raster']]
  }

  ## load other rasters if need be
  if(!is.null(raster_dict)){
    message("Loading non-population rasters")
    ##raster_list <- NULL
    for(rr in 1:nrow(raster_dict)){
      raster_list[[as.character(raster_dict[rr, 1])]] <- raster(as.character(raster_dict[rr, 2]))
    }
  }

  ## get points in parallel -----------------
  message('Running getPoints() on each polygon in parallel.\n')
  if(is.null(raster_dict)){
    getPointsWrapper <- function(x){
      poly_name <- paste(pars[x,c('shapefile', 'location_code')], collapse = '__')
      poly <- polys[[poly_name]]

      ## NOTE: simply resampling on 2010 population data for now for convenience. Could redo by year, but wouldnt make much of difference anyway.
      if (is.null(poly)) {
        ## make a dummy, empty points dataframe
        points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
      } else {
        points <- try(
          getPoints(shape    = poly,
                    raster   = popraster[[3]],
                    n        = density,
                    perpixel = perpixel,
                    prob     = prob) )
        if(inherits(points,'try-error')){
          points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
        } else {
          colnames(points) <- c('longitude', 'latitude', 'weight')
        }
      }
      return(points)
    }
    chunk_points <- mclapply(1:n_chunk,getPointsWrapper,mc.cores=cores)
  }else{
    getPointsWrapper.w.dict <- function(x){
      which.rast <- as.character(x[2])
      x <- as.numeric(x[1])
      poly_name <- paste(pars[x,c('shapefile', 'location_code')], collapse = '__')
      poly <- polys[[poly_name]]

      if (is.null(poly)) {
        ## make a dummy, empty points dataframe
        points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
      } else {
        rast <-
          points <- try(
            getPoints(shape    = poly,
                      raster = (if(length(raster_list[[which.rast]]) == 1) raster_list[[which.rast]][[1]] else raster_list[[which.rast]][[1]]),
                      n        = density,
                      perpixel = perpixel,
                      prob     = prob) )
        if(inherits(points,'try-error')){
          points <- data.frame(longitude = NA, latitude  = NA, weight    = NA)[0, ]
        } else {
          colnames(points) <- c('longitude', 'latitude', 'weight')
        }
      }
      return(points)
    }
    ## match raster_dict to n_chunk and make a matrix to apply over
    n.rast <- data.frame(n = 1:n_chunk, rast = pars[, raster_col, with=FALSE])
    chunk_points <- apply(n.rast, 1, getPointsWrapper.w.dict)
  }

  # duplicate records, and add their integration points
  message('Duplicating polygon records for each of their new integration points.\n')
  getNewChunksWrapper <-function(x){
    chunk  <- chunks[[x]]
    points <- chunk_points[[x]]

    if (nrow(points) == 0) {
      new_chunk <- chunk[0, ]
      warning(paste('Chunk',x,'missing spatial info and will be dropped.'))
    } else {

      dupRecords <- function(j) {
        record <- chunk[j, , drop = FALSE]       # pull out the record
        # drop columns also in "points"
        record <- record[, !(colnames(record) %in% colnames(points)), with=FALSE]
        # faster than rbind / replicate
        record_dup <- cbind(record, points)

        return(record_dup)
      }

      duped_records <- lapply(1:nrow(chunk), dupRecords)

      # append the duplicated data
      new_chunk <- rbindlist(duped_records)
    }
    return(new_chunk)
  }

  new_chunks <- mclapply(1:n_chunk,getNewChunksWrapper,mc.cores=cores)

  # remove old chunks from data
  message('Finishing up..\n')
  idx_all <- unlist(indices)
  data <- data[-idx_all, ]

  new_chunks_all <- rbindlist(new_chunks, fill = TRUE)
  new_chunks_all$pseudocluster = TRUE
  #if(length(data$start_year!=0)) {
  data$pseudocluster= FALSE
  data <- rbind(data, new_chunks_all, fill=TRUE)
 
  # count rows with no lat/longs and remove them
  no_ll <- which(is.na(data$latitude) | is.na(data$longitude))
  message(paste('Dropping',length(no_ll),'rows with missing spatial information.'))
  if(length(no_ll) != 0) data <- data[-no_ll, ]

  stopCluster(cl)

  data <- data.table(data)
  return(data)

}
