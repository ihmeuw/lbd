##########################################################################
## Functions for setting up an MBG model run
##########################################################################

## Make time stamp in standardized format.
  make_time_stamp <- function(time_stamp) {

    run_date <- gsub("-","_",Sys.time())
    run_date <- gsub(":","_",run_date)
    run_date <- gsub(" ","_",run_date)

    if(time_stamp==FALSE) run_date <- 'scratch'

    return(run_date)

  }
  
## Create directory structure
#   Arguments:
#     indicator_group = Category of indicator, i.e. "education"
#     indicator       = Specific outcome to be modeled within indicator category, i.e. "edu_0"
  create_dirs <- function(indicator_group, indicator) {

    dir.create(paste0('<<<< FILEPATH REDACTED >>>>'))
    dir.create(paste0('<<<< FILEPATH REDACTED >>>>'))

    indicator_dir <- paste0('<<<< FILEPATH REDACTED >>>>')

    for(dir in c('output','model_image_history')) {
      dir.create(paste0('<<<< FILEPATH REDACTED >>>>'), showWarnings = FALSE)
    }

  }

## Make template rasters (pull from central analysis folders managed by Lucas based off the area to model specified in config)
#   Arguments:
#     simple = Single polygon that defines boundaries of the entire area you want to model over.
#   Returns: Empty raster over modeling area. To be used for cropping covariates quickly and projecting model.
  get_template_raster <- function(simple) {

    message('Creating rasters of admin units')

    # Centrally controlled folder of analysis shapefiles
    analysis_raster_dir <- '<<<< FILEPATH REDACTED >>>>'

    # Load empty raster for largest analysis area, mask/crop to selected analysis area
    template_raster <- raster(paste0(analysis_raster_dir, 'stage2_analysis.tif'))
    template_raster <- mask(crop(template_raster,simple),simple)

    return(template_raster)

  }


## Load input data from required location
#   Arguments:
#     indicator = Specific outcome to be modeled within indicator category, i.e. "edu_0"
#     simple    = Single polygon that defines boundaries of the entire area you want to model over.
#   Returns: Input data subset to modeling area.
   load_input_data <- function(indicator, simple = NULL, agebin = 0, removeyemen = FALSE, pathaddin = "",
                               withdate=FALSE, date='', years='five_year',range=5, update_run_date = FALSE,
                               withtag=FALSE, datatag='', use_share=FALSE, yl = year_list) {

    # Ensure str_match loaded
    str_match <- stringr::str_match

    if(withdate){
      if(date=='')
        rd=run_date
      if(date!='')
        rd=date
    } else {
      rd = run_date
    }

    # Load input data by indicator
    if(use_share==FALSE) load_dir <- '<<<< FILEPATH REDACTED >>>>'
    if(use_share==TRUE) load_dir  <- '<<<< FILEPATH REDACTED >>>>'

    if(!withdate & !withtag) filename <- paste0(load_dir, indicator)
    if(withtag)              filename <- paste0(load_dir, indicator, datatag)
    if(withdate)             filename <- paste0('<<<< FILEPATH REDACTED >>>>', indicator)

    # try to see if an RDS exists, if so use that, if not use a csv
    if(file.exists(paste0(filename,'.RDS'))){
      message('READING INPUT DATA FROM RDS FILE')
      d <- readRDS(paste0(filename,'.RDS'))
    } else {
      message('READING INPUT DATA FROM CSV FILE')
      d <- read.csv(paste0(filename,'.csv'))
    }

    d$latitude  <- as.numeric(as.character(d$latitude))
    d$longitude <- as.numeric(as.character(d$longitude))
    message(nrow(d))
    d=d[d$latitude<=90,]
    d=d[d$latitude>=-90,]
    d=d[d$longitude<=180,]
    d=d[d$longitude>=-180,]
    d <- subset(d, !is.na(latitude))
    d <- subset(d, latitude!=0)
    message(nrow(d))

    # Check for necessary columns
    if(!(indicator %in% names(d))) stop(paste0("Your input data does not contain a column for your indicator: ", indicator))

    # Subset to within modeling area
    if(!is.null(simple)){
      coordinates(d) <- c("longitude", "latitude")
      proj4string(d) <- proj4string(simple)
      d$keep <- !is.na(over(d, as(simple, "SpatialPolygons")))
      message(paste0(round(mean(d$keep), 2)*100, '% of input data in specified template'))
      d <- d[d$keep==TRUE,]
    }
    d <- as.data.table(d)

    if(agebin!=0)   d = d[age%in%agebin,]
    if(removeyemen) d = d[country!='Yemen' & country!='YEM',]

    # remap any years as needed
    if(years=='five_year') {
      d <- d[year >= 1998 & year <= 2002, year := 2000]
      d <- d[year >= 2003 & year <= 2007, year := 2005]
      d <- d[year >= 2008 & year <= 2012, year := 2010]
      d <- d[year >= 2013 & year <= 2017, year := 2015]
    }

    if (nrow(subset(d, year < min(yl))) > 0) {
      warning(paste0("Dropping all data before min(year_list) = ", min(yl), "..."))
      d <- subset(d, year >= min(yl))
    }
    if (nrow(subset(d, year > max(yl))) > 0) {
      warning(paste0("Dropping all data after max(year_list) = ", max(yl), "..."))
      d <- subset(d, year <= max(yl))
    }

    # Change all "country" assignments to national level (in case subnational in the input data)
    if (nrow(d[grepl("[A-Z]*_[.]*", country),]) > 0) {
      subnat_countries <- unique(d[grepl("[A-Z]*_[.]*", country), country])
      warning(paste0("Changing subnational to national country codes for the following: ",
                     paste0(subnat_countries, collapse = ",")))
      d[grepl("[A-Z]*_[.]*", country), country := str_match(country,"([A-Z]*)_[.]*")[,2]]
    }

    # creaste a weighted SS to base QTs on
    if(sum(c('N','weight') %in% colnames(d)) == 2) d[,weighted_n := N*weight]

    # Save a copy
    if(update_run_date == TRUE) {
      if(dir.exists('<<<< FILEPATH REDACTED >>>>') == TRUE) {
        existing_dir <- '<<<< FILEPATH REDACTED >>>>'
        new_try <- existing_dir
        index <- 0
        while(dir.exists(new_try)) {
          index <- index + 1
          new_try <- paste0(existing_dir, '_', index)
        }
        run_date <- paste0(run_date, '_', index)
        dir.create(new_try, showWarnings = FALSE)
        run_date_dir <- new_try
      }
      if(dir.exists('<<<< FILEPATH REDACTED >>>>') == FALSE) {
        run_date_dir <- '<<<< FILEPATH REDACTED >>>>'
        dir.create(run_date_dir, showWarnings = FALSE)
      }
      write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
      return(list(d, run_date))
    }

    if(update_run_date == FALSE) {
      if(agebin==0){
        run_date_dir <- '<<<< FILEPATH REDACTED >>>>'
        dir.create(run_date_dir, showWarnings = FALSE)
        write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
      } else {
        run_date_dir <- '<<<< FILEPATH REDACTED >>>>'
        dir.create(run_date_dir, showWarnings = FALSE)
        write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
      }
      return(d)
    }

  }

## Read in shapefile and dissolve to single polygon for creating one big mesh (no admin1 boundaries)
#   gaul_list = any
#   buffer = distance to buffer around object
#   tolerance = tol in the gSimplify function. Higher # = coarser polygon
#   use_premade = shouold premade versions of regions / Africa be used?
#   raking = pull raking shapefile
#   shapefile_version = string specifying which shapefile version to pull

load_simple_polygon <- function(gaul_list, buffer, tolerance = 0.2,
                                subset_only = F, makeplots = F, use_premade = F,
                                custom_shapefile_path = NULL, custom_shapefile = NULL,
                                raking = F, shapefile_version = 'current') {

  # Logic check
  if (!is.null(custom_shapefile_path) & !is.null(custom_shapefile)) stop("You cannot specify both a custom shapefile and a custom shapefile path")

  # If using custom shapefiles, don't use premade polys
  if (!is.null(custom_shapefile_path) | !is.null(custom_shapefile)) use_premade <- F

  if (use_premade == T) {
    # Check for common gaul lists so we can query premade simple polygons
    if(identical(gaul_list, c(29, 42, 45, 47, 50, 66, 90, 94, 106, 105, 144, 155, 159, 181, 182, 214, 217, 221, 243))) {
      message("Your GAUL list matches WSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>/wssa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(8,49,59,68,76,89))) {
      message("Your GAUL list matches CSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>/cssa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(43,58,70,77,79,133,150,152,170,205,226,74,257,253,270))) {
      message("Your GAUL list matches ESSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>/essa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(4,40762,40765,145,169,6,248))) {
      message("Your GAUL list matches NAME, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>/name.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(35,142,172,227,235,271))) {
      message("Your GAUL list matches SSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>/sssa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(4,6,8,29,35,42,43,45,47,49,50,58,59,66,68,70,
                              74,76,77,79,89,90,94,95,105,106,142,144,145,150,
                              152,155,159,169,170,172,181,182,205,214,217,221,
                              226,235,243,248,253,268,270,271,40762,40765,
                              227,257,133,269))) {
      message("Your GAUL list matches AFRICA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>/africa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
  }

  # Otherwise, make a new simple_poly
  # count the vertices of a SpatialPolygons object, with one feature
  vertices <- function(x) sum(sapply(x@polygons[[1]]@Polygons, function(y) nrow(y@coords)))

  # ~~~~~~~~~~~~~~~~~
  # load data

  if (is.null(custom_shapefile_path) & is.null(custom_shapefile)) {
    message("Opening master shapefile...")
    master_shape <- readOGR(get_admin_shapefile(admin_level = 0, raking = raking,
                                                version = shapefile_version))
    master_shape@data$ADM0_CODE <- as.numeric(as.character(master_shape@data$ADM0_CODE))
    subset_shape <- master_shape[master_shape@data$ADM0_CODE %in% gaul_list, ]
  } else if (!is.null(custom_shapefile_path) & is.null(custom_shapefile)) {
    message("Opening custom shapefile...")
    master_shape <- readOGR(custom_shapefile_path)
    subset_shape <- master_shape
  } else if (is.null(custom_shapefile_path) & !is.null(custom_shapefile)) {
    master_shape <- custom_shapefile
    subset_shape <- master_shape
  }

  if(subset_only==TRUE) {
    return(list(subset_shape=subset_shape,spoly_spdf=NULL))
  }

  if(subset_only==FALSE) {

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    message('Making a super-low-complexity map outline for INLA mesh creation')

    message(paste0('Full polygon vertices: ', vertices(subset_shape)))

    # Merge everything together (where intersecting)
    af <- gUnaryUnion(subset_shape)

    # Initial simplification
    af_simple <- gSimplify(af, tol = tolerance, topologyPreserve = TRUE)

    # Remove tiny features
    # Get all sub polygons and their areas
    polys <- af_simple@polygons[[1]]@Polygons
    areas <- sapply(polys, function(x) x@area)

    # If there is more than one sub polygon, remove the ditzels (many single-country subsets are a single polygon,
    #   like Uganda, which would break these few lines)
    if(length(areas)>1) {
      # find top 5% by area
      big_idx <- which(areas > quantile(areas, 0.95))

      # convert back into a spatialPolygons object
      spoly <- SpatialPolygons(list(Polygons(polys[big_idx], ID = 1)))
    }
    if(length(areas)==1) {
      spoly <- af_simple
      big_idx <- 1
    }

    # Buffer slightly
    spoly <- gBuffer(spoly, width = buffer)

    # simplify again to reduce vertex count
    spoly <- gSimplify(spoly, tol = tolerance, topologyPreserve = TRUE)

    # Get list of original polygons
    polys2 <- af@polygons[[1]]@Polygons

    # Check if all are within the simple polygon
    check_if_in_spoly <- function(the_poly, compare_to, the_proj = projection(master_shape)) {
      the_poly <- SpatialPolygons(list(Polygons(list(the_poly), ID = 1)))
      projection(the_poly) <- the_proj
      projection(compare_to) <- the_proj

      if(suppressWarnings(gIsValid(the_poly)) == F) return(TRUE) #Ignore invalid polygons

      poly_intersect <- rgeos::gIntersection(the_poly, compare_to)

      if(is.null(poly_intersect)) {
        return(FALSE)
      } else {
        return(ifelse((raster::area(poly_intersect) == raster::area(the_poly)), TRUE, FALSE))
      }
    }

    over_list <- sapply(polys2, function(x) check_if_in_spoly(x, compare_to = spoly))

    if (all(over_list) == FALSE) {
      # Add back in polygons if missed by above procedure (e.g. islands dropped)
      big_idx <- unique(c(big_idx, which(over_list == F)))
      spoly <- SpatialPolygons(list(Polygons(polys[big_idx], ID = 1)))
      spoly <- gBuffer(spoly, width = buffer)
      spoly <- gSimplify(spoly, tol = tolerance, topologyPreserve = TRUE)
    }

    # Now check again with new spoly
    over_list <- sapply(polys2, function(x) check_if_in_spoly(x, compare_to = spoly))

    # If still not all enclosed, tolerance probably too high. Return warning
    if (all(over_list) == FALSE) {
      number_false = length(over_list[over_list == F])
      number_total = length(over_list)
      warning(paste0(number_false, " of ", number_total, " polygons are NOT enclosed within your simple polygon. \n",
                     "Adjust your buffer and tolerance values."))
    }

    # Return results
    message(paste0('Simplified vertices: ', vertices(spoly)))

    # plot to check it encloses all of the important bits
    if(makeplots) plot(spoly)
    if(makeplots) plot(af, add = TRUE, border = grey(0.5))

    # turn into an SPDF
    spoly_spdf <- SpatialPolygonsDataFrame(spoly,
                                           data = data.frame(ID = 1),
                                           match.ID = FALSE)

    # add projection information
    projection(spoly_spdf) <- projection(master_shape)

    return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
  }
}


## build_space_mesh() ------------------------------------------------------------
#'
#' @title Buils spatial FEM mesh
#'
#' @description Build a finite-elements triangulation mesh to use in
#'   SPDE GP approximation
#'
#' @param d data.frame or data.table of observational data with
#'   'longitude' and 'latitude' columns
#'
#' @param simple simple_polygon defining modeling domain
#'
#' @param max_edge string of 2 element numeric vector in R
#'   notation. e.g. "c(0.25, 5)". First entry decribes the maximum
#'   allowed triangle edge length within the simple_polygon modeling
#'   boundary. Second entry describes the maximum allowed triaqngle
#'   edge length in the extended buffer region. Units are in
#'   lat-long distances. Used if s2mesh=FALSE.
#'
#' @param mesh_offset string of 2 element numeric vector in R
#'   notation. e.g. "c(1, 5)". Describes the automatic extension
#'   distance from the 'simple' boundary. The entries are the inner
#'   and outer extension distances, respectively. Units are in
#'   lat-long distance. Used if s2mesh=FALSE.
#'
#' @param plot_mesh Logical. Should a plot of the mesh be generated?
#'   Not currently implemented for s2mesh
#'
#' @param s2mesh Logical. Should the mesh be created on the surface of
#'   a sphere? If TRUE, s2params is used to specify mesh parameters
#'   instead of max_edge and mesh_offset
#'
#' @param s2params string of 3 element numeric vector in R
#'   notation. e.g. "c(25, 500, 1000)". The entries describe the
#'   minimum triangle edge length allowed, hos far to extend the mesh
#'   beyond the 'simple' boundary, and the maximum allowed triangle
#'   edge length, respectively. Units are in kilometers. Used only if
#'   s2mesh=TRUE.
#'
#' @return an 'inla.mesh' object
#'

build_space_mesh <- function(d, simple, max_edge, mesh_offset,
                             plot_mesh = F, s2mesh = FALSE,
                             s2params = NULL){

  if(!as.logical(s2mesh)){ ## build mesh on R2 plane
    message(paste0('Creating spatial mesh, max edge parameter: ', max_edge))
    max.edge <- eval(parse(text=max_edge))
    mesh_offset <- eval(parse(text=mesh_offset))
    mesh_s <- inla.mesh.2d(
      boundary = inla.sp2segment(simple),
      loc = cbind(d$longitude,d$latitude),
      max.edge = max.edge,
      offset = mesh_offset,
      cutoff = max.edge[1]
    )

    if (plot_mesh) {
      plot(mesh_s, asp=1)
      points(d$longitude, d$latitude, col=d$year)
    }

  } else{ ## build mesh on sphere surface

    if(is.null(s2params)){
      stop("You've chosen to build an s2 mesh but haven't provided parameters to do so (i.e. s2params=NULL)!")
    }

    message(paste0('Creating spatial SPHERICAL mesh, min edge, max edge, extension kms are: ', s2params))
    s2params <- eval(parse(text = s2params))

    ## convert data locs to 3d coords on the unit-sphere
    true.radius.of.earth = 6371
    s3 <- lonlat3D(d$longitude, d$latitude)

    ## convert boundary of simple_polygon to 3d coords on unit-sphere
    boundary <- inla.sp2segment(simple)
    boundary.loc <- lonlat3D(boundary$loc[, 1], boundary$loc[, 2])

    ## make a s2 domain mesh using data locs and boundary locs
    all.loc <- rbind(s3, boundary.loc)
    mesh_s <- inla.mesh.create(loc=all.loc,
                               cutoff = s2params[1] / true.radius.of.earth, ## minimum triangle edge allowed
                               extend = list(offset = s2params[2] / true.radius.of.earth), ## how far to extend mesh
                               refine=list(max.edge = s2params[3] / true.radius.of.earth)) ## max triangle edge allowed

    ## TODO add 3d mesh plotting
    ## an example of how to do this can be found in the code examples here:
    ## http://www.r-inla.org/examples/case-studies/simpson2011
    if(plot_mesh){
      ## plot 3d mesh
    }
  }

  return(mesh_s)

}

## Create temporal mesh (defaulting to the four period U5M approach for now, come back and make more flexible later)
build_time_mesh <- function(periods=1:4) {
  mesh_t <- inla.mesh.1d(
    loc = c(periods),
    degree = 1,
    boundary = rep("free", 2)
  )
  return(mesh_t)
}


## Alternative method of rasterizing:
##### We go from shapefile to SpatialPointsDataFrame, which will give us
#### an SPDF of boundaries. Then we merge this on top of the result that
#### comes out of rasterize(), and that will fix the admin border issues

#'
#' @title Shapefile to SpatialPointsDataFrame
#'
#' @description Converts a shapefile into a SpatialPointsDataFrame for a given field
#'
#' @param sp.df SpatialPolygonDataFrame. Input admin shapefile
#' @param field String. The field describing the shapes, for e.g. ADM0_CODE
#'
#' @return A SpatialPointsDataFrame with coords and value of field
#'
#' @importFrom sp SpatialPointsDataFrame
#' @export
#'
#'
shape_to_spointsdf <- function(sp.df, field) {
  
  ## Loop over 'field' units
  each_country <- lapply(c(1:length(sp.df@polygons)), function(i) {
    ## Loop over polygons for field 'i'
    rbindlist(lapply(c(1:length(sp.df@polygons[[i]]@Polygons)), function(j) {
      data.table(sp.df@polygons[[i]]@Polygons[[j]]@coords, V3 = sp.df@data[[field]][i])
    }))  
  })
  
  
  results <- data.frame(rbindlist(each_country))
  results <- sp::SpatialPointsDataFrame(coords=results[,1:2], data=data.frame(field=results[,3]))
  
  
  ## Fix field name to the string supplied
  names(results) <- paste0(field)

  return(results)
}


#' @title Rasterize with border checks
#'
#' @description Rasterizing using a shapefile and a template raster, such that
#' we account for any pixels that are on the border of \code{field} units, which
#' could have been lost due to raster::rasterize only evaluating on centroids
#'
#' @param shapes SpatialPolygonDataFrame.. Input shapefile
#'
#' @param template_raster SpatialPolygonDataFrame.. The reference raster (usually WorldPop)
#'
#' @param field String The field with appropriate administrative unit information (usually ADM0_CODE)
#'
#' @param link_table String or data.table. If data.table it is used as-is. If String: either an absolute
#'   file path to an RDS file OR a short name for the administrative shape file e.g., "2019_02_27" or "current".
#'
#' @return A raster with border and insides properly filled
#'
#' @details rasterize_check_coverage has three distinct use cases based off of the value of link_table
#'
#' 1. \code{link_table} is NULL. In this case rasterize_check_coverage will behave identically to raster::rasterize
#'
#' 2. \code{link_table} is a String referring to relase of admin shapefiles ("current" or e.g., "2019_02_27"). In this case
#'    \code{field} should be "ADM0_CODE", "ADM1_CODE" or "ADM2_CODE". This will load the lbd_standard_link.rds file,
#'    from the related admin shapefile directory, aggregate area_fraction as necessary to match the level of \code{field},
#'    and then apply those values to pixels in the space defined by \code{shapes}.
#'
#' 3. \link{link_table} is a data.table OR a String absolute path to a RDS file containing a data.table. This will use the
#'    provided \code{link_table} to assign values to the result raster similarly to use case #2.
#'
#' Note that for both use cases 2 and 3 all pixel_id coordinates must be in the same raster space. This is currently the
#' area defined by cropping the world raster to the pixels occupied by stage 1 and stage 2 countries.
#'
#' @export
#'
rasterize_check_coverage <- function(shapes, template_raster, field, ..., link_table = modeling_shapefile_version) {
  # backwards-compatible behavior - just call rasterize()
  if (is.null(link_table)) return(raster::rasterize(shapes, template_raster, field = field, ...))

  # Validate arguments
  is_admin_link_table <- FALSE
  if (is.data.table(link_table)) {
    is_admin_link_table <- TRUE
    # nothing to do - already a link table loaded in memory
  } else if (R.utils::isAbsolutePath(link_table)) {
    link_table <- readRDS(link_table)
  } else if (is_admin_shapefile_string(link_table)) {
    is_admin_link_table <- TRUE
    # load link table with pre-computed ownership percentages for each pixel cell
    link_table_file <- paste0(get_admin_shape_dir(link_table), "lbd_standard_link.rds")
    link_table <- readRDS(link_table_file)
  } else {
    stop("link_table argument was neither a data.table, an admin shapefile string, or an absolute path to a RDS file.")
  }

  if (! field %in% names(link_table)) {
    msg <- paste("WARNING: rasterize_check_coverage called with field", field,
                 "which is not present in link_table. Defaulting to raster::rasterize()")
    message(msg)
    return(raster::rasterize(shapes, template_raster, field = field, ...))
  }

  # aggregate link table generically for admin 0/1/2
  # Note: we need `with=FALSE` because `field` is passed as a parameter (not a hard-coded string)
  table <- link_table[,c("pixel_id", field, "area_fraction"), with=FALSE]
  if (is_admin_link_table && field != "ADM2_CODE") {
    # sum rows; area_fraction now represents the total area coverage by ADM0/1_CODE instead of ADM2_CODE
    table <- table[, .(area_fraction = sum(area_fraction)), by = c("pixel_id", field)]
  }
  # subset table so that we have 1 entry per pixel_id - the value of `field` with the maximum
  # area_fraction value for that pixel_id
  # https://stackoverflow.com/a/24558696
  pixel_owner <- table[table[, .I[which.max(area_fraction)], by=pixel_id]$V1]
  pixel_owner <- pixel_owner[order(pixel_id)]
  
  # generate world raster with pixel values for `field`
  world_pixel_owner <- suppressWarnings(empty_world_raster())
  # subset to only those pixels owned by a shape we're interested in
  owned_pixels <- pixel_owner[pixel_owner[[field]] %in% shapes[[field]]]
  world_pixel_owner[owned_pixels$pixel_id] <- owned_pixels[[field]]

  result <- raster::crop(world_pixel_owner, template_raster, snap="near")
  if (raster::ncell(result) != raster::ncell(template_raster)) {
    message <- paste("Error in creating result raster. Should have created a raster of shape",
                     paste(dim(result), collapse=","),
                     "but instead created a raster of shape",
                     paste(dim(template_raster), collapse=","))
    stop(message)
  }
  return(result)
}

#' @title Create raster representing world
#'
#' @description Creates a raster file consistent with LBD's world raster.
#'
#' @param whole_world Logical if TRUE return raster for entire world. If FALSE,
#'   return raster for stage1+stage2
#'
empty_world_raster <- function(whole_world = FALSE) {
  if (whole_world) {
    result <- raster::raster(nrow=4320, ncol=8640,
                             crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                             xmn=-180, xmx=180, ymn=-90, ymx=90,
                             vals=NA)
  } else {  # stage 1+2
    result <- raster::raster(nrow=2123, ncol=6610,
                             crs=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"),
                             xmn=-118.375, xmx=157.0417, ymn=-34.875, ymx=53.58333,
                             vals=NA)
  }
  return(result)
}


## Load list of raster inputs (pop and simple)
build_simple_raster_pop <- function(subset_shape, u5m=FALSE, field=NULL, raking=F, link_table=modeling_shapefile_version) {

  if (is.null(field)) {
    if ('GAUL_CODE' %in% names(subset_shape@data)) field <- 'GAUL_CODE'
    if ('ADM0_CODE' %in% names(subset_shape@data)) field <- 'ADM0_CODE'
  }

  if(raking){
    field <- 'loc_id'
    # no 'loc_id' field in the link table, so we can't use it
    link_table <- NULL
  }

  if(u5m==FALSE){
    master_pop <- brick('<<<< FILEPATH REDACTED >>>>/WorldPop_total_global_stack.tif')
  } else {
    master_pop <- brick(raster('/<<<< FILEPATH REDACTED >>>>/worldpop_a0004t_5y_2000_00_00.tif'),
    raster('/<<<< FILEPATH REDACTED >>>>/worldpop_a0004t_5y_2005_00_00.tif'),
    raster('/<<<< FILEPATH REDACTED >>>>/worldpop_a0004t_5y_2010_00_00.tif'),
    raster('/<<<< FILEPATH REDACTED >>>>/worldpop_a0004t_5y_2015_00_00.tif'))
  }

  cropped_pop <- crop(master_pop, extent(subset_shape), snap="out")
  ## Fix rasterize
  initial_raster <- rasterize_check_coverage(subset_shape, cropped_pop, field = field, link_table = link_table)
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))!=0) {
    rasterized_shape <- 
      raster::merge(
        rasterize_check_coverage(subset(subset_shape, !(get(field) %in% unique(initial_raster))),
                                 cropped_pop,
                                 field = field,
                                 link_table = link_table),
        initial_raster)
  }
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))==0) {
    rasterized_shape <- initial_raster
  }
  masked_pop <- raster::mask(x=cropped_pop, mask=rasterized_shape)

  raster_list <- list()
  raster_list[['simple_raster']] <- rasterized_shape
  raster_list[['pop_raster']] <- masked_pop

  return(raster_list)

}

## #############################################################################
## GET ADMIN CODES AND RELATED FUNCTIONS
## #############################################################################


## load_adm0_lookup_table() ----------------------------------------------------
#'
#' @title Load the GAUL lookup table
#'
#' @description Loads the most recent version of the lookup table that links
#'   ADM0 codes with other identifiers such as GBD location IDs, MBG modeling
#'   regions, and ISO codes, among others.
#'
#' @return Returns a data.frame of the ADM0 lookup table. If package data.table
#'   is loaded, also converts the lookup table to a data.table
#'
load_adm0_lookup_table <- function(){
  # Define consistent path to the ADM0 lookup table
  lookup_table_filepath <- paste0('/<<<< FILEPATH REDACTED >>>>/',
                                  'stage_master_list.csv')
  # If data.table is loaded, read in as a data.table
  if ('package:data.table' %in% search()){
    lookup_table <- fread(lookup_table_filepath)
  } else {
     # Otherwise, read in as a data.frame
     lookup_table <- read.csv(lookup_table_filepath)
  }
  # Set ISO codes as lowercase for easy lookup
  lookup_table$iso3 <- tolower(lookup_table$iso3)
  return (lookup_table)
}



## pull_custom_modeling_regions() ----------------------------------------------
#'
#' @title Pull custom modeling regions
#'
#' @description Define modeling regions that are not simple combinations of the
#'   default MBG regions (in other words, regions that are not combinations of
#'   four-letter MBG regions such as "wssa" or "seas+ocea" or ISO-3 codes such
#'   as 'ZAF' or 'CHN').
#'
#' @param custom_region character vector of custom modeling regions
#'
#' @return Returns a named list of custom modeling regions with associated
#'    "standard" (non-custom) modeling regions that can be directly interpreted
#'    by get_adm0_codes().
#'
pull_custom_modeling_regions <- function(custom_regions){
  custom_regions <- tolower(custom_regions)

  # FULL LIST OF ALL REFERENCE REGIONS
  # If you need to add a new custom region, add it to this list
  ref_reg_list <- list(
    'africa'          = 'noaf+essa+wssa+cssa+sssa-yem',
    'middle_east'     = 'mide+stan-pak',
    'eastern_europe'  = "blr+est+lva+ltu+mda+ukr",
    'latin_america'   = 'caca+trsa+ansa',
    'south_asia'      = 'soas+chn_d2+pak',
    'central_america' = 'caca',
    'south_america'   = 'ansa+trsa',
    # se_asia was historically inclusive of East Asia + SE Asia
    'se_asia'         = 'eaas+seas+ocea+png',

    #data coverage regions
    'africa_dcp' = 'noaf+essa+wssa+cssa+sssa+yem',
    'middle_east_dcp' = 'mide+stan-yem-pak',
    'latin_america_dcp' = 'latin_america+cub',
    'south_asia_dcp' = 'south_asia-mdv-syc',
    'se_asia_dcp' = 'eaas+seas+png+idn+phl+tls+mys+twn',

    'stage1' = 'noaf+essa+wssa+cssa+sssa-yem',
    # ONLY stage 2 countries (not inclusive of Stage 1)
    'stage2' = 'ansa+caca+stan+eaas+mide+ocea+soas+seas+trsa+yem',
    'stage3' = 'all-stage1-stage2',

    # India + all disputed territories it claims: Jammu & Kashmir, Arunachal
    #   Pradesh, Aksay Chin, Demjok, Tirpani, Bara Hotii, & Samdu
    'ind_all_territories' = 'ind+ind_d1+ind_d2+ind_d3+chn_d2',

    # Getting into modeler-defined regions
    'cssa_diarrhea'  = 'cssa-ago',
    'sssa_diarrhea'  = 'sssa+ago-swz',
    'essa_diarrhea'  = 'essa+swz-ken_d1-yem',
    'cssa_diarrhea2' = 'cssa+cmr+tcd-ago-cod',
    'essa_diarrhea2' = 'essa+sdn+swz-ken_d1-yem',
    'name_diarrhea2' = 'noaf-esh-egy-egy_d1-sdn-sdn_d2',
    'sssa_diarrhea2' = 'sssa+ago-swz',
    'wssa_diarrhea2' = 'wssa+cmr+tcd-cmr-tcd',

    'essa_edu' = 'dji+eri+eth+som+ssd+sdn',
    'cssa_edu' = 'cssa+essa+nga+swz+cmr-dji-eri-eth-som-ssd-sdn-ken_d1-yem',
    'name_edu' = 'noaf-egy-egy_d1-sdn-sdn_d2-esh',
    'sssa_edu' = 'sssa-swz',
    'wssa_edu' = 'wssa-cmr-nga',

    'essa_sdn' = 'essa+sdn-ken_d1-yem',

    'cessa' = 'cssa+essa-ken_d1-yem',
    'cwssa' = 'cssa+wssa',

    'namelite' = 'noaf-dza-lby-tun-sdn_d1-sdn_d2-egy_d1-esh',

    # Region 'cessa2' originally contained American Samoa--this has been dropped
    'cessa2' = 'cssa+essa-ago-zmb-mwi-moz-ken_d1-yem',
    'sssa2'  = 'sssa+ago+zmb+mwi+moz',

    'cssa_cam'   = 'cssa+cmr',
    'wssa_nocam' = 'wssa-cmr',

    'name_hi' = 'noaf-sdn-sdn_d2-egy_d1-esh',
    'essa_hi' = 'zmb+ken',
    'essa_lo' = 'essa+sdn-zmb-ken-ken_d1-yem',
    'cssa_hi' = 'gnq+cog+gab',
    'cssa_lo' = 'cssa-gnq-cog-gab',
    'wssa_hi' = 'gha',
    'wssa_lo' = 'wssa-gha-mrt',
    'sssa_hi' = 'bwa+nam+zaf',
    'sssa_lo' = 'sssa-bwa-nam-zaf',

    'essa_hilo' = 'essa+lso+sdn+swz+zwe-ken_d1-yem',

    # Vaccines
    'vax_soas' = 'soas+pak-syc',
    'vax_seas' = 'seas+ocea-asm-fji-kir-wsm-ton',
    'vax_eaas' = 'eaas',

    # Need to move these back to vax_seas when 0long issue fixed
    # Not using this for now
    'vax_seas_0long' = 'asm+fji+kir+wsm+ton',

    'vax_caeu' = 'arm+aze+geo+kgz+mda+tjk+tkm+ukr+uzb',

    'vax_crbn' = 'cub+dma+dom+grd+hti+jam+lca+vct+vir',
    'vax_ctam' = 'blz+col+cri+slv+gtm+hnd+mex+nic+pan+ven',
    'vax_ansa' = 'ansa-col-ven',
    'vax_trsa' = 'trsa',

    'vax_name' = 'noaf+mide+afg+omn',
    'vax_cssa' = 'cssa',
    'vax_essa' = 'essa+syc',
    'vax_sssa' = 'sssa',
    'vax_wssa' = 'wssa',

    # For modelers who are still using the defunct NAME region
    'name_historic' = 'noaf-esh',

    # Custom regions for Diarrhea/ORT/WASH/LRI/HAP
    'dia_afr_horn' = 'dji+eri+eth+sdn+som+ssd+yem',
    'dia_cssa' = 'ago+caf+cod+cog+gab+gnq+stp',
    'dia_wssa' = 'ben+bfa+civ+cmr+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+ner+nga+sen+sle+tcd+tgo',
    'dia_name' = 'dza+egy+lby+mar+tun',
    'dia_sssa' = 'bwa+nam+zaf',
    'dia_mcaca' = 'blz+cri+cub+dma+dom+grd+gtm+hnd+hti+jam+lca+mex+nic+pan+slv+vct',
    'dia_s_america' = 'bol+bra+col+ecu+guy+per+pry+sur+tto+ven',
    'dia_central_asia' = 'kgz+tjk+tkm+uzb',
    'dia_chn_mng' = 'chn+mng',
    'dia_se_asia' = 'khm+lao+mmr+mys+tha+vnm',
    'dia_malay' = 'idn+phl+png+tls',
    'dia_south_asia' = 'bgd+btn+ind+lka+npl+pak',
    'dia_mid_east' = 'afg+irn+irq+jor+pse+syr',
    'dia_essa' = 'bdi+com+ken+lso+mdg+moz+mwi+rwa+swz+syc+tza+uga+zmb+zwe',
    'dia_oceania' = 'asm+fji+fsm+kir+mhl+slb+ton+vut+wsm',
    
    # Additional custom regions for ORT
    'dia_s_america_n' = 'guy+col+sur+tto',
    'dia_s_america_s' = 'bol+ecu+per+pry',
    'dia_wssa_e' = 'cmr+ner+tcd',
    'dia_wssa_w' = 'ben+bfa+civ+cpv+gha+gin+gmb+gnb+lbr+mli+mrt+sen+sle+tgo'

  )
  # Warn if there are any custom regions not in the reference list
  missing_regions <- custom_regions[ !(custom_regions %in% names(ref_reg_list)) ]
  if( length(missing_regions) > 0 ){
    message(paste0('WARNING: The following custom regions are not defined: ',
                   paste(missing_regions, collapse=','))
            )
  }
  # Return a named list of all custom regions
  custom_regions_list <- ref_reg_list[ names(ref_reg_list) %in% custom_regions ]
  return(custom_regions_list)
}



## get_adm0_codes() ------------------------------------------------------------
#'
#' @title Get Admin0 (Country) Codes
#'
#' @description Pull Admin0 Codes for the specified countries or regions,
#'   optionally excluding certain Admin0 codes.
#'
#' @param adm0s Vector of ISO-3 codes, four-letter modeling region names, region
#'   group names as defined in get_region_groups, or 'all' (returns all Admin0
#'   codes)
#' @param strict (default `FALSE`) Causes the function to fail if an ISO code or
#'   region name is not included in the full list.
#' @param lookup_table (default `NULL`) Sets a data.frame or data.table to be
#'   used as the lookup_table. If `NULL`, the lookup table will be loaded using
#'   the `load_adm0_lookup_table()` function.
#' @param core_repo (default NULL) THIS ARGUMENT IS DEPRECATED AND WILL BE
#'   REMOVED IN THE FUTURE. Please remove it from your function calls.
#' @param adm0_type (default 'detect') Which class of admin0 codes
#'   should be pulled? Must be one of 'gaul', 'gadm', or 'detect'. If
#'   'gaul' or 'gadm', it will simply use what was specified. If
#'   'detect', the function will detect the adm0_type directly from an
#'   admin_shapefile. It will check the adm0_type specified by the
#'   shapefile_version argument if it is specified, otherwise if will
#'   detect the adm0_type using 'current'.
#' @param shapefile_version string specifying shapefile_version to be
#'   used in adm0_type determination if adm0_type is set to 'detect'
#' @param subnational_raking Logical. set to true if you want to use raking version of shapefile
#' @return Vector of numeric ADM0 codes. Warns if not all regions align with a
#'   ADM0 code.
#'
get_adm0_codes <- function(adm0s,
                           strict = FALSE,
                           lookup_table = NULL,
                           core_repo = NULL,
                           adm0_type = 'detect',
                           shapefile_version = NULL,
                           subnational_raking = FALSE
                           ){

  # Check adm0_type input
  if(!(adm0_type %in% c('gaul', 'gadm', 'detect'))){
    stop("You must select either 'gaul', 'gadm', or 'detect' for adm0_type!")
  }
  # Check that data.table has been loaded
  if(!('data.table' %in% .packages())){
    stop("Please load the data.table package before running get_adm0_codes!")
  }
  # core_repo deprecation message
  if( !is.null(core_repo) ){
    message(paste0(
      "WARNING: The 'core_repo' argument has been deprecated in get_adm0_codes()",
      " -- please remove this argument."
    ))
  }

  # Determine adm0_type
  if(adm0_type == 'detect'){
    if(is.null(shapefile_version)){
      adm0_type <- detect_adm_shapefile_date_type(
        get_admin_shapefile(version = 'current', raking = subnational_raking)
      )$shpfile_type
    } else{
      adm0_type <- detect_adm_shapefile_date_type(
        get_admin_shapefile(version = shapefile_version, raking = subnational_raking)
      )$shpfile_type
    }
  }
  message(sprintf('Using ADM codes from %s', adm0_type))

  # Read in the ADM0 code data.table once at the beginning of the function
  if(is.null(lookup_table)) lookup_table <- load_adm0_lookup_table()

  # Check that a valid adm0_type has been defined
  # If so, define the field that will be pulled from the lookup table
  adm_type <- tolower(adm0_type)
  adm0_field_reference = list(
    gadm = 'gadm_geoid',
    gaul = 'GAUL_CODE'
  )
  if(!(adm0_type %in% names(adm0_field_reference))){
    stop(paste0(
      "Invalid adm0_type value. Valid values are: ",
      paste(names(adm0_field_reference,collapse=', '))
    ))
  }
  pull_field <- adm0_field_reference[[adm0_type]]

  ## PROCESS INCOMING ADM0 TEXT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # In case multiple character strings were passed in, concatenate with '+'
  adm0s <- paste(adm0s, collapse='+')
  # Set all as lowercase
  adm0s <- tolower(adm0s)

  # Split on - and + characters
  adm0s <- gsub(';','',adm0s)
  adm0s <- gsub('\\+',';;\\+',adm0s)
  adm0s <- gsub('-',';;-',adm0s)
  adm0s_vec <- unlist(base::strsplit(adm0s, ';;'))

  # Making the function past-compatible, for modelers who are still using the
  #  original 'name' modeling region for North Africa (now 'noaf').
  #  'name_historic' is now a custom region, equivalent to 'noaf-esh', in the
  #  pull_custom_modeling_regions().
  adm0s_vec <- gsub('^([+-])?name$', '\\1name_historic', adm0s_vec)

  # If the string was empty, end now
  if (length(adm0s_vec)==0){
    message("No GAUL codes were returned; exiting")
    return(numeric(0))
  }

  # The first region and all regions beginning in '+' are included
  # Strip the beginning + signs
  include <- gsub('\\+','', c(adm0s_vec[1], adm0s_vec[ grepl('^\\+',adm0s_vec) ]) )
  # All regions beginning with '-' are treated as exclusions
  # Strip any minus signs from them
  exclude <- gsub('-','', adm0s_vec[ grepl('^-',adm0s_vec) ] )

  ## CHECK FOR CUSTOM REGIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define regex for all non-custom regions
  standard_code_regex <- "^([a-z]{3,4}|[a-z]{3}_d[0-9]{1,5})$"
  # Split both included and excluded data into custom and non-custom regions
  standard_inc <- include[ grepl(standard_code_regex, include) ]
  custom_inc   <- include[ !grepl(standard_code_regex, include) ]
  # Check in excluded data
  standard_exc <- exclude[ grepl(standard_code_regex, exclude) ]
  custom_exc   <- exclude[ !grepl(standard_code_regex, exclude) ]

  # Get definitions and ADM0 codes for included and excluded GAULs, if needed
  # `custom_inc_codes` and `custom_exc_codes` will be added/subtracted at the end
  # Check out that functional recursion
  recursive_pull <- function(x) get_adm0_codes(
    x, lookup_table=lookup_table, adm0_type=adm0_type
  )

  if ( length(custom_inc) > 0 ){
    custom_inc_defs  <- pull_custom_modeling_regions(custom_inc)
    custom_inc_codes <- unlist(lapply(custom_inc_defs, recursive_pull))
  } else {
    custom_inc_codes <- integer(0)
  }
  if ( length(custom_exc) > 0 ){
    custom_exc_defs  <- pull_custom_modeling_regions(custom_exc)
    custom_exc_codes <- unlist(lapply(custom_exc_defs, recursive_pull))
  } else {
    custom_exc_codes <- integer(0)
  }

  ## PROCESS NON-CUSTOM REGIONS (MBG MODELING REGIONS + ISOS) ~~~~~~~~~~~~~~~~~~
  # Helper function to pull all GAUL codes for a set of modeling regions or
  #  ISO codes
  utility_pull_adm0_codes <- function(region_codes, standard_regex){
    # Return none if the length of the inputs is 0
    if (length(region_codes)==0) return(integer(0))
    # Input data assertions
    if( !all(grepl(standard_regex, region_codes)) ){
      stop('All region codes must match the MBG or ISO code formats.')
    }
    isos     <- region_codes[nchar(region_codes)!=4]
    mbg_regs <- region_codes[nchar(region_codes)==4]
    if ('all' %in% isos){
      # 'all' is a special case where all GAUL codes are pulled
      pulled_adm0s <- unique(lookup_table[get(pull_field)>=0, get(pull_field)])
    } else (
      # Pull all GAUL codes matching the ISOs or regions
      pulled_adm0s <- unique(
        lookup_table[((iso3 %in% isos)|(mbg_reg %in% mbg_regs)) & (get(pull_field)>=0),
                     get(pull_field)]
      )
    )
    # Warn if any iso codes or MBG regions aren't in the lookup table
    missing_isos <- isos[ !(isos %in% lookup_table[,iso3]) & (isos!='all') ]
    if (length(missing_isos) > 0) message(paste0("WARNING: Missing these ISOs: ",
                                                 paste(missing_isos,collapse=',')))
    missing_regs <- mbg_regs[ !(mbg_regs %in% lookup_table[,mbg_reg]) ]
    if (length(missing_regs) > 0) message(paste0("WARNING: Missing these MBG regions: ",
                                                 paste(missing_regs,collapse=',')))
    return(as.integer(pulled_adm0s))
  }

  # Pull GAUL codes for standard include and exclude regions
  standard_inc_codes <- utility_pull_adm0_codes(standard_inc, standard_code_regex)
  standard_exc_codes <- utility_pull_adm0_codes(standard_exc, standard_code_regex)

  ## COMBINE TO CREATE FINAL ADM0 CODE SET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get all (standard + custom) codes to include
  all_include_codes <- unique( c(standard_inc_codes, custom_inc_codes) )
  # Get all (standard + custom) codes to exclude
  all_exclude_codes <- unique( c(standard_exc_codes, custom_exc_codes) )
  # Final codes = codes to include - codes to exclude
  return_codes <- all_include_codes[ !(all_include_codes %in% all_exclude_codes) ]

  return(return_codes)
}


## get_gaul_codes() ------------------------------------------------------------
#'
#' @title Get GAUL Codes (TO BE DEPRECATED)
#'
#' @description Pull Admin0 codes for the specified countries or regions,
#'   optionally excluding certain Admin0 codes. The only difference from the
#'   standard get_adm0_codes() function is that the default adm0_type is 'gaul'
#'   rather than 'gadm'.
#'
#' NOTE: THIS FUNCTION IS BEING DEPRECATED - ALWAYS USE get_adm0_codes() INSTEAD.
#'
get_gaul_codes <- function(
                           adm0s,
                           strict = FALSE,
                           lookup_table = NULL,
                           core_repo = NULL,
                           adm0_type = 'gaul',
                           shapefile_version = NULL
                           ){
  get_adm0_codes(
    adm0s        = adm0s,
    strict       = strict,
    lookup_table = lookup_table,
    core_repo    = core_repo,
    adm0_type    = adm0_type,
    shapefile_version = shapefile_version
  )
}


## get_adm_codes_subnat() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Get subnational administrative codes (adm1 and adm2)
#'
#' @description Pulls subnational administrative codes worldwide for a given
#'   administrative level (1 or 2) and a given set of country-level geographies,
#'   represented by a vector of standard Admin0 codes pulled using
#'   get_adm0_codes(). This function generalizes get_gaul_codes_subnat(), which
#'   will be deprecated.
#'
#' @param adm0_list A numeric vector of Admin0 codes pulled using
#'   get_adm0_codes().
#' @param admin_level A subnational administrative level: 1 and 2 are currently
#'   supported.
#' @param shapefile_version A character string indicating which shapefile version to pull
#'
get_adm_codes_subnat <- function(adm0_list, admin_level, shapefile_version = 'current'){
  # Read the table data for the current administrative shapefile
  hierarchy <- read.dbf(get_admin_shapefile(admin_level, suffix=".dbf",
                                            version = shapefile_version))
  hierarchy <- as.data.table(hierarchy)
  # Rename the admin0 field to ADM0_CODE
  adm0_name <- grep("ADM0_CODE", names(hierarchy), value = T)
  names(hierarchy)[names(hierarchy)==adm0_name] <- "ADM0_CODE"
  # Keep only rows within the given Admin0s, and keep only the column for the
  #  relevant administrative level codes
  hierarchy <- hierarchy[ADM0_CODE %in% adm0_list,]
  adm_specific_name <- grep(
    paste0('ADM', admin_level, '_CODE'), names(hierarchy),
    value = T
  )
  adm_list <- unique(hierarchy[[adm_specific_name]])
  return(adm_list)
}


## get_gaul_codes_subnat() ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#'
#' @title Get subnational GAUL codes
#'
#' @description Pulls subnational administrative codes worldwide for a given
#'   administrative level (1 or 2) and a given set of country-level geographies,
#'   represented by a vector of standard Admin0 codes pulled using
#'   get_adm0_codes(). This function is a simple wrapper for
#'   get_adm_codes_subnat() and has been kept for backwards-compatibility.
#'
#' NOTE: THIS FUNCTION IS BEING DEPRECATED. ALWAYS USE get_adm_codes_subnat()
#'   INSTEAD.
#'
get_gaul_codes_subnat <- function(gaul_list, admin_level, shapefile_version) {
  subnat_list <- get_adm_codes_subnat(
    adm0_list   = gaul_list,
    admin_level = admin_level,
    shapefile_version = shapefile_version
  )
  return(subnat_list)
}



add_gauls_regions <- function(df, simple_raster, shapefile_version = 'current') {

  # Extract GAUL_CODE from simple_raster using lat/longs
  df$ADM0_CODE <- raster::extract(simple_raster, df[ , c('longitude', 'latitude'), with=F])

  # Add names of regions by GAUL_CODE
  for(r in Regions) {
    df <- df[ADM0_CODE %in% get_adm0_codes(r, shapefile_version = shapefile_version), region := r]
  }

  ##Check that merge of GAUL_CODE worked properly
  df <- df[!is.na(region), good_records := 1]
  message(paste0(length(df$good_records) - length(df[good_records==1, N]), ' out of ', length(df$good_records), ' could not have GAUL/region extracted properly. Probably coastal points, need to fix.'))
  df <- df[good_records==1, ]

}


gaul_convert <- function(countries, from = "iso3", verbose = F, shapefile_version = 'current') {

  # Purpose: Convert a vector of countries (ihme_loc_id format) to vector of GAUL codes
  # Inputs:
  #         countries: vector of countries in ihme_loc_id format
  #         from: format of input
  #               options: "iso3" = "ihme_loc_id" = "ihme_lc_id"
  #                        "name" (the loc_name or loc_nm_short)
  #               ("iso3" is treated as "ihme_loc_id" for backwards compatability)
  #
  # Outputs: a vector of gaul codes

  # load reference table
  str_match <- stringr::str_match

  # Catch if already passed gaul codes
  if(class(countries) =="numeric") return (countries)
  if(all(grepl("^[[:digit:]]+$", countries))) return(countries)

  gaul_table <- get_location_code_mapping(shapefile_version = shapefile_version)

  # convert input & output to lower case for easier matching
  #lowercase_cols <- c("short_name", "official_name", "iso3", "iso2", "uni", "undp")
  #gaul_table[, (lowercase_cols) := lapply(.SD, tolower), .SDcols = lowercase_cols,]

  ## ## convert input & output to lower case for easier matching
  if (verbose == T) message("\nlowercasing columns. the columns that get lowered are:")
  for(i in 1:ncol(gaul_table)){
    if(any(class(gaul_table[[i]]) %in% c("factor", "character"))) {
      if (verbose == T) message(sprintf("On column: %s", colnames(gaul_table)[i]))
      gaul_table[[i]] <- tolower(gaul_table[[i]])
    }
  }

  # Lowercase & ensure character for input
  countries <- tolower(countries)
  countries <- as.character(countries)

  ## Catch if a subnational in XXX_##### IHME syntax
  ## This returns national-level gaul codes only
  if(any(grepl("_", countries))) {
    countries[grepl("_", countries)] <- str_match(countries[grepl("_", countries)], "(.*)_.*")[,2]
  }

  if (from == "iso3" | from == "ihme_loc_id" | from == "ihme_lc_id") {

    gaul_code <- sapply(countries, function(x) gaul_table[ihme_lc_id == x, GAUL_CODE]) %>% as.numeric

  } else if(from == "name") {

    # Matching only national-level for now
    # Drop undefined & subnational rows
    gaul_table_nat <- subset(gaul_table, GAUL_CODE != -1)
    gaul_table_nat <- subset(gaul_table_nat, level == 3)

    gaul_code <- sapply(countries, function(x) gaul_table_nat[loc_nm_sh == x, GAUL_CODE]) %>% as.numeric

    #check to see if this matched all of the provided items in the vector; use partial / fuzzy matching if not
    if(length(gaul_code[is.na(gaul_code)]) > 0) {

      # Create a table to fill in
      table_matching <- cbind(countries, gaul_code) %>% as.data.table
      names(table_matching) <- c("country", "gaul_code")
      table_matching$gaul_code <- as.numeric(table_matching$gaul_code)

      approx_matched <- table_matching[is.na(gaul_code), country]

      # Indicate that approximate matching took place

      message("\nNot all country names provided were found in the lookup table.")
      message("Attempting to match names provided with those in lookup table.")
      message(paste0("Approximate matching attempted for: ", paste(approx_matched, collapse = ', '), "\n"))

      approx_match <- function(country) {
        # First, try matching to long form of name
        gaul_code <- gaul_table_nat[grep(country, gaul_table_nat$loc_name),]$GAUL_CODE

        # If that doesn't work, grep within the short name
        if (length(gaul_code) == 0) gaul_code <- gaul_table_nat[grep(country, gaul_table_nat$loc_nm_sh),]$GAUL_CODE

        # If that doesn't work, grep within the long name
        if (length(gaul_code) == 0) gaul_code <- gaul_table_nat[grep(country, gaul_table_nat$loc_name),]$GAUL_CODE

        # Could fill in other matching here if desired

        # Warn if nonspecific
        if (length(gaul_code) > 1) warning(paste0("\"", country, "\" matches multiple country names in the lookup table. Please be more specific."))

        # Finally, if no matches, return NA
        if (length(gaul_code) != 1) gaul_code <- NA

        return(as.numeric(gaul_code))
      }

      # Try approximate matching
      table_matching[is.na(gaul_code)]$gaul_code <- sapply(table_matching[is.na(gaul_code)]$country,approx_match)

      not_matched <- table_matching[is.na(gaul_code)]$country

      # Error checking
      if(length(not_matched) > 0) {
        warning(paste0("Some countries could not be matched:\n", paste(not_matched, collapse=', ')))
      }
      gaul_code <- table_matching$gaul_code %>% as.numeric
    }

  } else {
    # Error catching for non-supported country type
    stop("\nPlease enter a valid country code type")
  }

  if(length(gaul_code[is.na(gaul_code)]) > 0){
    # Error catching for failure to match all country codes
    message(paste0("CAUTION! Returning NA values.\nMatches not found for all country codes in input list.\n",
                "Please check your input values"))
  }
  return(gaul_code)
}

## Make map of period indices to run any time periods in your data (like annual instead of 5-year)
  make_period_map <- function(modeling_periods) {
    data_period <- sort(modeling_periods)
    period_ids <- seq(data_period)
    period_map <- as.data.table(data_period)
    period_map <- period_map[, period_id := period_ids]
    return(period_map)
  }

  interpolate_gbd <- function(gbd) {
    new_gbd <- list()
    for(this_year in c(2000:2015)) {
      if(this_year %in% 2000:2004) copied_data <- gbd[year == 2000, ]
      if(this_year %in% 2005:2009) copied_data <- gbd[year == 2005, ]
      if(this_year %in% 2010:2014) copied_data <- gbd[year == 2010, ]
      if(this_year %in% 2015:2015) copied_data <- gbd[year == 2015, ]
      copied_data <- copied_data[, year := this_year]
      new_gbd[[as.character(this_year)]] <- copied_data
    }
    new_gbd <- rbindlist(new_gbd)
    new_gbd <- new_gbd[order(name, year)]
    new_gbd <- new_gbd[!(year %in% c(2000,2005,2010,2015)), mean := NA]
    library(zoo)
    for(country in unique(new_gbd[, name])) {
      new_gbd <- new_gbd[name == country, mean := na.approx(new_gbd[name == country, mean])]
    }
    return(new_gbd)
  }

make_regions_map <- function(regions, subset_shape, shapefile_version = 'current') {
  col_list <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
  i <- 1
  plot(subset_shape)
  for(reg in regions) {
    shapes <- subset_shape[subset_shape$GAUL_CODE %in% get_adm0_codes(reg, shapefile_version = shapefile_version), ]
    plot(shapes, add=TRUE, col=col_list[i])
    i <- i + 1
  }
}


merge_with_ihme_loc <-function(d, re=Regions, shapefile_version = 'current'){
  message(nrow(d))
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  d <- d[, ihme_lc_id := as.character(country)]

  d <- merge(d, gaul_to_loc_id, by='ihme_lc_id',all.x=T)

  message(nrow(d))
  for(r in re) {
    d <- d[GAUL_CODE %in% get_adm0_codes(r, shapefile_version = shapefile_version), region := r]
  }
  return(d)
}


## check_custom_regions ################################################

#' Check custom region list against base MBG list and see what's added
#' or missing in comparison
#'
#' @param region_list vector of region names (defined in `pull_custom_modeling_regions()`)
#' @return list of two tables
#'         [[1]] gaul_codes in custom, but not default regions,
#'           [[2]] gaul_codes in default, but not custom regions, and
#'       [[3]] duplicated gaul_codes
#' @examples
#'
#' region_list <- c('vax_soas','vax_seas','vax_eaas','vax_caeu',
#'                 'vax_crbn','vax_ctam','vax_ansa','vax_trsa',
#'                 'vax_name','vax_cssa','vax_essa','vax_sssa','vax_wssa')
#' check_list <- check_custom_regions(region_list)

check_custom_regions <- function(region_list, shapefile_version = 'current') {

  # Function to pull a table of gauls by region
  make_gaul_table <- function(rr) {
    data.table(region = rr,
             gaul_code = get_adm0_codes(rr, shapefile_version = shapefile_version))
  }

  # Pull gauls for custom regions
  custom_reg_table <- lapply(region_list, make_gaul_table) %>% rbindlist
  setnames(custom_reg_table, "region", "custom_region")

  # Pull gauls for default regions
  default_reg_table <- lapply(c('stage1', 'stage2'), make_gaul_table) %>% rbindlist
  setnames(default_reg_table, "region", "default_region")

  # Merge tables together
  reg_table <- merge(custom_reg_table, default_reg_table, all.x = T, all.y = T)

  # Grab the rest of the information about these gaul codes
  info_table <- fread("/<<<< FILEPATH REDACTED >>>>/stage_master_list.csv")
  reg_table <- merge(reg_table, info_table, all.x = T, all.y = F, by.x = "gaul_code", by.y = "GAUL_CODE")

  custom_not_default <- subset(reg_table, !is.na(custom_region) & is.na(default_region))
  default_not_custom <- subset(reg_table, is.na(custom_region) & !is.na(default_region))

  # Check for duplicated gauls
  dups <- reg_table$gaul_code[duplicated(reg_table$gaul_code)] %>% unique
  duplicated_table <- subset(reg_table, gaul_code %in% dups)

  message(paste0("There are ", nrow(custom_not_default), " gaul code(s) in the CUSTOM list not included in the default MBG regions."))
  message(paste0("There are ", nrow(default_not_custom), " gaul code(s) in the DEFAULT mbg regions not included in your custom list."))
  message(paste0("There are ", length(dups), " duplicated gaul code(s) in the CUSTOM list"))
  message(paste0("\nReturning a list of three tables:\n",
           "  1) gaul codes in custom, but not default regions\n",
           "  2) gaul codes in default, but not custom regions\n",
           "  3) duplicated gaul codes"))

    return(list(custom_not_default = custom_not_default,
                default_not_custom = default_not_custom,
                duplicated = duplicated_table))
}

## plot_custom_regions() ################################################

#' Make a world plot of custom region lists
#'
#' @param region_list vector of region names (defined in `pull_custom_modeling_regions()`)
#' @return ggplot object of the map of the custom regions
#' @examples
#' png(file = "/path/to/file.png",
#'     width = 14,
#'     height = 6,
#'     units = "in",
#'     res = 300)
#'
#' plot_custom_regions(region_list, plot_title = "Vaccine Regions")
#'
#' dev.off()

plot_custom_regions <- function(region_list,
                                plot_title = NULL,
                                verbose = T,
                                shapefile_version = 'current') {

  # Replace some functions to be quieter
  geom_polygon_quiet <- function(...) {suppressMessages(ggplot2::geom_polygon(...))}
  geom_path_quiet     <- function(...) {suppressMessages(ggplot2::geom_path(...))}

  # Function to pull a table of gauls by region
  make_gaul_table <- function(rr) {
    data.table(region = rr,
             gaul_code = get_adm0_codes(rr, shapefile_version = shapefile_version))
  }

  reg_table <- lapply(region_list, make_gaul_table) %>% rbindlist

  # Add in stage 1/2 countries that aren't in custom list
  default_reg_table <- data.table(in_s1s2 = T,
                                  gaul_code = get_adm0_codes(c('stage1', 'stage2'), shapfile_version = shapefile_version))

  reg_table <- merge(reg_table, default_reg_table, all.x = T, all.y = T)
  reg_table[is.na(region) & in_s1s2 == T, region := "missing"]

  # Load gaul shapefile
  if (verbose == T) message("Loading shapefiles...")
  world_shape <- readRDS('/<<<< FILEPATH REDACTED >>>>/g2015_2014_0_simp_tol_0.1.rds')
  world_shape <- subset(world_shape, ADM0_NAME != "Antarctica")
  world_shape <- merge(world_shape, unique(reg_table), all.x = T, all.y = F, by.x = "ADM0_CODE", by.y = "gaul_code")

  # Make the map

  # A custom, mostly blank theme to use for mapping
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank())

  # Discrete color palette from Carto colors
  carto_discrete <- c("#7F3C8D","#11A579","#F2B701","#E73F74",
                      "#3969AC","#80BA5A","#E68310","#008695",
                      "#CF1C90","#f97b72","#4b4b8f","#A5AA99",
                      "#66C5CC","#F6CF71","#F89C74","#DCB0F2",
                      "#87C55F","#9EB9F3","#FE88B1","#C9DB74",
                      "#8BE0A4","#B497E7","#D3B484","#B3B3B3")

  color_scheme <- carto_discrete[1:length(region_list)]
  names(color_scheme) <- region_list

  # Fortify once
  if (verbose == T) message("Fortifying shapefile...")
  world_shape@data$id <- rownames(world_shape@data)
  world_shape_df <- suppressMessages(as.data.table(fortify(world_shape, region="id")))
  world_shape_df <- merge(world_shape_df, as.data.table(world_shape@data))

  # Get centroids for small areas in order to overlay circles
  if (verbose == T) message("Calculating centroids for small countries...")
  small_shapes <- subset(world_shape, area(world_shape) < 1e11)

  cents <- gCentroid(small_shapes, byid = T, id = small_shapes$id) %>%
    as.data.frame() %>%
    cbind(rownames(.), .) %>%
    as.data.table %>%
    setnames(., names(.), c("id", "x", "y"))

  cents <- merge(cents, small_shapes@data, by="id") %>% as.data.table

  cents <- subset(cents, !is.na(region) | (in_s1s2 == TRUE))

  # Fix Kiribati if present -- spans 0 long so centroid wraps to close to 0,0
  cents[ADM0_CODE == 135, x := -154.75583729] # http://www.marineregions.org/gazetteer.php?p=details&id=8441
  cents[ADM0_CODE == 135, y := -3.81449579] # http://www.marineregions.org/gazetteer.php?p=details&id=8441

  # Make map and add on options
  if (verbose == T) message("Drawing map...")
  gg_map <- ggplot() +
    geom_polygon_quiet(data = world_shape_df,
                       aes(x = long,
                           y = lat,
                           group = group),
                           fill = "gray")

  if (length(subset(world_shape_df, region == "missing")) > 0) {
    gg_map <- gg_map + geom_polygon_quiet(data = subset(world_shape_df, region == "missing"),
                       aes(x = long,
                           y = lat,
                           group = group),
                           fill = "white")
  }

  # Plot regions by color and overlay borders
  gg_map <- gg_map +
    geom_polygon_quiet(data = subset(world_shape_df, region != "missing" & !is.na(region)),
                       aes(x = long,
                           y = lat,
                           group = group,
                           fill = region)) +
    geom_path_quiet(data = world_shape_df,
                    aes(x = long,
                        y = lat,
                        group = group),
                    size = 0.2,
                    color = "black")

  # If small geographies present, overlay a circle for easier visualization
  if (nrow(cents) > 0) {
    gg_map <- gg_map +
      geom_point(data = cents,
                 aes(x = x, y = y, color = region),
                 size = 4,
                 alpha = 0.5)
  }

  # Add on plot options including title if needed
  gg_map <- gg_map +
    theme_empty +
    coord_equal(ratio = 1) +
    scale_fill_manual(name = "Region",values = color_scheme) +
    scale_color_manual(name = "Region",values = color_scheme) +
    guides(color = FALSE)

  if (!is.null(plot_title)) gg_map <- gg_map + labs(title = plot_title)

  # Return the final map
  return(gg_map)
}


#' load_libs
#'
#' @description A function to load libraries based on whether the
#' code is being run on cluster prod or the geos nodes. Tries
#' to load each package, and returns errors or warnings if the package fails to
#' load rather than breaking the loop.
#'
#' @param packages A vector or list of packages
#' @param stop Boolean; if T will stop code if not all packages load. if F (default),
#' the code will contiue and will throw a warning.
load_libs<-function(packages,stop=F){
  package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                        paste0('<<<< FILEPATH REDACTED >>>>/geos_packages'),
                        paste0('<<<< FILEPATH REDACTED >>>>/packages'))
  .libPaths(package_lib)
  message(paste0('Loading packages from ',package_lib))

  if(Sys.info()[1] == "Windows"){stop("STOP! you will overwrite these packages if you run from windows")}

  # Try loading each of the packages in the list provided
  for(package in c(packages)) {
    try({library(package, lib.loc = package_lib, character.only=TRUE)},silent=T)
  }

  failed_packages<-packages[!packages %in% (.packages())]

  if(length(failed_packages>0)){
    message(paste0("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n",
                   "The following packages failed to load: \n",
                   "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n"))
    message(paste0(" ",failed_packages))
    message("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n")
    if(stop){
      stop("Code halted; important packages failed to load!")
    }else{
      warning("Some packages did not load properly (see summary above), but we are going to proceed anyway.")
    }
  }else{
    message(paste0("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n",
                   "Packages have been loaded! \n",
                   "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n"))
  }
}


#' crop_set_mask
#'
#' @description  Crop the raster to the extent of a simple_raster object,
#' set the extent, and then mask it based on the simple_raster object.
#'
#'
#' @param raster_object A raster object (raster, Brick, etc) of a specific region
#' @param simple_raster A simple raster object that serves as the template for that region
#' @return A raster/Brick/Stack, with the extents and masks of the simple raster.
crop_set_mask<-function(raster_object,template_raster){
  raster_object<-raster::crop(raster_object, extent(template_raster))
  raster_object<-raster::setExtent(raster_object, template_raster)
  raster_object<-raster::mask(raster_object,template_raster)
  return(raster_object)
}


#' @title Uploading config info to SQLite database
#' @description This function will upload the MBG config info into a SQLite database. There are no primary keys set, so any number of duplicate entries can be recorded.
#'
#' @param config The config \code{data.table} prepped apriori.
#' @param user User name
#' @param core_repo Path to core repository
#' @param indicator_group Indicator group
#' @param indicator Indicator
#' @param run_date Run date datestamp
#' @param verbose Print config before exiting?
#'
#' @importFrom DBI dbConnect dbWriteTable dbDisconnect
#' @importFrom RSQLite SQLite
#'
#' @export
#'
upload_config_to_db <- function(config, user, core_repo, indicator_group, indicator, run_date, verbose = FALSE) {

  ## Prep config file (transpose and add fields)
  config_t <- transpose(config)[2,]
  colnames(config_t) <- config$V1

  ## Add missing columns
  missing_metadata <- data.table(user = user,
                                 core_repo = core_repo,
                                 indicator_group = indicator_group,
                                 indicator = indicator,
                                 run_date = run_date)
  config_binded <- cbind(missing_metadata, config_t)

  ## Add datestamp of adding the config
  config_binded[, datestamp:= paste0(Sys.time())]

  ## Open connection to database
  dbpath <- "/<<<< FILEPATH REDACTED >>>>/run_db_20190216.sqlite"
  configdb <- DBI::dbConnect(RSQLite::SQLite(), dbpath)

  ## Upload config
  DBI::dbWriteTable(configdb, "config_upload", config_binded, append = TRUE)

  ## Disconnect db
  DBI::dbDisconnect(configdb)

  ## Print config?
  if(verbose) print(config_binded)


  print("Config data uploaded")

  return(NULL)

}
