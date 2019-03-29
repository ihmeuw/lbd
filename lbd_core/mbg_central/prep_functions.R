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

## Load parameters from config file into memory
#   Arguments:
#     repo            = Location where you've cloned "mbg" repository.
#     indicator_group = Category of indicator
#     indicator       = Specific outcome to be modeled within indicator categor
  load_config <- function(repo, indicator_group, indicator, config_name=NULL, covs_name = NULL, post_est_only=FALSE, run_date = '') {

    # Pull from config .csv file
    if (is.null(config_name)) {
      ## If new model run, pull config from repo
      if(post_est_only==FALSE) config <- fread(paste0(repo, '/', indicator_group, '/config_', indicator, '.csv'), header=FALSE)
      ## If running analysis on existing model, use config from that model's outputs folder
      if(post_est_only==TRUE) config <- fread('<<<< FILEPATH REDACTED >>>>/config.csv')
    } else {
      config <- fread(paste0(repo, '/', indicator_group, '/', config_name, '.csv'), header=FALSE)
    }

    # If a covariate .csv file exists, use that instead
    if (!is.null(covs_name)) {

      # Grab fixed effects & measures (and gbd fixed effects & measures) from CSV if present
      covs <- fread(paste0(repo, '/', indicator_group, '/', covs_name, '.csv'), header = TRUE)
      # After update to data.table 1.11.4, "T" and "F" are not read in as
      # logical, but as characters, which we need to remedy here. We are
      # assuming that the "covs.csv" has an "include" and "gbd" column here
      covs[, gbd := as.logical(gbd)]
      covs[, include := as.logical(include)]
      covs <- subset(covs, include == T) # Use only those where include flag set to true
      fe <- subset(covs, gbd == F)
      gbd <- subset(covs, gbd == T)
      gbd[measure != "output", measure := "covariate"] # For backwards compatability -- basically it assumes you meant 'covariate' if you specified anything other than 'outcome' (eg, mean or NA)
      fixed_effects <- paste0(fe$covariate, collapse = " + ")
      fixed_effects_measures <- paste0(fe$measure, collapse = " + ")
      gbd_fixed_effects <- paste0(gbd$covariate, collapse = " + ")
      gbd_fixed_effects_measures <- paste0(gbd$measure, collapse = " + ")

      # Remove any other versions from original config
      config <- subset(config, !(V1 %in% c("fixed_effects", "fixed_effects_measures", "gbd_fixed_effects", "gbd_fixed_effects_measures")))
      config <- config %>%
        rbind(., list("fixed_effects", fixed_effects)) %>%
        rbind(., list("fixed_effects_measures", fixed_effects_measures)) %>%
        rbind(., list("gbd_fixed_effects", gbd_fixed_effects)) %>%
        rbind(., list("gbd_fixed_effects_measures", gbd_fixed_effects_measures))
    }

    # Assign all the covariates to the environment
    for (param in config[, V1]) {
      assign(param, config[V1==param, V2], envir=globalenv())
    }

    return(config)
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
    if(use_share==FALSE) load_dir <- "<<<< FILEPATH REDACTED >>>>"
    if(use_share==TRUE) load_dir  <- "<<<< FILEPATH REDACTED >>>>"

    if(!withdate & !withtag) filename <- paste0(load_dir, indicator)
    if(withtag)              filename <- paste0(load_dir, indicator, datatag)
    if(withdate)             filename <- "<<<< FILEPATH REDACTED >>>>"

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
      if(dir.exists("<<<< FILEPATH REDACTED >>>>") == TRUE) {
        existing_dir <- "<<<< FILEPATH REDACTED >>>>"
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
      if(dir.exists("<<<< FILEPATH REDACTED >>>>" == FALSE)) {
        run_date_dir <- "<<<< FILEPATH REDACTED >>>>"
        dir.create(run_date_dir, showWarnings = FALSE)
      }
      write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
      return(list(d, run_date))
    }

    if(update_run_date == FALSE) {
      if(agebin==0){
        run_date_dir <- "<<<< FILEPATH REDACTED >>>>"
        dir.create(run_date_dir, showWarnings = FALSE)
        write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
      } else {
        run_date_dir <- "<<<< FILEPATH REDACTED >>>>"
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

load_simple_polygon <- function(gaul_list, buffer, tolerance = 0.2,
                                subset_only = F, makeplots = F, use_premade = F,
                                custom_shapefile_path = NULL, custom_shapefile = NULL, raking = F) {

  # Logic check
  if (!is.null(custom_shapefile_path) & !is.null(custom_shapefile)) stop("You cannot specify both a custom shapefile and a custom shapefile path")

  # If using custom shapefiles, don't use premade polys
  if (!is.null(custom_shapefile_path) | !is.null(custom_shapefile)) use_premade <- F

  if (use_premade == T) {
    # Check for common gaul lists so we can query premade simple polygons
    if(identical(gaul_list, c(29, 42, 45, 47, 50, 66, 90, 94, 106, 105, 144, 155, 159, 181, 182, 214, 217, 221, 243))) {
      message("Your GAUL list matches WSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load("<<<< FILEPATH REDACTED >>>>")
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(8,49,59,68,76,89))) {
      message("Your GAUL list matches CSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load("<<<< FILEPATH REDACTED >>>>")
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(43,58,70,77,79,133,150,152,170,205,226,74,257,253,270))) {
      message("Your GAUL list matches ESSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load("<<<< FILEPATH REDACTED >>>>")
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(4,40762,40765,145,169,6,248))) {
      message("Your GAUL list matches NAME, loading premade simple_polygon and assigning subset_shape to global environment...")
      load("<<<< FILEPATH REDACTED >>>>")
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(35,142,172,227,235,271))) {
      message("Your GAUL list matches SSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load("<<<< FILEPATH REDACTED >>>>")
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
      load("<<<< FILEPATH REDACTED >>>>")
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
    master_shape <- readOGR(get_admin_shapefile(admin_level = 0, raking = raking))
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

## Create spatial mesh
#   d = Data table with columns "latitude" and "longitude".
#   simple = Single polygon to create mesh over.
#   max_edge = Largest allowed triangle edge length (see INLA's inla.mesh.2d documentation).
#   mesh_offset = One or two values, for an inner and an optional outer extension distance (see INLA's inla.mesh.2d documentation).
build_space_mesh <- function(d, simple, max_edge, mesh_offset, plot_mesh = F) {
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

    return(mesh_s)
}

## Create temporal mesh
build_time_mesh <- function(periods=1:4) {
  mesh_t <- inla.mesh.1d(
    loc = c(periods),
    degree = 1,
    boundary = rep("free", 2)
  )
  return(mesh_t)
}

## Load list of raster inputs (pop and simple)
build_simple_raster_pop <- function(subset_shape, u5m=FALSE, field=NULL, raking=F) {

  if (is.null(field)) {
    if ('GAUL_CODE' %in% names(subset_shape@data)) field <- 'GAUL_CODE'
    if ('ADM0_CODE' %in% names(subset_shape@data)) field <- 'ADM0_CODE'
  }

  if(raking){
    field <- 'loc_id'
  }

  if(u5m==FALSE){
    master_pop <- brick("<<<< FILEPATH REDACTED >>>>")
  } else {
    master_pop <- brick(raster("<<<< FILEPATH REDACTED >>>>"),
    raster("<<<< FILEPATH REDACTED >>>>"),
    raster("<<<< FILEPATH REDACTED >>>>"),
    raster("<<<< FILEPATH REDACTED >>>>"))
  }

  cropped_pop <- crop(master_pop, extent(subset_shape), snap="out")
  ## Fix rasterize
  initial_raster <- rasterize(subset_shape, cropped_pop, field = field)
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))!=0) {
    rasterized_shape <- merge(rasterize(subset(subset_shape, !(get(field) %in% unique(initial_raster))), cropped_pop, field = field), initial_raster)
  }
  if(length(subset(subset_shape, !(get(field) %in% unique(initial_raster))))==0) {
    rasterized_shape <- initial_raster
  }
  #rasterized_shape <- rasterize(subset_shape, cropped_pop, field='GAUL_CODE')
  masked_pop <- raster::mask(x=cropped_pop, mask=rasterized_shape)

  raster_list <- list()
  raster_list[['simple_raster']] <- rasterized_shape
  raster_list[['pop_raster']] <- masked_pop

  return(raster_list)

}

## #############################################################################
## GET GAUL CODES AND RELATED FUNCTIONS
## #############################################################################


## load_gaul_lookup_table() ----------------------------------------------------
#'
#' @title Load the GAUL lookup table
#'
#' @description Loads the most recent version of the lookup table that links
#'   GAUL codes with other identifiers such as GBD location IDs, MBG modeling
#'   regions, and ISO codes, among others.
#'
#' @return Returns a data.frame of the GAUL lookup table. If package data.table
#'   is loaded, also converts the lookup table to a data.table
#'
load_gaul_lookup_table <- function(){
  # Define consistent path to the GAUL lookup table
  lookup_table_filepath <- "<<<< FILEPATH REDACTED >>>>"
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
#'    by get_gaul_codes().
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
    'middle_east_dcp' = 'mide+stan-yem-pak+tur+isr+leb',
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

    # For modelers who are still using the defunct NAME region
    'name_historic' = 'noaf-esh'
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



## get_gaul_codes() ------------------------------------------------------------
#'
#' @title Get GAUL Codes
#'
#' @description Pull GAUL Codes for the specified countries or regions,
#'   optionally excluding certain GAUL codes.
#'
#' @param gaul Vector of ISO-3 codes, four-letter modeling region names, region
#'   group names as defined in get_region_groups, or 'all' (returns all GAUL
#'   codes)
#' @param strict (default `FALSE`) Causes the function to fail if an ISO code or
#'   region name is not included in the full list.
#' @param lookup_table (default `NULL`) Sets a data.frame or data.table to be
#'   used as the lookup_table. If `NULL`, the lookup table will be loaded using
#'   the `load_gaul_lookup_table()` function.
#' @param core_repo (default `FALSE`) Causes the function to fail if an ISO code or
#'   region name is not included in the full list.
#'
#' @return Vector of numeric GAUL codes. Warns if not all regions align with a
#'   GAUL code.
#'
get_gaul_codes <- function(gaul,
                           strict = FALSE,
                           lookup_table = NULL) {
  # Load data.table
  suppressMessages(load_R_packages(c('data.table','assertthat')))

  # Read in the GAUL code data.table once at the beginning of the function
  if(is.null(lookup_table)) lookup_table <- load_gaul_lookup_table()

  ## PROCESS INCOMING GAUL TEXT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # In case multiple character strings were passed in, concatenate with '+'
  gaul <- paste(gaul, collapse='+')
  # Set all as lowercase
  gaul <- tolower(gaul)

  # Split on - and + characters
  gaul <- gsub(';','',gaul)
  gaul <- gsub('\\+',';;\\+',gaul)
  gaul <- gsub('-',';;-',gaul)
  gaul_vec <- unlist(base::strsplit(gaul, ';;'))

  # Making the function past-compatible, for modelers who are still using the
  #  original 'name' modeling region for North Africa (now 'noaf').
  #  'name_historic' is now a custom region, equivalent to 'noaf-esh', in the
  #  pull_custom_modeling_regions().
  gsub('^([+-])?name$', '\\1name_historic', gaul_vec)

  # If the string was empty, end now
  if (length(gaul_vec)==0){
    message("No GAUL codes were returned; exiting")
    return(numeric(0))
  }

  # The first region and all regions beginning in '+' are included
  # Strip the beginning + signs
  include <- gsub('\\+','', c(gaul_vec[1], gaul_vec[ grepl('^\\+',gaul_vec) ]) )
  # All regions beginning with '-' are treated as exclusions
  # Strip any minus signs from them
  exclude <- gsub('-','', gaul_vec[ grepl('^-',gaul_vec) ] )

  ## CHECK FOR CUSTOM REGIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define regex for all non-custom regions
  standard_code_regex <- "^([a-z]{3,4}|[a-z]{3}_d[0-9]{1,5})$"
  # Split both included and excluded data into custom and non-custom regions
  standard_inc <- include[ grepl(standard_code_regex, include) ]
  custom_inc <- include[ !grepl(standard_code_regex, include) ]
  # Check in excluded data
  standard_exc <- exclude[ grepl(standard_code_regex, exclude) ]
  custom_exc <- exclude[ !grepl(standard_code_regex, exclude) ]

  # Get definitions and GAUL codes for included and excluded GAULs, if needed
  # `custom_inc_codes` and `custom_exc_codes` will be added/subtracted at the end
  # Check out that functional recursion
  if ( length(custom_inc) > 0 ){
    custom_inc_defs  <- pull_custom_modeling_regions(custom_inc)
    custom_inc_codes <- unlist(lapply(custom_inc_defs,function(x) get_gaul_codes(x,
                                                         lookup_table=lookup_table)
                               ))
  } else {
    custom_inc_codes <- integer(0)
  }
  if ( length(custom_exc) > 0 ){
    custom_exc_defs  <- pull_custom_modeling_regions(custom_exc)
    custom_exc_codes <- unlist(lapply(custom_exc_defs,function(x) get_gaul_codes(x,
                                                         lookup_table=lookup_table)
                               ))
  } else {
    custom_exc_codes <- integer(0)
  }

  ## PROCESS NON-CUSTOM REGIONS (MBG MODELING REGIONS + ISOS) ~~~~~~~~~~~~~~~~~~
  # Helper function to pull all GAUL codes for a set of modeling regions or
  #  ISO codes
  utility_pull_gauls <- function(region_codes, standard_regex){
    # Return none if the length of the inputs is 0
    if (length(region_codes)==0) return(integer(0))
    # Input data assertions
    assert_that(all(grepl(standard_regex, region_codes)),
                msg=paste0('All region codes must match the MBG or ISO code formats.'))
    isos     <- region_codes[nchar(region_codes)!=4]
    mbg_regs <- region_codes[nchar(region_codes)==4]
    if ('all' %in% isos){
      # 'all' is a special case where all GAUL codes are pulled
      pulled_gauls <- unique(lookup_table[,GAUL_CODE])
    } else (
      # Pull all GAUL codes matching the ISOs or regions
      pulled_gauls <- unique(lookup_table[(iso3 %in% isos) | (mbg_reg %in% mbg_regs),
                                          GAUL_CODE])
    )
    # Warn if any iso codes or MBG regions aren't in the lookup table
    missing_isos <- isos[ !(isos %in% lookup_table[,iso3]) ]
    if (length(missing_isos) > 0) message(paste0("WARNING: Missing these ISOs: ",
                                                 paste(missing_isos,collapse=',')))
    missing_regs <- mbg_regs[ !(mbg_regs %in% lookup_table[,mbg_reg]) ]
    if (length(missing_regs) > 0) message(paste0("WARNING: Missing these MBG regions: ",
                                                 paste(missing_regs,collapse=',')))
    return(as.integer(pulled_gauls))
  }

  # Pull GAUL codes for standard include and exclude regions
  standard_inc_codes <- utility_pull_gauls(standard_inc, standard_code_regex)
  standard_exc_codes <- utility_pull_gauls(standard_exc, standard_code_regex)

  ## COMBINE TO CREATE FINAL GAUL CODE SET ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get all (standard + custom) codes to include
  all_include_codes <- unique( c(standard_inc_codes, custom_inc_codes) )
  # Get all (standard + custom) codes to exclude
  all_exclude_codes <- unique( c(standard_exc_codes, custom_exc_codes) )
  # Final codes = codes to include - codes to exclude
  return_codes <- all_include_codes[ !(all_include_codes %in% all_exclude_codes) ]

  return(return_codes)
}



gaul_convert <- function(countries, from = "iso3", verbose = F) {

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

  gaul_table <- get_location_code_mapping()

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

make_regions_map <- function(regions, subset_shape) {
  col_list <- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999')
  i <- 1
  plot(subset_shape)
  for(reg in regions) {
    shapes <- subset_shape[subset_shape$GAUL_CODE %in% get_gaul_codes(reg), ]
    plot(shapes, add=TRUE, col=col_list[i])
    i <- i + 1
  }
}


merge_with_ihme_loc <-function(d,re=Regions){
  message(nrow(d))
  gaul_to_loc_id <- get_location_code_mapping()
  d <- d[, ihme_lc_id := as.character(country)]

  d <- merge(d, gaul_to_loc_id, by='ihme_lc_id',all.x=T)

  message(nrow(d))
  for(r in re) {
    d <- d[GAUL_CODE %in% get_gaul_codes(r), region := r]
  }
  return(d)
}

