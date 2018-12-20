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
#     indicator_group = Category of indicator, i.e. "education"
#     indicator       = Specific outcome to be modeled within indicator category, i.e. "edu_0"
  load_config <- function(repo, indicator_group, indicator, config_name=NULL, covs_name = NULL, post_est_only=FALSE, run_date = '') {

    # Pull from config .csv file
    if(is.null(config_name)) {
      ## If new model run, pull config from repo
      if(post_est_only==FALSE) config <- fread(paste0(repo, indicator_group, '/config_', indicator, '.csv'), header=FALSE)
      ## If running analysis on existing model, use config from that model's outputs folder
      if(post_est_only==TRUE) config <- fread(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/config.csv'))
    }
    if(!is.null(config_name)) config <- fread(paste0(repo, indicator_group, '/', config_name, '.csv'), header=FALSE)

    # If a covariate .csv file exists, use that instead
    if(!is.null(covs_name)) {
      # Grab fixed effects & measures (and gbd fixed effects) from CSV if present
      covs <- fread(paste0(repo, indicator_group, '/', covs_name, '.csv'), header = TRUE)
      covs <- subset(covs, include == T) # Use only those where include flag set to true
      fe <- subset(covs, gbd == F)
      gbd <- subset(covs, gbd == T)
      fixed_effects <- paste0(fe$covariate, collapse = " + ")
      measures <- paste0(fe$measure, collapse = " + ")
      gbd_fixed_effects <- paste0(gbd$covariate, collapse = " + ")

      # remove any other versions from original config
      config <- subset(config, !(V1 %in% c("fixed_effects", "fixed_effects_measures", "gbd_fixed_effects")))
      config <- config %>%
                rbind(., list("fixed_effects", fixed_effects)) %>%
                rbind(., list("fixed_effects_measures", measures)) %>%
                rbind(., list("gbd_fixed_effects", gbd_fixed_effects))
    }

    # Assign all the covariates to the environment
    for(param in config[, V1]) {
      assign(param, config[V1==param, V2], envir=globalenv())
    }

    return(config)

  }

## Create directory structure
#   Arguments:
#     indicator_group = Category of indicator, i.e. "education"
#     indicator       = Specific outcome to be modeled within indicator category, i.e. "edu_0"
  create_dirs <- function(indicator_group, indicator) {

    dir.create(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group))
    dir.create(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator))

    indicator_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator)

    for(dir in c('output','model_image_history')) {
      dir.create(paste0(indicator_dir,'/',dir), showWarnings = FALSE)
    }

  }

## Make template rasters (pull from central analysis folders managed by Lucas based off the area to model specified in config)
#   Arguments:
#     simple = Single polygon that defines boundaries of the entire area you want to model over.
#   Returns: Empty raster over modeling area. To be used for cropping covariates quickly and projecting model.
  get_template_raster <- function(simple) {

    message('Creating rasters of admin units')
    root <- ifelse(Sys.info()[1]=="Windows", "<<<<< FILEPATH REDACTED >>>>>", "<<<<< FILEPATH REDACTED >>>>>")

    # Centrally controlled folder of analysis shapefiles
    analysis_raster_dir <- paste0(root,'<<<< FILEPATH REDACTED >>>>>')

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
   load_input_data <- function(indicator, simple = NULL, agebin = 0, removeyemen = FALSE, pathaddin = "", withdate=FALSE, date='', years='five_year',range=5, update_run_date = FALSE, withtag=FALSE, datatag='', use_share=FALSE, year_list = c(2000:2015)) {

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
    root <- ifelse(Sys.info()[1]=="Windows", "<<<< FILEPATH REDACTED >>>>>", "<<<< FILEPATH REDACTED >>>>>")
    if(use_share==FALSE) load_dir <- paste0(root,'<<<< FILEPATH REDACTED >>>>>')
    if(use_share==TRUE) load_dir <- '<<<< FILEPATH REDACTED >>>>>'
    if(!withdate & !withtag) d = read.csv(paste0(load_dir, indicator, '.csv'))
    if(withtag) d = read.csv(paste0(load_dir, indicator, datatag, '.csv'))
    if(withdate) d = read.csv(paste0(root,'<<<< FILEPATH REDACTED >>>>>',rd,'/', indicator, '.csv'))
    d$latitude  = as.numeric(as.character(d$latitude))
    d$longitude = as.numeric(as.character(d$longitude))
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

    if(agebin!=0) d = d[age==agebin,]
    if(removeyemen) d = d[country!='Yemen' & country!='YEM',]

    if(years=='five_year') {
      d <- d[year >= 1998 & year <= 2002, year := 2000]
      d <- d[year >= 2003 & year <= 2007, year := 2005]
      d <- d[year >= 2008 & year <= 2012, year := 2010]
      d <- d[year >= 2013 & year <= 2017, year := 2015]
    }
    if(years=='annual') {
      d <- d[year >= 1998 & year <= 2000, year := 2000]
    }

    d <- subset(d, year >= 1998)
    if (nrow(subset(d, year > max(year_list))) > 0) {
      warning(paste0("Dropping all data after max(year_list) = ", max(year_list), "..."))
      d <- subset(d, year <= max(year_list))
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
      if(dir.exists(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, "/output/", run_date)) == TRUE) {
        existing_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, "/output/", run_date)
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
      if(dir.exists(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, "/output/", run_date)) == FALSE) {
        run_date_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, "/output/", run_date)
        dir.create(run_date_dir, showWarnings = FALSE)
      }
      write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
      return(list(d, run_date))
    }

    if(update_run_date == FALSE) {
      if(agebin==0){
        run_date_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, "/output/", run_date)
        dir.create(run_date_dir, showWarnings = FALSE)
        write.csv(d, file=paste0(run_date_dir, "/input_data", pathaddin, ".csv"))
      } else {
        run_date_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, "_age",agebin,"/output/", run_date)
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

load_simple_polygon <- function(gaul_list, buffer, tolerance = 0.4,
                                subset_only = FALSE,returnsubsetshape=FALSE,
                                makeplots = F, use_premade = T) {
  if (use_premade == T) {
    # Check for common gaul lists so we can query premade simple polygons
    if(identical(gaul_list, c(29, 42, 45, 47, 50, 66, 90, 94, 106, 105, 144, 155, 159, 181, 182, 214, 217, 221, 243))) {
      message("Your GAUL list matches WSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>>/simple_polygons/wssa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(8,49,59,68,76,89))) {
      message("Your GAUL list matches CSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>>/simple_polygons/cssa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(43,58,70,77,79,133,150,152,170,205,226,74,257,253,270))) {
      message("Your GAUL list matches ESSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>>/simple_polygons/essa.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(4,40762,40765,145,169,6,248))) {
      message("Your GAUL list matches NAME, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>>/simple_polygons/name.RData')
      assign("subset_shape", subset_shape, envir=globalenv())
      if(makeplots) plot(spoly_spdf)
      if(makeplots) plot(subset_shape, add = TRUE, border = grey(0.5))
      return(list(subset_shape=subset_shape,spoly_spdf=spoly_spdf))
    }
    if(identical(gaul_list, c(35,142,172,227,235,271))) {
      message("Your GAUL list matches SSSA, loading premade simple_polygon and assigning subset_shape to global environment...")
      load('<<<< FILEPATH REDACTED >>>>>/simple_polygons/sssa.RData')
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
      load('<<<< FILEPATH REDACTED >>>>>/simple_polygons/africa.RData')
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
  message("Opening master shapefile...")
  master_shape <- readRDS('<<<< FILEPATH REDACTED >>>>>/master_shape_all.rds')
  subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% gaul_list, ]

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

      intersect <- gIntersection(the_poly, compare_to)

      if(is.null(intersect)) {
        return(FALSE)
      } else {
        return(ifelse((area(intersect) == area(the_poly)), TRUE, FALSE))
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

## Create spatial mesh (pull from central analysis folders managed by Lucas based off the area to model specified in config)
#   d = Data table with columns "latitude" and "longitude".
#   simple = Single polygon to create mesh over.
#   max_edge = Largest allowed triangle edge length (see INLA's inla.mesh.2d documentation).
#   mesh_offset = One or two values, for an inner and an optional outer extension distance (see INLA's inla.mesh.2d documentation).
build_space_mesh <- function(d, simple, max_edge, mesh_offset) {
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

    plot(mesh_s,asp=1);points(d$longitude,d$latitude,col=d$year)

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

## Load list of raster inputs (pop and simple)
build_simple_raster_pop <- function(subset_shape,u5m=FALSE) {

  if(u5m==FALSE){
    master_pop <- brick('<<<< FILEPATH REDACTED >>>>>/WorldPop_total_global_stack.tif') #WorldPop_allStages_stack.tif')
  } else {
    master_pop <- brick(raster('<<<< FILEPATH REDACTED >>>>>/worldpop_a0004t_5y_2000_00_00.tif'),
    raster('<<<< FILEPATH REDACTED >>>>>/worldpop_a0004t_5y_2005_00_00.tif'),
    raster('<<<< FILEPATH REDACTED >>>>>/worldpop_a0004t_5y_2010_00_00.tif'),
    raster('<<<< FILEPATH REDACTED >>>>>/worldpop_a0004t_5y_2015_00_00.tif'))
  }

  cropped_pop <- crop(master_pop, extent(subset_shape), snap="out")
  ## Fix rasterize
  initial_raster <- rasterize(subset_shape, cropped_pop, field = 'GAUL_CODE')
  if(length(subset_shape[!subset_shape$GAUL_CODE%in%unique(initial_raster),])!=0) {
    rasterized_shape <- merge(rasterize(subset_shape[!subset_shape$GAUL_CODE%in%unique(initial_raster),], cropped_pop, field = 'GAUL_CODE'), initial_raster)
  }
  if(length(subset_shape[!subset_shape$GAUL_CODE%in%unique(initial_raster),])==0) {
    rasterized_shape <- initial_raster
  }
  #rasterized_shape <- rasterize(subset_shape, cropped_pop, field='GAUL_CODE')
  masked_pop <- mask(x=cropped_pop, mask=rasterized_shape)

  raster_list <- list()
  raster_list[['simple_raster']] <- rasterized_shape
  raster_list[['pop_raster']] <- masked_pop

  return(raster_list)

}

## Grab lists of GAUL codes
get_gaul_codes <- function(gauls) {

  # Inputs: gauls = vector of continents, regions, stages, or countries
  #                 **now also accepting combinations**
  # Outputs: vector of gaul codes

  gaul_list <- NULL

  # Set up list of gaul_codes
   gaul_ref <- list()
   gaul_ref[["eastern_europe"]] <- c(165,254,26,147,140,78)
   gaul_ref[["middle_east"]] <- c(117,1,187,137,267,141,201,215,249,118,130,238,21,255)
   gaul_ref[["latin_america"]] <- c(71,72,11,63,99,123,246,20,209,211,24,108,200,258,191,75,
                                    180,162,111,103,61,28,30,57,107,73,194,233,263,195,37,33)
   gaul_ref[["se_asia"]] <- c(147296,67,147295,171,44,264,139,153,240,196,116,242)
   gaul_ref[["south_asia"]] <- c(23,31,115,175,188)
   gaul_ref[["africa"]] <- c(4,6,8,29,35,42,43,45,47,49,50,58,59,66,68,70,
                             74,76,77,79,89,90,94,95,105,106,142,144,145,150,
                             152,155,159,169,170,172,181,182,205,214,217,221,
                             226,235,243,248,253,268,270,271,40762,40765,
                             227,257,133,269)
   gaul_ref[["central_america"]] <- c(11,20,28,30,24,61,63,71,72,99,103,111,108,
                                      123,209,162,180,191,200,75,246,211,258)
   gaul_ref[["south_america"]] <- c(33,37,57,73,107,195,194,233,263)
   gaul_ref[["stage1"]] <- c(170,152,205,226,270,133,150,74,257,253,271,79,77,70,
                             58,43,89,76,50,214,59,68,45,49,8,169,6,40765,4,172,235,
                             142,227,35,243,221,217,182,181,159,155,144,105,90,106,94,
                             47,66,42,29,72,108,180,111,103,195,33,239,132,167,67,171,
                             44,264,139,153,240,196,116,242,175,188,115,117,231,31,23,
                             1,267,269,91,118,130,238,3,192,262,225,157,135)
   gaul_ref[["stage2"]] <- c(145,75,162,28,57,107,73,261,138,147295,154,187,215,249,19,
                             255,254,34,170,152,205,226,270,133,150,74,257,253,271,79,77,
                             70,58,43,89,76,50,214,59,68,45,49,8,169,6,40765,4,172,235,
                             142,227,35,243,221,217,182,181,159,155,144,105,90,106,94,
                             47,66,42,29,72,108,180,111,103,195,33,239,132,167,67,171,
                             44,264,139,153,240,196,116,242,175,188,115,117,231,31,23,
                             1,267,269,91,118,130,238,3,192,262,225,157,135)
   gaul_ref[["cssa"]]   <- c(8,49,59,68,76,89)
   gaul_ref[["cssa_diarrhea"]]   <- c(49,59,68,76,89) # Removed Angola
   gaul_ref[["sssa_diarrhea"]]   <- c(8,35,142,172,227,271) # Added Angola, removed Swaziland
   gaul_ref[["essa_diarrhea"]]   <- c(43,58,70,77,79,133,150,152,170,205,226,74,235,257,253,270) # Added Swaziland

   gaul_ref[["cssa_diarrhea2"]]   <- c(49,45,50,59,89,76)
   gaul_ref[["essa_diarrhea2"]]   <- c(43,58,70,77,79,133,150,152,170,205,226,74,235,257,253,270,6) # Added Swaziland
   gaul_ref[["name_diarrhea2"]]   <- c(4,40762,145,169,248) # no yemen for now,269)
   gaul_ref[["sssa_diarrhea2"]]   <- c(8,35,142,172,227,271) # Added Angola, removed Swaziland
   gaul_ref[["wssa_diarrhea2"]]   <- c(29,42,47,66,90,94,106,105,144,155,159,181,182,214,217,221,243)

   gaul_ref[["essa_edu"]]   <- c(226,79,74,6,77,70)
   gaul_ref[["cssa_edu"]]   <- c(182,43,58,133,150,152,170,205,235,257,253,270,49,45,59,89,76,68,8) # Added Swaziland
   gaul_ref[["name_edu"]]   <- c(4,40762,145,169,248) # no yemen for now,269)
   gaul_ref[["sssa_edu"]]   <- c(35,142,172,227,271) # Added Angola, removed Swaziland
   gaul_ref[["wssa_edu"]]   <- c(29,42,47,66,90,94,106,105,144,155,159,181,214,217,221,243,50)

   gaul_ref[["essa_edu2"]]   <- c(43,58,70,77,79,133,150,152,170,205,226,74,257,253,270)
   gaul_ref[["name_edu2"]]   <- c(4,40762,145,169,6,248) # no yemen for now,269)
   gaul_ref[["sssa_edu2"]]   <- c(35,142,172,227,235,271)
   gaul_ref[["wssa_edu2"]]   <- c(29,42,45,47,50,66,90,94,106,105,144,155,159,181,182,214,217,221,243)
   gaul_ref[["cssa_edu2"]]   <- c(8,49,59,68,76,89)

   gaul_ref[["essa"]]   <- c(43,58,70,77,79,133,150,152,170,205,226,74,257,253,270)
   gaul_ref[["name"]]   <- c(4,40762,40765,145,169,6,248) # no yemen for now,269)
   gaul_ref[["namelite"]] <- c(169,6,40765) # remove tunisia, algeria and libya from mesh
   gaul_ref[["sssa"]]   <- c(35,142,172,227,235,271)
   gaul_ref[["wssa"]]   <- c(29,42,45,47,50,66,90,94,106,105,144,155,159,181,182,214,217,221,243)
   gaul_ref[["cessa"]]  <- c(8,49,59,68,76,89,43,58,70,77,79,133,150,152,170,205,226,74,257,253,270)
   gaul_ref[["cwssa"]]  <- c(8,49,59,68,76,89,29,42,45,47,50,66,90,94,106,105,144,155,159,181,182,214,217,221,243)

   gaul_ref[["cessa2"]]  <- c(49,59,68,76,89,43,58,70,77,79,133,150,205,226,74,257,253,5) # remove AGO, ZMB, MWI, MOZ, add cameroon
   gaul_ref[["sssa2"]]    <- c(35,142,172,227,235,271,8,270,152,170) # add AGO, ZMB, MWI, MOZ

   gaul_ref[["cssa_cam"]]   <- c(8,49,59,68,76,89,45) # added cameroon
   gaul_ref[["wssa_nocam"]]   <- c(29,42,47,50,66,90,94,106,105,144,155,159,181,182,214,217,221,243) # removed cameroon

   gaul_ref[["name_hi"]] <- c(4,40765,145,169,248,40762) # No Sudan
   gaul_ref[["essa_hi"]] <- c(133,270) # Zambia and Kenya
   gaul_ref[["essa_lo"]] <- c(43,58,70,77,79,150,152,170,205,226,74,257,253,6) # Contains Sudan
   gaul_ref[["cssa_hi"]] <- c(59,76,89) # GNQ, COG, and Gabon
   gaul_ref[["cssa_lo"]] <- c(8,49,68)
   gaul_ref[["wssa_hi"]] <- c(94) #Ghana
   gaul_ref[["wssa_lo"]] <- c(29,42,45,47,50,66,90,106,105,144,155,181,182,214,217,221,243) # no Ghana
   gaul_ref[["sssa_hi"]] <- c(35,172,227) # Botswana, namibia, south africa
   gaul_ref[["sssa_lo"]] <- c(142,235,271) # LSO, Zimbabwe, Swaziland
   gaul_ref[["mwi"]] <- 152
   gaul_ref[["nga"]] <- 182
   gaul_ref[["egy"]] <- 40765
   gaul_ref[["gha"]] <- 94
   gaul_ref[["cod"]] <- 68
   gaul_ref[["zaf"]] <- 227
   gaul_ref[["ago"]] <- 8
   gaul_ref[["caf"]] <- 49
   gaul_ref[["cog"]] <- 59
   gaul_ref[["cod"]] <- 68
   gaul_ref[["gnq"]] <- 76
   gaul_ref[["gab"]] <- 89
   gaul_ref[["bwa"]] <- 35
   gaul_ref[["lso"]] <- 142
   gaul_ref[["nam"]] <- 172
   gaul_ref[["zaf"]] <- 227
   gaul_ref[["swz"]] <- 235
   gaul_ref[["zwe"]] <- 271
   gaul_ref[["ben"]] <- 29
   gaul_ref[["bfa"]] <- 42
   gaul_ref[["cmr"]] <- 45
   gaul_ref[["cpv"]] <- 47
   gaul_ref[["tcd"]] <- 50
   gaul_ref[["civ"]] <- 66
   gaul_ref[["gmb"]] <- 90
   gaul_ref[["gha"]] <- 94
   gaul_ref[["gin"]] <- 106
   gaul_ref[["gnb"]] <- 105
   gaul_ref[["lbr"]] <- 144
   gaul_ref[["mli"]] <- 155
   gaul_ref[["mrt"]] <- 159
   gaul_ref[["ner"]] <- 181
   gaul_ref[["nga"]] <- 182
   gaul_ref[["stp"]] <- 214
   gaul_ref[["sen"]] <- 217
   gaul_ref[["sle"]] <- 221
   gaul_ref[["tgo"]] <- 243
   gaul_ref[["bdi"]] <- 43
   gaul_ref[["com"]] <- 58
   gaul_ref[["dji"]] <- 70
   gaul_ref[["eri"]] <- 77
   gaul_ref[["eth"]] <- 79
   gaul_ref[["ken"]] <- 133
   gaul_ref[["mdg"]] <- 150
   gaul_ref[["mwi"]] <- 152
   gaul_ref[["moz"]] <- 170
   gaul_ref[["rwa"]] <- 205
   gaul_ref[["som"]] <- 226
   gaul_ref[["ssd"]] <- 74
   gaul_ref[["tza"]] <- 257
   gaul_ref[["uga"]] <- 253
   gaul_ref[["zmb"]] <- 270
   gaul_ref[["dza"]] <- 4
   gaul_ref[["egy"]] <- 40765
   gaul_ref[["lby"]] <- 145
   gaul_ref[["mar"]] <- 169
   gaul_ref[["sdn"]] <- 6
   gaul_ref[["tun"]] <- 248
   gaul_ref[["essa_hilo"]] <- c(133,270,43,58,70,77,79,150,152,170,205,226,74,257,253,6,142,235,271) # Added LSO, ZWE, swaziland, and SDN

  # Start to build gaul list by matching reference list above
   gaul_list <- lapply(gauls, function(x) {gaul_ref[[x]]})
   names(gaul_list) <- gauls

  # Try to bring in additional gaul codes with gaul_convert
   gaul_list <- lapply(names(gaul_list), function(x) {

     gl <- c()
     if (is.null(gaul_list[[x]]) == T) {
       try(gl <- gaul_convert(x, from = "iso3"))
     } else {
       gl <- gaul_list[[x]]
     }

     # Message if this didn't work
     if (is.null(gl) == T) {
       message(paste0("/nUnable to find a match for the following: ", x))
       message("Check your input vector to ensure valid region names / country codes")
     }
     return(gl)
   })

  # Rearrange gaul_list to be a single vector
   gaul_list <- unlist(gaul_list)

  # Check for duplicates & drop if needed
   if (length(gaul_list) != length(unique(gaul_list))) {
     message("Duplicate gaul codes found - your input regions may overlap.")
     message("Dropping duplicates...")
   }
  gaul_list <- unique(gaul_list)
  return(gaul_list)

}

get_gaul_codes_subnat <- function(gaul_list, admin_level) {

  if(admin_level %in% c(0,1)) hierarchy <- read.dbf(paste0("<<<< FILEPATH REDACTED >>>>>admin", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, ".dbf"))
  if(admin_level == 2) hierarchy <- read.dbf(paste0("<<<< FILEPATH REDACTED >>>>>/admin", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, "_modified.dbf"))
  hierarchy <- as.data.table(hierarchy)
  adm0_name <- grep("ADM0_CODE", names(hierarchy), value = T)
  names(hierarchy)[names(hierarchy)==adm0_name] <- "ADM0_CODE"
  hierarchy <- hierarchy[ADM0_CODE %in% gaul_list,]
  adm_specific_name <- grep(paste0('ADM', admin_level, '_CODE'), names(hierarchy), value = T)
  adm_list <- unique(hierarchy[[adm_specific_name]])
  return(adm_list)

}


add_gauls_regions <- function(df, simple_raster) {

  # Extract GAUL_CODE from simple_raster using lat/longs
  df$GAUL_CODE <- raster::extract(simple_raster, df[ , c('longitude', 'latitude'), with=F])

  # Add names of regions by GAUL_CODE
  for(r in Regions) {
    df <- df[GAUL_CODE %in% get_gaul_codes(r), region := r]
  }

  ##Check that merge of GAUL_CODE worked properly
  df <- df[!is.na(region), good_records := 1]
  message(paste0(length(df$good_records) - length(df[good_records==1, N]), ' out of ', length(df$good_records), ' could not have GAUL/region extracted properly. Probably coastal points, need to fix.'))
  df <- df[good_records==1, ]

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

  if (Sys.info()["sysname"] == "Linux") {
    j_root <- "<<<< FILEPATH REDACTED >>>>>"
  } else {
    j_root <- "<<<< FILEPATH REDACTED >>>>>"
  }

  str_match <- stringr::str_match

  # Catch if already passed gaul codes
  if(class(countries) =="numeric") return (countries)
  if(all(grepl("^[[:digit:]]+$", countries))) return(countries)

  table_file <- paste0(j_root, "<<<< FILEPATH REDACTED >>>>>/gaul_to_loc_id.csv")
  gaul_table <- read.csv(table_file) %>% data.table

  # convert input & output to lower case for easier matching
  #lowercase_cols <- c("short_name", "official_name", "iso3", "iso2", "uni", "undp")
  #gaul_table[, (lowercase_cols) := lapply(.SD, tolower), .SDcols = lowercase_cols,]

  ## ## convert input & output to lower case for easier matching
  if (verbose == T) message("\nlowercasing columns. the columns that get lowered are:")
  for(i in 1:ncol(gaul_table)){
    if(class(gaul_table[[i]]) == "factor"){
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
  gaul_to_loc_id <- fread("<<<< FILEPATH REDACTED >>>>>/gaul_to_loc_id.csv")
  d <- d[, ihme_lc_id := as.character(country)]

  d <- merge(d, gaul_to_loc_id, by='ihme_lc_id',all.x=T)

  message(nrow(d))
  for(r in re) {
    d <- d[GAUL_CODE %in% get_gaul_codes(r), region := r]
  }
  return(d)
}

#' set_root
#'
#' @description A function to set the "root" based on the OS used.
#' @author Rebecca Stubbs
#'
set_root<-function(){
  root<<-ifelse(Sys.info()[1]=="Windows", "<<<< FILEPATH REDACTED >>>>>", "<<<< FILEPATH REDACTED >>>>>") # Setting "root" as global variable
}

#' load_libs
#'
#' @description A function to load libraries based on whether the
#' code is being run on cluster prod or the geos nodes. Tries
#' to load each package, and returns errors or warnings if the package fails to
#' load rather than breaking the loop.
#' @author Rebecca Stubbs
#'
#' @param packages A vector or list of packages
#' @param stop Boolean; if T will stop code if not all packages load. if F (default),
#' the code will contiue and will throw a warning.
load_libs<-function(packages,stop=F){
  set_root() # Seting the "root" based on OS
  package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                        paste0(root,'<<<< FILEPATH REDACTED >>>>>/geos_packages'),
                        paste0(root,'<<<< FILEPATH REDACTED >>>>>/packages'))
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
#' @author Rebecca Stubbs
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
