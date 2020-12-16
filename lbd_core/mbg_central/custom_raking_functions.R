#################################################################################
### Custom raking function
## builds a new simple polygon and raster, crosswalks the cell pred object to match the new raster. Generates population weights based on the shapefile and field passed in, calculates raking factors using the rake_to df
## Inputs: cell_pred - cell pred object to be raked
##         shapefile path - path to shapefile that will be used for raking
##         field - field in shapefile that has admin identifiers that match with rake_to
##         rake_to - df with name, year, and mean columns. values in name must match up with values in field
##         region - region used to produce cell pred object
##         year_list - years
## Outputs: List with new cell pred, simple raster, raking factors, summary rasters


#' @title Custom Raking Function
#' @description Builds a new simple polygon and raster, crosswalks the cell pred object to match the new raster. 
#'  Generates population weights based on the shapefile and field passed in, 
#'  calculates raking factors using the rake_to df
#'
#' @param cell_pred Cell pred object to be raked.
#' @param shapefile_path Path to shapefile that will be used for
#'   raking.
#' @param field Field in shapefile that has admin identifiers that
#'   match with rake_to.
#' @param rake_to Df with name, year, and mean columns. values in name
#'   must match up with values in field.
#' @param reg Region used to produce cell pred object.
#' @param year_list List of years
#' @param modeling_shapefile_version string indicating which version of
#'   shapefile to use for matching cell_pred pixels to adm codes
#'
#' @return Returns a new cell pred object, simple raster, raking factors, pre-raking aggregate numbers, and rasters of mean, lower, upper, and cirange for years in year list
#' @export

custom_rake <- function(cell_pred, shapefile_path, field, rake_to, reg, year_list,
                        modeling_shapefile_version = 'current') {

  #initialize results list
  outputlist <- list()

  #Calculate interval month, number of months between years in year list
  year_diff <- diff(year_list)
  if(length(unique(year_diff))!=1) {
    stop("Please use annual or 5-year intervals exclusively in year_list")
  } else {
    interval_mo <- year_diff[[1]] * 12
  }

  ## get simple polygon and simple raster used to produce cell pred
  message('Loading simple polygon')
  simple_polygon <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = modeling_shapefile_version),
                                        buffer = 0.4,
                                        shapefile_version = modeling_shapefile_version)
  subset_shape   <- simple_polygon[['subset_shape']]
  simple_polygon <- simple_polygon[['spoly_spdf']]

  message('Loading simple raster')
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

  #get new simple polygon and simple raster that will be raked to
  message('Loading new simple polygon to be raked to')
  new_simple_polygon <- load_simple_polygon(gaul_list = get_adm0_codes(reg), buffer = 0.4, custom_shapefile_path = shapefile_path)
  new_subset_shape   <- new_simple_polygon[['subset_shape']]
  new_simple_polygon <- new_simple_polygon[['spoly_spdf']]

  message('Loading new simple raster to be raked to')
  new_raster_list    <- build_simple_raster_pop(new_subset_shape)
  new_simple_raster  <- new_raster_list[['simple_raster']]
  new_pop_raster     <- new_raster_list[['pop_raster']]

  #get extents of original and simple raster to line up - extend and crop just in case
  new_simple_raster <- extend(new_simple_raster, simple_raster,values=NA)
  new_simple_raster <- crop(new_simple_raster, extent(simple_raster))

  #check original and new simple rasters match in extent and resolution
  if(extent(new_simple_raster) != extent(simple_raster)){
    stop("new simple raster extent does not match original simple raster")
  }
  if(any(res(new_simple_raster) != res(simple_raster))){
    stop("new simple raster resolution does not match original simple raster")
  }

  #crosswalk cell_pred to new raster and rename outputs
  message('Crosswalking cell pred object to new raster')
  new_pred_object <- crosswalk_cell_pred_add_NA(simple_raster, new_simple_raster, cell_pred, year_list)
  cell_pred <- new_pred_object[[1]]
  simple_raster <- new_pred_object[[2]]

  message('Getting population raster')
  ## Pull 2000-2015 annual population brick using new covariates function
  if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
  pop_raster_annual <- load_worldpop_covariate(template_raster = simple_polygon,
                                               pop_measure = pop_measure,
                                               pop_release = pop_release,
                                               start_year = min(year_list),
                                               end_year = max(year_list),
                                               interval = interval_mo)[[1]]

  ## extend and crop pop raster to ensure it matches the simple raster
  pop_raster_annual  <- extend(pop_raster_annual, simple_raster, values=NA)
  pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
  pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
  pop_raster_annual  <- mask(pop_raster_annual, simple_raster)

  ## check to ensure the pop raster matches the simple raster in extent and resolution
  if(extent(pop_raster_annual) != extent(simple_raster)){
    stop("population raster extent does not match simple raster")
  }
  if(any(res(pop_raster_annual) != res(simple_raster))){
    stop("population raster resolution does not match simple raster")
  }

  #get custom admin raster
  custom_admin_raster <- load_custom_admin_raster(shapefile_path, field, simple_raster)

  ## Create population weights using the annual brick and feed custom year argument to aggregation function
  message('Building population weights object')
  pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                          simple_raster = simple_raster,
                                          pop_raster    = pop_raster_annual,
                                          gaul_list     = get_adm0_codes(reg),
                                          custom_admin_raster = custom_admin_raster)

  message('Making condSim')
  cond_sim_draw_adm0 <- make_condSim(admin_level    = 0,
                                     pop_wts_object = pop_wts_adm0,
                                     cell_pred      = cell_pred,
                                     gaul_list      = get_adm0_codes(reg),
                                     summarize      = FALSE,
                                     years          = year_list)

  # save some raw country estimates to compare with GBD in a plot later on
  cond_sim_raw_adm0 <- apply(cond_sim_draw_adm0   , 1, mean)
  adm0_geo          <- cbind(mean=cond_sim_raw_adm0,
                             lower=apply(cond_sim_draw_adm0, 1, quantile, probs=.025),
                             upper=apply(cond_sim_draw_adm0, 1, quantile, probs=.975))
  outputlist[['adm0_geo']] <- data.table(split_geo_names(adm0_geo),adm0_geo)

  message('Calculating raking factors')
  rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
                              rake_to     = rake_to)

  # rake the cell preds
  message('Raking cell pred object')
  raked_cell_pred <- rake_predictions(raking_factors = rf,
                                      pop_wts_object = pop_wts_adm0,
                                      cell_pred      = cell_pred)

  outputlist[["raked_cell_pred"]] <- raked_cell_pred
  outputlist[["simple_raster"]] <- simple_raster
  outputlist[["raking_factors"]] <- rf

  message('Summarizing raked cell preds')
  for(summeasure in c('mean','cirange','lower','upper')){
    message(sprintf('    %s',summeasure))
    outputlist[[sprintf('%s_raked_raster', summeasure)]] <-
      make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                              mask                 = simple_raster,
                              return_as_raster     = TRUE,
                              summary_stat         = summeasure)
  }

  return(outputlist)
}


#' @title Save outputs from custom raking function
#' @description This function saves outputs from FILEPATH
#'
#' @param custom_rake_output Output list from custom raking function
#' @param outdir Directory for files to be saved to
#' @param indicator Name of indicator being modeled
#' @param age_group Name of age group
#' @param prefix Character string to be added to files to avoid overwriting non-custom files in folder
#' @param reg Name of region that was modeled
#' @return NULL
#' @export
save_custom_raking_outputs <- function(custom_rake_output,
                                       outdir,
                                       indicator,
                                       prefix,
                                       reg) {

  ol <- custom_rake_output
  raked_cell_pred <- ol[["raked_cell_pred"]]
  save(raked_cell_pred, file=sprintf('%s/%s_%s_raked_cell_draws_eb_bin0_%s_0.RData',outdir, prefix, indicator, reg))
  simple_raster <- ol[["simple_raster"]]
  writeRaster(simple_raster, file= sprintf('%s/%s_%s_simple_raster', outdir, prefix, indicator), format = "GTiff")
  raking_factors <- ol[["raking_factors"]]
  write.csv(raking_factors, file=sprintf('%s/%s_%s_%s_rf.csv',outdir, prefix, indicator, reg))
  adm0_geo <- ol[["adm0_geo"]]
  write.csv(adm0_geo, file = sprintf('%s/%s_%s_adm0_geo.csv', outdir, prefix, indicator))
  mean_raster <- ol[["mean_raked_raster"]]
  writeRaster(mean_raster, file= sprintf('%s/%s_%s_mean_raked_2000_2015', outdir, prefix, indicator), format = "GTiff")
  cirange_raster <- ol[["cirange_raked_raster"]]
  writeRaster(cirange_raster, file= sprintf('%s/%s_%s_cirange_raked_2000_2015', outdir, prefix, indicator), format = "GTiff")
  upper_raster <- ol[["upper_raked_raster"]]
  writeRaster(upper_raster, file= sprintf('%s/%s_%s_upper_raked_2000_2015', outdir, prefix, indicator), format = "GTiff")
  lower_raster <- ol[["lower_raked_raster"]]
  writeRaster(lower_raster, file= sprintf('%s/%s_%s_lower_raked_2000_2015', outdir, prefix, indicator), format = "GTiff")
}

#' Custom Aggregation Function
#' @title Custom Aggregation Function
#' @description Create custom aggregates
#'
#' @param cell_pred Cell pred object to be aggregated
#' @param custom_shapefile_path Path to shapefile that will be used for aggregation
#' @param custom_shapefile Alternatively can specify a shapefile object directly
#' @param field Field in shapefile that should be used for aggregation (ADM2_CODE, etc).
#' @param reg Region used to produce cell pred object.
#' @param yl List of years
#' @param ss Summary statistics
#' @param return_shapefile Boolean, should the function return the shapefile?
#' @param verbose Should the function be talkative?
#' @param modeling_shapefile_version string used to specify which shapefile version to use for merging adm codes to pixels. should be the version used to generate the cell_preds
#' 
#' @return Returns a list containing:
#'          \code{draws}: data.table with rows for each aggregated area & year,
#'                   and columns representing [field], year, and draws[V1...Vn]
#'          \code{shapefile_data}: contents of the shapefile defined by \code{shapefile_path},
#'                            which can be merged on to "draws" to add on admin names, etc
#'          \code{summarized}: data table with row for each aggregated area & year,
#'                        and columns for [field], year, and all summary stats in \code{ss}
#'                        (if \code{ss}=NULL, then this is also NULL)
#'          \code{shapefile}: the shapefile itself (if return_shapefile = T; otherwise will be NULL)
#' @export
custom_aggregate <- function(cell_pred,
                             shapefile_path,
                             field,
                             reg,
                             yl,
                             pop_meas = pop_measure,
                             ss = NULL,
                             return_shapefile = T,
                             verbose = T,
                             modeling_shapefile_version = 'current') {

  #initialize results list
  outputlist <- list()

  ## get simple polygon and simple raster used to produce cell pred
  if (verbose) message('Loading simple polygon')
  simple_polygon <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapfile_version = modeling_shapefile_version),
                                        buffer = 0.4, shapefile_version = modeling_shapefile_version)
  subset_shape   <- simple_polygon[['subset_shape']]
  simple_polygon <- simple_polygon[['spoly_spdf']]

  if (verbose) message('Loading simple raster')
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

  #get new simple polygon and simple raster that will be raked to
  if (verbose) message('Loading new simple polygon to be raked to')

  new_subset_shape   <- new_simple_polygon[['subset_shape']]
  new_simple_polygon <- new_simple_polygon[['spoly_spdf']]

  if (verbose)  message('Loading new simple raster to be raked to')
  new_raster_list    <- build_simple_raster_pop(new_subset_shape, field=field)
  new_simple_raster  <- new_raster_list[['simple_raster']]
  new_pop_raster     <- new_raster_list[['pop_raster']]

  #get extents of original and simple raster to line up - extend and crop just in case
  new_simple_raster <- extend(new_simple_raster, simple_raster,values=NA)
  new_simple_raster <- crop(new_simple_raster, extent(simple_raster))

  #crosswalk cell_pred to new raster and rename outputs
  if (verbose) message('Crosswalking cell pred object to new raster')
  new_pred_object <- crosswalk_cell_pred_add_NA(simple_raster, new_simple_raster, cell_pred, length(yl))
  cell_pred <- new_pred_object[[1]]
  simple_raster <- new_pred_object[[2]]

  # if interval month doesn't exist, use year list to determine 5-year or annual
  if(!exists("interval_mo")){
    if((yl[[2]] - yl[[1]]) %% 5 == 0) {
      interval_mo <- 60
    } else {
      interval_mo <- 12
    }
  }

  if (verbose) message('Getting population raster')
  ## Pull annual population brick using new covariates function
  if (class(yl) == "character") yl <- eval(parse(text=yl))
  pop_raster_annual <- load_worldpop_covariate(template_raster = simple_polygon,
                                               pop_measure = pop_meas,
                                               pop_release = pop_release,
                                               start_year = min(yl),
                                               end_year = max(yl),
                                               interval = as.numeric(interval_mo))[[1]]

  pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
  pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
  pop_raster_annual  <- mask(pop_raster_annual, simple_raster)

  #get custom admin raster
  custom_admin_raster <- load_custom_admin_raster(shapefile_path, field, simple_raster)

  ## Create population weights using the annual brick and feed custom year argument to aggregation function
  if (verbose) message('Building population weights object')
  pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                          simple_raster = simple_raster,
                                          pop_raster    = pop_raster_annual,
                                          gaul_list     = get_adm0_codes(reg),
                                          custom_admin_raster = custom_admin_raster)

  if (verbose) message('Making condSim')
  cond_sim_draw_adm0 <- make_condSim(admin_level    = 0,
                                     pop_wts_object = pop_wts_adm0,
                                     cell_pred      = cell_pred,
                                     gaul_list      = get_adm0_codes(reg),
                                     summarize      = FALSE,
                                     years          = yl)

  # Add ID and year on to this admin level draw object
  admin_draws <- data.table(split_geo_names(as.matrix(cond_sim_draw_adm0)),cond_sim_draw_adm0)
  admin_draws$year <- as.numeric(admin_draws$year)
  setnames(admin_draws, "name", field)

  if (verbose) message("Storing admin draws")
  outputlist[["draws"]] <- admin_draws

  if (verbose) message("Storing shapefile information")

  if (!is.null(custom_shapefile_path)) the_shapefile <- readOGR(custom_shapefile_path)
  if (!is.null(custom_shapefile)) the_shapefile <- custom_shapefile

  outputlist[["shapefile_data"]] <- as.data.table(the_shapefile)

  # Calculate summary statistics
  if (!is.null(ss)) {
    if (verbose) message("Calculating and storing summary statistics")
    ss_df <- lapply(ss, function(summstat) {
      draw_cols <- names(admin_draws)[grepl("V[0-9]+", names(admin_draws))]
      out <- admin_draws[, .(stat = apply(.SD, 1, summstat)), by = list(get(field), year), .SDcols = draw_cols]
      setnames(out, c("stat", "get"), c(summstat, field))
      return(out)
    })

    outputlist[["summarized"]] <- Reduce(merge,ss_df)
  } else {
    outputlist[["summarized"]] <- NULL
  }

  if (return_shapefile == T) {
    outputlist[["shapefile"]] <- the_shapefile
  } else {
    outputlist[["shapefile"]] <- NULL
  }

  return(outputlist)
}


#' @title Load Custom Admin Raster
#'
#' @description takes in a shapefile path or an spdf object and converts it into a admin raster using the field variable
#' 
#' @param shapefile_path character string - Full file path to shapefile to be used to make admin raster.
#' @param field character string - Field in shapefile to be used as value in raster.
#' @param simple_raster Simple raster to set extent and resolution of admin raster.
#' @param shapefile option to pass in shapefile rather than path to save some time. Default `NULL`
#'
#' @return Raster with resolution and extent the same as simple_raster, with shapefile field as value
#' @export
load_custom_admin_raster  <- function(shapefile_path, field, simple_raster, shapefile=NULL){
  
  sr = simple_raster
  
  if(is.null(shapefile)){
    shapes <- rgdal::readOGR(shapefile_path)
  } else {
    shapes <- shapefile
  }
  
  # crop
  cropped_shapes <- crop(shapes, extent(sr), snap="out")
  
  #ensure field used for raster is numeric and rename to raster_id
  names(cropped_shapes)[names(cropped_shapes) == field] <- 'raster_id'
  
  #get the number of unique factors in shapefile field of interest
  factor_count <- length(unique(cropped_shapes$raster_id))
  #convert from factor to numeric
  cropped_shapes$raster_id <- as.numeric(as.character(cropped_shapes$raster_id))
  #get the number of unique numbers in shapefile field of interest
  number_count <- length(unique(cropped_shapes$raster_id))
  
  #ensure that no values were dropped in the conversion from factor to numeric
  if(factor_count != number_count) {
    stop("In load_custom_admin_raster(), conversion of shapefile field from factor to number resulted in a different number of unique values")
  }
  # rasterize shapefile with our custom function - make sure field is numeric, not factor
  initial_raster <- rasterize_check_coverage(cropped_shapes, sr, field = "raster_id")
  
  masked_shapes <- raster::mask(x=initial_raster, mask=sr)
  
  return(masked_shapes)
}
