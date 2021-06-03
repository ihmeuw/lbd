#' Directory where central computation keeps R environments for the shared functions
CC_ENV_DIR <- "<<<< FILEPATH REDACTED >>>>"
 

#' @title Load gbd data
#' @description Pull GBD estimates from database and return as a data table
#' @param gbd_type "covariate" or "output", depending which database you need to pull from
#' @param gbd_name GBD cov_name_short if "covariate" and GBD cause_id if "output"
#' @param gaul_list list of GAUL or GADM codes you want to pull
#' @param measure_id if "output", which measure you want to pull (defaults to incidence), Default: 6
#' @param age_group_id if "output", which age group you want to pull (defaults to Under 5), Default: 1
#' @param metric_id if "output", which metric you want to pull (defaults to rate), Default: 3
#' @param year_ids numeric vector, years of gbd data to load, Default: c(2000:2017)
#' @param return_by_age_sex if 'yes', returns sex_id and age_group_id columns in output table, Default: 'no'
#' @param collapse_age_sex if T, collapses age/sex/location/year to location/year via weighted mean, Default: FALSE
#' @param shapefile_version LBD standard shapefile version, Default: 'current'
#' @param gbd_round_id Which round of GBD to pull from, 5=GBD2017, 6=GBD2019 Default: 5
#' @param named_location_field string specifying which name of numerical location identifier that should be returned in 'name' column. could be 'GAUL_CODE' or 'location_id', Default: 'GAUL_CODE'
#' @param loc_ids vector of GBD location IDs used to pull GBD estimates
#' @param ... any other fields to be passed to all of the get_* functions from GBD CC
#' @return Returns 3-column data.table where "name" = (usually) GAUL code, "year" = year, and "mean" = value. Consistent with expected input for calc_raking_factors.
load_gbd_data     <- function(gbd_type,
                              gbd_name,
                              gaul_list,
                              measure_id = 6,
                              age_group_id = 1,
                              metric_id = 3,
                              year_ids = c(2000:2017),
                              return_by_age_sex = "no",
                              collapse_age_sex = FALSE,
                              shapefile_version = 'current',
                              gbd_round_id = 5,
                              named_location_field = "GAUL_CODE",
                              loc_ids = use_global_if_missing(loc_ids),
                              ...) {

  # get GAUL to location_id mapping
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  gaul_to_loc_id <- gaul_to_loc_id[, list(location_id = loc_id, GAUL_CODE)]

  if(is.null(loc_ids)){
    loc_ids <- gaul_to_loc_id[GAUL_CODE %in% gaul_list, location_id]
  }
  
  if (gbd_type == "covariate") { # load covariate data

    # get covariate metadata from covariate_name_short
    metadata <- get_covariate_metadata()
    covariate_id <- metadata[covariate_name_short == tolower(gbd_name), covariate_id]
    covariate_by_age <- metadata[covariate_name_short == tolower(gbd_name), by_age]
    covariate_by_sex <- metadata[covariate_name_short == tolower(gbd_name), by_sex]

    # get covariate data
    source(path_join(CC_ENV_DIR, 'get_covariate_estimates.R'))
    if (covariate_by_age) {
        gbd_estimates <- get_covariate_estimates(covariate_id = covariate_id, location_id = loc_ids, year_id = year_ids, age_group_id = age_group_id, gbd_round_id = gbd_round_id, ...)
    } else {
        gbd_estimates <- get_covariate_estimates(covariate_id = covariate_id, location_id = loc_ids, year_id = year_ids, gbd_round_id = gbd_round_id, ...)
    }

    # collapse to all age, both sexes (if specified, and if there are multiple age and/or sex groups)
    if (collapse_age_sex & gbd_estimates[, uniqueN(age_group_name) > 1 | uniqueN(sex_id) > 1]) {

      # get population data
      source(path_join(CC_ENV_DIR, 'get_population.R'))
      gbd_pops <- get_population(age_group_id = gbd_estimates[, unique(age_group_id)],
                                 location_id = gbd_estimates[, unique(location_id)],
                                 year_id = gbd_estimates[, unique(year_id)],
                                 sex_id = gbd_estimates[, unique(sex_id)],
                                 gbd_round_id = gbd_round_id,
                                 ...)

      # population-weight the covariate data
      gbd_estimates <- merge(gbd_estimates, gbd_pops, by=c('location_id', 'sex_id', 'age_group_id', 'year_id'))
      gbd_estimates <- gbd_estimates[, list(mean_value = weighted.mean(mean_value, population, na.rm=T)), by='location_id,year_id']
    }

    # format and return
    gbd_estimates <- merge(gbd_estimates, gaul_to_loc_id, by="location_id")
    setnames(gbd_estimates, c(named_location_field, "year_id", "mean_value"), c("name", "year", "mean"))

    if(return_by_age_sex=='no') gbd_estimates <- gbd_estimates[, list(name, year, mean)]
    if(return_by_age_sex=='yes') gbd_estimates <- gbd_estimates[, list(name, year, mean, sex_id, age_group_id)]

    return(gbd_estimates)
  } else if (gbd_type == "output") { # load cause data

    # get cause metadata
    source(path_join(CC_ENV_DIR, 'get_cause_metadata.R'))
    metadata <- get_cause_metadata(cause_set_id = 2, gbd_round_id = gbd_round_id)
    cause_id <- suppressWarnings(as.numeric(gbd_name))
    if (is.na(cause_id)) cause_id <- metadata[acause == gbd_name, cause_id]

    # get cause data
    source(path_join(CC_ENV_DIR, 'get_outputs.R'))
    gbd_estimates <- get_outputs(topic = "cause",
                                 version = "best",
                                 gbd_round_id = gbd_round_id,
                                 cause_id = cause_id,
                                 measure_id = measure_id,
                                 metric_id = metric_id,
                                 age_group_id = age_group_id,
                                 location_id = loc_ids,
                                 year_id = year_ids,
                                 ...)

    all_data <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
    setnames(all_data, c(named_location_field, "year_id", "val"), c("name", "year", "mean"))
    all_data <- all_data[, list(name, year, mean)]

    return(all_data)

  } else {
    stop(sprintf("Invalid `gbd_type`: choose 'covariate' or 'output' not '%s'", gbd_type))
  }
}


#' @title Get GBD locations
#' @description Loads in the shapefiles to get a complete set of codes for GBD
#' locations. this is particularly useful if you want to get
#' subnational GBD location codes which are a combination of different
#' administrative levels
#'
#' @param reg region you want to pull codes for. e.g. 'africa'
#' @param rake_subnational logical - do you want to pull location codes for subnational raking or just national raking
#' @param shapefile_version string specifying version date of shapefile to use
#'
#' @return returns a dataframe containing, at a minimum, an ADM_CODE column and an loc_id (ihme loc id) column. If rake_subnational=T, it will also contain ADM0*, ADM1* and rak_level (admin level for raking) columns
get_gbd_locs <- function(reg,
                         rake_subnational = T,
                         shapefile_version = raking_shapefile_version){
  require(sf)
  if(rake_subnational == T) {
    connector <-
      st_read(get_admin_shapefile(admin_level = 0, raking = T, version = shapefile_version), quiet = T) %>%
      st_set_geometry(NULL) %>%
      mutate(ADM0_CODE = as.numeric(as.character(ADM0_CODE))) %>%
      mutate(ADM1_CODE = as.numeric(as.character(ADM1_CODE))) %>%
      filter(ADM0_CODE %in% get_adm0_codes(reg, shapefile_version = shapefile_version, subnational_raking = TRUE))

    ## get the lowest raking level and add it onto connector
    if('ad_level' %in% colnames(connector)){
      connector$rak_level = connector$ad_level ## gadm (ad_level) and gaul (rak_level) use different names for this!
    }

    connector <- connector %>%
      dplyr::select(ADM0_CODE, ADM1_CODE, loc_id, rak_level) %>%
      mutate(rak_level = as.character(rak_level)) %>%
      dplyr::rename(location_id = loc_id) %>%
      mutate(location_id = as.numeric(as.character(location_id)))

    ADM_CODE <- connector$ADM0_CODE
    for(rl in unique(connector$rak_level)){
      ADM_CODE[which(rl == connector$rak_level)] <- connector[[paste0('ADM', rl, '_CODE')]][which(rl == connector$rak_level)]
    }

    connector <- cbind(connector, ADM_CODE)

  } else {
    connector <-
    get_location_code_mapping(shapefile_version = shapefile_version)[ADM_CODE %in% get_adm0_codes(reg, shapefile_version = shapefile_version),
                                                                     list(location_id = loc_id, ADM_CODE)]
  }
  return(data.table(connector))
}


#' @title Get GBD estimates
#' @description Loads national or subnational gbd estimates to be used for raking.
#'
#' @param gbd_name = GBD cov_name_short if "covariate" and GBD cause_id if "output"
#' @param region = name of region you want to pull results for. e.g. 'africa'
#' @param measure_id = if "output", which measure you want to pull (defaults to incidence)
#' @param age_group_id = if "output", which age group you want to pull (defaults to Under 5)
#' @param metric_id = if "output", which metric you want to pull (defaults to rate)
#' @param year_ids = numeric vector of years to pull
#' @param shapefile_version string of dated shapefile version to use when generating list of ihme loc ids and ADM code
#' @param rake_subnational Logical. do you want subnational estimates or just national ones
#' @param gbd_round_id numeric gbd round id to pull from
#' @param decomp_step character decomp_step to pull from, only relevant for gbd_round_id 5 and above
#'
#' @return Returns 3-column data.table where "name" = ihme loc id, "year" = year, and
#'      "mean" = value. Consistent with expected input for rake_cell_pred and calculate_raking_factors.
get_gbd_estimates <- function(gbd_type,
                              gbd_name,
                              region,
                              measure_id = 6,
                              age_group_id = 1,
                              metric_id = 3,
                              year_ids = c(2000:2017),
                              shapefile_version = 'current',
                              rake_subnational = TRUE,
                              gbd_round_id = 6,
                              decomp_step = NULL) {

  ## get GAUL to location_id mapping
  gaul_to_loc_id <- get_gbd_locs(
    reg = region,
    rake_subnational = rake_subnational,
    shapefile_version = shapefile_version
  )

  ## If we had subnational raking on, then we additionally pull in the national estimates
  ## because fractional rates raking needs it
  if(rake_subnational) {
    gaul_to_loc_id_nats <- get_gbd_locs(
      reg = region,
      rake_subnational = FALSE,
      shapefile_version = shapefile_version
    )

    gaul_to_loc_id <- rbindlist(list(as.data.table(gaul_to_loc_id)[, list(location_id, ADM_CODE)],
                                     as.data.table(gaul_to_loc_id_nats)), use.names = TRUE)

  } else {
    gaul_to_loc_id <- as.data.table(gaul_to_loc_id)[, list(location_id, ADM_CODE)]
  }

  gaul_to_loc_id <- unique(gaul_to_loc_id)
  loc_ids <- gaul_to_loc_id[, location_id]

  ## get cause metadata
  source(path_join(CC_ENV_DIR, 'get_cause_metadata.R'))
  metadata <- get_cause_metadata(cause_set_id = 2, gbd_round_id = gbd_round_id)
  cause_id <- suppressWarnings(as.numeric(gbd_name))
  if (is.na(cause_id)) cause_id <- metadata[acause == gbd_name, cause_id]

  ## get cause data
  source(path_join(CC_ENV_DIR, 'get_outputs.R'))
  gbd_estimates <- get_outputs(topic = "cause",
                               version = "latest",
                               gbd_round_id = gbd_round_id,
                               cause_id = cause_id,
                               measure_id = measure_id,
                               metric_id = metric_id,
                               age_group_id = age_group_id,
                               location_id = loc_ids,
                               year_id = year_ids,
                               decomp_step = decomp_step)

  all_data <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
  setnames(all_data, c("location_id", "year_id", "val"), c("name", "year", "mean"))
  all_data <- all_data[, list(name, year, mean)]

  return(all_data)

}


#' @title Get Africa population raster
#' @description DEPRECATED. Given a simple raster, read in a hard-coded population raster brick and crop to simple raster
#'
#' @param simple_raster simple raster as made by \code{\link{build_simple_raster_pop}}
#'
#' @return a masked and cropped population raster to input simple raster
get_population_data     <- function(simple_raster){

  # load population raster
  # RB 27SEPT: Note this is a temporary location, and is only Africa so updates will be necessary
  pop<- brick(file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/raw/covariates/new_20160421/pop_stack.tif'))

  # make sure population is cropped and extented to same as simple_raster
  # this is important, otherwise IDX wont work.
  if(!is.null(simple_raster)){
    pop <- mask(crop(pop,simple_raster),simple_raster)
    extent(pop)=extent(simple_raster)
  }
  return(pop)
}


#' @title Load admin raster
#' @description Load in LBD standard admin shapefile, crops it based on simple_raster, then runs it through \code{\link{rasterize_check_coverage}} to create an admin raster.
#'
#' @param admin_level 0, 1, or 2 - administrative level of shapefile to use
#' @param simple_raster simple raster as made by \code{\link{build_simple_raster_pop}}, used for resolution and extent
#' @param disag_fact disaggregation factor for raster. Each raster gets divided into this number^2 e.g. disag_fact of 2 will result in 4x the number of pixels, as the pixel is divided into 2 both rowwise and columnwise. default NULL/
#' @param shapefile_version LBD standard shapefile version
#'
#' @return a raster with values corresponding to the adm codes at the given level 
load_admin_raster  <- function(admin_level, simple_raster, disag_fact=NULL,
                               shapefile_version = 'current'){

  if(!admin_level %in% c(0,1,2)) stop("admin_level must be either 0, 1, or 2")

  # load admin raster
  if(!is.null(disag_fact)){
    sr = disaggregate(simple_raster,fact=disag_fact)
  } else {
    sr = simple_raster
  }

  # UPDATED: master gaul admin shapefiles
  shapes <- shapefile(get_admin_shapefile(admin_level, version = shapefile_version))

  # The variable we rasterize on must be numeric.
  shapes@data[[paste0('ADM', admin_level,'_CODE')]] <- as.numeric(as.character(shapes@data[[paste0('ADM', admin_level,'_CODE')]]))

  # crop
  cropped_shapes <- crop(shapes, extent(sr), snap="out")

  ## Fix rasterize
  initial_raster <- rasterize_check_coverage(cropped_shapes, sr, field = paste0('ADM', admin_level,'_CODE'))
  if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])!=0) {
    rasterized_shape <- merge(rasterize_check_coverage(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),], sr, field = paste0('ADM', admin_level,'_CODE')), initial_raster)
  }
  if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])==0) {
    rasterized_shape <- initial_raster
  }
  #rasterized_shape <- rasterize(cropped_shapes, sr, field=paste0('ADM', admin_level,'_CODE'))
  masked_shapes <- mask(x=rasterized_shape, mask=sr)

  return(masked_shapes)

}


#################################################################################
### CI Range estimate given draws
#################################################################################
#' @title Cirange
#' @description calculate confidence interval range for a numeric vector. Used for summary stats
#' 
#' @param x numeric vector
#' 
#' @family Summary statistics
#' 
#' @seealso
#' \code{\link{upper}} 
#' \code{\link{lower}} 
#' 
#' @return numeric value for confidence interval range
cirange = function(x){
  z=quantile(x,probs=c(.025,.975),na.rm=T)
  return(z[2]-z[1])
}

#' @title Lower confidence interval
#' @description calculate lower confidence interval (0.025). Used for summary stats
#' 
#' @param x numeric vector
#' 
#' @family Summary statistics
#' 
#' @seealso 
#' \code{\link{upper}} 
#' \code{\link{cirange}} 
#' 
#' @return numeric value for lower confidence interval
lower = function(x) quantile(x,probs=.025,na.rm=T)


#' @title Upper confidence interval
#' @description calculate upper confidence interval (0.975). Used for summary stats
#' 
#' @param x numeric vector
#' 
#' @family Summary statistics
#' 
#' @seealso
#' \code{\link{lower}} 
#' \code{\link{cirange}} 
#' 
#' @return numeric value for lower confidence interval
upper = function(x) quantile(x,probs=.975,na.rm=T)


#' @title Quantile coefficient of dispersion
#' @description calculate Quantile coefficient of dispersion
#' 
#' @param x numeric vector
#' @param quantile_low default 5, numeric between 0-100. 
#' @param quantile_high default 95, numeric between 0-100. 
#' 
#' @family Summary statistics
#' 
#' @return numeric value for Quantile coefficient of dispersion
qcd <- function(x, quantile_low = 5, quantile_high = 95) {

  low <- quantile_low/100
  high <- quantile_high/100
  quantiles <- quantile(x, c(low, high), na.rm = T)
  q_low <- as.numeric(quantiles[1])
  q_high <- as.numeric(quantiles[2])

  output <- (q_high - q_low) / (q_high + q_low)

  return(output)
}


#' @title Quantile coefficient of dispersion odds
#' @description calculate quantile coefficient of dispersion odds
#' 
#' @param x numeric vector
#' @param quantile_low default 5, numeric between 0-100. 
#' @param quantile_high default 95, numeric between 0-100. 
#' 
#' @family Summary statistics
#' 
#' @return numeric value for Quantile coefficient of dispersion odds
odds_qcd <- function(x, quantile_low = 5, quantile_high = 95) {

  low <- quantile_low/100
  high <- quantile_high/100
  quantiles <- quantile(x, c(low, high), na.rm = T)
  q_low <- as.numeric(quantiles[1])
  q_high <- as.numeric(quantiles[2])

  q_low <- q_low / (1-q_low)
  q_high <- q_high / (1-q_high)

  output <- (q_high - q_low) / (q_high + q_low)

  return(output)
}


#' @title Interquartile range
#' @description calculate interquartile range
#' 
#' @param x numeric vector
#' @param quantile_low default 25, numeric between 0-100. 
#' @param quantile_high default 75, numeric between 0-100. 
#' 
#' @family Summary statistics
#' 
#' @return numeric value for interquartile range
iqr <- function(x, quantile_low = 25, quantile_high = 75) {

  low <- quantile_low / 100
  high <- quantile_high / 100
  quantiles <- quantile(x, c(low, high), na.rm = T)
  q_low <- as.numeric(quantiles[1])
  q_high <- as.numeric(quantiles[2])

  output <- (q_high - q_low)

  return(output)

}


#' @title Get percentile
#' @description calculate percentile from a numeric vector
#' 
#' @param x numeric vector
#' @param percentile numeric between 0-100, percentile to get
#' 
#' @family Summary statistics
#' 
#' @return numeric value for given percentile
get_percentile <- function(x, percentile) {

  # Simply get and return a percentile
  percentile <- as.numeric(percentile)
  output <- quantile(x, percentile / 100, na.rm = T)
  return(output)

}

#' @title Get probability below
#' @description calculate probability of getting a number from a vector (`x`) below target value (`value`)
#' 
#' @param x numeric vector
#' @param value target value
#' @param equal_to default F. if F use `<`, if T use `<=`
#' 
#' @family Summary statistics
#' 
#' @return numeric value for probability of getting a number from a vector below target value
p_below <- function(x, value, equal_to = F) {

  value <- as.numeric(value)
  if (equal_to == T) output <- sum(x <= value)
  if (equal_to == F) output <- sum(x < value)

  output <- output / length(x)
  return(output)

}

#' @title Get probability above
#' @description calculate probability of getting a number from a vector (`x`) above target value (`value`)
#' 
#' @param x numeric vector
#' @param value target value
#' @param equal_to default F. if F use `<`, if T use `<=`
#' 
#' @family Summary statistics
#' 
#' @return numeric value for probability of getting a number from a vector above target value
p_above <- function(x, value, equal_to = F) {

  value <- as.numeric(value)
  if (equal_to == T) output <- sum(x >= value)
  if (equal_to == F) output <- sum(x > value)

  output <- output / length(x)
  return(output)

}


#' @title Calculate Coffey-Feingold-Bromberg metric
#' @description Calculate the Coffey-Feingold-Bromberg metric for a numeric vector `v` - useful as a normed measure of variability for a set of proportions. This implementation assumes equal weights. 
#' 
#' Reference:
#' Coffey, M. P., Feingold, M., & Bromberg, J. (1988).
#' A normed measures of variability among proportions.
#' Computational Statistics & Data Analysis, 7(2), 127-141.
#' https://doi.org/10.1016/0167-9473(88)90088-6
#' 
#' @param v numeric vector
#' 
#' @family Summary statistics
#' 
#' @return Coffey-Feingold-Bromberg metric
cfb <- function(v) {

  # calculate mean (u) and sample size (n)
  u <- mean(v, na.rm = T)
  if (is.na(u)) return(NA)
  n = length(v)

  # define numerator (h)
  h = var(v)

  # calculate denominator (max h)
  # (assuming equal weights)
  r <- n*u - floor(n*u)
  h_max <- u*(1-u) - r*(1-r)/n

  # return H statistic (sqrt(h/h_max))
  return(sqrt(h / h_max))

}


#' @title Make cell pred summary
#' @description Takes in raked or raw draw-level estimates and makes a stat summary vector or raster brick as a result. Summarization occurs across columns (e.g. a row with values 1 and 3 will be summarized as 2).
#'
#' @param draw_level_cell_pred raked or raw draw-level estimates
#' @param mask default simple_raster, if creating a raster brick (`return_as_raster=T`), masks output using this raster. Number of rows in draw_level_cell_pred must be a multiple of number of non-na pixels in the mask. 
#' @param return_as_raster default T. If T, returns a raster brick split into years determined by taking the number of rows in the `draw_level_cell_pred` and dividing by the number of non-na pixels in the `mask`. If F, returns a vector of summary statistics of length `nrow(draw_level_cell_pred)`
#' @param summary_stat summary statistic to apply. See `Summary statistics` family of functions.
#' @param ... any other arguments to pass to the `summary_stats` functions
#'            (note that currently will be passed to all functions)
#'            
#' @return a stat summary vector or raster brick
make_cell_pred_summary    <- function(draw_level_cell_pred,
                                      mask                 = simple_raster,
                                      return_as_raster     = TRUE,
                                      summary_stat         = 'mean',
                                      ...){

  # make summary
  summ <- apply(draw_level_cell_pred, 1, summary_stat, ...)

  # put it in a raster
  if(return_as_raster){
    yrs = dim(draw_level_cell_pred)[1]/length(cellIdx(mask))
    message(sprintf('Making a RasterBrick with %i layers',yrs))
    summ <- insertRaster(mask,  matrix(summ,  ncol = yrs))
  }


  return(summ)

}


#' @title Make admin prediction summary
#' @description Take an admin pred object and an sp_hierarchy list and generate a sorted,
#' cleaned table for a given set of summary statistics
#'
#' @param admin_pred the admin draws object (e.g. admin_1, etc)
#' @param sp_hierarchy_list spatial hierarchy object from the admin draws RData files
#' @param summary_stat summary statistic to apply
#' @param ... any other arguments to pass to the `summary_stats` functions
#'            (note that currently will be passed to all functions)
#' @return data table of summary stats and admin info
#' @examples
#' # load("indicator_raked_admin_draws_eb_bin0.RData") # Replace with your admin draws
#' make_admin_pred_summary(admin_2,
#'                         sp_hierarchy_list,
#'                         c("mean", "cirange", "upper", "lower"))
make_admin_pred_summary    <- function(admin_pred,
                                       sp_hierarchy_list,
                                       summary_stats = 'mean',
                                       ...){

  ### Get set up
  str_match <- stringr::str_match

  # Split up your input data
  ad_code <- subset(admin_pred, select = grep("ADM[0-9]_CODE", names(admin_pred)))
  year <- subset(admin_pred, select = "year")
  pop <- subset(admin_pred, select = "pop")
  draws <- subset(admin_pred, select = grep("V[0-9]*", names(admin_pred)))

  # Get the admin level
  ad_code_name <- names(ad_code)
  ad_level <- as.numeric(str_match(ad_code_name,"ADM([0-9])_CODE")[,2])

  # Get all admin levels
  all_admin_levels <- as.numeric(str_match(names(sp_hierarchy_list), "ADM([0-9])_CODE")[,2])
  all_admin_levels <- unique(all_admin_levels)[!is.na(unique(all_admin_levels))]

  ### Make summary
  summ_list <- lapply(summary_stats, function(ss) {
    apply(draws, 1, ss, ...)})

  if (length(summ_list) > 1) {
    output_df <- as.data.table(cbind(year, ad_code, do.call(cbind, summ_list)))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  } else if (length(summ_list) == 1) {
    output_df <- as.data.table(cbind(year, ad_code, summ_list[[1]]))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  }

  ### Add on identifiers

  # Drop smaller admin levels
  drop_levels <- all_admin_levels[all_admin_levels > ad_level]

  if (length(drop_levels) > 0) {
    for (d in drop_levels) {
      drop_cols <- names(sp_hierarchy_list)[grepl(paste0("ADM",d), names(sp_hierarchy_list))]
      sp_hierarchy_list <- subset(sp_hierarchy_list, select = !(names(sp_hierarchy_list) %in% drop_cols))
    }
  }
  sp_hierarchy_list <- unique(sp_hierarchy_list)

  output_df <- merge(sp_hierarchy_list, output_df, all.y = T, all.x = F)

  # Clean up & sort
  order_vars <- c(paste0("ADM", 0:ad_level, "_NAME"), "year")
  setorderv(output_df, order_vars)

  # Get col order
  ad_col_order <- sort(names(output_df)[grep("AD*_*", names(output_df))])
  cols_to_order <- c(ad_col_order, "region", "year", summary_stats)
  cols_to_order <- cols_to_order[cols_to_order %in% names(output_df)]
  if (length(cols_to_order) < length(names(output_df))) {
    other_cols <- names(output_df)[!(names(output_df) %in% cols_to_order)]
    cols_to_order <- c(cols_to_order, other_cols)
  }

  setcolorder(output_df, cols_to_order)

  return(output_df)

}

#' @title summarize_admins
#' @description Function to summarize admin_pred objects
#'
#' This is a wrapper for `make_admin_pred_summary()`
#'
#' @param ind indicator
#' @param ig indicator_group
#' @param summstats Summary statistics (functions) to compute.
#'                  Order will be the order they are written in csv
#'                  This is passed to `make_admin_pred_summary()`
#' @param raked Raked (T), unraked (F), or both (`c(T,F)`)?
#' @param ad_levels Admin levels to summarize (0, 1, and 2 by default)
#' @return Writes csv files to `sharedir/pred_derivatives/admin_summaries/`
#' @examples
#' summarize_admins(summstats = c("mean", "lower", "upper", "cirange"),
#'                  ad_levels = c(0,1,2),
#'                  raked     = c(T,F))
#' @rdname summarize_admins
summarize_admins <- function(ind = indicator,
                             ig = indicator_group,
                             summstats = c("mean", "lower", "upper", "cirange"),
                             raked = c(T,F),
                             ad_levels = c(0,1,2),
                             file_addin = NULL,
                             ...) {

  sharedir       <- file.path(fp_list['mbg_root'], ig, ind)
  input_dir <- paste0(sharedir, "/output/", run_date, "/")
  output_dir <- paste0(input_dir, "/pred_derivatives/admin_summaries/")
  dir.create(output_dir, recursive = T, showWarnings = F)

  # Convert raked to character
  rr <- character()
  if (T %in% raked) rr <- c(rr, "raked")
  if (F %in% raked) rr <- c(rr, "unraked")

  # If file_addin present, use it
  if (!is.null(file_addin)) file_addin <- paste0("_", file_addin)
  if (is.null(file_addin)) file_addin <- ""

  # Summarize and save admin preds
  for (rake in rr) {
    load(paste0(input_dir, ind, "_", rake, "_admin_draws_eb_bin0_0.RData"))
    sp_hierarchy_list <- mutate_if(sp_hierarchy_list, is.factor, as.character)
    sp_hierarchy_list <- mutate_at(sp_hierarchy_list, grep('_CODE', names(sp_hierarchy_list), value = T), as.numeric)

    for (ad in ad_levels) {
      message(paste0("Summarizing ", ind, ": admin ", ad, " (", rake, ")"))
      ad_summary_table <- make_admin_pred_summary(admin_pred = get(paste0("admin_", ad)),
                                                  sp_hierarchy_list,
                                                  summary_stats = summstats,
                                                  ...)
      fwrite(ad_summary_table,
             file = paste0(output_dir, ind, "_admin_", ad, "_", rake, file_addin, "_summary.csv"))
    }
  }
}


#' @title Save Postestimation objects
#' @description Takes in objects created during postestimation and saves them in a standard format in accordance with their object type.
#'
#' @param x any R object
#' @param filetype either 'rdata', 'raster', or 'csv'. The type of file to save `x` as. 
#' @param filename The name of the file to be saved. Is prefixed by `indic_`
#' @param indic indicator name, default indicator.
#'            
#' @return None
save_post_est   <- function(x,
                            filetype,
                            filename,
                            indic = indicator){

  output_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indic, '/output/', run_date)


  dir.create(output_dir, showWarnings = FALSE)

  filetype = tolower(filetype)
  if(!filetype %in% c('rdata','raster','csv'))
    stop('filetype argument has to be either rdata or raster or csv')


  if(filetype=='raster')
    writeRaster(
      x,
      file = paste0(output_dir, '/', indic,'_',filename),
      format='GTiff',
      overwrite = TRUE
    )

  if(filetype=='rdata')
    save(
      x,
      file = paste0(output_dir, '/', indic,'_',filename,'.RData'),
      compress = TRUE
    )

  if(filetype=='csv')
    write.csv(
      x,
      file = paste0(output_dir, '/', indic,'_',filename,'.csv')
    )


}

#' @title Load Cell Pred
#' @description Loads in RData file containing cell pred object from standard location. 
#' 
#' filenames:
#' u5m = T, ageasindic = T : indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'NA.RData'
#' u5m = T, ageasindic = F : indicator,'_age',agebin,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'.RData
#' u5m = F : indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0.RData'
#'
#' @param indicator_group indicator_group of model run
#' @param indicator indicator of model run
#' @param rd run_date of model run
#' @param region region of model run
#' @param agebin agebin of model run, usually 0 if not using ages
#' @param u5m Boolean default F, is under-5 mortality team? Appends `other` to end of filename if T
#' @param other String, default ''. Appended to end of filename.
#' @param ageasindic Boolean default T, changes the filename to load in if `u5m` is T
#'            
#' @return Cell pred object
load_cell_preds <- function(indicator_group,
                            indicator,
                            rd = run_date,
                            region,
                            agebin,
                            u5m=FALSE,
                            other='',
                            ageasindic=TRUE){


  if(u5m){
    if(ageasindic==FALSE){
      load(paste0(fp_list['mbg_root'],indicator_group,'/',indicator,'/output/',rd,'/',
                  indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'NA.RData')) # the 0 are no holdout
    } else {
      load(paste0(fp_list['mbg_root'],indicator_group,'/',indicator,'_age',agebin,'/output/',rd,'/',
                  indicator,'_age',agebin,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'.RData'))
    }
    cell_pred <- cptmp
  } else {
    load(paste0(fp_list['mbg_root'],indicator_group,'/',indicator,'/output/',rd,'/',
                indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0.RData')) # the 0 are no holdouts
  }
  return(cell_pred)
}


#' @title Load Cell Pred Stack
#' @description Loads in RData file containing stacked results cell draws in model run date folder from filename: indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0_stacked_results.RData'. Not implemented for holdouts
#'
#' @param indicator_group indicator_group of model run
#' @param indicator indicator of model run
#' @param rd run_date of model run, default `run_date` object
#' @param region region of model run
#' @param agebin agebin of model run, usually 0 if not using ages
#'            
#' @return Cell pred object
load_cell_preds_stack <- function(indicator_group,
                                  indicator,
                                  rd = run_date,
                                  region,
                                  agebin){

  load(paste0(fp_list['mbg_root'],indicator_group,'/',indicator,'/output/',rd,'/',
              indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0_stacked_results.RData')) # the 0 are no holdouts

  return(cell_pred)
}


#' @title Generate Model Fit Statistics
#' @description Takes in in and out of sample data, loads in cell draws, then finds the error, MAE, RMSE, and coverage for in and out of sample predictions. Generates a new subfolder in the model run folder `/fit_stats` and saves a csv there
#'
#' @param is_data in sample data included in the model run
#' @param oos_data out of sample data left out of the model run for validation
#' @param indicator indicator of model run
#' @param indicator_group indicator_group of model run
#' @param run_date run_date of model run
#' @param pathaddin typically `paste0('_bin',age,'_',reg,'_',holdout)`, used at the end of model output files
#'            
#' @return None
fit_stats <- function(is_data,
                      oos_data,
                      indicator,
                      indicator_group,
                      run_date,
                      pathaddin) {

  ### Fit statistics
  model_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/', run_date)
  image_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/model_image_history')

  ## Append in- and out-of-sample coords
  is_data <- is_data[, c('latitude','longitude','period','N',indicator), with=FALSE]
  is_data <- is_data[, oos := 0]
  oos_data <- oos_data[, c('latitude','longitude','period', 'N',indicator), with=FALSE]
  oos_data <- oos_data[, oos := 1]
  all_data <- rbind(is_data, oos_data)
  if(indicator_family=='binomial') all_data <- all_data[, obs := get(indicator) / N]
  if(indicator_family!='binomial') all_data <- all_data[, obs := get(indicator)]

  ## Extract predictions at input data and bind to input datatable
  message("Loading cell_preds and inputs...")
  load(paste0(image_dir,'/', run_date, pathaddin, '.RData'))
  mean_preds <- brick(paste0(model_dir, '/', indicator, '_prediction_eb', pathaddin))
  load(paste0(model_dir, '/', indicator, '_cell_draws_eb', pathaddin, '.RData'))

  cell_pred.dt <- as.data.table(cell_pred)
  cols <- names(cell_pred.dt)
  message("Making upper and lower credible interval rasters...")
  cell_pred.dt <- cell_pred.dt[ , upper := apply(.SD, 1, quantile, c(.975), na.rm=TRUE), .SDcols=cols]
  cell_pred.dt <- cell_pred.dt[ , lower := apply(.SD, 1, quantile, c(.025), na.rm=TRUE), .SDcols=cols]
  upper <- cell_pred.dt[, upper]
  lower <- cell_pred.dt[, lower]

  upper_raster <- insertRaster(simple_raster,
                               matrix(upper,
                                      ncol = length(unique(all_data$period))))
  lower_raster <- insertRaster(simple_raster,
                               matrix(lower,
                                      ncol = length(unique(all_data$period))))

  message("Extracting upper and lower values at all coordinates...")
  # extract cluster covariates
  all_data$longitude<-as.numeric(all_data$longitude)
  all_data$latitude <-as.numeric(all_data$latitude )

  # add dummy column for time-varying covariates
  all_data$pred <- NA
  all_data$upper <- NA
  all_data$lower <- NA

  # loop through time varying predictions to insert them
  for (period in sort(unique(all_data$period))) {

    # find matching rows
    message(paste0('Extracting predictions for period ',period))
    idx_tv<- which(all_data$period == period)

    # get period prediction raster
    craster  <- mean_preds[[period]]
    crs(craster)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

    uraster  <- upper_raster[[period]]
    crs(uraster)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

    lraster  <- lower_raster[[period]]
    crs(lraster)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

    # extract values
    extr <- extract(craster, all_data[idx_tv, c('longitude', 'latitude'),with=F])
    upper_extr <- extract(uraster, all_data[idx_tv, c('longitude', 'latitude'),with=F])
    lower_extr <- extract(lraster, all_data[idx_tv, c('longitude', 'latitude'),with=F])

    # add into df
    all_data[['pred']][idx_tv]=extr
    all_data[['upper']][idx_tv]=upper_extr
    all_data[['lower']][idx_tv]=lower_extr

  }

  ## Make fit statistics
  ## Add to df to combine later however we want
  message("Calculating all fit statistics at each coordinate...")

  # Drop missing predictions
  message(paste0('Predictions missing for ', length(all_data[is.na(pred), pred]), '/', length(all_data[ , pred]), ' observations'))
  message('If you have missings, this either due to missing covariate values (bad) or points that were in the buffer zone (fine).')
  message('Dropping missing rows...')
  df_no_nas <- all_data[!is.na(pred), ]

  # Coverage
  df_no_nas <- df_no_nas[obs >= lower & obs <= upper, covered := 1]
  df_no_nas <- df_no_nas[obs < lower | obs > upper, covered := 0]

  # Error
  df_no_nas <- df_no_nas[, error := obs - pred]

  # Make statistics over entire IS dataset
  mean_error <- mean(df_no_nas[oos==0, error])
  mean_absolute_error <- mean(df_no_nas[oos==0, abs(error)])
  rmse <- sqrt(mean(df_no_nas[oos==0, error]^2))
  coverage <- mean(df_no_nas[oos==0, covered])

  # Make statistics over entire OOS dataset
  oos_mean_error <- mean(df_no_nas[oos==1, error])
  oos_mean_absolute_error <- mean(df_no_nas[oos==1, abs(error)])
  oos_rmse <- sqrt(mean(df_no_nas[oos==1, error]^2))
  oos_coverage <- mean(df_no_nas[oos==1, covered])

  # Average over all observations
  message("IN-SAMPLE:")
  message(paste0('                       Average error:    ', round(mean_error, 3)))
  message(paste0('                       Average MAE:      ', round(mean_absolute_error, 3)))
  message(paste0('                       Average RMSE:     ', round(rmse, 3)))
  message(paste0('                       Average coverage: ', round(coverage, 3)))

  message("OUT-OF-SAMPLE:")
  message(paste0('                       Average error:    ', round(oos_mean_error, 3)))
  message(paste0('                       Average MAE:      ', round(oos_mean_absolute_error, 3)))
  message(paste0('                       Average RMSE:     ', round(oos_rmse, 3)))
  message(paste0('                       Average coverage: ', round(oos_coverage, 3)))

  fit_folder <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/', run_date, '/fit_stats')
  message(paste0('Saving fit statistics in ', fit_folder))
  dir.create(fit_folder, showWarnings = FALSE)
  message(pathaddin)
  write.csv(df_no_nas, file = paste0(fit_folder, '/fit_stats_', pathaddin, '.csv'))

}


#' @title Generate Model Fit Statistics for draws
#' @description Takes model draws and calculates fit statistics, adding them on to the draws as additional columns. Adds on mean, upper, lower, covered, and error columns. Prints out statistics over the full OOS dataset: error, MAE, RMSE, coverage
#'
#' @param draws cell_pred object in data.table form with draws names having a prefix
#' @param draw_prefix string, the prefix in the name for draw columns used as an identifier, e.g. for draws V1-V100, this would be "V".
#' @param observed column named for observed value to be compared against for calculating coverage and error?
#'            
#' @return cell_pred object with new columns for fit statistics
fit_stats_ho_id <- function(draws,
                            draw_prefix,
                            observed) {

  message('Your draws must contain columns oos and ho_id.')
  message('Summarizing draws called phat_*...')
  cols <- grep(draw_prefix, names(draws), value=TRUE)
  draws <- draws[, mean := .(mean = rowMeans(.SD)), by = ho_id, .SDcols=cols]
  draws <- draws[ , upper := apply(.SD, 1, quantile, c(.975), na.rm=TRUE), .SDcols=cols]
  draws <- draws[ , lower := apply(.SD, 1, quantile, c(.025), na.rm=TRUE), .SDcols=cols]

  message('Calculating fit statistics at each ho_id...')
  # Coverage
  draws <- draws[get(observed) >= lower & get(observed) <= upper, covered := 1]
  draws <- draws[get(observed) < lower | get(observed) > upper, covered := 0]

  # Error
  draws <- draws[, error := get(observed) - mean]

  # Make statistics over entire OOS dataset
  oos_mean_error <- mean(draws[oos==1, error])
  oos_mean_absolute_error <- mean(draws[oos==1, abs(error)])
  oos_rmse <- sqrt(mean(draws[oos==1, error]^2))
  oos_coverage <- mean(draws[oos==1, covered])

  message("OUT-OF-SAMPLE:")
  message(paste0('                       Average error:    ', round(oos_mean_error, 3)))
  message(paste0('                       Average MAE:      ', round(oos_mean_absolute_error, 3)))
  message(paste0('                       Average RMSE:     ', round(oos_rmse, 3)))
  message(paste0('                       Average coverage: ', round(oos_coverage, 3)))

  message('Returning draws with fit stat columns added.')
  return(draws)

}


#' @title Weighted Vector Aggregation
#' @description Weighted aggregation of a vector or matrix
#'
#' @param p_i vector or matrix to aggregate
#' @param N_i sample size or weights to aggregate with
#'            
#' @return aggregated value
Agg <- function(p_i, N_i){
  return(sum(p_i * N_i)/sum(N_i))
}


#' @title Save Model Run as MBG Covariate
#' @description DEPRECATED? Created objects do not follow new dated folder structure for covariates. Reads in the raked and unraked mean rasters and saves out each year layer into the covariate directory in a new folder `{indicator}_un/raked`.
#' 
#' @note TODO update to use new folder structure in covariates folder
#'
#' @param indicator_group indicator_group of model run
#' @param indicator indicator of model run
#' @param run_date run_date of model run
#' @param measure measure to use in folder name (should be `"mean"` since the mean raster is read in)
#'            
#' @return None
save_best_model <- function(indicator_group, indicator, run_date, measure) {

  ######################################################################
  ######################################################################

  ## Set repo location and indicator group
  core_repo <- file.path(fp_list$code_root, 'lbd_core/')


  ## Load libraries and miscellaneous MBG project functions.
  mbg_functions <- c('mbg_functions.R', 'prep_functions.R',
                     'covariate_functions.R', 'misc_functions.R',
                     'post_estimation_functions.R', 'gbd_functions.R',
                     'shiny_functions.R', 'holdout_functions.R',
                     'categorical_variable_functions.R',
                     'validation_functions.R',
                     'seegMBG_transform_functions.R')
  source(paste0(core_repo, '/mbg_central/setup.R'))
  source_functions(paste(core_repo, 'mbg_central', mbg_functions, sep = '/'))
  load_R_packages(c('foreign', 'rgeos', 'data.table','raster','rgdal','INLA',
                    'seegSDM','seegMBG','plyr','dplyr'))

  ######################################################################
  ######################################################################

  ## load recentmodel runs
  raked <- brick(paste0(fp_list['mbg_root'], indicator_group, "/", indicator, "/output/", run_date, "/", indicator, "_mean_raked_raster.tif"))
  unraked <- brick(paste0(fp_list['mbg_root'], indicator_group, "/", indicator, "/output/", run_date, "/", indicator, "_mean_raster.tif"))

  ## extend my africa rasters
  raked <- extend(raked, extent(-180, 180, -90, 90), keepres=TRUE)
  unraked <- extend(unraked, extent(-180, 180, -90, 90), keepres=TRUE)

  ## load in a covariate to match my new global rasters to the same pixels
  ## Hard-code central directories
  central_cov_dir <- paste0(fp_list$covariate_root, '09_MBG_covariates/') # central folder per Lucas
  central_tv_covs <- c('evi','lights_new','LST_day','total_pop','rates','malaria','fertility','urban_rural', 'land_cover', 'LST_avg', 'gpcp_precip', 'aridity_cruts')
  central_ntv_covs <- c('access','irrigation','LF','LF_vector','reservoirs','aridity','elevation','annual_precip','PET','dist_rivers_lakes','dist_rivers_only','lat','lon','latlon')
  evi <- brick(paste0(central_cov_dir, 'EVI_stack.tif'))

  ## make sure we match other covariates
  raked <- setExtent(raked, evi, keepres = TRUE, snap = TRUE)
  unraked <- setExtent(unraked, evi, keepres = TRUE, snap = TRUE)

  ## Make all dirs
  main_dir <- file.path(fp_list$covariate_root, '00_MBG_STANDARD/')
  dir.create(paste0(main_dir, indicator, '_raked/'))
  dir.create(paste0(main_dir, indicator, '_raked/', measure))
  dir.create(paste0(main_dir, indicator, '_raked/', measure, '/1y/'))
  dir.create(paste0(main_dir, indicator, '_unraked/'))
  dir.create(paste0(main_dir, indicator, '_unraked/', measure))
  dir.create(paste0(main_dir, indicator, '_unraked/', measure, '/1y/'))

  ## Write raster layer for each year, raked and unraked
  for(i in 1:length(names(unraked))) {
    subset_raked <- raked[[i]]
    subset_unraked <- unraked[[i]]
    year <- (i-1) + 2000
    message(paste0("Saving ", year))
    writeRaster(subset_raked,
                filename = paste0(fp_list$covariate_root, '00_MBG_STANDARD/', indicator, '_raked/', measure, '/1y/', indicator, '_raked_', measure, '_1y_', year, "_00_00.tif"),
                format = "GTiff",
                overwrite = TRUE)
    writeRaster(subset_unraked,
                filename = paste0(fp_list$covariate_root, '00_MBG_STANDARD/', indicator, '_unraked/', measure, '/1y/', indicator, '_unraked_', measure, '_1y_', year, "_00_00.tif"),
                format = "GTiff",
                overwrite = TRUE)
  }

  message(sprintf('All layers successfully saved in %s', main_dir))
  message(paste0('Indicator ', indicator, ' saved as ', indicator, '_raked and ', indicator, '_unraked to call in MBG config files in the fixed_effects parameter.'))

}


#' @title Calculate (very rough) covariate importance/weights in a stacking INLA model
#' @description This function attempts to provide some quantitative assessment of how important (in a predictive sense) each raw covariate (i) works out to be when running a stacked generalization ensemble through INLA. For each of the child models (j), it produces a covariate weight importance specific to that child stacker (w_ij), and then it combines the raw covariate weights across stackers using the scaled child model coefficient fit (Bscale_j) in INLA: wt_raw_cov_i = sum_j (w_ij * Bscale_j)
#'
#' @param ind indicator name
#' @param ind_gp indicator group name
#' @param reg region name
#' @param age agebin of model run, usually 0 if not using ages
#' @param holdout holdout of model run. 0 if running an in-sample model without holdouts
#'            
#' @return matrix of raw covariate "importance" weights. Each row contains a model, and each column contains the importance weights for a raw covariate. The final row contains the raw covariates weight in the final INLA model, and the final column contains the (scaled) weights of each child model in the final INLA model fit.
get.cov.wts <- function(rd,
                        ind,
                        ind_gp,
                        reg,
                        age = 0,
                        holdout = 0) {
  ## ##########################################
  ## load the workspaces and objects we need ##
  ## ##########################################

  pathaddin <- paste0('_bin', age, '_', reg, '_', holdout)

  this_config <- fread(paste0(fp_list['mbg_root'], ind_gp, '/', ind, '/output/', rd, '/config.csv'))
  stacker_list <- this_config[V1 == 'stacked_fixed_effects', V2]
  stackers_used <- strsplit(stacker_list," ")
  stackers_used <- stackers_used[[1]][stackers_used[[1]] != "+"]

  ## load the data used to fit the stackers and reconstruct the design matrix
  load(paste0(fp_list['mbg_root'], ind_gp, '/', ind, '/model_image_history/', rd, pathaddin, '.RData'))
  fit.data <- as.data.frame(df)
  fit.data <- fit.data[, -(1:max(match(stackers_used, colnames(fit.data))))] # this drops all the ID and post-stacking variables, leaving just the pre-stacking ones
  if("gaul_code" %in% colnames(fit.data)) fit.data$gaul_code <- NULL
  if("period" %in% colnames(fit.data)) fit.data$period <- NULL
  X <- fit.data
  rc <- colnames(X)

  ## load the stacker fits
  load(paste0(fp_list['mbg_root'], ind_gp, '/', ind, '/output/', rd,
              sprintf('/child_model_list_%s_%i.RData', reg, holdout)))
  smo <- child_models
  rm(child_models)

  ## load the INLA/TMB fit
  load(paste0(fp_list['mbg_root'], ind_gp, '/', ind, '/output/', rd,
              sprintf('/%s_model_eb_bin%i_%s_%i.RData', ind, age, reg, holdout)))
  fit <- res_fit
  rm(res_fit)

  ## make a matrix to hold the p-values/importances
  imp.mat <- as.data.frame(matrix(ncol = length(rc), nrow = length(smo), dimnames = list(names(smo), rc)))

  ## #####################################################################################
  ## now, for each of the models, add the pvalues to the corresponding rows and columns ##
  ## #####################################################################################

  ## ~~~~~
  ## gam ~
  ## ~~~~~
  gam <- smo[['gam']]
  if(!is.null(gam) & class(gam)[1] != 'try-error'){
    smoothed   <- summary(gam)$s.table[, 4] ## smoothed table
    names(smoothed) <- substr(names(smoothed), 3, nchar(names(smoothed))-1)
    unsmoothed <- summary(gam)$p.table[, 4] ## parametric table
    all.p <- c(smoothed, unsmoothed)

    ## now match names and stick into our pval.matrix
    imp.mat["gam",] <- all.p[rc]

    ##  convert from p-val to `importance`
    imp.mat["gam", ] <- -log(imp.mat["gam", ])
    imp.mat["gam", ][is.na(imp.mat["gam", ])] <- 0
    if (sum(is.infinite(unlist(imp.mat["gam",]))) > 0) { # if a covariate has 'infinite' importance, assign it 1 and everything else 0
      imp.mat["gam",][!is.infinite(unlist(imp.mat["gam",]))] <- 0
      imp.mat["gam",][is.infinite(unlist(imp.mat["gam",]))] <- 1
    }
    imp.mat["gam", ] <- imp.mat["gam", ] / sum(imp.mat["gam", ])

  }

  ## ~~~~~
  ## gbm ~
  ## ~~~~~
  gbm <- smo[['gbm']]
  if(!is.null(gbm) & class(gbm)[1] != 'try-error'){
    rel.inf <- gbm$contributions
    imp.mat["gbm",] <- rel.inf[rc, "rel.inf"]

    ## convert to scaled importance
    imp.mat["gbm", ] <- imp.mat["gbm", ] / sum(imp.mat["gbm", ])

  }

  ## all the penalized regression DO NOT give SD or p-vals...
  ## I 'standardize' the coefs as per a suggestion in this thread:
  ## https://stats.stackexchange.com/questions/14853/variable-importance-from-glmnet

  ## ~~~~~~~
  ## lasso ~
  ## ~~~~~~~
  lasso <- smo[['lasso']]
  if(!is.null(lasso) & class(lasso)[1] != 'try-error'){
    l <- lasso$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(lasso$lambda == l)

    sds <- apply(X, 2, sd, na.rm = T)
    unscaled.coefs <- lasso$beta[, l.idx]
    unscaled.coefs <- unscaled.coefs[rc]

    ## scaled
    scaled.coefs   <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    imp.mat["lasso", ] <- scaled.coefs[rc]

    ## convert to scaled importance
    imp.mat["lasso", ] <- imp.mat["lasso", ] / sum(imp.mat["lasso", ])

  }

  ## ~~~~~~~
  ## ridge ~
  ## ~~~~~~~
  ridge <- smo[['ridge']]
  if(!is.null(ridge) & class(ridge)[1] != 'try-error'){
    l <- ridge$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(ridge$lambda == l)

    sds <- apply(X, 2, sd, na.rm = T)
    unscaled.coefs <- ridge$beta[, l.idx]
    unscaled.coefs <- unscaled.coefs[rc]

    ## scaled
    scaled.coefs   <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    imp.mat["ridge", ] <- scaled.coefs[rc]

    ## convert to scaled importance
    imp.mat["ridge", ] <- imp.mat["ridge", ] / sum(imp.mat["ridge", ])

  }

  ## ~~~~~~
  ## enet ~
  ## ~~~~~~
  enet <- smo[['enet']]
  if(!is.null(enet) & class(enet)[1] != 'try-error'){
    l <- enet$cv_1se_lambda ## this is the CV lambda from the child fit
    l.idx <- which(enet$lambda == l)

    sds <- apply(X, 2, sd, na.rm = T)
    unscaled.coefs <- enet$beta[, l.idx]
    unscaled.coefs <- unscaled.coefs[rc]

    ## scaled
    scaled.coefs   <- abs(unscaled.coefs) * sds

    ## put them in the matrix
    imp.mat["enet", ] <- scaled.coefs[rc]

    ## convert to scaled importance
    imp.mat["enet", ] <- imp.mat["enet", ] / sum(imp.mat["enet", ])

  }

  ## ~~~~~~~~~
  ## xgboost ~
  ## ~~~~~~~~~
  xgboost <- smo[['xgboost']]
  if(!is.null(xgboost) & class(xgboost)[1] != 'try-error'){
    load_R_packages("caret")
    # Extract row names to then assign to column names
    xg_names <- rownames(varImp(xgboost, scale = FALSE)$importance)


    # Extract coefficients, already correctly scaled
    scaled.coefs <- varImp(xgboost, scale = FALSE)$importance[[1]]


    # Assign column names
    names(scaled.coefs) <- xg_names

    # put them in the matrix
    imp.mat["xgboost", ] <- scaled.coefs[rc]
  }

  ## ##########################################
  ## Now we propagate through the INLA coefs ##
  ## ##########################################

  ## to account for different scaling in INLA we also create SDs of the covariates that go into INLA
  inla.X <- df[, paste0(names(smo), "_cv_pred"), with=F]
  if (this_config[V1 == "indicator_family", V2] == "binomial" & this_config[V1 == "stackers_in_transform_space", V2] == T) {
    inla.X <- logit(inla.X)
  }
  inla.sds   <- apply(inla.X, 2, sd, na.rm = T)

  ## get the coefficients
  if ("sdrep" %in% names(fit)) { # TMB, w/ stackers as fixed effects
    inla.coefs <- fit$sdrep$par.fixed[names(fit$sdrep$par.fixed) == "alpha_j"]
    names(inla.coefs) <- fit$fenames
    inla.coefs <- inla.coefs[stackers_used]
  } else if ("covar" %in% names(fit$summary.random)) { # INLA, w/ stackers as random effects (to sum to one)
    inla.coefs <- fit$summary.random$covar$"mean"
    names(inla.coefs) <- stackers_used
  } else { # INLA, w/ stackers as fixed effects
    inla.coefs <- fit$summary.fixed[stackers_used, "mean"]
    names(inla.coefs) <- stackers_used
  }

  ## get the scaled coefficients
  inla.coefs <- inla.coefs[names(smo)]
  scaled.inla.coefs <- inla.coefs * inla.sds

  ## make and scale the INLA weighted relative importance and add as a row
  inla.rel.imp <- apply(imp.mat, 2, function(x) sum(x*scaled.inla.coefs))
  inla.rel.imp <- abs(inla.rel.imp)
  inla.rel.imp <- inla.rel.imp/sum(inla.rel.imp)
  imp.mat <- rbind(imp.mat, "INLA COMBINED" = inla.rel.imp)

  ## add the scaled coefs as a column
  imp.mat <- cbind(imp.mat, "SCALED.INLA.COEFS" = c(scaled.inla.coefs, NA))

  ## save to general output directory
  save(imp.mat,
       file = paste0(fp_list['mbg_root'], ind_gp, '/', ind, '/output/', rd,
                     sprintf('/cov_wts_%s_holdout_%i.RData', reg, holdout)))

  return(imp.mat)
}


#' @title Plot Covariate Weights
#' @description makes a 2d heatmap/color plot of the matrix produced by \code{\link{get.cov.wts}}
#'
#' @param rd run date
#' @param ind indicator name
#' @param ind_gp indicator group name
#' @param reg region name
#' @param plot.inla.col boolean default TRUE, should scaled inla coefficients for the child models should be included in the plot?
#' @param age agebin of model run, usually 0 if not using ages
#' @param holdout holdout of model run. 0 if running an in-sample model without holdouts
#'            
#' @return ggplot object of heatmap
plot.cov.wts <- function(rd,
                         ind,
                         ind_gp,
                         reg,
                         plot.inla.col = TRUE,
                         age = 0,
                         holdout = 0) {

  ## load a cov.wts object
  load(paste0(fp_list['mbg_root'], ind_gp, '/', ind, '/output/', rd,
              sprintf('/cov_wts_%s_holdout_%i.RData', reg, holdout)))
  cov.wts <- imp.mat

  ## remove the final column which has inla weight
  if(!plot.inla.col){
    cov.wts <- cov.wts[, -ncol(cov.wts)]
  }else{
    ## rescale to be -1:1
    cov.wts[, ncol(cov.wts)] <- cov.wts[, ncol(cov.wts)] / sum(cov.wts[, ncol(cov.wts)], na.rm = TRUE)
  }

  ## melt
  cw.m <- as.data.table(reshape2::melt(as.matrix(cov.wts)))
  colnames(cw.m)[1:3] <- c('Model', "Covar", "Imp")
  ## cw.m <- na.omit(cw.m)

  ## reorder factors
  cw.m$Model <- factor(cw.m$Model, levels = c("INLA COMBINED", sort(setdiff(unique(cw.m$Model), "INLA COMBINED"))))

  ## setup the plot

  ## pdf("~/test2.pdf", width = 10, height = 3)
  base_size <- 9
  p <- ggplot(cw.m, aes(Covar, Model)) +
    geom_tile(aes(fill = Imp), colour = "white") +
    scale_fill_gradient2(low = "red", mid = 'white', high = "steelblue", midpoint = 0) +
    theme_grey(base_size = base_size) +
    labs(x = "Covariate", y = "Model") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme(axis.text.x = element_text(size = base_size * 0.8,
                                     angle = 270,
                                     hjust = 0,
                                     colour = "grey50"))
  ## print(p)
  ## dev.off()
  return(p)

}


#' @title Make AROC raster
#' @description provides weighted aroc taken across years of prediction. Uses the input raster to calculate the aroc for each pixel and returns it as a single layer raster.
#'
#' @param pred.rast takes a stacked raster object (raster brick), different layers are different years of prediction by default it assumes that spacing between layers is 1 year and that the top layer is the earliest year layers
#' @param year.map Default NULL, if the above default is not correct, this vector (of length == raster depth) defines which years are in the raster.
#' @param pow default 1, exponential weight used in determining year_wt_i:
#' if pow==0, get uniform wts
#' if pow==1, get linear wts
#' if pow > 1, get exponential wts, ...
#' @param uselogit Default FALSE. If set to TRUE, then it converts rasters to logit space before calculations
#'            
#' @return returns a raster of AROC where each pixel in the 1 layer raster the weighted AROC for that pixel determined from the pred.rast
get.aroc <- function(pred.rast, year.map = NULL, pow = 1, uselogit = FALSE){

  ## make the year.map if not supplied
  if(is.null(year.map)){
    year.map <- 1:nlayers(pred.rast)
  }

  ## setup lengths of years and pixels
  n.yr <- length(year.map)
  n.px <- length(as.vector(pred.rast[[1]]))

  ## vectorize the pred.rast and make a matrix of RC between years
  pix.mat <- matrix(ncol = n.yr, nrow = n.px)

  ## convert to logit space?
  if(uselogit == TRUE){
    pix.mat <- log(pix.mat / (1 - pix.mat))
  }

  ## calculate RC, substitute RC between yr_{i+1} and yr_i into pix.mat_i
  message("calculating rates of change")
  for(i in 1:(n.yr - 1)){
    pix.mat[, i] <- log(values(pred.rast[[i + 1]])) - log(values(pred.rast[[i]])) / (year.map[i + 1] - year.map[i])
  }

  ## remove last column of matrix
  pix.mat <- pix.mat[, -n.yr]
  ## now we have a matrix of AROC between years for all pixels in pred.rast

  ## get unstandardized year weights
  year.wt <- ( year.map - min(year.map) ) ^ pow
  ## leave out the earliest year which now has weight zero
  year.wt <- year.wt[-1]
  ## standardize to sum to 1
  year.wt <- year.wt / sum(year.wt)

  ## perform the linear combo of year wts and rate of changes by year
  ## to get weighted AROC pixel estimates
  message("making weighted rates of change")
  pix.aroc <- pix.mat %*% year.wt

  ## and put back into a raster
  message("preparing AROC raster")
  aroc.rast <- pred.rast[[1]]
  values(aroc.rast) <- pix.aroc

  return(aroc.rast)
}


#' @title Forecast Raster with AROC
#' @description Uses an AROC raster (see \code{\link{get.aroc}}) to forecast values out from a starting raster.
#'
#' @param start.rast Raster to start the forecast with. Should be the most recent layer of the raster.
#' @param aroc.rast Raster whose pixels are AROC values
#' @param n.yrs Default 15, number of years to forecast
#' @param uselogit Default FALSE. if TRUE, convert to logit space before forecasting and inverts back to make final maps
#'            
#' @return A forecast raster layer n.yrs after the year of the start.rast
make.forecast <- function(start.rast, aroc.rast, n.yrs = 15, uselogit = FALSE){

  ## we use p*e^{r*t}

  fc.rast <- start.rast

  if(uselogit) fc.rast <- log(fc.rast / (1 - fc.rast))

  values(fc.rast) <- values(start.rast) * exp(values(aroc.rast) * n.yrs)

  if(uselogit) fc.rast <- exp(fc.rast) / (1 + fc.rast)

  return(fc.rast)

}


#' @title Compile Results Table
#' @description DEPRECATED. Seems to be limited to Africa and is tied to some other function that creates a `/table_{year}` subfolder in the run_date folder. Creates a new subdirectory within the run date folder, `/summary_tables` and saves a number of objects there. Seems to rely on files with `_percent2010` in the name to create admin summaries for the probability of meeting a certain threshold or SDG value.
#'
#' @param indicator_group indicator_group of model run
#' @param indicator indicator of model run
#' @param run_date run_date of model run
#' @param measure Default 'mortality'. Used as a regex pattern for finding files in the `table_{year}` run_date subfolder
#' @param baseline_year default 2000. Year in `table_{year}`
#' @param year_to_map default 2015. Year that gets mapped to Africa raster saved out by function
#' @param goal_threshold default 0.001. Replaces mean/upper/lower values in read in files where NA
#' @param metric default 'sdgprob'. Can also be 'percent'.
#' @param shapefile_version default 'current'. Standard shapefile version string. 
#'            
#' @return None
compile_results_table <- function(indicator_group,
                                  indicator,
                                  run_date,
                                  measure = 'mortality',
                                  baseline_year = 2000,
                                  year_to_map = 2015,
                                  goal_threshold = 0.001,
                                  metric = 'sdgprob',
                                  shapefile_version = 'current') {

  results_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/', run_date, '/table_', baseline_year)
  dir.create(paste0(results_dir, '/summary_tables'))
  all_files <- list.files(results_dir, pattern = measure, full.names = TRUE)
  if(metric=='percent') all_files <- all_files[grepl('_percent2010', all_files)]
  if(metric=='sdgprob') all_files <- all_files[!grepl('_percent2010', all_files)]
  all_admins <- rbindlist(lapply(all_files, fread))
  all_admins <- all_admins[is.na(mean), mean := goal_threshold]
  all_admins <- all_admins[is.na(upper), upper := goal_threshold]
  all_admins <- all_admins[is.na(lower), lower := goal_threshold]
  setnames(all_admins, 'admin2', 'ADM2_CODE')
  setnames(all_admins, 'admin1', 'ADM1_CODE')
  setnames(all_admins, 'admin0', 'ADM0_CODE')
  all_gauls <- unique(all_admins[, ADM0_CODE])
  for(admin_level in c(0,1,2)) {
    admin_names <- read.dbf(get_admin_shapefile(admin_level, suffix='.dbf',
                                                version = shapefile_version))
    admin_names <- as.data.table(admin_names)
    admin_names <- admin_names[, c(paste0('ADM', admin_level, '_CODE'), paste0('ADM', admin_level, '_NAME')), with=FALSE]
    all_admins <- merge(all_admins, admin_names, by=paste0('ADM', admin_level, '_CODE'), all.x=TRUE)
  }

  write.csv(all_admins,
            paste0(results_dir, '/summary_tables/', indicator, '_', measure, '_', metric, '_summary_table.csv'))

  ## Save copy of all_admins (by level) for Lucas
  for(adm in c(0,1,2)) {
    if(adm == 0) {
      this_admin_data <- all_admins[is.na(ADM1_CODE) & is.na(ADM2_CODE), ]
      this_admin_data <- this_admin_data[year == year_to_map, ]
    }
    if(adm == 1) {
      this_admin_data <- all_admins[!is.na(ADM1_CODE) & is.na(ADM2_CODE), ]
      this_admin_data <- this_admin_data[year == year_to_map, ]
    }
    if(adm == 2) {
      this_admin_data <- all_admins[!is.na(ADM1_CODE) & !is.na(ADM2_CODE), ]
      this_admin_data <- this_admin_data[year == year_to_map, ]
    }
    write.csv(this_admin_data, paste0(results_dir, '/pixels/adm', adm, '_', indicator, '_', measure, '_', metric, year_to_map, '_summary_table.csv'))
  }

  make_africa_raster <- function(gaul) {

    # Convert raster to SpatialPointsDataFrame
    message(paste0('rasterizing ', gaul, '...'))
    pixel_probs <- fread(paste0(results_dir, '/pixels/', measure, '_', gaul, '_probs.csv'))
    simple_raster <- raster(paste0(results_dir, '/simple/', gaul, '.tif'))
    probs_raster <- insertRaster(simple_raster, matrix(pixel_probs[, p_goal], ncol = 1))
    return(probs_raster)

  }

  africa_raster <- lapply(unique(all_admins[, ADM0_CODE]), make_africa_raster)
  final_africa_raster = do.call(raster::merge, africa_raster)
  writeRaster(final_africa_raster,
              file = paste0(results_dir, '/pixels/', indicator, '_', measure, '_', metric, year_to_map),
              format='GTiff',
              overwrite = TRUE)

}


#' @title Extract from brick
#'
#' @description Extract values from a raster or rasterbrick to lat/lon points by relevant year of data.
#'
#' @author Rebecca Stubbs
#'
#' @param data A data.table or data.frame (that will be converted to a data.table) with
#'             columns for latitude, longitude, and year.
#' @param yearcol A string name of a column that defines the year of the data within the data object
#' @param varname A string name of what you want the column of extracted values to be
#' @param raster_years A vector that describes the years that serve as bands in the rasterbrick, in order.
#' @param rb A rasterbrick with the bands that correspond to the raster_years.
#' @return Returns the data object, sorted based on year, with a new column with the extracted values (name
#' as specified based on the varname parameter).
extract_from_brick<-function(data,
                             yearcol="original_year",
                             rb,
                             varname="extracted_values",
                             raster_years=2000:2015){

  # Making sure you've passed the right arguments to the function:
  if(sum(c("longitude","latitude",yearcol) %in% names(data))!=3){
    stop("You need to have the columns longitude, latitude, and a specified column that describes year (as the yearcol parameter) in the data object passed to this function.")
  }
  if(varname=="extracted_values"){warning("The variable added to the data.table will be named 'extracted_values' unless you supply a different name to the varname parameter.")}

  if(min(data[[yearcol]])<min(raster_years)){stop("You have years of raw input data that extend beyond the min or max of the raster years supplied to the function. Consider setting years outside those bounds to the min/max of the raster years.")}

  # Ordering the data object by year
  data<-data[order(data[[yearcol]])]

  # Renaming the rasterbrick layers to be descriptive of year
  names(rb)<-paste0("year_",as.character(raster_years))

  # Define a vector to store the results
  extracted_values<-c()

  for (yr in raster_years){
    message(paste0("Extracting values for ",yr))
    # Select only the data points from that year
    singleyr<-data[data[[yearcol]]==yr,]

    # Convert to spatial object
    coordinates(singleyr) <- ~longitude + latitude

    # setting the proj4string (the projection information) of this layer to be the
    # same as the results raster
    proj4string(singleyr)<-proj4string(rb)

    # Extract raster values to points
    values<-raster::extract(rb[[paste0("year_",as.character(yr))]], singleyr)
    extracted_values<-c(extracted_values,values)
  }

  data[[varname]]<-extracted_values
  return(data)
}

#' @title Get weight correlation between raw data and unraked estimates
#' @description Returns the weighted correlation between the raw data points input into a model run,
#' and the unraked estimates.
#' @author Rebecca Stubbs
#'
#' @param indicator_group The string of the indicator group.
#' @param indicator The string name of the indicator.
#' @param rundate The string name of the rundate.
#'
#' @return Returns the output from a call to the weights::wtd.cors() function
#' between the raw data points and the estimates.
get_weighted_correlation_raw_estimated<-function(indicator_group,
                                                 indicator,
                                                 rundate){

  source(file.path(fp_list$code_root, 'lbd_core/mbg_central/prep_functions.R'))
  source(file.path(fp_list$code_root, 'lbd_core/mbg_central/setup.R'))
  load_R_packages(c('data.table','raster','sp','weights'))

  input_folder<-file.path(fp_list$mbg_store, 'input_data/')

  message("Reading in the raw data")
  # Read in CSV, reduce table size by eliminating irrelevant columns
  raw_data<-fread(paste0(input_folder,indicator,".csv"))
  raw_data[,raw:=raw_data[[indicator]]/N]
  raw_data<-raw_data[,list(year=original_year,latitude,longitude,raw,N,weight)]

  message("Loading in the mean raster (Unraked)")
  # Load in Raster of Values
  results_tif<-raster::stack(paste0(fp_list['mbg_root'],indicator_group,"/",indicator,"/output/",rundate,"/",indicator,"_mean_raster.tif"))
  names(results_tif)<-as.character(seq(2000,2015))

  raw_estimate_comparison<-list()

  message("Extracting values to points")
  for(yr in seq(2000,2015)){
    print(yr)
    raw_singleyr<-raw_data[year==yr,]

    # Convert to spatial object
    coordinates(raw_singleyr) <- ~longitude + latitude

    # setting the proj4string (the projection information) of this layer to be the
    # same as the results raster
    proj4string(raw_singleyr)<-proj4string(results_tif)

    # Extract raster values to points
    values<-raster::extract(results_tif[[paste0("X",as.character(yr))]], raw_singleyr)

    # Going back into table-space, and adding the data.table to a list to rbind later
    raw_singleyr<-as.data.table(raw_singleyr)
    raw_singleyr[,estimate:=values]
    raw_estimate_comparison[[as.character(yr)]]<-raw_singleyr
  }

  message("Calculating correlation coefficient using the weights::wtd.cors function")
  # Combining together the raw estimates into 1 data.table
  raw_estimate_comparison<-rbindlist(raw_estimate_comparison)
  cor_results<-weights::wtd.cors(x=raw_estimate_comparison$raw,
                                 y=raw_estimate_comparison$estimate,
                                 weight=raw_estimate_comparison$N * raw_estimate_comparison$weight)

  return(cor_results)
}


#' @title Prep Objects for Postestimation
#'
#' @description Save out a temporary .RData file that gets read in during postestimation. Primarily used to pass through GBD raking targets from the model launch script to the postestimation script. File gets saved in the model run folder in `/temp_post_est/post_est_temp_objs.RData`.
#'
#' @param indicator indicator for model run
#' @param indicator_group indicator_group for model run
#' @param run_date run_date for model run
#' @param save_objs vector of names of objects to save e.g. c("gbd"). Can pass any objects that exist in the global environment to this arg. Recommended not to pass through config arguments, because it can overwrite values accessed through the config. 
#' 
#' @return None
prep_postest <- function(indicator,
                         indicator_group,
                         run_date,
                         save_objs) {

  # Save a list of objects in a standard location for parallel scripts to pull from
  main_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/', run_date, '/')
  temp_dir <- paste0(main_dir, "temp_post_est/")
  temp_file <- paste0(temp_dir, "post_est_temp_objs.RData")
  dir.create(temp_dir, showWarnings = F)
  save(list = save_objs, file = temp_file)
}


#' @title Load, Combine Regions, and Save Post-Estimation Objects
#'
#' @description Run after Raking/Aggregation to combine objects across modeling regions. Function specifically combines rasters and raking factors.
#'
#' @param regions default `strata`, character vector of modelling regions to combine
#' @param summstats default c("mean", "cirange", "upper", "lower"). Which summary rasters to combine
#' @param raked default c("raked", "unraked"). Combine raked/unraked rasters?
#' @param rf_table default TRUE. Combine raking factors?
#' @param run_summ defualt TRUE. Make a run summary graph? See \code{\link{graph_run_summary}}
#' @param indic default `indicator`. Indicator for model run
#' @param ig default `indicator_group`. Indicator_group for model run
#' @param sdir default `sharedir`
#' @param proj default FALSE. If T, Combine projected rasters
#' @param proj_folder default NULL. Unused. 
#' 
#' @return None
post_load_combine_save <- function(regions = strata,
                                   summstats = c("mean", "cirange", "upper", "lower"),
                                   raked = c("raked", "unraked"),
                                   rf_table = TRUE,
                                   run_summ = TRUE,
                                   indic = indicator,
                                   ig = indicator_group,
                                   sdir = sharedir,
                                   proj = FALSE,
                                   proj_folder = NULL) {

  message(paste0("indic: ", indic))
  message(paste0("ig: ", ig))

  rake_addin <- character()
  if ("unraked" %in% raked) {
    lookup_dir <- paste0(sprintf("share/geospatial/mbg/%s/%s/output/%s/", ig, indic, run_date))
    ur <- length(grep(paste0(indic, ".*unraked.*raster.tif"), list.files(lookup_dir)))
    if(proj) ur <- length(grep(paste0(indic, ".*unraked_PROJ.*raster_PROJ.tif"), list.files(lookup_dir)))
    if (ur > 0) rake_addin <- c(rake_addin, unraked = "_unraked")
    if (ur == 0) rake_addin <- c(rake_addin, unraked = "")
  }

  if ("raked" %in% raked) {
    rake_addin <- c(rake_addin, raked = "_raked")
  }

  # loop through and combine all rasters
  message("\nCombining rasters...")
  for (rake in rake_addin) {
    message(names(rake_addin)[which(rake_addin == rake)])
    rr <- rake
    for (ss in summstats) {
      message(paste0("  ", ss))
      rlist <- list()
      for (reg in regions) {
        message(paste0("    ", reg))
        rlist[[reg]] <-
          brick(ifelse(proj,
                       file.path(fp_list['mbg_root'], ig, indic, 'output', run_date,
                                 sprintf("%s_%s%s_%s_raster_PROJ.tif", indic, reg, rake, ss)),
                       file.path(fp_list['mbg_root'], ig, indic, 'output', run_date,
                                 sprintf("%s_%s%s_%s_raster.tif", indic, reg, rake, ss))
          ))
      }
      if (length(rlist) > 1) rlist <- do.call(raster::merge, unname(rlist)) else rlist <- rlist[[1]]
      if (ss == "cirange") ssname <- "range" else ssname <- ss # naming convention
      save_post_est(
        rlist, "raster",
        ifelse(!proj,
               paste0(ssname, rr, "_raster"),
               paste0(ssname, rr, "_raster_PROJ")),
        indic
      )
    }
  }

  # do rf also
  if (rf_table) {
    message("RF table")
    rflist <- list()
    for (reg in regions) {
      rflist[[reg]] <-
        if(proj) {
          read.csv(file.path(fp_list['mbg_root'], ig, indic, 'output', run_date,
                             sprintf("%s_%s_rf_PROJ.csv", indic, reg)))
        } else {
          read.csv(file.path(fp_list['mbg_root'], ig, indic, 'output', run_date,
                             sprintf("%s_%s_rf.csv", indic, reg)))
        }
    }
    if(!proj) {
      save_post_est(do.call(rbind.fill, rflist), "csv", "rf", indic)
    } else {
      save_post_est(do.call(rbind.fill, rflist), "csv", "rf_PROJ", indic)
    }

  }

  # make a run summary graph
  if (run_summ) {
    graph_run_summary(
      run_date = run_date,
      indicator_group = ig,
      indicator = indic
    )
  }
}


#' @title Cleanup After Postestimation
#'
#' @description Remove `temp_post_est` directory and region rasters after \code{\link{post_load_combine_save}} runs. 
#'
#' @param indicator indicator of model run
#' @param indicator_group indicator_group of model run
#' @param run_date run_date of model run
#' @param strata if `delete_region_rasters == T`, a vector of region names to delete.
#' @param delete_region_rasters default F. If T, remove region rasters.
#' 
#' @return None
clean_after_postest <- function(indicator,
                                indicator_group,
                                run_date,
                                strata,
                                delete_region_rasters = F) {

  # Delete intermediate files that we no longer need
  main_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/', run_date, '/')
  temp_dir <- paste0(main_dir,'temp_post_est/')

  # Deleting - be careful!
  unlink(temp_dir, recursive = T)

  # If desired, insert code here to delete other temporary objects (eg. region-specific rasters)
  grep_string <- paste0(indicator, "_(", paste(strata, collapse = "|"), ").*_raster.tif")
  region_rasters <- grep(grep_string,list.files(main_dir), value=T) %>%
    paste0(main_dir, .)
  if(delete_region_rasters == T) unlink(region_rasters)
}



#' @title Get GBD age group for the given WorldPop measure
#' @description This function will return the age group from the GBD age group ID table
#' which corresponds to the WorldPop measures used
#'
#' @note Currently only taking in measures starting with "a" and for 12 month interval (ending with "t"). And \code{total}.
#'
#' @param pop_measure The WorldPop measure, e.g. \code{a0004t}. Default: \code{a0004t}
#'
#' @return The \code{age_group_id} for the associated \code{pop_measure} input.
#'
#' @export
get_age_group_from_worldpop <- function(pop_measure) {

  ## Strsplit out the 4 digits from the middle (always fixed length),
  ## only if its not total
  age_measures <- pop_measure
  if (!pop_measure %in% c("total", "wocba")) {
    age_measures <- substr(x = pop_measure, start = 2, stop = 5)
  }

  ## NOT ALL the ages have perfect 1:1 measure, and so we may
  ## have to aggregate out some age groups

  ## Here's a dictionary of age mappings
  age_dict <- list(
    "0004" = c(2, 3, 4, 5),
    "0514" = c(23),
    "1014" = c(7),
    "1519" = c(8),
    "1549" = c(24),
    "2024" = c(9),
    "2529" = c(10),
    "3034" = c(11),
    "3539" = c(12),
    "4044" = c(13),
    "4549" = c(14),
    "5054" = c(15),
    "5559" = c(16),
    "6064" = c(17),
    "total" = c(22),
    "wocba" = c(24)
  )

  ## Lookup the GBD age group for the given WorldPop measure and return it
  return(unlist(age_dict[paste0(age_measures)]))
}


#' @title Get GBD sex_id for the given WorldPop measure
#' @description This function will return the sex_id by parsing the final character of standard worldpop measures
#' which end in either "m", "f", or "t", representing male, female, or both sexes. Also accepts "total" and "wocba" measures.
#'
#' @param pop_measure The WorldPop measure, e.g. \code{a0004t}.
#'
#' @return The \code{sex_id} for the associated \code{pop_measure} input.
get_sex_id_from_worldpop <- function(pop_measure) {
  
  #define sex group for worldpop measures that don't match the standard format
  sex_dict <- list(
    "total" = 3,
    "wocba" = 2
  )
  
  if(pop_measure %in% names(sex_dict)) {
    return(unname(unlist(sex_dict[pop_measure])))
  }
  
  if(endsWith(pop_measure, "m")) {
    return(1)
  } else if (endsWith(pop_measure, "f")) {
    return(2)
  } else if (endsWith(pop_measure, "t")) {
    return(3)
  } else {
    message("could not identify GBD sex_id from worldpop, using 3 (both sexes) as default")
    return(3)
  }
}

#' @title Extrapolate GBD measure to future
#' @description Use a simple loess model to extrapolate GBD input data to desired year
#'
#' @param gbd The data.table with name, year, and mean values (output from \code{get_gbd_estimates})
#' @param year_list List of in sample years (config argument)
#' @param year_forecast_end Year to forecast to. Default: 2030
#'
#' @return A data.table with name, year and mean columns extrapolated to \code{year_forecast_end}
#'
#' @importFrom stats loess
#' @export
loess_extrap_gbd_input_by_loc <- function(gbd, year_list, year_forecast_end = 2030) {
  gbd_forecast <- data.table(expand.grid(name = unique(gbd$name),
                                         year = union(year_list, c(max(year_list):year_forecast_end))))
  gbd_forecast <- merge(gbd_forecast, gbd, c('name', 'year'), all.x = TRUE)

  ## Loop over all locations
  lapply(unique(gbd$name), function(loc) {

    ## Run a simple loess and extrapolate in log space
    gbd_fc_vector <- predict(stats::loess(formula = log(mean) ~ year,
                                          data = gbd_forecast[name == loc],
                                          control = loess.control(surface = "direct")),
                             newdata = gbd_forecast[name == loc])
    gbd_forecast[name == loc, mean_forecast:= gbd_fc_vector]

    ## Intercept shift at last year of in-sample data
    gbd_forecast[year == max(year_list), int_shift:= log(mean) - mean_forecast]
    gbd_forecast[name == loc, int_shift:= mean(int_shift, na.rm = TRUE), by = 'name']
    gbd_forecast[name == loc & year >= max(year_list), mean:= exp(mean_forecast + int_shift)]
    gbd_forecast[, c('mean_forecast', 'int_shift'):= NULL]
    return(0)
  })

  return(gbd_forecast)
}


#' @title Get pre-saved FHS population outputs
#' @description Access pre-saved RDS population outputs from FHS (saved out in \code{<geospatial_root>/fhs-outputs/population})
#'
#' @param population_version A version tag found from looking into: \code{<geospatial_root>/fhs-outputs/population}
#' Default: \code{"20190403_test_new_cluster_1000d_rerun_fix_draw_squeezed_agg_ordered"}
#' @param pop_measure WorldPop measure; only one allowed for now! Default: a0004t
#' @param sex_id Sex ID. Default: 3
#' @param scenario FHS scenario. Default: 0 (reference).
#' @param year_ids. Year_id to query. Default: NULL (get all years)
#' @param gbd_regions Regions to query (GBD Location IDs). Default: NULL (get all regions)
#'
#' @return A data.table with location_id, year_id, age_group_id = pop_measure, sex_id, run_id, population
#'
#' @export
get_fhs_population <- function(population_version = "20190403_test_new_cluster_1000d_rerun_fix_draw_squeezed_agg_ordered",
                               pop_measure = "a0004t",
                               sex_ids = 3,
                               scenarios = 0,
                               year_ids = NULL,
                               gbd_regions = NULL) {

  ## Get pop data
  pop_data <- readRDS(paste0(fp_list$geospatial_root, "fhs-outputs/population/", population_version, ".rds"))

  ## Keep on subsettin'
  pop_data <- pop_data[(age_group_id %in% get_age_group_from_worldpop(pop_measure = pop_measure)) & (scenario %in% scenarios) & (sex_id %in% sex_ids)]

  if(!is.null(year_ids)) {
    pop_data <- pop_data[year_id %in% year_ids]
  }
  if(!is.null(gbd_regions)) {
    pop_data <- pop_data[location_id %in% gbd_regions]
  }

  ## Return data with alignment of our choice, keeping unique summed up population entries
  pop_data[, age_group_id:= NULL]
  pop_data[, age_group_id:= pop_measure]
  pop_data <- pop_data[, .(population = sum(population)), by = c('location_id', 'year_id', 'age_group_id', 'sex_id')]

  ## Make fake run_id column
  pop_data[, run_id:= paste0(population_version)]

  setkeyv(pop_data, c('location_id', 'year_id', 'age_group_id', 'sex_id'))
  return(pop_data)

}