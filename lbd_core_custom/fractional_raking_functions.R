################################################################################
# these functions to accomodate an alternative workflow for a rates cell_pred
# which also makes admin aggregations of rates and counts.
# 
# NB these functions work for National and Sub-National linear raking of rates only.  
# They do not directly support a counts cell_pred (see Michael Collison's work on 
# fractional raking above) or logit raking.  It is also untested with the GADM 
# shapefiles. (the HIV team will be doing that in the future but it hasn't happend yet)
################################################################################

#' Rakes a cell_pred containing rates to GBD targets
#' @title fractional_rake_rates
#' @description  takes a rates cell_pred %>% 
#'               links it to a fractional raking and aggregation link table %>%
#'               adds population per fractional cell %>%
#'               rakes that population to ensure matching to GBD %>%
#'               saves those raking factors %>%
#'               uses that population to create a linked counts cell_pred %>%
#'               aggregates the cell_pred to the GBD geography level %>%
#'               converts the aggregations back to rates and rakes to GBD tagets %>%
#'               saves those raking fators.
#'
#' @param cell_pred cell_pred object from mbg models.  Each cell must have a rate in it
#' @param simple_raster the simple raster that the cell_pred is based on
#' @param simple_polygon the simple polygon that the cell_pred is based on
#' @param pixel_id list of the pixels in the simple raster that have non na values
#' @param shapefile_version which shapefile geographies are being used 
#' @param reg the modeling region
#' @param pop_measure the worldpop agegroup on which the model is built
#' @param year_list the modeled years
#' @param use_intermediate_years Boolean to indicate whether or not to rake to intermediate years. Default: TRUE
#' @param interval_mo the time in months between the modeled years
#' @param rake_subnational a logical value indicating the use of subnational raking targets or not
#' @param age_group the gbd age group that the model is built on
#' @param sex_id the gbd sex group that the model is built on
#' @param sharedir sharedir       <- '<<<< FILEPATH REDACTED >>>>'
#' @param run_date model run date
#' @param indicator modeled indicator
#' @param gbd gbd object prepared containing the raking targets
#' @param gbd_pops output from central code "get_population" function
#' @param countries_not_to_rake countries (vector or summed string) to not rake to GBD (we set rake factor to 1 for those)
#' @param countries_not_to_subnat_rake character vector of iso3 codes for countries not to subnationally rake
#' @param custom_output_folder Output the rake factors and outputs to custom folder path if specified. Default: NULL
#' @param rake_method if set to "logit" creates raking factors in logit space, otherwise assumes linear raking
#'
#' @return automatically saves out two raking factors tables, one for population and one for the indicator
#' @export


fractional_rake_rates <- function(cell_pred = cell_pred,
                                  simple_raster = simple_raster,
                                  simple_polygon = simple_polygon,
                                  pixel_id = pixel_id,
                                  shapefile_version = shapefile_version,
                                  reg = reg,
                                  pop_measure = pop_measure,
                                  year_list = year_list,
                                  use_intermediate_years = TRUE,
                                  interval_mo = interval_mo,
                                  rake_subnational = rake_subnational,
                                  age_group = age_group,
                                  sex_id = sex_id,
                                  sharedir = sharedir,
                                  run_date = run_date,
                                  indicator = indicator,
                                  gbd = gbd,
                                  rake_method = "linear",
                                  gbd_pops = gbd_pops,
                                  countries_not_to_rake = NULL,
                                  countries_not_to_subnat_rake = NULL,
                                  custom_output_folder = NULL,
                                  hold = holdout,
                                  meas = 'prevalence',
                                  etiology = 'none'){
  
  # setting a reference for the number of draws
  ndraws = ncol(cell_pred)
  
  #####################################################################
  #load the cell id to admin units link
  link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)
  
  #####################################################################
  # collect and load the population data from the WorldPop rasters
  covdt <- load_populations_cov(reg, pop_measure, measure = 'count', simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
  
  if(!use_intermediate_years) {
    print("Subsetting to supplied years only")
    covdt <- covdt[year %in% year_list]
  }
  
  #####################################################################
  # Prepping the cell_pred and link table to be linked by making sure they have the appropriate identifiers.  Also performs a 
  # zippering at the region boundary where cells that have a portion of their area outside of the modeling region are reintegrated
  # as a whole cell and considered to be only in one region.  This works becasue a given cell is only modeled in one region.
  link <- prep_link_table(link_table = link_table,
                          simple_raster = simple_raster,
                          pixel_id = pixel_id)
  
  cell_ids <- link_table[[2]]
  
  # getting the connector for sub-national or national raking, This connector gets the IHME location_code for our
  # gbd targets and connects that to the ADM0_CODE or ADM1_CODE as nessecary 
  
  connector <- get_gbd_locs(rake_subnational = rake_subnational,
                            reg = reg,
                            shapefile_version = shapefile_version)
  
  # getting the connector for sub-national raking - used to implement countries_not_to_subnat_rake
  nat_connector <- get_gbd_locs(rake_subnational = F,
                                reg = reg,
                                shapefile_version = shapefile_version)
  
  # merge the connectors on to the link table
  link <- sub_nat_link_merge(rake_subnational,
                             link,
                             connector,
                             nat_connector,
                             countries_not_to_subnat_rake)
  
  # set cell pred as a data table, and rename things
  cell_pred <- prep_cell_pred(cell_pred = cell_pred,
                              cell_ids  = cell_ids,
                              pixel_id  = pixel_id,
                              covdt     = covdt)
  
  
  # merge cell_pred on the link
  cell_pred = merge(link, cell_pred, by.x = 'ID', by.y = 'cell_id',allow.cartesian = TRUE)
  
  # space
  link <- NULL
  
  ## Raking Population ###################################################################
  # This is done to ensure that the total pop in each raking geography is the same as GBD
  message("raking population")
  
  #convert to fractional population 
  cell_pred = cell_pred[,pop := pop * area_fraction] 
  
  #NA out population where the pixel value is NA (to prevent weirdness with denominators)
  cell_pred = cell_pred[is.na(V1), pop := NA]
  
  scalars <- calculate_pop_scalars(cell_pred = cell_pred,
                                   age_group = age_group,
                                   connector = connector,
                                   sex       = sex_id,
                                   sharedir  = sharedir,
                                   run_date  = run_date,
                                   indicator = indicator,
                                   stratum   = reg,
                                   gbd_pops  = gbd_pops,
                                   custom_output_folder = custom_output_folder)
  # add back to the cell_pred as a population rf
  cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))
  
  # rake fractional populations
  cell_pred$pop_raked <- 0
  cell_pred = cell_pred[,pop_raked := pop * pop_scalar]
  cell_pred$pop <- NULL
  
  ## Raking actual data ###################################################################
  message("raking rates")
  # Calculate Fractional Raking Factors
  fractional_rf <- calculate_fractional_rfs(ndraws    = ndraws,
                                            cell_pred = cell_pred,
                                            gbd       = gbd,
                                            sharedir  = sharedir,
                                            run_date  = run_date,
                                            indicator = indicator,
                                            shapefile_version = shapefile_version,
                                            stratum   = reg,
                                            countries_not_to_rake = countries_not_to_rake,
                                            custom_output_folder = custom_output_folder,
                                            countries_not_to_subnat_rake = countries_not_to_subnat_rake,
                                            rake_method = rake_method,
                                            h = hold,
                                            m = meas,
                                            e = etiology)
}


#' Fractional aggregation of a rates cell pred
#' @title Fractional aggregation of a rates cell pred
#' @description takes a rates cell_pred =>
#'               links it to a fractional raking and aggregation link table =>
#'               rakes the rates based on the raking factors above =>
#'               saves raked rates cell_pred =>
#'               adds population per fractional cell =>
#'               rakes that population to ensure matching to GBD using pop raking factors =>
#'               uses that population to create a counts cell_pred =>
#'               saves that raked counts cell_pred =>
#'               aggregates the cell_pred to ADM 0, 1, 2 levels =>
#'               saves aggregations at ADM 0,1,2 levels, raked and unraked, rates and counts
#' @param cell_pred cell_pred object from mbg models.  Each cell must have a rate in it.  This needs to be unraked
#' @param simple_raster the simple raster that the cell_pred is based on
#' @param simple_polygon the simple polygon that the cell_pred is based on
#' @param pixel_id list of the pixels in the simple raster that have non na values
#' @param shapefile_version which shapefile geographies are being used
#' @param reg the modeling region
#' @param pop_measure the worldpop agegroup on which the model is built
#' @param year_list the modeled years
#' @param use_intermediate_years Boolean to indicate whether or not to rake to intermediate years. Default: TRUE
#' @param interval_mo the time in months between the modeled years
#' @param rake_subnational a logical value indicating the use of subnational raking targets or not
#' @param sharedir Path to sharedir e.g.: \code{sharedir <- '<<<< FILEPATH REDACTED >>>>'}
#' @param run_date model run date string
#' @param indicator modeled indicator
#' @param main_dir Path to main_dir e.g.: \code{main_dir <- '<<<< FILEPATH REDACTED >>>>'}
#' @param age age holdout group, assumed to be 0 but left flexible
#' @param holdout n fold hold out group, assumed to be 0 but left flexible
#' @param countries_not_to_subnat_rake character vector of iso3 codes for countries not to subnationally rake
#' @param return_objects Return cell pred and raking factors into environments? Default: FALSE
#' @param custom_output_folder Get the rake factors from this folder, and also save the Rdata files here. Default: NULL
#' @param rake_method if set to "logit" creates raking factors in logit space, otherwise assumes linear raking
#'
#' @return automatically saves unraked and raked versions admin aggregations of both rates and counts,
#' as well as raked cell_pred objects of both rates and counts. If desired, one can also return raked_cell_pred
#' and the raking factors into environment.
#'
#' @note the raked cell_pred objects cannot be used to create raked aggregation
#' @note you need to start with an unraked cell_pred, link it, rake it, then aggregate it.
#' @export
fractional_agg_rates <- function(cell_pred = cell_pred,
                                 simple_raster = simple_raster,
                                 simple_polygon = simple_polygon,
                                 pixel_id = pixel_id,
                                 shapefile_version = shapefile_version,
                                 reg = reg,
                                 pop_measure = pop_measure,
                                 year_list = year_list,
                                 use_intermediate_years = TRUE,
                                 interval_mo = interval_mo,
                                 rake_subnational = rake_subnational,
                                 sharedir = sharedir,
                                 run_date = run_date,
                                 indicator = indicator,
                                 main_dir = main_dir,
                                 rake_method = "linear",
                                 age = 0,
                                 holdout = 0,
                                 return_objects = FALSE,
                                 countries_not_to_subnat_rake = countries_not_to_subnat_rake,
                                 custom_output_folder = NULL,
                                 meas = 'prevalence',
                                 etiology = 'none') {
  # setting a reference for the number of draws
  ndraws <- ncol(cell_pred)
  region <- reg
  #####################################################################
  # load the cell id to admin units link
  
  link_table <- get_link_table(simple_raster, shapefile_version = shapefile_version)
  
  #####################################################################
  # collect and load the population data
  covdt <- load_populations_cov(reg, pop_measure, measure = "count", simple_polygon, simple_raster, year_list, interval_mo, pixel_id = pixel_id)
  
  if(!use_intermediate_years) {
    print("Subsetting to supplied years only")
    covdt <- covdt[year %in% year_list]
  }
  
  #####################################################################
  # Prepping the cell_pred and link table to be linked and then merging them
  link <- prep_link_table(
    link_table = link_table,
    simple_raster = simple_raster,
    pixel_id = pixel_id
  )
  
  cell_ids <- link_table[[2]]
  
  # getting the connector for sub-national raking
  connector <- get_gbd_locs(
    rake_subnational = rake_subnational,
    reg = reg,
    shapefile_version = shapefile_version
  )
  
  # getting the connector for sub-national raking
  nat_connector <- get_gbd_locs(
    rake_subnational = F,
    reg = reg,
    shapefile_version = shapefile_version
  )
  
  # merge the connector on to the link table
  link <- sub_nat_link_merge(
    rake_subnational,
    link,
    connector,
    nat_connector,
    countries_not_to_subnat_rake
  )
  
  # set cell pred as a data table, and rename things
  cell_pred <- prep_cell_pred(
    cell_pred = cell_pred,
    cell_ids = cell_ids,
    pixel_id = pixel_id,
    covdt = covdt
  )
  
  # merge on the link
  cell_pred <- merge(link, cell_pred, by.x = "ID", by.y = "cell_id", allow.cartesian = TRUE)
  
  
  ############################################################
  # adding the raking factors and scaling the populations
  
  message("adding raking factors")
  # convert to fractional population
  cell_pred <- cell_pred[, pop := pop * area_fraction]
  
  scalars <- read.csv(file = ifelse(is.null(custom_output_folder),
                                    paste0(sharedir, "/output/", run_date, "/", indicator, "_", reg, "_pop_rf.csv"),
                                    paste0(custom_output_folder, "/", indicator, "_", reg, "_pop_rf.csv")
  ))
  
  # add back to the cell_pred as a population rf
  cell_pred <- merge(cell_pred, scalars, by = c("location_id", "year"))
  
  #set filename add-in if using etiology
  if (!(etiology == 'none')){
    e_string <- paste0('_',etiology,'_')
  } else (e_string <- '_')
  
  ## load the raking factors
  fractional_rf <- read.csv(file = ifelse(is.null(custom_output_folder),
                                          paste0(sharedir, "/output/", run_date, "/", indicator, e_string, reg, "_", meas, "_rf_", holdout, ".csv"),
                                          paste0(custom_output_folder, "/", indicator, e_string, reg, "_", meas, "_rf_", holdout, ".csv")
  ))
  
  ## merge them onto the cell_pred
  cell_pred <- merge(cell_pred, fractional_rf, by = c("location_id", "year"))
  
  
  ############################################################
  # creating a raked rates cell_pred (this happens first becasue once we go to counts in the cell_pred we can't do back without loading a fresh copy)
  message("creating a raked rates cell_pred")
  # rake rates
  overs <- paste0("V", 1:ndraws)
  
  if(rake_method == "linear"){
    cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * rf)]
  }else{
    cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) invlogit(logit(get(x)) + rf))]
  }
  
  # multiply the cell_pred by the area fraction for the dedupe function (so that each cell will add to 1 and the constituent rates are weighted by area)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * area_fraction) ]
  
  raked_cell_pred <- dedupe_linked_cell_pred(cell_pred, overs)
  
  # Save raked rates cell_pred object
  save(raked_cell_pred, file = ifelse(is.null(custom_output_folder),
                                      paste0(
                                        sharedir, "/output/", run_date, "/",
                                        indicator, e_string, "cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                      ),
                                      paste0(
                                        custom_output_folder, "/",
                                        indicator, e_string, "cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                      )
  ))
  
  
  # SPACE!!!!!
  if (!return_objects) {
    raked_cell_pred <- NULL
    gc()
  }
  
  # un do the area fraction (so that you can use this cell pred again)
  overs <- paste0("V", 1:ndraws)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) / area_fraction) ]
  
  ############################################################
  # creating a raked counts cell_pred
  message("creating a raked counts cell_pred")
  # rake fractional populations
  cell_pred$pop_raked <- 0
  cell_pred <- cell_pred[, pop_raked := pop * pop_scalar]
  cell_pred$pop <- NULL
  
  message("converting from prevalence to counts")
  # set the variables to aggregate
  overs <- paste0("V", 1:ndraws)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * pop_raked)]
  
  # NA out population where the pixel value is NA (to prevent weirdness with denominators)
  cell_pred <- cell_pred[is.na(V1), pop_raked := NA]
  raked_cell_pred_c <- dedupe_linked_cell_pred(cell_pred, overs)
  save(raked_cell_pred_c, file = ifelse(is.null(custom_output_folder),
                                        paste0(
                                          sharedir, "/output/", run_date, "/",
                                          indicator, e_string, "c_cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                        ), paste0(
                                          custom_output_folder, "/",
                                          indicator, e_string, "c_cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"
                                        )
  ))
  
  # SPACE
  raked_cell_pred_c <- NULL
  
  ############################################################
  # creating a raked counts aggregations
  message("creating a raked counts aggregations")
  
  # do calculations!
  admin_2 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM2_CODE", "ADM0_CODE")]
  admin_1 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM1_CODE", "ADM0_CODE")]
  admin_0 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE")]
  
  
  setnames(admin_2, grep("V[0-9]", names(admin_2), value = T), c(overs, "pop_raked"))
  setnames(admin_1, grep("V[0-9]", names(admin_1), value = T), c(overs, "pop_raked"))
  setnames(admin_0, grep("V[0-9]", names(admin_0), value = T), c(overs, "pop_raked"))
  
  # create the spatial hierarchy
  sp_hierarchy_list <- unique(link[ADM0_CODE %in% unique(admin_0[, ADM0_CODE]), .(ADM0_CODE, ADM1_CODE, ADM2_CODE, ADM0_NAME, ADM1_NAME, ADM2_NAME, region)])
  sp_hierarchy_list[, region := reg]
  
  # cleaning raked admin draws in count space
  admin_0$pop <- admin_0$pop_raked
  admin_0$pop_raked <- NULL
  admin_1$ADM0_CODE <- NULL
  admin_1$pop <- admin_1$pop_raked
  admin_1$pop_raked <- NULL
  admin_2$ADM0_CODE <- NULL
  admin_2$pop <- admin_2$pop_raked
  admin_2$pop_raked <- NULL
  
  ## save raked counts aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     paste0(main_dir, "/", indicator, e_string, "raked", "_", meas, "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, e_string, "raked", "_", meas, "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
       )
  )
  
  ############################################################
  # creating raked rates aggregations (you can work back at the admin level because there are no admin's with pop = 0)
  message("creating a raked rates aggregations")
  
  # convert back to rates/prevalence
  overs <- paste0("V", 1:ndraws)
  admin_0 <- admin_0[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_1 <- admin_1[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_2 <- admin_2[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  
  ## save raked rates aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     paste0(main_dir, "/", indicator, e_string, "raked", "_", meas, "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, e_string, "raked", "_", meas, "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
       )
  )
  
  ############################################################
  # creating unraked counts aggregations (can be done two ways, un raking the counts cell_pred or reloading the unraked cell_pred and multiplying by the pop.  I chose this to avoid reloading cell_preds)
  message("creating a unraked counts aggregations")
  
  # unrake all draws
  overs <- paste0("V", 1:ndraws)
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) / rf) ]
  
  ## Create unraked counts agregation
  admin_2 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM2_CODE", "ADM0_CODE")]
  admin_1 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM1_CODE", "ADM0_CODE")]
  admin_0 <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "ADM0_CODE")]
  
  setnames(admin_2, grep("V[0-9]", names(admin_2), value = T), c(overs, "pop_raked"))
  setnames(admin_1, grep("V[0-9]", names(admin_1), value = T), c(overs, "pop_raked"))
  setnames(admin_0, grep("V[0-9]", names(admin_0), value = T), c(overs, "pop_raked"))
  
  admin_0$pop <- admin_0$pop_raked
  admin_0$pop_raked <- NULL
  admin_1$pop <- admin_1$pop_raked
  admin_1$pop_raked <- NULL
  admin_1$ADM0_CODE <- NULL
  admin_2$pop <- admin_2$pop_raked
  admin_2$pop_raked <- NULL
  admin_2$ADM0_CODE <- NULL
  gc()
  
  ## save unraked counts aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     paste0(main_dir, "/", indicator, e_string, "unraked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, e_string, "unraked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
       )
  )
  
  
  ############################################################
  # creating unraked rates aggregations
  message("creating a unraked rates aggregations")
  
  # convert back to rates
  overs <- paste0("V", 1:ndraws)
  admin_0 <- admin_0[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_1 <- admin_1[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  admin_2 <- admin_2[, (overs) := lapply(overs, function(x) get(x) / pop) ]
  
  ## save unraked rates aggregations
  save(admin_0, admin_1, admin_2, sp_hierarchy_list,
       file = ifelse(is.null(custom_output_folder),
                     paste0(main_dir, indicator, e_string, "unraked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, e_string, "unraked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
       )
  )
  
  ## Return cell pred and raking factors if desired
  if (return_objects) {
    output_list <- list()
    output_list[["rf"]] <- data.table(fractional_rf)
    output_list[["raked_cell_pred"]] <- raked_cell_pred
    return(output_list)
  }
}

#' @title Calculate fractional raking factor
#' @description This calculates the raking factors, it takes a rates cell_pred which is the main difference between it and the
#' get_fractional_rf function above (see Michael Collison's work).
#'
#' @param ndraws the number of draws used in this model run, can be calculated from the number of columns in a raw cell_pred
#' @param cell_pred the preped and linked cell_pred that has also had its population raked.
#' @param gbd a data frame of the gbd targets that we are raking to
#' @param sharedir the share directory for the indicator group, used to save the raking factors
#' @param run_date your run date
#' @param indicator the modeled indicator
#' @param shapefile_version the shapefile version
#' @param stratum  the region you are modeling over, used for saving the raking factors
#' @param countries_not_to_rake countries (vector or summed string) to not rake to GBD (we set rake factor to 1 for those)
#' @param custom_output_folder Output the rake factors to custom folder path if specified. Default: NULL
#' @param countries_not_to_subnat_rake countries (vector or summed string) to not rake to subnationals (we set rake factor to 1 for those)
#' @param MaxJump default `10`. Maximum size of a jump to the answer (for logit raking) - this will likely never need to be changed
#' @param MaxIter default `80`. Number of jumps towards the solution (for logit raking) - this will likely never need to be changed
#' @param FunTol default `1e-5`. Maximum allowed difference between the raking target and raked results (for logit raking) - this will likely never need to be changed
#' @param iterate default `F`. If logit raking for a location-year fails, try again with `MaxJump` and `MaxIter` times 10. If that fails, try again times 100. For circumstances where raking target is very far from estimate and raking does not converge.
#' @param zero_heuristic default `F`.  If logit raking, this will automatically set rf = -9999 for any country-year with a target value of 0.  This produces a raked value of 0 for all pixels.  Raking to a target of zero in logit space is very time-consuming and the algorithm
#'                       can only approach zero asymptotically.  For situations where the target is truly zero (most useful for, say, an intervention pre-introduction) this will both speed up the process and ensure that zeros are returned.
#' @param approx_0_1 default `F`. If logit raking, any values of zero will be replaced with 1e-10 and values of 1 will be replaced with (1-(1e-10)).  Otherwise, logit transformations will fail in `NewFindK`. Useful if some areas have very low or high predicted values in `cell_pred`, 
#'                   such that some draws are either 0 or 1
#' @param rake_method if set to "logit" creates raking factors in logit space, otherwise assumes linear raking
#'
#' @return a table of raking facotrs for each of the raking geographies used.  This is so that
#' the fractionally aggregated mbg estimates rasters for a given geography year will equal
#' the GBD estimates for that geography in that year.  This also assumes that the GBD is internally consistent
#' in converting rates to counts, namely rate*pop=count
#'
#' @export

calculate_fractional_rfs <- function(ndraws = ndraws,
                                     cell_pred = cell_pred,
                                     gbd = gbd,
                                     sharedir = sharedir,
                                     run_date = run_date,
                                     indicator = indicator,
                                     shapefile_version = shapefile_version,
                                     stratum = stratum,
                                     countries_not_to_rake = NULL,
                                     custom_output_folder = NULL,
                                     countries_not_to_subnat_rake = NULL,
                                     MaxJump = 10,
                                     MaxIter = 80, 
                                     FunTol = 1e-5,
                                     approx_0_1 = F,
                                     zero_heuristic = F,
                                     iterate = F,
                                     rake_method = "linear",
                                     h = hold,
                                     m = 'prevalence',
                                     e = NULL) {
  message("converting from prevalence to counts")
  # set the variables to aggregate
  overs <- paste0("V", 1:ndraws)
  
  # convert to counts
  cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) * pop_raked)]
  
  # do calculations!
  rake_geo <- cell_pred[, lapply(c(overs, "pop_raked"), function(x) sum(get(x), na.rm = T)), by = c("year", "location_id")]
  setnames(rake_geo, grep("V[0-9]", names(rake_geo), value = T), c(overs, "pop_raked"))
  
  # convert back to rates/prevalence
  rake_geo <- rake_geo[, (overs) := lapply(overs, function(x) get(x) / pop_raked) ]
  
  # merge to admin estimates
  rake_geo <- merge(rake_geo, gbd, by.x = c("location_id", "year"), by.y = c("name", "year"), all.x = TRUE)
  
  # finding the mean of the draws at the raking geographies
  rake_geo$gbd_prev <- rake_geo$mean
  rake_geo[, mbg_rate:= rowMeans(.SD), .SDcols = paste0('V', c(1:ndraws))]
  rake_geo$rf <- rake_geo$gbd_prev / rake_geo$mbg_rate
  
  # Clear Out
  message("creating fractional raking factors table")
  fractional_rf <- rake_geo
  fractional_rf$mbg_prev <- fractional_rf$mbg_rate
  fractional_rf <- fractional_rf[, c("location_id", "year", "mbg_prev", "gbd_prev", "rf")]
  
  
  if(rake_method == "logit"){
    cell_pred <- cell_pred[, (overs) := lapply(overs, function(x) get(x) / pop_raked)]
    #for each country year, find the logit raking factor
    #redefine cys
    cys <- unique(rake_geo[,.(location_id, year)])
    rf <- lapply(seq(nrow(cys)), function(x) {
      
      #year and location
      theloc = cys[x,location_id]
      theyear = cys[x,year]
      
      if (nrow(rake_geo[location_id == theloc & year == theyear]) == 0) {
        return(0)
      } else if (rake_geo[location_id == theloc & year == theyear, .(mean)] == 0 & zero_heuristic == T) {
        # catch true zeroes (i.e. pre-introduction of an intervention) and return -9999. This will speed things up & will replace later with 0
        return(-9999)        
      } else {
        ret <- try(LogitFindK(gbdval     = rake_geo[location_id == theloc & year == theyear,mean],
                              pixelval   = as.matrix(cell_pred[location_id == theloc & year == theyear & !is.na(V1), paste0('V', 1:ndraws), with = F]), #pass the cell pred rows that corrospond to this country year
                              weightval  = as.numeric(cell_pred[location_id == theloc & year == theyear & !is.na(V1), pop_raked]),
                              MaxJump    = MaxJump,
                              MaxIter    = MaxIter,
                              FunTol     = FunTol,
                              approx_0_1 = approx_0_1))
        
        
        if(iterate) {
          # Iterate over larger values of MaxJump and NumIter if needed
          if (is.null(ret) | "try-error" %in% class(ret)) {
            message(paste0("Your GBD and MBG estimates are quite far apart for location ", theloc, " | year ", theyear))
            message("Increasing MaxJump and NumIter by 10-fold, but you should investigate this...")
            
            ret <- try(LogitFindK(gbdval     = rake_geo[location_id == theloc & year == theyear,mean],
                                  pixelval   = as.matrix(cell_pred[location_id == theloc & year == theyear & !is.na(V1), paste0('V', 1:ndraws), with = F]), #pass the cell pred rows that corrospond to this country year
                                  weightval  = as.numeric(cell_pred[location_id == theloc & year == theyear & !is.na(V1), pop_raked]),
                                  MaxJump    = MaxJump*10,
                                  MaxIter    = MaxIter*10,
                                  FunTol     = FunTol,
                                  approx_0_1 = approx_0_1))
          }
          
          # If we get this far, the estimates are generally *really* far apart
          if (is.null(ret) | "try-error" %in% class(ret)) {
            message(paste0("Your GBD and MBG estimates are REALLY far apart for location ", theloc, " | year ", theyear))
            message("Increasing MaxJump and NumIter by 100-fold, but you should investigate this...")
            
            ret <- try(LogitFindK(gbdval     = rake_geo[location_id == theloc & year == theyear,mean],
                                  pixelval   = as.matrix(cell_pred[location_id == theloc & year == theyear & !is.na(V1), paste0('V', 1:ndraws), with = F]), #pass the cell pred rows that corrospond to this country year
                                  weightval  = as.numeric(cell_pred[location_id == theloc & year == theyear & !is.na(V1), pop_raked]),
                                  MaxJump    = MaxJump*100,
                                  MaxIter    = MaxIter*100,
                                  FunTol     = FunTol,
                                  approx_0_1 = approx_0_1))
          }
          
          # Throw error if previous two didn't work
          if (is.null(ret) | "try-error" %in% class(ret)) {
            stop(paste0("Your GBD and MBG estimates are WAY TOO far apart for location ", theloc, " | year ", theyear, " - stopping."))
          }
        }
        
        return(ret)
      }
    })
    rf <- unlist(rf)
    fractional_rf$rf <- rf
  }
  
  # Don't rake if countries_not_to_rake is provided
  if (!is.null(countries_not_to_rake)) {
    
    # Get the GBD location IDs from the ISO codes to set raking
    # factors to 1
    nonrake_table <- get_gbd_locs(
      reg = countries_not_to_rake,
      rake_subnational = FALSE, shapefile_version = shapefile_version
    )
    if (nrow(nonrake_table) > 0) {
      fractional_rf[location_id %in% nonrake_table$location_id, rf := 1]
    }
    rm(nonrake_table)
  }
  
  # Don't subnational rake if countries_not_to_subnat_rake is provided
  if (!is.null(countries_not_to_subnat_rake)) {
    
    # Get the GBD location IDs from the ISO codes to set raking
    # factors to 1. 
    # The main difference here from the above case is that we are pulling
    # in the subnational location IDs if they exist, and also filtering for the
    # case where rak_level (raking level) is >=1 (anything more detailed than ADM0)
    nonrake_table <- get_gbd_locs(
      reg = countries_not_to_subnat_rake,
      rake_subnational = TRUE,
      shapefile_version = shapefile_version
    )[rak_level >= 1]
    if (nrow(nonrake_table) > 0) {
      fractional_rf[location_id %in% get_gbd_locs(
        reg = countries_not_to_subnat_rake,
        rake_subnational = TRUE, shapefile_version = shapefile_version
      )$location_id, rf := 1]
    }
    rm(nonrake_table)
  }
  
  # Fill NAs with 1 BUT WITH A SOLEMN WARNING
  if(nrow(fractional_rf[is.na(rf)]) >= 1 ) {
    warning("You have NAs in your output raking factor. Please check your GBD outputs to make sure you didn't miss any locations.")
    warning("Forcing those missing raking factors to 1")
    fractional_rf[is.na(rf), rf:= 1]
  }
  
  #set filename add-in if using etiology
  if (!(e == 'none')){
    e_string <- paste0('_',e,'_')
  } else (e_string <- '_')
  
  # saving the raking factors
  if (!is.null(custom_output_folder)) {
    write.csv(fractional_rf, file = paste0(custom_output_folder, "/", indicator, e_string, stratum, "_", m, "_rf_", h, ".csv"))
  } else {
    write.csv(fractional_rf, file = paste0(sharedir, "/output/", run_date, "/", indicator, e_string, stratum, "_", m, "_rf_", h, ".csv"))
  }
  
  return(fractional_rf)
}



