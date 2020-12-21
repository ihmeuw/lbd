#################### FRACTIONAL RAKING #############################
#
# In using 5x5 pixels, some pixels are split between several admin
# units. Fractional raking is a way of raking counts which takes
# this into consideration. For any pixel that is split between
# admin units, the fractional area of that admin unit within the
# pixel is multiplied by the count so that counts are properly
# attributed to their respective admin units.
#
# The key to fractional raking is a [link table] which has a row
# for each pixel-admin unit, and includes the gaul codes for admins
# 0-2, as well as the fractional area of the admin unit in that pixel
# (a pixel split between two admin 2 units would have two rows in the
# link table). There is a global link table generated from the stage 1
# and stage 2 simple raster, where every pixel gets a unique id.
# Practically, when using a link table, the global link table and raster
# are subsetted based on the simple raster of the region being raked
# using get_link_table()
#
# Build_link_polygon() and build_link_table() are functions used to
# create the link table, which will only need to be used when
# there is a change in the shapefile used to generate simple rasters.
#
# Many of the functions require the parameter overs, which is a vector
# of the column names of the draws, typically named V1 to Vndraws
# when the cell pred is converted to a data.frame or data.table.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fractional Raking Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @title Fractionally Rake Counts
#' 
#' @description A function used to rake mbg count ouputs to GBD estimates using a fractional methodolgy. Counts in each pixel are assigned to admin units and split by area fraction for pixels with multiple admin units. Counts are then summed to ADM0 (except for those specified as subnational) and raking factors for those admin units are calculated. From there, the raking factors are applied to the linked cell pred, which is then processed in 2 ways: aggregated up to adm0/1/2, and deduplicated (the linked cell pred has 1 row for each pixel-adm unit-year which is reduced to 1 row per each pixel-year). N.B. If you rake subnationally and are missing raking targets for one or more subnationals, the final sum of counts will not match GBD at ADM0.
#' 
#' @param count_cell_pred Cell pred matrix object with counts to be raked.
#' @param rake_to Df with `name`, `year`, and `mean` columns. Values in name must be gbd loc_ids.
#' @param reg character string - Region used to produce cell pred object.
#' @param year_list integer vector - vector of years
#' @param rake_subnational  If true, rakes to ADM1 units not in `countries_not_to_subnat_rake`
#' @param countries_not_to_subnat_rake default set in config, character vector of subnational country iso3 codes NOT to rake to
#' @param countries_not_to_rake default c("GUF", "ESH), french guiana and western sahara. Any country iso3s passed here will have raking factor set to 1. 
#' @param simple_raster default `NULL`, option to pass in simple raster to function if it's been loaded already. 
#' @param modeling_shapefile_version default `current`, string identifying version of of shapefile used in modeling
#' @param raking_shapefile_version default `current`, string identifying version of of shapefile to use in raking
#' 
#'
#' @return Returns a named list with a raked cell pred object, raking factors, and adm0, 1, and 2 aggregated draws raked and unraked
#' 
fractionally_rake_counts <- function(count_cell_pred,
                                     rake_to,
                                     reg,
                                     year_list,
                                     rake_subnational,
                                     countries_not_to_subnat_rake,
                                     countries_not_to_rake = c("ESH+GUF"),
                                     simple_raster = NULL,
                                     modeling_shapefile_version = 'current',
                                     raking_shapefile_version = 'current') {
  #get draw column names
  ndraws <- ncol(count_cell_pred)
  overs <- paste0('V',1:ndraws)
  
  ## load objects used in modeling
  if(is.null(simple_raster)) {
    ## get simple polygon and simple raster used to produce cell pred
    message('Loading Simple Polygon')
    simple_polygon <- load_simple_polygon(gaul_list = get_adm0_codes(reg,
                                                                     shapefile_version = modeling_shapefile_version),
                                          buffer = 0.4,
                                          shapefile_version = modeling_shapefile_version)
    subset_shape   <- simple_polygon[['subset_shape']]
    simple_polygon <- simple_polygon[['spoly_spdf']]
    
    message('Loading Simple Raster\n')
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    
  } else {
    message('Using Supplied Simple Raster')
  }
  
  message("Loading Link Table and Prepping Cell Pred for Raking")
  
  #get link table
  link_table_output <- get_link_table(simple_raster, raking_shapefile_version)
  
  #format get_link_table outputs and add pixel ids to cell pred
  link <- link_table_output[["link_table"]]
  pixel_id <- link_table_output[["pixel_ids"]]
  pixel_id <- rep(pixel_id, times=length(year_list))
  year <- rep(year_list, each=(nrow(count_cell_pred)/length(year_list)))
  
  if(length(pixel_id) != nrow(count_cell_pred)){
    message("Something has gone wrong with the matchup of pixels between the cell pred and link table-")
    message("Length of pixel vector: ", length(pixel_id))
    message("Rows in cell pred: ", nrow(count_cell_pred))
    message("This may be due to a mismatch between the cell pred and simple raster, or pixels in the simple raster not being in the link table.")
    stop()
  }
  
  id_count_cell <- cbind(count_cell_pred, pixel_id, year)
  
  # subset link table to countries in region -- fix for border issue with fractional raking
  link <- link[ADM0_CODE %in% unique(simple_raster),]
  
  # build connector object between ADM codes and GBD loc ids
  connector <- data.table(get_gbd_locs(rake_subnational = rake_subnational,
                                       reg = reg,
                                       shapefile_version = raking_shapefile_version))
  
  if(rake_subnational){
    connector[, ADM_CODE := ADM0_CODE]
    connector[ADM1_CODE != 0, ADM_CODE := ADM1_CODE]
  }
  connector <- connector[,c("ADM_CODE", "location_id")]
  
  # add on ADM_CODE to raking factors so everything can be done in GADM codes rather than GBD loc ids
  rake_to <- merge(rake_to, connector, by.x = "name", by.y = "location_id", all.x = T)
  
  # set ADM_CODE which is used as the adm_level for raking
  # when subnational countries are provided, sets ADM_CODE to the ADM1_CODE
  link[, ADM_CODE := ADM0_CODE]
  if(rake_subnational){
    
    # Get the GBD location IDs from the ISO codes to set raking
    # factors to 1. We are pulling
    # in the subnational location IDs if they exist, and also filtering for the
    # case where rak_level (raking level) is >=1 (anything more detailed than ADM0)
    subnational_country_code <- get_gbd_locs(
      reg = countries_not_to_subnat_rake,
      rake_subnational = TRUE,
      shapefile_version = raking_shapefile_version
    )[rak_level >= 1, ADM0_CODE]
    link[ADM0_CODE %in% subnational_country_code, ADM_CODE := ADM1_CODE]
  }
  
  # get raking factors 
  message("\nCalculating Raking Factors")
  fractional_rf <- get_fractional_rf(id_count_cell, link, rake_to, reg, overs, raking_shapefile_version)
  
  if(nrow(fractional_rf[is.na(rf)]) > 0) {
    message("Warning: there are some admin units that are missing raking factors. This is likely due to a lack of raking targets. If these are subnational locations, your counts will not add up at the ADM0 level.")
    message("These locations (GADM code) are missing raking factors: ", paste(unique(fractional_rf[is.na(rf), ADM_CODE]), collapse = ", "))
  }
  
  # getting unraked draws - merge cell pred and link table
  unraked_linked_cell_pred <- merge(link, id_count_cell, by.x = "ID", by.y = "pixel_id", all.y = T, allow.cartesian=T)
  # multiply counts by area fraction
  unraked_linked_cell_pred[, (overs) := .SD * as.numeric(area_fraction), .SDcols = overs]
  # aggregate unraked draws to adm0/1/2
  unraked_aggregate_count_list <- aggregate_counts(unraked_linked_cell_pred, overs, raked = F)
  
  # setting rf for non-raking countries
  fractional_rf[ADM0_CODE %in% get_adm0_codes(countries_not_to_rake, shapefile_version = raking_shapefile_version), rf := 1]
  
  # apply raking factors
  message("\nApplying Raking Factors")
  raked_counts_linked_cell_pred <- apply_fractional_rf(id_count_cell, link, fractional_rf, overs)
  
  # deduplicate linked cell pred that has multiple rows per pixel
  message("Aggregating...\n")
  deduped_counts_cell_pred <- dedupe_linked_cell_pred(raked_counts_linked_cell_pred, overs)
  
  # aggregate deaths to Adm 0-2
  raked_aggregate_count_list <- aggregate_counts(raked_counts_linked_cell_pred, overs, raked = T)
  
  # compile outputs
  output_list <- list()
  output_list[["raked_cell_pred"]] <- deduped_counts_cell_pred
  output_list <- c(output_list, unraked_aggregate_count_list, raked_aggregate_count_list)
  output_list[["raking_factors"]] <- fractional_rf
  
  return(output_list)
}

#' @title Save outputs from fractional counts raking
#'
#' @description Save outputs from fractional count raking
#' coming out of \code{fractionally_rake_counts} call.
#'
#' @param output_list The list of outputs from \code{fractionally_rake_counts} call.
#' @param sharedir Output to directory
#' @param indicator Indicator
#' @param age Age
#' @param reg Region
#'
#' @param holdout Holdout
#'
#' @examples
#' \dontrun{
#' raked_frax_counts_save(
#'   output_list = outputs,
#'   sharedir = '<<<< FILEPATH REDACTED >>>>',
#'   indicator = "tr_had_diarrhea",
#'   age = 0,
#'   reg = "tza",
#'   holdout = 0
#' )
#' }
#' 
#' @rdname raked_frax_counts_save
#' @export
raked_frax_counts_save <- function(output_list, sharedir, indicator, age, reg, holdout) {
  
  ## Save raking factors
  fwrite(output_list[["raking_factors"]],
         file = paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_", reg, "_rf.csv")
  )
  
  # Save raked rates cell_pred object (which is equal to counts values!!)
  
  #### Need to unpack it, otherwise R gets sad ####
  raked_cell_pred <- output_list[["raked_cell_pred"]]
  
  save(raked_cell_pred,
       file = paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_raked_cell_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
  )
  
  ## Save raked counts aggregations
  raked_adm0_draws <- output_list[["raked_adm0_draws"]]
  raked_adm1_draws <- output_list[["raked_adm1_draws"]]
  raked_adm2_draws <- output_list[["raked_adm2_draws"]]
  save(raked_adm0_draws, raked_adm1_draws, raked_adm2_draws,
       file = paste0(main_dir, indicator, "_", "raked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
  )
  
  ## save unraked counts aggregations
  unraked_adm0_draws <- output_list[["unraked_adm0_draws"]]
  unraked_adm1_draws <- output_list[["unraked_adm1_draws"]]
  unraked_adm2_draws <- output_list[["unraked_adm2_draws"]]
  
  save(unraked_adm0_draws, unraked_adm1_draws, unraked_adm2_draws,
       file = paste0(main_dir, indicator, "_", "unraked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
  )
  
  
  return(0)
}





#' Calculate Fractional Raking Factors
#'
#' @Description: Calculate raking factors considering cases where a pixel is in
#' multiple admin units for deaths or cases
#'
#' @param cell_pred a cell pred with deaths/cases in data.table format with a pixel_id column from `get_link_table()`
#' @param link link table gotten from `get_link_table()`
#' @param gbd data.table with a columns - `gaul`, `year_id`, and `gbd_field (see below)`
#' with estimates for the indicator being raked to for each year-admin0
#' @param region the region used to make the simple raster
#' @param overs a vector of the column names in the cell pred corresponding to draws
#' (typically V1:Vndraws)
#' @param shapefile_version string indicating which shapefile version to match ADM0 codes against
#'
#' @return a data.table with raking factor by adm0 gaul code and year
#'
#' @import data.table
#' @import dplyr
#' 
get_fractional_rf <- function(cell_pred, link, gbd, region, overs, shapefile_version) {
  
  #merge link table and cell pred on pixel id
  linked_cell_pred <- merge(link, cell_pred, by.x = "ID", by.y = "pixel_id", all.y = T, allow.cartesian=T)
  
  if(nrow(linked_cell_pred[is.na(ADM_CODE)]) > 0){
    message("There were ", nrow(linked_cell_pred[is.na(ADM_CODE)]), " rows in the dataset that did not have a corresponding pixel in the link table" )
    message("The following pixels (pixel_id from link table) were affected: ")
    message(paste(linked_cell_pred[is.na(ADM_CODE)]$pixel_id, collapse = ", "))
  }
  
  #multiply death draws by area fractions and aggregate draws to ADM0 and year
  linked_cell_pred <- linked_cell_pred[, lapply(overs, function(x) sum(get(x) * as.numeric(area_fraction), na.rm = T)), by = c('year','ADM_CODE', 'ADM0_CODE')]
  
  #get mean of draws
  linked_cell_pred[,mbg_rf := rowMeans(linked_cell_pred[,overs, with=F])]
  
  #merge with gbd and calculate population raking factor
  rf_df <- merge(linked_cell_pred, gbd, by= c("ADM_CODE", "year"), all.x=T)
  rf_df <- rf_df[, rf := mean / mbg_rf]
  
  #drop unecessary columns
  rf <- rf_df[,c("ADM_CODE", "year", "rf", "ADM0_CODE")]
  
  #set raking factors for countries not in region to 1 to avoid overraking in countries
  #partially included in the simple raster. Not doing this causes the values of border
  #pixels to be much higher than they should be.
  rf[!(ADM0_CODE %in% get_adm0_codes(region, shapefile_version = shapefile_version)), rf := 1]
  return(rf)
}


#' Apply Fractional Raking Factors
#'
#' @Description: Apply raking factors considering cases where a pixel is in
#' multiple admin units for deaths or cases
#'
#' @param cell_pred a cell pred with deaths/cases
#' @param link link table generated by `build_link_table()`
#' @param fractional_rf a data.table generated by `get_fractional_rf`
#' @param overs a vector of the column names in the cell pred corresponding to draws
#' (typically V1:Vndraws)
#'
#' @return a linked cell pred with the rf applied
#'
#' @import data.table
#'
apply_fractional_rf <- function(cell_pred, link, fractional_rf, overs) {
  #generate id for rows in cell pred (merges below mess up ordering)
  cell_pred[, cell_pred_id := .I]
  #merge cell pred and link file
  linked_cell_pred <- merge(link, cell_pred, by.x = "ID", by.y = "pixel_id", all.y = T, allow.cartesian=T)
  #merge linked cell pred with raking factors
  linked_cell_pred <- merge(linked_cell_pred, fractional_rf, by.x = c("ADM_CODE", "year", "ADM0_CODE"), by.y = c("ADM_CODE", "year", "ADM0_CODE"), all.x =T)
  #converting from "units" type to double, truncates decimal otherwise
  linked_cell_pred$area_fraction <- as.double(linked_cell_pred$area_fraction)
  #calculate deaths on all draws
  raked_cell_pred <- linked_cell_pred[, (overs) := lapply(overs, function(x) get(x)* rf * area_fraction)]
  
  #sort on cell_pred_id
  setorder(raked_cell_pred, cell_pred_id)
  
  return(raked_cell_pred)
}


#' Deduplicate linked cell preds
#'
#' @Description: Collapse rows of cell preds with multiple admin units in one pixel back into
#' single rows
#'
#' @param linked_cell_pred a cell pred merged to a link table generated by `get_link_table()`
#' @param overs a vector of the column names in the cell pred corresponding to draws
#' (typically V1:Vndraws)
#'
#' @return a cell pred object in original matrix form which corresponds to simple raster
#'
#' @import data.table
#'
dedupe_linked_cell_pred <- function(linked_cell_pred, overs) {
  
  #sum deaths for each duplicated pixel
  dedup_dt <- linked_cell_pred[, lapply(.SD, sum), by=cell_pred_id, .SDcols=overs]
  
  #reorder on cell_pred_id to match original order
  dedup_dt$cell_pred_id <- as.numeric(as.character(dedup_dt$cell_pred_id))
  setorder(dedup_dt, cell_pred_id)
  
  #drop columns and convert back into matrix
  deduped_linked_cell_pred <- dedup_dt[,overs, with = FALSE]
  deduped_linked_cell_pred <- as.matrix(deduped_linked_cell_pred)
  colnames(deduped_linked_cell_pred) <- NULL
  return(deduped_linked_cell_pred)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Post-Raking Processing
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Aggregate raked linked cell pred
#'
#' @Description: rake linked cell pred to admins0-2, maintaining draws
#'
#' @param linked_cell_pred a raked and linked cell pred object
#' @param year_list list of years in cell_pred
#' @param raked Boolean, raked or unraked cell draws
#'
#' @return a list with 3 data.tables for admins0-2 with draws intact
#'
aggregate_counts <- function(linked_cell_pred, overs, raked) {
  admin_2 <- linked_cell_pred[, lapply(overs, function(x) sum(get(x), na.rm = T)), by = c('year', 'ADM2_CODE', 'ADM0_CODE')]
  admin_1 = linked_cell_pred[, lapply(overs, function(x) sum(get(x), na.rm = T)), by = c('year','ADM1_CODE', 'ADM0_CODE')]
  admin_0 = linked_cell_pred[, lapply(overs, function(x) sum(get(x), na.rm = T)), by = c('year','ADM0_CODE')]
  if (raked) {
    return(list("raked_adm0_draws" = admin_0, "raked_adm1_draws" = admin_1,"raked_adm2_draws" = admin_2))
  } else {
    return(list("unraked_adm0_draws" = admin_0, "unraked_adm1_draws" = admin_1,"unraked_adm2_draws" = admin_2))
  }
}


#' Convert cell pred into raster brick
#'
#' @Description: Insert cell pred object into a series of raster layers for each year
#' then combine into rasterBrick
#'
#' @param cell_pred a cell pred object in matrix format
#' @param simple_raster the simple raster corresponding to the cell_pred_object
#' @param year_list list of years in cell_pred
#' @param func function to summarize cell pred. `mean`, `upper`, and `lower` are acceptable
#'
#' @return a rasterBrick with a layer for each year in the cell_pred
#'
#' @import matrixStats
#'
cell_pred_to_raster_brick <- function(cell_pred, simple_raster, year_list, func) {
  require(matrixStats)
  
  raster_list <- vector("list", length = length(year_list))
  
  #apply summary function
  if(func == "mean") {
    message("summarizing draws to mean")
    death_vec <- rowMeans(cell_pred)
  } else if (func == "upper") {
    message("summarizing draws to upper quantile - this can be slow for large cell_preds")
    death_vec <- matrixStats::rowQuantiles(cell_pred, probs = 97.5 / 100)
  } else if (func == "lower") {
    message("summarizing draws to lower quantile - this can be slow for large cell_preds")
    death_vec <- matrixStats::rowQuantiles(cell_pred, probs = 2.5 / 100)
  } else {
    stop("Not a valid function, please choose 'mean', 'upper', or 'lower'.")
  }
  
  #split years into separate vectors in a list
  death_vec_list <- split(death_vec, cut(seq_along(death_vec), length(year_list), labels = FALSE))
  #get raster pixels with data
  pixel_id <- which(!is.na(getValues(simple_raster)))
  message("Building raster brick \n")
  #loop through year vectors and insert into rasters
  for(i in 1:length(death_vec_list)){
    death_ras <- simple_raster
    raster_list[[i]] <- insertRaster(death_ras, cbind(death_vec_list[[i]]))
  }
  #set years as raster layer names
  names(raster_list) <- year_list
  
  #convert to rasterBrick
  death_raster <- brick(raster_list)
  return(death_raster)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Link Table Functions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Subset global link table to regional table
#'
#' @description: Reads in the global link table and id raster and subsets using
#' a regional simple raster.
#'
#' @param simple_raster Simple raster used for modelling region
#'
#' @return list of 2 objects:
#'    - a link table
#'    - a vector of pixel ids used to link one year in a cell pred to the link table.
#'
#' @import data.table
#' @import raster
#'
get_link_table <- function(simple_raster, shapefile_version) {
  #read in link_table
  global_link <- data.table(readRDS(sprintf("<<<< FILEPATH REDACTED >>>>/lbd_standard_link.rds", shapefile_version)))
  #read in id_raster
  global_raster <- readRDS(sprintf("<<<< FILEPATH REDACTED >>>>/lbd_standard_id_raster.rds", shapefile_version))
  #crop id_raster to simple_raster
  cropped_raster <- raster::crop(global_raster, simple_raster)
  #mask id_raster with simple_raster
  masked_raster <- raster::mask(cropped_raster, simple_raster)
  #extract id_raster
  pixel_ids <- raster::extract(masked_raster, extent(masked_raster), na.rm = T)
  pixel_ids <- pixel_ids[!is.na(pixel_ids)]
  #subset link table by pixel ids from extract
  link_table <- global_link[ID %in% pixel_ids,]
  #make sure list of pixel ids is same number of nonna cells in simple raster
  if(length(pixel_ids) != length(which(!is.na(getValues(simple_raster))))) {
    stop("The number of pixel_ids does not match the number of non-na values in simple raster. \nPart of the simple raster may be outside the extent of the global id raster.")
  }
  #return subsetted link table and list of pixel ids
  return(list("link_table" = link_table, "pixel_ids" = pixel_ids))
}


#' @title Link Table Generation
#'
#' @description Builds a link table where each row represents an admin2-pixel, and includes
#' the admin0-2 levels, the id of the pixel (based on simple raster), and the percent of the
#' pixel covered by a country. Used for fractional raking where a pixel is in multiple
#' administrative units.
#'
#' @param shapefile_version version of admin file to base link table on
#' @param cores number of cores to run foreach over
#' @param region the iso3 region(s) used to generate the simple raster
#' @param sample_ids whether or not to sample pixel ids from the stage1+stage2 raster
#' @param custom_shapefile_path A non-admin shapefile to build a link table for
#' @param custom_shapefile_field the field in custom_shapefile with admin identifiers
#' 
#' @return returns the link data.table described above, and the id raster and polygon used in the calculation
#' @export

build_link_table <- function(shapefile_version,
                              cores,
                              region = "stage1+stage2",
                              sample_ids = FALSE,
                              custom_shapefile_path = NULL,
                              custom_shapefile_field = 'ADM0_CODE') {
  
  # load shapefile, subset
  if (is.null(custom_shapefile_path)) {
    polys <- sf::st_read(get_admin_shapefile(admin_level = 2, version = shapefile_version), stringsAsFactors = FALSE)
    codes <- get_adm0_codes(region, shapefile_version = shapefile_version)
    polys <- polys[polys$ADM0_CODE %in% codes,]
  } else {
    polys <- sf::st_read(custom_shapefile_path, stringsAsFactors = FALSE, quiet = TRUE)
  }
  
  # get simple_raster (replaces build_simple_raster_pop)
  master_pop <- raster::brick('<<<< FILEPATH REDACTED >>>>/WorldPop_total_global_stack.tif')
  cropped_pop <- raster::crop(master_pop, polys, snap = 'out')
  simple_raster <- raster::rasterize(polys, cropped_pop, field = custom_shapefile_field)
  
  # resample from stage1+stage2 raster to get proper pixel ids
  if (sample_ids) {
    world_raster <- empty_world_raster()
    world_raster[] <- 1:length(world_raster)
    simple_raster <- raster::resample(world_raster, simple_raster, method = 'ngb')
  } else {
    simple_raster[] <- 1:length(simple_raster)
  }
  
  # pixel polygons (replaces build_link_polygon)
  pixels <- spex::polygonize(simple_raster)
  
  # change 'layer' to 'pixel_id'
  setnames(pixels, 'layer', 'pixel_id')

  # assign crs
  sf::st_crs(pixels) <- 4326
  
  # intersection function to be passed into foreach, once for each record in polys
  get_intersect <- function(poly) {
    # validate geometry
    if (!sf::st_is_valid(poly)) {
      poly <- lwgeom::st_make_valid(poly)
    }
    
    # crop the baseline raster, get values
    px_crop <- raster::values(raster::crop(simple_raster, poly, snap = 'out'))

    # subset possible pixels
    target <- pixels[pixels$pixel_id %in% px_crop,]
    
    # get a list of target polys (e.g. zip code) that intersect with the poly (e.g. critical habitat)
    intersects <- sf::st_intersects(poly, target) # returns sparse matrix
    
    # compute start area of the poly
    target <- target[intersects[[1]],]  # sparse matrix index
    target$start_area <- sf::st_area(target)
    
    # compute the intersection
    intersection <- sf::st_intersection(poly, target)
    
    # no intersection
    if (nrow(intersection) == 0) {
      sf::st_geometry(poly) <- NULL
      return(poly)
    }
    
    # compute end area and area fraction
    intersection$end_area <- sf::st_area(intersection)
    intersection$area_fraction <- intersection$end_area / intersection$start_area
    
    # remove the geometry to save space
    intersection$geometry <- NULL
    
    return(intersection)
  }
  
  # go through id polygons in parallel and compute intersections with ADM2 polygon
  registerDoParallel(cores = cores)
  link_table <- foreach(i = 1:nrow(polys), .final = rbindlist) %dopar% get_intersect(polys[i,])
  
  # create ID variable (copy of pixel_id) for backwards compatibility
  link_table$ID <- link_table$pixel_id
  
  # converting from "unit" to numeric
  link_table$area_fraction <- as.numeric(link_table$area_fraction)
  
  # fixes area_fractions for pixels that are partially in water
  link_table <- fix_link(link_table)
  
  return(list("link_table" = link_table, "id_raster" = simple_raster, "id_poly" = pixels))
}

#' Build Link Polygon
#'
#' @description: Builds a polygon where there is a layer for each pixel in a raster,
#' with an ID field referring to the id of the pixel in the raster
#'
#' @param region character string of region to be pulled from `get_adm0_codes()`
#' @param simple_raster default `NULL`. provides option to pass in simple raster rather
#' than making it using a region name to save time.
#' @param shapefile_version string to identify version of shapefile to build simple_raster
#'
#' @return returns an sf polygon
#'
#' @example idpoly <- build_link_polygon(region)
#'
#' @import raster
#' @import sp
#' @import rgeos
#' @import spex
#' @import sf
#'
build_link_polygon <- function(region, simple_raster = NULL, shapefile_version = 'current') {
  require(raster)
  require(sp)
  require(rgeos)
  require(spex)
  require(sf)
  
  #load simple raster (pass in as parameter for speed up)
  if (is.null(simple_raster)) {
    simple_polygon <- load_simple_polygon(gaul_list = get_adm0_codes(region, shapefile_version = shapefile_version),
                                          buffer = 0.4,
                                          shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon[['subset_shape']]
    simple_polygon <- simple_polygon[['spoly_spdf']]
    
    message('Loading simple raster')
    raster_list    <- build_simple_raster_pop(subset_shape) #,u5m=TRUE)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
  }
  
  #make a copy of simple raster and set all values to 0
  new_ras <- simple_raster
  values(new_ras) <- 0
  #give each pixel a unique id
  new_ras[1:length(simple_raster)] <- 1:length(simple_raster)
  
  #convert raster to sf polygon and give it an ID field
  message("converting ID raster to polygon")
  idpoly <- spex::polygonize(new_ras)
  idpoly$ID <- 1:nrow(idpoly)
  
  return(idpoly)
}


#' Fix area fractions in link table
#'
#' @description: Scales area fractions for pixels with water. Pixels with 1 admin unit and
#' water have the area fraction set to 1, as all values will occur on the land. Pixels
#' containing multiple admin units and water are scaled proportionally. Internal function
#' used in `build_link_table()`
#'
#' @param link link table
#'
#' @return a link table with adjusted area fractions
#'
fix_link <- function(link) {
  #figure out which pixels have multiple rows
  link_fix <- link %>%
    group_by(ID) %>%
    dplyr::summarise(total_area = sum(area_fraction),
              n = n()) %>%
    data.table()
  
  #merge link, link_fix and cast some variables to numeric
  link_fix2 <- merge(link, link_fix, by = "ID")
  link_fix2[, area_fraction := as.numeric(area_fraction)]
  link_fix2[, total_area := as.numeric(total_area)]
  #set area to 1 for all pixels in a single admin unit and water
  link_fix2[n == 1 & area_fraction < 1, area_fraction := 1]
  #recalculate area for pixels in multiple admin units and water proportionally
  link_fix2[n > 1 & total_area < 1, area_fraction := area_fraction / total_area]
  return(link_fix2)
}


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
#' @param sharedir sharedir
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
                                  meas = 'prevalence'){
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
                                            m = meas)
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
#' @param sharedir Path to shared directory
#' @param run_date model run date string
#' @param indicator modeled indicator
#' @param main_dir Path to main directory
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
                                 meas = 'prevalence') {
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
  
  ## load the raking factors
  fractional_rf <- read.csv(file = ifelse(is.null(custom_output_folder),
                                          paste0(sharedir, "<<<< FILEPATH REDACTED >>>>", indicator, "_", reg, "_", meas, "_rf_", holdout, ".csv"),
                                          paste0(custom_output_folder, "/", indicator, "_", reg, "_", meas, "_rf_", holdout, ".csv")
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
                                        sharedir, "<<<< FILEPATH REDACTED >>>>",
                                        indicator, "_cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"                                      
                                        ),
                                      paste0(
                                        custom_output_folder, "/",
                                        indicator,  "_cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"                                      
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
                                          sharedir, "<<<< FILEPATH REDACTED >>>>",
                                          indicator, "_c_cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"                                        
                                        ), paste0(
                                          custom_output_folder, "/",
                                          indicator, "_c_cell_draws_raked_", meas, "_eb_bin", age, "_", reg, "_", holdout, ".RData"
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
                     paste0(main_dir, "/", indicator, "_", "raked", "_", meas, "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, "_", "raked", "_", meas, "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
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
                     paste0(main_dir, "/", indicator, "_", "raked", "_", meas, "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, "_", "raked", "_", meas, "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
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
                     paste0(main_dir, "/", indicator, "_", "unraked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, "_", "unraked", "_c", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
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
                     paste0(main_dir, indicator, "_", "unraked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData"),
                     paste0(custom_output_folder, "/", indicator, "_", "unraked", "_admin_draws_eb_bin", age, "_", reg, "_", holdout, ".RData")
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




# Functions called internally to the above two functions

#' @title fetch_from_rdata
#' @description: this pulls out pieces of an r data object
#'               
#' @param file_location where to look for the Rdata object
#' @param item_name the item to pull from the Rdata object
#' @param use_grep Weather or not to use grep to look for similar file names
#'
#' @return loads the specific item from an Rdata object
#' @export

fetch_from_rdata <- function(file_location, item_name, use_grep = F) {
  load(file_location)
  if (use_grep) {
    return(mget(grep(item_name, ls(), value = T)))
  }else{
    return(get(item_name))
  }
}



#' @title load_populations_cov
#' @description: This gets the population in each cell which you need to convert the rates to counts for fractional raking
#' it can also get the stackers if you would like
#'               
#' @param reg modeling region over which you are operating
#' @param pop_measure the world pop covariate measure used in your model
#' @param measure needs to be depricated, currently a flag that if set to 'prevalence' this function also grabs teh covatiate stakers
#' @param simple_polygon the simple polygon for your region
#' @param simple_raster the simple raster for your model
#' @param year_list The years you are modeling over, will be used to pull populations
#' @param interval_mo the number of months between time steps
#' @param outputdir currently not used, there is a commented out line of code which can save the population raster if you want
#'
#' @return used internally above, returns a data frame with the population in each cell-year so rates can be turned into counts
#' @return for fractional raking and aggregation
#' @export


load_populations_cov <- function(reg,
                                 pop_measure,
                                 measure = 'prevalence',
                                 simple_polygon,
                                 simple_raster,
                                 year_list,
                                 interval_mo,
                                 outputdir,
                                 pixel_id){
  message("Loading the world pop rasters for this region")
  ## Pull 2000-2015 annual population brick using new covariates function
  pop_raster_annual <- load_worldpop_covariate(template_raster = simple_polygon,
                                               pop_measure = pop_measure,
                                               pop_release = pop_release,
                                               start_year = min(year_list),
                                               end_year = max(year_list),
                                               interval = as.numeric(interval_mo))
  
  ## extend and crop pop raster to ensure it matches the simple raster #not convinced this is totally needed
  pop <- pop_raster_annual[[1]]
  pop  <- extend(pop, simple_raster, values = NA)
  pop  <- crop(pop, extent(simple_raster))
  pop  <- setExtent(pop, simple_raster)
  pop  <- raster::mask(pop, simple_raster)
  
  ## check to ensure the pop raster matches the simple raster in extent and resolution
  if (extent(pop) != extent(simple_raster)) {
    stop("population raster extent does not match simple raster")
  }
  if (any(res(pop) != res(simple_raster))) {
    stop("population raster resolution does not match simple raster")
  }
  
  message("loading the covariate stakers for this model")
  if (measure == 'prevalence') {
    covs = fetch_from_rdata(paste0('<<<< FILEPATH REDACTED >>>>', run_date, pathaddin, '.RData'), 'cov_list')
    fes = fetch_from_rdata(paste0('<<<< FILEPATH REDACTED >>>>', run_date, pathaddin, '.RData'), 'all_fixed_effects')
    
    submodels = trimws(strsplit(fes, '+', fixed = T)[[1]])
    covs = covs[submodels]
    #make sure spatial extent is the same
    covs = lapply(covs, function(x) invlogit(crop(x, simple_raster)))
  }else{
    covs = list()
  }
  
  #bring everything into one place
  covs$pop = crop(pop, simple_raster)
  covnames = names(covs)
  
  #ensure the dimensions are the same
  for (ccc in covs) {
    stopifnot(dim(ccc)[1:2] == dim(simple_raster)[1:2])
  }
  
  message("converting the covariate stackers and pop data in to a data table")
  #convert to datables, reshape and stuff
  brick_to_dt = function(bbb, pixel_id = pixel_id){
    dt = setDT(as.data.frame(bbb))
    dt[, pxid := .I] #probably uncessary
    
    #drop rows now in cellIdx
    dt = dt[pixel_id,]
    
    dt = melt(dt, id.vars = 'pxid', variable.factor = F)
    dt = dt[,.(value)]
    return(dt)
  }
  covdt <- covnames
  for (iii in 1:length(covs)) {
    dt <- brick_to_dt(bbb = covs[[iii]], pixel_id = pixel_id)
    covdt[[iii]] <- dt
  }
  #covdt = lapply(covs, brick_to_dt(bbb = covs, pixel_id = pixel_id))
  covdt = do.call(what = cbind, covdt)
  setnames(covdt, names(covs))
  
  # Add pixel_id, but make sure that its recycled explicitly as per data.table 1.12.2 guidelines
  covdt[, pixel_id := rep(pixel_id, times = nrow(covdt) / length(pixel_id))]
  
  #set pop to 0 when pop is na
  covdt[is.na(pop), pop := 0]
  
  #add year to covdt
  yyy = as.vector(unlist(lapply(min(year_list):max(year_list), function(x) rep.int(x, times = length(pixel_id)))))
  covdt[,year := yyy]
  
  #free up a bit of space
  rm(covs)
  
  return(covdt)
}



#' @title sub_nat_link_merge
#' @description: This merges the link table and the list of GBD geographies so that we know what aggregations to use for raking.
#'               
#' @param rake_subnational T/F are we using the subnational targets for raking.  What targets get used is in the raking shp
#' @param link the link table
#' @param connector the table identifying what geographies go to what raking targets
#' @param nat_connector a connector object from `get_gbd_locs()` with `rake_subnational = F`. Used to replace subnationals specified in `countries_not_to_subnat_rake`
#' @param countries_not_to_subnat_rake a character vector of iso3 codes gotten from the config. All countries in this list will be raked nationally.
#'
#' @return a link table that has not only the spatial hierarchy in it but also a column for raking geographies.
#'
#' @export
sub_nat_link_merge <- function(rake_subnational,
                               link,
                               connector,
                               nat_connector = NULL,
                               countries_not_to_subnat_rake = NULL){
  
  if (rake_subnational == T) {
    connector0 <- connector[rak_level == 0, ]
    connector1 <- connector[rak_level == 1, ]
    
    #move countries from subnat connector to national connector
    if(!is.null(countries_not_to_subnat_rake)){
      if(is.null(nat_connector)) {
        stop("A connector table built with get_gbd_locs() with rake_subnational = F must be passed in when countries_not_to_subnat_rake is not null.")
      } else {
        nat_connector <- nat_connector[ADM_CODE %in% get_adm0_codes(countries_not_to_subnat_rake, shapefile_version = modeling_shapefile_version),]
        nat_connector[, rak_level := 0]
        nat_connector[, ADM1_CODE := 0]
        nat_connector[, ADM0_CODE := ADM_CODE]
        
        connector0 <- rbind(connector0, nat_connector)
        connector1 <- connector1[!(ADM0_CODE %in% get_adm0_codes(countries_not_to_subnat_rake, shapefile_version = modeling_shapefile_version)),]
        
      }
    }
    
    # Remove any duplicated rows
    connector0 <- unique(connector0)
    connector1 <- unique(connector1)
    
    # Ensure that all subnationals are properly assigned to connectors
    # When the raking shapefile version does not have all subnationals, 
    # there are possible scenarios where requested subnationals will not be present.
    subnational_countries <- c("POL+RUS+UKR+NZL+JPN+USA+ITA+NOR+SWE+GBR+MEX+BRA+IRN+IND+PAK+CHN+IDN+PHL+ETH+KEN+ZAF+NGA")
    
    #invert countries_not_to_subnat_rake
    subnat_countries_to_rake <- setdiff(get_adm0_codes(subnational_countries, shapefile_version = modeling_shapefile_version), 
                                        get_adm0_codes(countries_not_to_subnat_rake, shapefile_version = modeling_shapefile_version))
    
    #check to see if any requested subnationals are in the national connector
    if(any(unique(connector0[, ADM0_CODE]) %in% subnat_countries_to_rake)){
      stage_list <- fread("<<<< FILEPATH REDACTED >>>>/stage_master_list.csv")
      stage_list <- stage_list[, c("iso3", "loc_id")]
      con_error <- connector0[connector0[, ADM0_CODE] %in% subnat_countries_to_rake,]
      
      con_error <- merge(con_error, stage_list, by.x = "location_id", by.y = "loc_id")
      
      message("The following countries are not present subnationally in the raking shapefile and will be raked nationally: ", paste(unique(con_error$iso3), collapse = ", "))
      message("All subnationals are included in admin shapefiles from 4/2019 or later")
    }
    
    connector0$ADM1_CODE <- NULL
    connector1$ADM0_CODE <- NULL
    
    link0 <- merge(link, connector0, by = c("ADM0_CODE"))
    link1 <- merge(link, connector1, by = c("ADM1_CODE"))
    link <- rbind(link0, link1)
  } else {
    link <- merge(link, connector, by.x = c("ADM0_CODE"), by.y = c("ADM_CODE"))
  }
  return(link)
}


#' @title prep_link_table
#' @description: This prepares the link table for merging and does the zippering operation on the region boundary.
#' The zippering means that a given cell is only in one modeled reagion.
#'               
#' @param link_table the link table for your region
#' @param simple_raster the simple raster for your region
#' @param pixel_ids a list of all valif pixel ids in your simple raster
#'
#' @return a link table such that cell fractions that fall either in the ocean or in a country outside of your region are dropped
#' and the areas of the remain cell fractions for those cells are rescaled to equal 1 so no data is lost.  The result is cells that
#' lie on the border of the region are zippered into one region or the other.
#' @export
#' 
prep_link_table <- function(link_table = link_table,
                            simple_raster = simple_raster,
                            pixel_id = pixel_id){
  #####################################################################
  # Prepping the cell_pred and link table to be linked
  
  
  #keep only locations in the region
  link <- link_table[[1]]
  link = link[ADM0_CODE %in% unique(simple_raster),]# This means only cells and cell fractions that are in the region are included in the link table
  
  #keep only cells where simple raster identifies it as belonging to the country or region
  link_map = data.table(af_id = link_table[[2]], reg_id = pixel_id, simp_ras_val = simple_raster[pixel_id])#creates table relating pixel number in the africa simple raster (id in link table) to pixel number in simple raster (id in simple raster/cell pred), to the value in the simple raster(ADM0_CODE).
  link = merge(link, link_map, by.x = 'ID', by.y = 'af_id')#merges that relational table to the link table
  
  #scale up such that the sum of area fractions for each pixel is equal to 1 so pieces are not lost
  link[, total_frac := sum(area_fraction), by = "ID"]
  link[, c('area_fraction','old_frac') := list(area_fraction/total_frac, area_fraction)]
  
  return(link)
}


#' @title prep_cell_pred
#' @description: This converts the cell pred to a data.table and sets up the variables needed for merging
#'               
#' @param cell_pred the cell pred object that will be raked and aggregated
#' @param cell_ids the cell ids from the global simple raster, used to connect the data to the link table
#' @param pixel_ids the cell ids from the simple raster
#' @param covdt the data table containing the populations and possibly the staker covariates as well
#'
#' @return returns the cell_pred as a dt with additional columns for cell_id and pixel_id wich allows for
#' merging to the link table.  Also adds the population and covariate columns if you would like.
#'
#' @export
prep_cell_pred <- function(cell_pred=cell_pred,
                           cell_ids = cell_ids,
                           pixel_id = pixel_id,
                           covdt = covdt){
  #set cell pred as a data table, and rename things
  if (!inherits(x = cell_pred, 'data.table')) {
    cell_pred = as.data.table(cell_pred)
    names(cell_pred) = paste0('V',1:ncol(cell_pred))
  }
  
  cell_pred[, cell_pred_id := .I] #cell_pred object ID
  cell_pred[,cell_id := rep(cell_ids, times = nrow(cell_pred) / length(cell_ids))]  #cell id references the africa map
  cell_pred[,pixel_id := rep(pixel_id, times = nrow(cell_pred) / length(pixel_id))] #pixel id references the regional map  
  
  #add population, year and potentially the stackers
  cell_pred = cbind(cell_pred, covdt[,c('year', 'pop'),with = F])
  
  #make sure it behaved
  stopifnot(any(!(cell_pred[,pixel_id] != rep.int(covdt[,pixel_id], 18))))
  
  return(cell_pred)
}

#' @title Calculate population scalars
#' @description this calculates the population raking factors needed so that prevalence rates and counts will align with GBD
#'
#' @param cell_pred the cell pred object that has been prepped and linked
#' @param age_group the GBD age_group you are modeling
#' @param connector the list of raking geographies
#' @param sex the GBD sex_id you are modeling
#' @param sharedir the share directory for the indicator group, used to save the population raking factors
#' @param run_date your run date
#' @param indicator the modeled indicator
#' @param stratum  the region you are modeling over, used for saving the population raking factors
#' @param gbd_pops output from central code "get_population" function, passed from main function
#' @param custom_output_folder Output the rake factors to custom folder path if specified. Default: NULL
#' 
#' @return a table of population raking facotrs for each of the raking geographies used.  This is so that
#' the fractionally aggregated worldpop rasters for a given age and sex group in a geography year will equal
#' the GBD population for that age and sex group in that geography in that year.
#'
#' @export
calculate_pop_scalars <- function(cell_pred = cell_pred,
                                  age_group = age_group,
                                  connector = connector,
                                  sex = sex_id,
                                  sharedir = sharedir,
                                  run_date = run_date,
                                  indicator = indicator,
                                  stratum = stratum,
                                  gbd_pops = gbd_pops,
                                  custom_output_folder = NULL){
  
  # do calculations!
  rake_geo_pop <- cell_pred[, lapply(c("pop"), function(x) sum(get(x), na.rm = T)), by = c("year", "location_id")]
  rake_geo_pop$world_pop_unraked <- rake_geo_pop$V1
  loc_ids <- unique(connector$location_id)
  
  # adjust to be GBD populations
  scalars <- merge(rake_geo_pop, gbd_pops, by.x = c("location_id", "year"), by.y = c("location_id", "year_id"), all.x = T)
  # Where there are no GBD populations (e.g. ESH, GUF), populate population column with worldpop
  scalars[is.na(population), population := world_pop_unraked]
  scalars[, pop_scalar := population / world_pop_unraked]
  
  # for records
  scalars$gbd_pop <- scalars$population
  scalars$population <- NULL
  
  # trim Scalars object
  scalars$age_group_id <- NULL
  scalars$sex_id <- NULL
  scalars$run_id <- NULL
  scalars$V1 <- NULL
  
  # Fill NAs with 1 BUT WITH A SOLEMN WARNING
  if(nrow(scalars[is.na(pop_scalar)]) >= 1 ) {
    warning("You have NAs in your population raking factor. Please check your GBD populations to make sure you didn't miss any locations.")
    warning("Forcing those missing raking factors to 1")
    scalars[is.na(pop_scalar), pop_scalar:= 1]
  }

  if (!is.null(custom_output_folder)) {
    write.csv(scalars, file = paste0(custom_output_folder, "/", indicator, "_", stratum, "_pop_rf.csv"))
  } else {
    write.csv(scalars, file = paste0(sharedir, "<<<< FILEPATH REDACTED >>>>", indicator, "_", stratum, "_pop_rf.csv"))
  }

  return(scalars)
}


#' @title Calculate fractional raking factor
#' @description This calculates the raking factors, it takes a rates cell_pred which is the main difference between it and the
#' get_fractional_rf function above.
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
                                     h = 0,
                                     m = 'prevalence') {
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
  
  # saving the raking factors
  if (!is.null(custom_output_folder)) {
  	write.csv(fractional_rf, file = paste0(custom_output_folder, "/", indicator, "_", stratum, "_", m, "_rf_", h, ".csv"))
  } else {
  	write.csv(fractional_rf, file = paste0(sharedir, "<<<< FILEPATH REDACTED >>>>", indicator, "_", stratum, "_", m, "_rf_", h, ".csv"))  
  }


  
  return(fractional_rf)
}





#' @title Get GBD population for use in fractional raking
#' @description Get GBD population for use in fractional raking
#'
#' @param pop_measure WorldPop measure, default: "a0004t"
#' @param reg Region
#' @param year_list the years to get population for
#'
#' @return A data.table with \code{c('location_id', 'year_id', 'sex_id', 'run_id', 'age_group_id', 'population')}
#' @note the \code{age_group_id} returned is equal to the value of \code{pop_measure}
#'
#' @export
prep_gbd_pops_for_fraxrake <- function(pop_measure = "a0004t", reg, year_list, ...) {
  if (class(year_list) == "character") year_list <- eval(parse(text = year_list))
  
  ## Get age group ID from pop_measure
  age_GBD <- get_age_group_from_worldpop(pop_measure)
  
  ## Get location ID from GBD
  loc_GBD <- unique(c(get_gbd_locs(reg)$location_id, get_gbd_locs(reg, rake_subnational = FALSE)$location_id))
  
  ## We need to use GBD's get_populations, NOT the LBD one, because of deprecated DB stuff
  source(paste0(CC_ENV_DIR, "/get_population.R"))
  gbd_pops <- get_population(age_group_id = age_GBD, location_id = loc_GBD, sex_id = 3, year_id = year_list, ...)
  
  ## Reduce to a single age group:
  gbd_pops <- gbd_pops[, .(population = sum(population)), by = c("location_id", "year_id", "sex_id", "run_id")]
  gbd_pops[, age_group_id := pop_measure]
  
  return(gbd_pops)
}


#' @title Create rasters for raking and aggregation
#' @description This snippet of code is taken from \code{rake_cell_pred} so that
#' we can feed this into making summary statistics later on in the fractional raking code
#'
#' @param reg Region
#' @param modeling_shapefile_version Modeling shapefile version
#' @param raking_shapefile_version Raking shapefile version
#' @param field Location field. Default: "loc_id"
#'
#' @return A list with objects: "simple_raster", "simple_polygon", "subset_shape", "new_simple_raster", "new_simple_polygon", "new_subset_shape"
#'
#' @export
prep_shapes_for_raking <- function(reg, modeling_shapefile_version, raking_shapefile_version, field = "loc_id") {
  
  ## Load simple polygon template to model over (in GADM SPACE)
  gaul_list <- get_adm0_codes(reg, shapefile_version = raking_shapefile_version)
  simple_polygon_list <- load_simple_polygon(
    gaul_list = gaul_list, buffer = 0.4,
    shapefile_version = raking_shapefile_version
  )
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  
  ## Load list of raster inputs (pop and simple)
  raster_list <- build_simple_raster_pop(subset_shape)
  simple_raster <- raster_list[["simple_raster"]]
  pop_raster <- raster_list[["pop_raster"]]
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  
  
  ### Craete raked_simple_raster template
  shapefile_path <- get_admin_shapefile(
    admin_level = 0,
    raking = FALSE,
    version = raking_shapefile_version
  )
  
  ## Get GADM codes
  gaul_list <- get_adm0_codes(reg, shapefile_version = raking_shapefile_version)
  
  # If not raking subnationally, use the admin0 shapefile
  message("Loading national raking shapefile")
  
  location_metadata <- get_location_code_mapping(shapefile_version = raking_shapefile_version)
  location_metadata <- location_metadata[, c("GAUL_CODE", "loc_id")]
  
  shapefile <- readOGR(shapefile_path, stringsAsFactors = FALSE, GDAL1_integer64_policy = TRUE)
  shapefile <- shapefile[shapefile$ADM0_CODE %in% gaul_list, ]
  
  # merge on loc_ids for making simple raster
  shapefile <- sp::merge(shapefile, location_metadata, by.x = "ADM0_CODE", by.y = "GAUL_CODE", all.x = T)
  shapefile@data[is.na(shapefile@data$loc_id), "loc_id"] <- -1
  
  
  # get simple raster from new gbd shapefile
  message("Loading raking raster\n")
  new_simple_polygon <- load_simple_polygon(
    gaul_list = NULL, ## doesn't matter since custom_shapefile specified
    buffer = 0.4, custom_shapefile = shapefile
  )
  new_subset_shape <- new_simple_polygon[["subset_shape"]]
  new_simple_polygon <- new_simple_polygon[["spoly_spdf"]]
  
  message("Loading simple raster\n")
  raking_link_table <- build_raking_link_table(shapefile_version = raking_shapefile_version)
  new_raster_list <- build_simple_raster_pop(new_subset_shape, field = field, link_table = raking_link_table)
  new_simple_raster <- new_raster_list[["simple_raster"]]
  
  # get extents of original and simple raster to line up - extend and crop just in case
  new_simple_raster <- extend(new_simple_raster, simple_raster, values = NA)
  new_simple_raster <- crop(new_simple_raster, extent(simple_raster))
  
  
  ## Return a named list of things we want
  return(list(
    "raster_list" = raster_list,
    "simple_raster" = simple_raster,
    "simple_polygon" = simple_polygon,
    "subset_shape" = subset_shape,
    "new_simple_raster" = new_simple_raster,
    "new_simple_polygon" = new_simple_polygon,
    "new_subset_shape" = new_subset_shape,
    "pixel_id" = pixel_id
  ))
}