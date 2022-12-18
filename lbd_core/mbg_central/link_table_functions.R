#' @title Assemble and return metadata for LBD standard raking.
#'
#' @description Assembles metadata from the LBD standard raking shapefile for
#' use in building raster layers based upon pixel coverage and not centroid
#' ownership.
#' 
#' @param shapefile_version Shapefile version. Default: 'current'
#' @param countries_not_to_subnat_rake as it sounds. Replaces subnational metadata with national metadata. Default: NULL
#' @return A list with admin_codes and code_locid
#' \item{admin_codes}{A named list with keys for all administrative levels
#' we rake to, likely "ADM0_CODE" and "ADM1_CODE". The values associated are
#' vectors containing all the admin codes we rake to.}
#' \item{code_locid}{A data.table with two columns: admin_code and loc_id.
#' This is a simplification of the loc_id, ADM0_CODE, and ADM1_CODE fields.
#' E.g., for geographies raked to admin 0 the admin_code column will contain
#' relevant ADM0_CODE corresponding to loc_id values.}
standard_raking_link_pieces <- function(shapefile_version = "current", countries_not_to_subnat_rake = NULL) {
  raking_dbf_path <- get_admin_shapefile(type = 'raking', suffix='.dbf', version = shapefile_version)
  raking_records <- data.table(foreign::read.dbf(raking_dbf_path))
  
  #replace metadata for countries that are specified as not subnational
  if(!is.null(countries_not_to_subnat_rake)){
    subnat_change_codes <- get_adm0_codes(countries_not_to_subnat_rake, shapefile_version = shapefile_version)
    raking_records[ADM0_CODE %in% subnat_change_codes, `:=` (ad_level=0, ADM1_CODE=0, ADM1_NAME=NA)]
    
    code_replacement <- get_gbd_locs("all", rake_subnational = F)
    
    for(code in subnat_change_codes){
      raking_records[ADM0_CODE %in% code, loc_id := code_replacement[ADM_CODE == code]$location_id]
    }
    raking_records <- unique(raking_records)
  }
  
  # key by "ADM0_CODE" or "ADM1_CODE"; returns those codes we rake to
  admin_codes <- list()
  # used to hold sub-data.table values for each admin level we rake to.
  code_locid_pieces <- list()
  i <- 1

  for (level in unique(raking_records$ad_level)) {
    field <- paste0("ADM", level, "_CODE")
    records_at_level <- raking_records[raking_records$ad_level == level, ]

    raking_targets <- records_at_level[[field]]
    admin_codes[[field]] <- raking_targets

    code_locid_pieces[[i]] <- data.table(admin_code = raking_targets,
                                         loc_id = records_at_level$loc_id)
    i <- i + 1
  }
  result <- list(admin_codes = admin_codes,
                 code_locid = do.call(rbind, code_locid_pieces))
  return(result)
}

#' @title Aggregates link_table pixel ownership for geographies at distinct admin levels.
#'
#' @description Normally aggregation of link table values from e.g., admin 2 to
#' admin 0 is straight-forward. This function handles aggregating link table values
#' to heterogeneous administrative levels by providing a information about which
#' geographies need to be aggregated to a particular level. In addition, it also
#' translates from Local Burden of Disease administrative codes to GBD location_ids.
#' Finally, if desired, because this is only used for creating the raking raster
#' at ADM0 level, we sub out all ADM1 and ADM2 values to use their parent 
#' location ID by mapping with the GBD location hierarchy
#'
#' @param shapefile_version Shapefile version. Default: "current"
#' 
#' @param force_adm0 Force link table aggregation to ADM0? This produces output consistent 
#' with \code{build_simple_raster_pop}. It is necessary to ensure pixels are consistently 
#' assigned between different rasters. Set this to FALSE with caution. Default: TRUE
#' 
#' @param countries_not_to_subnat_rake as it sounds. Replaces subnational metadata with national metadata. Default: NULL
#'
#' @return a link table with names "pixel_id", "area_fraction", and "loc_id".
#'
#' @export
#'
build_raking_link_table <- function(shapefile_version = "current", force_adm0 = TRUE, countries_not_to_subnat_rake = NULL) {
  
  ## Get some standard link pieces for raking
  link_pieces <- standard_raking_link_pieces(shapefile_version, countries_not_to_subnat_rake)
  
  ## Get the raking targets
  rake_targets <- link_pieces[["admin_codes"]]
  
  ## Get link table to aggregate
  code_locid <- link_pieces[["code_locid"]]
  
  ## Get our full link table
  link_table <- standard_link_table(shapefile_version)
  
  sub_tables <- lapply(names(rake_targets), function(admin_key){
    codes <- rake_targets[[admin_key]]
    # subset link_table
    links_for_codes <- link_table[link_table[[admin_key]] %in% codes, 
                                  c("pixel_id", "area_fraction", admin_key), with=FALSE]
    # aggregate to that admin level
    aggregated <- links_for_codes[, .(area_fraction = sum(area_fraction)), by = c("pixel_id", admin_key)]
    # add "loc_id" field
    indices <- match(aggregated[[admin_key]], code_locid$admin_code)
    aggregated[, "loc_id"] <- code_locid[indices, "loc_id"]
    # delete `admin_key` field
    aggregated[, (admin_key):=NULL]
    aggregated
  })
  
  # merge results from different admin levels
  raking_table <- rbindlist(sub_tables)
  
  if(force_adm0) {
    # Get the GBD location hierarchies, where
    # we map the subnational countries to their parent ADM0 countries
    ### This is for avoiding corner cases where we have no raking targets ###
    source(paste0(CC_ENV_DIR, "/get_location_metadata.R"))
    gbd_locset_copy <- get_location_metadata(22, decomp_step = "step4", gbd_round_id = 6) #This is GBD2019 - final step
    gbd_locset <- gbd_locset_copy[, .(location_id, ihme_loc_id)]
    
    # GBD subnational locations are identified in ihme_loc_id as: NATIONAL_SUBNATIONALCODE
    # where the first 3 characters of an ihme_loc_id with "_" in it will be the location's 
    # parent ihme_loc_id.
    # We take the parent location out of those.
    # e.g.: "USA_501" -> "USA"
    gbd_locset[ihme_loc_id %like% "_", parent_to_subnat := substr(ihme_loc_id, 1, 3)]
    gbd_locset_submap <- gbd_locset[ihme_loc_id %in% unique(na.omit(gbd_locset$parent_to_subnat)), .(parent_to_subnat = ihme_loc_id, subnat_map = location_id)]
    gbd_locset_submap <- unique(merge(gbd_locset, gbd_locset_submap, "parent_to_subnat")[, .(location_id, subnat_map)])
    
    ## Merge that on to raking_link_table and override all subnats with nats IDs where fixing is necessary
    raking_table <- merge(raking_table,
                          gbd_locset_submap[, .(loc_id = location_id, subnat_map)],
                          by = "loc_id", all.x = TRUE
    )
    raking_table[!is.na(subnat_map), loc_id := subnat_map]
    raking_table[, subnat_map := NULL]
  }
  
  
  return(raking_table)
}


#' @title Load the standard link table and return it
#'
#' @description Loads the lbd_standard_link.rds file associated with an admin shapefile.
#'
#' @param admin_version the administrative shapefile version. Either "current" or
#' a valid YYYY_MM_DD string.
#'
#' @return the data.table stored in the RDS file.
standard_link_table <- function(admin_version) {
  path <- paste0(get_admin_shape_dir(admin_version), "lbd_standard_link.rds")
  return(readRDS(path))
}
