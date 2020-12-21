# ---------------------------------------------------------------------------------------------
# aggregate_child_stackers()
#
# Function that aggregates child stacker results in a population-weighted manner
#
# Inputs:
# Indicator, indicator group, run date, reg (which can be a vector of regions),
# shapefile version, population measure for population-weighted aggregation, age, and holdout
#
# Outputs:
# List of admin 0 (ad0) and admin 1 (ad1) aggregated child stacker results
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
aggregate_child_stackers <- function(indicator,
                                     indicator_group,
                                     run_date,
                                     regions,
                                     modeling_shapefile_version,
                                     pop_measure = 'a0004t',
                                     age = 0,
                                     holdout = 0){
                                 
  # -----------------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------------------
  # Set-up
  
  # Load required functions
  library('data.table')
  library('raster')
  source('<<<< FILEPATH REDACTED >>>>/mbg_central/fractional_raking_functions.R')
  source('<<<< FILEPATH REDACTED >>>>/mbg_central/covariate_functions.R')
  source('<<<< FILEPATH REDACTED >>>>/mbg_central/prep_functions.R')
  source('<<<< FILEPATH REDACTED >>>>/mbg_central/shapefile_functions.R')
  
  # Create empty list to write to
  reg_outputs <- list()
  # --------------------------------------------------------------------------------------
  
  
  # ----------------------------------
  # Begin loop over regions
  
  for (reg in regions) {
    message(paste0('\nLoading stackers and population estimates for: ', reg))
    
    # ---------------------------------------------------------------------------------------------------------------
    # Load rasters and stackers
    
    # Get the simple and new_simple rasters prepped up for us
    raster_outputs <- prep_shapes_for_raking(
      reg = reg,
      modeling_shapefile_version = modeling_shapefile_version,
      raking_shapefile_version = modeling_shapefile_version,
      field = 'loc_id'
    )
    
    # Take out the objects from the list that actually matters to us:
    simple_raster <- raster_outputs[['simple_raster']]
    simple_polygon <- raster_outputs[['simple_polygon']]
    pixel_id <- raster_outputs[['pixel_id']]
    
    # Add stacking children to this process
    covs = fetch_from_rdata(paste0('<<<< FILEPATH REDACTED >>>>', 
                                   run_date, '_bin', age,'_', reg, '_', holdout, '.RData'), 'cov_list')
    fes = fetch_from_rdata(paste0('<<<< FILEPATH REDACTED >>>>', 
                                  run_date, '_bin', age,'_', reg, '_', holdout, '.RData'), 'all_fixed_effects')
    submodels = trimws(strsplit(fes, '+', fixed = T)[[1]])
    covs = covs[submodels]
    
    # Set stacker names
    covnames = names(covs)
    
    # Ensure the dimensions are the same
    for(ccc in covs){
      message(names(ccc))
      message(dim(ccc)[1:2])
      message(dim(simple_raster)[1:2])
      stopifnot(dim(ccc)[1:2] == dim(simple_raster)[1:2])
    }
  
    # Convert raster to data table and clean
    raster_to_dt = function(ras){
      dt <- crop(ras, extent(simple_raster))
      dt <- data.table(extract(dt, pixel_id))
      dt[, pixel_id := pixel_id]
      dt <- melt(dt, id.vars='pixel_id')
      dt[, c('pixel_id', 'variable') := NULL]
      return(dt)
    }
    
    covdt = lapply(covs, raster_to_dt)
    covdt = do.call(what = cbind, covdt)
    setnames(covdt, names(covs))
  
    # Add pixel_id, but make sure that its recycled explicitly as per data.table 1.12.2 guidelines
    covdt[, pixel_id := rep(pixel_id, times = nrow(covdt) / length(pixel_id))]
    
    # Add year to covdt
    yyy = as.vector(unlist(lapply(min(year_list):max(year_list), function(x) rep.int(x, times = length(pixel_id)))))
    covdt[, year := yyy]
    
    # Free up a bit of space
    rm(covs)
    # ---------------------------------------------------------------------------------------------------------------
    
    
    # ----------------------------------------------------------------------------------------
    # Add in population
    
    # Pull worldpop estimates from covariate database
    pop <- load_and_crop_covariates_annual(covs = 'worldpop',
                                           measures = pop_measure,
                                           simple_polygon = simple_polygon,
                                           start_year  = min(year_list),
                                           end_year    = max(year_list),
                                           interval_mo = 12,
                                           agebin=1)
    
    # Convert to data table
    pop <- raster_to_dt(pop$worldpop)
    
    # Add to stacker results
    covdt <- cbind(covdt, pop)
    # ----------------------------------------------------------------------------------------
    
    
    # -------------------------------------------------------------------------------------------------------
    # Use link table to assign admin codes to stacker results
    
    # Load the cell id to admin units link
    link_table <- get_link_table(simple_raster, shapefile_version = modeling_shapefile_version)
    
    # Prep the cell_pred and link table to be linked by making sure they have the appropriate identifiers
    link <- prep_link_table(
      link_table = link_table,
      simple_raster = simple_raster,
      pixel_id = pixel_id
    )
    
    # Add link table version of pixel id to stacker data table
    covdt$pixel_id <- link_table$pixel_ids
    
    # Merge identifiers with stacker data frame
    children <- merge(covdt, link, by = 'pixel_id', allow.cartesian = T)
    
    # Save regional results in list
    reg_outputs[[reg]] <- children
    
    # Free up a bit of space
    rm(covdt, link, link_table, pop, children)
    # -------------------------------------------------------------------------------------------------------
    
  }
  
  # End loop over regions
  # ----------------------------------
  
  
  # -------------------------------------------------------------------
  # Aggregate
  
  # Convert list to data table
  children <- rbindlist(reg_outputs, use.names = T)
  
  # Get population for calculating weighted average
  children[, pop := value*area_fraction]
  
  # Weighted mean by admin 0
  child0 <- children[, lapply(.SD, weighted.mean, w = pop, na.rm = T),
                     by = c('ADM0_CODE', 'year'),
                     .SDcols = covnames]
  
  # Weighted mean by admin 1
  child1 <- children[, lapply(.SD, weighted.mean, w = pop, na.rm = T),
                     by = c('ADM0_CODE', 'ADM1_CODE', 'year'),
                     .SDcols = covnames]
  
  # Create output list
  children <- list(child0, child1)
  names(children) <- c('ad0', 'ad1')
  # -------------------------------------------------------------------

  
  # -------------------------------------------------------------------------------------------------------
  # End function
  return(children)
  
}
# -----------------------------------------------------------------------------------
