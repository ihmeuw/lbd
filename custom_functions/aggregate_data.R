# ---------------------------------------------------------------------------------------------
# aggregate_input_data()
#
# Function that aggregates model input data in a population-weighted manner
#
# Inputs:
# Indicator, indicator group, run date, reg (which can be a vector of regions),
# and shapefile version
#
# Outputs:
# List of admin 0 (ad0) and admin 1 (ad1) aggregated input data
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
aggregate_input_data <- function(indicator,
                                 indicator_group,
                                 run_date,
                                 regions,
                                 modeling_shapefile_version) {
  # -----------------------------------------------------------------------------------
  

  # ---------------------------------------------------------------
  # Set-up
  
  # Load packages
  libs <- c('data.table', 'dplyr', 'raster', 'sf', 'fasterize')
  lapply(libs, library, character.only = TRUE)
  
  # Set arguments
  code <- '<<<< FILEPATH REDACTED >>>>'
  outputdir <- '<<<< FILEPATH REDACTED >>>>'
  
  # Create empty lists to write to
  reg_ad0_inputs <- list()
  reg_ad1_inputs <- list()
  # ---------------------------------------------------------------
  
  
  # ----------------------------------
  # Begin loop over regions
  
  for (reg in regions) {
    message(paste0('Aggregating input data for: ', reg))
  
    # ------------------------------------------------------------------------------
    # Load data and rasters
    
    # Load input data
    setwd('<<<< FILEPATH REDACTED >>>>')
    input_data <- fread(paste0(indicator, '.csv'))
    input_data[, prop := get(indicator)/N]
    
    # Reading in shapefiles & population raster
    file_dir <- paste0('<<<< FILEPATH REDACTED >>>>',
                       modeling_shapefile_version, '/')
    a0_shp <- readRDS(paste0(file_dir,
                             'lbd_standard_admin_0.rds'))
    a1_shp <- readRDS(paste0(file_dir,
                             'lbd_standard_admin_1.rds'))
    pop_ras <- raster(
      paste0('<<<< FILEPATH REDACTED >>>>',
             'worldpop_total_1y_2010_00_00.tif'))
    # ------------------------------------------------------------------------------
    
    
    # ----------------------------------------------------------------------------------------------------------------
    # Aggregate data
    
    # Input data and admin_draws formatting function
    format_id <- function(df, indi = indicator, repo = code, region = reg,
                          shp_version = modeling_shapefile_version) {
      source(paste0(code, '/mbg_central/prep_functions.R'))
      source(paste0(code, '/mbg_central/shapefile_functions.R'))
      gaul_list <- get_adm0_codes(region, shapefile_version = shp_version)
      
      lookup <- fread('<<<< FILEPATH REDACTED >>>>/stage_master_list.csv')
      isos <- lookup[gadm_geoid %in% gaul_list,iso3]
      df <- df[df$country %in% isos,]
      
      return(df)
    }
    
    # Format shps as rasters
    a0_shp$ADM0_CODE <- as.numeric(as.character(a0_shp$ADM0_CODE))
    a1_shp$ADM1_CODE <- as.numeric(as.character(a1_shp$ADM1_CODE))
    a0_raster <- fasterize(st_as_sf(a0_shp), pop_ras, field = 'ADM0_CODE')
    a1_raster <- fasterize(st_as_sf(a1_shp), pop_ras, field = 'ADM1_CODE')
    
    # Format input data & extract adminIDs
    input_data <- format_id(input_data)
    input_data$ADM1_CODE <- raster::extract(a1_raster, dplyr::select(input_data, 
                                                                     longitude, latitude))
    input_data$ADM0_CODE <- raster::extract(a0_raster, dplyr::select(input_data, 
                                                                     longitude, latitude))
    
    # Summarize input data by admin 0 IDs
    a0_input <- as.data.table(input_data)
    a0_input[, input_mean := weighted.mean(prop, sum_of_sample_weights*weight),
             by = c('nid', 'year', 'ADM0_CODE')]
    a0_input[, input_ss := sum(weight*N),
             by = c('nid', 'year', 'ADM0_CODE')]
    
    # Summarize input data by admin 1 IDs                                  
    a1_input <- as.data.table(input_data)
    a1_input[, input_mean := weighted.mean(prop, sum_of_sample_weights*weight),
             by = c('nid', 'year', 'ADM0_CODE', 'ADM1_CODE')]
    a1_input[, input_ss := sum(weight*N),
             by = c('nid', 'year', 'ADM0_CODE', 'ADM1_CODE')]
    
    # Add regional data to parent lists
    reg_ad0_inputs[[reg]] <- unique(a0_input[, c('nid', 'year', 'ADM0_CODE', 'input_mean', 'input_ss')])
    reg_ad1_inputs[[reg]] <- unique(a1_input[, c('nid', 'year', 'ADM0_CODE', 'ADM1_CODE', 'input_mean', 'input_ss')])
    # ----------------------------------------------------------------------------------------------------------------
    
  }
  
  # End loop over regions
  # ----------------------------------
  

  # -------------------------------------------------
  # Clean up
  
  # Bind regional data tables to create ad0/ad1 list
  reg_inputs <- list(rbindlist(reg_ad0_inputs),
                     rbindlist(reg_ad1_inputs))
  names(reg_inputs) <- c('ad0', 'ad1')
  # -------------------------------------------------
  
  
  # -------------------------------------------------------------------------------------------------------
  # End function
  return(reg_inputs)
  
}
# -----------------------------------------------------------------------------------
