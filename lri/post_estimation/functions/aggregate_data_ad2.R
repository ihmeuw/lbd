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
# List of admin 2 (ad2) aggregated input data
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
aggregate_input_data_ad2 <- function(indicator,
                                     indicator_group,
                                     run_date,
                                     regions,
                                     modeling_shapefile_version,
                                     file_name = indicator) {
  
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
  reg_ad2_inputs <- list()
 
  # ---------------------------------------------------------------
  
  
  # ----------------------------------
  # Begin loop over regions
  
  for (reg in regions) {
    message(paste0('Aggregating input data for: ', reg))
    
    # ------------------------------------------------------------------------------
    # Load data and rasters
    
    # Load input data
    setwd('<<<< FILEPATH REDACTED >>>>')
    input_data <- fread(paste0(file_name, '.csv'))
    input_data[, prop := get(indicator_group)/N]
    
    # Reading in shapefiles & population raster
    file_dir <- '<<<< FILEPATH REDACTED >>>>'
    a2_shp <- readRDS(paste0(file_dir,
                             'lbd_standard_admin_2.rds'))
    pop_ras <- raster('<<<< FILEPATH REDACTED >>>>')
    # ------------------------------------------------------------------------------
    
    
    # ----------------------------------------------------------------------------------------------------------------
    # Aggregate data
    
    # Input data and admin_draws formatting function
    format_id <- function(df, indi = indicator, repo = code, region = reg,
                          shp_version = modeling_shapefile_version) {
      source(paste0(code, '/mbg_central/prep_functions.R'))
      source(paste0(code, '/mbg_central/shapefile_functions.R'))
      gaul_list <- get_adm0_codes(region, shapefile_version = shp_version)
      
      lookup <- fread('<<<< FILEPATH REDACTED >>>>')
      isos <- lookup[gadm_geoid %in% gaul_list,iso3]
      df <- df[df$country %in% isos,]
      
      return(df)
    }
    
    # Format shps as rasters
    a2_shp$ADM2_CODE <- as.numeric(as.character(a2_shp$ADM2_CODE))
    a2_raster <- fasterize(st_as_sf(a2_shp), pop_ras, field = 'ADM2_CODE')

    # Format input data & extract adminIDs
    input_data <- format_id(input_data)
    input_data$ADM2_CODE <- raster::extract(a2_raster, dplyr::select(input_data, 
                                                                     longitude, latitude))
    
    # Summarize input data by admin 2 IDs                                  
    a2_input <- as.data.table(input_data)
    a2_input[, input_mean := weighted.mean(prop, sum_of_sample_weights*weight),
             by = c('nid', 'year', 'ADM2_CODE')]
    a2_input[, input_ss := sum(weight*N),
             by = c('nid', 'year', 'ADM2_CODE')]
    
    # Add regional data to parent lists
    reg_ad2_inputs[[reg]] <- unique(a2_input[, c('nid', 'year', 'ADM2_CODE', 'input_mean', 'input_ss')])
    # ----------------------------------------------------------------------------------------------------------------
    
  }
  
  # End loop over regions
  # ----------------------------------
  
  
  # -------------------------------------------------
  # Clean up
  
  # Bind regional data tables to create ad2 list
  reg_inputs <- list(rbindlist(reg_ad2_inputs))
  names(reg_inputs) <- c('ad2')
  
  # -------------------------------------------------
  
  
  # -------------------------------------------------------------------------------------------------------
  # End function
  return(reg_inputs)
  
}
# -----------------------------------------------------------------------------------