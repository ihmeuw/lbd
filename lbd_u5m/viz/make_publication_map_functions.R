## #############################################################################
## 
## FUNCTIONS FOR MAKING PUBLICATION-READY MAPS
## 
## Purpose: Pull output data and create a large series of publication-ready maps
## NOTE: Please run this on the geospatial R Singularity image!
## 
## #############################################################################

library(colorspace)
library(data.table)
library(dplyr)
library(fasterize)
library(ggplot2)
library(raster)
library(RColorBrewer)
library(stringr)
library(viridis)


## Program constants (DO NOT ALTER) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
U5M_BASE_PATH <- '<<<< FILEPATH REDACTED >>>>'
U5M_ADMIN_LEVELS  <- c(0, 1, 2)
U5M_AGE_GROUPS    <- list('under5', 'infant', 'neonatal')
U5M_RAKING_LABELS <- c('raked','unraked')
U5M_STATS <- c(
  'mean','lower','upper','deathcounts_mean'
)
## Stage 2 mapping bounds
S2_BOUNDS <- list(
  'lat_min'  = -57,
  'lat_max'  = 63,
  'long_min' = -129,
  'long_max' = 165
)

## pull_prep_admin_data() ----------------------------------------------------->
#' 
#' @title Pull and prep administrative data for mapping
#' 
#' @description Pulls all data for a given administrative unit and preps it for
#'   mapping:
#'   - Pulls mean, lower, upper data from file (mortality and death counts)
#'   - Reshapes wide
#'   - Calculates AROC of the mean maps for q and deathcounts
#'   - Calculates NN/U5 and Inf/U5 for all years
#'   
#' @param run_date Underscored model run date, used to point to model run folder
#' @param adm_level Integer: which admin level should be prepped and pulled?
#' 
#' @return Prepped data.table with all inputs joined and reshaped wide
#' 
pull_prep_admin_data <- function(run_date, adm_level){
  message(sprintf("Pulling all data aggregated to the admin%s level...",adm_level))
  # Template filepath for all aggregated CSVs (and column names!)
  template_dt <- CJ(
    ag = U5M_AGE_GROUPS,
    rake = U5M_RAKING_LABELS,
    stat = U5M_STATS
  )
  template_dt[,
    fp := sprintf('<<<< FILEPATH REDACTED >>>>', U5M_BASE_PATH, 
                  ag, run_date, ag, stat, rake, adm_level)
  ]
  template_dt[, col_name := paste(ag, stat, rake, sep='_')]
  # Helper function to determine AROC given one of the input tables
  get_agg_aroc <- function(in_data){
    first_year_data <- in_data[year==min(year), .(admin_code, value)]
    last_year_data <- in_data[year==max(year), .(admin_code, value)]
    num_years <- max(in_data[,year]) - min(in_data[,year]) + 1
    aroc_data <- merge(
      x        = first_year_data,
      y        = last_year_data,
      by       = c('admin_code'),
      suffixes = c('_first','_last')
    )
    aroc_data[value_last == 0, value_last := 1E-8]
    aroc_data[value_first == 0, value_first := 1E-8]
    aroc_data[, value := log(value_last / value_first)/(num_years)]
    aroc_data <- aroc_data[, .(value, admin_code)]
    aroc_data[, year := 'aroc']
    return(aroc_data)
  }
  # Pull data for all age groups, measures, and raking statuses
  # Bind together into a single data.table
  data_long_list <- lapply(
    as.list(1:nrow(template_dt)),
    function(row_i){
      if(file.exists(template_dt[row_i, fp])){
        # Read in data for this ag/measure/raked combo
        in_data <- fread(template_dt[row_i, fp])
        setnames(in_data, paste0('ADM',adm_level,"_CODE"), 'admin_code')
        # Get the annualized rate of change and rbind it to the data
        in_data <- rbind(in_data, get_agg_aroc(in_data), fill=TRUE)
        # Set the 'year' column as a character
        in_data[, year := as.character(year)]
        # Create a col_name column based on year and ag/measure/raked combo
        in_data[, wide_name := paste0(template_dt[row_i, col_name], '_', year)]
        # Keep only needed fields
        in_data <- in_data[, .(admin_code, value, wide_name, year)]
        setkeyv(in_data, 'admin_code')
        return(in_data)
      } else {
        return(NULL)
      }
    }
  )
  full_data_long <- data_long_list[ 
    sapply(data_long_list, function(i) 'data.table' %in% class(i))
  ] %>% rbindlist
  # Remove any "arocs" that are not means
  full_data_long <- full_data_long[!grepl('aroc',wide_name)|grepl('mean',wide_name)]
  # Get the start and end years
  start_year <- suppressWarnings(min(as.numeric(full_data_long[,year]), na.rm=TRUE))
  end_year   <- suppressWarnings(max(as.numeric(full_data_long[,year]), na.rm=TRUE))
  full_data_long[,year := NULL]
  # Reshape wide using the column name
  full_data_wide <- dcast(
    data=full_data_long,
    formula = 'admin_code ~ wide_name',
    value.var='value',
    fun.aggregate=mean 
  )
  for (yr in start_year:end_year){
    nn_col <- paste0('neonatal_mean_raked_',yr)
    inf_col <- paste0('infant_mean_raked_',yr)
    u5_col <- paste0('under5_mean_raked_',yr)
    # Get NN / U5M mortality rate ratios
    nn_u5_col <- paste0('nn_u5_ratio_',yr)
    if ( all(c(nn_col,u5_col) %in% names(full_data_wide)) ){
      full_data_wide[, (nn_u5_col) := get(nn_col) / get(u5_col) ]
    }
    # Get infant / U5M mortality rate ratios
    inf_u5_col <- paste0('inf_u5_ratio_',yr)
    if ( all(c(inf_col,u5_col) %in% names(full_data_wide)) ){
      full_data_wide[, (inf_u5_col) := get(inf_col) / get(u5_col) ]
    }
  }
  # There should now be a huge data.table with (num gauls) rows and columns for
  #  every year/measure/age group/raked combo
  full_data_wide[, adm_level := adm_level]
  message("Done pulling and reshaping administrative data.")
  return(full_data_wide)
}



## pull_prep_bra_mex_data() --------------------------------------------------->
#' 
#' @title Pull and prep Mexico and Brazil admin1 data for mapping
#' 
#' @description Pulls all data for Brazilian and Mexican states and preps it for
#'   mapping:
#'   - GBD 2017 q data (mean, lower, upper) from mortality tables
#'   - GBD 2017 D data (mean, lower, upper) from get_outputs
#'   - Reshapes wide
#'   - Calculates AROC of the mean maps for q and deathcounts
#'   - Calculates NN/U5 and Inf/U5 for all years
#'   
#' @param shp_version [default 'current'] Character vector telling which shp 
#'   date should be pulled (to merge on state data). 
#'   
#' @return Prepped data.table with all inputs joined and reshaped wide
#' 
pull_prep_bra_mex_data <- function(shp_version){
  message("Pulling all admin1 data for Brazil and Mexico...")
  ## Load GBD shared functions
  source('<<<< FILEPATH REDACTED >>>>')
  source('<<<< FILEPATH REDACTED >>>>')
  source('<<<< FILEPATH REDACTED >>>>')
  ## Pull location IDs and names for all Brazilian and Mexican states
  locs <- get_location_metadata(
    location_set_id = 21,
    gbd_round_id = 5
  )
  mb <- locs[ location_name %in% c("Brazil","Mexico"), .(location_id,location_name)]
  setnames(mb, c('location_id','location_name'), c('parent_id','parent_name'))
  mb_state <- merge(
    x = locs,
    y = mb,
    by = c('parent_id')
  )[, .(location_id,location_ascii_name,parent_name,ihme_loc_id)]
  ## Pull GAUL/GADM codes for Brazil and Mexico states
  mb_admins <- as.data.table(foreign::read.dbf(
    gsub("shp$","dbf",get_admin_shapefile(admin_level=1, version=shp_version))
  ))[ ADM0_NAME %in% c('Brazil', 'Mexico'), .(ADM1_NAME, ADM1_CODE, ADM0_NAME)]
  setnames(
    mb_admins, 
    c('ADM0_NAME','ADM1_NAME','ADM1_CODE'), 
    c('parent_name','location_ascii_name','admin_code')
  )
  # Convert to ASCII for merge
  mb_admins$location_ascii_name <- stringi::stri_trans_general(
    mb_admins$location_ascii_name, 
    'Latin-ASCII'
  )
  # Run a few more manual conversions
  mb_admins[location_ascii_name=='Michoacan De Ocampo',location_ascii_name:='Michoacan de Ocampo']
  mb_admins[location_ascii_name=='Queretaro Arteaga',location_ascii_name:='Queretaro']
  mb_admins[location_ascii_name=='Coahuila De Zaragoza',location_ascii_name:='Coahuila']
  mb_admins[
    location_ascii_name=='Veracruz De Ignacio De La Llave',
    location_ascii_name:='Veracruz de Ignacio de la Llave'
  ]
  mb_admins[
    location_ascii_name=='Distrito Federal' & parent_name=='Mexico',
    location_ascii_name:='Mexico City'
  ]
  mb_state <- merge(
    x = mb_admins,
    y = mb_state,
    by = c('location_ascii_name','parent_name'),
    all = TRUE
  )
  if(nrow(mb_state) != nrow(mb_admins)) stop("Issue with admin merge")
  ## Pull GBD 2017 Q
  mb_locs <- mb_state$location_id
  gbd_q_inf_u5 <- get_life_table(
    age_group_id = c(1, 28),
    location_id = mb_locs,
    sex_id=3,
    year_id=2000:2017,
    life_table_parameter_id=3, # nQx
    gbd_round_id=5
  )[, .(location_id, year_id, age_group_id, mean)]
  gbd_q_nn <- fread(
    "<<<< FILEPATH REDACTED >>>>"
  )[
    (year_id >= 2000) & (year_id <= 2017) & (age_group_id == 42) & (sex_id==3) ,
  ][
    location_id %in% mb_locs, 
    .(location_id, year_id, age_group_id, mean, lower, upper)
  ]
  ## Pull GBD 2017 D
  gbd_d <- get_outputs(
    'cause',
    cause_id = 294, # All causes
    age_group_id = c(1, 28, 42),
    sex_id = 3, 
    location_id = mb_locs,
    year_id = 2000:2017,
    measure_id = 1, # Deaths
    metric_id = 1, # Counts
    gbd_round_id = 5
  )[, .(location_id, year_id, age_group_id, val)]
  setnames(gbd_d, 'val', 'value')
  gbd_d[, variable := 'deathcounts_mean']
  ## Combine all and reshape wide
  gbd_q_nn_long <- melt(
    gbd_q_nn,
    id.vars=c('location_id','year_id','age_group_id'),
    measure.vars=c('mean','lower','upper')
  )
  gbd_q_inf_u5_long <- melt(
    gbd_q_inf_u5,
    id.vars=c('location_id','year_id','age_group_id'),
    measure.vars=c('mean')
  )
  all_long <- rbindlist(
    list(gbd_d, gbd_q_nn_long, gbd_q_inf_u5_long),
    use.names=TRUE
  )
  all_long <- merge(
    x = all_long,
    y = data.table(
      age_group_id = c(1, 28, 42),
      age_group_name = c('under5','infant','neonatal')
    ),
    by = c('age_group_id')
  )
  all_long[,wide_name:=paste(age_group_name,variable,'raked',year_id,sep='_')]
  all_long <- all_long[, .(location_id, value, wide_name)]
  data_wide <- dcast(
    data=all_long,
    formula = 'location_id ~ wide_name',
    value.var='value',
    fun.aggregate=mean # TODO: Make sure this is not returning multiple values
  )
  ## Calculate AROC
  all_topics <- gsub('_2017$', '', grep('_2017$',names(data_wide), value=T))
  for(tp in all_topics){
    fy <- paste0(tp,'_2000')
    ly <- paste0(tp,'_2017')
    arc <- paste0(tp,'_aroc')
    data_wide[, (arc) := log( get(ly) / get(fy) ) / 18 ]      
  }
  ## Calculate NN/U5 and INF/U5
  for(yr in 2000:2017){
    nn_col <- paste0('neonatal_mean_raked_',yr)
    inf_col <- paste0('infant_mean_raked_',yr)
    u5_col <- paste0('under5_mean_raked_',yr)
    inf_u5_col <- paste0('inf_u5_ratio_',yr)
    nn_u5_col <- paste0('nn_u5_ratio_',yr)
    data_wide[, (nn_u5_col) := get(nn_col) / get(u5_col) ]
    data_wide[, (inf_u5_col) := get(inf_col) / get(u5_col) ]
  }
  ## Merge onto admin codes 
  full_data_wide <- merge(
    x = data_wide,
    y = mb_state[, .(location_id, admin_code, parent_name)],
    by = 'location_id'
  )
  # Return full data
  message("Done pulling admin1 data for Brazil and Mexico.")
  return(full_data_wide)  
}


## pull_format_admin_shp() ---------------------------------------------------->
#' 
#' @title Pull and format administrative shapefiles
#' 
#' @description Pre-process administrative shapefiles at a given level
#' 
#' @param adm_level Which administrative divisions (0/1/2) should be pulled?
#' @param shp_version [default 'current'] Character vector telling which shp 
#'   date should be pulled. 
#' @param simp_prop [default NULL] What proportion of points should be kept in 
#'   the simplified shapefile? The default, NULL, does not simplify.
#' @param stage2_only [default FALSE] Limit to results from Stage 2?
#' @param include_non_gbd [default TRUE] Should non-GBD locations (namely ESH
#'   Western Sahara and GUF French Guiana in Stage 2) be included in the output
#'   shapefile?
#' 
#' @return a prepped shapefile ready for ggplotting
#' 
pull_format_admin_shp <- function(
  adm_level,
  shp_version = 'current',
  simp_prop   = NULL,
  stage2_only = FALSE,
  include_non_gbd = TRUE
  ){
  # Load the shapefile, simplifying Admin 0 shapefiles a bit by default
  if(is.null(simp_prop) & (adm_level==0)) simp_prop <- 0.05
  shp <- fast_load_shapefile(
    shp_path = get_admin_shapefile(admin_level=adm_level, version=shp_version),
    simplify_tol=simp_prop
  )
  # Pull the admin0 lookup table 
  lookup_table <- load_adm0_lookup_table()
  adm_field <- ifelse(
    detect_adm_shapefile_date_type(
      shpfile_path=get_admin_shapefile(admin_level=0, version=shp_version)
    )$shpfile_type=='gaul',
    'GAUL_CODE',
    'gadm_geoid'
  )
  setnames(lookup_table, adm_field, 'ADM0_CODE')
  if(include_non_gbd==FALSE){
    # Drop all non-GBD locations, which will not have a valid location ID
    lookup_table <- lookup_table[ loc_id > 0, ]
  }
  if(stage2_only){
    lookup_table <- lookup_table[ Stage %in% c('1','2a','2b'), ]
  }
  # Fortify for ggplotting, but don't add fields
  shp$admin_code <- shp@data[, paste0('ADM',adm_level,'_CODE')]
  shp@data       <- shp@data[, c('ADM0_CODE','admin_code')]
  shp_fort <- prep_shp_data_for_mapping(
    shp       = shp,
    dataset   = lookup_table[, .(ADM0_CODE, location_name)],
    merge_var ='ADM0_CODE'
  )
  # Subset to the Stage 2 boundaries
  shp_fort[lat  < S2_BOUNDS$lat_min,  lat  := S2_BOUNDS$lat_min  ]
  shp_fort[lat  > S2_BOUNDS$lat_max,  lat  := S2_BOUNDS$lat_max  ]
  shp_fort[long < S2_BOUNDS$long_min, long := S2_BOUNDS$long_min ]
  shp_fort[long > S2_BOUNDS$long_max, long := S2_BOUNDS$long_max ]
  message(sprintf("Done prepping admin%s data.\n",adm_level))
  return(shp_fort)
}


## pull_format_disputed_shp() ------------------------------------------------->
#' 
#' @title Pull and format the the disputed territories shapefile
#' 
#' @description Pre-process the shapefile that will be used to represent the
#'   outlines of disputed territories
#' 
#' @param shp_version [default 'current'] Character vector telling which shp 
#'   date should be pulled. 
#' @param simp_prop [default NULL] What proportion of points should be kept in 
#'   the simplified shapefile? The default, NULL, does not simplify.
#' 
#' @return a prepped border shapefile ready for ggplotting
#' 
pull_format_disputed_shp <- function(
  shp_version='current', 
  simp_prop=NULL
  ){
  # Pull from file
  shp <- fast_load_shapefile(
    shp_path = get_admin_shapefile(type='disputed_mask', version=shp_version),
    simplify_tol=simp_prop
  )
  # Format for ggplotting
  shp_fort <- prep_shp_data_for_mapping(
    shp       = shp,
    dataset   = shp@data[, c("ADM0_CODE","ADM0_NAME")],
    merge_var ='ADM0_CODE'
  )
  message("Done pulling and prepping disputed border shapefile.")
  return(shp_fort)
}


## make_raster_aroc() --------------------------------------------------------->
#' 
#' @title Make annualized rate of change raster
#' 
#' @description Given a raster brick, create an annualized rate of change raster
#'   using only the first raster, final raster, and number of layers
#' 
#' @param r_brick Raster brick to create an AROC raster with
#' 
#' @return A RasterLayer containing the AROC raster
#' 
make_raster_aroc <- function(r_brick){
  # Stop if the input is not a RasterBrick
  if(!('RasterBrick' %in% class(r_brick))) stop("The input is not a RasterBrick.")
  # Denominator: Number of layers (years)
  n_years <- dim(r_brick)[3]
  # Remove zeroes, which would create errors
  first_yr_data <- r_brick[[1]]
  first_yr_data[ first_yr_data <= 0 ] <- 1E-8
  last_yr_data  <- r_brick[[ dim(r_brick)[3] ]]
  last_yr_data[ last_yr_data <= 0 ] <- 1E-8
  # Calculate annualized rate of change
  aroc_raster <- log(last_yr_data / first_yr_data) / n_years
  return(aroc_raster)
}


## pull_prep_raster_data() ---------------------------------------------------->
#' 
#' @title Pull and prep raster data for mapping 
#' 
#' @description Pulls raster data for a given U5M run date
#' 
#' @param run_date Underscored model run date, used to point to model run folder
#' @param start_year [default 2000] first year of model analysis
#' @param end_year [default 2017] final year of model analysis
#' 
#' @return Named list with two items, "annual" and "aroc", containing all annual
#'   raster bricks and annualized rate-of-change RasterLayers, respectively.
#' 
pull_prep_raster_data <- function(
  run_date, 
  start_year = 2000,
  end_year   = 2017
  ){
  message("Pulling annual rasters...")
  template_dt <- CJ(
    ag = U5M_AGE_GROUPS,
    rake = U5M_RAKING_LABELS,
    stat = U5M_STATS
  )
  # Create filepaths and raster names
  template_dt[,
    fp := sprintf('<<<< FILEPATH REDACTED >>>>', U5M_BASE_PATH, 
                  ag, run_date, ag, stat, rake, start_year, end_year)
  ]
  template_dt[, r_name := paste(ag, stat, rake, sep="_")]

  # Pull all raster bricks, dropping raster bricks that don't actually exist
  r_list <- as.list(template_dt[,r_name])
  names(r_list) <- template_dt[,r_name]
  rasters_annual <- lapply(
    r_list,
    function(raster_title){
      in_fp <- template_dt[r_name==raster_title, fp]
      if(file.exists(in_fp)) return( brick(in_fp) ) else return(NULL)
    }
  )
  rasters_annual <- rasters_annual[ 
    sapply(rasters_annual, function(x) 'RasterBrick' %in% class(x))
  ]
  # Make NN/U5 and INF/U5 rasters
  message("Calculating NN/under-5 and infant/under-5 ratio rasters...")
  if( all(c('neonatal_mean_raked','under5_mean_raked') %in% names(rasters_annual)) ){
    rasters_annual$nn_u5_ratio <- (
      rasters_annual$neonatal_mean_raked / rasters_annual$under5_mean_raked
    )
  }
  if( all(c('infant_mean_raked','under5_mean_raked') %in% names(rasters_annual)) ){
    rasters_annual$inf_u5_ratio <- (
      rasters_annual$infant_mean_raked / rasters_annual$under5_mean_raked
    )
  }
  # Get annualized rate of change rasters
  message("Making annualized rate-of-change rasters...")
  rasters_aroc <- lapply(rasters_annual, make_raster_aroc)
  names(rasters_aroc) <- paste0(names(rasters_aroc),'_aroc')
  # Return full dataset
  out_list <- list(
    'annual' = rasters_annual,
    'aroc'   = rasters_aroc
  )
  message("Done pulling rasters.")
  message(paste(
    "Pulled", length(rasters_annual), "annual rasters with",
    dim(rasters_annual[[1]])[3], "years each and", length(rasters_aroc),
    "AROC rasters."
  ))
  return(out_list)
}


## pull_raster_masks() -------------------------------------------------------->
#' 
#' @title Pull Raster Masks
#' 
#' @description Pulls lake data and the population mask for MBG stage 2 models
#' 
#' @param shp_version Which shapefile version should be used to subset to
#'   Stage 2?
#' 
#' @return Named list with two RasterLayer items, 'lakes' and 'pop'
#' 
pull_raster_masks <- function(shp_version){
  # Pull data from specified file paths as RasterLayers
  mask_folder <- '<<<< FILEPATH REDACTED >>>>'
  lakes <- raster( paste0(mask_folder,'<<<< FILEPATH REDACTED >>>>') )
  pop   <- raster( paste0(mask_folder,'<<<< FILEPATH REDACTED >>>>') )
  # Make the Stage 2 mask from a fasterized shapefile
  sf_adm0 <- sf::st_read(
    dsn = get_admin_shapefile(admin_level=0, version=shp_version)
  )
  lookup_table <- load_adm0_lookup_table()
  stage2 <- lookup_table[ Stage %in% c('1', '2a', '2b') & iso3!='chn', ]
  adm_field <- ifelse(
    detect_adm_shapefile_date_type(
      shpfile_path=get_admin_shapefile(admin_level=0, version=shp_version)
    )$shpfile_type=='gaul',
    'GAUL_CODE',
    'gadm_geoid'
  )
  setnames(stage2, adm_field, 'ADM0_CODE')
  sf_stage2 <- sf_adm0[ sf_adm0[['ADM0_CODE']] %in% stage2[,ADM0_CODE], ]
  sf_stage2$dummy_var <- 1
  mask_for_lakes <- fasterize::fasterize(
    sf     = sf_stage2,
    raster = lakes,
    field  = 'dummy_var',
    fun    = 'first'
  )
  mask_for_pop <- fasterize::fasterize(
    sf     = sf_stage2,
    raster = pop,
    field  = 'dummy_var',
    fun    = 'first'
  )
  # Subset to stage2 bounds
  lakes <- lakes * mask_for_lakes
  pop   <- pop   * mask_for_pop
  # Convert to data.tables with 'lat' and 'long' columns
  lakes_dt <- as.data.frame(lakes, xy=TRUE, na.rm=TRUE) %>% as.data.table %>% 
    setnames(., c('long','lat','x'))
  pop_dt   <- as.data.frame(pop,   xy=TRUE, na.rm=TRUE) %>% as.data.table %>% 
    setnames(., c('long','lat','x'))
  out_list <- list(
    lakes    = lakes,
    pop      = pop,
    lakes_dt = lakes_dt,
    pop_dt   = pop_dt
  )
  return(out_list)
}



## prep_u5m_data_for_mapping() ------------------------------------------------>
#' 
#' @title Prep all U5M run data for mapping 
#' 
#' @description Pull all spatial data needed for mapping in U5M
#' 
#' @param run_date U5M model run date
#' @param core_repo Path to the lbd_core repository, for sourcing shapefile funcs
#' @param u5m_repo Path to the U5M repository, used to source mapping utilities
#' @param shp_version [default '2018_08_28'] What version of the admin
#'   shapefiles should be pulled and joined to the admin data? Default is the 
#'   most recent GAUL shapefile.
#' @param start_year [default 2000] first year of model analysis
#' @param end_year [default 2017] final year of model analysis
#' @param resume [default TRUE] Searches the outputs folder for a prepped object
#'   containing the outputs of this function and loads if it has already been
#    created
#' @param cache_results [default TRUE] Save a copy of the prepped dataset in the
#'   run directory
#' @param simp_prop [default NULL] What proportion of points should be kept by
#'   the simplification algorithm? NOTE: SET THIS TO NULL FOR NOW -- I HAVE NOT
#'   FOUND A GOOD SOLUTION TO THE POLY SIMPLIFICATION ISSUE YET
#' @param include_non_gbd [default TRUE] Should countries not included in GBD
#'   (ESH Western Sahara and GUF French Guiana) be included in the map?
#'   
#' @return Named list containing the following:
#'   - adm*_dt:  Wide prepped data for each admin* level
#'   - adm*_shp: All adm* data joined with its fortified shapefile
#'   - rast:     All outputs of `pull_prep_raster_data`, including:
#'     - annual: RasterBricks containing annual estimates of model outputs
#'     - aroc: Annualized rate-of-change rasters for each of the annual bricks
#'   - masks: List containing two RasterLayer masks, 'lakes', and 'pop'
#' 
prep_u5m_data_for_mapping <- function(
  run_date,
  core_repo,
  u5m_repo,
  shp_version     = '<<<< REDACTED >>>>', 
  start_year      = 2000,
  end_year        = 2017,
  resume          = TRUE,
  cache_results   = TRUE,
  simp_prop       = NULL,
  prep_mex_bra    = TRUE,
  include_non_gbd = TRUE
  ){
  # Source mapping functions
  source(paste0(core_repo,'<<<< FILEPATH REDACTED >>>>'))
  source(paste0(core_repo,'<<<< FILEPATH REDACTED >>>>'))
  source(paste0(u5m_repo,'<<<< FILEPATH REDACTED >>>>'))
  # Get the filepath of what would be the cached object
  cache_fp <- sprintf(
    '<<<< FILEPATH REDACTED >>>>',
    U5M_BASE_PATH, run_date
  )
  # Read and return the cached data, if it exists
  if(resume & file.exists(cache_fp)){
    message("Pulling cached source data from file...")
    mapping_data_full <- get(load(cache_fp))
    message("... source data pulled.")
    return(mapping_data_full)
  }

  ## ** The rest of the function assumes that the cached data was not pulled **
  # Iterate through all admin units
  adm_list <- as.list(U5M_ADMIN_LEVELS)
  names(adm_list) <- paste0('adm',U5M_ADMIN_LEVELS)
  # PULL ALL ADMINISTRATIVE DATA
  admin_dts <- lapply(
    adm_list,
    function(lev) pull_prep_admin_data(run_date=run_date, adm_level=lev)
  )
  # PULL ADMINISTRATIVE DATA FOR BRAZIL AND MEXICO
  if(prep_mex_bra==TRUE){
    mb_admin_dt <- pull_prep_bra_mex_data(shp_version=shp_version)
  } else {
    mb_admin_dt <- data.table(admin_code = integer(0))
  }
  # PULL ALL ADMINISTRATIVE SHAPEFILES AND FORTIFY
  admin_shps <- lapply(
    adm_list,
    function(lev) pull_format_admin_shp(
      adm_level       = lev,
      shp_version     = shp_version,
      simp_prop       = simp_prop,
      stage2_only     = TRUE,
      include_non_gbd = include_non_gbd
  ))
  bg_shp <- pull_format_admin_shp(
    adm_level   = 0,
    shp_version = shp_version,
    simp_prop   = simp_prop,
    stage2_only = FALSE,
    include_non_gbd = TRUE
  )
  # Pull disputed area boundaries
  disputed <- pull_format_disputed_shp(
    shp_version = shp_version,
    simp_prop   = simp_prop
  )
  # PULL ALL RASTER DATA
  rast <- pull_prep_raster_data(
    run_date   = run_date, 
    start_year = start_year,
    end_year   = end_year
  )
  masks <- pull_raster_masks(shp_version = shp_version)

  ## Combine all data into an output list, cache (if specified), and return
  mapping_data_full <- list(
    dt       = admin_dts,
    mb_dt    = mb_admin_dt,
    shp      = admin_shps,
    bg_shp   = bg_shp,
    disputed = disputed,
    rast     = rast,
    masks    = masks
  )
  if(cache_results) save(mapping_data_full, file=cache_fp)
  message("All results formatted and pulled.")
  return(mapping_data_full)
}


## autogen_map_title() -------------------------------------------------------->
#' 
#' @title Automatically generate a U5M map title
#' 
#' @param topic The same as the `topic` argument in `make_admin_map()`
#' @param geo_level The same as the `geo_level` argument to make_admin_map()
#' @param year Either a numeric year or 'aroc'
#' 
#' @return A title for the map
#' 
autogen_map_title <- function(topic, geo_level, year){
  ## The title format will be ""
  ## Example: "Mortality Rate among Children Under 5 by Admin2 Unit (Mean, 2000)"
  ## OR: "Annualized Rate of Change in Mortality Rate among Children Under 5 (Mean)"
  # Format sentence based on whether the 'year' is AROC or not
  if (topic=='nn_u5_ratio'){
    return(sprintf('Ratio of Neonatal to Under-5 Mortality in %s',year))
  }
  if (topic=='inf_u5_ratio'){
    return(sprintf('Ratio of Infant to Under-5 Mortality in %s',year))
  }
  if (year=='aroc'){
    title_template <- 'Annualized Rate of Change in %s among %s%s (%s)'
  } else {
    title_template <- paste0('%s among %s%s (%s, ',year,')')
  }
  # Determine phenomenon (mortality rate or death counts)
  if( grepl('deathcounts',topic) ){
    phenom <- 'Number of Deaths'
  } else {
    phenom <- 'Mortality Rate'
  }
  # Determine age group
  ag_dt <- data.table(
    keys = c('under5', 'infant', 'neonatal'),
    vals = c('Children Under 5', 'Infants', 'Neonates')
  )
  # Determine geographic level
  geo_level_list <- list(
    'adm0' = ' by Country',
    'adm1' = ' by Admin1 Unit',
    'adm2' = ' by Admin2 Unit',
    'rast' = ''
  )
  # Determine statistic
  stat_dt <- data.table(
    keys = c('mean', 'lower', 'upper'),
    vals = c('Mean', '95% CI Upper Limit', '95% CI Lower Limit')
  )
  ## Compose the title
  full_title <- sprintf(
    title_template,
    phenom,
    ag_dt[ str_detect(topic,keys), vals ],
    geo_level_list[[geo_level]],
    stat_dt[ str_detect(topic,keys), vals ]
  )
  return(full_title)
}

## detect_color_scale() ------------------------------------------------------->
#' 
#' @title Detect color scale for U5M map
#' 
#' @description Auto-generates color breakpoints, labels for those breakpoints,
#'   and a set of colors that will serve as inputs to scale_fill_gradientn()
#' 
#' @param topic The same as the `topic` argument to make_admin_map()
#' @param geo_level The same as the `geo_level` argument to make_admin_map()
#' @param year The same as the `year` argument to make_admin_map()
#' 
#' @return named list with four items: 'legend_title', col_breaks', 'col_labs',
#'   grad_vals'
#' 
detect_color_scale <- function( topic, geo_level, year ){

  if( year == 'aroc'){
    # This is annualized rate-of-change data
    legend_title <- 'Annualized\nrate of\ndecline'
    col_breaks <- c(-0.1, -0.05, -0.02, 0.0)
    col_labs <- c('10%','5%', '2%', 'No decline')
    grad_vals <- c(
      '#002080', '#003866', '#00504C', '#006833', '#008019', '#009900', 
      '#53B639', '#A6D373', '#FAF0AD', '#EDA259', '#E15505'
    )
  } else if ( str_detect(topic,'deathcounts') ){
    # This is death counts data
    legend_title <- 'Number\nof deaths'
    if(geo_level=='rast') legend_title <- 'Number\nof deaths\nper grid\ncell'
    # Determine scaling factor based on geo_level
    geo_level_scale_list <- list(
      'adm0' = c(5000, 50000, 200000, 500000),
      'adm1' = c(500,  5000,  20000,  50000 ),
      'adm2' = c(50,   500,   2000,   5000  ),
      'rast' = c(0.1,  1,     4,      10    )
    )
    geo_level_labs_list <- list(
      'adm0' = c('5k',   '50k', '200k', '>500k'),
      'adm1' = c('500',  '5k',  '20k',  '>50k' ),
      'adm2' = c('50',   '500', '2k',   '>5k'  ),
      'rast' = c('',     '1',   '4',    '>10'  )
    )
    if( !(geo_level %in% names(geo_level_scale_list)) ) stop("Improper geo_level.")
    col_breaks <- geo_level_scale_list[[ geo_level ]]
    col_labs   <- geo_level_labs_list[[ geo_level ]]
    # Only four break points for now
    grad_vals <- c("#f7f7c0","#ec5f5f","#6d2a80","#400040")
  } else if ( endsWith(topic, 'ratio') ){
    # This is Infant/U5M or NN/U5M ratio data
    # Color pallette is blue - orange
    grad_vals <- brewer.pal("BuPu", n=9)[2:9]
    if(startsWith(topic, 'nn')){
      legend_title <- 'Neonatal/\nUnder-5\nMortality\nRatio'
      col_breaks <- c(.25, .5, .75)
      col_labs <- c('<0.25','0.5','>0.75')
    } else {
      legend_title <- 'Infant/\nUnder-5\nMortality\nRatio'
      col_breaks <- c(.5, .7, .9)
      col_labs <- c('<0.5','0.7','>0.9')
    }
  } else {
    # This is mortality rate data
    legend_title <- 'Mortality rate\nper 1000\nlive births'
    # Breaks and gradients vary by age group
    if( str_detect(topic, 'under5') ){
      # Handle under-5 data
      col_breaks <- c(0, 5, 25, 100, 200)
      col_labs   <- c("", "5", "25", "100", ">200")
      grad_vals <- c(
        '#009999','#5aada1','#88c1a9','#b1d6b0','#d8eab8','#ffffbf',"#fff2b6",
        "#fee5ae","#fed8a5","#fecb9c","#fdbe93","#fdb18b","#fda482","#fc9779",
        "#fc8a70","#fc7d68","#fb705f","#fb6356","#fb564d","#fa4945","#fa3c3c",
        "#f53b3a","#f03938","#eb3837","#e63635","#e13533","#db3431","#d6322f",
        "#d1312e","#cc2f2c","#c72e2a","#c22d28","#bd2b26","#b82a25","#b32823",
        "#ae2721","#a8261f","#a3241d","#9e231c","#99211a","#942018"
      )
    } else if ( str_detect(topic, 'infant') ){
      # Handle infant data
      col_breaks <- c(0, 2.5, 12.5, 75, 150)
      col_labs   <- c("", "2.5", "12.5", "75", ">150")
      grad_vals <- c(
        '#009999','#5aada1','#88c1a9','#b1d6b0','#d8eab8','#ffffbf',"#fff7ba",
        "#ffefb5","#fee8af","#fee0aa","#fed8a5","#fed0a0","#fec89a","#fdc195",
        "#fdb990","#fdb18b","#fda985","#fda180","#fc9a7b","#fc9276","#fc8a70",
        "#fc826b","#fc7a66","#fb7361","#fb6b5b","#fb6356","#fb5b51","#fb534c",
        "#fa4c46","#fa4441","#fa3c3c","#f73b3b","#f33a3a","#f03938","#ec3837",
        "#e93736","#e63635","#e23534","#df3532","#db3431","#d83330","#d5322f",
        "#d1312e","#ce302c","#ca2f2b","#c72e2a","#c42d29","#c02c28","#bd2b26",
        "#b92a25","#b62924","#b32823","#af2722","#ac2720","#a8261f","#a5251e",
        "#a2241d","#9e231c","#9b221a","#972119","#942018"
      )
    } else {
      # Handle neonatal data
      col_breaks <- c(0, 2.5, 12.5, 50, 100)
      col_labs   <- c("", "2.5", "12.5", "50", ">100")
      grad_vals <- c(
        '#009999','#5aada1','#88c1a9','#b1d6b0','#d8eab8','#ffffbf',"#fff2b6",
        "#fee5ae","#fed8a5","#fecb9c","#fdbe93","#fdb18b","#fda482","#fc9779",
        "#fc8a70","#fc7d68","#fb705f","#fb6356","#fb564d","#fa4945","#fa3c3c",
        "#f53b3a","#f03938","#eb3837","#e63635","#e13533","#db3431","#d6322f",
        "#d1312e","#cc2f2c","#c72e2a","#c22d28","#bd2b26","#b82a25","#b32823",
        "#ae2721","#a8261f","#a3241d","#9e231c","#99211a","#942018"
      )
    }
  }
  ## COLLECT AND RETURN VALUES
  out_list <- list(
    'col_title'  = legend_title,
    'col_breaks' = col_breaks,
    'col_labs'   = col_labs,
    'grad_vals'  = grad_vals
  )
  return(out_list)
}

## make_publication_map() ----------------------------------------------------->
#' 
#' @title Make a map for publication
#' 
#' @description Given input administrative-level data, prepare a map
#' 
#' @param full_data The output of prep_u5m_data_for_mapping
#' @param topic Either the name of a field in the admin dataset, or the name of 
#'   raster brick or layer in the `rast` list item.
#' @param year Either a year from 2000 to 2017 or 'aroc'
#' @param geo_level [default 'adm2'] Determines geographic level of mapping. One 
#'   of 'adm0', 'adm1', 'adm2', or 'rast'
#' @param title [default 'detect'] Title of the plot
#' @param subtitle [default NULL] Subtitle of the plot
#' @param adm0_restriction [default NULL] which adm0 codes should be plotted? If
#'   NULL, the default, all countries will be plotted
#' @param mapping_bounds What bounds should the map be plotted to? Should be a 
#'   named vector with names 'W','E','N','S'.
#' @param col_scale_title [default 'detect'] Title of the color scale legend
#' @param col_breaks [default 'detect'] Breaks to show in the color scale 
#'   legend. The first and last of these will become the minimum and maximum 
#'   plotting value
#' @param col_labels [default 'detect'] Labels accompanying the color scale breaks
#' @param col_scale [default rev(brewer.pal("Spectral", n=11))] Color scale that
#'    will be used as the basis for the gradient in `scale_color_gradientn()`
#' @param col_transformation [default 'identity'] Use this argument to add any
#'    transformation to the color scale (for example, plotting on a log scale
#'    rather than linear). Valid options are any arguments to the `trans` arg
#'    in ggplot `scale_*` functions.
#' @param add_projection [default TRUE] Should the Mollweide projection be 
#'    applied to this map? If not, the map will be in the standard Mercator
#'    projection.
#' @param add_masks [default 'detect'] Should the population and lake masks be
#'    added? If 'detect', will add for raster data but not for admin data
#' @param point_size What size should raster cells be plotted on the map?
#' 
#' @return ggplot objection containing the prepped map. Ggplot arguments can
#'    still be added to this map using `fig <- fig + [...]`
#' 
make_publication_map <- function(
  full_data,
  topic,
  year,
  geo_level        = 'adm2',
  title            = 'detect',
  subtitle         = NULL,
  adm0_restriction = NULL,
  mapping_bounds   = c('S'=-35.5, 'N'=53, 'W'=-108, 'E'=157.5),
  col_scale_title  = 'detect',
  col_breaks       = 'detect',
  col_labels       = 'detect',
  col_scale        = 'detect',
  col_transformation = 'identity',
  add_projection   = TRUE,
  add_masks        = 'detect',
  add_disputed     = TRUE,
  add_bra_mex      = TRUE,
  point_size       = 0.08
  ){
  ## Auto-detect title and color settings, then apply as needed
  default_title       <- autogen_map_title(topic, geo_level, year)
  default_color_scale <- detect_color_scale(topic, geo_level, year)
  mask_default        <- (geo_level == 'rast')
  # Apply to defaults only
  if(title=='detect') title <- default_title
  if(col_scale_title=='detect') col_scale_title <- default_color_scale$col_title
  if(col_breaks=='detect') col_breaks <- default_color_scale$col_breaks
  if(col_labels=='detect') col_labels <- default_color_scale$col_labs
  if(col_scale=='detect') col_scale <- default_color_scale$grad_vals
  if(add_masks=='detect') add_masks <- mask_default
  # Set upper and lower limits for color scale
  col_min <- min(col_breaks)
  col_max <- max(col_breaks)
  # Load background borders
  ad0_sf <- suppressWarnings(suppressMessages(
    sf::st_read(get_admin_shapefile(admin_level=0))
  ))

  ## Pull data for administrative units and raster-level data differently
  if (startsWith(geo_level, 'adm')){
    ## Format spatial data for administrative units
    topic_field <- paste0(topic,'_',year)
    in_data <- full_data$dt[[geo_level]][, 
                    c('admin_code', topic_field), with=F]
    setnames(in_data, topic_field, 'map_var')
    if (!str_detect(topic,'deathcounts') & (year!='aroc') & !str_detect(topic,'ratio')){
      in_data[, map_var := map_var * 1000 ]
    }
    # If there is an admin0 code restriction, subset the shapefile
    shp_for_merge <- full_data$shp[[geo_level]]
    if(!is.null(adm0_restriction)){
      shp_for_merge <- shp_for_merge[ADM0_CODE %in% adm0_restriction]
    }
    # Merge the input data with the underlying fortified shapefile
    map_data <- merge(
      x  = shp_for_merge,
      y  = in_data,
      by = 'admin_code'
    )
    map_data <- map_data[order(shp_order)]
    # Truncate the data to the limits of the color scale
    map_data[ map_var < col_min, map_var := col_min ]
    map_data[ map_var > col_max, map_var := col_max ]
  } else {
    ## Format spatial data for raster layers
    if(year=='aroc'){
      # Pull the AROC data for the topic as a rasterLayer
      map_rast <- full_data$rast$aroc[[paste0(topic, "_aroc")]]
    } else {
      # Pull the year of data as a rasterLayer
      map_rast <- full_data$rast$annual[[topic]][[ as.numeric(year) - 1999 ]]
    }
    # If an admin0 restriction exists, mask the raster
    if(!is.null(adm0_restriction)){
      ad0_sf_sub <- ad0_sf[ad0_sf$ADM0_CODE %in% adm0_restriction,]
      ad0_sf_sub$dummy <- 1
      country_mask <- fasterize::fasterize(
        sf = ad0_sf_sub,
        raster = map_rast,
        field = 'dummy'
      )
      map_rast <- map_rast * country_mask
    }
    # Convert the raster into a data.table for ggplotting
    map_data <- as.data.frame(map_rast, xy=TRUE, na.rm=TRUE) %>% as.data.table
    names(map_data) <- c('long', 'lat', 'mapvar')
    if (!str_detect(topic,'deathcounts') & (year!='aroc') & !str_detect(topic,'ratio')){
      map_data[, mapvar := mapvar * 1000 ]
    }
    map_data[ mapvar < col_min, mapvar := col_min ]
    map_data[ mapvar > col_max, mapvar := col_max ]
    # Subset out any data outside of the stage 2 boundaries
    map_data <- map_data[ (lat  < S2_BOUNDS$lat_max ) &
                          (lat  > S2_BOUNDS$lat_min ) &
                          (long < S2_BOUNDS$long_max) &
                          (long > S2_BOUNDS$long_min),
    ]
  }

  ## IF THE DATA EXISTS FOR BRAZIL AND MEXICO, INCLUDE IT
  topic_field <- paste0(topic,'_',year)
  map_mb <- (
    ( add_bra_mex==TRUE ) & 
    ( topic_field %in% names(full_data$mb_dt) ) &
    ( !str_detect(topic,'deathcounts') | (geo_level=='adm1') )  # Death counts at adm1 only
  )
  if( map_mb == TRUE ){
    mb_data <- full_data$mb_dt[, c('admin_code',topic_field), with=FALSE]
    setnames(mb_data, topic_field, 'map_var')
    if (!str_detect(topic,'deathcounts') & (year!='aroc') & !str_detect(topic,'ratio')){
      mb_data[, map_var := map_var * 1000 ]
    }
    # Merge the input data with the underlying fortified shapefile
    mb_data <- merge(
      x  = full_data$shp[['adm1']],
      y  = mb_data,
      by = 'admin_code'
    )
    mb_data <- mb_data[order(shp_order)]
    # Truncate the data to the limits of the color scale
    mb_data[ map_var < col_min, map_var := col_min ]
    mb_data[ map_var > col_max, map_var := col_max ]
  } else {
    mb_data <- data.table()
  }

  ## CREATE PLOT
  # Add all parts of the plot that are included for both rasters and polygons
  point_stroke_size <- 0
  na_color <- '#E6E6E6'
  fig <- ggplot() + 
    labs(
      title    = title,
      subtitle = subtitle,
      fill     = col_scale_title,
      color    = col_scale_title
    ) +
    theme_map() + 
    theme(
      legend.position = c(0.95, 0.80),
      legend.title = element_text(size=9, color='#222222'),
      legend.text  = element_text(size=6, color='#222222')
    )
  # Add different layers for rasters versus polygons
  if(geo_level=='rast'){
    # Plotting behavior for rasters
    fig <- fig +
      geom_point(
        data   = map_data,
        aes(x=long, y=lat, color=mapvar),
        shape  = 15,
        size   = point_size,
        stroke = point_stroke_size
      ) +
      scale_colour_gradientn(
        limits = c(col_min, col_max),
        colors = col_scale,
        breaks = col_breaks,
        labels = col_labels,
        trans  = col_transformation,
        na.value   = na_color
      )
  } else {
    # Plotting behavior for polygons
    fig <- fig + 
      geom_polygon(
        data = map_data,
        aes(x=long, y=lat, group=group, fill=map_var)
      ) +
      scale_fill_gradientn(
        limits = c(col_min, col_max),
        colors = col_scale,
        breaks = col_breaks,
        labels = col_labels,
        trans  = col_transformation,
        na.value   = na_color
      )
  }
  # Add Brazil and Mexico data, if desired
  if(map_mb==TRUE){
    fig <- fig + 
      geom_polygon(
        data = mb_data,
        aes(x=long, y=lat, group=group, fill=map_var)
      )
    if(geo_level=='rast'){
      fig <- fig + 
        scale_fill_gradientn(
          limits = c(col_min, col_max),
          colors = col_scale,
          breaks = col_breaks,
          labels = col_labels,
          trans  = col_transformation,
          na.value   = na_color
        )
    }
  }

  # Add masks, if desired
  if(add_masks==TRUE){
    fig <- fig +
      geom_point(
        data   = full_data$masks$lakes_dt,
        aes(x=long, y=lat),
        color  = '#D8ECF3',
        shape  = 15,
        size   = point_size,
        stroke = point_stroke_size
      ) + 
      geom_point(
        data   = full_data$masks$pop_dt,
        aes(x=long, y=lat),
        color  = na_color,
        shape  = 15,
        size   = point_size,
        stroke = point_stroke_size
      )
  }
  # Add country borders
  fig <- fig + geom_sf(
    data = ad0_sf,
    fill=NA,
    color='#222222',
    lwd=.2,
    size=.2
  )
  # Add disputed territories as dashed borders, if specified
  if(add_disputed==TRUE){
    fig <- fig + geom_polygon(
      data = full_data$disputed,
      aes(x=long, y=lat, group=group),
      fill=NA,
      linetype=3,
      color='#222222',
      lwd=.1
    )
  }
  # Add projections to the plot and crop to the smaller dataset
  if(add_projection==TRUE){
    fig <- fig + coord_sf(
      crs = 'mollweide',
      xlim = c(mapping_bounds['W'], mapping_bounds['E']),
      ylim = c(mapping_bounds['S'], mapping_bounds['N'])
    )
  } else {
    fig <- fig + coord_sf(
      xlim = c(mapping_bounds['W'], mapping_bounds['E']),
      ylim = c(mapping_bounds['S'], mapping_bounds['N'])
    )       
  }
  # Return the figure
  return(fig)
}
