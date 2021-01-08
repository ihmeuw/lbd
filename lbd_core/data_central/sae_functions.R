################################################################################
# This function is used to fractionally aggregate any covariate raster (MBG standard or custom)
# to areas defined by any shapefile. 
################################################################################

#' Function used to fractionally aggregates a rasterized covariate to custom areal units (used in SAE modeling)
#' @title frac_agg_covs
#' @description: This function builds a custom link table defined by the custom shapefile (or uses predefined link table);
#' load in cropped mbg covariates or a custom covariate to the shapefile;
#' and fractionally aggregating based on the custom link table and an aggregation method. 
#' This function uses frac_agg_covs_raster for the final aggregation step
#' 
#'
#' @param cov_config Covariate config that specifies four things: The MBG standard covariates names (list found here: <<<< FILEPATH REDACTED >>>>); 
#'                   the statistical measure of the covariate (usually mean or median, also found in link above); the agg_method (aggregation method, for now limited to sum, mean or population weighted mean); 
#'                   and the covariate release date e.g., "2019_06_10" (if this is an empty or missing column it will substitute most current covariate release). 
#'                   If using a custom covariate, this need only contain the covariate (what you want to call the covariate) and the agg_method
#' @param years Years you wish to aggregate covariates
#' @param shapefile_path The filepath to the shapefile
#' @param shapefile_field The field in the shapefile that indicates the areal units
#' @param core_repo LBD core repo to source central functions
#' @param worldpop_age_sex_group Age/sex group for the worldpop raster (i.e. "a1549t" for 15-49 both sexes), needed if agg_method = "pop_weight", default is total pop
#' @param worldpop_age_sex_release If agg_method = "pop_weight", defines which release of worldpop_age_sex_group should be used 
#' @param cores How many cores should be used in parallel processing
#' @param link_table You can pass a Link table if you want to provide your own, or else it will build a custom link table from the shapefile using build_link_table()
#' @param id_raster Cropped raster used to build link table, it is an output of build_link_table()
#' @param custom_covariate Covariate raster that is the same number of layers as years that be used instead of mbg standard covariates 
#' 
#' @return Returns data table with shapefile field, year, and weighted covariate value for each covariate as columns
#' 
#' @examples
#' Aggregating total population in Pakistan at the admin 1 level from worldpop
#' cov_config <- data.table(covariate = "worldpop", measure = "total", agg_method = "sum", release = "2017_04_27")
#' frac_agg_covs(cov_config, years = 2000:2017, shapefile_path = path_to_shapefile,
#'               shapefile_field = "ADM1_CODE")
#' 
#' Aggregating total population in Pakistan from custom population raster at the admin 1 level 
#' cov_config <- data.table(covariate = "custom_pop", agg_method = "sum")
#' frac_agg_covs(years = 2000:2017, shapefile_path = path_to_shapefile,
#'               shapefile_field = "ADM1_CODE", custom_covariate = custom_population_raster)
#' 

frac_agg_covs <- function(cov_config, 
                          years, 
                          shapefile_path,
                          shapefile_field, 
                          core_repo = paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/"),
                          worldpop_age_sex_group = "total",
                          worldpop_age_sex_release = "2020_03_20",
                          cores = 1, 
                          link_table = NULL,
                          id_raster = NULL,
                          custom_covariate = NULL) {
  
  # Load necessary packages 
  library(tidyr)
  library(purrr)
  library(dplyr)
  library(data.table)
  library(rgdal)
  library(assertthat)
  library(rgeos)
  library(raster)
  library(parallel)
  library(doParallel)

  # Load mbg functions
  source(paste0(core_repo, '/mbg_central/setup.R'))
  load_mbg_functions(core_repo)
  
  # Run argument check to fail early if incorrectly specified -----------------------------------------
  if (!is.numeric(years)) {
    stop("Years must be a numeric vector")
  }
  
  if (!all(map_lgl(c("covariate", "agg_method"), ~ . %in% names(cov_config)))) {
    stop("At minimum, the covariate config must contain the name of the covariate and agg_method")
  }
  
  # Ensure the aggregation method is in pop_weight, sum or unweighted mean 
  if (!all(map_lgl(cov_config$agg_method, ~ . %in% c("pop_weight", "sum", "mean")))) { 
    stop("Aggregation method must be either population weighted (pop_weight), sum (sum), or an unweighted mean (mean)")
  }
  
  # Additional checks if loading MBG standard covariates 
  if (is.null(custom_covariate)) {
    cov_config <- data.table(cov_config)
    
    # Update any missing "release" values with newest release; populates "release" column if not present
    update_fixed_effect_config_with_missing_release(cov_config)
    
    # Ensure that covariate config is correctly specified
    if (!all(map_lgl(c("covariate", "measure", "agg_method", "release"), ~ . %in% names(cov_config)))) {
      stop("Covariate config must contain the name of the covariate, the measure, agg_method, and release")
    }
  }
  
  # Load link table ---------------------------------------------------------------------------------
  # If no link table is provided, either load link table from VR filepath or create custom link table 
  message("Loading link table")
  if (is.null(link_table)) {
    link_filepaths <- vr_convert_link_filepath(shapefile_path)
    if (file.exists(link_filepaths$link_path)) {
      
      message("Grabbing related link table and id raster")
      link_table <- readRDS(link_filepaths$link_path)
      id_raster  <- readRDS(link_filepaths$id_path)
      
    } else {
      
      message("Building custom link table")
      link_table_list <- build_link_table(shapefile_version = NULL,
                                          cores = cores,
                                          region = NULL,
                                          custom_shapefile_path = shapefile_path,
                                          custom_shapefile_field = shapefile_field)
      link_table <- link_table_list$link_table
      id_raster  <- link_table_list$id_raster
    }
  }
  
  # Link table shapefile field must be numeric 
  link_table <- link_table %>% dplyr::mutate_at(vars(shapefile_field), funs(as.numeric(as.character(.)))) %>% data.table()
  
  # Ensure both link table and id_raster are provided, they are necessary for creating simple raster 
  if (!is.null(link_table)) {
    if (is.null(id_raster)) {
      stop("If you are providing a link table, you must also provide the id_raster which is an output of build_link_table()")
    }
  }
  
  # Create simple raster using link table and ID raster
  message("Creating simple raster")
  subset_shape <- st_read(shapefile_path, quiet = T) %>% dplyr::select(shapefile_field)
  simple_raster <- build_simple_raster_pop(subset_shape, 
                                           field = shapefile_field, 
                                           link_table = link_table, 
                                           id_raster = id_raster,
                                           pop_release = "2019_11_25")[['simple_raster']]
  
  # Grab covariates --------------------------------------------------------------------
  if (is.null(custom_covariate)) {
    message("Loading MBG standard covariates")
    # Load covariates from mbg central
    source(paste0(core_repo, "mbg_central/covariate_functions.R"))
    cov_rasters <- 
      mclapply(split(cov_config, 1:nrow(cov_config)), mc.cores = cores, function(cov) {
        loader <- MbgStandardCovariateLoader$new(start_year = min(years),
                                                 end_year = max(years),
                                                 interval = 12,
                                                 covariate_config = cov)
        cov_rasters <- loader$get_covariates(simple_raster)}) %>% unlist()
  } else {
    message("Cropping custom covariate")
    cov_rasters <- 
      mclapply(split(cov_config, 1:nrow(cov_config)), mc.cores = cores, function(cov) {
        custom_cov <- custom_covariate[[cov$covariate]]
        
        # The covariate raster must have the same number of layers as years
        if (nlayers(custom_cov) != length(years)) {
          stop("The custom covariate must have the same number of layers as years you wish to aggregate over")
        }
        if (!compareCRS(custom_cov, simple_raster)) {
          message("Rasters do not share same coordinate reference system, check in on this")
        }
        
        # Crop to simple raster 
        custom_cov <- raster::crop(custom_cov, raster::extent(simple_raster))
      })
  }

  # When aggregating population, measure is now the covariate name to differentiate between different age groupings 
  if (all(cov_config$covariate == "worldpop")) {
    cov_config$covariate <- cov_config$measure
  }
  
  # Name cov raster by covariate
  names(cov_rasters) <- cov_config$covariate
  
  # load population raster if population weighting for any covariate
  if ("pop_weight" %in% cov_config$agg_method) {
    pop_covariate_config <- 
      data.table(covariate = 'worldpop', 
                 measure = worldpop_age_sex_group,
                 release = worldpop_age_sex_release)
    pop_loader <- 
      MbgStandardCovariateLoader$new(start_year = min(years), 
                                     end_year = max(years), 
                                     interval = 12,
                                     covariate_config = pop_covariate_config)
    pop_raster <- pop_loader$get_covariates(simple_raster)$worldpop
  }
  
  # Fractionally aggregate covariates and return data table ----------------------------------
  message("Fractionally aggregating covariates")                        
  aggregated_covs <-
    mclapply(split(cov_config, 1:nrow(cov_config)), mc.cores = cores, function(cov) {
      frac_agg_cov_raster(cov_name = cov$covariate,
                          cov_raster = cov_rasters[[cov$covariate]],
                          years = years,
                          link_table = link_table,
                          shapefile_field = shapefile_field,
                          agg_method = cov$agg_method,
                          pop_raster = if (cov$agg_method == "pop_weight") pop_raster else NULL)}) %>% 
    purrr::reduce(full_join, by = c(shapefile_field, "year")) %>% 
    data.table()
}

#' Fractionally aggregates a cropped covariate raster given an aggregation method and link table
#' @title frac_agg_cov_raster
#' @description: This function is often used in conjunction with frac_agg_cov, and it
#' aggregates a covariate raster given a link table and aggregation method
#' frac_agg_covs function
#' 
#'
#' @param cov_name Name of the covariates
#' @param cov_raster raster of covariate, cropped to population raster if population weighting 
#' @param years years included in covariate raster
#' @param link_table Link table for fractional aggregation
#' @param shapfile_field The field in the shapefile that indicates the areal units
#' @param agg_method Aggregation method. Currently limited to population weight, sum, or mean (this can be expanded)
#' @param pop_raster Population raster if using populatoin weighted mean 
#' 
#' @return Returns data table with shapefile field, year, and weighted covariate value
#' 

frac_agg_cov_raster = function(cov_name,
                               cov_raster,
                               years,
                               link_table,
                               shapefile_field,
                               agg_method,
                               pop_raster) {
  
  message("Working on ", cov_name)
  # Check model inputs
  if (agg_method == "pop_weight") {
    if (extent(cov_raster) != extent(pop_raster)) stop("Covariate raster and population raster are not aligned")
  }
  
  # Assemble covariate values by pixel_id ---------------------------
  if (nlayers(cov_raster) == length(years)) {
    cov_values <-
      rbindlist(lapply(1:nlayers(cov_raster), function(i) {
        data.table(pixel_id = 1:length(values(cov_raster[[i]])),
                   cov = values(cov_raster[[i]]),
                   year = years[i])
      }))
  } else if (nlayers(cov_raster) == 1) {
    # If synoptic, repeat over all years specified by covariate
    message(paste(cov_name, "is a synoptic covariate, copying data for all years"))
    pixel_years <- 
      data.table(expand.grid(1:length(values(cov_raster)), years)) %>% 
      setnames(c("pixel_id", "year"))
    cov_values <- 
      data.table(pixel_id = 1:length(values(cov_raster)),
                 cov = values(cov_raster)) %>% 
      left_join(pixel_years, by = "pixel_id")
  } else {
    stop("Covariate raster does not have correct amount of layers representing each year")
  }
  
  # Document what % of pixels found in link table are missing in covariate raster
  cov_values <- tidyr::drop_na(cov_values, cov)
  
  pct_missing_pixels <- 
    setdiff(unique(link_table$ID), unique(cov_values$pixel_id)) %>% 
    length() %>% 
    {100*(. / length(unique(link_table$pixel_id)))} %>% round(1)
  
  message(paste0(pct_missing_pixels, "% of pixels found in link table missing ", cov_name, " data"))
  
  # Make sure each pixel has estimates for every year
  complete_pixel_years <- 
    cov_values %>% 
    dplyr::count(pixel_id) %>% 
    dplyr::filter(n != length(years)) %>% 
    nrow()
  
  if (complete_pixel_years != 0) {
    # Drop zero value and see if this helps find complete pixels
    cov_values <- filter(cov_values, cov != 0)
    complete_pixel_years <-
      cov_values %>%
      dplyr::count(pixel_id) %>%
      dplyr::filter(n != length(years)) %>%
      nrow()
    if (complete_pixel_years != 0) {
      incomplete_pixels <- 
        cov_values %>%
        dplyr::count(pixel_id) %>%
        dplyr::filter(n != length(years)) %>% 
        pull(pixel_id)
      
      message ("Some pixels are missing covariate values in certain years, check on covariate raster data")
    }
  }
  
  if (agg_method == "pop_weight") {
    # Assemble weighted covariate raster by pixel_id if population weighting
    pop_values <- 
      rbindlist(lapply(1:nlayers(pop_raster), function(i) {
        data.table(pixel_id = 1:length(values(pop_raster[[i]])),
                   pop = values(pop_raster[[i]]),
                   year = years[i])
      }))
    
    # Document what % of pixels found in link table are missing in weighted raster 
    pop_values <- drop_na(pop_values, pop)
    
    pct_missing_pixels <- 
      setdiff(unique(link_table$ID), unique(pop_values$pixel_id)) %>% 
      length() %>% 
      {100*(. / length(unique(link_table$pixel_id)))} %>% round(1)
    
    message(paste0(pct_missing_pixels, "% of pixels found in link table missing population data"))
    
    # Merge covariate values by pixel id and year with weighted raster
    cov_values <- dplyr::full_join(cov_values, pop_values, by = c("pixel_id", "year"))
    
    # See if there are pixels where we have population but are missing covariate 
    missing <- 
      cov_values %>% 
      dplyr::summarize(cov = round(100 * mean(is.na(cov)), 1),
                       pop = round(100 * mean(is.na(pop)), 1))
    
    if (missing$pop != 0 | missing$pop != 0) {
      message(paste0(missing$cov, "% of pixels found in link table missing covariate data where there is population data"))
      message(paste0(missing$pop, "% of pixels found in link table missing population data where there is covariate data"))
      message("Make sure this is not worrying")
    }
  }
  
  # Return fractionally aggregated covariates by field using link table --------------------
  cov_aggregated <- 
    link_table %>% 
    dplyr::rename(field = shapefile_field) %>% 
    dplyr::select(pixel_id, field, area_fraction) %>%
    inner_join(cov_values, by = "pixel_id") %>% # Join only places that are in the link table and in the covariate raster 
    group_by(field, year)
  
  # Choose aggregation method
  if (agg_method == "pop_weight") {
    cov_aggregated <- dplyr::summarize(cov_aggregated, covariate = weighted.mean(cov, pop*area_fraction, na.rm = T))
  } else if (agg_method == "sum") {
    cov_aggregated <- dplyr::summarize(cov_aggregated, covariate = sum(cov*area_fraction, na.rm = T))
  } else if (agg_method == "mean") {
    cov_aggregated <- dplyr::summarize(cov_aggregated, covariate = mean(cov*area_fraction, na.rm = T))
  } else {
    stop("Aggregation method not population weighted, sum, or mean")
  }
  
  # Return collapsed data table
  cov_aggregated <- data.table(ungroup(cov_aggregated))
  cov_aggregated <- setnames(cov_aggregated, c("covariate", "field"), c(cov_name, shapefile_field))
  
  # Make sure it is square, there should be values for each area in each year 
  assert_that(cov_aggregated %>% 
    dplyr::select(-cov_name) %>% 
    dplyr::summarise_all(funs(length(unique(.)))) %>%
    prod() == nrow(cov_aggregated))
    
  return(cov_aggregated)
}

#' @title create_adj_matrix
#' @description: This function creates an adjacency matrix given a shapefile path to be used in SAE modeling 
#' 
#' @param shapefile_path The filepath to the custom shapefile

#' 
#' @return Adjacency matrix used in small area estimation modeling, with neighbors defined by shared borders

create_adj_matrix <- function(shapefile_path) {
  
  library(sf)
  library(spdep)
  library(rgdal)
  
  shapefile <- rgdal::readOGR(shapefile_path)
  neighbors <- spdep::poly2nb(shapefile)
  adjmat <- as(spdep::nb2mat(neighbors, style = "B", zero.policy = T), "dgTMatrix")
}
