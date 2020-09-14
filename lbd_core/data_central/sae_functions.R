################################################################################
# This function is used to fractionally aggregate any MBG standard covariate raster
# to areas defined by a custom shapefile.
#
################################################################################

#' Function used to fractionally aggregates a rasterized covariate to custom areal units (used in SAE modeling)
#' @title frac_agg_covs
#' @description: This function builds a custom link table defined by the custom shapefile (or uses predefined link table);
#' load in cropped mbg covariates to the shapefile;
#' Load in population raster also cropped to the shapefile;
#' and fractionally aggregating based on the custom link table and an aggregation method.
#' This function uses frac_agg_covs_raster for the final aggregation step
#'
#'
#' @param covs MBG standard covariates you wish to aggregate to shapefile areas
#' @param years Years you wish to aggregate covariates
#' @param measures Statistical measures of aggregated covariates (mean, median, sum). Must be same length as covariates
#' @param releases the release date of each covariate e.g., "2019_06_10". Must be same length as covariates.
#' @param shapefile_path The filepath to the custom shapefile
#' @param shapefile_field The field in the shapefile that indicates the areal units
#' @param core_repo LBD core repo to source central functions
#' @param agg_method Method of aggregation, default to population-weighted mean
#' @param worldpop_age_sex_group Age/sex group for the worldpop raster (i.e. "a1549t" for 15-49 both sexes), needed if agg_method = "pop_weight", default is total pop
#' @param cores How many cores should be used in parallel processing
#' @param link_table Link table if you want to provide your own, or else it will build a custom link table from the shapefile
#'
#' @return Returns data table with shapefile field, year, and weighted covariate value for each covoaraite as its own column

frac_agg_covs <- function(covs, years, measures, releases, shapefile_path,
                          shapefile_field,
                          core_repo = "<<<< FILEPATH REDACTED >>>>",
                          agg_method = "pop_weight",
                          worldpop_age_sex_group = "total",
                          worldpop_age_sex_release = "2017_04_27",
                          cores = 1,
                          link_table = NULL) {

  # Load necessary packages
  library(tidyr)
  library(purrr)
  library(dplyr)
  library(data.table)
  library(rgdal)
  library(rgeos)
  library(raster)
  library(parallel)

  # Ensure no masking of these functions
  rename <- dplyr::rename
  select <- dplyr::select
  summarize <- dplyr::summarize

  # Ensure the aggregation method is in pop_weight, sum or unweighted mean
  if (!agg_method %in% c("pop_weight", "sum", "mean")) {
    stop("Aggregation method must be either population weighted (pop_weight), sum (sum), or an unweighted mean (mean)")
  }

  # Make sure measures is same length as covs
  if (length(covs) != length(measures)) {
    stop("Statistical measure must be passed to this function for each covariate selected, lengths do not match")
  }

  if (length(covs) != length(releases)) {
    stop("Release date must be passed to this function for each covariate measure selected, lengths do not match")
  }

  # Read-in or build custom link table
  if (is.null(link_table)) {
    message("Building custom link table")
    source(paste0(core_repo, "mbg_central/fractional_raking_functions.R"))
    source(paste0(core_repo, "mbg_central/prep_functions.R"))
    link_table <- build_link_table(shapefile_version = NULL,
                                   cores = cores,
                                   region = NULL,
                                   custom_shapefile_path = shapefile_path,
                                   custom_shapefile_field = shapefile_field)$link_table
  }

  # Build simple raster
  message("Building simple raster")
  subset_shape  <- load_simple_polygon(gaul_list = NULL, buffer = 0.4, subset_only = FALSE, custom_shapefile_path = shapefile_path)[[1]]
  simple_raster <- build_simple_raster_pop(subset_shape, field = shapefile_field, link_table = NULL)[['simple_raster']]

  # Load covariate rasters
  message("Loading covariates and population raster")
  source(paste0(core_repo, "mbg_central/covariate_functions.R"))

  covariate_config <- data.table(covariate = covs,
                                 measure = measures,
                                 release = releases)
  loader <- MbgStandardCovariateLoader$new(start_year = min(years),
                                           end_year = max(years),
                                           interval = 12,
                                           covariate_config = covariate_config)
  cov_rasters <- loader$get_covariates(simple_raster)


  # If aggregating population, rename to be measures
  if (all(covs == "worldpop")) {
    covs <- measures
  }

  # load population raster if population weighting
  if (agg_method == "pop_weight") {
    pop_raster <- load_worldpop_covariate(template_raster = simple_raster,
                                          measure = worldpop_age_sex_group,
                                          release = worldpop_age_sex_release,
                                          start_year = min(years),
                                          end_year = max(years),
                                          interval = 12)[[1]]
  }

  # Fractionally aggregate covariates and return data table
  message("Fractionally aggregating covariates")
  aggregated_covs <-
    mclapply(1:length(covs), mc.cores = cores, function(i) {
      frac_agg_cov_raster(cov_name = covs[i],
                          cov_raster = cov_rasters[[i]],
                          years = years,
                          link_table = link_table,
                          shapefile_field = shapefile_field,
                          agg_method = agg_method,
                          pop_raster = if (agg_method == "pop_weight") pop_raster else NULL)}) %>%
    reduce(full_join, by = c(shapefile_field, "year")) %>%
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
  } else if (nlayers(cov_raster == 1)) {
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

  # Document what % of pixels found in link table are missing in covariate raster,
  # This is a rough check
  cov_values <-
    cov_values %>%
    drop_na(cov)

  pct_missing_pixels <-
    setdiff(unique(link_table$ID), unique(cov_values$pixel_id)) %>%
    length() %>%
    {100*(. / length(unique(link_table$pixel_id)))} %>% round(1)

  message(paste0(pct_missing_pixels, "% of pixels found in link table missing ", cov_name, " data"))

  if (agg_method == "pop_weight") {
    # Assemble weighted covariate raster by pixel_id if population weighting
    pop_values <-
      rbindlist(lapply(1:nlayers(pop_raster), function(i) {
        data.table(pixel_id = 1:length(values(pop_raster[[i]])),
                   pop = values(pop_raster[[i]]),
                   year = years[i])
      }))

    # Document what % of pixels found in link table are missing in weighted raster
    pop_values <-
      pop_values %>% drop_na(pop)

    ct_missing_pixels <-
      setdiff(unique(link_table$ID), unique(pop_values$pixel_id)) %>%
      length() %>%
      {100*(. / length(unique(link_table$pixel_id)))} %>% round(1)

    message(paste0(pct_missing_pixels, "% of pixels found in link table missing population data"))

    # Merge covariate values by pixel id and year with weighted raster
    cov_values <- full_join(cov_values, pop_values, by = c("pixel_id", "year"))
    if (any(is.na(cov_values))) {
      message("Unexpected missingness when joining covariates and population raster")
    }
  }

  # Make sure each pixel has estimates for every year
  complete_pixel_years <-
    cov_values %>%
    count(pixel_id) %>%
    filter(n != length(years)) %>%
    nrow()
  if (complete_pixel_years != 0) {
    stop ("Some pixels are missing covariate values in certain years, check on covariate raster data")
  }

  # Return fractionally aggregated covariates by field using link table --------------------
  cov_aggregated <-
    link_table %>%
    dplyr::rename(field = shapefile_field) %>%
    dplyr::select(pixel_id, field, area_fraction) %>%
    right_join(cov_values, by = "pixel_id") %>%
    group_by(field, year)

  # Choose aggregation method
  if (agg_method == "pop_weight") {
    cov_aggregated <-
      cov_aggregated %>%
      summarize(covariate = weighted.mean(cov, pop*area_fraction))
  } else if (agg_method == "sum") {
    cov_aggregated <-
      cov_aggregated %>%
      summarize(covariate = sum(cov*area_fraction))
  } else if (agg_method == "mean") {
    cov_aggregated <-
      cov_aggregated %>%
      summarize(covariate = mean(cov*area_fraction))
  } else {
    stop("Aggregation method not population weighted, sum, or mean")
  }

  # Return collapsed data table
  cov_aggregated <-
    cov_aggregated %>%
    ungroup() %>%
    data.table() %>%
    setnames(c("covariate", "field"), c(cov_name, shapefile_field))

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
  adjmat <- as(spdep::nb2mat(neighbors, style="B"), "dgTMatrix")
}
