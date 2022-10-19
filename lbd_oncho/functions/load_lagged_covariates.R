load_lagged_covariates <- function(covariate_config, template, raster_agg = raster_agg_factor, start_year = min(year_list), end_year = max(year_list)) {
  raster_list <- list()
  for (i in 1:nrow(covariate_config)) {
    cov_raster <- load_worldpop_covariate(template_raster = template,
                                          covariate = covariate_config[i, covariate],
                                          pop_measure = covariate_config[i, measure],
                                          pop_release = covariate_config[i, release],
                                          start_year = start_year - covariate_config[i, year_lag],
                                          end_year = end_year - covariate_config[i, year_lag],
                                          interval = as.integer(interval_mo))
    
    if (as.integer(raster_agg) != 1) { ## Aggregate covariate raster by raster_agg factor. Using temp raster to accommodate resampling of categorical covariates.
      temp_raster <- suppressWarnings(raster::aggregate(cov_raster[[1]], as.integer(raster_agg))) ## Create temp raster for resampling.
      if ((dataType(cov_raster[[1]]) == "LOG1S") | grepl("INT", dataType(cov_raster[[1]])) | (isTRUE(all.equal(sort(unique(as.data.table(raster::freq(cov_raster[[1]], useNA = "no", merge = TRUE, digits = 1)))$value), c(0, 1))))) { # Categorical covariates
        cov_raster <- raster::resample(cov_raster[[1]], temp_raster, method = "ngb")
      } else if (names(cov_raster) %in% c("worldpop")) { # Covariates requiring sum rather than mean during aggregation
        values(cov_raster[[1]])[is.na(values(cov_raster[[1]]))] <- 0 # raster::aggregate introduces NaN values for worldpop, for an unknown reason; this line avoids the problem
        cov_raster[[1]] <- suppressWarnings(raster::aggregate(cov_raster[[1]], fun=sum, as.integer(raster_agg)))
      } else { # Continuous covariates
        cov_raster[[1]] <- raster::resample(cov_raster[[1]], temp_raster)
      }
    }
    
    raster_list <- c(raster_list, cov_raster)
  }
  return(raster_list)
}