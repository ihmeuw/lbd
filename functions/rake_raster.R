# Aniruddha Deshpande, adesh@uw.edu
rake_raster <- function(
                 showProgress = FALSE,
                 subset,
                 weights,
                 values,
                 gbd_table,
                 output_rf = FALSE
  ) {
  counter <- 0
  rf_list <- list()
  brick_list <- list()
  for (j in 1:nlayers(values)) {
  rf_yr_list <- list()
  ras_list <- list()  
    for (i in unique(subset)) {
      counter <- counter + 1
      if (showProgress) {
        print(paste("loc_id", i, "pct", counter/(length(unique(subset)) * nlayers(values))))
      }
      # Create a subsetRaster raster which will be used to subset to geography
      subsetRaster <- subset
      subsetRaster[subsetRaster != i] <- NA
      subsetRaster[!is.na(subsetRaster)] <- 1
      
      # Create a weights raster for subset geography
      subsetRaster_weights <- subsetRaster*weights[[j]]
      
      # Calculate Sum(value * weightsulation)
      lbd_ras <- values[[j]]*subsetRaster
      lbd_num <- lbd_ras*subsetRaster_weights
      lbd_num <- sum(getValues(lbd_num), na.rm = TRUE)

      # Calculate total weightsulation in the subset geography
      lbd_denom <- sum(getValues(subsetRaster_weights), na.rm = TRUE)
      
      # Calculate weighted mean
      lbd_agg <- lbd_num/lbd_denom
      
      # Grab GBD value & calculate raking factor
      gbd_values <- as.numeric(filter(gbd_table, year_id == j+1999 &
                         location_id == i)$val)

      if (length(gbd_values)== 0) {
          rf = 0
          lbd_ras <- lbd_ras * rf 
        } else {
          if (lbd_denom == 0) {
            rf <- NA
            lbd_ras[!is.na(lbd_ras)] <- gbd_values/ncell(lbd_ras[!is.na(lbd_ras)])
          } else {
            rf <- gbd_values/lbd_agg
            lbd_ras <- rf * lbd_ras
          }
      }

      if (output_rf) {
        rf_ras <- lbd_ras
        rf_ras[!is.na(rf_ras)] <- rf
        rf_yr_list[[length(rf_yr_list) + 1]] <- rf_ras  
      }
      ras_list[[length(ras_list) + 1]] <- lbd_ras
    }

    if (output_rf) {
      rf_list[[j]] <- do.call(raster::merge, rf_yr_list)  
    } else {
      rf_list <- NA
    }
    brick_list[[j]] <- do.call(raster::merge, ras_list)
  }
  raked_brick <- do.call(stack, brick_list)
  output_list <- list(raked_brick, stack(rf_list))
  names(output_list) <- c("raked_brick", "rf_list")
  return(output_list)
}

