### Nudge points with NA covariate extractions to nearest grid cell with covariate data, if within a threshold distance (in meters).
nudge_coordinates <- function(temp_data = the_data, the_covs = the_covs, covs_cs_df = covs_cs_df, threshold = 10000, save_fixed_coordinate = F) {
  temp_data$true_lat <- temp_data$latitude
  temp_data$true_long <- temp_data$longitude

  message("Nudging points with NA covariate extractions to nearest grid cell with covariate data, if within a threshold distance...")

  the_covs <- the_covs[(grep("gaul_code", the_covs, invert = TRUE))]
  for (a in 1:length(the_covs)) {
    current <- the_covs[[a]]
    has_na_data <- which(is.na(temp_data[, get(current)]))
    
    message(paste0(current, ": ", length(has_na_data), " rows with NA covariate data"))
    
    if (length(has_na_data) == 0) {
      next
    }

     Identify nearest data point with covariate data, and temporarily use those coordinates
    template <- all_cov_layers[[(as.character(current))]][[1]]

    for (b in 1:length(has_na_data)) {
      message(b)
      dist <- replace(distanceFromPoints(template, temp_data[has_na_data[[b]], c("longitude", "latitude")]), is.na(template), NA)
      dist[dist > threshold] <- NA
      if (cellStats(dist, "sum") == 0) {
        message("No cell found within threshold distance. Moving to next data point.")
      } else {
        coords <- as.data.table(xyFromCell(template, which.min(dist)[1]))
        suppressWarnings(temp_data[has_na_data[b], "longitude"] <- coords$x)
        suppressWarnings(temp_data[has_na_data[b], "latitude"] <- coords$y)
      }
    }

    # Extract data for all data rows
    cs_covs_new <- extract_covariates(temp_data,
      all_cov_layers[[current]],
      id_col = "a_rowid",
      return_only_results = TRUE,
      centre_scale = TRUE,
      period_var = "year",
      period_map = period_map
    )

    setnames(cs_covs_new$covs, current, "var")
    temp_data[, (current) := cs_covs_new$covs$var]
    covs_cs_df[covs_cs_df$name == (current), ] <- cs_covs_new[[2]][cs_covs_new[[2]]$name == (current), ]
    temp_data <- as.data.table(temp_data)
    temp_data$longitude <- temp_data$true_long
    temp_data$latitude <- temp_data$true_lat
  }
  
  # Return results
  if (save_fixed_coordinate == T) {
    return(list("temp_data" = temp_data, "covs_cs_df" = covs_cs_df, "was_na_data" = was_na_data))
  } else {
    return(list("temp_data" = temp_data, "covs_cs_df" = covs_cs_df))
  }
}
