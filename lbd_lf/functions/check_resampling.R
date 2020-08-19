check_resampling <- function(check_resamp_rd, indicator_group = "lf", indicator_poly = "had_lf_poly", shp_path = <<<< FILEPATH REDACTED >>>>, ignore_warnings = FALSE, return_list = FALSE) {
  ### Checks for common problems that arise during polygon resampling.
  ### 1) NA location_code in data set
  ### 2) shapefile missing from shapefile directory
  ### 3) NA or missing referenced location_code in shapefiles
  ### 4) resampling results in a very low number of clusters

  ### setup ########################################################
  require(doParallel)
  require(plyr)
  require(snow)
  require(doSNOW)
  resamp_input <- fread(<<<< FILEPATH REDACTED >>>>)
  resamp_output <- fread(<<<< FILEPATH REDACTED >>>>)

  input_count <- table(resamp_input[, .(shapefile, location_code)]) %>%
    as.data.table() %>%
    .[N > 0]
  output_count <- table(resamp_output[, .(shapefile, location_code)]) %>%
    as.data.table() %>%
    .[N > 0]
  setnames(output_count, "N", "output_N")
  setnames(input_count, "N", "input_N")
  check_counts <- merge(input_count, output_count, by = c("shapefile", "location_code"), all = T)

  if (nrow(check_counts[output_N == input_N, ]) > 0) {
    message(paste0(
      "Warning! ", nrow(check_counts[output_N == input_N, ]), " of ", nrow(check_counts),
      " (", round(nrow(check_counts[output_N == input_N, ]) / nrow(check_counts), 3) * 100, "%) shapefile/location_code combos resampled to only 1 cluster"
    ))
  }

  if (nrow(check_counts[is.na(output_N) == T, ]) > 0) {
    message(paste0(
      "Warning! ", nrow(check_counts[is.na(output_N) == T, ]), " of ", nrow(check_counts),
      " (", round(nrow(check_counts[is.na(output_N) == T, ]) / nrow(check_counts), 3) * 100, "%) shapefile/location_code combos went unsampled. ",
      sum(check_counts[is.na(output_N) == T, input_N]), " polygons total."
    ))
  }

  print(check_counts[is.na(output_N) == T, ])

  if (nrow(resamp_input[is.na(location_code) == T | location_code == "NA", ]) > 0) {
    message(paste0(
      "Warning! ", nrow(resamp_input[is.na(location_code) == T | location_code == "NA", ]), " of ", nrow(check_counts) + length(resamp_input[is.na(location_code) == T | location_code == "NA", shapefile]),
      " (", round(nrow(resamp_input[is.na(location_code) == T | location_code == "NA", ]) / (nrow(check_counts) + length(resamp_input[is.na(location_code) == T | location_code == "NA", shapefile])), 3) * 100, "%) referenced shapefiles are missing location_codes"
    ))
  }

  message("Checking referenced shapefiles")
  data <- copy(resamp_input)

  # Confirm necessary columns in data
  nec_cols <- c("shapefile", "location_code")

  if (any(!(nec_cols %in% names(data)))) {
    message("Missing the following columns in your data:")
    print(nec_cols[!nec_cols %in% names(data)])
  }

  # change lat long to latiude longitude if needed
  names(data)[names(data) == "lat"] <- "latitude"
  names(data)[names(data) == "long"] <- "longitude"

  # find records with shapefiles
  shp_idx <- which(!is.na(data$shapefile))
  message(paste(length(shp_idx), "of", nrow(data), "rows of data are polygon assigned
"))

  # identify all shapefiles from the dataset
  all_shapes <- unique(data$shapefile[shp_idx])
  message(paste(length(all_shapes), "unique shapefiles found in data.
"))

  # check they're in the directory
  message(paste0("Checking shapefile directory (", shp_path, ") for matches.. "))
  if (!all(paste0(all_shapes, ".shp") %in% list.files(shp_path))) {
    message("Missing the following shapefiles:")
    print(all_shapes[!(paste0(all_shapes, ".shp") %in% list.files(shp_path))])
    if (!ignore_warnings) message("Function breaking because not all shapefiles are a match. Please check on the above shapefile.")
  } else {
    message("All shapefiles in data match a shapefile by name in the directory.
")
  }

  # report out what countries
  check_counts[, combined_name := paste0(shapefile, "_", location_code)]

  if (nrow(check_counts[output_N == input_N, ]) > 0) {
    message("Summary of polygons resampled to one cluster")
    resamp_input[paste0(shapefile, "_", location_code) %in% check_counts[output_N == input_N, combined_name], .(country, year)] %>%
      table() %>%
      print()
  }
  if (nrow(check_counts[is.na(output_N) == T, ]) > 0) {
    message("Summary of polygons dropped during resampling")
    resamp_input[paste0(shapefile, "_", location_code) %in% check_counts[is.na(output_N) == T, combined_name], .(country, year)] %>%
      table() %>%
      print()
  }


  if (return_list == TRUE) {
    return(list(check_counts[output_N == input_N, ], check_counts[is.na(output_N) == T, ], input_count[is.na(location_code) == T | location_code == "NA", ]))
  }
}
