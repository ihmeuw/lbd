# HEADER ------------------------------------------------------------------
# 04_pointpoly.R
# Purpose: Run polygon resampling algorithm to create candidate points
#          from polygon data           
#**************************************************************************

# Note: run from 00_master.R, which sets up all directories & files
# Outputs: pointpoly-positioned rds files

missing_shapefiles <- list()

for (vax in vaccines) {

  # Define objects 
  vax_prefix <- vaccine_list[[vax]][["prefix"]]
  vax_title  <- vaccine_list[[vax]][["title"]]
  vax_doses  <- vaccine_list[[vax]][["all_doses"]]
  age_range  <- vaccine_list[[vax]][["age_range"]]
  age_min <- as.numeric(min(age_range))
  age_max <- as.numeric(max(age_range))

  # Load data
  df_input <- readRDS(paste0(vaccine_cleaned_file_prefix, vax_prefix, ".rds"))

  ## Crop input data to Africa -----------------------------------------
  gaul_list <- get_gaul_codes('africa') 
  gaul_list <- gaul_list[gaul_list != 269] # Remove Yemen
  gaul_list <- c(gaul_list, 40764) # Add SDN

  # create a matching variable that excludes everything after "_" (subnationals)
  df_input[, match_country := country]
  df_input[grep("_", df_input$country), match_country :=  str_match(country, "([^_]+)_.*")[,2]]

  t_lookup <- unique(df_input$match_country) %>% as.character %>% as.data.table
  names(t_lookup) <- "match_country" 
  t_lookup[, gaul := gaul_convert(match_country)]

  # Merge lookup table
  df_input <- merge(df_input, t_lookup, all.x = T, by = "match_country")
  df_input[, match_country := NULL]

  # Drop if no gaul code found
  if (length(unique(df_input[is.na(gaul)]$country)) > 0) {
    message(paste0("No gaul codes found for ", 
                   paste(unique(df_input[is.na(gaul)]$country), collapse = ", "),
                   " - dropping..."))
  }
  df_input <- df_input[!is.na(gaul)]

  # Drop if not in gaul list & remove gaul variable
  keep_countries <- unique(df_input[gaul %in% gaul_list]$country)
  drop_countries <- unique(df_input[!(gaul %in% gaul_list)]$country)
  message(paste0("Dropping the following countries (not in gaul_list): ", paste(drop_countries, collapse = ", ")))
  df_input <- df_input[gaul %in% gaul_list]
  df_input[, gaul := NULL]

  message("\n###################################################")
  message(paste0("Working on polygon->point for ", vax))
  message("###################################################\n")

  if (nrow(df_input[point == 0]) > 0) {

    # Run the time-intensive point/poly solution (now better in parallel!)
    df_pointpoly <- resample_polygons_fast(data = df_input, 
                                           cores = cores, 
                                           indic = vax, 
                                           density = 0.001)
  
  } else {
    # If no polys, don't resample!
    df_pointpoly <- df_input
  }

  ### Write to output file & list------------------------------------------
  
  message(paste0("Saving resampled data for ", vax_title, " here: "))
  message(paste0("  ", vaccine_resampled_file_prefix, vax_prefix, ".rds"))
  saveRDS(df_pointpoly, paste0(vaccine_resampled_file_prefix, vax_prefix, ".rds"))

}
