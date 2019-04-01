# HEADER ------------------------------------------------------------------
# 05_generate_vaccine_csvs.R
# Purpose: Generate final csv files in format ready for MBG models / 
#          ordinal regression
#**************************************************************************

# Note: run from 00_master.R, which sets up all directories & files
# Outputs: input data csvs (in "out_dir") ready for MBG models

for (vax in vaccines) {

  # Define objects 
  vax_prefix <- vaccine_list[[vax]][["prefix"]]
  vax_title  <- vaccine_list[[vax]][["title"]]
  vax_doses  <- vaccine_list[[vax]][["all_doses"]]
  cond_vax   <- vaccine_list[[vax]][["cond_vaccines"]]
  all_vax    <- vaccine_list[[vax]][["all_vaccines"]]

  # Load data
  df_vax <- readRDS(paste0(vaccine_resampled_file_prefix, vax_prefix, ".rds"))

  # Add a rowID common between all derivative data sets
  # (for use in creation of holdouts)
  df_vax$row_id <- seq.int(nrow(df_vax))

  # Define list of vaccines to save csvs for
  # Need _cov for last dose, _cond for intermediate doses, and _cov for
  # any doses that are going to be the ultimate models (for validation)

  last_vax <- all_vax[length(all_vax)]
  last_dose <- max(vax_doses)

  csv_vaccines <- unique(c(cond_vax, 
                           paste0(vax_prefix, last_dose, "_cov"), 
                           final_model_vaccines))

  message(paste0("\nSaving .csv files for ", vax_title, "..."))

  for (ind in csv_vaccines) {

    # Create csv for the last dose -----------------------------------------

    df_temp <- copy(df_vax)

    # Extract dose number
    dose <- str_match(ind, "^[a-z]*([0-9]?)_.*")[2] %>% as.numeric 

    if (grepl("_cov", ind)) {

      # Add up all of the dose_X columns to get P(doses >= dose)

      df_temp[, outcome := rowSums(.SD),
                .SDcols = paste0(vax, "_dose_", dose:last_dose)]

      # Drop extra columns
      drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
      df_temp <- subset(df_temp, select = !names(df_temp) %in% drop_vars)

    } else if (grepl("_cond", ind)) {

      setnames(df_temp, paste0(vax_prefix, "_dose_", 0:dose), paste0("d", 0:dose))

      # Drop other dose columns (ignoring: conditional)
      drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
    
      # Subset to just the columns of interest
      df_temp <- subset(df_temp, select = !(names(df_temp) %in% drop_vars))

      # Overwrite N with sum of remaining variables
      df_temp[, N := rowSums(.SD), .SDcols = paste0("d", 0:dose)]

      # Drop extra columns
      df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:(dose-1))))

      # Rename to "outcome" for standardization
      setnames(df_temp, paste0("d", dose), "outcome")

      # Drop if N = 0 (none in that row have [dose] or fewer doses)
      df_temp <- subset(df_temp, N > 0)

    } else if (grepl("1_3_abs_dropout", ind)) {

      setnames(df_temp, paste0(vax_prefix, "_dose_", 0:last_dose), paste0("d", 0:last_dose))
      df_temp[, outcome := d1 + d2]  

      # Drop extra columns
      df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:last_dose)))

    } else if (grepl("1_3_rel_dropout", ind)) {

      # Get the denominator (P doses > 1)
      df_temp[, denom := rowSums(.SD),
                .SDcols = paste0(vax, "_dose_", 1:last_dose)]

      setnames(df_temp, paste0(vax_prefix, "_dose_", 0:last_dose), paste0("d", 0:last_dose))

      # Calculate relative dropout
      df_temp[, outcome := (d1 + d2) / denom]  

      # Drop extra columns
      df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:last_dose)))
      df_temp[, denom := NULL]

    }

    # Rename outcome to indicator
    setnames(df_temp, "outcome", ind)

    # Write
    message(paste0("\nWriting csv for ", ind, ". File can be found at: ",
                   "\n   ", out_dir, ind, ".csv"))
    write.csv(df_temp, paste0(out_dir, "/", ind, ".csv"))

  }

}
