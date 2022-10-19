#### Multiple imputation for missing covariate values
impute_missing_covariate_values <- function(all_cov_layers = all_cov_layers, simple_raster = simple_raster, year_list = year_list, cores_to_use = 1) {
  ### Resample covariate rasters to match simple_raster
  message("Resampling covariate rasters to match simple_raster...")
  all_cov_layers_resampled <- copy(all_cov_layers)

  for (i in 1:length(all_cov_layers_resampled)) {
    if ((dataType(all_cov_layers_resampled[[i]]) == "LOG1S") | (isTRUE(all.equal(sort(unique(as.data.table(raster::freq(all_cov_layers_resampled[[i]], useNA = "no", merge = TRUE, digits = 1)))$value), c(0, 1))))) { # Binary covariates
      all_cov_layers_resampled[[i]] <- raster::resample(all_cov_layers_resampled[[i]], simple_raster, method = "ngb")
    } else { # Continuous covariates
      all_cov_layers_resampled[[i]] <- raster::resample(all_cov_layers_resampled[[i]], simple_raster)
    }
  }

  ### Impute missing raster cells
  raster_mods <- copy(all_cov_layers_resampled)
  complete_data_raster <- all_cov_layers_resampled[[1]][[1]]

  ### Set any cells that are missing data in at least one covariate to NA
  message("Creating imputation training set...")
  for (i in 1:length(all_cov_layers_resampled)) {
    for (j in 1:nlayers(all_cov_layers_resampled[[i]])) {
      current <- all_cov_layers_resampled[[i]][[j]]
      complete_data_raster <- overlay(current, complete_data_raster, fun = function(x, y) {
        y[is.na(x[])] <- NA
        return(y)
      })
    }
  }

  ### Create imputation training set
  cell_idx <- seegSDM:::notMissingIdx(complete_data_raster)
  imputation_training_IDs <- sample(cell_idx, min(1000, length(cell_idx)), replace = FALSE, prob = NULL)

  if (nlayers(all_cov_layers[[1]]) == 1) {
    index <- 1
  } else {
    index <- which(year_list == sample(year_list, 1))
  }

  imputation_training_data <- data.table(x = all_cov_layers_resampled[[1]][[index]][imputation_training_IDs])
  colnames(imputation_training_data)[1] <- gsub(".1", "", names(all_cov_layers_resampled[[1]]), fixed = TRUE)
  for (i in 2:length(all_cov_layers)) {
    if (nlayers(all_cov_layers[[i]]) == 1) {
      index <- 1
    } else {
      index <- which(year_list == sample(year_list, 1))
    }

    imputation_training_data <- cbind(imputation_training_data, all_cov_layers_resampled[[i]][[index]][imputation_training_IDs])
    colnames(imputation_training_data)[i] <- gsub(".1", "", names(all_cov_layers_resampled[[i]][[1]]), fixed = TRUE)
  }

  imputation_training_data <- cbind(imputation_training_IDs, imputation_training_data)
  colnames(imputation_training_data)[1] <- "cell_idx"

  ### Drop constant covariates from training set
  imputation_training_data <- Filter(var, imputation_training_data)

  ### Create data set for cells with missing data
  missing_data_raster <- overlay(simple_raster, complete_data_raster, fun = function(x, y) {
    x[!is.na(y[])] <- NA
    return(x)
  })

  cell_idx_missing <- seegSDM:::notMissingIdx(missing_data_raster)

  if (length(cell_idx_missing) == 0) {
    message("No cells with missing covariate values... exiting imputation.")
    return(all_cov_layers_resampled)
  }
  
  message(paste0("Imputing missing covariate values for ", length(cell_idx_missing), " cells."))
  
  ### Impute missing values and set raster cells to imputed values
  registerDoParallel(cores = cores_to_use)

  imputed <- foreach(j = 1:length(year_list), .packages = c("data.table", "Hmisc", "raster")) %dopar% {
    message(paste0("Starting imputations for year ", year_list[[j]]))
    all_cov_layers_resampled_imputed <- all_cov_layers_resampled
    missing_data <- data.table(x = all_cov_layers_resampled[[1]][[1]][cell_idx_missing])
    missing_data[] <- NA

    for (i in 1:length(all_cov_layers_resampled_imputed)) {
      if (nlayers(all_cov_layers_resampled_imputed[[i]]) == 1) {
        index <- 1
      } else {
        index <- j
      }
      missing_data <- cbind(missing_data, all_cov_layers_resampled[[i]][[index]][cell_idx_missing])
      colnames(missing_data)[i + 1] <- gsub(".1", "", names(all_cov_layers_resampled_imputed[[i]][[1]]), fixed = TRUE)
    }

    missing_data <- missing_data[, -1]
    missing_data <- cbind(cell_idx_missing, missing_data)
    colnames(missing_data)[1] <- "cell_idx"

    missing_data_subset <- as.data.table(as.data.frame(missing_data)[, colnames(imputation_training_data)]) # Match columns
    missing_data_subset <- missing_data_subset[apply(missing_data_subset, MARGIN = 1, function(x) sum(is.na(x))) < 6]
    missing_data_subset <- missing_data_subset[apply(missing_data_subset, MARGIN = 1, function(x) sum(is.na(x))) > 0]
    covars_to_impute <- data.table(t(apply(missing_data_subset, MARGIN = 2, function(x) sum(is.na(x)))))
    col_names <- colnames(covars_to_impute)[which(covars_to_impute > 0)]

    for (i in 1:nrow(missing_data_subset)) {
      message(paste0("Imputing row ", i, " of ", nrow(missing_data_subset)))
      current_training <- rbind(imputation_training_data, missing_data_subset[i])
      current_training <- current_training[, lapply(.SD, function(v) if (uniqueN(v, na.rm = FALSE) > 1) v)]
      impute_arg <- eval(parse(text = paste("aregImpute( ~", paste(colnames(current_training)[2:ncol(current_training)], collapse = "+"), ", nk=0, data=current_training, n.impute=5)")))
      for (k in 1:length(impute_arg$imputed)) {
        if (!is.null(impute_arg$imputed[k][[1]])) {
          all_cov_layers_resampled_imputed[[names(impute_arg$imputed[k])]][missing_data_subset[i, cell_idx]][j] <- median(impute_arg$imputed[k][[1]])
        }
      }
    }

    return(list("j" = j, "imputed" = all_cov_layers_resampled_imputed))
  }

  stopImplicitCluster()

  ### Reorganize imputed rasters to match all_cov_layers
  imputed_final <- all_cov_layers_resampled

  message("Wrapping up imputations...")
  for (i in 1:length(imputed)) { # Cycle through years
    for (j in 1:length(imputed[[1]]$imputed)) { # Cycle through covariates
      if (nlayers(imputed_final[[j]]) == 1) { # If the variable is synoptic...
        imputed_final[[j]] <- imputed[[i]]$imputed[[j]]
      } else { # ...the variable is not synoptic.
        imputed_final[[j]][[i]] <- imputed[[i]]$imputed[[j]][[i]]
      }
    }
  }

  return(imputed_final)
}
