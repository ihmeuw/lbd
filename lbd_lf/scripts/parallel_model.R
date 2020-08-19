#####################################################################
## Parallel script for running MBG models                          ##
#####################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(<<<< FILEPATH REDACTED >>>>)
load_from_parallelize()
message(paste0("Using ", core_repo))

user <- Sys.info()[["user"]]

## Set repo locations
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

path <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- paste0(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(matrixStats)
library(sf)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator,
                                config_file = paste0(<<<< FILEPATH REDACTED >>>>), run_tests = FALSE)

age <- 0
holdout <- as.integer(holdout)

outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)
dir.create(outputdir, showWarnings = FALSE)

## make a pathaddin that get used widely
pathaddin <- paste0('_bin', 0, '_', region, '_', holdout)

# print run options
message("options for this run:\n")
for (arg in c(
  "region", "age", "run_date", "test", "holdout",
  "indicator", "indicator_group", "pathaddin", "outputdir"
)) {
  message(paste0(arg, ":    ", get(arg), "     // type: ", class(get(arg))))
}

# We need to be in the singularity image, and specifically the LBD one if using TMB
if(!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if(as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Make sure this inla patch is implemented if running on geos
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

## cores to use
cores_to_use <- Sys.getenv("SGE_HGR_fthread")
message(paste("Model set to use", cores_to_use, "cores"))

# print out session info so we have it on record
sessionInfo()

source(paste0(<<<< FILEPATH REDACTED >>>>))

message(paste0("perform_initial_INLA_approx: ", perform_initial_INLA_approx))
message(paste0("initial_INLA_approx_levels: ", initial_INLA_approx_levels))

reg <- region

#reload spatial templates
load(paste0(<<<< FILEPATH REDACTED >>>>))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (class(year_list) == "character") year_list <- eval(parse(text = year_list))
if (class(predict_years) == "character") predict_years <- eval(parse(text = predict_years))

# determine indices based on year_list for pdf plotting later on.
y_start <- seq(1, length(year_list), 16)
if (length(year_list) > 16) {
  y_end <- c(seq(16, length(year_list), 16), length(year_list))
} else {
  y_end <- 16
}

y_start_pred <- seq(1, length(unlist(predict_years)), 16)
if (length(unlist(predict_years)) > 16) {
  y_end_pred <- c(seq(16, length(unlist(predict_years)), 16), length(unlist(predict_years)))
} else {
  y_end_pred <- 16
}

##############################################################
# skip a large chunk if requested in config
if (!as.logical(skiptoinla)) {
  message("You have chosen to not skip directly to inla.")
  
  ## some set up
  if (class(z_list) == "character") z_list <- eval(parse(text = z_list))
  
  ## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
  # For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
  if (holdout != 0) {
    message(paste0("Holdout != 0 so loading holdout data only from holdout ", holdout))
    message("Please be sure you have a list object called stratum_ho in your environment.")
    # if stratifies by age then make sure loads correctly
    if (age != 0) df <- as.data.table(stratum_ho[[paste("region", region, "_age", age, sep = "__")]])
    if (age == 0) df <- as.data.table(stratum_ho[[paste("region", region, sep = "__")]])
    
    df[, "indicator_original" := get(indicator)]
    
    if (oos_preds_inla) {
      df_holdout <- df[fold == holdout]
      df_holdout[, (indicator) := NA]
      
      df <- rbind(df[fold != holdout, ], df_holdout)
      
    } else {
      df <- df[fold != holdout, ]
    }
    df$first_entry <- 1
    df <- df[region == region,]
    run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
    write.csv(df, file = paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
  } else {
    run_date_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
    dir.create(run_date_dir, showWarnings = FALSE)
    df <- df[region == region,]
    write.csv(df, file = paste0(<<<< FILEPATH REDACTED >>>>), row.names = F)
  }
  
  
  if(as.logical(use_subnat_res)) {
    
    #get adm0 codes for countries to get subnat REs
    if("all" %in% subnat_country_to_get){
      countries_to_get_subnat_res <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
    } else {
      countries_to_get_subnat_res <- get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version)[get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version) %in% get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)]
    }
    
    #load and subset standard admin1 shape to countries with subnational REs 
    subnat_full_shp    <- readRDS(get_admin_shapefile( admin_level = 1, raking = F, suffix = '.rds', version = modeling_shapefile_version ))
    subnat_shapefile <- raster::subset(subnat_full_shp, 
                                       ADM0_CODE %in% countries_to_get_subnat_res)
    
    simple_polygon_list2 <- load_simple_polygon(gaul_list = NULL, 
                                                buffer = 1, tolerance = 0.4, 
                                                custom_shapefile = subnat_shapefile)
    subset_shape2        <- simple_polygon_list2[[1]]
    
    ## Load list of raster inputs (pop and simple)
    raster_list2        <- build_simple_raster_pop(subset_shape2, field = "ADM1_CODE")
    simple_raster2      <- raster_list2[['simple_raster']]
    
    #simple_raster2 has to be the same size as the simple raster for predict to work correctly
    simple_raster2 <- raster::extend(simple_raster2, extent(simple_raster))
    simple_raster2 <- raster::crop(simple_raster2, extent(simple_raster))
    
    ## Merge ADM0/1 codes to df
    adm1_subset_lox <- over(SpatialPoints(df[,.(long = longitude, lat = latitude)], 
                                          CRS(proj4string(subnat_shapefile))), subnat_shapefile)
    df[, subnat_re_ADM1_CODE := as.numeric(as.character(adm1_subset_lox$ADM1_CODE))]
    df[, subnat_re_ADM0_CODE := as.numeric(as.character(adm1_subset_lox$ADM0_CODE))]
    
    #create new ADM1 columns for each country in subnat_country_to_get so data can be fit separately
    for(i in 1:length(unique(na.omit(df$subnat_re_ADM0_CODE)))) {
      df[subnat_re_ADM0_CODE == unique(na.omit(df$subnat_re_ADM0_CODE))[i], (paste0("SUBNAT", i)) := subnat_re_ADM1_CODE]
    }
  } else {
    simple_raster2 <- NULL
  }
  
  # if there is another weight column, multiply it with weight now
  if (exists("other_weight")) {
    if (other_weight != "") {
      message(paste0("Multiplying weight and ", other_weight))
      df[["weight"]] <- df[["weight"]] * df[[other_weight]]
    }
  }
  
  ## Some built in data checks that cause known problems later on
  if (indicator_family %in% c("binomial") & any(df[!is.na(get(indicator)), get(indicator) / N] > 1)) {
    stop("You have binomial data where k > N. Check your data before proceeding")
  }
  
  if (any(df[["weight"]] %in% c(Inf, -Inf) | any(is.na(df[["weight"]])))) {
    stop("You have illegal weights (NA, Inf, -Inf). Check your data before proceeding")
  }
  
  ## Save distribution of data for this region
  png(paste0(<<<< FILEPATH REDACTED >>>>))
  if(indicator_family=='binomial') hist(df[df$first_entry==1, get(indicator)]/df$N[df$first_entry==1]) else hist(df[df$first_entry==1, get(indicator)])
  dev.off()
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Define modeling space. In years only for now.
  if (yearload=='annual') period_map <- make_period_map(modeling_periods = c(min(year_list):max(year_list)))
  if (yearload=='five-year') period_map <- make_period_map(modeling_periods = seq(min(year_list), max(year_list), by = 5))
  
  # make covariates conditional
  cov_layers <- gbd_cov_layers <- NULL
  
  ## Pull all covariate bricks/layers
  if (nchar(fixed_effects) > 0) {
    message("Grabbing raster covariate layers")
    
    cov_layers <- load_lagged_covariates(covariate_config = fixed_effects_config,
                                         template = simple_raster_raf1,
                                         start_year = min(year_list),
                                         end_year = max(year_list),
                                         raster_agg = as.integer(raster_agg_factor))
  }
  
  
  ## Combine all covariates
  all_cov_layers <- c(cov_layers, gbd_cov_layers)
  
  # regenerate all fixed effects equation from the cov layers
  all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")
  
  ## Make stacker-specific formulas where applicable
  all_fixed_effects_brt <- all_fixed_effects
  
  # Set Up Country Fixed Effects
  if (as.logical(use_child_country_fes) | as.logical(use_inla_country_fes) | as.logical(country_re_pp)) {
    message("Setting up country fixed effects")
    fe_gaul_list <- sort(unique(c(gaul_convert(unique(df[, country]), shapefile_version = modeling_shapefile_version), gaul_list)))
    fe_template <- cov_layers[[1]][[1]]
    
    fe_subset_shape <- subset_shape[subset_shape@data$ADM0_CODE %in% fe_gaul_list,]
    
    gaul_code <- build_simple_raster_pop(fe_subset_shape, field = 'ADM0_CODE')$simple_raster
    gaul_code <- raster::resample(gaul_code, cov_layers[[1]][[1]], method="ngb")
    gaul_code <- setNames(gaul_code,'gaul_code')
    gaul_code <- create_categorical_raster(gaul_code)
    
    ## update covlayers and add country fixed effects
    all_cov_layers <- update_cov_layers(all_cov_layers, gaul_code)
    all_fixed_effects_cfes <- paste(all_fixed_effects,paste(names(gaul_code)[1:length(names(gaul_code))], collapse = " + "), sep = " + ")
    
    ## update specific stacker formulas (for now we just want country effects in BRT)
    all_fixed_effects_brt <- all_fixed_effects_cfes
  }
  
  # Add these to the fixed effects if we want them in stacking
  if (as.logical(use_child_country_fes) | as.logical(country_re_pp)) {
    gaul_fes <- paste(names(gaul_code)[1:length(names(gaul_code))], collapse = " + ")
    all_fixed_effects <- paste(all_fixed_effects, gaul_fes, sep = " + ")
  }
  
  if (use_age_diagnostics_in_stackers & as.logical(use_stacking_covs)) {
    age_midpoint_raster <- create_age_midpoint_rasters(template = simple_raster, start_year = min(year_list), end_year = max(year_list), reg = region, shapefile_version = modeling_shapefile_version)
    all_cov_layers[["age_midpoint"]] <- age_midpoint_raster
    
    df$diagnostic <- as.factor(tolower(df$diagnostic))
    predict_diagnostic <- tolower(predict_diagnostic)
    diagnostic_raster <- copy(simple_raster)
    values(diagnostic_raster) <- which(levels(df$diagnostic) == predict_diagnostic)
    diagnostic_raster <- raster::mask(diagnostic_raster, simple_raster)
    all_cov_layers[["diagnostic"]] <- diagnostic_raster
    
    all_fixed_effects <- paste0(all_fixed_effects, " + age_midpoint", " + diagnostic")
    all_fixed_effects_brt <- paste0(all_fixed_effects_brt, " + age_midpoint", " + diagnostic")
    
    df$diagnostic <- as.integer(df$diagnostic)
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Stacking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  # Figure out which models we're going to use
  child_model_names <- stacked_fixed_effects %>%
    gsub(" ", "", .) %>%
    strsplit(., "+", fixed = T) %>%
    unlist()
  message(paste0("Child stackers included are: ", paste(child_model_names, collapse = " // ")))
  
  the_covs <- format_covariates(all_fixed_effects)
  
  # copy the dataset to avoid unintended namespace conflicts
  the_data <- copy(df) %>% as.data.table()
  
  # shuffle the data into folds for stacking
  the_data <- the_data[sample(nrow(the_data)), ]
  the_data[, fold_id := cut(seq(1, nrow(the_data)), breaks = as.numeric(n_stack_folds), labels = FALSE)]
  
  # add a row id column
  the_data[, a_rowid := seq(1:nrow(the_data))]
  
  # plot covariates as a simple diagnostic here
  pdf(sprintf(<<<< FILEPATH REDACTED >>>>), height = 12, width = 12)
  for (covname in names(all_cov_layers)) {
    c <- copy(all_cov_layers[[covname]])
    if (nlayers(c) > 1) {
      names(c) <- paste0(covname, "_", year_list)
      for (i in 1:length(y_start)) {
        plot(c[[y_start[i]:y_end[i]]], maxpixel = 1e6)
      }
    } else {
      names(c) <- paste0(covname)
      plot(c, main = covname, maxpixel = 1e6)
    }
  }
  dev.off()
  rm(c)
  
  if (as.logical(impute_missing_covariates)) {
    if (!as.logical(use_saved_imputations)) {
      #### Impute raster cells with missing covariate values
      all_cov_layers_imputed <- impute_missing_covariate_values(all_cov_layers = all_cov_layers, simple_raster = simple_raster, year_list = year_list, cores_to_use = min(7, as.integer(fthread_parallel)))
      saveRDS(all_cov_layers_imputed, file = paste0(<<<< FILEPATH REDACTED >>>>))
    } else {
      all_cov_layers_imputed <- readRDS(paste0(<<<< FILEPATH REDACTED >>>>))
    }
    
    # plot resampled and imputed covariates as a simple diagnostic here
    pdf(sprintf(<<<< FILEPATH REDACTED >>>>), height = 12, width = 12)
    for (covname in names(all_cov_layers_imputed)) {
      c <- copy(all_cov_layers_imputed[[covname]])
      if (nlayers(c) > 1) {
        names(c) <- paste0(covname, "_", year_list)
        for (i in 1:length(y_start)) {
          plot(c[[y_start[i]:y_end[i]]], maxpixel = 1e6)
        }
      } else {
        names(c) <- paste0(covname)
        plot(c, main = covname, maxpixel = 1e6)
      }
    }
    dev.off()
    rm(c)
    
    all_cov_layers <- all_cov_layers_imputed
  }
  
  ### Clamp covariates if desired
  if (clamp_covariates) {
    if (!is.null(all_cov_layers)) {
      cs_covs_no_cs <- extract_covariates(the_data,
                                          all_cov_layers,
                                          id_col = "a_rowid",
                                          return_only_results = TRUE,
                                          centre_scale = FALSE,
                                          period_var = "year",
                                          period_map = period_map
      )
      
      for (i in 1:length(all_cov_layers)) {
        print(i)
        current_cov <- attr(all_cov_layers[i], "names")
        range_min <- min(cs_covs_no_cs[, get(current_cov)], na.rm = TRUE)
        range_max <- max(cs_covs_no_cs[, get(current_cov)], na.rm = TRUE)
        for (j in 1:nlayers(all_cov_layers[[i]])) {
          print(j)
          vals <- values(all_cov_layers[i][[1]][[j]])
          if (sum(vals < range_min, na.rm=TRUE) > 0) {
            vals[vals < range_min] <- range_min
            all_cov_layers[i][[1]] <- setValues(all_cov_layers[i][[1]], values=vals, layer=j)
          }
          
          if (sum(vals > range_max, na.rm=TRUE) > 0) {
            vals[vals > range_max] <- range_max
            all_cov_layers[i][[1]] <- setValues(all_cov_layers[i][[1]], values=vals, layer=j)
          }
        }
      }
    }
    
    # plot resampled and imputed covariates as a simple diagnostic here
    pdf(sprintf(<<<< FILEPATH REDACTED >>>>), height = 12, width = 12)
    for (covname in names(all_cov_layers)) {
      c <- copy(all_cov_layers[[covname]])
      if (nlayers(c) > 1) {
        names(c) <- paste0(covname, "_", year_list)
        for (i in 1:length(y_start)) {
          plot(c[[y_start[i]:y_end[i]]], maxpixel = 1e6)
        }
      } else {
        names(c) <- paste0(covname)
        plot(c, main = covname, maxpixel = 1e6)
      }
    }
    dev.off()
    rm(c)
  }
  
  # extract covariates to the points and subset data where it's missing covariate values
  if (!is.null(all_cov_layers)) {
    cs_covs <- extract_covariates(the_data,
                                  all_cov_layers,
                                  id_col = "a_rowid",
                                  return_only_results = TRUE,
                                  centre_scale = TRUE,
                                  period_var = "year",
                                  period_map = period_map
    )
  } else {
    cs_covs <- NULL
  }
  
  # A check to see if any of the variables do not vary across the data. This could break model later so we check and update some objects
  covchecklist <- check_for_cov_issues(check_pixelcount = check_cov_pixelcount, check_pixelcount_thresh = ifelse(exists("pixelcount_thresh"), as.numeric(pixelcount_thresh), 0.95), drop_nonvarying_covs = FALSE)
  for (n in names(covchecklist)) {
    assign(n, covchecklist[[n]])
  }
  
  ## Check for data where covariate extraction failed
  if (is.null(cs_covs)) {
    rows_missing_covs <- 0
  } else {
    rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
  }
  
  if (rows_missing_covs > 0) {
    pct_missing_covs <- round((rows_missing_covs / nrow(the_data)) * 100, 2)
    warning(paste0(
      rows_missing_covs, " out of ", nrow(the_data), " rows of data ",
      "(", pct_missing_covs, "%) do not have corresponding ",
      "covariate values and will be dropped from child models..."
    ))
    if (rows_missing_covs / nrow(the_data) > 0.1) {
      stop(paste0(
        "Something has gone quite wrong: more than 10% of your data does not have ",
        "corresponding covariates.  You should investigate this before proceeding."
      ))
    }
  }
  
  if (!is.null(cs_covs[[1]])) {
    the_data <- merge(the_data, cs_covs[[1]][, c(the_covs[!(the_covs %in% c("age_midpoint", "diagnostic"))], "a_rowid"), with = F], by = "a_rowid", all.x = F, all.y = F) # added subset to make sure no weird columns are added
  }
  
  # store the centre scaling mapping
  covs_cs_df <- cs_covs[[2]]
  
  ### Nudge points with NA covariate extractions to nearest grid cell with covariate data, if within a threshold distance.
  if (!is.null(cs_covs[[1]])) {
    nudged_results <- nudge_coordinates(temp_data = the_data, the_covs = the_covs, covs_cs_df = covs_cs_df, threshold = 10000)
    the_data <- nudged_results[[1]]
    covs_cs_df <- nudged_results[[2]]
  }
  
  # this will drop rows with NA covariate values
  if (oos_preds_inla) { # Retain rows for holdout fold (i.e., rows with NA indicator value)
    the_data_na <- na.omit(the_data, c("N", the_covs), invert = TRUE)
    the_data <- na.omit(the_data, c("N", the_covs))
  } else {
    the_data_na <- na.omit(the_data, c(indicator, "N", the_covs), invert = TRUE)
    the_data <- na.omit(the_data, c(indicator, "N", the_covs))
  }
  
  message(paste0(nrow(the_data_na), " rows being dropped due to NA covariate values. These have been saved as dropped_data_NA_covariates.csv."))
  output_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  write.csv(the_data_na, paste0(<<<< FILEPATH REDACTED >>>>), row.names = FALSE)
  
  if (as.logical(use_stacking_covs)) {
    message('Fitting Stackers')
    # Run the child stacker models
    child_model_run <- run_child_stackers(models = child_model_names, input_data = the_data)
    
    # Bind the list of predictions into a data frame
    child_mods_df <- do.call(cbind, lapply(child_model_run, function(x) x[[1]]))
    
    ## Combine the children models with the_data
    the_data  <- cbind(the_data, child_mods_df)
    
    ## Rename the child model objects into a named list
    child_model_objs <- setNames(lapply(child_model_run, function(x) x[[2]]), child_model_names)
    
    ## return the stacked rasters
    stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                          period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                          child_models     = child_model_objs,
                                          indicator_family = indicator_family,
                                          centre_scale_df  = covs_cs_df)
    
    ## Plot stackers
    pdf(paste0(<<<< FILEPATH REDACTED >>>>))
    for (i in 1:length(stacked_rasters)) {
      s_r <- copy(stacked_rasters[[i]])
      for (i in 1:length(y_start_pred)) {
        plot(s_r[[y_start_pred[i]:y_end_pred[i]]])
      }
    }
    dev.off()
    rm(s_r)
    
    message('Stacking is complete')
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Set the fixed effects to use in INLA based on config args
  if (!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)) {
    all_fixed_effects <- ""
  } else if (as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)) {
    all_fixed_effects <- stacked_fixed_effects
  } else if (!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)) {
    all_fixed_effects <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
  } else if (as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)) {
    all_fixed_effects <- stacked_fixed_effects
  } else if (as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)) {
    all_fixed_effects <- paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
  } else if (!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)) {
    all_fixed_effects <- paste(all_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ") # Fixed effects names when using interaction (currently hard coded into this script)
  } else if (as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)) {
    all_fixed_effects <- paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ") # Fixed effects names when using interaction (currently hard coded into this script)
  }
  
  ## Copy things back over to df
  df <- copy(the_data)
  
  ## Remove the covariate columns so that there are no name conflicts when they get added back in
  df <- df[, paste0(the_covs) := rep(NULL, length(the_covs))]
  
  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
  df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]
  
  ## Create a full raster list to carry though to the shiny/next steps
  if (as.logical(use_stacking_covs)) {
    cov_list      <- c(unlist(stacked_rasters),unlist(all_cov_layers))
    child_mod_ras <- cov_list[child_model_names]
  } else {
    cov_list <- unlist(all_cov_layers)
    child_model_names <- ''
  }
  
  ## make sure this inla patch is implemented if running on geos
  if (grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
  
  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d = df,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset,
                             s2mesh = as.logical(use_s2_mesh),
                             s2params = s2_mesh_params,
                             cutoff = mesh_cutoff)
  
  ## Build temporal mesh (standard for now)
  if (length(unique(year_list)) == 1) {
    mesh_t <- NULL
  } else {
    mesh_t <- build_time_mesh(periods = eval(parse(text = mesh_t_knots)))
  }
  
  output_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  save(mesh_s, mesh_s_ppp, mesh_t, file = (paste0(<<<< FILEPATH REDACTED >>>>))) # Save mesh_s and mesh_t to file for later use
  
  ## ## For raw covs, don't want to center-scale (as that will happen in `build_mbg_data_stack()`)
  ##
  ## ## This is a bit weird, but when stacking covs are used the oos-stackers (used in `fit_mbg()`)
  ## ## do not get center scaled in save_mbg_input() - this just harmonizes the measures.  If this
  ## ## step isn't done, then the covs get double-center-scaled with odd results.
  ##
  ## ## For predict_mbg, non-center-scaled covs are pulled from cov_list (either stackes or raw) and
  ## ## center-scaled within the function.  So both fit and predict take place on center-scaled covs
  
  if (as.logical(use_raw_covs) == TRUE) {
    centre_scale_covs <- FALSE
  } else {
    centre_scale_covs <- TRUE
  }
  
  ## Save the full environment image
  save.image(paste0(<<<< FILEPATH REDACTED >>>>))
  
  ## Save all inputs for MBG model into correct location on /share
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 simple_raster_raf1= simple_raster_raf1,
                 mesh_s            = mesh_s,
                 mesh_s_ppp        = mesh_s_ppp,
                 mesh_t            = mesh_t,
                 cov_list          = cov_list,
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = centre_scale_covs)
  
} else { ## BEGIN SKIPTOINLA
  message("Now copying saved MBG inputs from that chosen run_date.")
  
  file.copy(
    from = paste0(<<<< FILEPATH REDACTED >>>>),
    to = paste0(<<<< FILEPATH REDACTED >>>>)
  )
}

## Reload data an prepare for MBG
load(paste0(<<<< FILEPATH REDACTED >>>>))

if (!exists("simple_raster2")) {
  simple_raster2 <- NULL
}

if (as.logical(skiptoinla) & !as.logical(skipinla)) { # rebuild space and time meshes if skiptoinla but not skipinla
  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d = df,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset,
                             s2mesh = as.logical(use_s2_mesh),
                             s2params = s2_mesh_params,
                             cutoff = mesh_cutoff)
  
  ## Build temporal mesh (standard for now)
  if (length(unique(year_list)) == 1) {
    mesh_t <- NULL
  } else {
    mesh_t <- build_time_mesh(periods = eval(parse(text = mesh_t_knots)))
  }
  
  output_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  save(mesh_s, mesh_s_ppp, mesh_t, file = (paste0(<<<< FILEPATH REDACTED >>>>))) # Save mesh_s and mesh_t to file for later use
}

load((paste0(<<<< FILEPATH REDACTED >>>>))) # Load mesh_s and mesh_t from file

# Bound GBM to 0-1 if desired
if (exists("gbm_bounded_0_1")) {
  if (as.logical(gbm_bounded_0_1) == T & !is.null(cov_list[["gbm"]])) {
    message("Truncating GBM values > 1 to 0.999")
    values(cov_list[["gbm"]])[values(cov_list[["gbm"]]) >= 1 & !is.na(values(cov_list[["gbm"]]))] <- 0.999
    gbm_cols <- grep(paste0("(gbm)(.*_pred)"), names(df), value=T)
    replace_one <- function(x) {
      x[x>=1 & !is.na(x)] <- 0.999
      return(x)
    }
    df[, (gbm_cols) := lapply(.SD, replace_one), .SDcols = gbm_cols]
  }
}

## Convert stackers to transform space, if desired
## NOTE: we do this here to ensure that the stacker rasters are saved in prevalence/untransformed space
## this is useful for diagnostics and other code that was built expecting the untransformed rasters
if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)) {
  message('Converting stackers to logit space')
  
  ## Transform the rasters
  for (ii in child_model_names) {
    
    ## Preserve variable names in the raster first
    tmp_rastvar <- names(cov_list[[ii]])
    
    ## Logit
    cov_list[[ii]] <- logit(cov_list[[ii]])
    
    ## Reassign names
    names(cov_list[[ii]]) <- tmp_rastvar
    rm(tmp_rastvar)
  }
  
  ## Transform the stacker values that are in df
  stacker_cols <- grep(paste0("(", paste(child_model_names, collapse="|"), ")(.*_pred)"), names(df), value=T)
  df[, (stacker_cols) := lapply(.SD, logit), .SDcols = stacker_cols]
  
}

if (yearload == "annual") period_map_predict <- make_period_map(modeling_periods = c(min(unlist(predict_years)):max(unlist(predict_years))))
if (yearload == "five-year") period_map_predict <- make_period_map(modeling_periods = seq(min(unlist(predict_years)), max(unlist(predict_years)), by = 5))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## for stacking, overwrite the columns matching the model_names so that we can trick INLA into being our stacker
if (as.logical(use_stacking_covs)) {
  df[, paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x, "_cv_pred")))]
}

## Generate MBG formula for INLA call (will run but not used by TMB)
t_mod_type <- paste0("'", t_mod_type, "'")

if (!as.logical(fit_with_tmb)) {
  mbg_formula <- build_mbg_formula_with_priors(fixed_effects               = all_fixed_effects,
                                               add_nugget                  = use_inla_nugget,
                                               nugget_prior                = nugget_prior,
                                               add_ctry_res                = use_inla_country_res,
                                               ctry_re_prior               = ctry_re_prior,
                                               temporal_model_type         = t_mod_type,
                                               temporal_model_theta_prior  = temporal_model_theta_prior,
                                               temporal_model_theta1_prior = rho_prior,
                                               temporal_model_pc_prec_prior = temporal_model_pc_prec_prior,
                                               temporal_model_pc_pacf1_prior = temporal_model_pc_pacf1_prior,
                                               no_gp                       = !as.logical(use_gp),
                                               coefs.sum1                  = coefs_sum1,
                                               use_space_only_gp           = use_space_only_gp,
                                               use_time_only_gmrf          = use_time_only_gmrf,
                                               use_pref_samp_pp            = as.logical(use_pref_samp_pp),
                                               pref_samp_pp_stages         = pref_samp_pp_stages,
                                               covars_in_pp                = covars_in_pp,
                                               country_re_pp               = country_re_pp,
                                               fit_crosswalk_inla          = fit_crosswalk_inla,
                                               timebycountry_RE            = use_timebyctry_res,
                                               adm0_list                   = gaul_list,
                                               crosswalk_training_data_set = crosswalk_training_data_set,
                                               use_rw1_crosswalk_covariates = use_rw1_crosswalk_covariates,
                                               pc_priors_ar_cor1            = pc_priors_ar_cor1,
                                               use_age_diagnostics_in_stackers = use_age_diagnostics_in_stackers,
                                               subnat_RE = use_subnat_res,
                                               subnat_country_to_get = subnat_country_to_get,
                                               subnat_re_prior = subnat_re_prior)
}

## For INLA we need to add data for missing time points to ensure we get predictions
##  for all relevant time points. The 0 observations do not contribute to the 
##  model fitting but they prevent INLA from auto-removing 
##  random effects that (conditionally) have no data impacting their fit
if (!as.logical(fit_with_tmb)) {
  if(use_timebyctry_res) {
    ## If we are using a time only effect by country then we need to make sure 
    ##  all year effects are estimated for each country.
    df$adm0code <- gaul_convert(df$country)
    for(adm0_code in gaul_list) {
      dfsub <- df[df$adm0code == adm0_code, ]
      missing_years <- setdiff(year_list, dfsub$year)
      if (length(missing_years) > 0) {
        fake_data <- dfsub[1:length(missing_years), ]
        fake_data[, year := missing_years]
        fake_data[, c(indicator, 'N', 'weight') := 0]
        fake_data[, period := NULL]
        fake_data <- merge(fake_data, period_map)
        df <- rbind(df, fake_data)
      }
    }
  } else {
    ## If not, we only need to make sure we have an observation for each missing
    ##  year (country irrelevant)
    missing_years <- setdiff(year_list, df$year)
    if (length(missing_years) > 0) {
      fake_data <- df[1:length(missing_years), ]
      fake_data[, year := missing_years]
      fake_data[, c(indicator, 'N', 'weight') := 0]
      fake_data[, period := NULL]
      fake_data <- merge(fake_data, period_map)
      df <- rbind(df, fake_data)
    }
  }
}

df[diagnostic == "ict", diagnostic := "ICT"] ### Standardize on "ICT"

## Create SPDE INLA stack
input_data <- build_mbg_data_stack(df            = df, # note that merge (if using TMB) will return data in a different (but internally consistent) order, just different than df
                                   fixed_effects = all_fixed_effects,
                                   spde_prior    = spde_prior,
                                   mesh_s        = mesh_s,
                                   mesh_s_ppp    = mesh_s_ppp,
                                   tmb = fit_with_tmb,
                                   mesh_t        = mesh_t, # not currently implemented with tmb
                                   use_ctry_res  = as.logical(use_inla_country_res),
                                   use_nugget    = use_inla_nugget, # implemented with tmb
                                   exclude_cs    = child_model_names, # raw covs will get center scaled here though (see notes above)
                                   coefs.sum1    = coefs_sum1, # not currently implemented tmb
                                   scale_gaussian_variance_N = scale_gaussian_variance_N,
                                   shapefile_version = modeling_shapefile_version,
                                   zl            = z_list, # if this is not zero and tmb==TRUE, it will trigger 3rd kronecker and fixed effects
                                   zcol          = zcol,   # must not be null if z_list is present
                                   use_gp = use_gp,
                                   use_space_only_gp = use_space_only_gp,
                                   use_time_only_gmrf = use_time_only_gmrf,
                                   use_timebyctry_res = use_timebyctry_res,
                                   adm0_list = gaul_list,
                                   use_pref_samp_pp = as.logical(use_pref_samp_pp),
                                   pref_samp_pp_stages = pref_samp_pp_stages,
                                   covars_in_pp = covars_in_pp,
                                   country_re_pp = country_re_pp,
                                   subset = subset_shape,
                                   subset_ppp = subset_shape_ppp,
                                   gaul_list_ppp = gaul_list_ppp,
                                   fit_crosswalk_inla = fit_crosswalk_inla,
                                   st_gp_int_zero = as.logical(st_gp_int_zero),
                                   s_gp_int_zero = as.logical(s_gp_int_zero),
                                   crosswalk_training_data_set = crosswalk_training_data_set,
                                   use_subnat_res = use_subnat_res)

## Generate other inputs necessary
outcome <- df[[indicator]] # N+_i - event obs in cluster
N <- df$N # N_i - total obs in cluster

if (!as.logical(fit_with_tmb)) {
  ## combine all the inputs
  stacked_input <- input_data[[1]]
  spde <- input_data[[2]]
  cs_df <- input_data[[3]]
  spde.sp <- input_data[[4]] ## used for space only (time stationary) gp
  spde2 <- input_data[[5]]
  
  weights <- inla.stack.data(stacked_input)$weight
} else { # TMB
  weights <- df$weight
}

# catch in case there is no weight column
if (is.null(weights)) {
  weights <- rep(1, nrow(df))
}

## Set the number of cores to be equal to input; f missing, then revert to cores_to_use value
if(Sys.getenv("OMP_NUM_THREADS") != "") {
  setompthreads(Sys.getenv("OMP_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input OpenMP threads.")
  setompthreads(cores_to_use)
}

if(Sys.getenv("MKL_NUM_THREADS") != "") {
  setmklthreads(Sys.getenv("MKL_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input MKL threads.")
  setmklthreads(cores_to_use)
}

## Fit MBG model
if (!as.logical(skipinla)) {
  if (!as.logical(fit_with_tmb)) {
    message('Fitting model with R-INLA')
    
    if (!exists("omp_strat")) {
      omp_strat <- "default"
      pardiso_license <- NULL
    } else if (grepl("pardiso", omp_strat)) {
      inla.setOption("pardiso.license", <<<< FILEPATH REDACTED >>>>)
      pardiso_license <- <<<< FILEPATH REDACTED >>>>
    } else {
      pardiso_license <- NULL
    }
    
    model_fit <- fit_mbg(indicator_family = indicator_family,
                         stack.obs        = stacked_input,
                         spde             = spde,
                         cov              = outcome,
                         N                = N,
                         int_prior_mn     = intercept_prior,
                         f_mbg            = mbg_formula,
                         run_date         = run_date,
                         keep_inla_files  = keep_inla_files,
                         cores            = as.integer(fthread_parallel),
                         blas_cores       = as.integer(fthread_parallel),
                         omp_strat        = omp_strat,
                         pardiso_license  = pardiso_license,
                         wgts             = weights,
                         intstrat         = intstrat,
                         fe_sd_prior      = 1 / 9,
                         verbose_output   = TRUE,
                         use_pref_samp_pp = as.logical(use_pref_samp_pp),
                         n = nrow(df),
                         nv = mesh_s$n,
                         pref_samp_pp_stages = pref_samp_pp_stages,
                         fit_crosswalk_inla = fit_crosswalk_inla,
                         perform_initial_INLA_approx = perform_initial_INLA_approx,
                         initial_INLA_approx_levels = initial_INLA_approx_levels)
    
  } else {
    message('Fitting model with TMB')
    message(sprintf('%i Data points and %i mesh nodes',nrow(df),length(input_data$Parameters$Epsilon_stz)))
    
    ## Save RDS file of input data for replication
    saveRDS(object = input_data, ## save this here in case predict dies
            file = sprintf(<<<< FILEPATH REDACTED >>>>))
    ## Run the model
    system.time(
      model_fit <- fit_mbg_tmb( lbdcorerepo     = core_repo,
                                cpp_template    = 'mbg_tmb_model',
                                tmb_input_stack = input_data,
                                control_list    = list(trace=1, eval.max=500, iter.max=300, abs.tol=1e-20),
                                optimizer       = 'nlminb', # TODO add optimx
                                ADmap_list      =  NULL,
                                sparse_ordering = as.logical(sparse_ordering))
    )
    
    ## Clamping
    clamp_covs <- TRUE # TODO CONFIG THIS
  }
  
  saveRDS(object = model_fit, ## save this here in case predict dies
          file = sprintf(<<<< FILEPATH REDACTED >>>>))
} else {
  ## Skipped fitting INLA so just load model and move to predict
  model_fit <- readRDS( file = sprintf(<<<< FILEPATH REDACTED >>>>))
}

# Run predict_mbg on chunks of X samples (to avoid memory issues)
if (exists("max_chunk")) {
  max_chunk <- as.integer(max_chunk)
} else {
  max_chunk <- 50
}

for (l in 1:length(cov_list)) {
  message(sprintf("On cov %i out of %i", l, length(cov_list)))
  cov_list[[l]] <- raster::crop(cov_list[[l]], raster::extent(simple_raster))
  cov_list[[l]] <- raster::setExtent(cov_list[[l]], simple_raster)
  cov_list[[l]] <- raster::mask(cov_list[[l]], simple_raster)
}

#### Save model and training data set
save(model_fit, file = (paste0(<<<< FILEPATH REDACTED >>>>)))
write.csv(df, file = (paste0(<<<< FILEPATH REDACTED >>>>)), row.names = FALSE)

samples <- as.numeric(samples)

if (!as.logical(fit_with_tmb)) {
  #### Sample posterior draws for linear predictor
  message(paste0("Saving draws from inla.posterior.sample()"))
  
  suppressWarnings(draws <- inla.posterior.sample(samples, model_fit, selection = list("APredictor" = 1:nrow(df))))
  input <- stacked_input$data$data
  
  for (i in 1:length(draws)) {
    draws[[i]] <- draws[[i]]$latent
  }
  
  pred_draws <- do.call(cbind, draws)
  
  ## Save pred_draws to disk
  outputdir_new <- sprintf(<<<< FILEPATH REDACTED >>>>) # new outputdir
  dir.create(outputdir_new)
  save(pred_draws, file = paste0(<<<< FILEPATH REDACTED >>>>))
  save(input, file = paste0(<<<< FILEPATH REDACTED >>>>))
  write.csv(df, file = paste0(<<<< FILEPATH REDACTED >>>>))
}

#### Run full pixel-level predictions unless only INLA predictions are requested (i.e., for OOS validation)
if ((holdout == 0) | !oos_preds_inla) {
  if (parallel_pred_agg & !as.logical(fit_with_tmb)) { # Run predict_mbg() in parallel (currently only for INLA)
    if (!file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
      ### Generate and save posterior draws from INLA model fit
      save_all_draws(model_fit = model_fit, samples = samples, rd = run_date, outputdir = outputdir)
    } else {
      message("Draws object already present on disk")
    }
    
    message(paste0("Making predictions in ", max_chunk, " draw chunks, in parallel by predict_years."))
    
    if (is.null(restart_predict_years)) {
      years <- predict_years
    } else {
      years <- restart_predict_years
    }
    
    if (!anyNA(restart_predict_years)) {
      ### Next, launch prediction script in parallel
      predict_job <- parallelize(slots            = fthread_predict,
                                 memory           = m_mem_free_predict,
                                 script           = predict_code,
                                 geo_nodes        = as.logical(use_geos_nodes),
                                 expand_vars      = list(predict_years = unlist(years)),
                                 save_objs        = c("region", "indicator", "indicator_group", "run_date", "config_file", "df", "simple_raster", "subset_shape", "simple_polygon", "pop_raster", "cov_list", "core_repo", "indic_repo", "holdout", "model_fit", "cs_df", "mesh_s", "mesh_t", "child_model_names", "all_fixed_effects", "period_map_predict", "stacked_input", "simple_raster2"),
                                 prefix           = "predict",
                                 log_location     = 'sharedir',
                                 script_dir       = paste0(<<<< FILEPATH REDACTED >>>>),
                                 run_time         = runtime_predict,
                                 threads          = fthread_predict,
                                 singularity_opts = list(SET_OMP_THREADS = fthread_predict, SET_MKL_THREADS = fthread_predict))
      
      monitor_jobs(predict_job)
    }
    
    save_mbg_preds_parallel(
      config = config,
      time_stamp = time_stamp,
      run_date = run_date,
      res_fit = model_fit,
      df = df,
      pathaddin = pathaddin
    )
    
  } else { # Classical approach to predict_mbg()
    
    message(paste0("Making predictions in ", max_chunk, " draw chunks."))
    
    ## Create vector of chunk sizes
    chunks <- rep(max_chunk, samples %/% max_chunk)
    if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
    c <- 0
    pm <- lapply(chunks, function(samp) {
      c <<- c + 1
      message(paste0("Starting prediction batch ", c, " of ", length(chunks), "..."))
      if(fit_with_tmb == FALSE){
        predict_mbg(res_fit       = model_fit,
                    cs_df         = cs_df,
                    mesh_s        = mesh_s,
                    mesh_t        = mesh_t,
                    cov_list      = cov_list,
                    samples       = samp,
                    simple_raster = simple_raster,
                    transform     = transform,
                    coefs.sum1    = coefs_sum1,
                    pred_gp       = as.logical(use_gp),
                    shapefile_version = modeling_shapefile_version,
                    predict_years = predict_years,
                    return_mean_sd = FALSE,
                    use_space_only_gp = use_space_only_gp,
                    use_time_only_gmrf = use_time_only_gmrf,
                    use_timebyctry_res = as.logical(use_timebyctry_res),
                    fit_crosswalk_inla = fit_crosswalk_inla,
                    predict_diagnostic = predict_diagnostic,
                    predict_age_start = predict_age_start,
                    predict_age_end = predict_age_end,
                    crosswalk_training_data_set = crosswalk_training_data_set,
                    simple_raster_subnats = simple_raster2)
      } else {
        predict_mbg_tmb(samples              = samp,
                        seed                 = NULL,
                        tmb_input_stack      = input_data,
                        model_fit_object     = model_fit,
                        fes                  = all_fixed_effects,
                        sr                   = simple_raster,
                        yl                   = year_list,
                        zl                   = z_list,
                        covs_list            = cov_list,
                        clamp_covs           = clamp_covs,
                        cov_constraints = covariate_constraint_vectorize(config),
                        use_full_interacting_effect = as.logical(use_gp),
                        use_space_only_gp = as.logical(use_space_only_gp),
                        use_time_only_gmrf = as.logical(use_time_only_gmrf),
                        use_age_only_gmrf = as.logical(use_age_only_gmrf),
                        coefs.sum1           = coefs_sum)
      }
    })
    
    if(fit_with_tmb == FALSE){
      preds <- vector("list", length(pm))
      spatial_random_effects <- vector("list", length(pm))
      for (i in 1:length(pm)) {
        preds[[i]] <- pm[[i]][[3]]
        spatial_random_effects[[i]] <- pm[[i]][[4]]
      }
      
      pm <- preds
    } else {
      preds <- vector("list", length(pm))
      spatial_random_effects <- vector("list", length(pm))
      for (i in 1:length(pm)) {
        preds[[i]] <- pm[[i]][[1]]
        spatial_random_effects[[i]] <- pm[[i]][[2]]
      }
      
      pm <- preds
    }
    
    
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # if z dimension has more than one level, then save each z as a different indicator
    if (length(z_list) > 1) {
      
      # reorder pm list, right now it's z within each chunk. rbind all z's together
      for (z in z_list) { # z_list must be integers starting with 1
        if (length(chunks) > 1) {
          for (ch in 2:length(chunks)) {
            pm[[1]][[z]] <- cbind(pm[[1]][[z]], pm[[ch]][[z]])
          }
        }
      }
      pm <- pm[[1]] # pm is now a list of cell_preds by z
      
      # loop over z and save as an indicator each one
      orig_indic <- indicator
      orig_paddin <- pathaddin
      orig_outdir <- outputdir
      
      message("Wrapping up")
      
      for (z in z_list) {
        indicator <- sprintf("%s_%s%i", orig_indic, zcol, z) # new indicator name
        pathaddin <- paste0("_bin", z, "_", reg, "_", holdout) # new pathaddin
        outputdir <- sprintf(<<<< FILEPATH REDACTED >>>>) # new outputdir
        dir.create(outputdir)
        message(sprintf("New indicator: %s", indicator))
        
        # make a mean raster
        mean_ras <- insertRaster(simple_raster, matrix(rowMeans(pm[[z]]), ncol = length(predict_years)))
        
        # save indicator stuff
        save_mbg_preds(
          config = config,
          time_stamp = time_stamp,
          run_date = run_date,
          mean_ras = mean_ras,
          sd_ras = NULL,
          res_fit = model_fit,
          cell_pred = pm[[z]],
          df = df,
          pathaddin = pathaddin
        )
        
        pdf(paste0(<<<< FILEPATH REDACTED >>>>))
        plot(mean_ras, maxpixel = 1e6)
        dev.off()
      }
      
      indicator <- orig_indic
      pathaddin <- orig_paddin
      outputdir <- orig_outdir
      
      message("done saving indicator-specific outputs by z")
    } else { # if no z colums (most peoples cases)
      
      ## Make cell preds and a mean raster
      cell_pred <- do.call(cbind, pm)
      mean_ras <- insertRaster(simple_raster, matrix(rowMeans(cell_pred), ncol = length(predict_years)))
      
      message("Wrapping up")
      save_mbg_preds(
        config = config,
        time_stamp = time_stamp,
        run_date = run_date,
        mean_ras = mean_ras,
        sd_ras = NULL,
        res_fit = model_fit,
        cell_pred = cell_pred,
        df = df,
        pathaddin = pathaddin
      )
      
      # plot the mean raster
      cols <- c("#810f7c", "#9ebcda", "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#800026")
      cut <- c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 1)
      cuts <- c(seq(0, .3, .005))
      cuts2 <- c(seq(0, 1, .05))
      
      names(mean_ras) <- paste0("mean_pred_year.", predict_years)
      
      outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)
      pdf(paste0(<<<< FILEPATH REDACTED >>>>))
      
      for (i in 1:length(y_start)) {
        plot(mean_ras[[y_start[i]:y_end[i]]], maxpixel = 1e6, breaks = cut, col = cols)
      }
    
      dev.off()
      
      ###### Output spatial random effects
      if(!exists("return_spatial_random_fields")) {
        return_spatial_random_fields <- TRUE
      }
      
      if(as.logical(return_spatial_random_fields)) {
        cell_pred <- do.call(cbind, spatial_random_effects)
        mean_ras <- insertRaster(simple_raster, matrix(rowMeans(cell_pred), ncol = length(predict_years)))
        
        outputdir <- sprintf(<<<< FILEPATH REDACTED >>>>)
        pdf(paste0(<<<< FILEPATH REDACTED >>>>))
        plot(mean_ras, maxpixel = 1e6)
        dev.off()
        
        writeRaster(mean_ras, file = (paste0(<<<< FILEPATH REDACTED >>>>)), overwrite = TRUE, format = "GTiff")
      }
      
    }
  }
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
