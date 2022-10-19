# --------------------------------------------------------------------------------------------------
# Script designed to facilitate covariate selection using variance inflation factor (VIF) analysis
# The purpose of the VIF analysis is to test for multicollinearity in covariates
#
# (1) load covariate data cropped to MBG input data by region
# (2) perform VIF analysis and drop covariates above a threshold by region
# (3) take results and summarize in a table that is human-readable
# (4) take results and summarize in MBG input format
#
#
# References:
# Faraway's Linear Models in R, Chapter 4: Diagnostics
# https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/

# NOTE: fmsb package was built under R version 3.6.0. The package will fail to load if using the latest singularity image.
# --------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
# Setup

user <- Sys.info()["user"]
repo <- <<<< FILEPATH REDACTED >>>>
indicator_group <- "oncho"
indicator <- "had_oncho_w_resamp"
config_par <- "config"
config_file <- config_par
cov_par <- "covs"
cov_file <- config_file
# indicate whether to use old run date
use_old_run_date <- FALSE
old_run_date_input <- ""
regions <- "oncho_endem_afr"
threshold_min <- 2
threshold_max <- 10
threshold_step <- 1
individual_countries <- FALSE
crop_covs <- TRUE
cropped_covs_file <- NULL

# set additional arguments
core_repo <- repo
age <- 0
holdout <- 0
threshold <- c(seq(as.numeric(threshold_min), as.numeric(threshold_max), by = as.numeric(threshold_step)))

# Load MBG packages and functions
message("Loading in required R packages and MBG functions")
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>)))
package_list <- c(package_list, "sf")

source(<<<< FILEPATH REDACTED >>>>)
mbg_setup(package_list = package_list, repos = core_repo)

### use newer version of raster package and load fasterize and sf
path <- <<<< FILEPATH REDACTED >>>>
library(fasterize, lib.loc = path)
library(sf, lib.loc = path)
library(fmsb, lib.loc = path)

indic_repo <- <<<< FILEPATH REDACTED >>>>
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

# Load and check config
config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator, 
                                config_file = <<<< FILEPATH REDACTED >>>>, run_tests = FALSE)
message(paste0("covs_name = ", covs_name))

# set run date
if (use_old_run_date == FALSE) {
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- old_run_date_input
}

# create directories to save results in
save_dir <- <<<< FILEPATH REDACTED >>>>
dir.create(save_dir)
mbg_input_dir <- <<<< FILEPATH REDACTED >>>>
dir.create(mbg_input_dir)

## Create proper year list object
if (class(year_list) == "character") year_list <- eval(parse(text = year_list))

## Throw a check for things that are going to be needed later
message("Looking for things in the config that will be needed for this script to run properly")
check_config()

# If running individual countries make sure all country FEs and REs off
if (individual_countries) {
  use_child_country_fes <- FALSE
  use_inla_country_fes <- FALSE
  use_inla_country_res <- FALSE
}
# ------------------------------------------------------------------------------------------------------------------------

if (!exists("rw1_raw_covar_list")) {
  rw1_raw_covar_list <- NULL
} else {
  rw1_raw_covar_list <- eval(parse(text = rw1_raw_covar_list))
}

# ------------------------------------------------------------------------------------
# Crop covariates to data

if (crop_covs == TRUE) {
  # set up blank list to store cropped covariate data
  all_cov_data <- list()
  
  # start loop over regions
  for (region in regions) {
    
    
    # --------------------------------------------------------------------------------------------------
    # Prep MBG inputs and load data
    
    # report region and indicator
    message(indicator)
    reg <- region
    message(reg)
    
    # make pathaddin
    pathaddin <- <<<< FILEPATH REDACTED >>>>
    
    ## some set up
    if (class(year_list) == "character") year_list <- eval(parse(text = year_list))
    if (class(z_list) == "character") z_list <- eval(parse(text = z_list))
    
    ## Load simple polygon template to model over
    gaul_list <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = use_premade, shapefile_version = modeling_shapefile_version)
    subset_shape <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
    
    ## Load list of raster inputs (pop and simple)
    raster_list <- build_simple_raster_pop(subset_shape)
    simple_raster <- raster_list[["simple_raster"]]
    pop_raster <- raster_list[["pop_raster"]]
    
    simple_raster <- raster::mask(simple_raster, pop_raster[[1]])
    
    ## Load input data and make sure it is cropped to modeling area
    df <- load_input_data(
      indicator = indicator,
      simple = simple_polygon,
      removeyemen = FALSE,
      years = "annual",
      update_run_date = FALSE,
      use_share = as.logical(use_share),
      withtag = as.logical(withtag),
      datatag = datatag,
      yl = year_list
    )
    
    # remove unecessary folder and data file
    unlink(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)
    
    ## Some built in data checks that cause known problems later on
    if (indicator_family == "binomial" & any(df[, get(indicator)] / df$N > 1)) {
      stop("You have binomial data where k > N. Check your data before proceeding")
    }
    if (any(df[["weight"]] %in% c(Inf, -Inf) | any(is.na(df[["weight"]])))) {
      stop("You have illegal weights (NA,Inf,-Inf). Check your data before proceeding")
    }
    # --------------------------------------------------------------------------------------------------
    
    
    # --------------------------------------------------------------------------------------------------------------
    # Pull covariates
    
    ## Define modeling space. In years only for now.
    if (yearload == "annual") {
      period_map <-
        make_period_map(modeling_periods = c(min(year_list):max(year_list)))
    }
    if (yearload == "five-year") {
      period_map <-
        make_period_map(modeling_periods = seq(min(year_list), max(year_list), by = 5))
    }
    
    ## Make placeholders for covariates
    cov_layers <- gbd_cov_layers <- NULL
    
    ## Pull all covariate bricks/layers
    if (nchar(fixed_effects) > 0) {
      message("Grabbing raster covariate layers")
      
      cov_layers <- load_lagged_covariates(covariate_config = fixed_effects_config,
                                           template = simple_raster,
                                           start_year = min(year_list),
                                           end_year = max(year_list),
                                           raster_agg = as.integer(raster_agg_factor))
      
    }
    
    if (as.logical(use_river_size)) {
      river_size <- raster(<<<< FILEPATH REDACTED >>>>)
      river_size <- projectRaster(river_size, crs = crs(cov_layers[[1]][[1]]))
      river_size <- raster::resample(river_size, cov_layers[[1]][[1]])
      river_size <- raster::crop(river_size, cov_layers[[1]][[1]])
      river_size[is.na(river_size[])] <- 0
      river_size <- raster::mask(river_size, cov_layers[[1]][[1]])
      
      cov_layers[["river_size"]] <- river_size
    }
    
    if (as.logical(use_oncho_suitability)) {
      oncho_suitability <- stack(<<<< FILEPATH REDACTED >>>>)
      oncho_suitability <- subset(oncho_suitability, 1)
      oncho_suitability <- projectRaster(oncho_suitability, crs = crs(cov_layers[[1]][[1]]))
      oncho_suitability <- raster::resample(oncho_suitability, cov_layers[[1]][[1]])
      oncho_suitability <- raster::crop(oncho_suitability, cov_layers[[1]][[1]])
      oncho_suitability[is.na(oncho_suitability[])] <- 0
      oncho_suitability <- raster::mask(oncho_suitability, cov_layers[[1]][[1]])
      
      cov_layers[["oncho_suitability"]] <- oncho_suitability
    }
    
    ## Pull country level gbd covariates
    if (nchar(gbd_fixed_effects) > 0) {
      message("Grabbing GBD covariates")
      
      effects <- trim(strsplit(gbd_fixed_effects, "\\+")[[1]])
      measures <- trim(strsplit(gbd_fixed_effects_measures, "\\+")[[1]])
      gbd_cov_layers <- load_gbd_covariates(
        covs = effects,
        measures = measures,
        year_ids = year_list,
        age_ids = gbd_fixed_effects_age,
        template = cov_layers[[1]][[1]],
        simple_polygon = simple_polygon,
        interval_mo = interval_mo
      )
    }
    
    ## Combine all covariates
    all_cov_layers <- c(cov_layers, gbd_cov_layers)
    
    ## regenerate all fixed effects equation from the cov layers
    all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")
    
    ## Make stacker-specific formulas where applicable
    all_fixed_effects_brt <- all_fixed_effects
    
    # get cov names
    the_covs <- format_covariates(all_fixed_effects)
    
    ## copy the dataset to avoid unintended namespace conflicts
    the_data <- copy(df)
    
    ## add a row id column
    the_data[, a_rowid := seq(1:nrow(the_data))]
    
    ## extract covariates to the points and subset data where its missing covariate values
    cs_covs <- extract_covariates(the_data,
                                  all_cov_layers,
                                  id_col = "a_rowid",
                                  return_only_results = TRUE,
                                  centre_scale = TRUE,
                                  period_var = "year",
                                  period_map = period_map
    )
    
    
    ## Check for data where covariate extraction failed
    rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
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
    
    the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)
    
    ## store the centre scaling mapping
    covs_cs_df <- cs_covs[[2]]
        
    # this will drop rows with NA covariate values
    the_data_na <- na.omit(the_data, c(indicator, "N", the_covs), invert = TRUE)
    the_data <- na.omit(the_data, c(indicator, "N", the_covs))
    
    ## stop if this na omit demolished the whole dataset
    if (nrow(the_data) == 0) stop("You have an empty df, make sure one of your covariates was not NA everywhere.")
    
    ## make sure we're working with a data table
    the_data <- data.table(the_data)
    
    ## create a data table with just covariates cropped by locations in region with data
    cov_cols <- names(all_cov_layers)
    cov_data <- the_data[, cov_cols, with = FALSE]
    cov_data <- cov_data[, lat := NULL]
    cov_data <- cov_data[, lon := NULL]
    
    
    ## save data to list
    all_cov_data[[region]] <- cov_data
    # --------------------------------------------------------------------------------------------------------------
  } # end loop over regions
  
  
  # save data (for easier loading in the future)
  cropped_covs_file <- paste0(save_dir, "covs_cropped_to_data_", Sys.Date(), ".RData")
  save(all_cov_data, file = cropped_covs_file)
  # ------------------------------------------------------------------------------------
} else {
  # load previously cropped covariate data
  load(cropped_covs_file)
}


# ---------------------------------------------------------------
# Stepwise VIF selection function

# Start function
stepwise_vif_selection <- function(thresh,
                                   covariate_data,
                                   reg,
                                   trace = TRUE) {
  # ---------------------------------------------------------------
  
  
  # ------------------------------------------------------------------------------------------------------
  # Data checks
  
  # make sure input data is a data table
  if (any(!"data.table" %in% class(covariate_data))) covariate_data <- data.table(covariate_data)
  
  # remove any covariate values that do not vary accross the region
  var_names <- names(covariate_data)
  for (val in var_names) {
    if (!is.na(val) & length(unique(covariate_data[[val]])) < 2) {
      message(paste0('\nRemoving covariate "', val, '" because it is constant across the region.'))
      covariate_data[, (val) := NULL]
    }
  }
  # ------------------------------------------------------------------------------------------------------
  
  
  # ------------------------------------------------------------------------------------------------------------------------
  # VIF selection
  
  # load vif analysis package
  user <- Sys.info()[["user"]]
  path <- <<<< FILEPATH REDACTED >>>>
  library(fmsb, lib.loc = path)
  
  # get initial vif value for all comparisons of variables
  vif_init <- NULL
  
  for (val in var_names) {
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = "+")
    form_in <- formula(paste(val, "~", form))
    vif_init <- rbind(vif_init, c(val, VIF(lm(form_in, data = covariate_data))))
  }
  
  vif_max <- max(as.numeric(vif_init[, 2]), na.rm = TRUE)
  
  if (vif_max < thresh) {
    if (trace == T) { # print output of each iteration
      prmatrix(vif_init, collab = c("var", "vif"), rowlab = rep("", nrow(vif_init)), quote = F)
      cat("
")
      
      message(paste("
All variables have VIF < ", thresh, ", max VIF ", round(vif_max, 2), sep = ""), "

")
    }
    
    cat(var_names)
  } else {
    
    # backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while (vif_max >= thresh) {
      vif_vals <- NULL
      var_names <- names(covariate_data)
      
      for (val in var_names) {
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = "+")
        form_in <- formula(paste(val, "~", form))
        vif_add <- VIF(lm(form_in, data = covariate_data))
        vif_vals <- rbind(vif_vals, c(val, vif_add))
      }
      
      max_row <- which(vif_vals[, 2] == max(as.numeric(vif_vals[, 2]), na.rm = TRUE))[1]
      
      vif_max <- as.numeric(vif_vals[max_row, 2])
      
      if (vif_max < thresh) break
      
      if (trace == T) { # print output of each iteration
        prmatrix(vif_vals, collab = c("var", "vif"), rowlab = rep("", nrow(vif_vals)), quote = F)
        cat("
")
        message(paste0("removed: ", vif_vals[max_row, 1], " ", vif_max, "

"))
        flush.console()
      }
      
      keep_cols <- var_names[-which(var_names == vif_vals[max_row, 1])]
      covariate_data <- covariate_data[, keep_cols, with = FALSE]
    }
    
    message(paste0("
Covariates selected for ", indicator, " in ", reg, " with a threshold of ", thresh, ":"))
    cat(names(covariate_data))
  }
  
  # vector of selected covariates
  selected_covs <- c(names(covariate_data))
  # ------------------------------------------------------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------------------
  # Add selection results to table
  
  # add 0/1 for whether or not covariate was selected to a data table
  d1 <- data.table(cov_names)
  d2 <- data.table(selected_covs)
  d2[, included := TRUE]
  dt <- merge(d1, d2, by.x = "cov_names", by.y = "selected_covs", all = T)
  dt[is.na(included), included := FALSE]
  dt <- dcast(melt(dt, id.vars = "cov_names"), variable ~ cov_names)
  dt[, variable := NULL]
  
  # add threshold, region, number of covariates selected, and run date that corresponds to data
  dt[, region := reg]
  dt[, vif_threshold := thresh]
  dt[, num_covs_selected := length(selected_covs)]
  dt[, data_run_date := run_date]
  
  # bind to original data table
  selection_results <- rbindlist(list(selection_results, dt), use.names = TRUE, fill = T)
  # --------------------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------------------------------------
  # Save mbg input files
  
  # create mbg input files
  selected_covs <- data.table(covariate = selected_covs)
  cov_config <- merge(covs, selected_covs, by = "covariate")
  
  # save mbg input files
  dir.create(paste0(mbg_input_dir, "vif_", i, "/"), showWarnings = FALSE)
  write.csv(cov_config, <<<< FILEPATH REDACTED >>>>)
  # ----------------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------
  # End function
  
  # return a vector of covariates selected
  return(selection_results)
}
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
# Stepwise covariate removal using VIF

# loop over thresholds
for (i in threshold) {
  
  # report threshold
  message(paste0("
Threshold: ", i))
  
  # load covariate data
  load(cropped_covs_file)
  
  # load input covariate names
  # covs <- fread(paste0(indic_repo, indicator_group, "/", cov_par, '.csv'))
  covs <- fread(paste0(indic_repo, "/covs/", cov_par, ".csv"))
  covs <- covs[!(covariate %in% c("lat", "lon"))]

  if (as.logical(use_river_size)) {
    covs <- rbind(covs, data.table("covariate" = "river_width", "measure" = "mean", "gbd" = "FALSE", "include" = "TRUE", "year_lag" = 0, "release" = ""))
  }
  
  cov_names <- covs[include == "TRUE", covariate]
  
  # cov_names <- cov_names[!(cov_names %in% c("lat", "lon"))]
  
  # set up table
  col_names <- c(cov_names, "region", "vif_threshold", "num_covs_selected", "data_run_date")
  selection_results <- setNames(data.table(matrix(
    nrow = 0,
    ncol = length(col_names)
  )), col_names)
  
  # run vif selection analysis
  for (region in regions) {
    message(region)
    
    # run stepwise vif selection function
    selection_results <- stepwise_vif_selection(
      thresh = i,
      covariate_data = all_cov_data[[region]],
      reg = region,
      trace = F
    )
  }
  
  # save selection results
  write.csv(selection_results, <<<< FILEPATH REDACTED >>>>)
  
  # reset for next threshold
  message("Selection complete. Results saved.")
}
# ---------------------------------------------------------------------------------------------------------------
