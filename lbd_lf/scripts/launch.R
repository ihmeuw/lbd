#####################################################################################################################################
## Launch Focal 3 MBG models
## This script preps data and spatial objects, launches parallel model scripts
## and performs post-prediction.
#####################################################################################################################################

#### Load objects from parallelize
source(<<<< FILEPATH REDACTED >>>>)
load_from_parallelize()
message(paste0("Using ", core_repo))

user <- Sys.info()[["user"]]

## Set repo locations (core_repo is passed from parent_master.R but is retained here in commented form to facilitate interactive debugging)
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>)

path <- paste0(<<<< FILEPATH REDACTED >>>>)

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- paste0(<<<< FILEPATH REDACTED >>>>)
package_list <- c(t(read.csv(sprintf(<<<< FILEPATH REDACTED >>>>), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(paste0(<<<< FILEPATH REDACTED >>>>))

mbg_setup(package_list = package_list, repos = core_repo)

## Load additional needed tables
library(data.table)
library(fasterize)
library(matrixStats)

## Focal 3-specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(paste0(<<<< FILEPATH REDACTED >>>>), recursive = TRUE)) {
  message(funk)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}

## Load config; set_up_config_focal_3() is a wrapper for set_up_config() which accommodates custom covariate options
config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator, 
                                config_file = paste0(<<<< FILEPATH REDACTED >>>>), run_tests = FALSE)

## Source Focal 3-specific setup script to prepare custom objects
source(paste0(<<<< FILEPATH REDACTED >>>>))

## Process year objects
if (class(year_list) == "character") year_list <- eval(parse(text = year_list))
if (class(predict_years) == "character") predict_years <- eval(parse(text = predict_years))
if (exists("aggregation_years")) {
  if (class(aggregation_years) == "character") aggregation_years <- eval(parse(text = aggregation_years))  
} else {
  aggregation_years <- copy(predict_years)
}

## Additional setup
sharedir <- sprintf(<<<< FILEPATH REDACTED >>>>)
region <- unlist(strsplit(region_list, " "))
Regions <- region

## Skip to post-estimation code or prepare for model launch
if (skip_to_aggregation) { # Skip model launching and go directly to post-estimation (e.g., if the model has already been fit but post-estimation failed)
  message("Skipping to post-estimation code")
  
  ## Make loopvars aka strata grid (format = regions, ages, holdouts)
  if (as.logical(makeholdouts)) loopvars <- expand.grid(region, 0, 1:n_ho_folds) else loopvars <- expand.grid(region, 0, 0)
  source(paste0(<<<< FILEPATH REDACTED >>>>))
  
} else {
  ## Make loopvars aka strata grid (format = regions, ages, holdouts)
  if (as.logical(makeholdouts)) loopvars <- expand.grid(region, 0, 1:n_ho_folds) else loopvars <- expand.grid(region, 0, 0)
  
  ## Skip launching the model(s) and go directly to creation of summary rasters? (e.g., if prediction finished but rasters were not produced)
  if (!skip_to_summary_rasters) {
    ## Loop through regions
    for (r in region) {
      ## Retrieve list of country-level administrative codes
      gaul_list <- get_adm0_codes(r, shapefile_version = modeling_shapefile_version)
      
      ## Check for existing spatial object file; load from disk if it is present, otherwise derive spatial objects
      if (!(file.exists(paste0(<<<< FILEPATH REDACTED >>>>)))) {
        ## Load simple polygon template to model over
        simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = use_premade, shapefile_version = modeling_shapefile_version)
        subset_shape <- simple_polygon_list[[1]]
        simple_polygon <- simple_polygon_list[[2]]
        
        ## Build administrative and population rasters
        raster_list <- build_simple_raster_pop(subset_shape)
        
        ## Create a full-resolution admin raster, masking by population raster (to remove pixels with missing population data)
        simple_raster_raf1 <- raster::mask(raster_list[["simple_raster"]], raster_list[["pop_raster"]][[1]])
        
        ## Aggregate population and administrative rasters by raster_agg_factor (default is 1, in which case aggregation is not performed)
        pop_raster <- suppressWarnings(raster::aggregate(raster_list[["pop_raster"]], fun = sum, fact = as.integer(raster_agg_factor)))
        simple_raster <- suppressWarnings(raster::resample(raster_list[["simple_raster"]], pop_raster[[1]], method="ngb"))
        
        ## Mask aggregated administrative raster by aggregated population raster
        simple_raster <- raster::mask(simple_raster, pop_raster[[1]])
        
        ## Set up spatial objects for preferential sampling, if requested
        if (as.logical(use_pref_samp_pp) & !is.null(ppp_region)) {
          ## Load simple polygon template to model over
          gaul_list_ppp <- get_adm0_codes(ppp_region, shapefile_version = modeling_shapefile_version)
          simple_polygon_list_ppp <- load_simple_polygon(gaul_list = gaul_list_ppp, buffer = 1, tolerance = 0.4, use_premade = use_premade, shapefile_version = modeling_shapefile_version)
          subset_shape_ppp <- simple_polygon_list_ppp[[1]]
        } else {
          gaul_list_ppp <- ""
          subset_shape_ppp <- subset_shape
        }
        
        ## Save spatial templates to disk
        save(simple_raster, simple_raster_raf1, subset_shape, simple_polygon, pop_raster, subset_shape_ppp, file = paste0(<<<< FILEPATH REDACTED >>>>))
      } else {
        ## Load existing spatial templates from disk
        load(paste0(<<<< FILEPATH REDACTED >>>>))
      }
    }
    
    ## Load input data
    df <- load_input_data(indicator   = indicator,
                          simple      = simple_polygon,
                          removeyemen = FALSE,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list,
                          region      = reg)
    
    ## Create df with pixel ids (and polygon centroids and age group weights if aggregation)
    df <- process_input_data(df, 
                             pop_release                = pop_release, 
                             interval_mo                = interval_mo, 
                             modeling_shapefile_version = modeling_shapefile_version,
                             poly_ag                    = as.logical(poly_ag), 
                             zcol                       = zcol, 
                             zcol_ag                    = zcol_ag, 
                             zcol_ag_id                 = zcol_ag_id, 
                             z_map                      = z_map, 
                             z_ag_mat                   = z_ag_mat)
    
    ### Check for incorrect ADM0
    # Load link table for GADM codes
    adm_link <- fread(<<<< FILEPATH REDACTED >>>>)
    
    ## Extract ADM0 using shapefile rather than raster
    adm0 <- over(SpatialPoints(df[, c("longitude", "latitude")], proj4string=CRS(proj4string(subset_shape))), subset_shape)
    setnames(adm0, "ADM0_CODE", "ADM0_CODE_over")
    df <- cbind(df, adm0)
    df <- merge(df, adm_link[, c("gadm_geoid", "iso3")], by.x = "ADM0_CODE_over", by.y = "gadm_geoid", all.x = TRUE)
    
    wrong_adm0 <- df[country != iso3]
    message(paste0(nrow(wrong_adm0), " rows dropped due to ADM0 mismatch. These have been saved in wrong_adm0.csv. Please check these."))
    write.csv(wrong_adm0, paste0(<<<< FILEPATH REDACTED >>>>))
    
    df <- df[country == iso3]
    
    ###############################################################################
    ## Make Holdouts
    ###############################################################################
    if (as.logical(makeholdouts)) {
      ## Load the full input data
      df <- load_input_data(
        indicator = indicator,
        simple = NULL,
        removeyemen = FALSE,
        years = yearload,
        withtag = as.logical(withtag),
        datatag = datatag,
        use_share = as.logical(use_share)
      )
      
      ## Hard-coded sensitivity analyses
      if (sensitivity_analysis == 1) {
        df <- df[country != "ETH"]
      } else if (sensitivity_analysis == 2) {
        df <- df[country != "BFA"]
      } else if (sensitivity_analysis == 3) {
        df <- df[!(country == "IND" & year %in% c(2011:2013))]
      }
      
      # add in location information
      df <- merge_with_ihme_loc(df, re = region)
      
      # temporary hack: remove Botswana
      df <- df[country != "BWA"]
      
      message("n_ho_folds: ")
      message(as.integer(n_ho_folds))
      
      # make a list of dfs for each region, with 5 qt folds identified in each
      stratum_ho <- make_folds(
        data = df,
        n_folds = as.integer(n_ho_folds),
        spat_strat = spat_strat,
        temp_strat = temp_strat,
        strat_cols = "region",
        ts = as.numeric(ho_ts),
        mb = as.numeric(ho_mb),
        yr_col = yr_col[1],
        seed = random_seed
      ) # extra arg passed into quadtree_folds() to save shapefile (should be passed in through ... argument)
    } else {
      stratum_ho <- NULL
    }
    
    ###############################################################################
    ## Launch Parallel Script
    ###############################################################################
    
    ## Make loopvars aka strata grid (format = regions, ages, holdouts)
    if (as.logical(makeholdouts)) loopvars <- expand.grid(region, 0, 1:n_ho_folds) else loopvars <- expand.grid(region, 0, 0)
    
    ## loop over them, save images and submit qsubs
    if (!exists("restart_folds") | (is.null(restart_folds) == T) | identical(restart_folds, "NULL")) {
      loop_ints <- 1:nrow(loopvars)
    } else {
      if (class(restart_folds) == "character") restart_folds <- eval(parse(text = restart_folds))
      restart_folds <- as.integer(restart_folds)
      if (max(restart_folds) > nrow(loopvars)) {
        warning("The specified restart_folds exceeds the number of possible folds.")
      } else {
        loop_ints <- restart_folds
      }
    }
    
    message("loop_ints")
    message(loop_ints)
    
    ## Save image of workspace to facilitate debugging and interactive sessions
    for (i in loop_ints) {
      save.image(paste0(<<<< FILEPATH REDACTED >>>>))
    }
    
    ## Launch parallel MBG job
    parallel_job <- parallelize(slots            = fthread_parallel,
                                memory           = m_mem_free_parallel,
                                script           = parallel_code,
                                geo_nodes        = as.logical(use_geos_nodes),
                                expand_vars      = list(region = region, holdout = loopvars[[3]]),
                                save_objs        = c("indicator", "indicator_group", "run_date", "config_file", "df", "gaul_list", "simple_raster", "subset_shape", "simple_polygon", "pop_raster", "stratum_ho", "core_repo"),
                                prefix           = "parallel",
                                log_location     = 'sharedir',
                                script_dir       = paste0(<<<< FILEPATH REDACTED >>>>),
                                run_time         = runtime_parallel,
                                threads          = fthread_parallel,
                                singularity_opts = list(SET_OMP_THREADS = fthread_parallel, SET_MKL_THREADS = fthread_parallel))
    
    monitor_jobs(parallel_job)
  }
  
  ###############################################################################
  ## Post-Estimation
  ###############################################################################
  
  if (!parallel_pred_agg) { # Produce summary rasters if predictions weren't run in parallel
    ## Save strata for Shiny to use in producing aggregated fit statistics
    strata <- as.character(loopvars[, 1])
    stratas <- paste0(as.character(loopvars[, 1]), "_", as.character(loopvars[, 3]))
    dir.create(paste0(<<<< FILEPATH REDACTED >>>>), showWarnings = FALSE)
    save(strata, file = paste0(<<<< FILEPATH REDACTED >>>>))
    
    ## Post estimation to be done by strata
    counter <- 0
    for (s in stratas) {
      counter <- counter + 1
      message(s)
      
      ## Load spatial templates
      load(paste0(<<<< FILEPATH REDACTED >>>>))
      
      ## Load session image
      if (skip_to_summary_rasters) {
        pathaddin <- paste0("_bin", loopvars[counter, 2], "_", as.character(loopvars[counter, 1]), "_", loopvars[counter, 3])
        load(paste0(<<<< FILEPATH REDACTED >>>>))
      }
      
      ## Load cell draws
      load(paste0(<<<< FILEPATH REDACTED >>>>))
      
      yl_tag <- paste0("_", min(predict_years), "_", max(predict_years))
      
      ## Summarize predictions for each cell
      message("Starting mean raster summary")
      mean_raster <- insertRaster(simple_raster, matrix(rowMeans(cell_pred), ncol = length(predict_years)))
      
      message("Saving mean raster summary")
      save_post_est(mean_raster, "raster", paste0(s, "_unraked_mean_raster"))
      
      message("Starting upper 95 raster summary")
      upper_95_raster <- insertRaster(simple_raster, matrix(rowQuantiles(cell_pred, probs = 0.975), ncol = length(predict_years)))
      
      message("Saving upper 95 raster summary")
      save_post_est(upper_95_raster, "raster", paste0(s, "_unraked_upper_raster"))
      
      message("Starting lower 95 raster summary")
      lower_95_raster <- insertRaster(simple_raster, matrix(rowQuantiles(cell_pred, probs = 0.025), ncol = length(predict_years)))
      
      message("Saving lower 95 raster summary")
      save_post_est(lower_95_raster, "raster", paste0(s, "_unraked_lower_raster"))
      
      message("Starting range raster summary")
      range_raster <- upper_95_raster - lower_95_raster
      
      message("Saving range raster summary")
      save_post_est(range_raster, "raster", paste0(s, "_unraked_range_raster"))
      
      ## Some cleanup
      rm(mean_raster)
      rm(cell_pred)
      rm(upper_95_raster)
      rm(lower_95_raster)
      rm(range_raster)
    }
  }
  
  if (length(region) > 1 & indicator_family == "oncho") { # Combine regional summary rasters
    for(summ in c("lower", "mean", "upper", "range")){
      ras_l <- grep(paste0(summ, "_raster.tif"), list.files(paste0(<<<< FILEPATH REDACTED >>>>)), value = T)
      ras_l <- lapply(ras_l, FUN = function(m) brick(paste0(<<<< FILEPATH REDACTED >>>>)))
      combined_ras <- do.call(ras_l, FUN = function(m) brick(paste0(<<<< FILEPATH REDACTED >>>>)))
      save_post_est(combined_ras, "raster", "oncho_endem_afr_unraked_", summ, "_raster")
    }
  }
  
  ## Now run post-estimation script
  source(paste0(<<<< FILEPATH REDACTED >>>>))
}
