### Post-estimation script
### Purpose: Sourced in Launch script to process model outputs for diagnostic purposes
##########################################################################################

# NOTE: there are many other diagnostic functions that could be incoporated. Check central post estimation script for reference.

##############################################################################
## Summarize model results
##############################################################################

### Load spatial templates
load(paste0(<<<< FILEPATH REDACTED >>>>))

# Duplicate child_model_list if skiptoinla from model run where child models were run
if (!(as.logical(makeholdouts))) {
  for (r in unique(as.character(loopvars[, 1]))) {
    for (h in unique(as.integer(loopvars[, 3]))) {
      if (as.logical(use_stacking_covs)) {
        if (as.logical(skiptoinla)) {
          file.copy(
            from = paste0(<<<< FILEPATH REDACTED >>>>),
            to = paste0(<<<< FILEPATH REDACTED >>>>)
          )
          
          file.copy(
            from = paste0(<<<< FILEPATH REDACTED >>>>),
            to = paste0(<<<< FILEPATH REDACTED >>>>)
          )
        }
        
        cov.wts <- get.cov.wts(
          rd = run_date, ## run_date
          ind = indicator, ## indicator
          ind_gp = indicator_group, ## indicator_group
          reg = r,
          age = 0,
          holdout = h
        )
        
        cov.plots <- plot.cov.wts(
          rd = run_date, ## run_date
          ind = indicator, ## indicator
          ind_gp = indicator_group, ## indicator_group
          reg = r,
          age = 0,
          holdout = h
        )
      }
    }
  }
}

###############################################################################
## Aggregate to admin2, admin1, and national levels
###############################################################################

if (!exists("skip_to_diagnostics")) {
  skip_to_diagnostics <- FALSE
} else {
  skip_to_diagnostics <- as.logical(skip_to_diagnostics)
}

if (!exists("skip_plot_stackers")) {
  skip_plot_stackers <- FALSE
} else {
  skip_plot_stackers <- as.logical(skip_plot_stackers)
}

pass_link_table <- modeling_shapefile_version

if (!skip_to_diagnostics) {
  if (as.logical(makeholdouts)) {
    holdouts <- 1:n_ho_folds
    
    if (as.logical(oos_preds_inla)) {
      outputdir <- sprintf(<<<< FILEPATH REDACTED >>>>) # new outputdir
      pathaddin <- paste0('_bin', 0, '_', region, '_0')
      if (exists("is_run_date")) {
        if (!is.null(is_run_date) & (is_run_date != "NULL")) {
          file.copy(
            from = paste0(<<<< FILEPATH REDACTED >>>>),
            to = paste0(<<<< FILEPATH REDACTED >>>>)
          )
          file.copy(
            from = paste0(<<<< FILEPATH REDACTED >>>>),
            to = paste0(<<<< FILEPATH REDACTED >>>>)
          )
        }
      }
      
      run_in_oos <- get_is_oos_draws_inla(
        ind_gp = indicator_group,
        ind = indicator,
        rd = run_date,
        age = 0,
        yrs = predict_years,
        get.oos = TRUE,
        write.to.file = TRUE,
        year_col = "year",
        shapefile_version = modeling_shapefile_version,
        pass_link_table = pass_link_table,
        holdouts = holdouts,
        region = region
      )
      
    } else {
      if (exists("is_run_date")) {
        if (!is.null(is_run_date) & (is_run_date != "NULL") & !as.logical(oos_preds_inla)) {
          pathaddin <- paste0(<<<< FILEPATH REDACTED >>>>)
          file.copy(
            from = paste0(<<<< FILEPATH REDACTED >>>>),
            to = paste0(<<<< FILEPATH REDACTED >>>>)
          )
        }
      }
      
      # Run oos metrics code if makeholdouts ################################################
      
      run_in_oos <- get_is_oos_draws(
        ind_gp = indicator_group,
        ind = indicator,
        rd = run_date,
        age = 0,
        nperiod = length(predict_years),
        yrs = predict_years,
        get.oos = TRUE,
        write.to.file = TRUE,
        year_col = "year",
        shapefile_version = modeling_shapefile_version,
        pass_link_table = pass_link_table
      )
    }
    
    ## set out_dir
    out_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
    dir.create(out_dir, recursive = T, showWarnings = F)
    
    ## for admin0
    draws.df <- fread(sprintf(<<<< FILEPATH REDACTED >>>>))
    
    country.pvtable <- get_pv_table(
      d = draws.df,
      indicator_group = indicator_group,
      rd = run_date,
      indicator = indicator,
      aggregate_on = "country",
      draws = as.numeric(samples),
      out.dir = out_dir,
      plot_ncol = 6
    )
    
    write.csv(<<<< FILEPATH REDACTED >>>>
    )
    
    ad1.pvtable <- get_pv_table(
      d = draws.df,
      indicator_group = indicator_group,
      rd = run_date,
      indicator = indicator,
      aggregate_on = "ad1",
      draws = as.numeric(samples),
      out.dir = out_dir,
      plot_ncol = 6
    )
    
    write.csv(ad1.pvtable,
              file = sprintf(<<<< FILEPATH REDACTED >>>>
              )
    )
    
    ad2.pvtable <- get_pv_table(
      d = draws.df,
      indicator_group = indicator_group,
      rd = run_date,
      indicator = indicator,
      aggregate_on = "ad2",
      draws = as.numeric(samples),
      out.dir = out_dir,
      plot_ncol = 6
    )
    
    write.csv(ad2.pvtable,
              file = sprintf(<<<< FILEPATH REDACTED >>>>
              )
    )
    
    if (!exists("validation_plot_by_country")) validation_plot_by_country <- TRUE
    if (is.null(validation_plot_by_country)) validation_plot_by_country <- TRUE
    
    if (validation_plot_by_country) {
      list_countries <- unique(draws.df$country)
      c_list <- lapply(list_countries, FUN = function(x) return(draws.df[country == x, ]))
      pvtable_list <- lapply(c_list, FUN = function(i) {
        if (nrow(i[!is.na(i$cell_pred_id)]) == 0) {
          message(paste0("No valid pixels for ", unique(i$country)))
        } else {
          if (nrow(i) == 1) {
            return(NULL)
          } else {
            c <- unique(i$country)
            out <- get_pv_table_lf(
              d = i,
              indicator_group = indicator_group,
              rd = run_date,
              indicator = indicator,
              aggregate_on = "Master_UID",
              draws = as.numeric(samples),
              out.dir = out_dir,
              page_per_year = FALSE,
              plot_title = paste0(c, ": OOS for ", indicator),
              plot = T,
              save_csv = F,
              plot_ci = T,
              plot_ncol = 6
            )
            return(out)
          }
        }
      })
      
      oos_plots <- lapply(pvtable_list, FUN = function(x) return(x[[2]]))
      names(oos_plots) <- list_countries
      save(oos_plots, file = paste0(<<<< FILEPATH REDACTED >>>>))
    }
    
    point.pvtable <- get_pv_table_lf(
      d = draws.df,
      indicator_group = indicator_group,
      rd = run_date,
      indicator = indicator,
      aggregate_on = "Master_UID",
      draws = as.numeric(samples),
      out.dir = out_dir,
      page_per_year = FALSE, plot_ci = T,
      plot_ncol = 6,
      result_agg_over = c("year", "oos")
    )
    
    point.pvtable.all <- get_pv_table_lf(
      d = draws.df,
      indicator_group = indicator_group,
      rd = run_date,
      indicator = indicator,
      aggregate_on = "Master_UID",
      draws = as.numeric(samples),
      out.dir = out_dir,
      page_per_year = FALSE, plot_ci = T,
      plot_ncol = 6,
      result_agg_over = c("oos"),
      plot = FALSE
    )
    
    point.pvtable.all <- cbind("Year" = paste0(min(point.pvtable[[1]]$Master_UID$Year), "-", max(point.pvtable[[1]]$Master_UID$Year)), point.pvtable.all[[1]])
    point.pvtable <- rbind(point.pvtable[[1]]$Master_UID, point.pvtable.all)
    
    write.csv(point.pvtable,
              file = sprintf(<<<< FILEPATH REDACTED >>>>
              )
    )
  } else {
    holdouts <- 0
    
    if (!exists("skip_submit_aggregation_script")) {
      skip_submit_aggregation_script <- FALSE
    }
    
    if (!as.logical(skip_submit_aggregation_script)) {
      overwrite <- FALSE
      raked <- FALSE
      age <- 0
      
      if (!parallel_pred_agg) {
        aggregation_job <- parallelize(slots            = fthread_agg,
                                       memory           = m_mem_free_agg,
                                       script           = "aggregate_results_lf.R",
                                       geo_nodes        = as.logical(use_geos_nodes),
                                       expand_vars      = list(region = Regions, holdout = holdouts),
                                       save_objs        = c("indicator", "indicator_group", "run_date", "config_file", "overwrite", "raked", "age", "core_repo"),
                                       prefix           = "aggregate",
                                       log_location     = 'sharedir',
                                       script_dir       = paste0(<<<< FILEPATH REDACTED >>>>),
                                       run_time         = runtime_agg,
                                       threads          = fthread_agg,
                                       singularity_opts = list(SET_OMP_THREADS = fthread_agg, SET_MKL_THREADS = fthread_agg))
        
        monitor_jobs(aggregation_job)
      } else {
        if (run_pop_threshold_analysis) {
          threshold_years <- c(2000, 2005, 2010, 2018)
          
          if (!is.null(pop_threshold_levels)) {
            threshold_indices <- pop_threshold_levels
            aggregation_job <- parallelize(slots            = fthread_agg,
                                           memory           = m_mem_free_agg,
                                           script           = "aggregate_results_lf_sensitivity_predict_year_parallel.R",
                                           geo_nodes        = as.logical(use_geos_nodes),
                                           expand_vars      = list(region = Regions, holdout = holdouts, year_predict = threshold_years, threshold_index = threshold_indices),
                                           save_objs        = c("indicator", "indicator_group", "run_date", "config_file", "overwrite", "raked", "age", "core_repo"),
                                           prefix           = "aggregate",
                                           log_location     = 'sharedir',
                                           script_dir       = paste0(<<<< FILEPATH REDACTED >>>>),
                                           run_time         = runtime_agg,
                                           threads          = fthread_agg,
                                           singularity_opts = list(SET_OMP_THREADS = fthread_agg, SET_MKL_THREADS = fthread_agg))
            
            monitor_jobs(aggregation_job)
          }
          
          thresholds <- list(c(0, 10000000, "_prev_0_pop_inf"), c(0, 750000, "_prev_0_pop_750k"), c(0, 250000, "_prev_0_pop_250k"), c(0, 100000, "_prev_0_pop_100k"), c(0, 50000, "_prev_0_pop_50k"), c(0, 25000, "_prev_0_pop_25k"), c(0, 10000, "_prev_0_pop_10k"), c(0, 5000, "_prev_0_pop_5k"))
          
          for (j in 1:8) {
            pop_threshold <- as.integer(thresholds[[j]][[1]])
            prev_threshold <- as.numeric(thresholds[[j]][[2]])
            file_suffix <- thresholds[[j]][[3]]
            
            admin_0_list <- vector("list", length(threshold_years))
            admin_1_list <- vector("list", length(threshold_years))
            admin_2_list <- vector("list", length(threshold_years))
            
            outputdir <- paste0(<<<< FILEPATH REDACTED >>>>)
            
            for (i in 1:length(threshold_years)) {
              load(paste0(<<<< FILEPATH REDACTED >>>>))
              admin_0_list[[i]] <- copy(admin_0)
              admin_1_list[[i]] <- copy(admin_1)
              admin_2_list[[i]] <- copy(admin_2)
              rm(admin_0); rm(admin_1); rm(admin_2);
            }
            
            admin_0 <- rbindlist(admin_0_list)
            admin_1 <- rbindlist(admin_1_list)
            admin_2 <- rbindlist(admin_2_list)
            
            save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = paste0(<<<< FILEPATH REDACTED >>>>))
          }
          
          quit()
        } else {
          threshold_index <- 1
          if (!skip_aggregation_parallelize_call) {
            aggregation_job <- parallelize(slots            = fthread_agg,
                                           memory           = m_mem_free_agg,
                                           script           = "aggregate_results_lf_sensitivity_predict_year_parallel.R",
                                           geo_nodes        = as.logical(use_geos_nodes),
                                           expand_vars      = list(region = Regions, holdout = holdouts, year_predict = predict_years),
                                           save_objs        = c("indicator", "indicator_group", "run_date", "config_file", "overwrite", "raked", "age", "threshold_index", "core_repo"),
                                           prefix           = "aggregate",
                                           log_location     = 'sharedir',
                                           script_dir       = paste0(<<<< FILEPATH REDACTED >>>>),
                                           run_time         = runtime_agg,
                                           threads          = fthread_agg,
                                           singularity_opts = list(SET_OMP_THREADS = fthread_agg, SET_MKL_THREADS = fthread_agg))
            
            monitor_jobs(aggregation_job)
          }
          
          admin_0_list <- vector("list", length(predict_years))
          admin_1_list <- vector("list", length(predict_years))
          admin_2_list <- vector("list", length(predict_years))
          
          outputdir <- paste0(<<<< FILEPATH REDACTED >>>>)
          
          for (i in 1:length(predict_years)) {
            load(paste0(<<<< FILEPATH REDACTED >>>>))
            admin_0_list[[i]] <- copy(admin_0)
            admin_1_list[[i]] <- copy(admin_1)
            admin_2_list[[i]] <- copy(admin_2)
            rm(admin_0); rm(admin_1); rm(admin_2);
          }
          
          admin_0 <- rbindlist(admin_0_list)
          admin_1 <- rbindlist(admin_1_list)
          admin_2 <- rbindlist(admin_2_list)
          
          save(admin_0, admin_1, admin_2, sp_hierarchy_list, file = paste0(<<<< FILEPATH REDACTED >>>>))
        }
      }
    }
    
    combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group, ages = 0, regions = region_list, holdouts = holdouts, raked = FALSE, delete_region_files = FALSE)
    
    summarize_admins(summstats = c("mean", "upper", "lower", "cirange"), ad_levels = c(0, 1, 2), raked = FALSE)
    
    ## Set run_date, indicator, indicator_group, out_dir per your preferences
    run_dates <- run_date
    run_label <- c(as.character(run_date))
    share_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
    in_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
    out_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
    dir.create(out_dir)
    
    in_file_ad0 <- paste0(<<<< FILEPATH REDACTED >>>>)
    in_file_ad1 <- paste0(<<<< FILEPATH REDACTED >>>>)
    in_file_ad2 <- paste0(<<<< FILEPATH REDACTED >>>>)
    
    # Prepare inputs ####################
    # Making a common list of regions input for below functions
    reg <- copy(region_list)
    if (reg[1] == "lf_endem_afr") {
      reg_l <- c("lf_endem_afr")
      reg_titles <- c(
        "LF Endemic Sub-Saharan Africa"
      )
    } else if (reg[1] == "IND") {
      reg_l <- c("IND")
      reg_titles <- c("India")
    } else if (reg[1] == "lf_s_se_asia") {
      reg_l <- c("lf_s_se_asia")
      reg_titles <- c("South/SE Asia")
    } else if (reg[1] == "lf_s_asia") {
      reg_l <- c("lf_s_asia")
      reg_titles <- c("South Asia")
    } else if (reg[1] == "lf_se_asia_pacific") {
      reg_l <- c("lf_se_asia_pacific")
      reg_titles <- c("SE Asia and Pacific")
    } else if (reg[1] == "lf_se_asia") {
      reg_l <- c("lf_se_asia")
      reg_titles <- c("Southeast Asia")
    } else if (reg[1] == "lf_hispaniola") {
      reg_l <- c("lf_hispaniola")
      reg_titles <- c("Hispaniola")
    } else if (reg[1] == "lf_vanuatu") {
      reg_l <- c("lf_vanuatu")
      reg_titles <- c("Vanuatu")
    } else if (reg[1] == "lf_wsm") {
      reg_l <- c("lf_wsm")
      reg_titles <- c("Samoa")
    } else if (reg[1] == "lf_bra") {
      reg_l <- c("lf_bra")
      reg_titles <- c("Brazil")
    } else if (reg[1] == "oncho_endem_afr") {
      reg_l <- c("oncho_endem_afr")
      reg_titles <- c("Oncho Endemic Africa")
    } else {
      message("Please specify reg_l and reg_titles yourself or add to the common list")
    }
    ###################################
    
    # Read in all mode admin aggregations, adding run label to each with add_run_label function
    ad0_df <- lapply(in_file_ad0, fread) %>% add_run_label(run_label)
    ad1_df <- lapply(in_file_ad1, fread) %>% add_run_label(run_label)
    ad2_df <- lapply(in_file_ad2, fread) %>% add_run_label(run_label)
    
    # Drop Ma'tan al-Sarra if present
    ad0_df <- subset(ad0_df, ADM0_CODE != 40762)
    ad1_df <- subset(ad1_df, ADM0_CODE != 40762)
    ad2_df <- subset(ad2_df, ADM0_CODE != 40762)
    
    admin_data <- input_aggregate_admin(indicator = indicator, indicator_group = indicator_group, regions = reg, indicator_family = "binomial", sample_column = "weighted_n", run_date = run_dates[length(run_dates)], svy_id = "Master_UID", shapefile_version = modeling_shapefile_version)
    
    ad0_data <- admin_data$ad0
    ad1_data <- admin_data$ad1
    ad2_data <- admin_data$ad2
    
    ad2_data[diagnostic == "ict", diagnostic := "ICT"]
    
    if (indicator_group == "lf") {
      mda_adm2 <- summarize_mda_admin2(cov = "lfmda", update = T, reg = reg, shapefile_version = modeling_shapefile_version, yl = predict_years)
    }
    
    # Run the plotting code ################################################
    
    if (indicator_group == "lf") {
      ind_title <- "LF Prevalence"
    } else {
      ind_title <- "Prevalence"
    }
    
    subnational_ts_plots(
      ad0_df = ad0_df,
      ad1_df = ad1_df,
      ad2_df = ad2_df,
      ad0_data = ad0_data,
      ad1_data = ad1_data,
      ad2_data = ad2_data,
      ind_title = ind_title,
      out_dir = out_dir,
      highisbad = T,
      val_range = c(0, 0.5),
      ad0_map_regions = reg_l,
      ad0_map_region_titles = reg_titles,
      plot_levels = c("ad0", "ad1", "ad2"),
      multiple_runs = F,
      verbose = T,
      plot_data = T,
      mda_data = mda_adm2,
      plot_stacker_trends = as.logical(use_stacking_covs),
      modeling_shapefile_version = modeling_shapefile_version
    )
    
  }
} else {
  message("Skipping to diagnostics")
}

##########################################################################
## Diagnostics ###########################################################
##########################################################################

reg <- copy(region_list)

## Folders & drive locations
commondir <- <<<< FILEPATH REDACTED >>>>
sharedir <- paste0(<<<< FILEPATH REDACTED >>>>)
out_dir <- paste0(<<<< FILEPATH REDACTED >>>>)

# Make sure that out_dir exists & create if not
dir.create(out_dir, recursive = T, showWarnings = F)

##########################################################################
## Plot model parameters #################################################

# Load model results
load(paste0(<<<< FILEPATH REDACTED >>>>))
if (!file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
  clean_model_results_table() # create a model results table
}

model_results <- read.csv(paste0(<<<< FILEPATH REDACTED >>>>), stringsAsFactors = F) %>% as.data.table()

names(model_results) <- c("region", "age", "parameter", "q0.025", "q0.5", "q0.975")

# Grab list of stackers
stackers <- model_results$parameter %>% unique()

stackers <- stackers[!(stackers %in% c(
  "int", "Nominal Range", "Nominal Variance",
  "GPRandom rho for time", "Precision for IID.ID",
  "Precision for CTRY.ID"
))]

# Plot stacker betas -----------------------------------------------------
#   Output of this part: plot of stacker betas & 95% CIs

gg_betas <- plot_stacker_betas(model_results = model_results, stackers = stackers)

png(filename = paste0(<<<< FILEPATH REDACTED >>>>),
    width = 10,
    height = 5,
    units = "in",
    type = "cairo-png",
    res = 200,
    pointsize = 10
)

print(gg_betas)
dev.off()

# Grab list of other parameters
other_params <- unique(model_results$parameter[!(model_results$parameter %in% c(stackers))])

# Plot other parameters --------------------------------------------------
#   Output of this part: plot of rho, int, variance, etc. from INLA

gg_other_params <- plot_other_params(model_results = model_results, other_params = other_params)

png(filename = paste0(<<<< FILEPATH REDACTED >>>>),
    width = 10,
    height = 5,
    units = "in",
    type = "cairo-png",
    res = 200,
    pointsize = 10
)

print(gg_other_params)
dev.off()

##########################################################################
## Plot maps of stackers #################################################

# Plot stackers vs true data for each region -----------------------------
#   Output of this part: maps of each stacker along with the mean raster
#   (unraked) and a plot of the data (with weight = alpha for resampled
#   polygons --> pseudoclusters), one for each year, for each region &
#   country.

if (!skip_plot_stackers) {
  plot_stackers(reg = reg, highisbad = F, individual_countries = T, shapefile_version = modeling_shapefile_version, yl = year_list, predict_years = predict_years)
}

##########################################################################
## Plot diagnostics for stackers

# Get a list of child models for each region

message("Loading child models for each region")
child_model_list <- lapply(Regions, function(reg) {
  child_model_file <- paste0(<<<< FILEPATH REDACTED >>>>)
  load(child_model_file, verbose = F)
  return(child_models)
})

names(child_model_list) <- Regions

# This list is now structured such that child_model_list[["cssa"]][["gam"]]
# will return the appropriate gam child object for the cssa region

# Plot GAM models ------------------------------------------------------
#   Output of this part: Plots of the component smooth functions that make
#   up the fitted gam model from stacking for each region's model.
#   Uses mgcv::plot.gam(). Limited to 12 plots per page, and will produce
#   multiple pages if needed to ensure that all covariates are graphed.

if ("gam" %in% stackers) {
  plot_gam_models(
    child_model_list = child_model_list,
    regions = Regions,
    o_dir = out_dir
  )
}

###############################################################################
# Create AROC objects & do projections
###############################################################################

make_aroc(
  ind_gp = indicator_group,
  ind = indicator,
  rd = run_date,
  matrix_pred_name = NULL,
  type = c("cell", "admin"),
  measure = "prevalence",
  year_list = predict_years,
  uselogit = FALSE,
  raked = FALSE,
  weighting_res = "domain",
  weighting_type = "exponential",
  pow = 1, ### SET TO ZERO TO USE FLAT WEIGHTING (pow = 1 weights more recent years more heavily)
  input_data = read.csv(sprintf(<<<< FILEPATH REDACTED >>>>)),
  mult_emp_exp = FALSE,
  extra_file_tag = "_exp_domain",
  shapefile_version = modeling_shapefile_version,
  pass_link_table = pass_link_table,
  parallel_pred_agg = parallel_pred_agg,
  regions = region
)

make_proj(
  ind_gp = indicator_group,
  ind = indicator,
  rd = run_date,
  type = c("cell", "admin"),
  proj_years = max(predict_years),
  measure = "prevalence",
  skip_cols = NULL,
  year_list = predict_years,
  uselogit = FALSE,
  raked = FALSE,
  extra_file_tag = "_exp_domain",
  shapefile_version = modeling_shapefile_version,
  parallel_pred_agg = parallel_pred_agg,
  regions = region
)

###############################################################################
# Look at performance against goals
###############################################################################
# Define goals: start by initializing goal object
goals <- add_goal(
  target_year = max(predict_years),
  target = 0.02,
  target_type = "less",
  abs_rel = "absolute",
  pred_type = c("cell", "admin")
)

# Add goals to existing goal object by specifying goal_obj
goals <- add_goal(
  goal_obj = goals,
  target_year = max(predict_years),
  target = 0.01,
  target_type = "less",
  abs_rel = "absolute",
  pred_type = c("cell", "admin")
)

# Run comparisons
compare_to_target(
  ind_gp = indicator_group,
  ind = indicator,
  rd = run_date,
  goal_obj = goals,
  measure = "prevalence",
  year_list = predict_years,
  shapefile_version = modeling_shapefile_version,
  uselogit = FALSE,
  raked = raking,
  pass_link_table = pass_link_table,
  regions = region
)
