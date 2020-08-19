### Process custom config options

if (exists("random_seed")) {
  if(identical(random_seed, "NULL")) {
    random_seed <- NULL
  } else {
    random_seed <- as.integer(random_seed)
  }
} else {
  random_seed <- NULL
}

if (!exists("use_premade")) {
  use_premade <- FALSE
}

if (!exists("restart_folds")) {
  restart_folds <- NULL
} else if (identical(restart_folds, "NULL")) {
  restart_folds <- NULL
}

if (!exists("skip_to_summary_rasters")) {
  skip_to_summary_rasters <- FALSE
} else {
  skip_to_summary_rasters <- as.logical(skip_to_summary_rasters)
}

if (!exists("sensitivity_analysis")) {
  sensitivity_analysis <- 0
} else {
  sensitivity_analysis <- as.integer(sensitivity_analysis)
}

if (!exists("run_pop_threshold_analysis")) {
  run_pop_threshold_analysis <- FALSE
} else {
  run_pop_threshold_analysis <- as.logical(run_pop_threshold_analysis)
}

if (!exists("pop_threshold_levels")) {
  pop_threshold_levels <- NULL
} else {
  if (!is.null(pop_threshold_levels)) {
    pop_threshold_levels <- eval(parse(text = pop_threshold_levels))
  }
}

# Impute missing covariate values?
if (!exists("impute_missing_covariates")) {
  impute_missing_covariates <- FALSE
}

# Reload saved imputated data set?
if (!exists("use_saved_imputations")) {
  use_saved_imputations <- FALSE
}

if (!exists("predict_years")) {
  predict_years <- copy(year_list)
} else if (is.null(predict_years) | identical(predict_years, "NULL")) {
  predict_years <- copy(year_list)
}

if (!exists("perform_initial_INLA_approx")) {
  perform_initial_INLA_approx <- FALSE
} else {
  perform_initial_INLA_approx <- as.logical(perform_initial_INLA_approx)
}

if (!exists("parallel_pred_agg")) {
  parallel_pred_agg <- FALSE
} else {
  parallel_pred_agg <- as.logical(parallel_pred_agg)
}

if (!exists("oos_preds_inla")) {
  oos_preds_inla <- FALSE
} else {
  oos_preds_inla <- as.logical(oos_preds_inla)
}

if (!exists("restart_predict_years")) {
  restart_predict_years <- NULL
} else {
  if (!is.null(restart_predict_years)) {
    restart_predict_years <- eval(parse(text = restart_predict_years))
  }
}

if (!exists("temporal_model_theta_prior")) {
  temporal_model_theta_prior <- "list(prior = 'loggamma', param = c(1, 0.00005))"
}

if (!exists("use_space_only_gp")) {
  use_space_only_gp <- FALSE
} else {
  use_space_only_gp <- as.logical(use_space_only_gp)
}

if (!exists("mesh_cutoff")) {
  mesh_cutoff <- eval(parse(text = mesh_s_max_edge))[1]
}

if (!exists("clamp_covariates")) {
  clamp_covariates <- FALSE
} else {
  clamp_covariates <- as.logical(clamp_covariates)
}

if (!exists("oos_preds_inla")) {
  oos_preds_inla <- FALSE
} else {
  oos_preds_inla <- as.logical(oos_preds_inla)
}

if (!exists("initial_INLA_approx_levels") | !perform_initial_INLA_approx) {
  initial_INLA_approx_levels <- 0
} else {
  initial_INLA_approx_levels <- as.integer(initial_INLA_approx_levels)
}

if (!exists("temporal_model_pc_prec_prior")) {
  temporal_model_pc_prec_prior <- "list(param = c(3, 0.01)"
}

if (!exists("temporal_model_pc_pacf1_prior")) {
  temporal_model_pc_pacf1_prior <- "list(param = c(0.5, 0.5)"
}

if (!exists("pc_priors_ar_cor1")) {
  pc_priors_ar_cor1 <- NULL
}

if (!exists("skip_aggregation_parallelize_call")) {
  skip_aggregation_parallelize_call <- FALSE
} else {
  skip_aggregation_parallelize_call <- as.logical(skip_aggregation_parallelize_call)
}

if (!as.logical(use_stacking_covs)) {
  use_age_diagnostics_in_stackers <- FALSE
} else {
  use_rw1_crosswalk_covariates <- FALSE
  fit_crosswalk_inla <- FALSE
}

if (exists("other_weight")) {
  if (other_weight == "" | is.null(other_weight)) {
    rm(other_weight)
  }
}
