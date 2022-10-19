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

if (!exists("assumed_programmatic_stage")) {
  assumed_programmatic_stage <- NULL
}

if (!exists("mesh_cutoff")) {
  mesh_cutoff <- eval(parse(text = mesh_s_max_edge))[1]
}

if (!exists("use_pref_samp_pp")) {
  use_pref_samp_pp <- FALSE
} else {
  use_pref_samp_pp <- as.logical(use_pref_samp_pp)
}

if (!exists("pref_samp_pp_stages")) {
  pref_samp_pp_stages <- "all"
}

if (!exists("covars_in_pp")) {
  covars_in_pp <- FALSE
} else {
  covars_in_pp <- as.logical(covars_in_pp)
}

if (!exists("country_re_pp")) {
  country_re_pp <- FALSE
} else {
  country_re_pp <- as.logical(country_re_pp)
}

if (!exists("use_spatially_varying_coefficients")) {
  use_spatially_varying_coefficients <- FALSE
} else {
  use_spatially_varying_coefficients <- as.logical(use_spatially_varying_coefficients)
}

if (!exists("fit_crosswalk_inla")) {
  fit_crosswalk_inla <- FALSE
} else {
  fit_crosswalk_inla <- as.logical(fit_crosswalk_inla)
}

if (!exists("clamp_covariates")) {
  clamp_covariates <- FALSE
} else {
  clamp_covariates <- as.logical(clamp_covariates)
}

if (!exists("perform_initial_INLA_approx")) {
  perform_initial_INLA_approx <- FALSE
} else {
  perform_initial_INLA_approx <- as.logical(perform_initial_INLA_approx)
}

if (!exists("oos_preds_inla")) {
  oos_preds_inla <- FALSE
} else {
  oos_preds_inla <- as.logical(oos_preds_inla)
}

if (!exists("oos_preds_inla_nudge_to_centroid")) {
  oos_preds_inla_nudge_to_centroid <- FALSE
} else {
  oos_preds_inla_nudge_to_centroid <- as.logical(oos_preds_inla_nudge_to_centroid)
}

if (!exists("predict_diagnostic")) {
  if (indicator_group == "lf") {
    predict_diagnostic <- "ict"
  } else if (indicator_group == "oncho") {
    predict_diagnostic <- "ss"
  }
}

if (!exists("predict_age_start")) {
  predict_age_start <- 0
} else {
  predict_age_start <- as.integer(predict_age_start)
}

if (!exists("predict_age_end")) {
  predict_age_end <- 94
} else {
  predict_age_end <- as.integer(predict_age_end)
}

if (!exists("ppp_region")) {
  ppp_region <- NULL
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

if (!exists("crosswalk_training_data_set")) {
  crosswalk_training_data_set <- NULL
} else if (identical(crosswalk_training_data_set, "NULL")) {
  crosswalk_training_data_set <- NULL
}

if (!exists("use_river_size")) {
  use_river_size <- FALSE
} else {
  use_river_size <- as.logical(use_river_size)
}

if (!exists("use_rw1_crosswalk_covariates")) {
  use_rw1_crosswalk_covariates <- FALSE
} else {
  use_rw1_crosswalk_covariates <- as.logical(use_rw1_crosswalk_covariates)
}

if (!exists("use_age_diagnostics_in_stackers")) {
  use_age_diagnostics_in_stackers <- FALSE
} else {
  use_age_diagnostics_in_stackers <- as.logical(use_age_diagnostics_in_stackers)
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

if (!exists("use_oncho_suitability")) {
  use_oncho_suitability <- FALSE
}

if (!exists("rw1_raw_covar_list")) {
  rw1_raw_covar_list <- NULL
} else {
  rw1_raw_covar_list <- eval(parse(text = rw1_raw_covar_list))
}

if (!exists("pp_covars")) {
  pp_covars <- NULL
} else if (!is.null(pp_covars) & !identical(pp_covars, "NULL")) {
  pp_covars <- eval(parse(text = pp_covars))
}

if (!exists("clamping_quantiles")) {
  clamping_quantiles <- c(0, 1)
} else {
  clamping_quantiles <- eval(parse(text = clamping_quantiles))
}

if (!exists("keep_Master_UID_list")) {
  keep_Master_UID_list <- NULL
} else {
  keep_Master_UID_list <- eval(parse(text = keep_Master_UID_list))
}

if (subnat_country_to_get != "all") {
  subnat_country_to_get <- eval(parse(text = subnat_country_to_get))
}

if (!exists("spde_alpha")) {
  spde_alpha <- 2
} else {
  spde_alpha <- as.numeric(spde_alpha)
}

if (!exists("rw_covariate_model")) {
  rw_covariate_model <- "rw1"
}

if (!exists("rw_covariate_n_groups")) {
  rw_covariate_n_groups <- 25
} else {
  rw_covariate_n_groups <- as.integer(rw_covariate_n_groups)
}

if (!exists("control.inla.list.h")) {
  control.inla.list.h <- 0.01
} else {
  control.inla.list.h <- as.numeric(control.inla.list.h)
}
