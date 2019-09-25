
if ('SINGULARITY_NAME' %in% names(Sys.getenv())) {
  stacking_cores_to_use <- 1 # MKL does not play nicely with mclapply using multiple cores, so if
                             # using the singularity image (which includes MKL), provide only one
                             # core to mclapply during stacking
} else {
  stacking_cores_to_use <- cores_to_use
}

#Fit a gam model
if ('gam' %in% child_model_names) {
  tic("Stacking - GAM")
  gam <- fit_gam_child_model(df                = the_data,
                             model_name        = 'gam',
                             fold_id_col       = 'fold_id',
                             covariates        = all_fixed_effects,
                             additional_terms  = NULL,
                             weight_column     = 'weight',
                             bam               = FALSE,
                             spline_args       = list(bs = 'ts', k = 3),
                             auto_model_select = TRUE,
                             indicator         = indicator,
                             indicator_family  = indicator_family,
                             cores             = stacking_cores_to_use)
  toc(log = T)
}

#Fit a BRT model with Xgboost (faster than GBM)
if ('xgboost' %in% child_model_names){
  tic("Stacking - Xgboost")
  
  if (!exists("xg_model_tune")){
    message("xg_model_tune not found in config. Will be set to true, and model will be tuned with default grid search")
    xg_model_tune = TRUE
    hyperparameter_filepath = NULL
  }
  
  if (xg_model_tune == T & !exists("hyperparameter_filepath")){
    message("Tuning xgboost on default grid search")
    hyperparameter_filepath = NULL
  }
  xgboost <- fit_xgboost_child_model(df = the_data,
                                     covariates = all_fixed_effects,
                                     weight_column = 'weight',
                                     indicator = indicator,
                                     indicator_family = indicator_family,
                                     outputdir = outputdir,
                                     region = reg,
                                     xg_model_tune = xg_model_tune,
                                     hyperparameter_filepath = hyperparameter_filepath)
  toc(log = T)
}

#fit some nets
#lasso
if ('lasso' %in% child_model_names) {
  tic("Stacking - lasso")
  lasso <- fit_glmnet_child_model(df               = the_data,
                                  model_name       = 'lasso',
                                  covariates       = all_fixed_effects,
                                  fold_id_col      = 'fold_id',
                                  additional_terms = NULL,
                                  indicator_family = indicator_family,
                                  indicator        = indicator,
                                  cores            = stacking_cores_to_use,
                                  alpha            = 0,
                                  weight_column    = 'weight')
  toc(log = T)
}