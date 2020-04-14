
# Parallelization in stackers is a little bit complicated. Some of the stackers
# rely on packages like "mgcv" and "glmnet" which have their own internal
# parallelization which we have to make sure plays nice with the MKL.
# `mclapply()` is used heavily in stackers as well which will hang with
# multithreaded operations like OpenMP or MKL.

# Those stacking functions which have a "cores" argument have an "auto" option
# (default) where the stacker function has code to properly determine how many
# cores to give each operation based on how it is parallelized, how many cores
# are available for the job based on slots, etc. It is highly recommended that
# the "auto" option be used. Otherwise, you may pass in an integer argument for
# the number of cores. The only value that is known to work generally on all
# stackers is 1, or serial operation.

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
                             cores             = 'auto')
  toc(log = T)
}

#Fit a GBM/BRT model
if ('gbm' %in% child_model_names) {
  tic("Stacking - GBM")
  gbm <- fit_gbm_child_model(df               = the_data,
                             model_name       = 'gbm',
                             fold_id_col      = 'fold_id',
                             covariates       = all_fixed_effects_brt,
                             weight_column    = 'weight',
                             tc               = as.numeric(gbm_tc),
                             lr               = as.numeric(gbm_lr),
                             bf               = as.numeric(gbm_bf),
                             indicator        = indicator,
                             indicator_family = indicator_family,
                             cores            = 'auto')
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
  
  if (xg_model_tune == T & exists("hyperparameter_filepath")){
    message("Tuning xgboost on pre-specified grid search")
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
                                  alpha            = 1,
                                  weight_column    = 'weight',
                                  cores            = 'auto')
  toc(log = T)
}

#ridge
if ('ridge' %in% child_model_names) {
  tic("Stacking - ridge")
  ridge <- fit_glmnet_child_model(df               = the_data,
                                  model_name       = 'ridge',
                                  covariates       = all_fixed_effects,
                                  fold_id_col      = 'fold_id',
                                  additional_terms = NULL,
                                  indicator_family = indicator_family,
                                  indicator        = indicator,
                                  alpha            = 0,
                                  weight_column    = 'weight',
                                  cores            = 'auto')
  toc(log = T)
}

#enet
if ('enet' %in% child_model_names) {
  tic("Stacking - enet")
  enet = fit_glmnet_child_model(df               = the_data,
                                model_name       = 'enet',
                                covariates       = all_fixed_effects,
                                fold_id_col      = 'fold_id',
                                additional_terms = NULL,
                                indicator_family = indicator_family,
                                indicator        = indicator,
                                alpha            = 0.5,
                                weight_column    = 'weight',
                                cores            = 'auto')
  toc(log = T)
}
