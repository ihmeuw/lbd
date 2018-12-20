
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
                             cores             = stacking_cores)
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
                             cores            = stacking_cores)
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
                                  cores            = stacking_cores,
                                  alpha            = 0,
                                  weight_column    = 'weight')
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
                                  cores            = stacking_cores,
                                  alpha            = 1,
                                  weight_column    = 'weight')
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
                                cores            = stacking_cores,
                                alpha            = 0.5,
                                weight_column    = 'weight')
  toc(log = T)
}
