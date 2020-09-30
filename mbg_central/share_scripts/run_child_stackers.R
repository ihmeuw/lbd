
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

#create directory to store any stacker outputs
file.path(outputdir, "stackers") %>% dir.create

#Fit a gam model
if ('gam' %in% child_model_names) {
  tic("Stacking - GAM")
  gam <- fit_gam_child_model(df                = the_data,
                             model_name        = 'gam',
                             fold_id_col       = 'fold_id',
                             covariates        = all_fixed_effects,
                             additional_terms  = NULL,
                             weight_column     = 'weight',
                             bam               = TRUE,
                             spline_args       = list(bs = gam_basis, k = gam_knots %>% as.integer),
                             auto_model_select = TRUE,
                             indicator         = indicator,
                             indicator_family  = indicator_family,
                             cores             = 'auto')
  toc(log = T)
  
  #save gam obj and visualizations
  saveRDS(gam$gam, paste0(outputdir, "stackers/", reg, '_gam_obj.RDS'))
  b <- getViz(gam$gam)
  
  pdf(file=paste0(outputdir, "stackers/", reg, '_gam_plot.pdf'))
    print(plot(b, allTerms = T), pages=1)
  dev.off()
         
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
if ('xgboost' %in% child_model_names) {
  tic("Stacking - Xgboost")

  if (!exists("xg_model_tune")){
    message("xg_model_tune not found in config. Will be set to true, and model will be tuned with default grid search")
    xg_model_tune = TRUE
    hyperparameter_filepath = NULL
  }

  if (as.logical(xg_model_tune) & as.logical(xg_grid_search)){
    message("Building sobol sequence grid to tune xgb hyperparameters")
    hyperparameter_filepath = NULL

    gen_sequence <- function(min, max, rounding=4, scalar) {
      
      round(min + (max - min) * scalar, rounding) %>% return
    }

    #use the sobol sequence to generate a quasi random grid
    #note, scrambling=3 means:  Owen+Faure-Tezuka type of scrambling is applied
    sobol_scalars <- data.frame(runif.sobol(n = 200, dim = 6, scrambling = 3, seed = 98118, init = TRUE))
    names(sobol_scalars) <- c("nrounds", 'max_depth', "eta", "colsample_bytree", "min_child_weight", 'subsample')
    #note that carat returns a bug if this is a dt instead of df
    xg_grid <- data.frame(nrounds = gen_sequence(min=100, max=500, rounding=0, scalar=sobol_scalars$nrounds), 
                          max_depth = gen_sequence(min=2, max=3, rounding=0, scalar=sobol_scalars$max_depth),
                          eta = gen_sequence(min=.02, max=.2, scalar=sobol_scalars$eta),
                          colsample_bytree = gen_sequence(min=.4, max=1, scalar=sobol_scalars$colsample_bytree),
                          min_child_weight = gen_sequence(min=1, max=5, scalar=sobol_scalars$min_child_weight),
                          subsample = gen_sequence(min=.1, max=1, scalar=sobol_scalars$subsample),
                          gamma = 0)
                               
    summary(xg_grid)
    
  } else {
    message("Tuning xgboost with adaptive random search")
    xg_grid <- NA
    hyperparameter_filepath = NULL
    
  }
  
  #fit xgboost model
  xgboost <- fit_xgboost_child_model(df = the_data,
                                     covariates = all_fixed_effects,
                                     weight_column = 'weight',
                                     indicator = indicator,
                                     indicator_family = indicator_family,
                                     outputdir = outputdir,
                                     region = reg,
                                     xg_model_tune = xg_model_tune,
                                     xg_grid = xg_grid,
                                     build_ensemble = as.logical(xg_ensemble),
                                     return_second_best = as.logical(xg_second_best),
                                     hyperparameter_filepath = hyperparameter_filepath,
                                     debug=F)
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
                                  alpha            = 0,
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
                                  alpha            = 1,
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
  
  #plot coefficients from cross validation
  library(coefplot, lib.loc = package_lib)
  coefs <- coefplot(enet[['enet']], lambda=enet[['enet']]$cv_1se_lambda, sort='magnitude') + theme_minimal()
  ggsave(filename=paste0(outputdir, "stackers/", reg, '_glmnet_coeffplot.png'),
         plot=coefs)
  toc(log = T)
}
