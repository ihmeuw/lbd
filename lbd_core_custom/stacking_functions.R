
# fit_gbm ---------------------------------------------------------------------------------------------------------------------------
#df: a data table (post extract covariates)
#covariates: rhs of formula specifying the covariates/columns to be used to help fit the model
#weight_column: column in the data table that specifies the weight
#tc: tree complexity
#lr: learning rate
#bf: bag fraction
#cv: cross-validation folds
#indicator: dependant variable.
#indicator_family: Binomial models assume N as the # of trials and are actually modelled with poission with N as the offset

fit_gbm= function(df, covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, tc = gbm_tc, lr = gbm_lr, bf = gbm_bf, cv = gbm_cv, nminobs = gbm_nminobs, ntrees = gbm_ntrees,
                  indicator, indicator_family = 'binomial', plot.main = F){
  
  library(dismo)
  
  # check to see if its a vector of characters or a psudo-formula
  covariates = format_covariates(add_additional_terms(covariates,additional_terms))
  
  df = copy(df)
  
  # format weights
  if(!is.null(weight_column)){
    df[,data_weight := get(weight_column)]
  } else{
    df[,data_weight := 1]
  }
  weight_column = 'data_weight' #specify the column
  
  # BRT
  message(paste('Fitting GBM/BRT with tc:',tc,'lr:',lr,'bf:',bf,'cv:',cv))
  
  # create gbm formula
  df <- as.data.table(df)
  df[, response := round(get(indicator), 0)]
  df[, log_n  := log(N)]
  cov_names <- paste(covariates, collapse = ' + ')
  gbm_formula <- formula(paste0('response ~ offset(log_n) + ', cov_names))
  
  # run gbm model a la Josh Osborne
  mod <- gbm(gbm_formula,
             distribution      = 'poisson',
             data              = df, 
             n.trees           = ntrees, 
             interaction.depth = tc,
             shrinkage         = lr,
             n.minobsinnode    = nminobs,
             weights           = df[,get(weight_column)],
             bag.fraction      = bf,
             train.fraction    = 1-(1/cv),
             verbose           = F)
  
  if(is.null(mod)) stop('BRT ATTEMPT FAILED')
  
  return(mod)
}

# fit_gbm_child_model ---------------------------------------------------------------------------------------------------
#
#a function to fit the GBMs for stacking
#basically a wrapper function for the fit_gam function (which is its own wrapper function-- hooray for rabbit holes)
#model_name = what do you want the full fit model to be called upon the return. Must sync with subsquent functions

fit_gbm_child_model = function(df, model_name = 'gbm', fold_id_col = 'fold_id', covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL,
                               tc = gbm_tc, lr = gbm_lr, bf = gbm_bf, cv = gbm_cv, nminobs = gbm_nminobs, ntrees = gbm_ntrees, indicator = indicator, indicator_family = indicator_family, cores = 'auto'){
  
  #######The function#######
  library(parallel)
  
  #prevent df scoping
  df = copy(df)
  
  #fit the baby trees in parallel
  folds = unique(df[,get(fold_id_col)])
  
  message('Fitting baby gbm models in parallel')
  # Set multithreading to serial for `mclapply()`:
  set_serial_threads()
  # Determine appropriate number of cores to use in `mclapply()`
  if(cores == 'auto') cores <- get_max_forked_threads(nobjs = length(folds))
  baby_models = mclapply(folds, function(fff)
    fit_gbm(df = df[get(fold_id_col) != fff,],
            covariates = covariates,
            additional_terms   = additional_terms,
            weight_column      = weight_column,
            tc                 = tc,
            lr                 = lr,
            bf                 = bf,
            cv                 = cv,
            nminobs            = nminobs,
            ntrees             = ntrees,
            indicator          = indicator,
            indicator_family   = indicator_family,
            plot.main = F), mc.cores = cores)
  # Return to multithreading (if any):
  set_original_threads()
  
  for(fff in folds){
    #use the data fit on K-1 of the folds to fit on the help out fold
    df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_models[[fff]], df[get(fold_id_col)==fff,], n.trees=baby_models[[fff]]$n.trees, type='response')]
  }
  
  
  #fit GBM
  message('Fitting Full GBM')
  full_model = fit_gbm(df = df,
                       covariates = covariates,
                       additional_terms   = additional_terms,
                       weight_column      = weight_column,
                       tc                 = tc,
                       lr                 = lr,
                       bf                 = bf,
                       cv                 = cv,
                       nminobs            = nminobs,
                       ntrees             = ntrees,
                       indicator          = indicator,
                       indicator_family   = indicator_family)
  
  
  #add a model name slot
  full_model$model_name = model_name
  
  #predict the main BRT
  print(summary(full_model))
  print(str(full_model))
  print(full_model)
  df[,paste0(model_name,'_full_pred') := predict(full_model, df, n.trees = full_model$n.trees, type = 'response')]
  
  suffixes = c('_full_pred', '_cv_pred')
  return_cols = paste0(model_name, suffixes)
  #print(return_cols)
  #set up with for the return call
  return(setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name))))
  
}

#fit gam: a function to fit a gam or bam
#df: data table with the outcome/indicator and some covariates already extracted. This is different than the gam_cov functions
#covariates: a  vector of covariate names and/or formula in the style of all_fixed_effects (e.g. cov1 + cov2 + cov3)
#additional terms: a vector or single character of column names to be included in the model fit
#weight_column: in df, is there a column that specifies the observation weights?
#bam: should bam be used rather than gam?
#spline args: additional arguments to be sent to the spline call of gam/bam
#auto_model_select: if true, it overwrites BAM instructions to fit a GAM on smaller (N<2000) datasets-- helps with convergence
#indicator: name of the column of the DV
#indicator_family: model family
#cores: # of cores are available for use
fit_gam = function(df, covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, bam = F, spline_args = list(), auto_model_select =F, indicator, indicator_family = 'binomial', cores = 'auto'){
  library(mgcv)
  library(data.table)
  #also requires seeg
  
  df = copy(df) #in case data table scoping gets wonky
  
  #check to see if the gam formula is prespecified. if not, make it
  #check to see if its a vector of characters or a psudo-formula
  covariates = format_covariates(covariates) #additional terms is handled below
  
  
  #remove binary vars from the splined versions and set as additional terms
  n_uniq_vals = unlist(setNames(df[,lapply(covariates, function(x) uniqueN(get(x)))], covariates)) #this needs to return a named vector
  drop_vars = names(n_uniq_vals[n_uniq_vals<=as.numeric(gam_knots)+1])
  covariates = covariates[!covariates %in% drop_vars]
  
  additional_terms = c(additional_terms,drop_vars)
  additional_terms = additional_terms[!(is.null(additional_terms) | is.na(additional_terms))]
  
  #set response variable (stolen from GAM trans)
  if(indicator_family=="binomial") response <- cbind(events = df[, get(indicator)], trials = df[, N] - df[, get(indicator)])
  if(indicator_family=="gaussian") response <- cbind(outcome = df[, get(indicator)])
  
  
  #sort out the formula body
  f_bod = paste(paste0('s(',covariates,', ',parseArgsS(spline_args),')'),collapse = " + ")
  gam_formula = paste0("response ~ 1 + ", f_bod)
  if(length(additional_terms)>0) gam_formula = paste0(gam_formula,' + ', paste(additional_terms, collapse = " + "))
  
  gam_formula = as.formula(gam_formula) #final model formula
  
  #sort out weights
  #format weights
  if(!is.null(weight_column)){
    df[,data_weight := get(weight_column)]
  } else{
    df[,data_weight := 1]
  }
  weight_column = 'data_weight'
  
  
  # The `gam` and `bam` functions in the mgcv package have their own internal
  # OpenMP parallelization which can also use the MKL. So, we'll use those cores
  # set aside for OpenMP here for the `gam` and `bam` internal parallelization,
  # namely OMP_NUM_THREADS:
  if(cores == 'auto') cores <- get_omp_threads()
  
  #fit the gam
  message(paste0('Fitting GAM/BAM with spline args of: ', names(spline_args)[1],'=',spline_args[1],' ', names(spline_args)[2],'=', spline_args[2]))
  if(bam){
    if(auto_model_select ==T & nrow(df)<2000){
      model = mgcv::gam(gam_formula, data = df, family = indicator_family, weights = df[,get(weight_column)], control = list(nthreads = as.numeric(cores)))
    } else{
      model = mgcv::bam(gam_formula, data = df, family = indicator_family, weights = df[,get(weight_column)], nthreads = as.numeric(cores), discrete = T)
    }
  } else{
    model = mgcv::gam(gam_formula, data = df, family = indicator_family, weights = df[,get(weight_column)], control = list(nthreads = as.numeric(cores)))
  }
  
  #return the gam object
  return(model)
  
}

#a wrapper function for fitting gam/bam models in the stacking framework. It'll run 1+k times internally
#arguments are the same as fit gam above.
#basically a wrapper function for the fit_gam function (which is its own wrapper function-- hooray for rabbit holes)
fit_gam_child_model = function(df, model_name = 'gam', fold_id_col = 'fold_id', covariates = all_fixed_effects,
                               additional_terms = NULL, weight_column = NULL, bam =F, spline_args = list(bs = 'ts', k = as.numeric(gam_knots)),
                               auto_model_select =T, indicator = indicator, indicator_family = 'binomial', cores = 'auto',
                               outdir = outputdir, region = reg){
  library(mgcv)
  #remove scoping surprises
  df = copy(df)
  
  #start by fitting the full gam
  message('Fitting the Full GAM model')
  full_model = fit_gam(df,
                       covariates = covariates,
                       additional_terms = additional_terms,
                       weight_column = weight_column,
                       bam = bam,
                       spline_args = spline_args,
                       auto_model_select =auto_model_select,
                       indicator = indicator,
                       indicator_family = indicator_family,
                       cores = cores)
  
  #add a name to the game object
  full_model$model_name = model_name
  
  #fit the child/cv gams
  message("Fitting baby gams")
  #fit the child/cv rfs
  folds = unique(df[,get(fold_id_col)])
  
  for(fff in folds){
    baby_model = fit_gam(df = df[get(fold_id_col) != fff,],
                         covariates=covariates,
                         additional_terms = additional_terms,
                         bam = bam,
                         weight_column = weight_column,
                         spline_args = spline_args,
                         auto_model_select = auto_model_select,
                         indicator = indicator,
                         indicator_family = indicator_family,
                         cores = cores)
    
    #fill in the data
    df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_model, df[get(fold_id_col)==fff,],type = 'response')]
  }
  
  #predict using full model fit earlier
  df[,paste0(model_name,'_full_pred') := predict(full_model,df,type = 'response')]
  
  #return a subset of the columns. Full pred denotes the fit from the full model. CV pred is the OOS stuff
  suffixes = c('_full_pred', '_cv_pred')
  return_cols = paste0(model_name, suffixes)
  
  #set up with for the return call
  gam_fit <- setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name)))
  
  # Create stacking directory to save results
  stack_dir <- paste0(outdir, 'stackers/')
  dir.create(stack_dir, showWarnings = F)
  
  # plot fit results
  pdf(paste0(stack_dir, region, '_gam_stacker_predictors.pdf'))
  plot(gam_fit$gam, pages = 1)
  dev.off()
  
  # return
  return(gam_fit)
}

## fit xgboost child model function ###################################################

#' This function fits xgboost, another implementation of boosted regression trees. Xgboost has a different fitting algorithm
#' as compared to gbm, and tends to fit much faster with fewer number of iterations required. This speeds up the run time considerably.
#'
#' This function takes a data frame with the extracted covariates at each lat/long and fits regression trees on your indicator.
#' The added bonus of xgboost for prevalence indicators is it can perform logistic regression with an outcome variable between
#' [0,1], so we no longer need to convert brt's to a poisson distribution.
#'
#' This function also has built in cross-validation that uses 5 fold cross validation repeated 5 times for model accuracy.
#' You can either pass the hyperparameter space as a file path to a csv file, or go with the default grid search by not providing a filepath.
#' Xgboost has a large number of parameters, but as of now the only 3 tunable hyperparameters. The number of rounds (nround) specifies the number of iterations.
#' Eta (learning rate) which specifies the shrinkage used in updating each tree fit to prevent overfitting. After each boosting step, we can
#' directly get the weights of new features, and eta shrinks the feature weights to make the boosting process more conservative. Finally, we also vary
#' the max depth, which specifies the maximum depth of the tree. Increasing this value will make the model more complex and more likely to overfit.
#'
#' @param df Data Frame. Data frame that contains your indicator and the covariate values extracted at those cluster locations
#' @param indicator Character. Model indicator. Make sure in your data frame that this exists
#' @param indicator_family Character. Indicator statistical family. Common inputs include Binomial and Gaussian
#' @param outputdir File path. Model output directory, defined early on in the parallel script
#' @param region Character. Region being modeled over
#' @param covariates Character. List of covariates being used in stacking. For example: access2 + aridity + fertility
#' @param weight_column Numeric. This column contains the weights form polygon resampling.
#' @param xg_model_tune Logical. If true, xgboost will perform grid search to pick best hyperparameters
#' @param hyperparameter_filepath File path. If true, xgboost will perform a grid search to find the optimal hyperparameter.
#'                                If false, it will read from config what the optimal parameters. You must specify nrounds, eta and max depth
#'
#' @return a list containing: 1) The out of sample and in sample predictions from the xgboost fit, labeled "dataset"
#'                            2) The child model object, labeled "xgboost"



fit_xgboost_child_model <- function(df,
                                    indicator = indicator,
                                    indicator_family = "binomial",
                                    outputdir,
                                    region,
                                    covariates = all_fixed_effects,
                                    weight_column = 'weight',
                                    k_folds = 5, #use 5-fold NID crossvalidation by default
                                    xg_model_tune = TRUE,
                                    xg_grid = NA,
                                    xg_ensemble_corr = .5,
                                    build_ensemble = F,
                                    hyperparameter_filepath = NULL,
                                    debug=F){
  
  load_R_packages(c("xgboost", "caret", 'ggrepel', 'caretEnsemble'))
  
  if (debug) browser()
  
  # Create stacking directory to save results
  stack_dir <- paste0(outputdir, "stackers/")
  dir.create(stack_dir, showWarnings = F)
  
  # Create model formula
  df <- as.data.table(df)
  setnames(df, indicator, "indicator")
  form <- as.formula(paste0('indicator ~ ', covariates))
  
  # Create custom weight column for xgboost, weight * sample size
  df[, xg_weight := get(weight_column) * N]
  
  # Make sure to model in prevalence space if binomial
  if (indicator_family == "binomial"){
    df[, indicator := indicator / N]
    objective_function = "reg:logistic"
  }
  
  # If gaussian indicator make objective function linear
  if (indicator_family == "gaussian") objective_function = "reg:linear"
  
  if(xg_model_tune == F & is.null(hyperparameter_filepath)){
    stop("If you are not tuning xgboost you must provide a filepath to chosen hyperparameters./n
         Look at the hyperparameter_filepath argument to this function")
  }
  
  if (xg_model_tune == T) {
    message("Model tuning xgboost")
    
    # Set grid search as default unless filepath is provided
    if(is.null(hyperparameter_filepath)) message("Tuning with default hyperparameter settings")
    else {
      message("Selecting pre-specified hyperparameter grid")
      hyperparam <- fread(hyperparameter_filepath)
      
      # Define grid search from the pre-specified hyperparameter settings
      xg_grid <- expand.grid(nrounds = hyperparam$nrounds,
                             max_depth = eval(parse(text=hyperparam$max_depth)),
                             eta = eval(parse(text=hyperparam$eta)),
                             colsample_bytree = .5,
                             min_child_weight = 1,
                             subsample = 1,
                             gamma = 0)
      
    }
    
    # Specify the training folds using NID to account for dependence and prevent overfitting
    k_folds <- uniqueN(df$nid) %>% ifelse(.<k_folds, ., k_folds) #make sure that nid count > k
    folds <- groupKFold(df$nid, k = k_folds) 
    
    # Set cross validation options, default to 5 times repeated 5-fold cross validation
    # Selection function is "oneSE" to pick simplest model within one standard error of minimum
    
    if (xg_grid %>% is.na) {
      message('Tuning XGBoost using adaptive random search')
      # Try using an adaptive random search instead of a cartesian grid, should be more efficient
      train_control <- trainControl(selectionFunction = "oneSE",
                                    method = "adaptive_cv",
                                    search = "random",
                                    adaptive = list(min = 5, alpha = 0.05, 
                                                    method = "gls", complete = TRUE),
                                    number = k_folds,
                                    repeats = 5,
                                    index = folds,
                                    returnResamp = 'all',
                                    savePredictions = 'all')
      
      
      # Fit model
      xg_fit <- train(form,
                      data = df,
                      trControl = train_control,
                      verbose = F,
                      tuneLength = 25, #tune for 25 runs 
                      metric = "RMSE",
                      method = "xgbTree",
                      objective = objective_function,
                      weights = df$xg_weight)
    } else {
      message('Tuning XGBoost using grid search')
      #use provided grid to search
      train_control <- trainControl(selectionFunction = "oneSE",
                                    method = "repeatedcv",
                                    number = k_folds,
                                    repeats = 5,
                                    index = folds,
                                    returnResamp = 'all',
                                    savePredictions = 'all')
      
      # Fit model
      xg_fit <- train(form,
                      data = df,
                      trControl = train_control,
                      verbose = F,
                      tuneGrid = xg_grid,
                      metric = "RMSE",
                      method = "xgbTree",
                      objective = objective_function,
                      weights = df$xg_weight)
    }
    
    # find the 2nd best model (balance between lowest correlation to best and highest RMSE)
    # first pull out the relevant objects from the tuning object
    best_tune <- xg_fit$bestTune %>% as.data.table %>% setkey
    message('Found best tune, with parameters: \n')
    summary(best_tune)
    
    tune_cols <- key(best_tune) #capture all the tuning params as a col vector
    preds <- xg_fit$pred %>% as.data.table %>% setkeyv(., tune_cols)
    samps <- xg_fit$resample %>% as.data.table %>% setkeyv(., c(tune_cols, 'Resample'))
    
    #pull out the best preds
    best_preds <- preds[best_tune] %>% setkeyv(., key(samps)) # isolate the predictions from best model
    preds <- preds[!best_tune] %>% setkeyv(., key(samps))  # all other predictions
    
    #add the samp scores to the preds
    best_preds <- merge(best_preds, samps, by=key(samps))
    preds <- merge(preds, samps, by=key(samps)) 
    
    #merge on the preds and calculate the correlation to the best model
    preds <- merge(preds, best_preds[, .(rowIndex, Resample, best_pred=pred, best_rmse=RMSE)], by=c('rowIndex', 'Resample'))
    preds[, delta_rmse := RMSE-best_rmse]
    preds[, model_id := .GRP, by=tune_cols]
    setkey(preds, model_id, Resample) #key on the model ID and fold
    preds[, cor := corSpearman(pred, best_pred), by=key(preds)]
    preds[, pooled_cor := corSpearman(pred, best_pred), by=model_id]
    preds[, pooled_rmse := RMSE(pred, obs), by=model_id]
    preds[, pooled_delta_rmse := pooled_rmse-RMSE(best_preds$pred, best_preds$obs), by=model_id]
    
    #return the results for each model/fold
    results <- preds[, .(Resample, model_id, cor, pooled_cor, RMSE, pooled_rmse, delta_rmse, pooled_delta_rmse)] %>% 
      unique(., by=key(.))
    
    #choose the best model based on having lowest pooled rmse below the pooled cor cutoff
    results[, second_best := 0]
    second_best_model_id <- NA
    
    #use instead, the best model that is below the 10% pooled correlation value
    xg_ensemble_corr <- quantile(results$pooled_cor, probs=.1)
    results[pooled_delta_rmse == results[pooled_cor < xg_ensemble_corr, pooled_delta_rmse %>% min], 
            `:=` (second_best=1, label=model_id)]
    second_best_model_id <- results[second_best==1, model_id %>% max]
    
    message('Found second best tune: #', second_best_model_id, '\ndelta RMSE=',
            results[second_best==1, pooled_delta_rmse %>% max] %>% round(2), 
            '\n-->model had ', xg_ensemble_corr %>% round(2), ' correlation w/ best')
    
    #extract the tune settings from the selected second best model
    second_best_tune <- preds[model_id==second_best_model_id, names(best_tune), with=F][1]
    message('Second best tune, has parameters: \n')
    summary(second_best_tune)
    
    #rbind the aggregated results too for graphing
    agg_results <- results[, .(Resample='pooled', model_id=model_id,
                               cor=pooled_cor, delta_rmse=pooled_delta_rmse, second_best=second_best, label=label)] %>% 
      unique(., by='model_id') %>% 
      list(., results) %>% 
      rbindlist(use.names = T, fill = T)
    
    #graph the results
    plot <- ggplot(agg_results, aes(x=cor, y=delta_rmse, label=label)) +
      geom_point() +
      geom_vline(xintercept = xg_ensemble_corr, color='red') +
      geom_text_repel(nudge_y = .25, nudge_x=.75) +
      facet_wrap(~Resample) +
      theme_minimal()
    
    ggsave(filename = paste0(stack_dir, region, '_xg_ensemble.png'),
           plot = plot)
    
    # Save model fit objects for future use
    saveRDS(xg_fit, paste0(outputdir, 'stackers/xg_fit_', region, ".RDS"))
    
    # Save the best parameters to csv file
    write.csv(best_tune, paste0(stack_dir, 'xgboost_best_tune_', region, '.csv'))
    write.csv(second_best_tune, paste0(stack_dir, 'xgboost_second_best_tune_', region, '.csv'))
    
  }
  
  #set up control for the ensemble model
  train_control_final <- trainControl(method = "cv",
                                      number = k_folds,
                                      savePredictions = "final")
  
  # if selected, build the ensemble using caretEnsemble
  if (build_ensemble) {
    
    message('building xgb ensemble')
    
    #build ensemble
    xg_ensemble <- caretList(form,
                             data = df,
                             trControl = train_control_final,
                             tuneList = list(xgbTree=caretModelSpec(method="xgbTree", tuneGrid=best_tune),
                                             xgbTree=caretModelSpec(method="xgbTree", tuneGrid=second_best_tune)),
                             metric = "RMSE",
                             objective = objective_function,
                             weights = df$xg_weight)
    
    greedy_ensemble <- caretEnsemble(
      xg_ensemble,
      metric="RMSE",
      trControl=train_control_final)
    
    summary(greedy_ensemble)
    
    #Plot variable importance ---------------------------------------------
    pdf(file=paste0(stack_dir, region, '_xgb_ensemble_covariate_importance.pdf'),
        height=8, width=12)
    lapply(greedy_ensemble$models, function(x) varImp(x) %>% plot)
    dev.off()
    
    # Extract out of sample and in sample predictions
    df[, 'xgboost_cv_pred'   := arrange(greedy_ensemble$ens_model$pred, rowIndex)[,"pred"]]
    df[, 'xgboost_full_pred' := predict(greedy_ensemble, df)]
    
    # Name model for later use in making stack rasters
    greedy_ensemble$model_name <- "xgboost_ensemble"
    
    
    xgboost <- list(dataset = df[, c('xgboost_cv_pred', 'xgboost_full_pred')],
                    xgboost = greedy_ensemble)
    
    #otherwise, return both the child models for ensembling in MBG
  } else {
    
    message("Fitting xgboost children on final best/second best tuned hyperparameters")
    
    xg_fit_final <- train(form,
                          data = df,
                          trControl = train_control_final,
                          verbose = F,
                          tuneGrid = best_tune,
                          metric = "RMSE",
                          method = "xgbTree",
                          objective = objective_function,
                          weights = df$xg_weight)
    xg_fit_final_2 <- train(form,
                            data = df,
                            trControl = train_control_final,
                            verbose = F,
                            tuneGrid = second_best_tune,
                            metric = "RMSE",
                            method = "xgbTree",
                            objective = objective_function,
                            weights = df$xg_weight)
    
    #Plot variable importance ---------------------------------------------
    pdf(file=paste0(stack_dir, region, '_covariate_importance.pdf'),
        height=8, width=12)
    lapply(list(xg_fit_final, xg_fit_final_2), function(x) varImp(x) %>% plot)
    dev.off()
    
    # Extract out of sample and in sample predictions
    df[, 'xgboost_cv_pred'   := arrange(xg_fit_final$pred, rowIndex)[,"pred"]]
    df[, 'xgboost_full_pred' := predict(xg_fit_final, df)]
    df[, 'xgboost2_cv_pred'   := arrange(xg_fit_final_2$pred, rowIndex)[,"pred"]]
    df[, 'xgboost2_full_pred' := predict(xg_fit_final_2, df)]
    
    # Name model for later use in making stack rasters
    xg_fit_final$model_name <- "xgboost"
    xg_fit_final_2$model_name <- "xgboost2"
    
    
    xgboost <- list(
      xgboost <- list(dataset = df[, .(xgboost_cv_pred, xgboost_full_pred)],
                      xgboost = xg_fit_final),
      xgboost2 <- list(dataset = df[, .(xgboost2_cv_pred, xgboost2_full_pred)],
                       xgboost2 = xg_fit_final_2)
    )
    
  }
  
  return(xgboost)
  
}
