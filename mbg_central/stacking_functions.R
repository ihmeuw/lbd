###functions for stacking

#Extract Covariates
#About: Reads in a data table (with latitude and longitude) and a list of rasters/bricks/stacks and extracts the values at the points
#df: a data table. Should have columns named latitude and longitude at a miniumum
#covariate list: a list object of rasters-like formats (e.g. layer, brick, stack) or just one of that raster-like object
#unlike some of the other MBG covariate functions, this extract_covariates function should be able to treat
#raster bricks as a de facto list-- only if 1 brick is passed the function doesn't automatically assume time varying)
#instead, that is implied from the layer names within the brick.
#For best use, provide a list OR a single brick when trying to analyze multiple sets of covariates
#reconcile_timevarying: If true, the functio n enforces period matching. For example, if evi is in 4 periods, there will be one EVI column
#where the values of the column are period/time/year specific. Period in this case is usually identified by
# a .#### at the end of the layer name (usually within a brick). As of 11/2, this is .1, .2, .3 ,.4 but it should
# work for something like .1990, .2000, etc. more testing is needed
#period_var: is there a column in the passed data frame designating the time period that matches the rasters
#return_only_results: return only the covariate columns? (better for cbinding)
extract_covariates = function(df, covariate_list, id_col = NULL, reconcile_timevarying = T, period_var = NULL, return_only_results = F, centre_scale = F, period_map = NULL){
  #Load libraryd packages. Should be redundent

  library('raster')
  library('data.table')
  df = copy(df) #data table does weird scoping things. Do this just in case

  #create an id column
  if(is.null(id_col)){
    df[,rrrr := 1:nrow(df)]
    id_col = 'rrrr'
  }

  # fill in standard 5-year period map if it's missing
  if(is.null(period_map)) period_map <- make_period_map(modeling_periods = c(2000,2005,2010,2015))


  #if covariate list is not actually a list, convert it to a list
  if(class(covariate_list) != 'list'){

    #if a brick or a stack, extract the basename
    if(class(covariate_list) %in% c('RasterBrick', 'RasterStack')){
      #get the list of brick names
      #find the prefixes, remove the .### at the end of the brick
      brick_names = unique(gsub("(\\.[0-9]+)$","",names(covariate_list)))

      #convert the brick into a list the rest of the function expects (e.g. time varying vs. single)
      #assuming I've done this properly, this will convert the raster brick into a list where items in the list are seperated by the prefix
      #(assumes a .### at the end of the name in the raster brick object refers to the period)

      covariate_list = setNames(lapply(brick_names, function(x) covariate_list[[grep(x, names(covariate_list), value =T )]]),brick_names)


    }

  }

  #rename the raster bricks
  #borrowed from some of the covariate functions
  #converts the time varying covariate names into something that places nicer with the extract function
  tv_cov_names = c()
  for(lll in names(covariate_list)){
    if(class(covariate_list[[lll]]) == 'RasterBrick'){
      tv_cov_names = append(tv_cov_names, lll)
      if (nrow(period_map) > 1)  names(covariate_list[[lll]]) = paste0('XXX',1:length(names(covariate_list[[lll]]))) #use XXX as a place holder
    }

  }



  #extract the rasters by points
  locations = SpatialPoints(coords = as.matrix(df[,.(longitude,latitude)]), proj4string = CRS(proj4string(covariate_list[[1]])))

  cov_values = as.data.frame(lapply(covariate_list, function(x) raster::extract(x, locations)))

  #fix the names of the time varying covariates
  if (nrow(period_map) > 1)  names(cov_values) = sub('XXX', '', names(cov_values))

  #right now 4 periods are assumed-- make more flexible later

  #reconcile time varying covariates
  #start by converting the year variable into periods
  #if period is empty -- infer
  if(is.null(period_var)){

    year_map = as.data.frame(sort(unique(df[,year])))
    year_map$period_hold = 1:nrow(year_map)
    names(year_map) = c('year','period_hold')

    df = merge(df, year_map, by = 'year', sort =F) #sorting screws up things relative to cov values
    period_var = 'period_hold'
  }

  #sort the dataset by rid

  setorderv(df, cols = c(id_col))


  if(return_only_results){
    df = df[, c(id_col, 'longitude', 'latitude', period_var), with = F]
  }

  #combine the dataset with the covariate extractions
  df = cbind(df, cov_values)

  #reshape to fill out the columns-- these just get the tv cov names in a way that is needed for two different efforts
  #make this less likey to cause a hiccup
  tv_cov_colnames = grep(paste(tv_cov_names, collapse = "|"), names(df), value = T) #unlisted
  tv_cov_colist = lapply(tv_cov_names, function(x) grep(x, names(df), value = T))

  #Keep only where period of the data matches the period of the covariate for time varying covariates
  if(reconcile_timevarying){
    #reshapes the time varying covariates long
    df = melt(df, id.vars = names(df)[!names(df) %in% tv_cov_colnames], measure = tv_cov_colist, value.name = tv_cov_names, variable.factor =F)

    #melt returns different values of variable based on if its reshaping 1 or 2+ columns.
    #enforce that it must end with the numeric after the period
    df <- df[,variable:= as.numeric(substring(variable,regexpr("(\\.[0-9]+)$", variable)[1]+1))]

    #keep only where the periods match
    df <- merge(df, period_map, by.x = period_var, by.y = 'data_period')
    df = df[period_id == variable, ]

    #clean up
    df = df[,variable := NULL]

  }

  #clean up
  if(period_var=='period_hold'){
    df = df[,period_hold := NULL]
  }

  #centre_scale the results
  if(centre_scale){
    design_matrix = data.frame(df[,names(covariate_list)[!grepl('gaul', names(covariate_list))], with =F])
    cs_df <- getCentreScale(design_matrix)
    design_matrix <- centreScale(design_matrix, df = cs_df)

    #replace the df columns with the design matrix
    df[, names(covariate_list)[!grepl('gaul', names(covariate_list))] := NULL]
    df = cbind(df, design_matrix)

  }

  df <- df[, year := NULL]

  #return the data frame with the year-location specific covariates
  if(return_only_results & reconcile_timevarying){
    df = df[, names(df)[!names(df) %in% c('latitude','longitude')], with = F]
  }

  #sort df just to be sure
  setorderv(df, cols = c(id_col))
  if(id_col == 'rrrr') df[,rrrr := NULL]

  if(centre_scale){
    ## add on rows to covs_cs_df for the gaul_codes
    to.add <- colnames(df)[grepl("gaul_code_", colnames(df))]
    m.add  <- rep(0, length(to.add))
    s.add  <- rep(1, length(to.add))
    to.add <- data.frame(name = to.add, mean = m.add, sd = s.add)
    cs_df <- rbind(cs_df, to.add)
    return(list(covs = df, cs_df = cs_df))
  } else{
    return(df)
  }
}




#for parsing gams string functions
parseArgsS <- function(l) {
  # parse a list of additional arguments to smoothers in gamTrans
  stopifnot(is.list(l))
  l_string <- paste(names(l),
                    lapply(l, addQuotes),
                    sep = ' = ',
                    collapse = ', ')
  return (l_string)
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

  #message(gam_formula)

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


# fit_gbm ---------------------------------------------------------------------------------------------------------------------------
#
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
#basically a wrapper function for the fit_gam function (which is its own wrapper function)
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
                                    xg_model_tune = TRUE,
                                    hyperparameter_filepath = NULL){

  load_R_packages(c("xgboost", "caret"))

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

  if (xg_model_tune == T){
    message("Model tuning xgboost")

    # Set grid search as default unless filepath is provided
    if(is.null(hyperparameter_filepath)){
      message("Tuning with default hyperparameter settings")
      xg_grid <- expand.grid(nrounds = 100,
                             max_depth = c(4, 6, 8, 10, 12),
                             eta = (3:8) / 100,
                             colsample_bytree = .5,
                             min_child_weight = 1,
                             subsample = 1,
                             gamma = 0)
    } else {
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
    # Set cross validation options, default to 5 times repeated 5-fold cross validation
    # Selection function is "oneSE" to pick simplest model within one standard error of minimum
    train_control <- trainControl(selectionFunction = "oneSE",
                                  method = "repeatedcv",
                                  number = 5,
                                  repeats = 5)
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

    # Save model fit object for future use
    saveRDS(xg_fit, paste0(outputdir, 'stackers/xg_fit_', region, ".RDS"))

    # Create plot showing which hyperparameters were selected
    cv_index <- length(xg_fit$control$index)
    xg_fit$results <-
      xg_fit$results %>%
      mutate(RMSE_low  = RMSE - (RMSESD / sqrt(cv_index)),
             RMSE_high = RMSE + (RMSESD / sqrt(cv_index)))

    gg1 <-
      ggplot(xg_fit$results,
             aes(x = eta,
                 y = RMSE,
                 color = factor(max_depth),
                 shape = factor(max_depth))) +
      geom_line() +
      geom_point() +
      xlim(range(xg_fit$results$eta)[1], range(xg_fit$results$eta)[2]) +
      labs(x = "Learning rate") +
      facet_wrap(~ nrounds) +
      theme_bw() +
      scale_color_discrete("Max Tree Depth") +
      scale_shape_discrete("Max Tree Depth")

    error_width <- diff(range(xg_fit$results$eta)) / 50

    gg1 <- gg1 +
      geom_errorbar(data = xg_fit$results[oneSE(xg_fit$results, "RMSE", cv_index, F),],
                    aes(x = eta, ymin = RMSE_low, ymax = RMSE_high),
                    alpha = 0.8,
                    width = error_width,
                    size = 0.5,
                    color = "black")

    ggsave(filename = paste0(stack_dir, region, '_hyperparameter.png'),
           plot = gg1)

    # Save the best parameters to csv file
    write.csv(xg_fit$bestTune, paste0(stack_dir, 'xgboost_best_tune_', region, '.csv'))
  }

  if (xg_model_tune == T) {
    # Extract best parameters
    xg_best_tune <- read.csv(paste0(stack_dir, 'xgboost_best_tune_', region, '.csv'))
  } else {
    # Extract best parameters from filepath
    xg_best_tune <- read.csv(hyperparameter_filepath)
  }

  # Define grid search and
  xg_grid_final <- expand.grid(nrounds = xg_best_tune$nrounds,
                               max_depth = xg_best_tune$max_depth,
                               eta = xg_best_tune$eta,
                               colsample_bytree = .5,
                               min_child_weight = 1,
                               subsample = 1,
                               gamma = 0)

  train_control_final <- trainControl(method = "cv",
                                      number = 5,
                                      savePredictions = "final")

  message("Fitting xgboost on final tuned hyperparameters")
  xg_fit_final <- train(form,
                        data = df,
                        trControl = train_control_final,
                        verbose = F,
                        tuneGrid = xg_grid_final,
                        metric = "RMSE",
                        method = "xgbTree",
                        objective = objective_function,
                        weights = df$xg_weight)

  # Plot the covariate importance of final model
  cov_plot <-
    ggplot(varImp(xg_fit_final, scale = FALSE)) +
    labs(x = "Covariate", y = "Relative Importance") +
    theme_bw()
  ggsave(filename = paste0(stack_dir, region, '_covariate_importance.png'),
         plot = cov_plot)

  # Extract out of sample and in sample predictions
  df[, 'xgboost_cv_pred'   := arrange(xg_fit_final$pred, rowIndex)[,"pred"]]
  df[, 'xgboost_full_pred' := predict(xg_fit_final, df)]

  # Name model for later use in making stack rasters
  xg_fit_final$model_name <- "xgboost"

  xgboost <- list(dataset = df[, c('xgboost_cv_pred', 'xgboost_full_pred')],
                  xgboost = xg_fit_final)
  return(xgboost)
}


#fetch covariate layer
#given a raster-like object and a period, returns the appropriate layer, assuming chronological order
fetch_covariate_layer = function(ras, period = 1){
  if(class(ras) == 'RasterBrick' | class(ras) == "RasterStack"){
    return(ras[[period]])
  } else{
    return(ras)
  }
}

#create a stacked prediction for each period and child model. Function returns a list of raster objects. First object is a brick of stacked results. Additional objects refer to the child models
#covariate_layers: list of raster objects that consititue the covs. Must share names with the DT that the models were fit on
#period: what period should be used. 1:N -- as of 11/2 we assume 4
#child_models: a list of model objects from the child models
#return children: return the prediction rasters from the child models as well
#todo: allow dataframes of constants to be passed
produce_stack_rasters = function(covariate_layers = all_cov_layers, #raster layers and bricks
                                 period = 1, #period of analysis
                                 child_models = list(),
                                 indicator_family = 'binomial',
                                 return_children = F,
                                 centre_scale_df = NULL){ #stacker model

  message(paste0('The Period is ', period))

  #fetch the covariates appropriate for the period
  period_covs =brick(lapply(covariate_layers, function(x) fetch_covariate_layer(x,period)))

  #create constants -- only flexible for year/period
  year = data.frame(year = (1995 +period*5))

  #predict the various models. This is super strict with variable names (and beware scoping issues with the names)
  #brick the results
  stacker_predictors <- brick(lapply(child_models, function(x) predict_model_raster(x,
                                                                                    period_covs,
                                                                                    constants = year,
                                                                                    indicator_family = indicator_family,
                                                                                    centre_scale_df = centre_scale_df
  )))

  return(stacker_predictors)
}


#Wrapper function for produce stacked rasters
#Runs all (or a subset of periods) and formats the results.
make_stack_rasters <- function(covariate_layers = all_cov_layers,
                               period = NULL,
                               child_models = list(),
                               indicator_family = 'binomial',
                               rd = run_date,
                               re = reg,
                               ind_gp = indicator_group,
                               ind = indicator,
                               ho = holdout,
                               centre_scale_df = NULL,
                               ...) {

  ## first, save child_models for use in get.cov.wts in post_estiamtion if not other places too
  save(child_models, file='<<<< FILEPATH REDACTED >>>>')

  if(is.null(period)){
    period = 1:4
  }

  res = lapply(period, function(the_period) produce_stack_rasters(covariate_layers = covariate_layers, #raster layers and bricks. Covariate rasters essentially
                                                                  period = the_period, #period of analysis. Fed in from the lapply
                                                                  child_models = child_models, #a list of model objects for the child learners
                                                                  indicator_family = indicator_family,
                                                                  centre_scale_df = centre_scale_df))
  
  ## Prep the rasters to be ordered in the following order:
  ## raster_brick[[child_models]][[period]]
  ret_obj <- lapply(names(child_models), function(x_cn) {
    raster::brick(lapply(period, function(x_t) res[[x_t]][[x_cn]]))
  })
  
  ## Set names of the output raster list
  names(ret_obj) <- names(child_models)

  # Make sure that the raster brick dimensions and names are all correct
  stopifnot(assertthat::are_equal(length(ret_obj), length(names(child_models))))
  
  j <- 0
  for (x_t in ret_obj) {
    j <- j + 1
    stopifnot(assertthat::are_equal(names(x_t), paste0(names(child_models)[j], ".", period)))
  }

  return(ret_obj)
}

# predict_model_raster -----------------------------------------------------------------------------------
#
#predict model raster: given a model object and some covariates, predict out the results in raster form
#model call: the model object you want to create rasters from
#covariate_layers: a list of raster-like objects named identically to the models fit on the tabular data
#constants: do any constants need to be fed to the prediction step?
#indicator_family: model family. Specifies what sorts of transformations need doing
#to do: add a transform switch

predict_model_raster = function(model_call,
                                covariate_layers,
                                constants = NULL,
                                indicator_family = 'binomial',
                                centre_scale_df = NULL){
  
  
  message(paste0('predicting out:', model_call$model_name))
  #convert the raster objects into a named matrix
  dm = as.data.frame(stack(covariate_layers),xy =T)
  
  #apply centre scaling
  if(!is.null(centre_scale_df)){
    cs_dm = centreScale(dm[,names(covariate_layers)], df = centre_scale_df)
    dm = cbind(dm[,c('x','y')], cs_dm)
  }
  
  orig_names = names(dm)
  
  #add constants(if present) to the raster
  if(!is.null(constants)){
    dm = cbind(dm, constants)
  }
  #create a template
  dm$rid = 1:nrow(dm)
  template = dm[,c('rid','x','y')]
  
  #drop rows with NA data
  dm = na.omit(dm)
  
  
  #if a gam or a bam
  #class(model_call)
  if(inherits(model_call, 'gam') | inherits(model_call, 'bam')){
    
    ret_obj = predict(model_call, newdata =dm, type = 'response')
    
  } else if(inherits(model_call, 'gbm')){
    
    ret_obj = predict(model_call, newdata=dm, n.trees = model_call$n.trees, type = 'response')
    
  } else if(inherits(model_call, 'glmnet')){
    #glmnet is dumb and wants a matrix for new data
    
    
    dm_n = names(dm)
    dm = as.matrix(dm)
    colnames(dm) = dm_n
    
    #predict(object, newx, s = NULL, type=c("link","response","coefficients","nonzero","class"), exact = FALSE, offset, ...)
    
    ret_obj = predict(model_call, newx = data.matrix(dm[,rownames(model_call$beta)]), s=model_call$cv_1se_lambda, type = 'link')
    
    #backtransform into probability/percentage space
    if( indicator_family=='binomial'){
      ret_obj = invlogit(ret_obj)
    }
    
    
    #return dm to its data frame form
    dm = data.frame(dm)
    colnames(ret_obj) = 'ret_obj'
  } else if(inherits(model_call, 'randomForest')){
    
    ret_obj = predict(model_call, newdata=dm, type = 'response')
    
    #back transform if binomial
    if( indicator_family=='binomial'){
      ret_obj = invlogit(ret_obj)
    }
    
  } else if(inherits(model_call, 'earth')){
    ret_obj = predict(model_call, newdata=dm, type = 'response')
    colnames(ret_obj) ='ret_obj'
    
  } else if (inherits(model_call, 'train')) {
    # All caret objects use this framework
    ret_obj = predict(model_call, newdata = dm)
  }
  
  #convert back to a raster
  
  #rasterFromXYZ(xyz, res=c(NA,NA), crs=NA, digits=5)
  ret_obj = cbind(data.frame(rid = dm[,c('rid')]), ret_obj = ret_obj)
  
  #restore to former glory
  ret_obj= merge(template, ret_obj, by = 'rid', all.x =T)
  setorder(ret_obj, rid)
  ret_obj = rasterFromXYZ(ret_obj[,c('x','y','ret_obj')], res = res(covariate_layers[[1]]), crs=crs(covariate_layers[[1]]))
  
  #return the object
  return(setNames(ret_obj,model_call$model_name))
}

#A function to stack the estimates from child modules using GLM
#df: dataset. It must have the *_cv_pred and *_full_pred columns from the child models
#indicator: name of the outcome/indicator column in the dataset
#indicator_family: family that the stacker should be fit with
#NOTE: This function is somewhat out of date. Use the gam stacker instead for now.
glm_stacker = function(df, #the dataset in data table
                       model_names=c('gam','gbm'), #prefixes of the models to be stacked
                       indicator = indicator, #the indicator of analysis
                       indicator_family = indicator_family){ #indicator family (e.g. binomial)
  ##start function##

  #copy dataset to avoid weird data table scoping issues
  df = copy(df)

  #format the outcome variable depending on the family
  if(indicator_family == 'binomial'){
    df[, failures := N-get(indicator)] #failures in the sense they didn't get sick
    outcome = df[, .(get(indicator),failures)]
    names(outcome)[1] = indicator
  } else{

    outcome = df[,.(get(indicator))]
    names(outcome)[1] = indicator
  }

  outcome = as.matrix(outcome)

  #format the child model results into the glm format
  #create new columns to hold the cv results, which are then replaced with the full results upon prediction
  df[, (model_names) := lapply(model_names, function(mn) get(paste0(mn,'_cv_pred')))]

  #collapse the model names to a basic formula
  glm_formula = as.formula(paste0('outcome~',paste(model_names, collapse = '+')))
  stacker = glm(glm_formula, family = indicator_family, data = df)

  #predict the results as fit from the crossvalidated stuff
  df[,stacked_cv_pred := predict(stacker, df, type = 'response')]

  #overwrite the columns to work on the full fit child modules
  df[, (model_names) := lapply(model_names, function(mn) get(paste0(mn,'_full_pred')))]
  df[,stacked_pred := predict(stacker, df, type = 'response')]

  #return the dataframe and the stacker model
  return(setNames(list(df[,stacked_pred], stacker),c('dataset','stacker_model'))) #[,.(stacked_pred)]

}

#Stack models using GAM/BAM
#df: the data frame
#indicator: indicator of analysis
#indicator_family: what is the analytical model of the family
#bam: should bam be used?
#spline_args: spline parameters passed to GAM/BAM function
gam_stacker = function(df, #the dataset in data table
                       model_names=c('gam','gbm'),
                       weight_column = NULL, #prefixes of the models to be stacked
                       bam = F, #whether or not bam should be used
                       spline_args = list(bs = 'ts', k = 3),
                       indicator = indicator, #the indicator of analysis
                       indicator_family = indicator_family,
                       cores = 'auto'){
  ##start function##
  
  #copy dataset to avoid weird data table scoping issues
  df = copy(df)
  
  #format the outcome variable depending on the family
  if(indicator_family == 'binomial'){
    df[, failures := N-get(indicator)] #failures in the sense they didn't get sick
    outcome = df[, .(get(indicator),failures)]
    names(outcome)[1] = indicator
  } else{
    
    outcome = df[,.(get(indicator))]
    names(outcome)[1] = indicator
  }
  
  outcome = as.matrix(outcome)
  
  #format the child model results into the glm format
  #create new columns to hold the cv results, which are then replaced with the full results upon prediction
  df[, (model_names) := lapply(model_names, function(mn) get(paste0(mn,'_cv_pred')))]
  
  #Fit the gam as the stacker
  stacker = fit_gam(df = df,
                    covariates = paste(model_names, collapse = ' + '),
                    additional_terms = NULL,
                    weight_column = NULL,
                    bam = bam,
                    spline_args = spline_args,
                    auto_model_select =T,
                    indicator = indicator,
                    indicator_family = indicator_family,
                    cores = cores)
  
  stacker$model_name = 'gam_stacker'
  
  #predict the results as fit from the crossvalidated stuff
  df[,stacked_cv_pred := predict(stacker, df, type = 'response')]
  
  #overwrite the columns to work on the full fit child modules
  df[, (model_names) := lapply(model_names, function(mn) get(paste0(mn,'_full_pred')))]
  df[,stacked_pred := predict(stacker, df, type = 'response')]
  
  #return the dataframe and the stacker model
  return(setNames(list(df[,stacked_pred], stacker),c('dataset','stacker_model'))) #[,.(stacked_pred)]
  
}

#fit random forests
fit_rf = function(df, covariates = all_fixed_effects, additional_terms = NULL, ntree = 1000, indicator, indicator_family = 'binomial'){
  df = copy(df)

  #add additional terms if requested
  the_covs = format_covariates(add_additional_terms(covariates,additional_terms))

  #random forest doesn't have a binomial or possion option. Use emperical logit (this also keeps it consistent with other methods that return logit probabilities)
  #create outcome variable
  if(indicator_family == 'binomial' | indicator_family == 'poisson'){
    #message('emplogit')
    df[, y := emplogit(get(indicator), N)]
    #message(names(df))
  } else{
    df[,y:= get(indicator)]
  }


  #fit a random forest
  message(paste0('Fitting Random Forest with ntree: ', ntree))
  model = randomForest(x = df[, the_covs, with = F] ,y = df[,y] , ntree = ntree)

  return(model)


}
#fit a random forest model for stacking
#df: data table
#model_name: what do you want to call this?
#covariates: formula of the fixed effects
#fold_id_col: what is the column iding the fold
#additional_terms: constants, other covarites. Usually only used for year and other non-raster covariates.
#indicator_family: analyitical family
#indicator: what indicator
#cores: how many cores can be used?
fit_rf_child_model = function(df,model_name = 'rf', fold_id_col = 'fold_id', covariates = all_fixed_effects, additional_terms = NULL, ntree = 1000, indicator = indicator, indicator_family = indicator_family, cores = 'auto'){
  library(randomForest)
  library(parallel)
  df = copy(df)
  message('Fitting the Full Random Forest')

  #format covariate string
  the_covs = format_covariates(add_additional_terms(covariates,additional_terms))

  full_model = fit_rf(df,
                      covariates = covariates,
                      additional_terms = additional_terms,
                      ntree = ntree,
                      indicator = indicator,
                      indicator_family = indicator_family)

  full_model$model_name = model_name

  #fit the child/cv rfs
  folds = unique(df[,get(fold_id_col)])

  message("Fitting baby forests in parallel")
  # Set multithreading to serial for `mclapply()`:
  set_serial_threads()
  # Determine appropriate number of cores to use in `mclapply()`
  if(cores == 'auto') cores <- get_max_forked_threads(nobjs = length(folds))
  #fit all the baby forests at once
  baby_models = mclapply(folds, function(fff) fit_rf(df[get(fold_id_col) != fff,],
                                                     covariates = covariates,
                                                     additional_terms = additional_terms,
                                                     ntree = ntree,
                                                     indicator = indicator,
                                                     indicator_family = indicator_family), mc.cores = cores)
  # Return to multithreading (if any):
  set_original_threads()
  #predict iteratively for simplification. I don't think RF likes changing variable names
  for(fff in folds){
    #message(paste0('Predicting Fold: ', fff))
    #sub_data = pred_filler[get(fold_id_col)==fff,]
    #le_preds = predict(baby_forests[[fff]], newdata = sub_data, type = 'response')
    df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_models[[fff]], newdata = df[get(fold_id_col)==fff,the_covs, with = F], type = 'response')]
  }

  #predict using full model fit earlier
  df[,paste0(model_name,'_full_pred') := predict(full_model,df[,the_covs, with = F],type = 'response')]

  #return a subset of the columns. Full pred denotes the fit from the full model. CV pred is the OOS stuff
  suffixes = c('_full_pred', '_cv_pred')
  return_cols = paste0(model_name, suffixes)

  #if binomial, take the invlogit
  if( indicator_family=='binomial'){
    df[,return_cols[1] := invlogit(get(return_cols[1]))]
    df[,return_cols[2] := invlogit(get(return_cols[2]))]
  }

  #print(return_cols)
  #set up with for the return call
  return(setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name))))
}


#format covariate string
format_covariates = function(covariates = covariates){
  if(length(covariates)==1){
    #split back into parts
    covariates = unlist(tstrsplit(covariates,"\\+"))
  }
  covariates = trimws(covariates)
  return(covariates)
}

#Take the emperical logit
emplogit = function(success, N, epsilon = NULL) {
  #http://stats.stackexchange.com/questions/109702/empirical-logit-transformation-on-percentage-data
  #squeeze in the edges
  tform = success/N

  #if epsilon is null, follow the instructions from the attached link
  if(is.null(epsilon)){
    epsilon = min(tform[tform >0 & tform <1])/2
  }

  tform[tform==0] = tform[tform == 0] + epsilon
  tform[tform==1] = tform[tform==1] - epsilon
  tform = log(tform/(1-tform))

  return(tform)
}

#fit penalized regressions(lasso, elastic net, ridge)
#df: the data, in data table format
#covariates: string in formula notation denoting the covariates
#additional_terms: other covariates to include
#alpha: alpha parameter for glmnet calls
#parallel: TRUE/FALSE to turn on/off parallelization
#offset: is there an offset, valid only for poission models. Not fully implemented
fit_glmnet = function(df, covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, alpha = 1, indicator, indicator_family = 'binomial', parallel = FALSE){

  library(glmnet)
  library(doParallel)
  df = copy(df)

  #add additional terms if requested
  the_covs = format_covariates(add_additional_terms(covariates,additional_terms))

  #format weights
  if(!is.null(weight_column)){
    data_weights = df[,get(weight_column)]
  } else {
    data_weights = rep(1,nrow(df))
  }

  #create outcome object
  #specifying a binomial object is annoyingly difficult in glmnet. Use emplogit instead.
  if(indicator_family == 'binomial'){
    not_indicator = df[,N] - df[,get(indicator)]
    response_var = cbind(not_indicator, outcome = df[,get(indicator)])
  } else if (indicator_family == 'poisson') {
    if(!is.null(offset)){
      stop('fit_glmnet does not currently work for offset poissons')
    } else {
      response_var = df[,get(indicator)]
    }
  } else {
    response_var = df[,get(indicator)]
  }


  #create design matrix
  dm = as.matrix(df[,the_covs, with = F])
  colnames(dm) = the_covs

  #search for lambda
  #these models are run as gaussian because of the prior transformation.
  message(paste0('Fitting glmnet with alpha: ', alpha))

  cv_res = cv.glmnet(x = dm , y= response_var, family = indicator_family, alpha = alpha, weights = data_weights, parallel = parallel)

  #fit full model using selected lambdas
  model = glmnet(x = dm , y= response_var, family = indicator_family, lambda = cv_res$lambda, alpha = alpha, weights = data_weights)

  #preserve the cv_1se_lambda
  model$cv_1se_lambda = cv_res$lambda.1se

  return(model)

}

#Fit a glmnet child model
#df: the data in data table format
#model_name: what do you want to call this?
#covariates: formula of the fixed effects
#fold_id_col: what is the column iding the fold
#additional_terms: constants, other covarites. Usually only used for year and other non-raster covariates.
#indicator_family: analyitical family
#indicator: what indicator
#parallel: TRUE/FALSE to turn on/off parallelization
#alpha: paramter for glmnet
fit_glmnet_child_model = function(df,model_name = 'glmnet', fold_id_col = 'fold_id', covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, alpha = 1, indicator = indicator, indicator_family = 'binomial', cores = 'auto'){
  library(glmnet)
  library(doParallel)
  df = copy(df)
  message('Fitting the Full GLMNET')

  the_covs = format_covariates(add_additional_terms(covariates, additional_terms))

  # glmnet has it's own internal parallelization that we don't have much control
  # over. It is either on or off and is parallel over each fold supplied by the
  # 'nfolds' argument to `cv.glmnet()`. Since we use the default of 10 for
  # 'nfolds' in `cv.glmnet()` in the `fit_glmnet()` function (i.e. we don't pass
  # in an 'nfolds' argument), we have to make sure we at least have 10 cores to
  # use before we turn the parallel option on. See this:
  # https://www.rdocumentation.org/packages/glmnet/versions/2.0-16/topics/cv.glmnet

  # Also, `cv.glmnet()` has caused multithreaded operations to hang similar to
  # `mclapply()`, so we set multithreaded operations to serial here
  set_serial_threads()
  if(cores == 'auto') cores <- get_total_threads()
  if(cores >= 10) {
    parallel <- TRUE
    # start the cluster if parallel
    # Apparently, with Linux systems we don't need the `cl <- makeCluster(cores)`
    # step. Doing it that way will also work, but it is much slower.
    # Do it here, so we don't have to spin up/down the cluster with each
    # `fit_glmnet()` call.
    registerDoParallel(cores = cores)
  } else parallel <- FALSE

  #fit the full model
  full_model = fit_glmnet(df,
                          covariates = covariates,
                          additional_terms = additional_terms,
                          weight_column = weight_column,
                          alpha = alpha,
                          indicator = indicator,
                          indicator_family = indicator_family,
                          parallel = parallel)

  full_model$model_name = model_name

  #fit the child/cv rfs
  folds = unique(df[,get(fold_id_col)])

  for(fff in folds){
    #message(paste0('Fitting and Predicting Fold: ', fff))
    baby_model = fit_glmnet(df[get(fold_id_col) != fff,],
                            covariates = covariates,
                            additional_terms = additional_terms,
                            weight_column = weight_column,
                            alpha = alpha,
                            indicator = indicator,
                            indicator_family = indicator_family,
                            parallel = parallel)

    new_data = df[get(fold_id_col)==fff,the_covs, with = F]

    n_nd = names(new_data)
    new_data = as.matrix(new_data)
    names(new_data) = n_nd

    df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_model, newx = new_data, s = baby_model$cv_1se_lambda, type = 'link')]

  }

  #stop the cluster just in case
  stopImplicitCluster()

  # Return to multithreading (if any):
  set_original_threads()

  #predict using full model fit earlier
  new_data = df[,the_covs, with = F]

  n_nd = names(new_data)
  new_data = as.matrix(new_data)
  names(new_data) = n_nd
  df[,paste0(model_name,'_full_pred') := predict(full_model,newx = new_data, s = full_model$cv_1se_lambda, type = 'link')]

  #return a subset of the columns. Full pred denotes the fit from the full model. CV pred is the OOS stuff
  suffixes = c('_full_pred', '_cv_pred')
  return_cols = paste0(model_name, suffixes)

  #if binomial, undo the logit
  if( indicator_family=='binomial'){
    df[,return_cols[1] := invlogit(get(return_cols[1]))]
    df[,return_cols[2] := invlogit(get(return_cols[2]))]
  }


  #set up with for the return call
  return(setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name))))
}


#fit earth: a function to fit a multivariate adaptive regression splines
#df: data table with the outcome/indicator and some covariates already extracted. This is different than the gam_cov functions
#covariates: a  vector of covariate names and/or formula in the style of all_fixed_effects (e.g. cov1 + cov2 + cov3)
#additional terms: a vector or single character of column names to be included in the model fit
#weight_column: in df, is there a column that specifies the observation weights?
#indicator: name of the column of the DV
#indicator_family: model family
fit_earth = function(df, covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, indicator, indicator_family = 'binomial'){
  library(earth)
  library(data.table)
  #also requires seeg

  df = copy(df) #in case data table scoping gets wonky

  the_covs = format_covariates(add_additional_terms(covariates, additional_terms))

  #set response variable
  if(indicator_family=="binomial") response <- cbind(success = df[, get(indicator)], failure = df[, N] - df[, get(indicator)])
  if(indicator_family=="gaussian") response <- cbind(outcome = df[, get(indicator)])

  #sort out weights
  #format weights
  if(!is.null(weight_column)){
    df[,data_weight := get(weight_column)]
  } else{
    df[,data_weight := 1]
  }
  weight_column = 'data_weight'

  #fit the earth
  message(paste0('Fitting earth'))

  model = earth(x = df[,the_covs,with = F], y = response, weights = df[,get(weight_column)], glm = list(family =indicator_family))

  #return the earth object
  return(model)

}

#fit_earth_child_model: a function to fit a multivariate adaptive regression splines, with the kfold crossval
#df: data table with the outcome/indicator and some covariates already extracted. This is different than the gam_cov functions
#model_name: model name
#fold_id_col: What column identifies the folds
#covariates: a  vector of covariate names and/or formula in the style of all_fixed_effects/right hand side of a formula (e.g. cov1 + cov2 + cov3)
#additional terms: a vector or single character of column names to be included in the model fit
#weight_column: in df, is there a column that specifies the observation weights?
#indicator: name of the column of the DV
#indicator_family: model family
#cores: # of cores are available for use
fit_earth_child_model = function(df, model_name = 'earth', fold_id_col = 'fold_id', covariates = all_fixed_effects,
                                 additional_terms = NULL, weight_column = NULL, indicator = indicator, indicator_family = 'binomial', cores = 'auto'){
  library(earth) #load the earth!
  library(parallel)

  #remove scoping surprises
  df = copy(df)

  the_covs = format_covariates(add_additional_terms(covariates, additional_terms))

  #start by fitting the full gam
  message('Fitting the Full earth model')
  full_model = fit_earth(df = df,
                         covariates = covariates,
                         additional_terms = additional_terms,
                         weight_column = weight_column,
                         indicator = indicator,
                         indicator_family = indicator_family)

  #add a name to the game object
  full_model$model_name = model_name

  #fit the child/cv gams
  message("Fitting baby earth")
  #fit the child/cv rfs
  folds = unique(df[,get(fold_id_col)])

  # Set multithreading to serial for `mclapply()`:
  set_serial_threads()
  # Determine appropriate number of cores to use in `mclapply()`
  if(cores == 'auto') cores <- get_max_forked_threads(nobjs = length(folds))
  baby_models = mclapply(folds, function(fff) fit_earth(df = df[get(fold_id_col) != fff,],
                                                        covariates = covariates,
                                                        additional_terms = additional_terms,
                                                        weight_column = weight_column,
                                                        indicator = indicator,
                                                        indicator_family = indicator_family), mc.cores = cores)
  # Return to multithreading (if any):
  set_original_threads()

  for(fff in folds){
    #use the data fit on K-1 of the folds to fit on the help out fold
    df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_models[[fff]], df[get(fold_id_col)==fff,the_covs, with = F],type='response')]
  }

  #predict using full model fit earlier
  df[,paste0(model_name,'_full_pred') := predict(full_model,df[,the_covs,with = F],type = 'response')]

  #return a subset of the columns. Full pred denotes the fit from the full model. CV pred is the OOS stuff
  suffixes = c('_full_pred', '_cv_pred')
  return_cols = paste0(model_name, suffixes)
  #print(return_cols)
  #set up with for the return call
  return(setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name))))
}

add_additional_terms= function(fes, add_terms = NULL){
  new_list = paste(unlist(sapply(c(fes, add_terms), function(x) format_covariates(x))), collapse = ' + ')
  return(new_list)
}

#' Create national estimates from stackers by population-weighting the pixel estimates.
#'
#' @param indicator indicator
#' @param indicator_group indicator_group
#' @param run_date model run date
#' @param age age group
#' @param holdout model holdout number
#' @param regions vector of model regions
#' @param year_list vector of modeled years
#' @param pop_measure population measure to use for population-weighting
#' @param stackers_logit_transform logical. Were stackers logit-transformed when fitting the final model?
#' @param predictor_logit_transform logical. Was the linear predictor logit-transformed (eg, as is
#'        typically the case in binomial models)?
#' @param results_file file path to a .csv file where results should be saved. Optional.
#' @param shapefile_version string specifies which shapefile version to pull
#'
#' @return A table of national-level estimates for each country, year, and stacker. If results_file
#'         is provided, this information is saved to the specified file. If results_file is not
#'         provided, this information is returned as a data.table.
#'

aggregate_stackers_admin0 <- function(indicator, indicator_group, run_date, age, holdout,
                                      regions, year_list, pop_measure, stackers_logit_transform,
                                      predictor_logit_transform, results_file = NULL,
                                      shapefile_version = 'current') {

  # Loop over regions
  all_regions <- lapply(regions, function(reg) {
    message(paste("Working on", reg, "\n"))
    pathaddin <- paste0('_bin', age, '_', reg, '_', holdout)

    # Get simple polygon and simple raster
    poly <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version),
                                buffer = 0.4, subset_only = FALSE,
                                shapefile_version = shapefile_version)
    subset_shape <- poly[[1]]
    simple_polygon <- poly[[2]]
    simple_raster <- build_simple_raster_pop(subset_shape)[["simple_raster"]]
    pixel_id <- seegSDM:::notMissingIdx(simple_raster)

    # Get population
    pop <- load_and_crop_covariates_annual(covs = 'worldpop', measures = pop_measure, simple_polygon = simple_polygon,
                                           start_year = min(year_list), end_year = max(year_list), interval_mo = 12, agebin=1)
    pop <- crop_set_mask(pop[[1]], simple_raster)
    pop <- data.table(raster::extract(pop, pixel_id))
    pop[, pixel_id := pixel_id]
    pop <- melt(pop, id.vars = "pixel_id", variable.name = "year", value.name = "pop")
    pop[, year := min(year_list) + as.numeric(year) - 1]
    pop[is.na(pop), pop:=0]

    # Get stackers
    load(paste0("<<<< FILEPATH REDACTED >>>>", run_date, pathaddin, ".RData"))
    cov_list <- cov_list[child_model_names]
    cov_list <- lapply(cov_list, function(x) raster::extract(x, pixel_id))
    cov_list <- data.table(do.call("cbind", cov_list))
    cov_list[, pixel_id := pixel_id]
    cov_list <- melt(cov_list, id.vars = "pixel_id", measure=patterns(child_model_names),
                     variable.name = "year", value.name = child_model_names)
    cov_list[, year := min(year_list) + as.numeric(year) - 1]

    # Get combined stackers result based on INLA coefficients
    load(paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_model_eb", pathaddin, ".RData"))
    if ("sdrep" %in% names(res_fit)) { # TMB, w/ stackers as fixed effects
      coef <- res_fit$sdrep$par.fixed[names(res_fit$sdrep$par.fixed) == "alpha_j"]
      names(coef) <- res_fit$fenames
      int <- coef["int"]
      coef <- coef[child_model_names]
    } else if ("covar" %in% names(res_fit$summary.random)) { # INLA, w/ stackers as random effects (to sum to one)
      coef <- res_fit$summary.random$covar$"mean"
      names(coef) <- child_model_names
      int <- res_fit$summary.fixed["int", "mean"]
    } else { # INLA, w/ stackers as fixed effects
      coef <- res_fit$summary.fixed[child_model_names, "mean"]
      names(coef) <- child_model_names
      int <- res_fit$summary.fixed["int", "mean"]
    }

    if (stackers_logit_transform) {
      cov_list[, combined := int + Reduce('+', lapply(child_model_names, function(x) logit(cov_list[[x]]) * coef[x]))]
    } else {
      cov_list[, combined := int + Reduce('+', lapply(child_model_names, function(x) cov_list[[x]] * coef[x]))]
    }

    if (predictor_logit_transform) cov_list[, combined := inv.logit(combined)]

    # Merge simple_raster, stackers, and population
    admin0 <- data.table(pixel_id = pixel_id, ADM0_CODE = raster::extract(simple_raster, pixel_id))
    admin0 <- merge(admin0, cov_list, by="pixel_id")
    admin0 <- merge(admin0, pop, by=c("pixel_id", "year"))
    if (nrow(admin0) != nrow(pop) | nrow(admin0) != nrow(cov_list)) stop("dimension mismatch")

    # Calculate population-weighted national averages
    admin0[, lapply(.SD, function(x) weighted.mean(x, pop, na.rm=T)), by='ADM0_CODE,year', .SDcols = c(child_model_names, "combined")]
  })

  # Save or return results for all regions
  all_regions <- rbindlist(all_regions)
  if (!is.null(results_file)) {
    write.csv(all_regions, file = results_file, row.names=F)
    return(paste("results saved to:", results_file))
  } else {
    return(all_regions)
  }
}
