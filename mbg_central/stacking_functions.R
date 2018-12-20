###functions for stacking

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

  #the periods thing is kind of dumb. for now, confirm that the year variable is in four periods to match the raster bricks
  # if(length(unique(df[,year]))!=4){
  #   stop("Your input data frame's year column is not in the 4 period mode to match the 4 periods of likely raster bricks \n Please fix and rerun")
  # }


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
      names(covariate_list[[lll]]) = paste0('XXX',1:length(names(covariate_list[[lll]]))) #use XXX as a place holder
    }

  }



  #extract the rasters by points
  locations = SpatialPoints(coords = as.matrix(df[,.(longitude,latitude)]), proj4string = CRS(proj4string(covariate_list[[1]])))

  cov_values = as.data.frame(lapply(covariate_list, function(x) extract(x, locations)))

  #fix the names of the time varying covariates
  names(cov_values) = sub('XXX', '', names(cov_values))

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

    #if(length(tv_cov_colist)==1 & (packageVersion("data.table") < package_version("1.10.0"))){
    #  setnames(df,paste0(tv_cov_names[[1]],'1'), tv_cov_names)
    #}

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
#stolen from the seeg stuff
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
fit_gam = function(df, covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, bam = F, spline_args = list(), auto_model_select =F, indicator, indicator_family = 'binomial', cores = 1){
  library(mgcv)
  library(data.table)
  #also requires seeg

  df = copy(df) #in case data table scoping gets wonky

  #check to see if the gam formula is prespecified. if not, make it
  #check to see if its a vector of characters or a psudo-formula
  covariates = format_covariates(covariates) #additional terms is handled below


  #remove binary vars from the splined versions and set as additional terms
  n_uniq_vals = unlist(setNames(df[,lapply(covariates, function(x) uniqueN(get(x)))], covariates)) #this needs to return a named vector
  binary_vars = names(n_uniq_vals[n_uniq_vals<=2])
  covariates = covariates[!covariates %in% binary_vars]

  additional_terms = c(additional_terms,binary_vars)
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

  #fit the gam
  message(paste0('Fitting GAM/BAM with spline args of: ', names(spline_args)[1],'=',spline_args[1],' ', names(spline_args)[2],'=', spline_args[2]))
  if(bam){
    if(auto_model_select ==T & nrow(df)<2000){
      model = mgcv::gam(gam_formula, data = df, family = indicator_family, weights = df[,get(weight_column)], control = list(nthreads = as.numeric(cores)))
    } else{
      model = mgcv::bam(gam_formula, data = df, family = indicator_family, weights = df[,get(weight_column)], nthreads = as.numeric(cores), discrete = T)
    }
  } else{
    model = mgcv::gam(gam_formula, data = df, family = indicator_family, weights = df[,get(weight_column)], control = list(nthreads = as.numeric(cores)),
                      method = "REML")
  }

  #return the gam object
  return(model)

}

#a wrapper function for fitting gam/bam models in the stacking framework. It'll run 1+k times internally
#arguments are the same as fit gam above.
#basically a wrapper function for the fit_gam function (which is its own wrapper function-- hooray for rabbit holes)
fit_gam_child_model = function(df, model_name = 'gam', fold_id_col = 'fold_id', covariates = all_fixed_effects,
                                additional_terms = NULL, weight_column = NULL, bam =F, spline_args = list(bs = 'ts', k = 3),
                                auto_model_select =T, indicator = indicator, indicator_family = 'binomial', cores = 1){
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
  #print(return_cols)
  #set up with for the return call
  return(setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name))))
}


#df: a data table (post extract covariates)
#covariates: rhs of formula specifying the covariates/columns to be used to help fit the model
#weight_column: column in the data table that specifies the weight
#tc: tree complexity
#lr: learning rate
#bf: bag fraction
#indicator: dependant variable.
#indicator_family: Binomial models assume N as the # of trials and are actually modelled with poission with N as the offset
fit_gbm= function(df, covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, tc = 4, lr = .005, bf = .75, indicator, indicator_family = 'binomial',
                  plot.main = F){


  library(dismo)

  #check to see if its a vector of characters or a psudo-formula
  covariates = format_covariates(add_additional_terms(covariates,additional_terms))

  df = copy(df)

  #format weights
  if(!is.null(weight_column)){
    df[,data_weight := get(weight_column)]
  } else{
    df[,data_weight := 1]
  }
  weight_column = 'data_weight' #specify the column


  # BRT function we use has no binomial. Use emperical logistic or poisson

  #set up poisson outcome structure for binomial and poisson data
  if(indicator_family %in% c('binomial','poisson'))  {
    indicator_family = 'poisson'
    offset=log(df[,N])
    df[,pre_round := get(indicator)]
    #message('WARNING: For Poisson to work, need to round decimals in the response')
    df[,paste0(indicator) := round(pre_round,0)] #round indicator to 0
  } else offset = NULL

  #run the brts. The logic is similiar to the BRT covs (e.g. copy pasted). NOTE, this will run a model for all periods at once. Brt_covs runs each year independantly.
  # learning brt

  message(paste('Fitting GBM/BRT with tc:',tc,'lr:',lr,'bf:',bf))

    
    mod <- try(
      gbm.step(data             = as.data.frame(df),
               gbm.y            = indicator,
               gbm.x            = covariates,
               offset           = offset,
               family           = indicator_family,
               site.weights     = df[,get(weight_column)],
               tree.complexity  = tc,
               learning.rate    = lr,
               bag.fraction     = bf,
               silent           = T,
               plot.main = F,
               plot.folds = F),silent=TRUE)
    
    if(is.null(mod)){
      message('First BRT attempt failed. Lowering Learning Rate by 1/10')
      mod <- try(
        gbm.step(data             = as.data.frame(df),
                 gbm.y            = indicator,
                 gbm.x            = covariates,
                 offset           = offset,
                 family           = indicator_family,
                 site.weights     = df[,get(weight_column)],
                 tree.complexity  = tc,
                 learning.rate    = lr*.1,
                 bag.fraction     = bf,
                 silent           = T,
                 plot.main = F,
                 plot.folds = F))
    }
    if(is.null(mod)){
      message('Second BRT attempt failed. Lowering Original Learning rate by 1/1000 AGAIN')
      mod <- try(
        gbm.step(data             = as.data.frame(df),
                 gbm.y            = indicator,
                 gbm.x            = covariates,
                 offset           = offset,
                 family           = indicator_family,
                 site.weights     = df[,get(weight_column)],
                 tree.complexity  = tc,
                 learning.rate    = lr*.001,
                 bag.fraction     = bf,
                 silent           = T,
                 plot.main = F,
                 plot.folds = F))
    }
    if(is.null(mod)){
      message('Third BRT attempt failed. Slow learn plus low tree complexity')
      mod <- try(
        gbm.step(data             = as.data.frame(df),
                 gbm.y            = indicator,
                 gbm.x            = covariates,
                 offset           = offset,
                 family           = indicator_family,
                 site.weights     = df[,get(weight_column)],
                 tree.complexity  = 2,
                 learning.rate    = lr*.001,
                 bag.fraction     = bf,
                 silent           = T,
                 plot.main = F,
                 plot.folds = F))
    }

    if(is.null(mod)) stop('ALL BRT ATTEMPTS FAILED')

    return(mod)
}

#a function to fit the GBMs for stacking
#basically a wrapper function for the fit_gam function (which is its own wrapper function-- hooray for rabbit holes)
#model_name = what do you want the full fit model to be called upon the return. Must sync with subsquent functions
fit_gbm_child_model = function(df, model_name = 'gbm', fold_id_col = 'fold_id', covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL,
                                tc = 4, lr = 0.005, bf = 0.75, indicator = indicator, indicator_family = indicator_family, cores = 1){

#######The function#######
  library(parallel)

  #prevent df scoping
  df = copy(df)

  #fit the baby trees in parallel
  folds = unique(df[,get(fold_id_col)])

  message('Fitting baby gbm models in parallel')
  baby_models = mclapply(folds, function(fff)
            fit_gbm(df = df[get(fold_id_col) != fff,],
            covariates = covariates,
            additional_terms   = additional_terms,
            weight_column      = weight_column,
            tc                 = tc,
            lr                 = lr,
            bf                 = bf,
            indicator          = indicator,
            indicator_family   = indicator_family,
            plot.main = F), mc.cores = cores)


  for(fff in folds){
    #use the data fit on K-1 of the folds to fit on the help out fold
    df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_models[[fff]], df[get(fold_id_col)==fff,], n.trees=baby_models[[fff]]$gbm.call$best.trees,type='response')]
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
           indicator          = indicator,
           indicator_family   = indicator_family)


  #add a model name slot
  full_model$model_name = model_name

  #predict the main BRT
  df[,paste0(model_name,'_full_pred') := predict(full_model,df,n.trees=full_model$gbm.call$best.trees,type = 'response')]

  suffixes = c('_full_pred', '_cv_pred')
  return_cols = paste0(model_name, suffixes)
  #print(return_cols)
  #set up with for the return call
  return(setNames(list(df[,return_cols,with = F], full_model),c('dataset',paste0(model_name))))

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
#stacker_model: the model object from the stacker
#return children: return the prediction rasters from the child models as well
#todo: allow dataframes of constants to be passed
produce_stack_rasters = function(covariate_layers = all_cov_layers, #raster layers and bricks
                                 period = 1, #period of analysis
                                 child_models = list(), #gam model fitted to full data
                                 stacker_model = stacker,
                                 indicator_family = 'binomial',
                                 return_children = F,
                                 centre_scale_df = NULL){ #stacker model

  message(paste0('The Period is ', period))

  #fetch the covariates appropriate for the period
  period_covs =brick(lapply(covariate_layers, function(x) {
    message(sprintf("on covariate %i: %s", which(names(covariate_layers) %in% names(x)), names(x)))
    fetch_covariate_layer(x,period)}
    ))

  #create constants -- only flexible for year/period
  year = data.frame(year = (1995 +period*5))

  #predict the various models. This is super strict with variable names (and beware scoping issues with the names)
  #brick the results
  stacker_predictors = brick(lapply(child_models, function(x) predict_model_raster(x, period_covs, constants = year, indicator_family = indicator_family, centre_scale_df = centre_scale_df)))

  #take the rasters of the stacker predictors and predict using stacker
  #to: reuse the predict_model_raster function and add glm functionality to it
  stacked_ras = predict_model_raster(model_call = stacker_model, covariate_layers = stacker_predictors, indicator_family = indicator_family)


  #return the results
  if(return_children){
    return(list(setNames(stacked_ras,'stacked_results'), stacker_predictors))
  }else{
    stacked_ras = setNames(stacked_ras,'stacked_results')
    return(stacked_ras)
  }
}


#Wrapper function for produce stacked rasters
#Runs all (or a subset of periods) and formats the results.
make_stack_rasters = function(covariate_layers = all_cov_layers, #raster layers and bricks
                              period = NULL, #period of analysis, NULL returns all periods
                              child_models = list(), #gam model fitted to full data
                              stacker_model = stacked_results[[2]],
                              indicator_family = 'binomial',
                              return_children = F,
                              rd = run_date,
                              re = reg,
                              ind_gp = indicator_group,
                              ind = indicator,
                              ho = holdout,
                              centre_scale_df = NULL){

  ## first, save child_models for use in get.cov.wts in post_estiamtion if not other places too
  save(child_models, file=sprintf('<<<<< FILEPATH REDACTED >>>>>/mbg/%s/%s/output/%s/child_model_list_%s_%i.RData',
                                  ind_gp, ind, rd, re, ho))

  if(is.null(period)){
    period = 1:4
  }

  res = lapply(period, function(the_period) produce_stack_rasters(
              covariate_layers = covariate_layers, #raster layers and bricks. Covariate rasters essentially
              period = the_period, #period of analysis. Fed in from the lapply
              child_models = child_models, #a list of model objects for the child learners
              stacker_model = stacker_model, #the model object of the stacker
              indicator_family= indicator_family,
              return_children = return_children,
              centre_scale_df = centre_scale_df))

  stacked_rasters = brick(lapply(1:length(res), function(x) res[[x]][[1]]))



  if(return_children){
    cms = lapply(1:length(res), function(x) res[[x]][[2]])
    cms = lapply(1:dim(cms[[1]])[[3]], function(x) brick(sapply(cms,'[[', x)))

    ret_obj = unlist(list(stacked_rasters,cms))
  } else{
    ret_obj = list(stacked_rasters)
  }

  #set the names of the return object
  #unique(gsub("(\\.[0-9]+)$","",names(testy[[1]])))
  ret_obj_names = unlist(lapply(ret_obj, function(x) unique(gsub("(\\.[0-9]+)$","",names(x)))))

  #return!
  return(setNames(ret_obj, ret_obj_names))
}

#predict model raster: given a model object and some covariates, predict out the results in raster form
#model call: the model object you want to create rasters from
#covariate_layers: a list of raster-like objects named identically to the models fit on the tabular data
#constants: do any constants need to be fed to the prediction step?
#indicator_family: model family. Specifies what sorts of transformations need doing
#to do: add a transform switch
predict_model_raster = function(model_call, covariate_layers, constants = NULL, indicator_family = 'binomial', centre_scale_df = NULL){
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

    ret_obj = predict(model_call, newdata=dm, n.trees = model_call$gbm.call$best.trees, type = 'response')

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
                       cores = 1){
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
fit_rf_child_model = function(df,model_name = 'rf', fold_id_col = 'fold_id', covariates = all_fixed_effects, additional_terms = NULL, ntree = 1000, indicator = indicator, indicator_family = indicator_family, cores = 1){
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
  #fit all the baby forests at once
  baby_models = mclapply(folds, function(fff) fit_rf(df[get(fold_id_col) != fff,],
                       covariates = covariates,
                       additional_terms = additional_terms,
                       ntree = ntree,
                       indicator = indicator,
                       indicator_family = indicator_family), mc.cores = cores)

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
#cores: how many cores are available
#offset: is there an offset, valid only for poission models. Not fully implemented
fit_glmnet = function(df, covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, alpha = 1, indicator, indicator_family = 'binomial', parallel = F){

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

  #print(head(dm))

  #register parallel
  #registerDoParallel(cores=cores)
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
#cores: how many cores can be used?
#alpha: paramter for glmnet
fit_glmnet_child_model = function(df,model_name = 'glmnet', fold_id_col = 'fold_id', covariates = all_fixed_effects, additional_terms = NULL, weight_column = NULL, alpha = 1, indicator = indicator, indicator_family = 'binomial', cores = 1){
  library(glmnet)
  library(doParallel)
  df = copy(df)
  message('Fitting the Full GLMNET')


  the_covs = format_covariates(add_additional_terms(covariates, additional_terms))

  #start the cluster
  if(cores>1) registerDoParallel(cores = cores)

  #fit the full model
  full_model = fit_glmnet(df,
                     covariates = covariates,
                     additional_terms = additional_terms,
                     weight_column = weight_column,
                     alpha = alpha,
                     indicator = indicator,
                     indicator_family = indicator_family,
                     parallel = cores>1)

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
                       parallel = cores>1)

    new_data = df[get(fold_id_col)==fff,the_covs, with = F]

    n_nd = names(new_data)
    new_data = as.matrix(new_data)
    names(new_data) = n_nd

    df[get(fold_id_col)==fff, paste0(model_name,'_cv_pred') := predict(baby_model, newx = new_data, s = baby_model$cv_1se_lambda, type = 'link')]

  }

  #stop the cluster just in case
  stopImplicitCluster()

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
                                additional_terms = NULL, weight_column = NULL, indicator = indicator, indicator_family = 'binomial', cores = 1){
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

  baby_models = mclapply(folds, function(fff) fit_earth(df = df[get(fold_id_col) != fff,],
                        covariates = covariates,
                        additional_terms = additional_terms,
                        weight_column = weight_column,
                        indicator = indicator,
                        indicator_family = indicator_family), mc.cores = cores)


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
