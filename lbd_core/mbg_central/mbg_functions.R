logit <- function(x) {
  log(x/(1-x))
}
invlogit <- function(x) {
  exp(x)/(1+exp(x))
}

save_mbg_input <- function(indicator = indicator, indicator_group = indicator_group, df = df , simple_raster = simple_raster, simple_raster2 = NULL, mesh_s = mesh_s,
                           mesh_t = mesh_t, cov_list = all_cov_layers, run_date = NULL, pathaddin = "",
                           child_model_names = NULL, all_fixed_effects = NULL, period_map = NULL, centre_scale = T) {

  period_map <- copy(period_map) ## to avoid global scoping issues
  just_covs <- extract_covariates(df, cov_list, return_only_results = T, centre_scale = centre_scale, period_var = 'year', period_map = period_map)
  if(centre_scale==TRUE){
    just_covs <- just_covs$covs
  }
  just_covs <- just_covs[, year := NULL]
  just_covs <- just_covs[, period_id := NULL]
  df = cbind(df, just_covs)

  #create a period column
  if(is.null(period_map)) {
    period_map <- make_period_map(c(2000,2005,2010,2015))
  }
  df[, period := NULL]
  setnames(period_map, 'data_period', 'year')
  setnames(period_map, 'period_id', 'period')
  df = merge(df, period_map, by = 'year', sort =F)

  ## Now that we've extracted covariate values to our data in the buffer zone, clip cov list to simple_raster instead of simple_polygon
  ##    (simple_raster is area we actually want to model over)
  for(l in 1:length(cov_list)) {
    cov_list[[l]]  <- crop(cov_list[[l]], extent(simple_raster))
    cov_list[[l]]  <- setExtent(cov_list[[l]], simple_raster)
    cov_list[[l]]  <- raster::mask(cov_list[[l]], simple_raster)
  }

  ## Save all inputs
  to_save <- c("df", "simple_raster", "mesh_s", "mesh_t", 'cov_list', 'child_model_names', 'all_fixed_effects','period_map')
  if(!is.null(simple_raster2)) to_save <- c(to_save, "simple_raster2")
  save(list = to_save, file = paste0("<<<< FILEPATH REDACTED >>>>", pathaddin, '.RData'))
  message(paste0('All inputs saved to "<<<< FILEPATH REDACTED >>>>', pathaddin,'.RData'))

}

transform_theta <- function(theta) {(exp(theta)-1) / (1+exp(theta))}

test_rho_priors <- function(temporal_model_theta1_prior) {
  temporal_model_theta1_prior <- eval(parse(text = temporal_model_theta1_prior))
  mean <- temporal_model_theta1_prior$param[1]
  prec <- temporal_model_theta1_prior$param[2]
  sd <- sqrt(1/prec)

  message(paste0('theta1_prior_prec: ', round(prec, 2)))
  lower <- transform_theta(mean - (1.96*sd))
  upper <- transform_theta(mean + (1.96*sd))
  message(paste0('95% range in Rho prior: ', round(transform_theta(mean), 2), ' (', round(lower, 2), ' - ', round(upper, 2), ')'))
}

build_mbg_formula_with_priors <- function(fixed_effects,
                                          positive_constrained_variables = NULL,
                                          interact_with_year = NULL,
                                          int = TRUE,
                                          nullmodel = FALSE,
                                          add_nugget = FALSE,
                                          nugget_prior = "list(prior = 'loggamma', param = c(2, 1))",
                                          add_ctry_res = FALSE,
                                          ctry_re_prior = "list(prior = 'loggamma', param = c(2, 1))",
                                          temporal_model_type="'ar1'",
                                          temporal_model_theta_prior = "list(prior = 'loggamma', param = c(1, 0.00005))",
                                          temporal_model_theta1_prior = "list(prior = 'normal', param = c(0, 1/(2.58^2)))",
                                          no_gp = FALSE,
                                          stacker_names = child_model_names,
                                          coefs.sum1 = FALSE,
                                          subnat_RE = FALSE,
                                          subnat_re_prior = "list(prior = 'loggamma', param = c(1, 5e-5))") {

  if(nchar(stacker_names[1]) == 0 & coefs.sum1 == TRUE){
    message("WARNING! You've chosen sum-to-1 but don't appear to be using any stackers. Unless you have a very good reason to do this, it probably doesn't make sense. As such, we're setting coefs.sum1 <- FALSE")
    coefs.sum1 <- FALSE
  }

  # Set up model equation
  if(int==TRUE) intercept='+int' else intercept =''
  f_null  <- formula(paste0('covered~-1',intercept))

  if(!is.null(interact_with_year)){
    for(iwy in interact_with_year){
      fixed_effects <- gsub(iwy,paste0(iwy,'* factor(period)'),fixed_effects)
    }
    if(!is.null(positive_constrained_variables))
      stop("Cannot Both Constrain and Interact, one must be NULL.")
  }


  if(!is.null(positive_constrained_variables)){
    if(!all(positive_constrained_variables %in% strsplit(fixed_effects,' \\+ ')[[1]]))
      stop('Not all your positive_constrained_variables match fixed_effects')
    for(pcv in positive_constrained_variables){
      v <- sprintf("f(%s, model='clinear',range=c(0,Inf),initial=0)",pcv)
      fixed_effects <- gsub(pcv,v,fixed_effects)
    }
  }

  f_nugget <- as.formula(paste0('~f(IID.ID, model = "iid", hyper = list(theta=', nugget_prior, '))'))
  f_res <- as.formula(paste0('~f(CTRY.ID, model = "iid", hyper = list(theta=', ctry_re_prior, '))'))
  f_subnat <- paste0('f(SUBNAT.ID, model = "iid", constr=TRUE, hyper = list(theta=', subnat_re_prior, "))")
  
  test_rho_priors(temporal_model_theta1_prior) ## Report how priors for theta1 (Rho) are being used in Rho space.
  f_space <- as.formula(paste0('~f(space,
                                   model = spde,
                                   group = space.group,
                                   control.group = list(model = ', temporal_model_type, ", ",
                                                       "hyper = list(theta = ", temporal_model_theta_prior, ", theta1 = ", temporal_model_theta1_prior, ")))"))

  f.e.v <- base::strsplit(all_fixed_effects, ' + ', fixed=T)[[1]] ## fixed effect vector
  f_sum1 <- as.formula(paste0('~f(covar,
                                  model = \'iid\',
                                  extraconstr = list(A = matrix(1, 1, ', length(f.e.v), '), e = 1),
                                    hyper=list(theta=list(initial=log(inla.set.control.fixed.default()$prec),
                                               fixed=TRUE)))'))

  if(nchar(fixed_effects) <= 1){ 
    nullmodel <- TRUE
  }

  ## build formula starting with f_null

  f_mbg <- f_null

  if(!nullmodel){
    if(coefs.sum1){
      f_lin <- f_sum1 
    } else{
      f_lin <- reformulate(fixed_effects)
    }

    f_mbg <- f_mbg + f_lin
  }

  if(!no_gp) f_mbg <- f_mbg + f_space
  if(add_nugget==TRUE) f_mbg <- f_mbg + f_nugget
  if(add_ctry_res == TRUE) f_mbg <- f_mbg + f_res
  if (subnat_RE == TRUE) f_mbg <- paste(f_mbg, f_subnat, sep = "+")
  
  message(f_mbg)
  return(f_mbg)
}

build_mbg_formula <- function(fixed_effects,
                              positive_constrained_variables=NULL,
                              interact_with_year=NULL,
                              int=TRUE,
                              nullmodel=FALSE,
                              add_nugget=FALSE,
                              add_ctry_res = FALSE) {

  # Set up model equation
  if(int==TRUE) intercept='+int' else intercept =''
  f_null  <- formula(paste0('covered~-1',intercept))

  if(!is.null(interact_with_year)){
    for(iwy in interact_with_year){
      fixed_effects <- gsub(iwy,paste0(iwy,'* factor(period)'),fixed_effects)
    }
    if(!is.null(positive_constrained_variables))
      stop("Cannot Both Constrain and Interact, one must be NULL.")
  }


  if(!is.null(positive_constrained_variables)){
    if(!all(positive_constrained_variables %in% strsplit(fixed_effects,' \\+ ')[[1]]))
      stop('Not all your positive_constrained_variables match fixed_effects')
    for(pcv in positive_constrained_variables){
      v <- sprintf("f(%s, model='clinear',range=c(0,Inf),initial=0)",pcv)
      fixed_effects <- gsub(pcv,v,fixed_effects)
    }
  }

  f_nugget <- ~ f(IID.ID, model='iid',
               hyper=list(theta=list(prior="loggamma",param=c(2, 1))))

  f_res <- ~ f(CTRY.ID, model = 'iid',
               hyper=list(theta=list(prior="loggamma",param=c(2, 1))))

  f_space <- ~ f(
    space,
    model = spde,
    group = space.group,
    control.group = list(model = 'ar1'))

  if(!nullmodel){
    f_lin <- reformulate(fixed_effects)
    f_mbg = f_null + f_lin + f_space
    if(add_nugget==TRUE) f_mbg = f_null + f_lin + f_space + f_nugget
    if(add_ctry_res == TRUE) f_mbg = f_null + f_lin + f_space + f_res
    if(add_nugget == TRUE & add_ctry_res == TRUE) f_mbg = f_null + f_lin + f_space + f_nugget + f_res
  } else {
    f_mbg = f_null + f_space
  }
  message(f_mbg)
  return(f_mbg)
}


# spde priors function
local.inla.spde2.matern.new = function(mesh, alpha=2, prior.pc.rho,
                                       prior.pc.sig)
{
  # Call inla.spde2.matern with range and standard deviationparametrization
  d = INLA:::inla.ifelse(inherits(mesh, "inla.mesh"), 2, 1)
  nu = alpha-d/2
  kappa0 = log(8*nu)/2
  tau0   = 0.5*(lgamma(nu)-lgamma(nu+d/2)-d/2*log(4*pi))-nu*kappa0
  spde   = inla.spde2.matern(mesh = mesh,
                             B.tau   = cbind(tau0,   nu,  -1),
                             B.kappa = cbind(kappa0, -1, 0))

  # Change prior information
  param = c(prior.pc.rho, prior.pc.sig)
  spde$f$hyper.default$theta1$prior = "pcspdega"
  spde$f$hyper.default$theta1$param = param
  spde$f$hyper.default$theta1$initial = log(prior.pc.rho[1])+1
  spde$f$hyper.default$theta2$initial = log(prior.pc.sig[1])-1

  # End and return
  return(invisible(spde))
}

build_mbg_data_stack <- function(df, fixed_effects, mesh_s, mesh_t,
                                 exclude_cs='',usematernnew=F,
                                 sig0=0.5,rho0=0.3,
                                 coefs.sum1 = FALSE,
                                 use_ctry_res = FALSE,
                                 use_subnat_res = FALSE,
                                 remove_non_subnats = FALSE,
                                 use_nugget = FALSE,
                                 stacker_names = child_model_names,
                                 yl   = year_list,
                                 zl   = z_list,
                                 zcol = zcol,
                                 scale_gaussian_variance_N = TRUE,
                                 tmb  = FALSE,
                                 shapefile_version = 'current') {

  if(nchar(stacker_names[1]) == 0 & coefs.sum1 == TRUE){
    message("WARNING! You've chosen sum-to-1 but don't appear to be using any stackers. Unless you have a very good reason to do this, it probably doesn't make sense. As such, we're setting coefs.sum1 <- FALSE")
    coefs.sum1 <- FALSE
  }

  if(tmb == TRUE){
    message('Returning a TMB model stack.')
    return(
      build_mbg_data_stack_tmb(d          = df,
                               yl         = yl,                 
                               fes        = fixed_effects,  
                               indic      = indicator,      
                               exclude_cs = exclude_cs,
                               nugget     = use_nugget,
                               country_re = use_ctry_res,
                               nid_re     = use_nid_res,
                               zl         = zl,
                               zcol       = zcol,
                               shapefile_version         = shapefile_version, 
                               scale_gaussian_variance_N = scale_gaussian_variance_N,
                               mesh       = mesh_s))            
  } else {
    # construct an SPDE model with a Matern kernel
    message('Building SPDE...')
    if(usematernnew){
      # rho0 is typical range, sig0 typical sd
      spde = local.inla.spde2.matern.new(mesh         = mesh_s,
                                         prior.pc.rho = c(rho0, 0.5),
                                         prior.pc.sig = c(sig0, 0.5),
                                         alpha        = 2)
    } else {
      spde <- inla.spde2.matern(mesh = mesh_s, alpha = 2)
    }

    ## Build projector matrix between data locs and spatial mesh
    data.locs <- as.matrix(df[, c('longitude', 'latitude'), with=F])
    if(mesh_s$manifold == "S2"){
      data.locs <- lonlat3D(data.locs[, 1], data.locs[, 2])
    }
    A <- inla.spde.make.A(
      mesh = mesh_s,
      loc = data.locs,
      group = df$period,
      group.mesh = mesh_t
    )
    
    if(coefs.sum1){
      f.e.v <- stacker_names ## fixed eff. vec.
      A.covar <- as.matrix(df[, f.e.v, with = FALSE])
    }

    if (is.null(mesh_t)){
      space = inla.spde.make.index("space",
                                   n.spde = spde$n.spde)
    } else {
    space = inla.spde.make.index("space",
                                 n.spde = spde$n.spde,
                                 n.group = mesh_t$m)
    }


    # find cov indices
    if(fixed_effects!="NONE" & nchar(fixed_effects) > 0){
      f_lin <- reformulate(fixed_effects)
      message('Indexing covariates...')
      covs_indices <- unique(c(match(all.vars(f_lin), colnames(df))))

      # make design matrix, center the scaling
      design_matrix <- data.frame(int = 1,
                                  df[,covs_indices,with=F])

      cs_df <- getCentreScale(design_matrix, exclude = c('int',exclude_cs))

      design_matrix <- centreScale(design_matrix,
                                   df = cs_df)
    } else{
      design_matrix <- data.frame(int = rep(1,nrow(df)))
      cs_df <- getCentreScale(design_matrix, exclude = c('int','rates'))
    }

    # construct a 'stack' object for observed data
    cov=df[[indicator]] # N+_i
    N=df$N                 # N_i

    if(use_ctry_res){
      ## add an numeric gaul code to use in random effects
      design_matrix$CTRY.ID <- gaul_convert(df$country, shapefile_version = shapefile_version)
    }
    
    if (use_subnat_res) {
      ## add subnat ID to use in random effects
      design_matrix$SUBNAT.ID <- df$ADM1_CODE
      
      if(remove_non_subnats) {
        ## NA out subnational rows which are 0
        design_matrix$SUBNAT.ID[design_matrix$SUBNAT.ID == 0] <- NA
      }
    }

    if(use_nugget){
      design_matrix$IID.ID <- 1:nrow(design_matrix)
    }

    message('Stacking data...')
    if(coefs.sum1 == TRUE & nchar(fixed_effects) > 1){
      ## build in covar effect to be used in sum1 constraint
      stack.obs <- inla.stack(
        data = list(covered = cov),
        A = list(A, 1, A.covar),
        effects = list(space,
                       design_matrix,
                       list(covar = 1:ncol(A.covar))),
        tag = 'est'
      )
    }else{
      stack.obs <- inla.stack(
        data = list(covered = cov),
        A = list(A, 1),
        effects = list(space,
                       design_matrix),
        tag = 'est'
      )
    }

    return_list <- list(stack.obs, spde, cs_df)

    return(return_list)
  }
}


fit_mbg <- function(indicator_family, stack.obs, spde, cov, N, int_prior_mn, int_prior_prec=1, f_mbg, run_date, keep_inla_files, cores, 
                    verbose_output = FALSE, wgts = 0, intstrat = 'eb', fe_mean_prior = 0, fe_sd_prior = 1) {

  if(wgts[1] == 0) wgts = rep(1, length(N))

  inla.setOption("enable.inla.argument.weights", TRUE)

  message('Fitting INLA model')

  inla_working_dir <- "<<<< FILEPATH REDACTED >>>>"
  dir.create(inla_working_dir, showWarnings = FALSE)
  if(indicator_family=='binomial') {
    system.time(
      res_fit <- inla(f_mbg,
                      data = inla.stack.data(stack.obs),
                      control.predictor = list(A = inla.stack.A(stack.obs),
                                               link = 1,
                                               compute = FALSE),
                      control.fixed = list(expand.factor.strategy = 'inla',
                                           prec.intercept = int_prior_prec,
                                           mean.intercept = int_prior_mn,
                                           mean = fe_mean_prior,
                                           prec = fe_sd_prior),
                      control.compute = list(dic = TRUE,
                                             cpo = TRUE,
                                             config = TRUE),
                      control.inla = list(int.strategy = intstrat, h = 1e-3, tolerance = 1e-6),
                      family = 'binomial',
                      num.threads = cores,
                      Ntrials = N,
                      verbose = verbose_output,
                      working.directory = inla_working_dir,
                      weights = wgts,
                      keep = TRUE)
    )
  }
  if(indicator_family=='gaussian') {
    system.time(
      res_fit <- inla(f_mbg,
                      data = inla.stack.data(stack.obs),
                      control.predictor = list(A = inla.stack.A(stack.obs),
                                               compute = FALSE),
                      control.fixed = list(expand.factor.strategy = 'inla',
                                           prec.intercept = 1,
                                           mean.intercept = int_prior_mn,
                                           mean = 0,
                                           prec = 2),
                      control.compute = list(dic = TRUE,
                                             cpo = TRUE,
                                             config = TRUE),
                      control.inla = list(int.strategy = intstrat, h = 1e-3, tolerance = 1e-6),
                      family = 'Gaussian',
                      num.threads = cores,
                      verbose = TRUE,
                      working.directory = inla_working_dir,
                      weights = wgts,
                      scale = N,
                      keep = TRUE)
    )
  }
  if(indicator_family=='beta') {
    system.time(
      res_fit <- inla(f_mbg,
                      data = inla.stack.data(stack.obs),
                      control.predictor = list(A = inla.stack.A(stack.obs),
                                               compute = FALSE),
                      control.fixed = list(expand.factor.strategy = 'inla',
                                           prec.intercept = 1,
                                           mean.intercept = int_prior_mn,
                                           mean = 0,
                                           prec = 2),
                      control.compute = list(dic = TRUE,
                                             cpo = TRUE,
                                             config = TRUE),
                      control.inla = list(int.strategy = intstrat, h = 1e-3, tolerance = 1e-6),
                      family = 'Beta',
                      num.threads = cores,
                      verbose = TRUE,
                      working.directory = inla_working_dir,
                      weights = wgts,
                      scale = N,
                      keep = TRUE)
    )
  }

  return(res_fit)

}

predict_mbg <- function(res_fit, cs_df, mesh_s, mesh_t, cov_list,
                        samples, simple_raster, transform,
                        ind=indicator, indg=indicator_group,
                        rd=run_date, region = reg, full_df = df,
                        coefs.sum1 = FALSE, pred_gp = TRUE,
                        fixed_effects = all_fixed_effects,
                        nperiod = ifelse(is.character(year_list),
                                         length(eval(parse(text = year_list))),
                                         length(year_list)),
                        stacker_names = child_model_names,
                        shapefile_version = 'current',
                        simple_raster_subnats = NULL,
                        subnat_country_to_get = NULL) {

  if(nchar(stacker_names[1]) == 0 & coefs.sum1 == TRUE){
    message("WARNING! You've chosen sum-to-1 but don't appear to be using any stackers. Unless you have a very good reason to do this, it probably doesn't make sense. As such, we're setting coefs.sum1 <- FALSE")
    coefs.sum1 <- FALSE
  }

  for(l in 1:length(cov_list)) {
    message(sprintf("On cov %i out of %i", l, length(cov_list)))
    cov_list[[l]]  <- crop(cov_list[[l]], extent(simple_raster))
    cov_list[[l]]  <- setExtent(cov_list[[l]], simple_raster)
    cov_list[[l]]  <- mask(cov_list[[l]], simple_raster)
  }

  message('Making predictions')

  # number of samples
  n_draws <- samples

  # dummy raster
  template <- simple_raster

  cell_idx <- seegSDM:::notMissingIdx(template)

  # sample from posterior over latents
  suppressWarnings(draws <- inla.posterior.sample(n_draws, res_fit))

  # Save full inla draws object above to get hyperpar draws later (tau for Gaussian, etc.)
  run_dir <- paste0('<<<< FILEPATH REDACTED >>>>', '/inla_draws')
  dir.create(run_dir, showWarnings = FALSE)
  save(list='draws', file = paste0(run_dir, '/inla_draws_', region, '.RData'))

  # get parameter names
  par_names <- rownames(draws[[1]]$latent)

  # index to spatial field and linear coefficient samples
  s_idx <- grep('^space.*', par_names) 
  if(coefs.sum1){
    l_idx <- grep('covar', par_names)
    l_idx_int <- match(sprintf('%s.1', res_fit$names.fixed), ## main effects
                       par_names)
    if(mean(is.na(l_idx_int)) == 1) l_idx_int <- match(sprintf('%s', res_fit$names.fixed), par_names)
    if(mean(is.na(l_idx_int)) == 1) l_idx_int <- match(sprintf('%s:1', res_fit$names.fixed), par_names)##
    l_idx = c(l_idx, 
              l_idx_int) 

    if(length(l_idx) < 1){ 
      stop('Looks like you may have requested sum.to.1 in predict_mbg on a model that was not run using sum.to.1 constraint.')
    }

  }else{
    l_idx <- match(sprintf('%s.1', res_fit$names.fixed), 
                   par_names)
    if(mean(is.na(l_idx)) == 1) l_idx <- match(sprintf('%s', res_fit$names.fixed), par_names) 

    if(length(l_idx) < 1){ ## 1 b/c of the intercept
      stop('Looks like you may have set coefs.sum.1 to FALSE in predict_mbg using a fitted model that WAS run using sum.to.1 constraint')
    }
  }

  if(pred_gp){
    pred_s <- sapply(draws, function (x)
      x$latent[s_idx])
  }
  pred_l <- sapply(draws, function (x)
    x$latent[l_idx])
  if(length(l_idx)==1) pred_l=t(as.matrix(pred_l))

  if(coefs.sum1){
    rownames(pred_l) <- c(stacker_names, 
                          res_fit$names.fixed) 
  } else {
    rownames(pred_l) <- res_fit$names.fixed
  }


  ## if we fit with a nugget, we also need to take draws of the nugget precision
  if(length(grep('^IID.ID.*', par_names)) > 0){
    pred_n <- sapply(draws, function(x) {
      nug.idx <- which(grepl('IID.ID', names(draws[[1]]$hyper)))
      x$hyperpar[[nug.idx]]})
  }else{
    pred_n <- NULL
  }

  if(sum(grepl('CTRY.ID', par_names)) > 0){

    ctry.res.idx <- grep('CTRY.ID*', par_names)
    pred_ctry_res  <-  sapply(draws, function (x) x$latent[ctry.res.idx])
    ctry.res.names <- par_names[ctry.res.idx]
    if (length(ctry.res.names) == 1) {
      message('WARNING:ONLY ONE COUNTRY PRESENT IN DATA AND RE FIT')
      print('WARNING:ONLY ONE COUNTRY PRESENT IN DATA AND RE FIT')
      pred_ctry_res <- matrix(pred_ctry_res, nrow = 1)
    }
    rownames(pred_ctry_res) <- res_fit$summary.random$CTRY.ID$ID

    pred_ctry_prec <- sapply(draws, function(x){
      ctry.prec.idx <- which(grepl('CTRY.ID', names(draws[[1]]$hyper)))
      x$hyperpar[[ctry.prec.idx]]}) ## this gets the precision for country REs
  }else{
    pred_ctry_res <- NULL
  }

  if (sum(grepl("SUBNAT.ID", par_names)) > 0) {
    
    subnat.res.idx <- grep("SUBNAT.ID*", par_names)
    pred_subnat_res <- sapply(draws, function(x) x$latent[subnat.res.idx])
    subnat.res.names <- par_names[subnat.res.idx]
    if (length(subnat.res.names) == 1) {
      message("WARNING:ONLY ONE SUBNAT UNIT PRESENT IN DATA AND RE FIT")
      print("WARNING:ONLY ONE SUBNAT UNIT PRESENT IN DATA AND RE FIT")
      pred_subnat_res <- matrix(pred_subnat_res, nrow = 1)
    }
    rownames(pred_subnat_res) <- res_fit$summary.random$SUBNAT.ID$ID
    
    pred_subnat_prec <- sapply(draws, function(x) {
      subnat.prec.idx <- which(grepl("SUBNAT.ID", names(draws[[1]]$hyper)))
      x$hyperpar[[subnat.prec.idx]]
    }) ## this gets the precision for country REs
  } else {
    pred_subnat_res <- NULL
  }
  
  ## get coordinates of cells to predict to
  ## these are are needed both for raster extraction and GP projection
  coords <- xyFromCell(template, seegSDM:::notMissingIdx(template))

  ## Create coords and stuff for the subnational simple raster
  if(!is.null(simple_raster_subnats)) {
    template_subnat <- simple_raster_subnats
    subnat_cell_idx <- seegSDM:::notMissingIdx(template_subnat)
    subnat_coords <- raster::xyFromCell(template_subnat, notMissingIdx(template_subnat))
  }
  
  ## setup coords for GP projection. convert coords to spherical if
  ## using spherical modeling mesh. used if you made a mesh on s2
  if(mesh_s$manifold == "S2"){

    gp_coords <- lonlat3D(coords[, 1], coords[, 2])
    
  } else {
    gp_coords <- coords
  }

  # ~~~~~~~~~~~~~~
  # project spatial effect
  # make predictions for all periods
  # get predictor matrix between spatial nodes and prediction locations
  ## replicate coordinates and years for GP projection

  gp_coords_periods <- do.call(rbind,
                               replicate(nperiod,
                                         gp_coords,
                                         simplify = FALSE))

  groups_periods <- rep(1:nperiod,
                        each = nrow(gp_coords))

  # Projector matrix
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = gp_coords_periods,
    group = groups_periods,
    group.mesh = mesh_t
  )

  ## get samples of s for all cells
  if(pred_gp){
    cell_s <- as.matrix(Matrix::crossprod(Matrix::t(A.pred), pred_s))
  }else{
    cell_s <- 0
  }

  tvnames = pars =c()
  for(i in 1:length(cov_list)){
    if(dim(cov_list[[i]])[3]!=1) tvnames[length(tvnames)+1]= names(cov_list)[i]
    if(dim(cov_list[[i]])[3]==1) pars[length(pars)+1]= names(cov_list)[i]
  }

  cov_list = setNames(lapply(cov_list, function(x) mask(resample(x,template),template)),names(cov_list))

  # extract the covariates
  vals = data.table(do.call(cbind, lapply(cov_list, function(x) raster::extract(x,coords))))

  # create the int column
  vals[,int:= 1]

  # reshape long
  vals[,id:= 1:nrow(vals)] #create an id to ensure that melt doesn't sort things

  if (!is.null(mesh_t)){
    tv_cov_colnames = grep(paste(tvnames, collapse = "|"), names(vals), value = T) #unlisted
    tv_cov_colist = lapply(tvnames, function(x) grep(paste0("^", x, "\\.[[:digit:]]*$"), names(vals), value = T))

    vals = melt(vals, id.vars = c('id','int',pars), measure = tv_cov_colist, value.name = tvnames, variable.factor =F)

    vals[,variable:= as.numeric(substring(variable,regexpr("(\\.[0-9]+)$", variable)[1]+1))]
  }

  vals = vals[,mget(rownames(pred_l))]

  vals = centreScale(vals,cs_df) #happens for all periods at once

  cell_l <- unname(as.matrix(as(data.matrix(vals), "dgeMatrix") %*% pred_l))

  if(!is.null(pred_n)){
    cell_n <- sapply(pred_n, function(x){
      rnorm(n = nrow(cell_l), sd = 1 / sqrt(x), mean = 0)
    })
  }

  if(!is.null(pred_ctry_res)){

    cell_ctry_res <- matrix(0, ncol = ncol(cell_l), nrow = nrow(cell_l))

    reg.gauls <- get_adm0_codes(reg, shapefile_version = shapefile_version)
    gauls.with.res <- as.numeric(rownames(pred_ctry_res))
    gaul.vec <- values(simple_raster)[cell_idx]
    gaul.t.vec <- rep(gaul.vec, nperiod) ## expand gaul vec in time
    for(gg in reg.gauls){
      gg.rows <- which(gaul.t.vec == gg)

      if(gg %in% gauls.with.res){
        message(paste0('adding RE for country ', gg))
        gg.idx  <- which(gauls.with.res == gg)
        cell_ctry_res[gg.rows, ] <- cell_ctry_res[gg.rows, ] +
          matrix(rep(pred_ctry_res[gg.idx, ], length(gg.rows)), ncol = ncol(cell_ctry_res), byrow = TRUE)
      }else{
        re.draws <- rnorm(n = ncol(cell_ctry_res), mean = 0, sd = 1 / sqrt(pred_ctry_prec))
        cell_ctry_res[gg.rows, ] <- cell_ctry_res[gg.rows, ] +
        matrix(rep(re.draws, length(gg.rows)), ncol = ncol(cell_ctry_res), byrow = TRUE)
      }

    }

  }
  
  cell_l_sd <-  apply(cell_l, 1, sd)

  ## project model uncertainty from node-level sd to pred grid level
  if(pred_gp){
    node_s_sd <- apply(pred_s, 1, sd)
    cell_s_sd <- as.matrix(A.pred %*% node_s_sd)
  }else{
    cell_s_sd <- 0
  }

  ## get uncertainty from nugget
  if(!is.null(pred_n)){
    cell_n_sd <- apply(cell_n, 1, sd)
  }

  ## get uncertainty from country res
  if(!is.null(pred_ctry_res)){
    cell_ctry_res_sd <- apply(cell_ctry_res, 1, sd)
  }

  ## get uncertainty from subnat res
  if (!is.null(pred_subnat_res)) {
    cell_subnat_res_sd <- apply(cell_subnat_res, 1, sd)
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## combine to get cell-level sd
  cell_sd <- cell_l_sd^2 + cell_s_sd^2
  
  # add nugget if present
  if (!is.null(pred_n)) {
    cell_sd <- cell_sd + cell_n_sd^2
  }
  
  # add country res if present
  if (!is.null(pred_ctry_res)) {
    cell_sd <- cell_sd + cell_ctry_res_sd^2
  }
  
  # add subnat res if present
  if (!is.null(pred_subnat_res)) {
    cell_sd <- cell_sd + cell_subnat_res_sd^2
  }
  
  cell_sd <- sqrt(cell_sd)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # combine to produce draws 
  cell_all <- cell_l + cell_s

  # add nugget if present
  if (!is.null(pred_n)) {
    cell_all <- cell_all + cell_n
  }

  # add country res if present
  if (!is.null(pred_ctry_res)) {
    cell_all <- cell_all + cell_ctry_res
  }
  
  # add subnat res if present
  if (!is.null(pred_subnat_res)) {
    cell_all <- cell_all + cell_subnat_res
  }

  # get predictive draws on probability scale
  if(transform=='inverse-logit') { cell_pred <- plogis(as.matrix(cell_all))
  } else cell_pred <- eval(parse(text=transform))

  # get prediction mean (integrated probability)
  pred_mean <- rowMeans(cell_pred)

  # create multi-band rasters for each of these metrics
  # each band giving the values for a different period
  mean_ras <- insertRaster(template,
                           matrix(pred_mean,
                                  ncol = nperiod))

  sd_ras <- insertRaster(template,
                         matrix(cell_sd,
                                ncol = nperiod))

  names(mean_ras) <-
    names(sd_ras) <-
    paste0('period_', 1:nperiod)

  return_list <- list(mean_ras, sd_ras, cell_pred)

  return(return_list)

}

save_mbg_preds <- function(config, time_stamp, run_date, mean_ras, sd_ras, res_fit, cell_pred, df, pathaddin = "") {

  if(time_stamp==TRUE) output_dir <- "<<<< FILEPATH REDACTED >>>>"
  if(time_stamp==FALSE) output_dir <- "<<<< FILEPATH REDACTED >>>>"
  dir.create(output_dir, showWarnings = FALSE)

  write.csv(config, paste0(output_dir, '/config.csv'), row.names = FALSE)

  if (!is.null(mean_ras)) {
    writeRaster(
      mean_ras,
      file = (paste0(output_dir, '/', indicator, '_prediction_eb', pathaddin)),
      overwrite = TRUE
    )
  }

  if (!is.null(sd_ras)) {

    writeRaster(
      sd_ras,
      file = (paste0(output_dir, '/', indicator, '_sd_eb', pathaddin)),
      overwrite = TRUE
    )
  }

  save(res_fit,
       file = (paste0(output_dir, '/', indicator, '_model_eb', pathaddin, '.RData')))
  
  save(
    cell_pred,
    file = (paste0(output_dir, '/', indicator,'_cell_draws_eb',pathaddin,'.RData')),
    compress = TRUE
  )

  write.csv(
    df,
    file = (paste0(output_dir, '/', indicator, '_trainingdata', pathaddin)),
    row.names = FALSE
  )

  write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

}
