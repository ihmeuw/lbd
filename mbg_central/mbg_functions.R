## Logit functions
logit <- function(x) {
  log(x/(1-x))
}
invlogit <- function(x) {
  exp(x)/(1+exp(x))
}

## Save mbg inputs
save_mbg_input <- function(indicator = indicator, indicator_group = indicator_group, df = df , simple_raster = simple_raster, simple_raster2 = NULL, mesh_s = mesh_s,
                           mesh_t = mesh_t, cov_list = all_cov_layers, run_date = NULL, pathaddin="",
                           child_model_names = NULL, all_fixed_effects = NULL, period_map = NULL,centre_scale = T) {

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
  save(list = to_save, file = paste0('<<<< FILEPATH REDACTED >>>>/', run_date, pathaddin, '.RData'))
  message(paste0('All inputs saved to /<<<< FILEPATH REDACTED >>>>/', run_date, pathaddin,'.RData'))

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

#' @title Build the spde and spde prior
#' 
#' @description construct the spde and relevant priors for Matern GP so that it 
#'   can be used with INLA or TMB
#'   
#' @param spde_prior String containing a list. Specifies the type of prior used
#'   for the Matern model parameters. If type = "nonpc" then use non pc 
#'   prior. User can specify nominal prior means: spde_prior$prior$range.nominal 
#'   and spde_prior$prior$variance.nominal. Defaults are INLA defaults, that is 
#'   20\% of the range of the mesh and 1, respectively. If type = "pc" then use 
#'   pc prior.  User can specify spde_prior$prior$range = (range0, Prange) i.e. 
#'   P(range < range0) = Prange and spde_prior$prior$sigma = (sigma0, Psigma) 
#'   i.e. P(sigma > sigma0) = Psigma. Defaults for Prange and Psigma are 0.05. 
#'   Default for range0 is 5\% of the range of the mesh and sigma0=3.
#'   
#' @param mesh_s An 'inla.mesh' object used to fit the spatial SPDE GP
#'   approximation
#'   
#' @param st_gp_int_zero Logical. Should the space-time GP be forced to interate to 0?
#' 
#' @return List containing the spde and the prior filled in with defaults if no
#'   other values are specified
build_spde_prior <- function(spde_prior, mesh_s, st_gp_int_zero) {
  if(spde_prior$type=="pc"){ # PC prior
    if(is.null(spde_prior$prior$sigma)) {
      spde_prior$prior$sigma <- c(3, 0.05) # P(sigma > 3) = 0.05
    }
    if(is.null(spde_prior$prior$range)) {
      mesh.range <- max(c(diff(range(mesh_s$loc[, 1])), 
                          diff(range(mesh_s$loc[, 2])), 
                          diff(range(mesh_s$loc[, 3]))))
      spde_prior$prior$range <- c(mesh.range*0.05, 0.05) # P(range < 5% max extent of mesh) = 0.05
    }
    message(paste("Building spde with pc prior,",
                  spde_prior$prior$range[2]*100,
                  "% probability that the range is lower than",
                  spde_prior$prior$range[1], 
                  "and a",
                  spde_prior$prior$sigma[2]*100,
                  "% probability that sigma is greater than",
                  spde_prior$prior$sigma[1]))
    spde <- inla.spde2.pcmatern(mesh = mesh_s, 
                                alpha = 2,
                                prior.range = spde_prior$prior$range,
                                prior.sigma = spde_prior$prior$sigma,
                                constr = st_gp_int_zero)
  } else { # Non PC prior
    if(is.null(spde_prior$prior$variance.nominal)) {
      spde_prior$prior$variance.nominal <- 1
    }
    spde <- inla.spde2.matern(mesh = mesh_s,  alpha = 2, constr = st_gp_int_zero,
                              prior.range.nominal = spde_prior$prior$range.nominal,
                              prior.variance.nominal = spde_prior$prior$variance.nominal)
  }
  return(list(spde=spde, spde_prior=spde_prior))
}

#' @title Build an INLA formula for mbg
#' @description construct a formula object for use with R-INLA::inla() model fitting in conjunction with the prepped data object from build_mbg_data_stack() 
#' @param fixed_effects PARAM_DESCRIPTION
#' @param positive_constrained_variables PARAM_DESCRIPTION, Default: NULL
#' @param interact_with_year PARAM_DESCRIPTION, Default: NULL
#' @param int PARAM_DESCRIPTION, Default: TRUE
#' @param nullmodel PARAM_DESCRIPTION, Default: FALSE
#' @param add_nugget PARAM_DESCRIPTION, Default: FALSE
#' @param nugget_prior PARAM_DESCRIPTION, Default: 'list(prior = 'loggamma', param = c(2, 1))'
#' @param add_ctry_res PARAM_DESCRIPTION, Default: FALSE
#' @param ctry_re_prior PARAM_DESCRIPTION, Default: 'list(prior = 'loggamma', param = c(2, 1))'
#' @param temporal_model_type PARAM_DESCRIPTION, Default: ''ar1''
#' @param temporal_model_theta_prior PARAM_DESCRIPTION, Default: 'list(prior = 'loggamma', param = c(1, 0.00005))'
#' @param temporal_model_theta1_prior PARAM_DESCRIPTION, Default: 'list(prior = 'normal', param = c(0, 1/(2.58^2)))'
#' @param no_gp PARAM_DESCRIPTION, Default: FALSE
#' @param stacker_names PARAM_DESCRIPTION, Default: child_model_names
#' @param coefs.sum1 PARAM_DESCRIPTION, Default: FALSE
#' @param subnat_RE Additional admin-1 random effect. Default: FALSE
#' @param use_space_only_gp Logical. include a space only (time stationary) gp. Default: FALSE
#' @param use_time_only_gmrf Logical. include a time only gp. Default: FALSE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @seealso
#'  \code{\link[base]{strsplit}}
#' @rdname build_mbg_formula_with_priors
#' @note The \code{formula} object was created at the very end of the function after dealing with everything before that in strings.
#'
#' @export
build_mbg_formula_with_priors <- function(fixed_effects,
                                          positive_constrained_variables = NULL,
                                          interact_with_year = NULL,
                                          int = TRUE,
                                          nullmodel = FALSE,
                                          add_nugget = FALSE,
                                          nugget_prior = "list(prior = 'loggamma', param = c(2, 1))",
                                          add_ctry_res = FALSE,
                                          ctry_re_prior = "list(prior = 'loggamma', param = c(2, 1))",
                                          ctry_re_sum0 = FALSE,
                                          spde_integrate0 = FALSE,
                                          temporal_model_type="'ar1'",
                                          temporal_model_theta_prior = "list(prior = 'loggamma', param = c(1, 0.00005))",
                                          temporal_model_theta1_prior = "list(prior = 'normal', param = c(0, 1/(2.58^2)))",
                                          no_gp = FALSE,
                                          stacker_names = child_model_names,
                                          coefs.sum1 = FALSE,
                                          subnat_RE = FALSE,
                                          subnat_re_prior = "list(prior = 'loggamma', param = c(1, 5e-5))",
                                          use_space_only_gp = FALSE,
                                          use_time_only_gmrf = FALSE,
                                          time_only_gmrf_type = "rw2",
                                          nid_RE = as.logical(use_nid_res)) {
  
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
  
  f_nugget <- as.formula(paste0('~f(IID.ID, model = "iid", hyper = list(theta=', nugget_prior, '), constr = TRUE)'))
  f_nid_res <- as.formula(paste0('~f(NID, model = "iid", hyper = list(theta=', nugget_prior, '), constr = TRUE)'))
  f_res <- as.formula(paste0('~f(CTRY.ID, model = "iid", hyper = list(theta=', ctry_re_prior, '), constr = ', as.logical(ctry_re_sum0), ')'))
  f_subnat <- as.formula(paste0('~f(SUBNAT.ID, model = "iid", hyper = list(theta=', subnat_re_prior, '), constr = TRUE)'))  
  
  ## spatial gps correlated across grouping (usually time)
  test_rho_priors(temporal_model_theta1_prior) ## Report how priors for theta1 (Rho) are being used in Rho space.
  f_space_time <- as.formula(paste0('~f(space,
                                    model = spde,
                                    group = space.group,
                                    constr = ', as.logical(spde_integrate0), ',
                                    control.group = list(model = ', temporal_model_type, ", ",
                                    "hyper = list(theta = ", temporal_model_theta_prior, ", theta1 = ", temporal_model_theta1_prior, ")))"))
  ## space only gp
  f_space <- as.formula('~f(sp.no.t, model = spde.sp)')
  
  ## time only gp
  f_time <- as.formula(paste0('~f(t.no.sp, model = "',time_only_gmrf_type,'")'))
  
  f.e.v <- base::strsplit(all_fixed_effects, ' + ', fixed=T)[[1]] ## fixed effect vector
  f_sum1 <- as.formula(paste0('~f(covar,
                              model = \'iid\',
                              extraconstr = list(A = matrix(1, 1, ', length(f.e.v), '), e = 1),
                              hyper=list(theta=list(initial=log(inla.set.control.fixed.default()$prec),
                              fixed=TRUE)))'))
  
  if(nchar(fixed_effects) <= 1){ ## then we don't have any fixed effects, so run the nullmodel
    nullmodel <- TRUE
  }
  
  ## build formula starting with f_null
  
  f_mbg <- f_null
  
  if(!nullmodel){
    if(coefs.sum1){
      f_lin <- f_sum1 ## use the sum1 formula instead of fixed effects.
      # the'fixed' effects may now be
      # found @ inla$fit$summary.random
    } else{
      f_lin <- reformulate(fixed_effects)
    }
    
    f_mbg <- f_mbg + f_lin
  }
  
  if(!no_gp) f_mbg <- f_mbg + f_space_time
  if(use_space_only_gp) f_mbg <- f_mbg + f_space
  if(use_time_only_gmrf) f_mbg <- f_mbg + f_time
  if(add_nugget==TRUE) f_mbg <- f_mbg + f_nugget
  if(nid_RE == TRUE) f_mbg <- f_mbg + f_nid_res
  if(add_ctry_res == TRUE) f_mbg <- f_mbg + f_res
  if(subnat_RE == TRUE) f_mbg <- f_mbg + f_subnat
  
  message(f_mbg)
  return(f_mbg)
}

## build_mbg_data_stack ################################################

#' @title Build an INLA data stack for modeling
#'
#' @description Generates an INLA stack that contains information
#'   mapping data to various INLA effects (GP, random effects, ...)
#'
#' @param df data.frame/table containing observational data. must have
#'   latitude and longitude columns as well as extracted values for
#'   all covariates that will be used in modeling
#'
#' @param fixed_effects string containing fixed effects used in
#'   res_fit separated by " + "
#' 
#' @param mesh_s An 'inla.mesh' object used to fit the spatial SPDE GP
#'   approximation in res_fit
#'
#' @param mesh_t An 'inla.mesh' object used to fit the temporal
#'   correlation structure in res_fit. if NULL, 
#'
#' @param exclude_cs vector of strings detailing covariates that are
#'   not center-scaled
#'
#' @param spde_prior String containing a list. Specifies the type of prior used
#'   for the Matern model parameters. If type = "nonpc" then use non pc 
#'   prior. User can specify nominal prior means: spde_prior$prior$range.nominal 
#'   and spde_prior$prior$variance.nominal. Defaults are INLA defaults, that is 
#'   20\% of the range of the mesh and 1, respectively. If type = "pc" then use 
#'   pc prior.  User can specify spde_prior$prior$range = (range0, Prange) i.e. 
#'   P(range < range0) = Prange and spde_prior$prior$sigma = (sigma0, Psigma) 
#'   i.e. P(sigma > sigma0) = Psigma. Defaults for Prange and Psigma are 0.05. 
#'   Default for range0 is 5\% of the range of the mesh and sigma0=3.
#'
#' @param coefs.sum1 Logical. If TRUE, add a constraint to ensure
#'   covariate coefs sum to 1 in fitted model
#'
#' @param use_ctry_res Logical. If TRUE, include country random
#'   effects
#'   
#' @param use_subnat_res Logical. If TRUE, include subnational 
#'   (admin-1) random effects
#'   
#' @param use_subnat_res Logical. If TRUE, include subnational (admin-1) random
#'   effects
#'   
#' @param remove_non_subnats Logical. If TRUE, NA out all points of the subnational
#'   random effect which are not associated with the country
#'   
#' @param use_nugget Logical. If TRUE, include a nugget effect
#' 
#' @param stacker_names string vector of child model names
#'
#' @param yl numeric vector of years in model
#'
#' @param zl numeric vector of third dimension in kronecker
#'   product. Only implemented in TMB
#'
#' @param zcol column name of z values associated with zl. Only
#'   implemented in TMB
#'
#' @param scale_gaussian_variance_N Logical. Do you want to scale
#'   gaussian variance by sample size? Only implemented in TMB.
#'
#' @param tmb Logical. use tmb?
#' 
#' @param shapefile_version character. Version of shape file to use
#' 
#' @param cov_constraints named int vector. integer vector indexed by covariate names.
#'
#' @param use_space_only_gp Logical. include a space only (time stationary) gp
#' 
#' @param use_time_only_gmrf Logical. include a time only (time stationary) gp
#'
#' @param st_gp_int_zero Logical. Should the space-time GP be forced to interate to 0? Default: FALSE
#' 
#' @param s_gp_int_zero Logical. Should the space GP be forced to interate to 0? Only used if use_space_only_gp=TRUE. Default: FALSE
#'
#' @return List containing 1) 'inla.stack' object, 2) 'inla.spde'
#'   object, 3) cs_df: a center-scale dataframe containing info on
#'   mean and SD used to scale covaraites in design matrix

build_mbg_data_stack <- function(df, fixed_effects, mesh_s, mesh_t,
                                 exclude_cs='',
                                 spde_prior="list(type='nonpc')",
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
                                 cov_constraints = covariate_constraint_vectorize(config),
                                 use_gp = TRUE, 
                                 use_space_only_gp = FALSE,
                                 use_time_only_gmrf = FALSE,
                                 shapefile_version = 'current',
                                 st_gp_int_zero = FALSE,
                                 s_gp_int_zero = FALSE,
                                 nid_RE = as.logical(use_nid_res)
) {
  
  if(nchar(stacker_names[1]) == 0 & coefs.sum1 == TRUE){
    message("WARNING! You've chosen sum-to-1 but don't appear to be using any stackers. Unless you have a very good reason to do this, it probably doesn't make sense. As such, we're setting coefs.sum1 <- FALSE")
    coefs.sum1 <- FALSE
  }
  
  if(use_gp & use_space_only_gp & !st_gp_int_zero & !s_gp_int_zero){
    message("WARNING! You've chosen to use a s-t and a s-only gp and have set both integration constraints to FALSE. This presents identifiability issues so we ase turning the space-only integrate-to-0 constraint to TRUE (i.e. s_gp_int_zero <- TRUE)")
    s_gp_int_zero <- TRUE
  }
  
  # if fitting with tmb, punt this over to
  if(tmb == TRUE){
    message('Returning a TMB model stack.')
    return(
      build_mbg_data_stack_tmb(d          = df,
                               yl         = yl,                 # year list
                               fes        = fixed_effects,  # fixed effects in the model
                               indic      = indicator,          # indicator
                               exclude_cs = exclude_cs,
                               nugget     = use_nugget,
                               country_re = use_ctry_res,
                               nid_re     = use_nid_res,
                               zl         = zl,
                               zcol       = zcol,
                               shapefile_version = shapefile_version,
                               scale_gaussian_variance_N = scale_gaussian_variance_N,
                               mesh       = mesh_s,            # spatial mesh
                               cov_constraints = cov_constraints))
    
    # else do the inla version
  } else {
    
    if(use_time_only_gmrf){
      # if we are using a time only effect then we need to make sure all 
      # year effects are estimated
      # the NA observations do not contribute to the model fitting but they 
      # prevent INLA from auto-removing random effects that (conditionally) 
      # have no data impacting their fit
      if(!all(yl %in% unique(df[["year"]]))){
        missyears <- yl[!(yl %in% unique(df[["year"]]))]
        missDF <- sapply(missyears, function(y){
          tempDF <- df[1,] # copy the first row of df
          tempDF[,indicator] <- NA # set the outcome to NA
          tempDF[,"year"] <- y # substitute the missing year
          tempDF
        })
        df <- cbind(df, missDF)
        
      }
    }
    
    # construct an SPDE model with a Matern kernel for the space-time GP
    message('Building SPDE...')
    
    spde_prior <- eval(parse(text=spde_prior)) # convert from string to list
    
    spde_list <- build_spde_prior(spde_prior, mesh_s, st_gp_int_zero)
    spde_prior <- spde_list$spde_prior
    spde <- spde_list$spde
    
    ## Build projector matrix between data locs and spatial mesh
    data.locs <- as.matrix(df[, c('longitude', 'latitude'),with=F])
    if(mesh_s$manifold == "S2"){
      ## then the mesh is on the sphere and we need to use 3d coords
      data.locs <- lonlat3D(data.locs[, 1], data.locs[, 2])
    }
    
    ## here we actually build the projector matrix, A
    ## it is grouped across periods in time
    A <- inla.spde.make.A(mesh = mesh_s,
                          loc = data.locs,
                          group = df$period,
                          group.mesh = mesh_t
    )
    
    if(coefs.sum1){
      ## make A matrix comprised of covariate fixed_effects column vectors
      f.e.v <- stacker_names ## fixed eff. vec.
      A.covar <- as.matrix(df[, f.e.v, with = FALSE])
    }
    
    ## confusingly, this 'space' index is actually space indices in time
    ## and it is a space index if only one time period is used
    if (is.null(mesh_t)){
      space = inla.spde.make.index("space",
                                   n.spde = spde$n.spde)
    } else {
      space = inla.spde.make.index("space",
                                   n.spde = spde$n.spde,
                                   n.group = mesh_t$m)
    }
    
    ## make another set of objects just for space only
    ## make another projection matrix that is fixed across all time
    A.sp <- inla.spde.make.A(mesh = mesh_s,
                             loc = data.locs)
    
    ## make another spde object in case you want one (st or s) to
    ## int to 0 but the other one to be unconstrained
    if(spde_prior$type=="pc"){ # PC prior
      spde.sp <- inla.spde2.pcmatern(mesh = mesh_s, 
                                     alpha = 2,
                                     prior.range = spde_prior$prior$range,
                                     prior.sigma = spde_prior$prior$sigma,
                                     constr = s_gp_int_zero)
    } else { # Non PC prior
      spde.sp <- inla.spde2.matern(mesh = mesh_s,  alpha = 2, constr = s_gp_int_zero,
                                   prior.range.nominal = spde_prior$prior$range.nominal,
                                   prior.variance.nominal = spde_prior$prior$variance.nominal)
    }
    
    ## here we make the actual space (time stationary) index
    ## which can be used even if there are multiple time points
    sp.no.t <- inla.spde.make.index("sp.no.t",
                                    n.spde = spde.sp$n.spde)
    
    #find cov indices
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
    
    if (nid_RE) {
      design_matrix$NID <- (df$nid)     
    }
    
    if(use_time_only_gmrf){
      design_matrix$t.no.sp <- df$year
    }
    
    ## initialize some lists for the stacking function
    ## this assumes you alway have the s-t gp and the covariates
    A.list <- list(A, 1) ## this list contains the objects that map between the effects and the data
    e.list <- list(space, design_matrix) ## this list contains the effects [or their indices for things in f()s]
    
    ## add on sum-to-1 for covars if selected
    if(coefs.sum1 == TRUE & nchar(fixed_effects) > 1){
      A.list <- c(A.list, list(A.covar = A.covar))
      e.list <- c(e.list, list(covar = 1:ncol(A.covar)))
    }
    
    ## add on space only selected
    if(use_space_only_gp){
      A.list <- c(A.list, list(A.sp = A.sp))
      e.list <- c(e.list, list(sp.no.t = sp.no.t))
    }
    
    message('Stacking data...')
    ## combine all the different pieces
    stack.obs <- inla.stack(data = list(covered = cov),
                            A = A.list,
                            effects = e.list,
                            tag = 'est'
    )
    
    return_list <- list(stack.obs, spde, cs_df, spde.sp)
    
    return(return_list)
  }
}


#' @title Fit MBG model in INLA
#' @description Fit MBG model in INLA
#' @param indicator_family PARAM_DESCRIPTION
#' @param stack.obs PARAM_DESCRIPTION
#' @param spde PARAM_DESCRIPTION
#' @param cov PARAM_DESCRIPTION
#' @param N PARAM_DESCRIPTION
#' @param int_prior_mn PARAM_DESCRIPTION
#' @param int_prior_prec PARAM_DESCRIPTION, Default: 1
#' @param f_mbg PARAM_DESCRIPTION
#' @param run_date PARAM_DESCRIPTION
#' @param keep_inla_files PARAM_DESCRIPTION
#' @param cores PARAM_DESCRIPTION
#' @param verbose_output PARAM_DESCRIPTION, Default: FALSE
#' @param wgts PARAM_DESCRIPTION, Default: 0
#' @param intstrat PARAM_DESCRIPTION, Default: 'eb'
#' @param fe_mean_prior PARAM_DESCRIPTION, Default: 0
#' @param fe_sd_prior PARAM_DESCRIPTION, Default: 1
#' @param omp_start OpenMP strategy to be used, one of \code{c('small', 'medium', 'large', 'huge', 'default', 'pardiso.serial', 'pardiso.parallel')}.
#' Default: \code{'huge'}, which is the INLA default for TAUCS solver.
#' @param blas_cores Number of BLAS threads to use. Default: 1.
#' @param pardiso_license Path to PARDISO license. Must be provided if using a Pardiso startegy in \code{omp_strat}. Default: NULL
#' @param sparse_ordering boolean: should the output precision matrix have a
#'   sparse (metis) ordering?
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # EXAMPLE1
#' }
#' }
#' @rdname fit_mbg
#' @export
fit_mbg <- function(indicator_family,
                    stack.obs,
                    spde,
                    cov,
                    N,
                    int_prior_mn,
                    int_prior_prec = 1,
                    f_mbg,
                    run_date,
                    keep_inla_files,
                    cores, 
                    verbose_output = FALSE,
                    wgts = 0,
                    intstrat = "eb",
                    fe_mean_prior = 0,
                    fe_sd_prior = 1,
                    omp_strat = "huge",
                    blas_cores = 1,
                    pardiso_license = NULL,
                    sparse_ordering = TRUE) {
  if (wgts[1] == 0) {
    wgts <- rep(1, length(N))
  }
  
  ## Assert that omp_strat is in the domains of what is possible
  stopifnot(omp_strat %in% c("small", "medium", "large", "huge", "default", "pardiso.serial", "pardiso.parallel"))
  
  
  message("Checking to see if PARDISO strategy is used")
  if (grepl("pardiso", omp_strat)) {
    
    ## Pardiso strategy NEEDS a license
    stopifnot(!is.null(pardiso_license))
    
    ## Set up and test PARDISO solver
    inla.setOption("pardiso.license", pardiso_license)
    (success_msg <- inla.pardiso.check())
    
    if (!grepl("SUCCESS", success_msg)) {
      warning("PARDISO solver was not set up properly. Reverting to 'huge' strategy")
      omp_strat <- "huge"
    } else {
      message(paste0("Pardiso solver being used with strategy ", omp_strat))
    }
  }
  
  
  # Determine control.inla settings
  control.inla.list <- list(int.strategy = intstrat, h = 1e-3, tolerance = 1e-6)
  if(sparse_ordering==TRUE) control.inla.list$reordering <- 'METIS'
  
  # ~~~~~~~~~~~~~~~~
  # fit the model
  # enable weights
  inla.setOption("enable.inla.argument.weights", TRUE)
  
  message("Fitting INLA model")
  
  # set a prior variance of 1.96 on the intercept as this is
  # roughly as flat as possible on logit scale without >1 inflection
  
  # code just to fit the model (not predict)
  inla_working_dir <- paste0("<<<< FILEPATH REDACTED >>>>/inla_", run_date)
  dir.create(inla_working_dir, showWarnings = FALSE)
  if (indicator_family == "binomial") {
    system.time(
      res_fit <- inla(f_mbg,
                      data = inla.stack.data(stack.obs),
                      control.predictor = list(
                        A = inla.stack.A(stack.obs),
                        link = 1,
                        compute = FALSE
                      ),
                      control.fixed = list(
                        expand.factor.strategy = "inla",
                        mean = list(int = int_prior_mn,
                                    default = fe_mean_prior),
                        prec = list(int = int_prior_prec,
                                    default = fe_sd_prior)
                      ),
                      control.compute = list(
                        dic = TRUE,
                        cpo = TRUE,
                        config = TRUE,
                        openmp.strategy = omp_strat
                      ),
                      control.inla = control.inla.list,
                      family = "binomial",
                      num.threads = cores,
                      blas.num.threads = blas_cores,
                      Ntrials = N,
                      verbose = verbose_output,
                      working.directory = inla_working_dir,
                      weights = wgts,
                      keep = TRUE
      )
    )
  }
  if (indicator_family == "gaussian") {
    system.time(
      res_fit <- inla(f_mbg,
                      data = inla.stack.data(stack.obs),
                      control.predictor = list(
                        A = inla.stack.A(stack.obs),
                        compute = FALSE
                      ),
                      control.fixed = list(
                        expand.factor.strategy = "inla",
                        prec.intercept = 1,
                        mean.intercept = int_prior_mn,
                        mean = 0,
                        prec = 2
                      ),
                      control.compute = list(
                        dic = TRUE,
                        cpo = TRUE,
                        config = TRUE,
                        openmp.strategy = omp_strat
                      ),
                      control.inla = control.inla.list,
                      family = "gaussian",
                      num.threads = cores,
                      blas.num.threads = blas_cores,
                      verbose = TRUE,
                      working.directory = inla_working_dir,
                      weights = wgts,
                      scale = N,
                      keep = TRUE
      )
    )
  }
  if (indicator_family == "beta") {
    system.time(
      res_fit <- inla(f_mbg,
                      data = inla.stack.data(stack.obs),
                      control.predictor = list(
                        A = inla.stack.A(stack.obs),
                        compute = FALSE
                      ),
                      control.fixed = list(
                        expand.factor.strategy = "inla",
                        prec.intercept = 1,
                        mean.intercept = int_prior_mn,
                        mean = 0,
                        prec = 2
                      ),
                      control.compute = list(
                        dic = TRUE,
                        cpo = TRUE,
                        config = TRUE,
                        openmp.strategy = omp_strat
                      ),
                      control.inla = control.inla.list,
                      family = "beta",
                      num.threads = cores,
                      blas.num.threads = blas_cores,
                      verbose = TRUE,
                      working.directory = inla_working_dir,
                      weights = wgts,
                      scale = N,
                      keep = TRUE
      )
    )
  }
  
  return(res_fit)
}




## predict_mbg ################################################

#' @title Predict MBG (INLA) draws
#'
#' @description Generate draws from a fitted INLA MBG model across a
#'   space-time modeling domain. Draws are generated at the same
#'   resolution as simple_raster (probably 5x5km).
#'
#' @param res_fit A fitted INLA model object
#'
#' @param cs_df A dataframe containing mean and standard deviations used to scale the covariates
#'
#' @param mesh_s An 'inla.mesh' object used to fit the spatial SPDE GP approximation in res_fit
#'
#' @param mesh_t An 'inla.mesh' object used to fit the temporal correlation structure in res_fit
#'
#' @param cov_list A named list containing covariate rasters used in the model fit
#'
#' @param samples Integer number of draws to take from the fitted posterior distribution
#'
#' @param simple_raster A simple_raster defining the modeling domain pixels
#'
#' @param transform String name of transformation to be applied to
#'   draws. Possible options are 'inverse-logit'. Any other string
#'   will result in the identity transform.
#'
#' @param ind indicator string 
#'
#' @param indg indicator_group string
#'
#' @param rd run_date string
#'
#' @param region region string
#'
#' @param full_df data.frame/table containing observational data used to fit model
#'
#' @param coefs.sum1 Logical. Were the coefficients constrained to sum-to-1 in the model fit
#'
#' @param pred_gp Logical. Should gaussian process draws be added to the draw objects.
#'
#' @param fixed_effects string containing fixed effects used in res_fit separated by " + "
#'
#' @param yl numeric vector of years in model
#'
#' @param nperiod number of temporal periods in res_fit
#'
#' @param stacker_names string vector containing names of child models used for stacking in res_fit
#'
#' @param shapefile_version Shapefile version
#' 
#' @param simple_raster_subnats simple raster for the subnational country
#' 
#' @param subnat_country_to_get ISO code for subnational country
#'
#' @param use_space_only_gp Logical. include a space only (time stationary) gp
#' 
#' @param use_time_only_gmrf Logical. include a time only gp
#' 
#' @return list containing mean raster, sd raster, and a matrix of
#'   space-time pixel cell predictions. The cell prediction matrix is
#'   wide on draws and long on space-time locations. The space-time
#'   rows are formatted such that all pixels in the first period take
#'   up rows 1:num.space.pixels, Then period 2 space pixels are next
#'   and so on.
#' 

predict_mbg <- function(res_fit, cs_df, mesh_s, mesh_t, cov_list,
                        samples, simple_raster, transform,
                        ind=indicator, indg=indicator_group,
                        rd=run_date, region = reg, full_df = df,
                        coefs.sum1 = FALSE, pred_gp = TRUE,
                        fixed_effects = all_fixed_effects, yl = year_list,
                        nperiod = ifelse(is.character(year_list),
                                         length(eval(parse(text = year_list))),
                                         length(year_list)),
                        stacker_names = child_model_names,
                        shapefile_version = 'current',
                        simple_raster_subnats = NULL,
                        subnat_country_to_get = NULL,
                        use_space_only_gp = FALSE,
                        use_time_only_gmrf = FALSE
) {
  
  if(nchar(stacker_names[1]) == 0 & coefs.sum1 == TRUE){
    message("WARNING! You've chosen sum-to-1 but don't appear to be using any stackers. Unless you have a very good reason to do this, it probably doesn't make sense. As such, we're setting coefs.sum1 <- FALSE")
    coefs.sum1 <- FALSE
  }
  
  ## Now that we've extracted covariate values to our data in the buffer zone, clip cov list to simple_raster instead of simple_polygon
  ##    (simple_raster is area we actually want to model over)
  for(l in 1:length(cov_list)) {
    message(sprintf("On cov %i out of %i", l, length(cov_list)))
    ## res(cov_list[[l]]) <- res(simple_raster) ## to fix strange issue with PERU covs having different x-resolution
    cov_list[[l]]  <- crop(cov_list[[l]], extent(simple_raster))
    cov_list[[l]]  <- setExtent(cov_list[[l]], simple_raster)
    cov_list[[l]]  <- mask(cov_list[[l]], simple_raster)
  }
  
  #updated by DCCASEY on 11-21-2016
  
  
  message('Making predictions')
  
  ## number of samples
  n_draws <- samples
  
  ## dummy raster
  template <- simple_raster
  
  cell_idx <- seegSDM:::notMissingIdx(template)
  
  ## sample from posterior over latents
  suppressWarnings(draws <- inla.posterior.sample(n_draws, res_fit))
  
  ## Save full inla draws object above to get hyperpar draws later (tau for Gaussian, etc.)
  run_dir <- paste0('<<<< FILEPATH REDACTED >>>>', indg, '/', ind, '/output/', rd, '/inla_draws')
  dir.create(run_dir, showWarnings = FALSE)
  save(list='draws', file = paste0(run_dir, '/inla_draws_', region, '.RData'))
  
  ## get parameter names
  par_names <- rownames(draws[[1]]$latent)
  
  ## index to spatial-temporal field and linear coefficient samples
  s_idx <- grep('^space.*', par_names) ## space-time random effects
  
  ## index for spatial field
  s_no_t_idx <- grep('^sp.no.t.*', par_names) ## space random effects
  
  if(coefs.sum1){
    l_idx <- grep('covar', par_names)
    l_idx_int <- match(sprintf('%s.1', res_fit$names.fixed), ## main effects
                       par_names)
    if(mean(is.na(l_idx_int)) == 1) l_idx_int <- match(sprintf('%s', res_fit$names.fixed), par_names)
    if(mean(is.na(l_idx_int)) == 1) l_idx_int <- match(sprintf('%s:1', res_fit$names.fixed), par_names)##
    l_idx = c(l_idx, ## covars in sum-to-1
              l_idx_int) ## intercept
    
    ## check
    if(length(l_idx) < 1){ ## 1 b/c intercept still will show up
      stop('Looks like you may have requested sum.to.1 in predict_mbg on a model that was not run using sum.to.1 constraint.')
    }
    
  }else{
    l_idx <- match(sprintf('%s.1', res_fit$names.fixed), ## main effects
                   par_names)
    if(mean(is.na(l_idx)) == 1) l_idx <- match(sprintf('%s', res_fit$names.fixed), par_names) ## fix for different INLA version naming conventions
    
    ## check
    if(length(l_idx) < 1){ ## 1 b/c of the intercept
      stop('Looks like you may have set coefs.sum.1 to FALSE in predict_mbg using a fitted model that WAS run using sum.to.1 constraint')
    }
  }
  
  ## get samples of space-time REs as matrices
  if(pred_gp){
    pred_s <- sapply(draws, function(x) {
      x$latent[s_idx]})
  }
  
  ## get samples of space REs as matrices
  if(use_space_only_gp){
    pred_s_no_t <- sapply(draws, function (x)
      x$latent[s_no_t_idx])
  }
  
  ## get samples of other effects/params as matrices
  pred_l <- sapply(draws, function(x) {
    x$latent[l_idx]})
  if(length(l_idx)==1) pred_l=t(as.matrix(pred_l))
  
  if(coefs.sum1){
    rownames(pred_l) <- c(stacker_names, ## covars in sum-to-1
                          res_fit$names.fixed) ## intercept
  } else {
    rownames(pred_l) <- res_fit$names.fixed
  }
  
  ## if we fit with a nugget and include it in prediction, we also need to take draws of the nugget precision
  if(length(grep('^IID.ID.*', par_names)) > 0){
    pred_n <- sapply(draws, function(x) {
      nug.idx <- which(grepl('IID.ID', names(draws[[1]]$hyper)))
      x$hyperpar[[nug.idx]]}) ## this gets the precision for the nugget
  }else{
    pred_n <- NULL
  }
  
  ## check to see if country random effects were used
  if(sum(grepl('CTRY.ID', par_names)) > 0){
    
    ## then country random effects were used - get the fitted values
    ctry.res.idx <- grep('CTRY.ID*', par_names)
    pred_ctry_res  <-  sapply(draws, function (x) x$latent[ctry.res.idx])
    ctry.res.names <- par_names[ctry.res.idx]
    if (length(ctry.res.names) == 1) {
      message('WARNING:ONLY ONE COUNTRY PRESENT IN DATA AND RE FIT')
      print('WARNING:ONLY ONE COUNTRY PRESENT IN DATA AND RE FIT')
      pred_ctry_res <- matrix(pred_ctry_res, nrow = 1)
    }
    rownames(pred_ctry_res) <- res_fit$summary.random$CTRY.ID$ID
    
    ## also get the hyperparam precision in case not all countries had data and we need to take draws
    pred_ctry_prec <- sapply(draws, function(x){
      ctry.prec.idx <- which(grepl('CTRY.ID', names(draws[[1]]$hyper)))
      x$hyperpar[[ctry.prec.idx]]}) ## this gets the precision for country REs
  }else{
    pred_ctry_res <- NULL
  }
  
  if(use_time_only_gmrf){
    
    ## then main temporal effects were used - get the fitted values
    t_no_s_idx <- grep('^t.no.sp.*', par_names) ## time random effects
    pred_time_res  <-  sapply(draws, function (x) x$latent[t_no_s_idx])
    time.res.names <- par_names[t_no_s_idx]
    if (length(time.res.names) == 1) {
      message('WARNING:ONLY ONE YEAR PRESENT IN DATA AND RE FIT')
      print('WARNING:ONLY ONE YEAR PRESENT IN DATA AND RE FIT')
      pred_time_res <- matrix(pred_time_res, nrow = 1)
    }
    rownames(pred_time_res) <- res_fit$summary.random$t.no.sp$ID
    
  }else{
    pred_time_res <- NULL
  }
  
  ## check to see if admin-1 (subnat) random effects were used
  if (sum(grepl("SUBNAT.ID", par_names)) > 0) {
    
    ## then country random effects were used - get the fitted values
    subnat.res.idx <- grep("SUBNAT.ID*", par_names)
    pred_subnat_res <- sapply(draws, function(x) x$latent[subnat.res.idx])
    subnat.res.names <- par_names[subnat.res.idx]
    if (length(subnat.res.names) == 1) {
      message("WARNING:ONLY ONE SUBNAT UNIT PRESENT IN DATA AND RE FIT")
      print("WARNING:ONLY ONE SUBNAT UNIT PRESENT IN DATA AND RE FIT")
      pred_subnat_res <- matrix(pred_subnat_res, nrow = 1)
    }
    rownames(pred_subnat_res) <- res_fit$summary.random$SUBNAT.ID$ID
    
    ## also get the hyperparam precision in case not all countries had data and we need to take draws
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
    subnat_coords <- raster::xyFromCell(template_subnat, seegSDM:::notMissingIdx(template_subnat))
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
  #nperiod <- mesh_t$n
  #nperiod <- max(mesh_t$loc)
  ## We now take this as a argument based on the length of the year_list vector
  
  ## replicate coordinates and years for GP projection
  gp_coords_periods <- do.call(rbind,
                               replicate(nperiod,
                                         gp_coords,
                                         simplify = FALSE))
  
  groups_periods <- rep(1:nperiod,
                        each = nrow(gp_coords))
  
  ## Projector matrix for space-time effects
  A.pred <- inla.spde.make.A(mesh = mesh_s,
                             loc = gp_coords_periods,
                             group = groups_periods,
                             group.mesh = mesh_t)
  
  A.pred.s.no.t <- inla.spde.make.A(mesh = mesh_s,
                                    loc = gp_coords)
  
  ## get samples of s-t effect for all cells and across time
  if(pred_gp){
    # Using `crossprod()` here for product of a large, sparse matrix and a dense matrix is ~2x
    # faster than using `%*%`
    # cell_s <- as.matrix(crossprod(t(A.pred), pred_s))
    cell_s <- as.matrix(Matrix::crossprod(Matrix::t(A.pred), pred_s))
  }else{
    cell_s <- 0
  }
  
  ## get samples of the spatial effect at all cells and replicate across time
  if(use_space_only_gp){
    ## get the stationary spatial effect
    cell_s_no_t <- as.matrix(Matrix::crossprod(Matrix::t(A.pred.s.no.t), pred_s_no_t))
    ## and replicate across time
    cell_s_no_t <- do.call('rbind', rep(list(cell_s_no_t), nperiod))
  }else{
    cell_s_no_t <- 0
  }
  
  # remove out temporally varying covariates for now, will deal with them later
  #split the effects names by varying or time varying
  tvnames = pars =c()
  for(i in 1:length(cov_list)){
    if(dim(cov_list[[i]])[3]!=1) tvnames[length(tvnames)+1]= names(cov_list)[i]
    if(dim(cov_list[[i]])[3]==1) pars[length(pars)+1]= names(cov_list)[i]
  }
  
  #now split the covariates themselves (not just the names)
  # process covariates
  
  ## start by ensuring the same mask and cell size
  ## I think we don't need this next line anymore since everything should be done at the top of this function...
  cov_list = setNames(lapply(cov_list, function(x) {
    mask(resample(x, template), template)
  }), names(cov_list))
  
  #extract the covariates
  vals = data.table(do.call(cbind, lapply(cov_list, function(x) raster::extract(x,coords))))
  
  #create the int column
  vals[,int:= 1]
  
  #reshape long
  vals[,id:= 1:nrow(vals)] #create an id to ensure that melt doesn't sort things
  
  if (!is.null(mesh_t)){
    #convert the time varying names into meltable values
    tv_cov_colnames = grep(paste(tvnames, collapse = "|"), names(vals), value = T) #unlisted
    tv_cov_colist = lapply(tvnames, function(x) grep(paste0("^", x, "\\.[[:digit:]]*$"), names(vals), value = T))
    
    vals = melt(vals, id.vars = c('id','int',pars), measure = tv_cov_colist, value.name = tvnames, variable.factor =F)
    
    #fix the names
    #melt returns different values of variable based on if its reshaping 1 or 2+ columns.
    #enforce that it must end with the numeric after the period
    vals[,variable:= as.numeric(substring(variable,regexpr("(\\.[0-9]+)$", variable)[1]+1))]
  }
  
  #Tests suggest this is not needed anymore/as a carry over. Keep just in case
  #if(length(tv_cov_colist)==1){
  #  setnames(vals,paste0(tvnames[[1]],'1'), tvnames)
  #}
  
  #keep only columns that were part of the regression
  vals = vals[,mget(rownames(pred_l))]
  
  #centreScale the values
  vals = centreScale(vals,cs_df) #happens for all periods at once
  
  #create the grid of fixed effects
  # Conversion to the "dgeMatrix" here will gives better performance for matrix
  # multiplication of two dense matrices. Conversion back to standard unnamed
  # 'matrix' object slows the operation down a bit, but ensures that the later
  # `apply` operations will be as fast as possible. Still much faster this way
  # and much faster when using the MKL.
  cell_l <- unname(as.matrix(as(data.matrix(vals), "dgeMatrix") %*% pred_l))
  
  ## add on nugget effect if applicable
  if(!is.null(pred_n)){
    cell_n <- sapply(pred_n, function(x){
      rnorm(n = nrow(cell_l), sd = 1 / sqrt(x), mean = 0)
    })
  }
  
  ## add on country random effects if applicable
  if(!is.null(pred_ctry_res)){
    
    cell_ctry_res <- matrix(0, ncol = ncol(cell_l), nrow = nrow(cell_l))
    
    ## loop through the countries in the region and add on random
    ## effects if we have a fitted value. otherwise add on iid values
    ## drawn with the precision for the random effects hyperparam
    reg.gauls <- get_adm0_codes(reg, shapefile_version = shapefile_version)
    gauls.with.res <- as.numeric(rownames(pred_ctry_res))
    gaul.vec <- values(simple_raster)[cell_idx]
    gaul.t.vec <- rep(gaul.vec, nperiod) ## expand gaul vec in time
    for(gg in reg.gauls){
      ## add on random effect to those pixels belonging to the country
      gg.rows <- which(gaul.t.vec == gg)
      
      if(gg %in% gauls.with.res){
        message(paste0('adding RE for country ', gg))
        gg.idx  <- which(gauls.with.res == gg)
        cell_ctry_res[gg.rows, ] <- cell_ctry_res[gg.rows, ] +
          matrix(rep(pred_ctry_res[gg.idx, ], length(gg.rows)), ncol = ncol(cell_ctry_res), byrow = TRUE)
      }else{
        ## this country had no data and we need to take a draw from the RE dist for draw in preds
        ## for draw, we draw one value of the random effect using pred_ctry_prec
        re.draws <- rnorm(n = ncol(cell_ctry_res), mean = 0, sd = 1 / sqrt(pred_ctry_prec))
        ## replicate by all rows in that country and add it on
        cell_ctry_res[gg.rows, ] <- cell_ctry_res[gg.rows, ] +
          matrix(rep(re.draws, length(gg.rows)), ncol = ncol(cell_ctry_res), byrow = TRUE)
      }
      
    }
    
  }
  
  # The final cell object should be arranged by space then time such that all
  # years are clustered together and in ascending order. 
  if(!is.null(pred_time_res)){
    cell_time_res <- matrix(0, ncol = ncol(cell_l), nrow = nrow(cell_l))
    gaul.vec.N <- length(values(simple_raster)[cell_idx])
    idx_hold <- 0
    for(y_ in sort(yl)){
      tidx <- ((idx_hold*gaul.vec.N)+1):((idx_hold+1)*gaul.vec.N)
      cell_time_res[tidx,] <- sapply(
        pred_time_res[as.character(y_),], 
        function(it) rep(it, gaul.vec.N))
      idx_hold <- idx_hold + 1
    }
    
  }
  
  ## add on subnational random effects if applicable
  if (!is.null(pred_subnat_res) & !is.null(subnat_country_to_get)) {
    cell_subnat_res <- matrix(0, ncol = ncol(cell_l), nrow = nrow(cell_l))
    
    ## We will identify pixels at where the parent country of the subnational is in
    ## and because simple_raster2 is perfectly aligned with a subset of simple_raster for that country,
    ## we can simply add on the REs for those 
    
    subnat_gauls.with.res <- as.numeric(rownames(pred_subnat_res))
    subnat_gaul.vec <- values(simple_raster2)[subnat_cell_idx]
    subnat_gaul.t.vec <- rep(subnat_gaul.vec, nperiod) ## expand gaul vec in time
    
    
    ## Get the rows where the simple_raster corresponds to the target
    ## subnat country
    # gg.rows <- which(values(simple_raster) == get_adm0_codes(subnat_country_to_get), arr.ind = T)
    # 
    # ## Get inverse of above (non subnat)
    # gg.not.rows <- which(values(simple_raster) != get_adm0_codes(subnat_country_to_get), arr.ind = T)
    
    for (gsubnat in subnat_gauls.with.res) {
      
      ## Find the rows in simple_raster2 corresponding to the iterating subnat
      gsubnat.rows <- which(subnat_gaul.t.vec == gsubnat)
      gg.idx <- which(subnat_gauls.with.res == gsubnat)
      
      if (gsubnat %in% subnat_gauls.with.res) {
        
        message(paste0("adding RE for subnat ", gsubnat))
        
        cell_subnat_res[gsubnat.rows, ] <- cell_subnat_res[gsubnat.rows, ] + 
          matrix(rep(pred_subnat_res[gg.idx, ], length(gsubnat.rows)), ncol = ncol(cell_subnat_res), byrow = TRUE)
        
      } else {
        
        ## this subnat unit had no data and we need to take a draw from the RE dist for draw in preds
        ## for draw, we draw one value of the random effect using pred_subnat_prec
        re.draws <- rnorm(n = ncol(cell_subnat_res), mean = 0, sd = 1 / sqrt(pred_subnat_prec))
        
        ## replicate by all rows in that country and add it on
        cell_subnat_res[gsubnat.rows, ] <- cell_subnat_res[gsubnat.rows, ] + 
          matrix(rep(re.draws, length(gsubnat.rows)), ncol = ncol(cell_subnat_res), byrow = TRUE)
      }
    }
  }
  
  
  ## get sd from fixed effects
  cell_l_sd <-  apply(cell_l, 1, sd)
  
  ## project model uncertainty from node-level sd to pred grid level for s-t REs
  if(pred_gp){
    node_s_sd <- apply(pred_s, 1, sd)
    # Large, sparse matrix-vector product faster as `%*%` than with `crossprod()`
    cell_s_sd <- as.matrix(A.pred %*% node_s_sd)
  }else{
    cell_s_sd <- 0
  }
  
  ## project model uncertainty from node-level sd to pred grid level for spatial REs
  if(use_space_only_gp){
    node_s_no_t_sd <- apply(pred_s_no_t, 1, sd)
    ## Large, sparse matrix-vector product faster as `%*%` than with `crossprod()`
    cell_s_no_t_sd <- as.matrix(A.pred.s.no.t %*% node_s_no_t_sd)
    ## replicate in time
    cell_s_no_t_sd <- rep(cell_s_no_t_sd, nperiod)
  }else{
    cell_s_no_t_sd <- 0
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
  
  ## get uncertainty from time only effects
  if(!is.null(pred_time_res)){
    cell_time_res_sd <- apply(cell_time_res, 1, sd)
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## combine to get cell-level sd
  ## ASSUMES INDEPENDENCE!
  
  # do all of this by adding squares
  cell_sd <- cell_l_sd ^ 2 + cell_s_sd ^ 2 + cell_s_no_t_sd ^ 2
  
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
  
  # add time independent if present
  if(!is.null(pred_time_res)){
    cell_sd <- cell_sd + cell_time_res_sd^2 
  }
  
  # finally, take the square root
  cell_sd <- sqrt(cell_sd)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # combine to produce draws (untransformed)
  cell_all <- cell_l + cell_s + cell_s_no_t
  
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
  
  # add time independent if present
  if(!is.null(pred_time_res)){
    cell_all <- cell_all + cell_time_res 
  }
  
  # get predictive draws on probability scale
  if(transform == 'inverse-logit') {
    cell_pred <- plogis(as.matrix(cell_all))
  } else {
    cell_pred <- eval(parse(text = transform))
  }
  
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
  
  names(mean_ras) <- paste0('period_', 1:nperiod)
  names(sd_ras)   <- paste0('period_', 1:nperiod)
  
  return_list <- list(mean_ras, sd_ras, cell_pred)
  
  return(return_list)
  
}

save_mbg_preds <- function(config, time_stamp, run_date, mean_ras, sd_ras, res_fit, cell_pred, df ,pathaddin="") {

  if(time_stamp==TRUE) output_dir <- '<<<< FILEPATH REDACTED >>>>'
  if(time_stamp==FALSE) output_dir <- '<<<< FILEPATH REDACTED >>>>'
  dir.create(output_dir, showWarnings = FALSE)

  # Save log of config file
  write.csv(config, paste0(output_dir,'/config.csv'), row.names = FALSE)

  if (!is.null(mean_ras)) {
    writeRaster(
      mean_ras,
      file = (paste0(output_dir, '/', indicator,'_prediction_eb',pathaddin)),
      overwrite = TRUE
    )
  }

  if (!is.null(sd_ras)) {
    # latent sd
    writeRaster(
      sd_ras,
      file = (paste0(output_dir, '/', indicator,'_sd_eb',pathaddin)),
      overwrite = TRUE
    )
  }

  # save model
  save(res_fit,
       file = (paste0(output_dir, '/', indicator,'_model_eb',pathaddin,'.RData')))
  # save draws (with compression) to recombine later
  save(
    cell_pred,
    file = (paste0(output_dir, '/', indicator,'_cell_draws_eb',pathaddin,'.RData')),
    compress = TRUE
  )
  # save training data
  write.csv(
    df,
    file = (paste0(output_dir, '/', indicator,'_trainingdata',pathaddin)),
    row.names = FALSE
  )

  # Write a an empty file to indicate done with this parallel script
  write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

}
