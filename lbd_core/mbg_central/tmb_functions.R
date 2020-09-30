


#' @title Build Stack for TMB MBG model
#' @description Organize Data and Parameter Stacks to fit a TMB MBG model
#'
#' @param d prepped model frame with no NAs
#' @param yl a vector of years for analysis (i.e. c(2001,2002,2003))
#' @param zl a vector of zcol for analysis (i.e. ages c(1,2,3,4,5). Must be integers starting with 1)
#' @param fes a string of model fixed effects of form 'var1 + var2 + var3' or as string vector
#'  corresponding to column names in d
#' @param indic indicator name, corresponding to the response variable in d
#' @param country_re TRUE/FALSE include country re. If true and there is a zcol then there will be random slope as well (todo)
#' @param nid_re TRUE/FALSE include re on survey. Need nid column in the data frame. defaults to use_nid_res from config.
#' @param exclude_cs character vector of covariates to exclude from centrescaling
#' @param mesh an inla mesh object
#' @param s2mesh Logical. Should the mesh be created on the surface of
#'   a sphere? If TRUE, s2params is used to specify mesh parameters
#'   instead of max_edge and mesh_offset
#' @param s2params string of 3 element numeric vector in R
#'   notation. e.g. "c(25, 500, 1000)". The entries describe the
#'   minimum triangle edge length allowed, hos far to extend the mesh
#'   beyond the 'simple' boundary, and the maximum allowed triangle
#'   edge length, respectively. Units are in kilometers. Used only if
#'   s2mesh=TRUE.
#' @param cov_constraints named int vector. integer vector indexed by covariate
#'   names, in the format returned by the function
#'   \code{\link{covariate_constraint_vectorize}}. NOTE: In the current
#'   implementation, priors for fixed effects are set to the same values
#'   regardless of whether the covariate fixed effect is constrained or not.
#'   Apply constraints with caution.
#'
#' @return returns a named list with Data and Parameters to be passed into fit_mbg_tmb()

build_mbg_data_stack_tmb <- function(d          = df,
                                     yl         = year_list,
                                     zl         = z_list,
                                     fes        = all_fixed_effects,
                                     indic      = indicator,
                                     country_re = use_country_res,
                                     nid_re     = use_nid_res,
                                     exclude_cs = '',
                                     nugget     = FALSE,
                                     zcol       = NULL,
                                     shapefile_version = 'current',
                                     scale_gaussian_variance_N = TRUE,
                                     mesh       = mesh_s,
                                     cov_constraints = covariate_constraint_vectorize(config)
                                     ){

  # ensure d is a dt
  d <- setDT(d)

  # zcol
  if(!zcol %in% colnames(d)){
    message('No Z column detected')
    d[[zcol]] <- 0
  } else if (is.null(zcol)) {
    message('Z column was set as null')
    d[[zcol]] <- 0
  }
  if( !all(unique(d[[zcol]]) %in% zl)) {
    message('WARNING: zl and d[[zcol]] do not completely match up.. ')
    # check that we there arent values in zcol not matching z_list
    d[, dropz := !get(zcol) %in% zl]
    if(any(d$dropz != FALSE)){
      message(sprintf('WARNING: Detected some z values in zcol (%s) which were not in the z_list',zcol))
      message(sprintf('WARNING: Due to this, dropping %i rows from the input data',sum(d$dropz==TRUE)))
      print(table(d$dropz,d$age))
      d <- subset(d, dropz == FALSE)
      d[, dropz := NULL]
    }
  }

  # make a fake data point with no weight for the max period and Z to fill out the GP
  d <- rbind(d, d[1,])
  d[[zcol]][nrow(d)]     <- max(zl)
  d[['period']][nrow(d)] <- length(yl)
  d[['weight']][nrow(d)] <- 0


  # look for z dimension
  num_z <- 1
  if(length(zl)>1){
    message(sprintf('More than one unique %s found, initiating Z in the GP', zcol))
    num_z <- length(zl)

    # set A proj grouping. The ordering here must match the ordering of epsilon_stz in the template
    grp <- setDT(expand.grid(1:length(yl), 1:max(zl)))
    setnames(grp,c('Var1','Var2'),c('period',zcol))
    grp[,group := 1:.N]
    d <- merge(d, grp, by = c('period',zcol), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
  } else {
    # set Aproj grouping to period if there is but one z value
    d$group <- d$period
  }

  # coordinates at data points. these are passed to TMB in long,lat so
  # we keep them and create another set for 3d coords
  coords   <- cbind(d$longitude,d$latitude)

  # if we have  mesh on s2, first convert coords to spherical to project to mesh
  data.locs <- coords ## long, lat
  if(mesh$manifold == "S2"){
    ## then the mesh is on the sphere and we need to use 3d coords
    data.locs <- lonlat3D(data.locs[, 1], data.locs[, 2])
  }

  # make a projection matrix from data to st mesh
  A.proj <- inla.spde.make.A(mesh  = mesh,
                             loc   = data.locs,
                             group = d$group)


  # make a clean design matrix. make sure all fes appear in d
  fes  <- unlist(strsplit(fes, ' \\+ '))
  if(!all(fes %in% names(d)))
    stop('Check your fes argument, not all covariate names appear in d.')
  if(length(fes)!=0) {
    X_xp <- as.matrix(cbind(int=1, d[,c(fes),with=FALSE]))
  } else {
    X_xp <- as.matrix(cbind(int=rep(1,nrow(d))))
  }

  # add in age fixed effects
  if(num_z > 1){
    message(sprintf('Adding fixed effects for levels of zcol (%s)',zcol))
    for(z in 2:num_z){
      X_xp <- cbind(X_xp, d[[zcol]] == z)
    }
    colnames(X_xp) <- c(colnames(X_xp)[colnames(X_xp)!=''],paste0('FE_z_level__',2:num_z))
    exclude_cs <- c(exclude_cs,paste0('FE_z_level__',2:num_z))
  }

  # cs_df. imports seegMBG
  cs_df <- getCentreScale(X_xp, exclude = c('int',exclude_cs))
  X_xp  <- centreScale(X_xp, df = cs_df)

  # get data range in case we want to clamp for prediction later
  clamper <- data.table(apply(X_xp,2,range))

  # sort nugget and RE indicators
  nugget     <- as.numeric(as.logical(nugget))
  country_re <- as.numeric(as.logical(country_re))
  nid_re     <- as.numeric(as.logical(nid_re)) # these do not get used in prediction


  # check there is more than one observed country or nid if REs for those are set
  if(length(unique(d$country)) == 1 & country_re == TRUE){
    message('WARNING: Only found one unique country in this data frame, so turning off country random effects.')
    country_re <- FALSE
  }
  if(length(unique(d$nid)) == 1 & nid_re == TRUE){
    message('WARNING: Only found one unique NID in this data frame, so turning off NID random effects.')
    nid_re <- FALSE
  }

  # make a table mappind the Admin0 code to a unique country random effect
  md    <- get_location_code_mapping(shapefile_version = shapefile_version)
  mdsub <- md[ihme_lc_id %in% unique(as.character(d$country)),]
  if(nrow(mdsub) != length(unique(as.character(d$country)))){
    message(sprintf('get_location_code_mapping() COUNTRY NAMES: %s',paste(sort(mdsub$ihme_lc_id),collapse=', ')))
    message(sprintf('IN DATA COUNTRY NAMES: %s',paste(sort(unique(as.character(d$country))),collapse=', ')))
    stop('get_location_code_mapping() and countries in data not of matching lengths')
  }
  cntry_re_map <- data.table(
                    country   = mdsub$ihme_lc_id,
                    adm_code = mdsub$ADM_CODE,
                    re_id     = 0:(nrow(mdsub)-1))
  cntry_re_vec <- cntry_re_map$re_id[match(as.character(d$country),cntry_re_map$country)]

  # make an nid_re mapping table
  nid_re_map <- unique(select(d, nid)) %>%
    mutate(re_id = 0:(n()-1))
  nidEFF <- select(d, nid) %>% left_join(nid_re_map, by="nid")
  nid_re_vec <- nidEFF$re_id
  countries_with1nid <- select(d, country, nid) %>%
     group_by(country) %>% # group by country to see...
     mutate(nidCount=n()) %>% # the number of nid units per country
     filter(nidCount==1) %>% ungroup %>% select(country) %>% unique
  if(country_re==TRUE & nid_re==TRUE & nrow(countries_with1nid) > 0) {
    message(paste("WARNING! The following countries", countries_with1nid$country,
                   "have only 1 NID. If you encounter convergence problems, you
                  may want to remove NID random effects."))
  }


  # set GP RE array, or matrix depending on if we have a z dimension
  if(num_z > 1) {
    Epsilon_stz <- array(0, dim=c(mesh$n,length(yl),num_z))
  } else {
    Epsilon_stz <- array(0, dim=c(mesh$n,length(yl)))
  }

  # set up vectors of model family
  # look for convention in the data of lik_fam_<<binom,gauss>>, if not there, default to config family
  lik_gaussian <- lik_binomial <- rep(0, nrow(d))
  if(('lik_fam_binom' %in% names(d)) & ('lik_fam_gauss' %in% names(d))) {
    lik_gaussian <- d$lik_fam_gauss
    lik_binomial <- d$lik_fam_binom
    message(sprintf('Found row specific data likelihood indicators, will use those. %i rows binom, %i rows gauss',
                    sum(lik_binomial),sum(lik_gaussian)))
  } else {
    if(indicator_family == 'binomial') {
      lik_binomial <- rep(1, nrow(d))
      message('Using indicator family binomial for all rows')
    } else if(indicator_family == 'gaussian') {
      lik_gaussian <- rep(1, nrow(d))
      message('Using indicator family gaussian for all rows')
    }
  }

  if(any(lik_gaussian+lik_binomial != 1))
    stop('Not all rows in your data have been assigned a model (binom or gauss), or some have been assigned multiple!')

  # also look for sd if already exists in the data for the gauss rows to use
  # This is useful for crosswalked values with some data uncertainty, convention is variable named sd_<<INDICATOR>>
  sd_i <- rep(0, nrow(d))
  if(paste0('sd_',indicator) %in%  names(d)){
    message('Found SD estimates to use for crosswalked values.')
    sd_i <- d[[paste0('sd_',indicator)]]
    sd_i[is.na(sd_i)] <- 0
  }


  # run a quick regression to get starting values for the fixed effects
  # This can speed up model fitting if iterations are slow.
  # Note this fit is currently assuming binomial
  if(all(lik_binomial==1)){
    message('LM for starting fe vals')
    y <- (d[[indic]][lik_binomial==1]+.0001)/d$N[lik_binomial==1]
    y[y<=0] <- 0.001
    y[y>=1] <- 0.999
    fe_start <- round( unname( lm(qlogis(y) ~ -1 + X_xp[lik_binomial==1,])$coefficients ), 4)
  } else {
    message('Default starting fe vals')
    fe_start <- rep(0,ncol(X_xp))
  }
  message(sprintf('starting values for fixed effects: %s',paste0(fe_start,collapse=', ')))


  # cannot allow a gaussian likelihood to also have a nugget in the linear term, it leads to issues
  if(nugget == 1 & all(lik_gaussian == 1)){
    message('WARNING:: Nugget in all gaussian model leads to identifiability issues. Removing nugget for you.')
    nugget <- 0
  }


  # check if user wants to scale gaussian variance by N, if not set them all to one in the gaussian rows
  n_i <- d$N
  if(scale_gaussian_variance_N == FALSE & sum(lik_gaussian)>1) {
    message('Not scaling gaussian error by N since scale_gaussian_variance_N == FALSE.')
    n_i[lik_gaussian==1] <- 1
  }

  # if there is only one country in the region, turn off country_re
  if(all(country_re == TRUE & length(get_adm0_codes(reg)) == 1)){
    message('WARNING: One country in this region, so turning off country random effects')
    country_re <- FALSE
  }

  # print some messages for random effects
  if(nugget == TRUE)     message('USING NUGGET (INDIVIDUAL OBSERVATION) RANDOM EFFECTS')
  if(country_re == TRUE) message('USING COUNTRY RANDOM EFFECTS')
  if(nid_re == TRUE)     message('USING NID RANDOM EFFECTS')


  # Build SPDE object (using INLA functions) and get prior in TMB readable format
  spde_list <- read_inla_prior_matern(spde_prior, mesh)
  spde <- spde_list$spde


  # Construct a list of all data necessary to TMB to fit
  Data <- list(
    num_i = nrow(d),       # Total number of observations
    num_s = mesh$n,        # Number of vertices in SPDE mesh
    num_t = length(yl),    # Number of periods
    num_z = num_z,         # 3rd dimension for GP,
    y_i = d[[indic]],      # Number of observed events in the cluster (N+ in binomial likelihood)
    n_i = d$N,             # Number of observed exposures in the cluster (N in binomial likelihood)
    t_i = d$period-1,      # Sample period ( starting at zero because C)
    c_re_i = cntry_re_vec, # vector of country ids, ( starting at zero because C)
    nid_re_i = nid_re_vec, # vector of survey ids, zero index added in cpp for null effect
    w_i = d$weight,        # Data weight for each row
    X_ij = X_xp,           # Covariate design matrix
    M0 = spde$param.inla$M0, # SPDE sparse matrix
    M1 = spde$param.inla$M1, # SPDE sparse matrix
    M2 = spde$param.inla$M2, # SPDE sparse matrix
    Aproj = A.proj,        # mesh to prediction point projection matrix
    lik_gaussian_i = lik_gaussian, # data likelihood for each row
    lik_binomial_i = lik_binomial, # data likelihood for each row
    sd_i           = sd_i, # crossalked standard deviation
    options = list(
      use_priors = 1,      # option1==1 use priors
      adreport_off = 1,    # option2==1 ADREPORT off
      nugget = nugget,     # option3==1 include nugget
      country_random = country_re, # option4==1 country random effects
      NID_random = nid_re, # option5==1 NID random effects
      useGP = as.numeric(as.logical(use_gp))
    ),
    prior_log_nugget_sigma = read_inla_prior_sigma(nugget_prior),
    prior_log_cre_sigma = read_inla_prior_sigma(ctry_re_prior),
    prior_log_nidre_sigma = read_inla_prior_sigma(nid_re_prior),
    prior_matern = spde_list$prior,
    fconstraints = tmb_cov_constraint(colnames(X_xp), cov_constraints)
  )

  # Set staring values for parameters
  Parameters <- list(alpha_j          = fe_start,  # FE parameters alphas
                     logtau           = spde$param.inla$theta.mu[1], # Matern/AR tau
                     logkappa         = spde$param.inla$theta.mu[2], # Matern Range
                     trho             = 0.95,                          # temporal rho
                     zrho             = 0.95,                          # 3rd dimension of GP rho (TODO)
                     log_nugget_sigma = -1,                            # log(SD) of the normal nugget term
                     log_cre_sigma    = -1,                            # log(SD) of the normal country intercept term (later add slopes as vector)
                     log_nidre_sigma  = -1,                            # log(SD) of the normal NID intercept term
                     log_gauss_sigma  = -1,                            # log(SD) of normal model
                     Epsilon_stz      = Epsilon_stz,                   # Random Effects: GP locations
                     nug_i            = rep(0,nrow(d)),                # Random Effects: Nugget Values
                     cntry_re         = rep(0,nrow(cntry_re_map)),     # Random Effects Values of country random effects (later add in slope stuff)
                     nid_re           = rep(0,max(nid_re_vec)+1))      # Random Effects Values of nid random effects (later add in slope stuff)


  # put bounds on parameters (Note, this is especially important for rhos)
  L  <- c(rep(-10,ncol(X_xp)),-10,-10,-99,-99,-10,-10,-10) ## updated the rho limits from .99999 since I transformed them in the cpp
  U  <- c(rep( 10,ncol(X_xp)), 10, 10, 99, 99, 10, 10, 10)
  pn <- c(rep('alpha_j',ncol(X_xp)),'logtau','logkappa','trho','zrho','log_nugget_sigma','log_cre_sigma','log_nidre_sigma')
  names(L) <- names(U) <- pn

  # return the list
  return(list(Data         = Data,
              Parameters   = Parameters,
              cs_df        = cs_df,
              clamper      = clamper,
              coords       = coords,
              mesh         = mesh,
              cntry_re_map = cntry_re_map,
              L            = L,
              U            = U))

}



#' @title Fit a TMB MBG model
#' @description Fit a TMB MBG model and pull jointPrecision report
#'
#' @param lbdcorerepo core repo location
#' @param cpp_template name of cpp template file within ./<lbdcorerepo>/mbg_central/
#' @param tmb_input_stack object that results from build_mbg_data_stack_tmb() or build_mbg_data_stack(...,tmb=TRUE)
#' @param ADmap_list map parameter for ignoring parameters
#' @param control_list pass control list to nlminb()
#' @param optimizer which software to use for optimization (optim, or nlminb)
#' @param sparse_ordering boolean: should the ADfun be adjusted before fitting
#'   so that the output SDreport object has a sparse ordering? This option
#'   requires the metis install for TMB.
#'
#' @return list of tmb model objects: ADfun, opt, sdrep
#'
#' @useDynLib mbg_tmb_model
#'
#' @note TODO support for other data models, sum to one constraint, country random effects
fit_mbg_tmb <- function(lbdcorerepo     = core_repo,
                        cpp_template    = 'mbg_tmb_model',
                        tmb_input_stack = input_data,
                        ADmap_list      = NULL,
                        control_list    = NULL,
                        optimizer       = 'nlminb',
                        sparse_ordering  = TRUE
                        ){

  message('WARNING: This TMB implementation does not support sum to one constraints yet.')

  # compile the cpp file and dynload it
  message('compiling template')
  TMB::compile(sprintf('%s/mbg_central/%s.cpp', lbdcorerepo, cpp_template))
  dyn.load( TMB::dynlib(sprintf('%s/mbg_central/%s', lbdcorerepo, cpp_template)) )

  # deal with parallelization
  threads <- system('echo $OMP_NUM_THREADS', intern = TRUE)
  if(threads != '') {
    message(sprintf('Detected %s threads in OMP environmental variable.',threads))
    TMB::openmp(as.numeric(threads))
  } else {
    message('Did not detect environmental OMP variable, defaulting to 4 cores. \n
             You can set this using OMP_NUM_THREADS or when launching singularity image.')
    TMB::openmp(4)
  }

  # set Data flag for Kaspers normalization fix
  tmb_input_stack$Data$flag <- 1
  # Initialize random effects
  randompars <- c()
  # Initialize Autodiff function map list (used to fix parameters)
  if(is.null(ADmap_list)) ADmap_list <- list()

  # if no z-col then set that ar1 to null
  if(length(dim(tmb_input_stack$Parameters$Epsilon_stz))==2){
    ADmap_list[['zrho']] <- factor(NA)
  }

  if(tmb_input_stack$Data$options$useGP==1){
    # if use_gp is on (options$useGP==1), Epsilon_stz should be random
    randompars <- c(randompars,"Epsilon_stz")
  } else {
    # if use_gp is off, map GP parameters out
    for(par in c('logtau','logkappa','trho','zrho')) ADmap_list[[par]] <- factor(NA)
    ADmap_list[['Epsilon_stz']] <- rep(factor(NA), length(tmb_input_stack$Parameters$Epsilon_stz))
  }

  if(tmb_input_stack$Data$options$nugget==1){
    # if nugget option is on add nug_i to randompars
    randompars <- c(randompars,'nug_i')
  } else {
    # if no nugget, add nug related parameters to ignorelist
    ADmap_list[['log_nugget_sigma']] <- factor(NA)
    ADmap_list[['nug_i']] <- rep(factor(NA),length(tmb_input_stack$Parameters$nug_i))
  }

  if(tmb_input_stack$Data$options$country_random==1) {
    # if country RE option is on add cntry_re to randompars
    randompars <- c(randompars,'cntry_re')
  } else {
    # if no country re, add cntry_re related parameters to ignorelist
    ADmap_list[['log_cre_sigma']] <- factor(NA)
    ADmap_list[['cntry_re']] <- rep(factor(NA),length(tmb_input_stack$Parameters$cntry_re))
  }

  if(tmb_input_stack$Data$options$NID_random==1) {
    # if NID RE option is on add nid_re to randompars
    randompars <- c(randompars,'nid_re')
  } else {
    # if no NID RE, add nid_re related parameters to ignorelist
    ADmap_list[['log_nidre_sigma']] <- factor(NA)
    ADmap_list[['nid_re']] <- rep(factor(NA),length(tmb_input_stack$Parameters$nid_re))
  }

  if(sum(tmb_input_stack$Data$lik_gaussian_i) == 0){
    # map out model sigma if no gaussian observations
    ADmap_list[['log_gauss_sigma']]    <- factor(NA)
  }

  # Print fixed parameters and random parameters to be used in TMB run
  message(paste0('ADMAP_LIST: ',  paste0(names(ADmap_list),collapse=', ')))
  message(paste0('Random Pars: ', paste0(randompars,collapse=', ')))


  # make the AD object
  message('Making AD object')
  obj <- TMB::MakeADFun(
    data       = tmb_input_stack$Data,
    parameters = tmb_input_stack$Parameters,
    map        = ADmap_list,
    random     = randompars,
    hessian    = TRUE,
    DLL        = cpp_template
  )

  # normalize
  obj <- TMB::normalize(obj, flag = "flag")

  # Reduce fill in of sparse Cholesky factorization (requires metis install of TMB)
  if(sparse_ordering){
    TMB::runSymbolicAnalysis(obj)
  }

  # Run optimizer
  message('Running MLE')
  if(optimizer == 'nlminb')
    opt0 <- do.call("nlminb",list(start       =    obj$par,
                                  objective   =    obj$fn,
                                  gradient    =    obj$gr,
                                  lower       =    tmb_input_stack$L,
                                  upper       =    tmb_input_stack$U,
                                  control     =    control_list))
  if(optimizer == 'optim')
    opt0 <- do.call("optim",list(par = obj$par, fn = obj$fn, control = control_list, gr = obj$gr, method = 'BFGS'))

  # run sdreport to get joint precision of all parameters
  for(i in 1:20) message('Getting Joint Precision')
  SD0 <- TMB::sdreport(obj, getJointPrecision = TRUE, bias.correct = TRUE)

  # return
  return(list(
    ADfun   = obj,
    opt     = opt0,
    sdrep   = SD0,
    fenames = colnames(tmb_input_stack$Data$X_ij)
  ))

}





#' Take multivariate normal draws given a mean vector and precision matrix
#'
#' @param mu vector of parameter means
#' @param prec joint precision matrix
#' @param n.sims number of draws
#'
#' @return length(mu) by n.sims matrix of parameter draws
#'
rmvnorm_prec <- function(mu, prec, n.sims) {
  z <- matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L <- Cholesky(prec, super = TRUE)
  z <- solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
  z <- solve(L, z, system = "Pt") ## z = Pt    %*% z
  z <- as.matrix(z)
  return(mu + z)
}









#' @title Predict MBG from a TMB Model
#'
#' @description Project out to a full sample space defined by the sr argument
#'
#' @param samples Number of draws to take
#' @param seed Seed to set for RNG
#' @param model_fit_object object output from function fit_mbg_tmb()
#' @param tmb_input_data input_data output of build_mbg_data_stack_tmb
#' @param fes a string of model fixed effects of form 'var1 + var2 + var3' or as string vector
#'  corresponding to column names in d
#' @param sr simple raster object
#' @param yl a vector of years for analysis (i.e. c(2001,2002,2003))
#' @param zl a vector of z for analysis (i.e. c(1,2,3,4,5)). Config item z_list. No z-dimension it should just be zero
#' @param covs_list named list of covariate bricks, each should be length(yl) long
#' @param clamp_covs Should covariates be 'clamped' to not predict outside of
#'   their observed range in the data?
#' @param cov_constraints named int vector. integer vector indexed by covariate
#'   names, in the format returned by the function
#'   \code{\link{covariate_constraint_vectorize}}. NOTE: In the current
#'   implementation, priors for fixed effects are set to the same values
#'   regardless of whether the covariate fixed effect is constrained or not.
#'   Apply constraints with caution.
#'
#' @return a cell_preds object
#'
predict_mbg_tmb <- function(samples,
                            seed             = NULL,
                            tmb_input_stack  = input_data,
                            model_fit_object = model_fit,
                            fes              = all_fixed_effects,
                            sr               = simple_raster,
                            yl               = year_list,
                            zl               = z_list,
                            use_space_only_gp    = as.logical(use_space_only_gp),
                            transform        = 'inverse-logit',
                            covs_list        = cov_list,
                            clamp_covs       = FALSE,
                            cov_constraints = covariate_constraint_vectorize(config)) {

  # Pull SD report object
  sdrep     <- model_fit_object$sdrep

  # pull a few useful things from the input data stack
  cs_transform         <- tmb_input_stack$cs_df
  mesh                 <- tmb_input_stack$mesh
  cntry_re_map         <- tmb_input_stack$cntry_re_map
  if(clamp_covs == TRUE) {
    clamper            <- tmb_input_stack$clamper
  } else {
    clamper            <- NULL
  }

  # set seed if it is requested
  if(!is.null(seed)) set.seed(seed)

  # vector of means
  mu    <- c(sdrep$par.fixed,sdrep$par.random)

  # simulate draws
  if(use_gp == TRUE){
    draws <- rmvnorm_prec(mu = mu , prec = sdrep$jointPrecision, n.sims = samples)

    ## separate out the draws
    parnames      <- c(names(sdrep$par.fixed), names(sdrep$par.random))
    epsilon_draws <- draws[parnames=='Epsilon_stz',]
    alpha_draws   <- draws[parnames=='alpha_j',]
  } else {

    non_eps_idx <- c(which(names(mu)!='Epsilon_stz'))

    draws <- rmvnorm_prec(mu = mu[non_eps_idx] , prec = sdrep$jointPrecision[non_eps_idx,non_eps_idx], n.sims = samples)

    ## separate out the draws
    parnames      <- names(mu)[non_eps_idx]
    epsilon_draws <- matrix(0, nrow=sum(names(mu)=='Epsilon_stz'),ncol=samples)
    alpha_draws   <- draws[parnames=='alpha_j',]
  }

  # seperate out Z FE draws
  FE_z_draws    <- NULL
  if(length(zl) > 1){
    FE_z_draws    <- alpha_draws[ grepl('FE_z_level__',model_fit_object$fenames),] # separate out z-level fixed effects from other FEs
    alpha_draws   <- alpha_draws[!grepl('FE_z_level__',model_fit_object$fenames),] # remove any z-level fixed effects from other FEs
  }


  if(length(zl) > 1)
    if(dim(FE_z_draws)[1] != (length(zl)-1) )
      stop('Incorrect number of fixed effects for levels of z in the z_list')

  # names of fes
  tmb_const <- tmb_cov_constraint(model_fit_object$fenames, cov_constraints)
  fes       <- unlist(strsplit(fes, ' \\+ '))

  # mask covariates to simple raster
  for(l in 1:length(covs_list)) {
    covs_list[[l]]  <- crop(covs_list[[l]], extent(sr))
    covs_list[[l]]  <- setExtent(covs_list[[l]], sr)
    covs_list[[l]]  <- mask(covs_list[[l]], sr)
  }

  # keep only covariates used in model, typically stacking
  covs_list <- covs_list[names(covs_list) %in% fes]

  # get coordinates of full projection space
  # Extract admin0 code
  f_orig <- data.table(cbind(xyFromCell(sr, seegSDM:::notMissingIdx(sr)), adm_code=as.vector(sr[seegSDM:::notMissingIdx(sr)])))
  f_orig$t <- f_orig$z <- 1 # set initial time and Z
  f_orig[,tmpord:=1:.N]

  # use the country code dt from input_data to map admin0 code to RE values
  f_orig <- merge(f_orig,cntry_re_map[,c('adm_code','re_id'),with=FALSE],by='adm_code',all.x=TRUE)
  f_orig <- f_orig[order(tmpord)] # make 100% sure everything is correctly ordered after the merge.
  f_orig[, re_id := re_id+1 ]  # to deal with indexing which started at 0 in the cpp
  f_orig$re_id[is.na(f_orig$re_id)] <- 0 # to deal with countries not in the data

  # add time periods and z periods as needed
  grp <- setDT(expand.grid(1:length(yl), 1:length(zl)))
  setnames(grp,c('Var1','Var2'),c('t','z'))
  grp[,group := 1:.N]
  fullsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,t  := grp$t[grp$group==g]]
    tmp[,gp := g]
    fullsamplespace <- rbind(fullsamplespace,tmp)
  }
  fullsamplespace[,idx := 1:.N]

  # pull out covariates in format we expect them
  # a list of length periods with a brick of named covariates inside
  new_cl <- list()
  if(length(covs_list)==0){
    message('No covariates detected, predicting using intercept only.')
  } else {
    message(sprintf('%i covariates detected.',length(covs_list)))
    for(p in 1:length(yl)){
      new_cl[[p]] <- list()
      for(n in names(covs_list)){
        if(dim(covs_list[[n]])[3]==1) { # synoptic mean covariates
          new_cl[[p]][[n]] <- covs_list[[n]]
        } else if (dim(covs_list[[n]])[3]==length(yl)) { # time varying covariates
          new_cl[[p]][[n]] <- covs_list[[n]][[p]]
        } else { # error if there is some other weird non-conforming year thing
          stop(sprintf('Covariate %n is a brick with %i layers, while year_list has %i years',
                       n,dim(covs_list[[n]])[3],length(yl)))
        }
      }
      new_cl[[p]] <- brick(new_cl[[p]])
    }
  }

  # get surface locs to project on to
  pcoords        <- cbind(x=fullsamplespace$x, y=fullsamplespace$y) ## used for cov raster extract

  ## setup coords for GP projection. convert coords to spherical if
  ## using spherical modeling mesh. used if you made a mesh on s2
  if(mesh_s$manifold == "S2"){
    gp_coords <- lonlat3D(pcoords[, 1], pcoords[, 2])
  } else {
    gp_coords <- pcoords
  }

  ## define grouping across periods
  groups_periods <- fullsamplespace$gp

  # extract cell values  from covariates, deal with timevarying covariates here
  cov_vals <- list()
  for(z in 1:length(zl)){
    cov_vals[[z]] <- list()
    for(p in 1:length(yl)){
      if(length(fes)>0) {

        # raster extract and keep only fes
        cov_vals[[z]][[p]] <- raster::extract(new_cl[[p]], pcoords[1:nrow(f_orig),])
        cov_vals[[z]][[p]] <- cov_vals[[z]][[p]][,colnames(cov_vals[[z]][[p]]) %in% c(fes)]
        # If there is only a single covariate, convert from vector to matrix
        if( (length(covs_list)==1) & !('matrix' %in% class(cov_vals[[z]][[p]]))){
          cov_vals[[z]][[p]] <- matrix(cov_vals[[z]][[p]], ncol=1)
        }

        # transform raw covariate values (center scaled) if needed (i.e. if cs_tranform is not 1 0 for that variable)
        cov_vals[[z]][[p]] <- centreScale(cov_vals[[z]][[p]],cs_transform)

        # clamp covariates if clamper is not null
        if(!is.null(clamper)){
          for(fe in fes){
            tmpvec <- cov_vals[[z]][[p]][,colnames(cov_vals[[z]][[p]])==fe]
            mn <- as.numeric(clamper[,fe,with=FALSE][1])
            mx <- as.numeric(clamper[,fe,with=FALSE][2])
            tmpvec[tmpvec<mn] <- mn
            tmpvec[tmpvec>mx] <- mx
            cov_vals[[z]][[p]][,colnames(cov_vals[[z]][[p]])==fe] <- tmpvec
          }
        }
        # add an intercept
        cov_vals[[z]][[p]] <- cbind(int = 1, cov_vals[[z]][[p]])
      } else {
        # if no covariates just do intercept only
        cov_vals[[z]][[p]] <- cbind(int = rep(1,nrow(f_orig)))
      }
      # if there is a z column, add on those fixed effects indicators
      if(length(zl) > 1){
        tmpzmat <- matrix(0,ncol = (length(zl)-1), nrow = nrow(f_orig))
        colnames(tmpzmat) <- paste0('FE_z_level__',2:length(zl))
        for(zz in 2:length(zl))
          if(z == zz)
            tmpzmat[,paste0('FE_z_level__',zz)] <- 1
        cov_vals[[z]][[p]] <- cbind(cov_vals[[z]][[p]], tmpzmat)
      }
    }
  }

  ## use inla helper functions to project the spatial effect.
  A.pred <- inla.spde.make.A(
    mesh  = mesh,
    loc   = gp_coords,
    group = groups_periods)

  ### values of GP ST surface at each cell (long by nperiods)
  # if we have multiple zs then do this by z since its possible to throw a SuiteSparse 'Problem too large' error here.
  if(length(zl) > 1){
    cell_s <- list()
    for(zz in zl){
      cell_s[[zz]] <- as.matrix(A.pred[(which(fullsamplespace$z==zz)),] %*% epsilon_draws)
    }
    cell_s <- do.call('rbind',cell_s)
  } else{
    cell_s <- as.matrix(A.pred %*% epsilon_draws)
  }

  # covariate values by alpha draws
  l_vals <- list()
  for(z in 1:length(zl)){
    l_vals[[z]] <- list()
    for(p in 1:length(yl))
      l_vals[[z]][[p]] <- cov_vals[[z]][[p]] %*% apply_constraints(tmb_const, rbind(alpha_draws,FE_z_draws))
  }
  cell_l <- do.call("rbind",unlist(l_vals, recursive = FALSE))


  # add the nugget if needed (i.e. nugget parameter was not mapped out)
  cell_nug <- matrix(0L, nrow = dim(cell_l)[1], ncol = dim(cell_l)[2])
  if('log_nugget_sigma' %in% parnames){
    if(exists('no_nugget_predict') & no_nugget_predict==TRUE){
      message('Nugget not included in predict')
    } else {
      message('adding nugget')
      for(s in 1:samples)
        cell_nug[,s] <- rnorm(dim(cell_nug)[1],0,exp(draws[parnames=='log_nugget_sigma',])[s])
    }
  }

  # add the country random intercept if it was estimated in the model
  cell_cre_int <- matrix(0L, nrow = dim(cell_l)[1], ncol = dim(cell_l)[2])
  if('log_cre_sigma' %in% parnames){
    message('adding country random intercept')
    cre <- data.table(draws[parnames=='cntry_re',])
    cre[, re_id := 1:.N]
    cre <- merge(fullsamplespace,cre,by='re_id',all.x=TRUE)
    cre <- cre[order(idx)]
    cre <- cre[,grep('V',colnames(cre)),with=FALSE]
    if(all(dim(cell_cre_int) == dim(cre))){
      cell_cre_int <- as.matrix(cre)
      rm(cre)
    } else {
      stop('CHECK COUNTRY RE DIMENSIONS')
    }

  }

  # add together linear and st components
  pred_tmb <- cell_l + cell_s + cell_nug + cell_cre_int

  # transform
  if(transform=='inverse-logit') {
    pred_tmb <- plogis(as.matrix(pred_tmb))
  } else {
    pred_tmb <- eval(parse(text=sprintf('%s(as.matrix(pred_tmb))',transform)))
  }

  # if there is more than one z, then return a list of length zl cell_preds
  if(length(zl) > 1){
    pred_tmb_list <- list()
    chunklength <- dim(pred_tmb)[1]/length(zl)
    for(z in 1:length(zl))
      pred_tmb_list[[z]] <- pred_tmb[((z-1)*chunklength+1):(chunklength*z),1:samples]
    pred_tmb <- pred_tmb_list
  }

  # return the predicted cell_pred object
  return(pred_tmb)

}




#' @title Get fitted model parameters from TMB object
#' @description Make a nice table of fitted model parameters from a geostat TMB object
#'
#' @param model_fit fitted TMB object that comes out of fit_mbg_tmb()
#' @param exp_fixed_effects Boolean, should the fixed effects be exponentiated (if model was fit with logit link). Defaults to TRUE
#' @param transform_hyperparams Boolean, should hyperparmeters be transformed from fitted to more natural space. Defaults to TRUE
#' @param draws Integer, number of draws to use. Defaults to 1000
#' @param calculate_range Boolean, should we calculate and report range using kappa in draws space. Defaults to TRUE
#' @param calculate_nom_var  Boolean, should we calculate and report nominal variance using kappa and tau in draws space. Defaults to TRUE
#'
#' @return formatted data.table object
fitted_param_table_tmb <- function(model_fit,
                                   exp_fixed_effects     = TRUE,
                                   transform_hyperparams = TRUE,
                                   draws                 = 1000,
                                   calculate_range       = TRUE,
                                   calculate_nom_var     = TRUE,
                                   cov_constraints = covariate_constraint_vectorize(config)) {

  # get draws of parameter values
  mu <- model_fit$sdrep$par.fixed
  pm <- model_fit$sdrep$jointPrecision[1:length(mu),1:length(mu)]
  draws <- rmvnorm_prec(mu,pm,draws)

  # deal with given names of parameters
  fn <- names(model_fit$sdrep$par.fixed)
  fn[fn == 'alpha_j'] <- model_fit$fenames

  # Apply constraint transformations
  tmb_const <- tmb_cov_constraint(model_fit$fenames, cov_constraints)
  draws[names(model_fit$sdrep$par.fixed) == "alpha_j",] <-
    apply_constraints(tmb_const, draws[names(model_fit$sdrep$par.fixed) == "alpha_j",])

  # Transform fixed effects
  if(exp_fixed_effects == TRUE){
    draws[which(fn %in% model_fit$fenames),] <- exp(draws[which(fn %in% model_fit$fenames),])
  }

  # Get the range parameter
  if(calculate_range == TRUE){
    # see equation 6.16 (pg 196) of Blangiardo and Cameletti Book
    ranger <- sqrt(8) / exp(draws[which(fn == 'logkappa'),])
    draws  <- rbind(draws, ranger)
    fn     <- c(fn, 'range')
  }

  # Get the nominal variance parameter
  if(calculate_nom_var == TRUE){
    # see equation 6.17 (pg 196) of Blangiardo and Cameletti Book
    nominal_variance <- 1 / (4 * pi * (exp(draws[which(fn == 'logkappa'),]))^2 * (exp(draws[which(fn == 'logtau'),]))^2)
    draws <- rbind(draws, nominal_variance)
    fn    <- c(fn, 'nominal_variance')
  }

  # transform hyperparmeters
  if(transform_hyperparams == TRUE){
    draws[which(fn == 'logtau'),]             <- exp(draws[which(fn == 'logtau'),])
    draws[which(fn == 'logkappa'),]           <- exp(draws[which(fn == 'logkappa'),])
    draws[which(fn == 'log_nugget_sigma'),]   <- exp(draws[which(fn == 'log_nugget_sigma'),])
    draws[which(fn == 'log_cre_sigma'),]      <- exp(draws[which(fn == 'log_cre_sigma'),])
    draws[which(fn == 'log_nidre_sigma'),]    <- exp(draws[which(fn == 'log_nidre_sigma'),])

    fn[fn == 'logtau']           <- 'tau'
    fn[fn == 'logkappa']         <- 'kappa'
    fn[fn == 'log_nugget_sigma'] <- 'nugget_SD'
    fn[fn == 'log_cre_sigma']    <- 'country_RE_SD'
    fn[fn == 'log_nidre_sigma']  <- 'NID_RE_SD'

    draws[which(fn == 'zrho'),] <- (exp( draws[which(fn == 'zrho'),] ) - 1) / (exp( draws[which(fn == 'zrho'),] ) + 1)
    draws[which(fn == 'trho'),] <- (exp( draws[which(fn == 'trho'),] ) - 1) / (exp( draws[which(fn == 'trho'),] ) + 1)
    fn[fn == 'zrho'] <- 'age_rho'
    fn[fn == 'trho'] <- 'year_rho'

  }

  # summarize draws and clean up data table
  su <- data.table(t(apply(draws,1,quantile,c(0.025,0.500,0.975))))
  su[, fn := fn]
  colnames(su) <- c('lower','median','upper','param_name')
  su <- su[,c('param_name','median','lower','upper'),with=FALSE]

  # return the final table
  return(su)
}

#' @title Read INLA SD prior for TMB
#'
#' @description Read in a prior specification that is suited for INLA and make
#'   it TMB readable.
#'
#' @param prior_string character, character vec of length 1 specifying priors
#'
#' @return List specifying a TMB prior, containing three elements:
#'   - type: Is the prior normal, loggamma, or pc.prec
#'   - par1: The first shape parameter. In the lognormal case, the mean
#'   - par2: The second shape parameter. In the lognormal case, the variance
#'
read_inla_prior_sigma <- function(prior_string){
  prior_list <- eval(parse(text=prior_string[1]))
  if(!(prior_list$prior %in% c("normal", "loggamma", "pc.prec"))){
    stop("TMB implementation only supports normal, loggamma, or PC priors for
         SD parameters.")
  }
  return(list(
    type = prior_list$prior,
    par1 = prior_list$param[1],
    par2 = prior_list$param[2]
  ))
}

#' @title Read INLA Matern GP priors for TMB
#'
#' @description Read in a prior specification from config and make
#'   it TMB readable.
#'
#' @param prior_string character, character vec of length 1 specifying priors
#'
#' @return List containing (1) spde object and (2) list specifying a TMB prior,
#'   containing three elements:
#'   - type: Is the prior pc or nonpc (i.e. normal)
#'   - par1: Vector of length 2 for the first parameter. In the nonpc case,
#'     corresponds to mean and precision for logtau. In the pc case,
#'     corresponds to range0 and prange.
#'   - par2: Vector of length 2 for the first parameter. In the nonpc case,
#'     corresponds to mean and precision for logkappa. In the pc case,
#'     corresponds to sigma0 and psigma.
read_inla_prior_matern <- function(prior_string, mesh_s){
  prior_list <- eval(parse(text=prior_string[1]))

  spde_list <- build_spde_prior(prior_list, mesh_s, st_gp_int_zero=FALSE)
  spde_prior <- spde_list$spde_prior
  spde <- spde_list$spde

  if(spde_prior$type == "nonpc") {
    par1 <- c(spde$param.inla$theta.mu[1], spde$param.inla$theta.Q[1,1])
    par2 <- c(spde$param.inla$theta.mu[2], spde$param.inla$theta.Q[2,2])
  } else {
    par1 <- spde_prior$prior$range
    par2 <- spde_prior$prior$sigma
  }


  return(list(spde=spde,
              prior=list(
                type = prior_list$type,
                par1 = par1,
                par2 = par2)
  ))
}

#' @title Modify constraints to fit TMB formatting
#'
#' @description Optionally add constraints to fixed effects in the TMB
#'   optimization model. In the TMB model, a constraint label of '0' indicates
#'   that the fixed effect is unconstrained, a label of '1' constrains the
#'   above zero, and and a label of '-1' constrains the variable below zero.
#'
#' @param fes character, fixed effects in model
#' @param constraints named int vector output from the
#'   `covariate_constraint_vectorize()` function. Specifies how to constrain
#'   each fixed effect.
#' @param zl int, additional constraints to pad on which will be 0
#'
tmb_cov_constraint <- function(
  model_fes,
  constraints = covariate_constraint_vectorize(config)
  ){
  tmb_const <- sapply(unlist(strsplit(model_fes, " \\+ ")), function(m){
    if(m %in% names(constraints)){
      x <- unname(constraints[m])
    }
    else{
      x <- 0
    }
    x
  })

  return(tmb_const)
}

#' @title Apply constraints to fitted value draws]\
#'
#' @description Given draws of fixed effect values generated from a fitted TMB
#'   model and a list of constraints applied to fixed effects in that model,
#'   apply transformations on constrained fixed effects to reproduce how they
#'   were incorporated in the model (constraining them above or below zero). If
#'   a fixed effect had a constraint label of 0 (unconstrained, the default),
#'   the untransformed draws will be returned.
#'
#' @param tmb_const int vector, tmb prepped constraints
#' @param alpha_draws matrix, fitted draws of coefficients
#'
#' @return matrix of transformed beta coefficients for fixed effect draws
#'
apply_constraints <- function(tmb_const, FE_draws){
  Nfe <- nrow(FE_draws)

  FEdraws_const <- sapply(1:Nfe, function(j){
    X <- FE_draws[j,]
    if(tmb_const[j] == 1){
      X <- exp(X)
    }
    else if(tmb_const[j] == -1){
      X <- -exp(X)
    }

    return(X)
  })

  return(t(FEdraws_const))
}

