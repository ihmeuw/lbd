


#' Organize Data and Parameter Stacks to fit a TMB MBG model
#' Author: RB
#'
#' @param d prepped model frame with no NAs
#' @param yl a vector of years for analysis (i.e. c(2001,2002,2003))
#' @param zl a vector of zcol for analysis (i.e. ages c(1,2,3,4,5). Must be integers starting with 1)
#' @param fes a string of model fixed effects of form 'var1 + var2 + var3' or as string vector
#'  corresponding to column names in d
#' @param indic indicator name
#' @param country_re TRUE/FALSE include country re. If true and there is a zcol then there will be random slope as well
#' @param nid_re TRUE/FALSE include re on survey. Need nid column in the data frame. defaults to use_nid_res from config.
#' @param exclude_cs character vector of covariates to exclude from centrescaling
#' @param indic indicator name, corresponding to the response variable in d
#' @param mesh an inla mesh object. NOTE! this gets remade inside this function... Roy said there was a good reason
#' @param s2mesh Logical. Should the mesh be created on the surface of
#'   a sphere? If TRUE, s2params is used to specify mesh parameters
#'   instead of max_edge and mesh_offset
#' @param s2params string of 3 element numeric vector in R
#'   notation. e.g. "c(25, 500, 1000)". The entries describe the
#'   minimum triangle edge length allowed, hos far to extend the mesh
#'   beyond the 'simple' boundary, and the maximum allowed triangle
#'   edge length, respectively. Units are in kilometers. Used only if
#'   s2mesh=TRUE.
#'
#' @return returns a named list with Data and Parameters to be passed into fit_mbg_tmb()
#'
#'
#'
build_mbg_data_stack_tmb <- function(d          = df,
                                     yl         = year_list,
                                     zl         = z_list,
                                     fes        = all_fixed_effects,
                                     indic      = indicator,
                                     country_re = use_inla_country_res,
                                     nid_re     = use_nid_res,
                                     exclude_cs = '',
                                     nugget     = FALSE,
                                     zcol       = NULL,
                                     shapefile_version = 'current',
                                     scale_gaussian_variance_N = TRUE,
                                     mesh       = mesh_s
                                     ){

  #
  # require some libraries
  require(INLA)
  require(data.table)

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
  num_z <- 1 #
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

  # remake mesh ##
  mesh   <- build_space_mesh(d           = d,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset,
                             s2mesh      = use_s2_mesh,
                             s2params   = s2_mesh_params)

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

  # Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary) to get matrices
  spde <- inla.spde2.matern(mesh, alpha = 2 )

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
  #
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

  # make a gaul/country/cntry_RE_idx mapping table
  md     <- get_location_code_mapping(shapefile_version = shapefile_version)
  cntrys <- unique(as.character(d$country))
  gaulz  <- md$GAUL_CODE[md$ihme_lc_id %in% cntrys]
  if(length(cntrys) != length(gaulz))
    stop("When making country RE IDs could not get gauls from metadata of each country ISO. ")
  cntry_re_map <- data.table(
                    country   = cntrys,
                    gaul_code = gaulz,
                    re_id     = 0:(length(cntrys)-1))
  cntry_re_vec <- cntry_re_map$re_id[match(as.character(d$country),cntry_re_map$country)]

  # make an nid_re mapping table NID NID
  if(use_nid_res){
    nidz   <- unique(as.numeric(d$nid))
    nid_re_map <- data.table(
      nid    = nidz,
      re_id  = 0:(length(nidz)-1))
    nid_re_vec <- nid_re_map$re_id[match(as.numeric(d$nid),nid_re_map$nid)]
  }else{
    nidz <- 0
    nid_re_map <- data.table(
      nid    = nidz,
      re_id  = 0:(length(nidz)-1))
    nid_re_vec <- rep(0, nrow(df))

  }
  # sort country RE, for now just intercept. later add Random slope by zcol


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
  print(sum(lik_binomial))
  print(sum(lik_gaussian))
  if(sum(lik_binomial)>1){
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

  # print some messages for random effectss
  if(nugget == TRUE)     message('USING NUGGET (INDIVIDUAL OBSERVATION) RANDOM EFFECTS')
  if(country_re == TRUE) message('USING COUNTRY RANDOM EFFECTS')
  if(nid_re == TRUE)     message('USING NID RANDOM EFFECTS')

  # NOTE (RB 2AUG2018) I have concerns about identifiability in situations where there is one survey in a
  #  country and country and nid REs are both turned on
  # identify these obervations and pass them to model so that the model can ignore nid random effects in these situation
  #  another option could be to use the map to set these values to zero in this case. Waiting first to see if this actually causes convergence
  #  issues before dealing with it

  # Construct a list of all Data necessary to TMB to fit
  Data <- list(num_i          = nrow(d),               # Total number of observations
               num_s          = mesh$n,                # Number of vertices in SPDE mesh
               num_t          = length(yl),            # Number of periods
               num_z          = num_z,                 # 3rd dimension for GP,
               y_i            = d[[indic]],            # Number of observed events in the cluster (N+ in binomial likelihood)
               n_i            = d$N,                   # Number of observed exposures in the cluster (N in binomial likelihood)
               t_i            = d$period-1,            # Sample period ( starting at zero because C)
               c_re_i         = cntry_re_vec,          # vector of country ids, ( starting at zero because C)
               nid_re_i       = nid_re_vec,            # vector of survey ids, ( starting at zero because C)
               w_i            = d$weight,              # Data weight for each row
               X_ij           = X_xp,                  # Covariate design matrix
               M0             = spde$param.inla$M0,    # SPDE sparse matrix
               M1             = spde$param.inla$M1,    # SPDE sparse matrix
               M2             = spde$param.inla$M2,    # SPDE sparse matrix
               Aproj          = A.proj,                # mesh to prediction point projection matrix
               lik_gaussian_i = lik_gaussian,          # data likelihood for each row
               lik_binomial_i = lik_binomial,          # data likelihood for each row
               sd_i           = sd_i,                  # crossalked standard deviation
               options        = c(1,                   # option1==1 use priors
                                  1,                   # option2==1 ADREPORT off
                                  nugget,              # option3==1 include nugget
                                  country_re,          # option4==1 country random effects
                                  nid_re))             # option5==1 NID random effects
                                                       # option6==1 covariates sum to 1 constraint
                                                       # option7==1 no intercept

                                                #

  # Set staring values for parameters
  Parameters <- list(alpha_j          = fe_start,  # FE parameters alphas
                     logtau           = -0.5,                           # Matern/AR tau
                     logkappa         = -0.5,	                       # Matern Range
                     trho             = 0.95,                          # temporal rho
                     zrho             = 0.95,                          # 3rd dimension of GP rho 
                     log_nugget_sigma = -1,                            # log(SD) of the normal nugget term
                     log_cre_sigma    = -1,                            # log(SD) of the normal country intercept term (later add slopes as vector)
                     log_nidre_sigma  = -1,                            # log(SD) of the normal NID intercept term
                     log_gauss_sigma  = -1,                            # log(SD) of normal model
                     Epsilon_stz      = Epsilon_stz,                   # Random Effects: GP locations
                     nug_i            = rep(0,nrow(d)),                # Random Effects: Nugget Values
                     cntry_re         = rep(0,nrow(cntry_re_map)),     # Random Effects Values of country random effects (later add in slope stuff)
                     nid_re           = rep(0,nrow(nid_re_map)))     # Random Effects Values of country random effects (later add in slope stuff)


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
              nid_re_map   = nid_re_map,
              L            = L,
              U            = U))

}





#' Fit a TMB MBG model and pull jointPrecision report
#' Author: RB
#'
#' @param lbdcorerepo core repo location
#' @param cpp_template name of cpp template file within ./<lbdcorerepo>/mbg_central/
#' @param tmb_input_stack object that results from build_mbg_data_stack_tmb() or build_mbg_data_stack(...,tmb=TRUE)
#' @param ADmap_list map parameter for ignoring parameters
#' @param control_list pass control list to nlminb()
#' @param optimizer which software to use for optimization (optim, or nlminb)
#'
#' @return list of tmb model objects: ADfun, opt, sdrep
#'
#'
fit_mbg_tmb <- function(lbdcorerepo     = core_repo,
                        cpp_template    = 'mbg_tmb_model',
                        tmb_input_stack = input_data,
                        ADmap_list      = NULL,
                        control_list    = NULL,
                        optimizer       = 'nlminb'
                        ){

  message('WARNING: This TMB implementation does not support sum to one constraints yet.')

  # make sure TMB is loaded
  require(TMB)

  # compile the cpp file and dynload it
  message('compiling template')
  TMB::compile(sprintf('%s/mbg_central/%s.cpp', lbdcorerepo, cpp_template))
  dyn.load( dynlib(sprintf('%s/mbg_central/%s', lbdcorerepo, cpp_template)) )

  # deal with parallelization
  threads <- system('echo $OMP_NUM_THREADS', intern = TRUE)
  if(threads != '') {
    message(sprintf('Detected %s threads in OMP environmental variable.',threads))
    openmp(as.numeric(threads))
  } else {
    message('Did not detect environmental OMP variable, defaulting to 4 cores. \n
             You can set this using OMP_NUM_TREADS or when launching singularity image.')
    openmp(4)
  }

  # set Data flag for Kaspers normalization fix
  tmb_input_stack$Data$flag <- 1

  # if no z-col then set that ar1 to null
  if(length(dim(tmb_input_stack$Parameters$Epsilon_stz))==2){
    if(is.null(ADmap_list))
      ADmap_list <- list()
    ADmap_list[['zrho']] <- factor(NA)
  }

  # REs
  randompars <- c("Epsilon_stz")

  # if nugget option is on add nug_i to randompars
  if(tmb_input_stack$Data$options[3]==1){
    randompars <- c(randompars,'nug_i')
  } else { # if no nugget, add nug related parameters to ignorelist
    if(is.null(ADmap_list))
      ADmap_list <- list()
    ADmap_list[['log_nugget_sigma']] <- factor(NA)
    ADmap_list[['nug_i']]            <- rep(factor(NA),length(tmb_input_stack$Parameters$nug_i))
  }

  # if country RE option is on add cntry_re to randompars
  if(tmb_input_stack$Data$options[4]==1) {
    randompars <- c(randompars,'cntry_re')
  } else { # if no nugget, add nug related parameters to ignorelist
    if(is.null(ADmap_list))
      ADmap_list <- list()
    ADmap_list[['log_cre_sigma']]    <- factor(NA)
    ADmap_list[['cntry_re']]         <- rep(factor(NA),length(tmb_input_stack$Parameters$cntry_re))
  }

  # if NID RE option is on add nid_re to randompars
  if(tmb_input_stack$Data$options[5]==1) {
    randompars <- c(randompars,'nid_re')
  } else { # if no nugget, add nug related parameters to ignorelist
    if(is.null(ADmap_list))
      ADmap_list <- list()
    ADmap_list[['log_nidre_sigma']]    <- factor(NA)
    ADmap_list[['nid_re']]         <- rep(factor(NA),length(tmb_input_stack$Parameters$nid_re))
  }


  # map out model sigma if no gaussian observations
  if(sum(tmb_input_stack$Data$lik_gaussian_i) == 0){
    if(is.null(ADmap_list))
      ADmap_list <- list()
    ADmap_list[['log_gauss_sigma']]    <- factor(NA)
  }

  # make the AD object
  message('Making AD object')
  obj <- MakeADFun(data       = tmb_input_stack$Data,
                   parameters = tmb_input_stack$Parameters,
                   map        = ADmap_list,
                   random     = randompars,
                   hessian    = TRUE,
                   DLL        = cpp_template) #sprintf('%s/mbg_central/%s', lbdcorerepo, tmpl))

  # normalize
  obj <- normalize(obj, flag = "flag")

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
  return( list ( ADfun   = obj,
                 opt     = opt0,
                 sdrep   = SD0,
                 fenames = colnames(tmb_input_stack$Data$X_ij)) )

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









#' Predict MBG from a TMB Model
#'
#' Project out to a full sample space defined by the sr argument
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
                            transform        = 'inverse-logit',
                            covs_list        = cov_list,
                            clamp_covs       = FALSE) {

  # libs
  require(raster)
  require(sp)

  #
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
  draws <- rmvnorm_prec(mu = mu , prec = sdrep$jointPrecision, n.sims = samples)

  ## separate out the draws
  parnames      <- c(names(sdrep$par.fixed), names(sdrep$par.random))
  epsilon_draws <- draws[parnames=='Epsilon_stz',]
  alpha_draws   <- draws[parnames=='alpha_j',]
  FE_z_draws    <- NULL
  if(length(zl) > 1){
    FE_z_draws    <- alpha_draws[ grepl('FE_z_level__',model_fit_object$fenames),] # separate out z-level fixed effects from other FEs
    alpha_draws   <- alpha_draws[!grepl('FE_z_level__',model_fit_object$fenames),] # remove any z-level fixed effects from other FEs
  }

  if(length(zl) > 1)
    if(dim(FE_z_draws)[1] != (length(zl)-1) )
      stop('Incorrect number of fixed effects for levels of z in the z_list')

  # names of fes
  fes       <- unlist(strsplit(fes, ' \\+ '))

  # mask covariates to simple raster
  for(l in 1:length(covs_list)) {
    covs_list[[l]]  <- crop(covs_list[[l]], extent(sr))
    covs_list[[l]]  <- setExtent(covs_list[[l]], sr)
    covs_list[[l]]  <- mask(covs_list[[l]], sr)
  }

  # keep only covariates used in model, typically stacking
  covs_list <- covs_list[names(covs_list) %in% fes]
  if((length(covs_list)+1+(length(zl)-1)) != dim(alpha_draws)[1])
    stop('Check that all covariates are accounted for, number of alphas does not conform with number of covariates in covs_list that match fes.')


  # get coordinates of full projection space
  #f_orig <- data.table(cbind(sp::coordinates(sr)), t=1) # sold version that missed 4 pixels in india
  f_orig <- data.table(cbind(xyFromCell(sr, seegSDM:::notMissingIdx(sr)), gaul_code=as.vector(sr[seegSDM:::notMissingIdx(sr)]))) # extract gaul
  f_orig$t <- f_orig$z <- 1 # set initial time and Z
  f_orig[,tmpord:=1:.N]

  # use the gaul dt from input_data to map gual to RE values
  f_orig <- merge(f_orig,cntry_re_map[,c('gaul_code','re_id'),with=FALSE],by='gaul_code',all.x=TRUE)
  f_orig <- f_orig[order(tmpord)] # make 100% sure everything is correctly ordered after the merge.
  f_orig[, re_id := re_id+1 ]  # to deal with indexing which started at 0 in the cpp

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

        # transform raw covariate values (center scaled) if needed (i.e. if cs_tranform is not 1 0 for that variable)
        cov_vals[[z]][[p]] <- centreScale(cov_vals[[z]][[p]],cs_transform)

        # clamp covariates if clamper is not null
        if(!is.null(clamper)){
         # message('Clamping')
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
      l_vals[[z]][[p]] <- cov_vals[[z]][[p]] %*% rbind(alpha_draws,FE_z_draws)
  }
  cell_l <- do.call("rbind",unlist(l_vals, recursive = FALSE))


  # add the nugget if needed (i.e. nugget parameter was not mapped out)
  cell_nug <- matrix(0L, nrow = dim(cell_l)[1], ncol = dim(cell_l)[2])
  if('log_nugget_sigma' %in% parnames){
    message('adding nugget')
    for(s in 1:samples)
      cell_nug[,s] <- rnorm(dim(cell_nug)[1],0,exp(draws[parnames=='log_nugget_sigma',])[s])
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
  pred_tmb_no_gp <- cell_l + cell_nug + cell_cre_int

  # transform
  if(transform=='inverse-logit') {
    pred_tmb <- plogis(as.matrix(pred_tmb))
    pred_tmb_no_gp <- plogis(as.matrix(pred_tmb_no_gp))
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




#' Make a nice table of fitted model parameters from a geostat TMB object
#' RB
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
                                   calculate_nom_var     = TRUE) {

  # get draws of parameter values
  mu <- model_fit$sdrep$par.fixed
  pm <- model_fit$sdrep$jointPrecision[1:length(mu),1:length(mu)]
  draws <- rmvnorm_prec(mu,pm,draws)

  # deal with given names of parameters
  fn <- names(model_fit$sdrep$par.fixed)
  fn[fn == 'alpha_j'] <- model_fit$fenames

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
