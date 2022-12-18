


#' @title Build Stack for TMB MBG model
#' @description Organize Data and Parameter Stacks to fit a TMB MBG model
#'
#' @param d prepped model frame with no NAs
#' @param yl a vector of years for analysis (i.e. c(2001,2002,2003))
#' @param zl a vector of zcol for analysis (i.e. ages c(1,2,3,4,5). Must be integers starting with 1)
#' @param fes a string of model fixed effects of form 'var1 + var2 + var3' or as string vector
#'  corresponding to column names in d
#' @param indic indicator name, corresponding to the response variable in d
#' @param country_re TRUE/FALSE include country re. If true and there is a zcol then there will be random slope as well
#' @param nid_re TRUE/FALSE include re on survey. Need nid column in the data frame. defaults to use_nid_res from config. 
#' @param exclude_cs character vector of covariates to exclude from centrescaling
#' @param mesh_s an inla mesh object 
#' @param nugget Logical. If TRUE, include a nugget effect
#' @param ycol column name of y values, defaults to year
#' @param zcol column name of z values associated with zl.
#' @param shapefile_version character. Version of shape file to use.
#' @param scale_gaussian_variance_N Logical. Do you want to scale gaussian variance by sample size?
#' @param mesh_int an inla mesh object.
#' @param cov_constraints named int vector. integer vector indexed by covariate 
#'   names, in the format returned by the function 
#'   \code{\link{covariate_constraint_vectorize}}. NOTE: In the current
#'   implementation, priors for fixed effects are set to the same values
#'   regardless of whether the covariate fixed effect is constrained or not. 
#'   Apply constraints with caution.
#' @param main_space_effect TRUE/FALSE to include a space-only GP
#' @param main_time_effect TRUE/FALSE to include a time-only GMRF
#' @param main_age_effect TRUE/FALSE to include an age-only GMRF
#' @param use_full_interacting_effect TRUE/FALSE to use an interacting GP - this will be space-time if no z-dimension is present, space-age if only one yar is presented, and space-time-age if all three dimensions exist
#' @param coefs.sum1 Logical. Were the coefficients constrained to sum-to-1 in the model fit
#' @param stacker_names string vector of child model names

#'
#' @return returns a named list with Data and Parameters to be passed into fit_mbg_tmb()
#' 
#' 
#' @export
build_mbg_data_stack_tmb <- function(d          = df,                 
                                     yl         = year_list,  
                                     zl         = z_list,
                                     xl         = sex_list,
                                     fes        = all_fixed_effects,  
                                     indic      = indicator, 
                                     country_re = use_country_res, 
                                     nid_re     = use_nid_res,
                                     exclude_cs = '', 
                                     nugget     = FALSE,
                                     zcol       = NULL,
                                     ycol       = "year",
                                     xcol       = NULL,
                                     shapefile_version = 'current', 
                                     scale_gaussian_variance_N = TRUE,
                                     mesh_int       = mesh_int,
                                     mesh_s         = mesh_s,
                                     cov_constraints = covariate_constraint_vectorize(config),
                                     main_space_effect = as.logical(use_space_only_gp),
                                     main_time_effect = as.logical(use_time_only_gmrf),
                                     main_age_effect = as.logical(use_age_only_gmrf),
                                     main_sex_effect = as.logical(use_sex_only_gmrf),
                                     int_gp_1_effect = as.logical(use_gp),
                                     int_gp_2_effect = as.logical(any(interacting_gp_2_effects != '')),
                                     use_covs   = any(as.logical(use_stacking_covs), as.logical(use_raw_covs)),
                                     use_sz_gp = as.logical(use_sz_gp),
                                     use_sx_gp = as.logical(use_sx_gp),
                                     use_tx_gp = as.logical(use_tx_gp),
                                     use_zx_gp = as.logical(use_zx_gp),
                                     use_cre_z_gp  = as.logical(use_cre_z_gp),
                                     use_cre_x_gp  = as.logical(use_cre_x_gp),
                                     use_anc_corrections = as.logical(use_anc_corrections),
                                     use_error_iid_re = as.logical(use_error_iid_re),
                                     use_observation_level_error_iid = as.logical(use_observation_level_error_iid),
                                     use_cyzx_error_iid = as.logical(use_cyzx_error_iid),
                                     use_z_fes  = FALSE,
                                     coefs.sum1 = coefs_sum1,
                                     stacker_names = child_model_names){   

  
  # ensure d is a dt
  d <- setDT(d)
  
  # order d by row_id
  d <- d[order(row_id), ]
  
  # fix row_id to increment by 1
  d <- unique(select(d, row_id)) %>% mutate(row_id_transf = 1:n()) %>% 
    right_join(d, by = "row_id")
  d <- setDT(d)
  
  # zcol
  if (!zcol %in% colnames(d)){
    message('No Z column detected')
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
  
  # xcol
  if(!xcol %in% colnames(d)){
    message('No Sex column detected')
    d[[xcol]] <- 0
  } else if (is.null(xcol)) {
    message('Sex column was set as null')
    d[[xcol]] <- 0
  }
  if( !all(unique(d[[xcol]]) %in% xl)) {
    message('WARNING: xl and d[[sexcol]] do not completely match up.. ')
    # check that we there arent values in zcol not matching z_list
    d[, dropx := !get(xcol) %in% xl]
    if(any(d$dropx != FALSE)){
      message(sprintf('WARNING: Detected some sex values in xcol (%s) which were not in the sex_list ',xcol))
      message(sprintf('WARNING: Due to this, dropping %i rows from the input data',sum(d$dropx==TRUE)))
      print(table(d$dropx,d[[xcol]]))
      d <- subset(d, dropx == FALSE)
      d[, dropx := NULL]
    }
  }
  
  # make a fake data point with no weight for the max period and Z to fill out the GP
  d <- rbind(d, d[1,])
  d[[zcol]][nrow(d)]     <- max(zl)
  d[[xcol]][nrow(d)] <- length(xl)
  d[['period']][nrow(d)] <- length(yl)
  d[['weight']][nrow(d)] <- 0
  d[['row_id_transf']][nrow(d)] <- max(d$row_id_transf) + 1
  
  message(paste0('interacting gp 1 effects: ', interacting_gp_1_effects))
  message(paste0('interacting gp 2 effects: ', interacting_gp_2_effects))
  
  # look for z & x dimensions 
  num_z <- length(zl)
  num_x <- length(xl)
  
  #Set up groups for interaction # 1
  if(length(zl)>1 & 'age' %in% interacting_gp_1_effects){ 
    message(sprintf('More than one unique %s found, initiating Z in GP 1', zcol))
    # look for sex dimension 
    if(length(xl)==1 | !('sex' %in% interacting_gp_1_effects)){ 
      # set A proj grouping. The ordering here must match the ordering of epsilon_int_1 in the template
      grp <- setDT(expand.grid(1:length(yl), 1:max(zl)))
      setnames(grp,c('Var1','Var2'),c('period',zcol))
      grp[,group_1 := 1:.N]
      d <- merge(d, grp, by = c('period',zcol), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
    } else if (length(xl)==2 & 'sex' %in% interacting_gp_1_effects){
      message(sprintf('More than one unique %s found, initiating Sex in GP 1', xcol))
      grp <- setDT(expand.grid(1:length(yl), 1:max(zl), 1:max(xl)))
      setnames(grp,c('Var1','Var2','Var3'),c('period',zcol, xcol))
      grp[,group_1 := 1:.N]
      d <- merge(d, grp, by = c('period',zcol,xcol), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
      
    }
   
    # reorder d by row_id
    d <- d[order(row_id_transf), ]
  } else {
    # set Aproj grouping to period if there is but one z value
    d$group_1 <- d$period
  }
  
  #Set up groups for interaction # 2
  if(length(zl)>1 & 'age' %in% interacting_gp_2_effects){ 
    message(sprintf('More than one unique %s found, initiating Z in GP 2', zcol))
    # look for sex dimension 
    if(length(xl)==1 | !('sex' %in% interacting_gp_2_effects)){ 
      # set A proj grouping. The ordering here must match the ordering of epsilon_int_1 in the template
      grp <- setDT(expand.grid(1:length(yl), 1:length(zl)))
      setnames(grp,c('Var1','Var2'),c('period',zcol))
      grp[,group_2 := 1:.N]
      d <- merge(d, grp, by = c('period',zcol), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
    } else if (length(xl)==2 & 'sex' %in% interacting_gp_2_effects){
      message(sprintf('More than one unique %s found, initiating Sex in GP 2', xcol))
      grp <- setDT(expand.grid(1:length(yl), 1:length(zl), 1:length(xl)))
      setnames(grp,c('Var1','Var2','Var3'),c('period',zcol, xcol))
      grp[,group_2 := 1:.N]
      d <- merge(d, grp, by = c('period',zcol,xcol), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
    }
    
    # reorder d by row_id
    d <- d[order(row_id_transf), ]
  } else {
    # set Aproj grouping to period if there is but one z value
    d$group_2 <- d$period
  }

  #Set up groups for SZ interaction
  if(use_sz_gp==TRUE){  
    d$group_sz <- d$agebin
  }
  
  #Set up groups for SZ interaction
  if(use_sx_gp==TRUE){  
    d$group_sx <- d$sex_id
  }
  
  
  #Set up groups for TX interaction
      # set A proj grouping. The ordering here must match the ordering of epsilon_tx in the template
      grp <- setDT(expand.grid(1:length(yl), 1:length(xl)))
      setnames(grp,c('Var1','Var2'),c('period',xcol))
      grp[,group_tx := 1:.N]
      d <- merge(d, grp, by = c('period',xcol), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
    # reorder d by row_id
    d <- d[order(row_id_transf), ]
  
    # set A proj grouping. The ordering here must match the ordering of epsilon_zx in the template
    grp <- setDT(expand.grid(1:length(zl), 1:length(xl)))
    setnames(grp,c('Var1','Var2'),c(zcol,xcol))
    grp[,group_zx := 1:.N]
    d <- merge(d, grp, by = c(zcol,xcol), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
    # reorder d by row_id
    d <- d[order(row_id_transf), ]
  
  # coordinates at data points. these are passed to TMB in long,lat so
  # we keep them and create another set for 3d coords
  coords   <- cbind(d$longitude,d$latitude)
  
  # if we have  mesh on s2, first convert coords to spherical to project to mesh
  data.locs <- coords ## long, lat
  if(mesh_int$manifold == "S2"){
    ## then the mesh is on the sphere and we need to use 3d coords
    data.locs <- lonlat3D(data.locs[, 1], data.locs[, 2])
  }
  
  # make a projection matrix from data to int_1 mesh            
  A.proj_int_1 <- inla.spde.make.A(mesh  = mesh_int,
                                  loc   = data.locs,
                                  group = d$group_1)
  
  # make a projection matrix from data to int_1 mesh            
  A.proj_sz <- inla.spde.make.A(mesh  = mesh_int,
                                   loc   = data.locs,
                                   group = d$group_sz)
  
  if(use_sz_gp==F) A.proj_sz=A.proj_int_1
  
  
  # make a projection matrix from data to int_1 mesh            
  A.proj_sx <- inla.spde.make.A(mesh  = mesh_int,
                                loc   = data.locs,
                                group = d$group_sx)
  
  if(use_sx_gp==F) A.proj_sx=A.proj_int_1
  
  #Make a projection vector from data to int_2
  A.proj_int_2 <- d[['group_2']] - 1
  if(use_tx_gp==TRUE) A.proj_tx <- d[['group_tx']] - 1 else A.proj_tx = A.proj_int_2
  if(use_zx_gp==TRUE) A.proj_zx <- d[['group_zx']] - 1 else A.proj_zx = A.proj_int_2
  
  # make a projection matrix from data to s mesh
  A.proj_s <- inla.spde.make.A(
    mesh = mesh_s,
    loc = data.locs
  )
  
  A.proj_t <- as.integer(d[[ycol]] - min(yl))
  A.proj_z <- as.integer(d[[zcol]] - min(zl))
  A.proj_x <- as.integer(d[[xcol]] - min(xl))
                                       
  # make a clean design matrix. make sure all fes appear in d
  fes  <- unlist(strsplit(fes, ' \\+ ')) # i dont think this is safe and we might want to consider changing to formulas
  if(!all(fes %in% names(d)))  
    stop('Check your fes argument, not all covariate names appear in d.')
  if(length(fes)!=0) {
    X_xp <- as.matrix(cbind(int=1, d[,c(fes),with=FALSE]))
  } else {
    X_xp <- as.matrix(cbind(int=rep(1,nrow(d))))
  }

  # add in age fixed effects
  if(num_z > 1 & use_z_fes == TRUE){ 
    message(sprintf('Adding fixed effects for levels of zcol (%s)',zcol))
    for(z in 2:num_z){
      X_xp <- cbind(X_xp, d[[zcol]] == z)
    } 
    colnames(X_xp) <- c(colnames(X_xp)[colnames(X_xp)!=''],paste0('FE_z_level__',2:num_z))
    exclude_cs <- c(exclude_cs,paste0('FE_z_level__',2:num_z))
  }
  
  # add in any additional fixed effects (not to be included in predict)
  if(exists('non_pred_fes') & !is.null('non_pred_fes')){
    message(sprintf('Adding in the following non predicted FEs: %s',paste(non_pred_fes,collapse=', ')))
    for(npf in non_pred_fes){
      if(!npf %in% names(d)){
        message(sprintf('Not including %s because it is not in the data table',npf))
        non_pred_fes <- non_pred_fes[non_pred_fes != npf]
      } else {
        if(min(d[[npf]]) == max(d[[npf]])){
          message(sprintf('Not including %s because it does not vary in the data',npf))
          non_pred_fes <- non_pred_fes[non_pred_fes != npf]          
        } else{
          X_xp <- cbind(X_xp,d[[npf]])
          colnames(X_xp) <- c(colnames(X_xp)[-ncol(X_xp)],npf)
          exclude_cs <- c(exclude_cs,npf)
        }
      }
    }
  } else {
    non_pred_fes <- c()
  }
  
  # cs_df. imports seegMBG
  cs_df <- getCentreScale(X_xp, exclude = c('int',exclude_cs))
  X_xp  <- centreScale(X_xp, df = cs_df)
  
  # get data range in case we want to clamp for prediction later
  clamper <- data.table(apply(X_xp,2,range))
  
  # create stacker_col_id to identify which columns in stacker have sum-to-1
  if(coefs.sum1 == TRUE) {
    stacker_col_id <- ifelse(colnames(X_xp) %in% stacker_names, 1, 0)
  } else {
    stacker_col_id <- rep(0, ncol(X_xp))
  }
  
  # sort nugget and RE indicators
  nugget     <- as.numeric(as.logical(nugget))
  country_re <- as.numeric(as.logical(country_re))
  nid_re     <- as.numeric(as.logical(nid_re)) # these do not get used in prediction
  
  # obtain d_i (unique observations)
  d_i <- d[first_entry==1, ]
  
  # check there is more than one observed country or nid if REs for those are set
  if(length(unique(d_i$country)) == 1 & country_re == TRUE){
    message('WARNING: Only found one unique country in this data frame, so turning off country random effects.')
    country_re <- FALSE
  }
  if(length(unique(d_i$nid)) == 1 & nid_re == TRUE){
    message('WARNING: Only found one unique NID in this data frame, so turning off NID random effects.')
    nid_re <- FALSE
  }
  
  # make a gaul/country/cntry_RE_idx mapping table
  md    <- get_location_code_mapping(shapefile_version = shapefile_version)
  mdsub <- md[ihme_lc_id %in% unique(as.character(d_i$country)),]
  if(nrow(mdsub) != length(unique(as.character(d$country)))){
    message(sprintf('get_location_code_mapping() COUNTRY NAMES: %s',paste(sort(mdsub$ihme_lc_id),collapse=', ')))
    message(sprintf('IN DATA COUNTRY NAMES: %s',paste(sort(unique(as.character(d$country))),collapse=', ')))
    stop('get_location_code_mapping() and countries in data not of matching lengths')
  }
  cntry_re_map <- data.table(
                    country   = mdsub$ihme_lc_id,
                    adm_code = mdsub$ADM_CODE,
                    re_id     = 0:(nrow(mdsub)-1))
  cntry_re_vec <- cntry_re_map$re_id[match(as.character(d_i$country),cntry_re_map$country)]
 
  # make an nid_re mapping table
  nid_re_map <- unique(select(d_i, nid)) %>% 
    mutate(re_id = 0:(n()-1))
  nidEFF <- select(d_i, nid) %>% left_join(nid_re_map, by="nid")
  nid_re_vec <- nidEFF$re_id
  countries_with1nid <- select(d_i, country, nid) %>% 
     group_by(country) %>% # group by country to see...
     mutate(nidCount=n()) %>% # the number of nid units per country
     filter(nidCount==1) %>% ungroup %>% select(country) %>% unique
  if(country_re==TRUE & nid_re==TRUE & nrow(countries_with1nid) > 0) {
    message(paste("WARNING! The following countries", countries_with1nid$country, 
                   "have only 1 NID. If you encounter convergence problems, you 
                  may want to remove NID random effects."))
  }

  # Pixel re
  pixelEFF <- unique(select(d, pixel_id)) %>% 
    mutate(pixel_id_transf = 1:n()) %>% 
    right_join(d, by = "pixel_id")
  
  pixel_re_vec <- pixelEFF$pixel_id_transf-1
  
  # set GP RE arrays, or matrix depending on if we have z and x dimensions
  if(num_z > 1 & 'age' %in% interacting_gp_1_effects) {
    if(num_x > 1 & 'sex' %in% interacting_gp_1_effects) {
      Epsilon_int_1 <- array(0, dim=c(mesh_int$n,length(yl),num_z, num_x))
    } else{
      Epsilon_int_1 <- array(0, dim=c(mesh_int$n,length(yl),num_z))
    }
  } else {
    Epsilon_int_1 <- array(0, dim=c(mesh_int$n,length(yl)))
  }
  
  if(num_x > 1 & 'sex' %in% interacting_gp_2_effects) {
    Epsilon_int_2 <- array(0, dim=c(length(yl),num_z, num_x))
    } else{
      Epsilon_int_2 <- array(0, dim=c(length(yl),num_z))
    }
  
  
  Epsilon_sz <- array(0, dim=c(mesh_int$n,length(zl)))
  Epsilon_sx <- array(0, dim=c(mesh_int$n,length(xl)))

  Epsilon_tx <- array(0, dim=c(length(yl),num_x))
  Epsilon_zx <- array(0, dim=c(length(zl),num_x))

  Epsilon_s <- rep(0, mesh_s$n)
  Epsilon_t <- rep(0, length(yl))
  Epsilon_z <- rep(0, num_z)
  Epsilon_x <- rep(0, num_x)

  # set up vectors of model family 
  # look for convention in the data of lik_fam_<<binom,gauss>>, if not there, default to config family
  lik_gaussian <- lik_binomial <- rep(0, nrow(d_i))
  if(('lik_fam_binom' %in% names(d_i)) & ('lik_fam_gauss' %in% names(d_i))) {
    lik_gaussian <- d_i$lik_fam_gauss
    lik_binomial <- d_i$lik_fam_binom
    message(sprintf('Found row specific data likelihood indicators, will use those. %i rows binom, %i rows gauss',
                    sum(lik_binomial),sum(lik_gaussian)))
  } else {
    if(indicator_family == 'binomial') {
      lik_binomial <- rep(1, nrow(d_i))
      message('Using indicator family binomial for all rows')
    } else if(indicator_family == 'gaussian') {
      lik_gaussian <- rep(1, nrow(d_i))
      message('Using indicator family gaussian for all rows')
    }
  }
  
  if(any(lik_gaussian+lik_binomial != 1))
    stop('Not all rows in your data have been assigned a model (binom or gauss), or some have been assigned multiple!')
  
  # also look for sd if already exists in the data for the gauss rows to use
  # This is useful for crosswalked values with some data uncertainty, convention is variable named sd_<<INDICATOR>>
  sd_i <- rep(0, nrow(d_i))
  if(paste0('sd_',indicator) %in%  names(d)){
    message('Found SD estimates to use for crosswalked values.')
    sd_i <- d_i[[paste0('sd_',indicator)]]
    sd_i[is.na(sd_i)] <- 0
  }
  
  
  # run a quick regression to get starting values for the fixed effects 
  # This can speed up model fitting if iterations are slow.
  # Note this fit is currently assuming binomial
  # Note sometimes the starting values can strongly impact results
  if(all(lik_binomial==1) & coefs.sum1 != TRUE){
    message('LM for starting fe vals')
    y <- (d_i[[indic]][lik_binomial==1]+.0001)/d_i$N[lik_binomial==1]
    y[y<=0] <- 0.001
    y[y>=1] <- 0.999
    fe_start <- round( unname( lm(qlogis(y) ~ -1 + X_xp[d$first_entry==1 & lik_binomial==1,])$coefficients ), 4)
    site_fe_start <- round( unname( lm(qlogis(y) ~ -1 + d_i[lik_binomial==1,]$ANC_01)$coefficients ), 4)
  } else {
    message('Default starting fe vals')
    fe_start <- rep(0,ncol(X_xp))
    site_fe_start <- 0
  }
  
   message(sprintf('starting values for fixed effects: %s',paste0(fe_start,collapse=', ')))
  
  
  # cannot allow a gaussian likelihood to also have a nugget in the linear term, it leads to issues
  if(nugget == 1 & all(lik_gaussian == 1)){
    message('WARNING:: Nugget in all gaussian model leads to identifiability issues. Removing nugget for you.')
    nugget <- 0
  }
  
  # check if user wants to scale gaussian variance by N, if not set them all to one in the gaussian rows
  n_i <- d_i$N
  if(scale_gaussian_variance_N == FALSE & sum(lik_gaussian)>1) {
    message('Not scaling gaussian error by N since scale_gaussian_variance_N == FALSE.')
    n_i[lik_gaussian==1] <- 1
  }

  # if there is only one country in the region, turn off country_re
  if(all(country_re == TRUE & length(get_adm0_codes(reg)) == 1)){
    message('WARNING: One country in this region, so turning off country random effects')
    country_re <- FALSE
  }
  

    #Set up ANC site IID effect
    site_index <-0:nrow(unique(d_i[first_entry==1 & ANC==TRUE, c('nid', 'cluster_id')])) #Number of unique ANC sites + 1 for survey data
    
    p <- unique(d_i[first_entry==1 & ANC==TRUE, c('nid', 'site')])
    d_i[ANC==FALSE,site_index:=1]   #Indexing surveys with a site index of 1, which is also applied to one ANC site. but this wont affect the effect for that site because the surveys always be zeroed out
    if(nrow(p)>0){
    for(i in 1:nrow(p)){
      n_id <- p$nid[i]
      site_id <- p$site[i]
      d_i[nid==n_id & site==site_id, site_index:=i-1]
    }
    }
    site_index_i <- d_i$site_index
    
    
    #Set up country-year IID effect
    #Set up ANC site IID effect
    error_index <-0:nrow(unique(d_i[first_entry==1, c('country', 'year')])) #Number of unique ANC sites + 1 for survey data
    
    p <- unique(d_i[first_entry==1, c('country', 'year')])
     if(nrow(p)>0){
      for(i in 1:nrow(p)){
        country_id <- p$country[i]
        year_id <- p$year[i]
        d_i[country==country_id & year==year_id, error_index:=i-1]
      }
    }

    
    if(use_observation_level_error_iid==TRUE){
      error_index <-0:nrow(d_i[first_entry==1]) #Number of unique ANC sites + 1 for survey data
      
      p <- unique(d_i[first_entry==1])
      if(nrow(p)>0){
        for(i in 1:nrow(p)){
          d_i[i, error_index:=i-1]
        }
      }
    }
    
    if(use_cyzx_error_iid==TRUE){
      error_index <-0:nrow(unique(d_i[first_entry==1, c('country', 'year', 'agebin', 'sex_id')])) #Number of unique ANC sites + 1 for survey data
      
      p <- unique(d_i[first_entry==1, c('country', 'year', 'agebin', 'sex_id')])
      if(nrow(p)>0){
        for(i in 1:nrow(p)){
          country_id <- p$country[i]
          year_id <- p$year[i]
          x_id <- p$sex_id[i]
          agebin_id <- p$agebin[i]
          d_i[country==country_id & year==year_id & sex_id==x_id & agebin==agebin_id, error_index:=i-1]
        }
      }
      
    }
    
    error_index_i <- d_i$error_index
    
    
    
    #Create Aproj for cre_z
    #Set up groups for cre_z
        # set A proj grouping. The ordering here must match the ordering of epsilon_int_1 in the template
        grp <- setDT(expand.grid(unique(cntry_re_map$country), 1:length(zl)))
        names(grp) <- c('country', 'agebin')
        grp[,group_cre_z := 1:.N]
        d <- merge(d, grp, by = c('country', 'agebin'), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
        d <- d[order(row_id_transf), ]
        
        A.proj_cre_z <- d$group_cre_z - 1
    
    Epsilon_cre_z <- array(0, dim=c(length(unique(cntry_re_map$country)), num_z))
    
    
    #Set up groups for cre_x
    # set A proj grouping. The ordering here must match the ordering of epsilon_int_1 in the template
    grp <- setDT(expand.grid(unique(cntry_re_map$country), 1:length(xl)))
    names(grp) <- c('country', 'sex_id')
    grp[,group_cre_x := 1:.N]
    d <- merge(d, grp, by = c('country', 'sex_id'), all.x = TRUE) # warning this may re-order things, so do not use d stuff from above
    d <- d[order(row_id_transf), ]
    
    A.proj_cre_x <- d$group_cre_x - 1
    
    Epsilon_cre_x <- array(0, dim=c(length(unique(cntry_re_map$country)), num_x))
    
    
    
    # print some messages for random effects
    if(nugget == TRUE)                message('USING PIXEL RANDOM EFFECTS')
    if(country_re == TRUE)            message('USING COUNTRY RANDOM EFFECTS')
    if(nid_re == TRUE)                message('USING NID RANDOM EFFECTS')
    if(int_gp_1_effect == TRUE)       message('USING int 1 GP')
    if(int_gp_2_effect == TRUE)       message('USING int 2 GP')
    if(main_space_effect  == TRUE)    message('USING S GP')
    if(main_time_effect == TRUE)      message('USING T GMRF')
    if(main_age_effect == TRUE)       message('USING Z GMRF')
    if(main_sex_effect == TRUE)       message('USING Sex GMRF')
    if(use_sz_gp == TRUE)             message('USING SZ GP')
    if(use_sx_gp == TRUE)             message('USING SX GP')
    if(use_tx_gp == TRUE)             message('USING TX GP')
    if(use_zx_gp == TRUE)             message('USING ZX GP')
    if(use_cre_z_gp == TRUE)             message('USING CRE-Z GP')
    if(use_cre_x_gp == TRUE)             message('USING CRE-X GP')
    if(use_covs == TRUE)              message('USING FEs')
    if(use_anc_corrections == TRUE)              message('USING IID SITE CORRECTION')
    
    
    # Build SPDE object (using INLA functions) and get prior in TMB readable format
    spde_int_list <- read_inla_prior_matern(spde_prior_int, mesh_int)
    spde_int <- spde_int_list$spde
    
    # Build SPDE object (using INLA functions) and get prior in TMB readable format
    spde_s_list <- read_inla_prior_matern(spde_prior_s, mesh_s)
    spde_s <- spde_s_list$spde    
    

  # Construct a list of all data necessary to TMB to fit
  Data <- list(
    num_z = num_z,                 # 3rd dimension for GP, 
    num_x = num_x,                 # 4th dimension for GP, 
    y_i = d_i[[indic]],      # Number of observed events in the cluster (N+ in binomial likelihood)
    n_i = d_i$N,  # Number of observed exposures in the cluster (N in binomial likelihood)
    c_re_i         = cntry_re_vec,          # vector of country ids, ( starting at zero because C)
    nid_re_i       = nid_re_vec,            # vector of survey ids, zero index added in cpp for null effect
    pixel_re_k = pixel_re_vec, # vector of pixel ids
    site_iid_i = site_index_i, # vector of site ids
    error_iid_i = error_index_i, # vector of site ids
    anc_01_i = d_i$ANC_01, # vector of pixel ids
    w_i = d_i$weight,        # Data weight for each row
    X_kj = X_xp,# Covariate design matrix
    M0_int  = spde_int$param.inla$M0,    # SPDE sparse matrix--interaction
    M1_int  = spde_int$param.inla$M1,    # SPDE sparse matrix--interaction
    M2_int  = spde_int$param.inla$M2,    # SPDE sparse matrix--interaction
    M0_s  = spde_s$param.inla$M0,    # SPDE sparse matrix--space-only
    M1_s  = spde_s$param.inla$M1,    # SPDE sparse matrix--space-only
    M2_s  = spde_s$param.inla$M2,    # SPDE sparse matrix--space-only
    Aproj_int_1        = A.proj_int_1,     # mesh to prediction point projection matrix for int_1
    Aproj_int_2        = A.proj_int_2,     # mesh to prediction point projection matrix for int_2
    Aproj_s          = A.proj_s,     # mesh to prediction point projection matrix for s
    Aproj_t          = A.proj_t,     # mesh to prediction point projection matrix for t
    Aproj_z          = A.proj_z,     # mesh to prediction point projection matrix for z
    Aproj_x          = A.proj_x,     # mesh to prediction point projection matrix for sex
    Aproj_sz          = A.proj_sz,     # mesh to prediction point projection matrix for space*age
    Aproj_sx          = A.proj_sx,     # mesh to prediction point projection matrix for space*sex
    Aproj_tx          = A.proj_tx,     # mesh to prediction point projection matrix for time*sex
    Aproj_zx          = A.proj_zx,     # mesh to prediction point projection matrix for age*sex
    Aproj_cre_z          = A.proj_cre_z,     # mesh to prediction point projection matrix for age*sex
    Aproj_cre_x          = A.proj_cre_x,     # mesh to prediction point projection matrix for age*sex
    ID_k = as.integer(d$row_id_trans-1),       # observation ID
    lik_gaussian_i = lik_gaussian,          # data likelihood for each row
    lik_binomial_i = lik_binomial,          # data likelihood for each row
    sd_i = sd_i,                  # crossalked standard deviation
    stacker_col_id = stacker_col_id, # vector indicating which col of X_xp is a stacker & sum-to-1 (1) or not (0)
    options = list(
      use_priors = 1,        # option1==1 use priors
      adreport_off = 1,        # option2==1 ADREPORT off
      pixel_random = nugget,   # option3==1 include nugget
      country_random = as.numeric(country_re),          # option4==1 country random effects
      NID_random = as.numeric(nid_re),   # option5==1 NID random effects
      main_space_effect = as.numeric(main_space_effect),
      main_time_effect  = as.numeric(main_time_effect),
      main_age_effect   = as.numeric(main_age_effect),
      main_sex_effect   = as.numeric(main_sex_effect),
      sz_gp_effect      = as.numeric(use_sz_gp),
      sx_gp_effect      = as.numeric(use_sx_gp),
      tx_gp_effect      = as.numeric(use_tx_gp),
      zx_gp_effect      = as.numeric(use_zx_gp),
      cre_z_gp_effect      = as.numeric(use_cre_z_gp),
      cre_x_gp_effect      = as.numeric(use_cre_x_gp),
      int_gp_1_effect      = as.numeric(int_gp_1_effect),
      int_gp_2_effect      = as.numeric(int_gp_2_effect),
      age_in_int_1_gp      = as.numeric('age' %in% interacting_gp_1_effects),
      sex_in_int_1_gp      = as.numeric('sex' %in% interacting_gp_1_effects),
      sex_in_int_2_gp      = as.numeric('sex' %in% interacting_gp_2_effects),
      use_covs           = as.numeric(use_covs),
      use_anc_corrections           = as.numeric(as.logical(use_anc_corrections)),
      use_error_iid_re       = as.numeric(as.logical(use_error_iid_re))),
    prior_log_pixelre_sigma = read_inla_prior_sigma(nugget_prior),
    prior_log_cre_sigma    = read_inla_prior_sigma(ctry_re_prior),
    prior_log_nidre_sigma  = read_inla_prior_sigma(nid_re_prior),
    prior_log_site_iid_sigma  = read_inla_prior_sigma(site_iid_prior),
    prior_log_error_iid_sigma  = read_inla_prior_sigma(error_iid_prior),
    prior_log_t_sigma  = read_inla_prior_sigma(t_sigma_prior),
    prior_log_z_sigma  = read_inla_prior_sigma(z_sigma_prior),
    prior_log_x_sigma  = read_inla_prior_sigma(x_sigma_prior),
    prior_log_int2_sigma  = read_inla_prior_sigma(int2_sigma_prior),
    prior_trho_mean              = read_inla_prior(trho_prior)$par1,
    prior_trho_sd                = read_inla_prior(trho_prior)$par2,
    prior_trho_me_mean              = read_inla_prior(trho_me_prior)$par1,
    prior_trho_me_sd                = read_inla_prior(trho_me_prior)$par2,
    prior_zrho_mean              = read_inla_prior(zrho_prior)$par1,
    prior_zrho_sd                = read_inla_prior(zrho_prior)$par2,
    prior_xrho_mean              = read_inla_prior(xrho_prior)$par1,
    prior_xrho_sd                = read_inla_prior(xrho_prior)$par2,
    prior_matern_int             = spde_int_list$prior,
    prior_matern_s               = spde_s_list$prior,
    prior_cre_zrho_mean              = read_inla_prior(cre_zrho_prior)$par1,
    prior_cre_zrho_sd                = read_inla_prior(cre_zrho_prior)$par2,
    prior_log_cre_z_sigma  = read_inla_prior_sigma(cre_z_sigma_prior),
    prior_cre_xrho_mean              = read_inla_prior(cre_xrho_prior)$par1,
    prior_cre_xrho_sd                = read_inla_prior(cre_xrho_prior)$par2,
    prior_log_cre_x_sigma  = read_inla_prior_sigma(cre_x_sigma_prior),
    fconstraints = tmb_cov_constraint(colnames(X_xp), cov_constraints),
    aggweight_k = d$agg_weight, #Weight for agg data ages
    frr = d$hiv_frr #fertility rate ratios for anc data
  )
  
  if(!exists('use_ar2')) use_ar2=F
  if(use_ar2==T) logit_rho=c(1,1) else logit_rho=2
  if(use_cre_specific_params==T) cre_gp_rho=rep(2, nrow(cntry_re_map)) else cre_gp_rho=2
  if(use_cre_specific_params==T) cre_gp_sigma=rep(-0.5, nrow(cntry_re_map)) else cre_gp_sigma=-0.5
  
  # Set staring values for parameters
  Parameters <- list(alpha_j          = fe_start,  # FE parameters alphas
                     site_fe          = 0.5, #site_fe_start,
                     logtau           = -4, #round(spde_int$param.inla$theta.mu[1],4), # Matern/AR tau
                     logkappa         = 4, #round(spde_int$param.inla$theta.mu[2],4), # Matern Range
                     logtau_sz           = -4, #round(spde_int$param.inla$theta.mu[1],4), # Matern/AR tau
                     logkappa_sz         = 4, #round(spde_int$param.inla$theta.mu[2],4), # Matern Range
                     logtau_sx           = -4, #round(spde_int$param.inla$theta.mu[1],4), # Matern/AR tau
                     logkappa_sx        = 4, #round(spde_int$param.inla$theta.mu[2],4), # Matern Range
                     trho             = logit_rho,                          # temporal rho
                     zrho             = 2,                          # 3rd dimension of GP rho 
                     xrho             = 2,                          # 3rd dimension of GP rho 
                     logtau_me           = -4, #round(spde_s$param.inla$theta.mu[1],4), # Matern/AR tau
                     logkappa_me         = 4, #round(spde_s$param.inla$theta.mu[2],4), # Matern Range
                     trho_me             = logit_rho,                          # temporal rho
                     zrho_me             = logit_rho,                          # 3rd dimension of GP rho
                     xrho_me             = 2,                          # 3rd dimension of GP rho 
                     trho_int2             = logit_rho,                          # temporal rho
                     zrho_int2             = logit_rho,                          # 3rd dimension of GP rho 
                     xrho_int2             = 1.5,                          # 3rd dimension of GP rho 
                     trho_tx             = 2,                     
                     xrho_tx             = 1.5,                     
                     zrho_zx             = logit_rho,                  
                     xrho_zx             = 1.5,                      
                     zrho_sz             = 1.5,                   
                     xrho_sx             = 1.5,                      
                     cre_zrho            = cre_gp_rho,
                     log_cre_z_sigma     = cre_gp_sigma,
                     cre_xrho            = cre_gp_rho,
                     log_cre_x_sigma     = cre_gp_sigma,
                     log_pixelre_sigma = -1,                            # log(SD) of the normal nugget term
                     log_site_iid_sigma = -1,                            # log(SD) of the normal site iid term
                     log_error_iid_sigma = -1,                            # log(SD) of the normal country year iid term
                     log_cre_sigma    = -1,                            # log(SD) of the normal country intercept term (later add slopes as vector)
                     log_nidre_sigma  = -1,                            # log(SD) of the normal NID intercept term  
                     log_gauss_sigma  = -1,                            # log(SD) of normal model
                     log_t_sigma  = -2,                            # log(SD) of the normal T-GMRF variance term 
                     log_z_sigma  = -0.5,                            # log(SD) of the normal Z-GMRF variance term 
                     log_x_sigma  = -0.2,                            # log(SD) of the normal X-GMRF variance term
                     log_int2_sigma  = -0.5,                            # log(SD) of the int2 variance term
                     log_tx_sigma  = -1,                            # log(SD) of the tx variance term
                     log_zx_sigma  = -1,                            # log(SD) of the zx variance term
                     Epsilon_int_1      = Epsilon_int_1,                   # Random Effects: GP locations
                     Epsilon_int_2      = Epsilon_int_2,                   # Random Effects: GP locations
                     Epsilon_s        = Epsilon_s,                   # Random Effects: GP locations
                     Epsilon_t        = Epsilon_t,                   # Random Effects: GP locations
                     Epsilon_z        = Epsilon_z,                   # Random Effects: GP locations
                     Epsilon_x        = Epsilon_x,                   # Random Effects: GP locations
                     Epsilon_sz        = Epsilon_sz,                   # Random Effects: GP locations
                     Epsilon_sx        = Epsilon_sx,                   # Random Effects: GP locations
                     Epsilon_tx        = Epsilon_tx,                   # Random Effects: GP locations
                     Epsilon_zx        = Epsilon_zx,                   # Random Effects: GP locations
                     Epsilon_cre_z     = Epsilon_cre_z,  
                     Epsilon_cre_x     = Epsilon_cre_x,  
                     pixel_re         = rep(0,max(pixel_re_vec)+1),    # Random Effects: Nugget Values
                     cntry_re         = rep(0,nrow(cntry_re_map)),     # Random Effects Values of country random effects (later add in slope stuff)
                     nid_re           = rep(0,max(nid_re_vec)+1),      # Random Effects Values of nid random effects (later add in slope stuff)
                     site_iid_re           = rep(0,max(site_index)+1),      # Random Effects Values of site iid effects (later add in slope stuff)
                     error_iid_re           = rep(0,max(error_index)+1)      # Random Effects Values of site iid effects (later add in slope stuff)
  )
  
  
  # put bounds on parameters (Note, this is especially important for rhos)
  L  <- c(rep(-10,ncol(X_xp)),-10,-10,-99,-99,-99,-10,-10,-10,-10,-10,-10,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-99,-10,-10,-10,-10,-10,-10,-10,-10,-10, -10, -99, -10, -99, -10) ## updated the rho limits from .99999 since I transformed them in the cpp 
  U  <- c(rep( 10,ncol(X_xp)), 10, 10, 99, 99, 99,10, 10,10, 10,10, 10,99, 99, 99,99, 99, 99, 99,  99, 99,99, 99, 99, 99, 99, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 99, 10, 99, 10)
  L[which(stacker_col_id==1)] <- -Inf
  U[which(stacker_col_id==1)] <- Inf
  pn <- c(rep('alpha_j',ncol(X_xp)),
          'logtau','logkappa',
          'trho','zrho', 'xrho',
          'logtau_me','logkappa_me',
          'logtau_sz','logkappa_sz',
          'logtau_sx','logkappa_sx',
          'trho_me','zrho_me', 'xrho_me',
          'trho_tx','xrho_tx', 'zrho_zx', 'xrho_zx','zrho_sz', 'xrho_sx',
          'trho_int2','zrho_int2', 'xrho_int2',
          'log_pixelre_sigma','log_cre_sigma','log_nidre_sigma',
          'log_t_sigma','log_z_sigma','log_x_sigma', 'log_int2_sigma','log_tx_sigma','log_zx_sigma',
          'log_site_iid_sigma', 'site_fe', 'log_cre_z_sigma', 'cre_zrho', 'log_cre_x_sigma', 'cre_xrho', 'log_error_iid_sigma')
  names(L) <- names(U) <- pn
  
  # return the list
  return(list(Data         = Data,
              Parameters   = Parameters,
              cs_df        = cs_df,
              clamper      = clamper,
              coords       = coords,
              mesh_int     = mesh_int,
              mesh_s       = mesh_s,
              cntry_re_map = cntry_re_map,
              non_pred_fes = non_pred_fes,
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
#'
#' @return list of tmb model objects: ADfun, opt, sdrep
#'
#' @useDynLib mbg_tmb_model
#'
#' @export
fit_mbg_tmb <- function(lbdcorerepo     = core_repo,
                        cpp_template    = 'mbg_tmb_model',
                        tmb_input_stack = input_data,
                        ADmap_list      = NULL,
                        control_list    = NULL,
                        optimizer       = 'nlminb',
                        sparse_ordering = TRUE,
                        int_gp_1_effs   = interacting_gp_1_effects,
                        int_gp_2_effs   = interacting_gp_2_effects  
){
  
  # compile the cpp file and dynload it
  #message('NOT compiling template!!! Did you make any changes to mbg_tmb_model.cpp before recompiling?')
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

  parlist <- NULL
  
  if(tmb_input_stack$Data$options$int_gp_1_effect==1){
    randompars <- c(randompars,"Epsilon_int_1")
  } else{   
      ADmap_list[['Epsilon_int_1']]   <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_int_1))
      parlist         <- c(parlist, 'logtau','logkappa')
      ADmap_list[['trho']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$trho))
  }
  
  if(tmb_input_stack$Data$options$int_gp_2_effect==1){
    randompars <- c(randompars,"Epsilon_int_2")
  } else{   
    ADmap_list[['Epsilon_int_2']]     <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_int_2))
    ADmap_list[['log_int2_sigma']]         <- factor(NA)
  }
  
  if(tmb_input_stack$Data$options$sz_gp_effect==1){
    randompars <- c(randompars,"Epsilon_sz")
  } else{   
    ADmap_list[['Epsilon_sz']]     <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_sz))
    parlist <- c(parlist, 'zrho_sz', 'logtau_sz','logkappa_sz')
  }
  
  if(tmb_input_stack$Data$options$sx_gp_effect==1){
    randompars <- c(randompars,"Epsilon_sx")
  } else{   
    ADmap_list[['Epsilon_sx']]     <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_sx))
    parlist <- c(parlist, 'xrho_sx', 'logtau_sx','logkappa_sx')
  }
  
  if(tmb_input_stack$Data$options$tx_gp_effect==1){
    randompars <- c(randompars,"Epsilon_tx")
  } else{   
    ADmap_list[['Epsilon_tx']]     <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_tx))
    ADmap_list[['log_tx_sigma']]         <- factor(NA)
    parlist <- c(parlist, 'trho_tx', 'xrho_tx')
  }
  
  if(tmb_input_stack$Data$options$zx_gp_effect==1){
    randompars <- c(randompars,"Epsilon_zx")
  } else{   
    ADmap_list[['Epsilon_zx']]     <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_zx))
    ADmap_list[['log_zx_sigma']]         <- factor(NA)
    parlist <- c(parlist, 'xrho_zx')
    ADmap_list[['zrho_zx']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$zrho_zx))
  }
  
  if(tmb_input_stack$Data$options$main_space_effect==1) {
    randompars <- c(randompars,"Epsilon_s")
  } else {
    ADmap_list[['Epsilon_s']]           <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_s))
    parlist         <- c(parlist, 'logtau_me','logkappa_me')
  }
  if(tmb_input_stack$Data$options$main_time_effect==1) {
    randompars <- c(randompars,"Epsilon_t")
  } else {
    ADmap_list[['Epsilon_t']]           <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_t))
    parlist                             <- c(parlist, 'log_t_sigma')
    ADmap_list[['trho_me']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$trho_me))
  }
  if(tmb_input_stack$Data$options$main_age_effect==1) {
    randompars <- c(randompars,"Epsilon_z")
  } else {
    ADmap_list[['Epsilon_z']]           <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_z))
    parlist                             <- c(parlist, 'log_z_sigma')
    ADmap_list[['zrho_me']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$zrho_me))
  }
  if(tmb_input_stack$Data$options$main_sex_effect==1) {
    randompars <- c(randompars,"Epsilon_x")
  } else {
    ADmap_list[['Epsilon_x']]           <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_x))
    parlist                             <- c(parlist, 'xrho_me', 'log_x_sigma')
  }
  
  if(tmb_input_stack$Data$options$cre_z_gp_effect==1){
    randompars <- c(randompars,"Epsilon_cre_z")
  } else{
    
    ADmap_list[['Epsilon_cre_z']]     <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_cre_z))
    if(use_cre_specific_params==T){
    ADmap_list[['log_cre_z_sigma']]         <- rep(factor(NA),length(tmb_input_stack$Parameters$log_cre_z_sigma))
    ADmap_list[['cre_zrho']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$cre_zrho))
    } else{
      ADmap_list[['log_cre_z_sigma']]         <- factor(NA)
      ADmap_list[['cre_zrho']]        <- factor(NA)
    }
  }
  
  if(tmb_input_stack$Data$options$cre_x_gp_effect==1){
    randompars <- c(randompars,"Epsilon_cre_x")
      } else{ 
    ADmap_list[['Epsilon_cre_x']]     <- rep(factor(NA),length(tmb_input_stack$Parameters$Epsilon_cre_x))
    if(use_cre_specific_params==T){
    ADmap_list[['log_cre_x_sigma']]         <- rep(factor(NA),length(tmb_input_stack$Parameters$log_cre_x_sigma))
    ADmap_list[['cre_xrho']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$cre_xrho)) 
    } else{
      ADmap_list[['log_cre_x_sigma']]         <- factor(NA)
      ADmap_list[['cre_xrho']]        <- factor(NA)
    }
  }
  
#  if(!('year' %in% int_gp_1_effs))  parlist <- c(parlist, 'trho')
  if(!('age' %in% int_gp_1_effs))   parlist <- c(parlist, 'zrho')
  if(!('sex' %in% int_gp_1_effs))   parlist <- c(parlist, 'xrho')
  
  if(!('year' %in% int_gp_2_effs))  ADmap_list[['trho_int2']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$trho_int2))
  if(!('age' %in% int_gp_2_effs))   ADmap_list[['zrho_int2']]        <- rep(factor(NA),length(tmb_input_stack$Parameters$zrho_int2))
  if(!('sex' %in% int_gp_2_effs))   parlist <- c(parlist, 'xrho_int2')
  
  if(tmb_input_stack$Data$options$use_anc_corrections==1){
    # if anc site option is on add site to randompars
    randompars <- c(randompars,'site_iid_re')
  } else {
    # if no site, add site related parameters to ignorelist
    ADmap_list[['site_iid_re']] <- rep(factor(NA),length(tmb_input_stack$Parameters$site_iid_re))
    parlist                     <- c(parlist, 'site_fe', 
                                     'log_site_iid_sigma')
  }
  
  if(tmb_input_stack$Data$options$use_error_iid_re==1){
    # if anc site option is on add site to randompars
    randompars <- c(randompars,'error_iid_re')
  } else {
    # if no site, add site related parameters to ignorelist
    ADmap_list[['error_iid_re']] <- rep(factor(NA),length(tmb_input_stack$Parameters$error_iid_re))
    parlist                     <- c(parlist,  'log_error_iid_sigma')
  }
  
if(tmb_input_stack$Data$options$pixel_random==1){
    # if pixel option is on add pixel to randompars
    randompars <- c(randompars,'pixel_re')
  } else {
    # if no pixel, add pixel related parameters to ignorelist
    ADmap_list[['pixel_re']] <- rep(factor(NA),length(tmb_input_stack$Parameters$pixel_re))
    parlist                     <- c(parlist, 'log_pixelre_sigma')
  }
  
  # if country RE option is on add cntry_re to randompars
  if(tmb_input_stack$Data$options$country_random==1) {
    randompars <- c(randompars,'cntry_re')
  } else { # if no nugget, add nug related parameters to ignorelist
    parlist                     <- c(parlist, 'log_cre_sigma')
    ADmap_list[['cntry_re']]         <- rep(factor(NA),length(tmb_input_stack$Parameters$cntry_re))
  }
  
  # if NID RE option is on add nid_re to randompars
  if(tmb_input_stack$Data$options$NID_random==1) {
    randompars <- c(randompars,'nid_re')
  } else { # if no nugget, add nug related parameters to ignorelist
    parlist                     <- c(parlist, 'log_nidre_sigma')
    ADmap_list[['nid_re']]         <- rep(factor(NA),length(tmb_input_stack$Parameters$nid_re))
  }
  
  # map out model sigma if no gaussian observations
  if(sum(tmb_input_stack$Data$lik_gaussian_i) == 0){
    ADmap_list[['log_gauss_sigma']]    <- factor(NA) 
  }
  
  
  for(par in parlist) ADmap_list[[par]] <- factor(NA)
  
  # print
  message(paste0('ADMAP_LIST: ',  paste0(names(ADmap_list),collapse=', ')))
  message(paste0('Random Pars: ', paste0(randompars,collapse=', ')))
  
  #If all of the random effects are used, set ADmap_list back to NULL
  if(length(ADmap_list) == 0) ADmap_list <- NULL
  
  # make the AD object
  message('Making AD object')
  obj <- MakeADFun(data       = tmb_input_stack$Data, 
                   parameters = tmb_input_stack$Parameters,  
                   map        = ADmap_list, 
                   random     = randompars, 
                   hessian    = TRUE, 
                   DLL        = cpp_template,
                   checkParameterOrder =F) #sprintf('%s/mbg_central/%s', lbdcorerepo, tmpl))
  
  # normalize
  obj <- normalize(obj, flag = "flag")
  
  # Reduce fill in of sparse Cholesky factorization (requires metis install of TMB)
  if(sparse_ordering){
    runSymbolicAnalysis(obj)
  }
  
  message(obj)
  message(obj$par)
  message(control_list)
  
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
  
  message('Successfully finished getting join precision. Now returning:')
  
  # return
  return(list(
    ADfun   = obj,
    opt     = opt0,
    sdrep   = SD0,
    fenames = colnames(tmb_input_stack$Data$X_kj)
  ))
  
}





#' @title Multivariate normal draws
#' @description Take multivariate normal draws given a mean vector and precision matrix
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
#' @param coefs.sum1 Logical. Were the coefficients constrained to sum-to-1 in the model fit
#'
#' @return a cell_preds object
#'
#' @export
predict_mbg_tmb <- function(samples,
                            seed             = NULL,
                            tmb_input_stack  = input_data,
                            model_fit_object = model_fit,
                            fes              = all_fixed_effects,
                            sr               = simple_raster,
                            yl               = year_list,
                            zl               = z_list,
                            xl               = sex_list,
                            transform        = 'inverse-logit',
                            covs_list        = cov_list,
                            clamp_covs       = FALSE,
                            cov_constraints = covariate_constraint_vectorize(config),
                            int_gp_1_effect           = as.logical(use_gp),
                            int_gp_2_effect           = as.logical(any(interacting_gp_2_effects != '')),
                            use_space_only_gp = as.logical(use_space_only_gp),
                            use_time_only_gmrf = as.logical(use_time_only_gmrf),
                            use_age_only_gmrf = as.logical(use_age_only_gmrf),
                            use_sex_only_gmrf = as.logical(use_sex_only_gmrf),
                            use_sz_gp = as.logical(use_sz_gp),
                            use_sx_gp = as.logical(use_sx_gp),
                            use_tx_gp = as.logical(use_tx_gp),
                            use_zx_gp = as.logical(use_zx_gp),
                            use_cre_z_gp = as.logical(use_cre_z_gp),
                            use_cre_x_gp = as.logical(use_cre_x_gp),
                            use_gbd_covariate = as.logical(use_gbd_covariate),
                            coefs.sum1 = coefs_sum1,
                            use_covs   = any(as.logical(use_stacking_covs), as.logical(use_raw_covs)),
                            int_gp_1_effs = interacting_gp_1_effects,
                            int_gp_2_effs = interacting_gp_2_effects,
                            mesh_s = mesh_s,
                            mesh_int = mesh_int) {

  sdrep     <- model_fit_object$sdrep
  
  # pull a few useful things from the input data stack
  cs_transform         <- tmb_input_stack$cs_df
  mesh_int                 <- tmb_input_stack$mesh_int
  mesh_s                 <- tmb_input_stack$mesh_s
  cntry_re_map         <- tmb_input_stack$cntry_re_map
  if(clamp_covs == TRUE) {
    clamper            <- tmb_input_stack$clamper
  } else {
    clamper            <- NULL
  }
  
  ## dummy raster
  cell_idx <- seegSDM:::notMissingIdx(sr)
  
  stacker_col_id <- tmb_input_stack$Data$stacker_col_id

  # set seed if it is requested
  if(!is.null(seed)) set.seed(seed)
  
  # vector of means
  mu    <- c(sdrep$par.fixed,sdrep$par.random)
  
  # simulate draws
  draws <- rmvnorm_prec(mu = mu , prec = sdrep$jointPrecision, n.sims = samples)
  message('Successfully sampled draws')
  ## separate out the draws
  parnames      <- c(names(sdrep$par.fixed), names(sdrep$par.random))
  
  if (int_gp_1_effect == T){
    epsilon_int_1_draws <- draws[parnames=='Epsilon_int_1',]
  }    
  if (int_gp_2_effect == T){
    epsilon_int_2_draws <- draws[parnames=='Epsilon_int_2',]
  }    
  if (use_space_only_gp == T){
    epsilon_s_draws <- draws[parnames=='Epsilon_s',]
  }
  if (use_time_only_gmrf == T){
    epsilon_t_draws <- draws[parnames=='Epsilon_t',]
  }
  if (use_age_only_gmrf == T){
    epsilon_z_draws <- draws[parnames=='Epsilon_z',]
  }
  if (use_sex_only_gmrf == T){
    epsilon_x_draws <- draws[parnames=='Epsilon_x',]
  }
  
  if (use_sz_gp == T){
    epsilon_sz_draws <- draws[parnames=='Epsilon_sz',]
  }
  
  if (use_sx_gp == T){
    epsilon_sx_draws <- draws[parnames=='Epsilon_sx',]
  }
  
  if (use_tx_gp == T){
    epsilon_tx_draws <- draws[parnames=='Epsilon_tx',]
  }
  
  if (use_zx_gp == T){
    epsilon_zx_draws <- draws[parnames=='Epsilon_zx',]
  }
  
  
  if (use_cre_z_gp == T){
    epsilon_cre_z_draws <- draws[parnames=='Epsilon_cre_z',]
  }
  
  if (use_cre_x_gp == T){
    epsilon_cre_x_draws <- draws[parnames=='Epsilon_cre_x',]
  }
  
  alpha_draws   <- draws[parnames=='alpha_j',]
  
  #Deal with some matrix issues if alpha_draws is 1
  if(length(alpha_draws)==samples) alpha_draws <- matrix(c(alpha_draws), nrow=1)
  
  # remove non-predicted FEs from draws so they are discareded from here out
  alpha_draws   <- alpha_draws[which(!model_fit_object$fenames %in% non_pred_fes),]
  model_fit_object$fenames <- model_fit_object$fenames[!model_fit_object$fenames %in% non_pred_fes]
  
  # separate out Z FE draws
  FE_z_draws    <- NULL
  
  # names of fes
  tmb_const <- tmb_cov_constraint(model_fit_object$fenames, cov_constraints)
  fes       <- unlist(strsplit(fes, ' \\+ '))
  
  # keep only covariates used in model, typically stacking
  covs_list <- covs_list[names(covs_list) %in% fes]
  
  # mask covariates to simple raster
  if(length(fes)>0 & length(covs_list)>0){
    for(l in 1:length(covs_list)) {
      covs_list[[l]]  <- crop(covs_list[[l]], extent(sr))
      covs_list[[l]]  <- setExtent(covs_list[[l]], sr)
      covs_list[[l]]  <- mask(covs_list[[l]], sr)
    }
  }
  
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
  
  ## add time periods and z periods as needed--full sample space
  message('Creating full sample space index')
  grp <- setDT(expand.grid(1:length(yl), 1:length(zl), 1:length(xl)))
  setnames(grp,c('Var1','Var2','Var3'),c('t','z','sx'))
  grp[,group := 1:.N]
  fullsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,t  := grp$t[grp$group==g]]
    tmp[,sx := grp$sx[grp$group==g]]
    tmp[,gp := g]
    fullsamplespace <- rbind(fullsamplespace,tmp)
  }
  fullsamplespace[,idx := 1:.N]
  
  ## add time periods and z periods as needed -- interaction #1 sample space
  message('Creating interacting gp 1 sample space index')
  if('year' %in% int_gp_1_effs) t <- length(yl) else if(int_gp_1_effect == TRUE) stop('interacting gp without year is not yet supported') else t <- 1
  if('age' %in% int_gp_1_effs)  z <- length(zl) else z <- 1 
  if('sex' %in% int_gp_1_effs)  x <- length(xl) else x <- 1 
  
  if(int_gp_1_effect == TRUE & !('space' %in% int_gp_1_effs)) stop('interacting gp without space is not yet supported')
  
  grp <- setDT(expand.grid(1:t, 1:z, 1:x))
  setnames(grp,c('Var1','Var2','Var3'),c('t','z','sx'))
  grp[,group := 1:.N]
  intsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,t  := grp$t[grp$group==g]]
    tmp[,sx := grp$sx[grp$group==g]]
    tmp[,gp := g]
    intsamplespace <- rbind(intsamplespace,tmp)
  }
  intsamplespace[,idx := 1:.N]
  
  #Sample space for SZ gp
  message('Creating sz gp sample space index')
  z=length(zl)
  grp <- setDT(expand.grid(1:z))
  setnames(grp,'z')
  grp[,group := 1:.N]
  szsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,z  := grp$z[grp$group==g]]
    tmp[,gp := g]
    szsamplespace <- rbind(szsamplespace,tmp)
  }
  szsamplespace[,idx := 1:.N]
  
  #Sample space for SX gp
  message('Creating sx gp sample space index')
  sex=length(xl)
  z=1
  t=1
  grp <- setDT(expand.grid(1:sex,z,t))
  setnames(grp,c('sex', 'z', 't'))
  grp[,group := 1:.N]
  sxsamplespace <- data.table()
  for(g in 1:max(grp$group)){
    tmp <- f_orig
    tmp[,sex  := grp$sex[grp$group==g]]
    tmp[,gp := g]
    sxsamplespace <- rbind(sxsamplespace,tmp)
  }
  sxsamplespace[,idx := 1:.N]
  
  # pull out covariates in format we expect them
  # a list of length periods with a brick of named covariates inside
  # now also including separation by age
  
  new_cl <- list()
  if(use_covs==FALSE){
    message('No covariates detected, predicting using intercept only.')
  } else {
    message(sprintf('%i covariates detected.',length(covs_list[which(!names(covs_list) == 'gbd_est')])))
    if(agesex_in_stacking==F){
      for(p in 1:length(yl)){
        new_cl[[p]] <- list()
        for(n in names(covs_list[which(!names(covs_list) == 'gbd_est')])){
          if(dim(covs_list[[n]])[3]==1) { # synoptic mean covariates
            new_cl[[p]][[n]] <- covs_list[[n]]
          } else if (dim(covs_list[[n]])[3]==length(yl)) { # time varying covariates
            new_cl[[p]][[n]] <- covs_list[which(!names(covs_list) == 'gbd_est')][[n]][[p]]
          } else { # error if there is a non-conforming year 
            stop(sprintf('Covariate %n is a brick with %i layers, while year_list has %i years',
                         n,dim(covs_list[[n]])[3],length(yl)))
          }
        }
        new_cl[[p]] <- brick(new_cl[[p]])
      }
    } else{
      for(x in xl){
        new_cl[[x]] <- list()
        for(z in zl){
          new_cl[[x]][[z]] <- list()
          for(y in yl-1999){
            new_cl[[x]][[z]][[y]] <- list()
            for(n in 1:length(unique(names(covs_list)))){
              new_cl[[x]][[z]][[y]][[n]] <- covs_list[[x*z + (n-1)]][[y]]
            }
            new_cl[[x]][[z]][[y]] <- brick(new_cl[[x]][[z]][[y]])
          }
        }
      }
      
    }
  }
  
  if(use_gbd_covariate==T){
    message('gbd covariate detected')
    new_cl_gbd<-list()
    for(x in xl){
      new_cl_gbd[[x]] <- list()
      for(z in zl){
        new_cl_gbd[[x]][[z]] <- list()
        for(y in yl-1999){
          new_cl_gbd[[x]][[z]][[y]] <- list()
          new_cl_gbd[[x]][[z]][[y]] <- covs_list[['gbd_est']][[paste('gbd_est',x,z,y+1999, sep='.')]]
        }
        new_cl_gbd[[x]][[z]][[y]] <- brick(new_cl_gbd[[x]][[z]][[y]])
        names(new_cl_gbd[[x]][[z]][[y]]) <-'gbd_est'
      }
    }
  }
  
  
  
  # get surface locs to project on to
  pcoords        <- cbind(x=fullsamplespace$x, y=fullsamplespace$y) ## used for cov raster extract
  intcoords        <- cbind(x=intsamplespace$x, y=intsamplespace$y) ## used for cov raster extract
  pcoords_s      <- cbind(x=fullsamplespace[gp == 1,]$x, y=fullsamplespace[gp == 1,]$y) ## used for cov raster extract
  szcoords        <- cbind(x=szsamplespace$x, y=szsamplespace$y) ## used for cov raster extract
  sxcoords        <- cbind(x=sxsamplespace$x, y=sxsamplespace$y) ## used for cov raster extract
  
  # extract cell values  from covariates, deal with timevarying covariates here
  cov_vals <- list()
  cov_vals_gbd <- list()
  for(x in 1:length(xl)) {
    cov_vals[[x]] <- list()
    cov_vals_gbd[[x]] <- list()
    for(z in 1:length(zl)){
      cov_vals[[x]][[z]] <- list()
      cov_vals_gbd[[x]][[z]] <- list()
      for(p in 1:length(yl)){
        if(length(fes)>0 & length(covs_list)>0) {
          # raster extract and keep only fes--extract only from rasters for specific age bin stackers!
          if(agesex_in_stacking==F) cov_vals[[x]][[z]][[p]] <- raster::extract(new_cl[[p]], pcoords[1:nrow(f_orig),])
          if(agesex_in_stacking==T) cov_vals[[x]][[z]][[p]] <- raster::extract(new_cl[[x]][[z]][[p]], pcoords[1:nrow(f_orig),])
          if(use_gbd_covariate==T) cov_vals_gbd[[x]][[z]][[p]] <- raster::extract(new_cl_gbd[[x]][[z]][[p]], pcoords[1:nrow(f_orig),])
          if(use_gbd_covariate==T) cov_vals[[x]][[z]][[p]]     <- cbind(cov_vals[[x]][[z]][[p]], cov_vals_gbd[[x]][[z]][[p]])
          
          #revert column names to the generic stacker names
          colnames(cov_vals[[x]][[z]][[p]])<-fes
          
          # If there is only a single covariate, convert from vector to matrix
          if( (length(covs_list)==1) & !('matrix' %in% class(cov_vals[[x]][[z]][[p]]))){
            cov_vals[[x]][[z]][[p]] <- matrix(cov_vals[[x]][[z]][[p]], ncol=1)
          }
          
          # transform raw covariate values (center scaled) if needed (i.e. if cs_tranform is not 1 0 for that variable)
          cov_vals[[x]][[z]][[p]] <- centreScale(cov_vals[[x]][[z]][[p]],cs_transform) 
          
          # clamp covariates if clamper is not null
          if(!is.null(clamper)){
            # message('Clamping')
            for(fe in fes){
              tmpvec <- cov_vals[[x]][[z]][[p]][,colnames(cov_vals[[x]][[z]][[p]])==fe]
              mn <- as.numeric(clamper[,fe,with=FALSE][1])
              mx <- as.numeric(clamper[,fe,with=FALSE][2])
              tmpvec[tmpvec<mn] <- mn
              tmpvec[tmpvec>mx] <- mx
              cov_vals[[x]][[z]][[p]][,colnames(cov_vals[[x]][[z]][[p]])==fe] <- tmpvec
            }
          }
          # add an intercept
          cov_vals[[x]][[z]][[p]] <- cbind(int = 1, cov_vals[[x]][[z]][[p]])
        } else {
          # if no covariates just do intercept only
          cov_vals[[x]][[z]][[p]] <- cbind(int = rep(1,nrow(f_orig)))
        }
        # if there is a z column, add on those fixed effects indicators
        if(length(zl) > 1 & use_z_fes == TRUE){
          tmpzmat <- matrix(0,ncol = (length(zl)-1), nrow = nrow(f_orig))
          colnames(tmpzmat) <- paste0('FE_z_level__',2:length(zl))
          for(zz in 2:length(zl))
            if(z == zz) 
              tmpzmat[,paste0('FE_z_level__',zz)] <- 1
          cov_vals[[x]][[z]][[p]] <- cbind(cov_vals[[x]][[z]][[p]], tmpzmat)
        }
      }
    }
  }
  # covariate values by alpha draws
  l_vals <- list()
  for(x in 1:length(xl)) {
    l_vals[[x]] <- list()
    for(z in 1:length(zl)){
      l_vals[[x]][[z]] <- list()
      for(p in 1:length(yl))  
      l_vals[[x]][[z]][[p]] <- cov_vals[[x]][[z]][[p]] %*% apply_sum_to_one(stacker_col_id,apply_constraints(tmb_const, rbind(alpha_draws,FE_z_draws)))
    }
    l_vals[[x]] <-  do.call("rbind",unlist(l_vals[[x]], recursive = FALSE))
  }
  cell_l <- do.call("rbind",l_vals)
  
  ## setup coords for GP projection. convert coords to spherical if
  ## using spherical modeling mesh. used if you made a mesh on s2
  if(mesh_int$manifold == "S2"){
    gp_coords <- lonlat3D(pcoords[, 1], pcoords[, 2])
    gp_coords_s <- lonlat3D(pcoords_s[, 1], pcoords_s[, 2])
    gp_coords_int <- lonlat3D(intcoords[, 1], intcoords[, 2])
    gp_coords_sz <- lonlat3D(szcoords[, 1], szcoords[, 2])
    gp_coords_sx <- lonlat3D(sxcoords[, 1], sxcoords[, 2])
  } else {
    gp_coords <- pcoords
    gp_coords_s <- pcoords_s
    gp_coords_int <- intcoords
    gp_coords_sz <- szcoords
    gp_coords_sx <- sxcoords
  }
  remove(pcoords)
  remove(pcoords_s)
  remove(intcoords)
  remove(szcoords)
  remove(sxcoords)
  gc()
  
  ## define grouping across periods
  groups_periods <- fullsamplespace$gp
  int_periods    <- intsamplespace$gp
  sz_periods    <- szsamplespace$z
  sx_periods    <- sxsamplespace$sex
  
  ### values of GP surfaces at each cell (long by nperiods)
  # if we have multiple zs then do this by z since its possible to throw a SuiteSparse 'Problem too large' error here. 
  message('Getting GP predictions')
  ##Interacting GP #1
  if(int_gp_1_effect){
    message('GP int 1')
    ## use inla helper functions to project the spatial effect.
    A.pred_int_1 <- inla.spde.make.A(
      mesh  = mesh_int,
      loc   = gp_coords_int,
      group = int_periods)

    if(length(zl) > 1 & 'age' %in% int_gp_1_effs){
      cell_int_1 <- list()
      if(length(xl) ==1 ) {
        for(zz in zl){
          cell_int_1[[zz]] <- as.matrix(A.pred_int_1[(which(intsamplespace$z==zz)),] %*% epsilon_int_1_draws)
        }
      } else if(length(xl) > 1 & !('sex' %in% int_gp_1_effs)) {
        for(xx in xl){
          cell_int_1[[xx]] <- list()
          for(zz in zl){
            cell_int_1[[xx]][[zz]] <- as.matrix(A.pred_int_1[(which(intsamplespace$z==zz)),] %*% epsilon_int_1_draws)
          }
          cell_int_1[[xx]] <- do.call('rbind',cell_int_1[[xx]])
        }
      } else if(length(xl) > 1  & 'sex' %in% int_gp_1_effs){
        for(xx in xl){
          cell_int_1[[xx]] <- list()
          for(zz in zl){
            cell_int_1[[xx]][[zz]] <- as.matrix(A.pred_int_1[(which(intsamplespace$z==zz & intsamplespace$sx==xx)),] %*% epsilon_int_1_draws)
          }
          cell_int_1[[xx]] <- do.call('rbind',cell_int_1[[xx]])
        }
      }
      cell_int_1 <- do.call('rbind',cell_int_1)
    } else{
      cell_int_1 <- as.matrix(A.pred_int_1 %*% epsilon_int_1_draws)
    }
    

  
    #restructure cell_int_1 based on which variables are intcluded in the interacting gp
    if(length(zl) > 1 & !('age' %in% int_gp_1_effs)){
      if(length(xl) > 1 & !('sex' %in% int_gp_1_effs)){
        cell_int_1 <- do.call('rbind', replicate(length(zl)*length(xl), cell_int_1, simplify=FALSE))
      } else if(length(xl) ==1){
        cell_int_1 <- do.call('rbind', replicate(length(zl), cell_int_1, simplify=FALSE))
      }
    }
    remove(A.pred_int_1)
    remove(epsilon_int_1_draws)
    gc()
  } 
  
  remove(intsamplespace)
  gc()
  
  ##SZ GP
  if(use_sz_gp){
    message('SZ GP')
    ## use inla helper functions to project the spatial effect.
    A.pred_sz <- inla.spde.make.A(
      mesh  = mesh_int,
      loc   = gp_coords_sz,
      group = sz_periods)
    
    cell_sz <- list()
    for(xx in xl){
      cell_sz[[xx]] <- list()
      for(zz in zl){
        cell_sz[[xx]][[zz]] <- list()
        for(yy in yl){
          cell_sz[[xx]][[zz]][[yy]] <- as.matrix(A.pred_sz[(which(szsamplespace$z==zz)),] %*% epsilon_sz_draws)
        }
        cell_sz[[xx]][[zz]] <- do.call('rbind',cell_sz[[xx]][[zz]])
      }
      cell_sz[[xx]] <- do.call('rbind',cell_sz[[xx]])
    }
    cell_sz <- do.call('rbind',cell_sz)
    
    remove(A.pred_sz)
    remove(epsilon_sz_draws)
    gc()  
  } 
  remove(szsamplespace)
  gc()

  
  ##SX GP
  if(use_sx_gp){
    message('SX GP')
    ## use inla helper functions to project the spatial effect.
    A.pred_sx <- inla.spde.make.A(
      mesh  = mesh_int,
      loc   = gp_coords_sx,
      group = sx_periods)
    
    cell_sx <- list()
    for(xx in xl){
      cell_sx[[xx]] <- list()
      for(zz in zl){
        cell_sx[[xx]][[zz]] <- list()
        for(yy in yl){
          cell_sx[[xx]][[zz]][[yy]] <- as.matrix(A.pred_sx[(which(sxsamplespace$sex==xx)),] %*% epsilon_sx_draws)
        }
        cell_sx[[xx]][[zz]] <- do.call('rbind',cell_sx[[xx]][[zz]])
      }
      cell_sx[[xx]] <- do.call('rbind',cell_sx[[xx]])
    }
    cell_sx <- do.call('rbind',cell_sx)
    
    remove(A.pred_sx)
    remove(epsilon_sx_draws)
    gc()
  } 
  remove(sxsamplespace)
  gc()
  
  ### values of GP S surface at each cell (long by nperiods)
  if(use_space_only_gp)  {
    
    message('Space-only GP')
    # make a projection matrix from data to s mesh
    A.pred_s <- inla.spde.make.A(
      mesh  = mesh_s,
      loc   = gp_coords_s)
    
      cell_s <- list()
        for(xx in 1:length(xl)){
          cell_s[[xx]] <- list()
          for(zz in 1:length(zl)){
            cell_s[[xx]][[zz]] <- list()
            for(yy in 1:length(yl)){
              cell_s[[xx]][[zz]][[yy]] <- as.matrix(A.pred_s %*% epsilon_s_draws)
            }
            cell_s[[xx]][[zz]] <- do.call('rbind',cell_s[[xx]][[zz]])
          }
          cell_s[[xx]] <- do.call('rbind',cell_s[[xx]])
        }
        cell_s <- do.call('rbind',cell_s)
        
        remove(A.pred_s)
        remove(epsilon_s_draws)
        gc()
  }
  
  #Get raster & raster*year vectors for use in GP-Z surface creation 
  s.vec <- 1:length(sr[cell_idx])  ## vector of spatial cell count
  st.vec <- rep(s.vec, length(yl)) ## expand s vec in time
  stz.vec <- rep(st.vec, length(zl)) ## expand s vec in time
  
  ##Interaction #2
  if(int_gp_2_effect) {
    message('GP int 2')
    cell_int_2 <- list()
    for(x in 1:length(xl)){
      cell_int_2[[x]] <- list()
      for(z in 1:length(zl)){
        cell_int_2[[x]][[z]] <- list()
        for(y in 1:length(yl)){
          cell_int_2[[x]][[z]][[y]] <- list()
          if(!('sex' %in% int_gp_2_effs)) {
            cell_int_2[[x]][[z]][[y]] <- sapply(
              epsilon_int_2_draws[y + length(yl)*(z-1), ],
              function(it) rep(it, length(s.vec))
            )
          } else{
            cell_int_2[[x]][[z]][[y]] <- sapply(
              epsilon_int_2_draws[y + length(yl)*(z-1) + length(yl)*length(zl)*(x-1), ],
              function(it) rep(it, length(s.vec))
            )
          }
        } 
        cell_int_2[[x]][[z]] <- do.call('rbind',cell_int_2[[x]][[z]])
      }
      
      cell_int_2[[x]] <- do.call('rbind',cell_int_2[[x]])
    }
    cell_int_2 <- do.call('rbind',cell_int_2)
  }
  
  ##TX interaction
  if(use_tx_gp) {
    message('GP TX')
    cell_tx <- list()
    for(x in 1:length(xl)){
      cell_tx[[x]] <- list()
      for(z in 1:length(zl)){
        cell_tx[[x]][[z]] <- list()
        for(y in 1:length(yl)){
          cell_tx[[x]][[z]][[y]] <- list()
          cell_tx[[x]][[z]][[y]] <- sapply(
            epsilon_tx_draws[y + length(yl)*(x-1), ],
            function(it) rep(it, length(s.vec))
          )
        } 
        cell_tx[[x]][[z]] <- do.call('rbind',cell_tx[[x]][[z]])
      }
      
      cell_tx[[x]] <- do.call('rbind',cell_tx[[x]])
    }
    cell_tx <- do.call('rbind',cell_tx)
  }
  
  ##ZX interaction
  if(use_zx_gp) {
    message('GP ZX')
    cell_zx <- list()
    for(x in 1:length(xl)){
      cell_zx[[x]] <- list()
      for(z in 1:length(zl)){
        cell_zx[[x]][[z]] <- list()
        for(y in 1:length(yl)){
          cell_zx[[x]][[z]][[y]] <- list()
          cell_zx[[x]][[z]][[y]] <- sapply(
            epsilon_zx_draws[z + length(zl)*(x-1), ],
            function(it) rep(it, length(s.vec))
          )
        } 
        cell_zx[[x]][[z]] <- do.call('rbind',cell_zx[[x]][[z]])
      }
      
      cell_zx[[x]] <- do.call('rbind',cell_zx[[x]])
    }
    cell_zx <- do.call('rbind',cell_zx)
  }
  
  ### values of GP T surface at each cell, across agebins
  if(use_time_only_gmrf) { 
    message('Time-only GP')
    cell_t<-list()
    for(x in 1:length(xl)){
      cell_t[[x]] <- list()
      for(z in 1:length(zl)) {
        cell_t[[x]][[z]]<-list()
        for(y in 1:length(yl)){
          cell_t[[x]][[z]][[y]] <- sapply(
            epsilon_t_draws[y, ],
            function(it) rep(it, length(s.vec))
          )
        }
        cell_t[[x]][[z]] <- do.call('rbind',cell_t[[x]][[z]])
      }
      cell_t[[x]] <- do.call('rbind',cell_t[[x]])
    }
    cell_t <- do.call('rbind',cell_t)
  }
  
  ### values of GP Z surface at each cell
  if(use_age_only_gmrf) {
    message('Age-only GP')
    cell_z<-list()
    for(x in 1:length(xl)){
      cell_z[[x]] <- list()
      for(z in 1:length(zl)) {
        cell_z[[x]][[z]]<-list()
        for(y in 1:length(yl)){
          cell_z[[x]][[z]][[y]] <- sapply(
            epsilon_z_draws[z, ],
            function(it) rep(it, length(s.vec))
          )
        }
        
        cell_z[[x]][[z]] <- do.call('rbind',cell_z[[x]][[z]])
      }
      
      cell_z[[x]] <- do.call('rbind',cell_z[[x]])
    }
    cell_z <- do.call('rbind',cell_z)
  }
  
  ### values of GP X surface at each cell
  if(use_sex_only_gmrf) {
    message('Sex-only GP')
    cell_x<-list()
    for(x in 1:length(xl)){
      cell_x[[x]] <- list()
      for(z in 1:length(zl)) {
        cell_x[[x]][[z]]<-list()
        for(y in 1:length(yl)){
          cell_x[[x]][[z]][[y]] <- sapply(
            epsilon_x_draws[x, ],
            function(it) rep(it, length(s.vec))
          )
        }
        cell_x[[x]][[z]] <- do.call('rbind',cell_x[[x]][[z]])
      }
      cell_x[[x]] <- do.call('rbind',cell_x[[x]])
    }
    cell_x <- do.call('rbind',cell_x)
  }
 
  if(use_inla_nugget==TRUE){ 
  cell_nug <- matrix(0L, nrow = dim(cell_l)[1], ncol = dim(cell_l)[2])
  if('log_pixelre_sigma' %in% parnames){
    if(exists('no_nugget_predict') & no_nugget_predict==TRUE){
      message('Pixel-level random effect not included in predict')
    } else {
      message('Adding pixel-level random effect')
      for(s in 1:samples)
        cell_nug[,s] <- rnorm(dim(cell_nug)[1],0,exp(draws[parnames=='log_pixelre_sigma',])[s])
    }
  }
  }
  
  # add the country random intercept if it was estimated in the model
  cell_cre_int <- matrix(0L, nrow = dim(cell_l)[1], ncol = dim(cell_l)[2])
  if('log_cre_sigma' %in% parnames){
    message('adding country random intercept')
    cre <- data.table(draws[parnames=='cntry_re',])
    cre[, re_id := 1:.N]
    cre_nodat <- c(rnorm(samples,0,exp(draws[parnames=='log_cre_sigma',])), 0) #zero on the end for a fake re id
    cre <- rbind(cre, unname(data.table(t(cre_nodat))))    # add a row of draws for countries with NO REs (random draws, could do zero?)
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
  
  if(use_cre_z_gp==T){
    message('adding country-level Z GP')
    cre_z <- data.table(draws[parnames=='Epsilon_cre_z',])
    cre_z[, group_cre_z := 1:.N]
    
    grp <- setDT(expand.grid(cntry_re_map$re_id+1, 1:length(zl)))
    names(grp) <- c('re_id', 'z')
    grp[,group_cre_z := 1:.N]
    
    cre_z <- merge(cre_z,grp,by='group_cre_z',all.x=TRUE)

    cre_z <- merge(fullsamplespace,cre_z,by=c('re_id', 'z'), all.x=TRUE)  #different 'by' value?
    cre_z <- cre_z[order(idx)]
    cre_z <- cre_z[,grep('V',colnames(cre_z)),with=FALSE]
    cell_cre_z <- as.matrix(cre_z)
    remove(cre_z)
    gc()
  }
  
  if(use_cre_x_gp==T){
    message('adding country-level x GP')
    cre_x <- data.table(draws[parnames=='Epsilon_cre_x',])
    cre_x[, group_cre_x := 1:.N]
    
    grp <- setDT(expand.grid(cntry_re_map$re_id+1, 1:length(xl)))
    names(grp) <- c('re_id', 'sx')
    grp[,group_cre_x := 1:.N]
    
    cre_x <- merge(cre_x,grp,by='group_cre_x',all.x=TRUE)
    
    cre_x <- merge(fullsamplespace,cre_x,by=c('re_id', 'sx'), all.x=TRUE)  #different 'by' value?
    cre_x <- cre_x[order(idx)]
    cre_x <- cre_x[,grep('V',colnames(cre_x)),with=FALSE]
    cell_cre_x <- as.matrix(cre_x)
    remove(cre_x)
    gc()
  }
  
  
  # add together linear and st components
  message('add all effects')
  pred_tmb <- cell_l  + cell_cre_int
  
  if(use_inla_nugget==TRUE) pred_tmb <- pred_tmb + cell_nug
  if(int_gp_1_effect)    pred_tmb <- pred_tmb + cell_int_1
  if(int_gp_2_effect)    pred_tmb <- pred_tmb + cell_int_2
  if(use_space_only_gp)  pred_tmb <- pred_tmb + cell_s
  if(use_time_only_gmrf) pred_tmb <- pred_tmb + cell_t
  if(use_age_only_gmrf)  pred_tmb <- pred_tmb + cell_z
  if(use_sex_only_gmrf)  pred_tmb <- pred_tmb + cell_x
  if(use_sz_gp)          pred_tmb <- pred_tmb + cell_sz
  if(use_sx_gp)          pred_tmb <- pred_tmb + cell_sx
  if(use_tx_gp)          pred_tmb <- pred_tmb + cell_tx
  if(use_zx_gp)          pred_tmb <- pred_tmb + cell_zx
  if(use_cre_z_gp)          pred_tmb <- pred_tmb + cell_cre_z
  if(use_cre_x_gp)          pred_tmb <- pred_tmb + cell_cre_x
  
  # transform
  if(transform=='inverse-logit') { 
    message('inverse-logit transforming')
    pred_tmb <- plogis(as.matrix(pred_tmb))
  } else {
    pred_tmb <- eval(parse(text=sprintf('%s(as.matrix(pred_tmb))',transform)))
  }
  
  # if there is more than one z and/or more than one x, then return a list of nested length xl by length zl cell_preds
  message('parse into list')
  if(length(zl) > 1){
    pred_tmb_list <- list()
    if(length(xl) == 1) {
      chunklength <- dim(pred_tmb)[1]/length(zl) 
      for(z in 1:length(zl)){
        pred_tmb_list[[z]] <- pred_tmb[((z-1)*chunklength+1):(chunklength*z),1:samples]
      }
    } else if (length(xl) > 1) {
      chunklength_x <- dim(pred_tmb)[1]/length(xl)
      chunklength_z <- chunklength_x/length(zl)
      for(x in 1:length(xl)){
        pred_tmb_list[[x]] <- list()
        add_x <- (x-1)*chunklength_x
        for(z in 1:length(zl)){
          pred_tmb_list[[x]][[z]] <- pred_tmb[((z-1)*chunklength_z+1+add_x):(chunklength_z*z + add_x),1:samples]
        }
      }
    }
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
    draws[which(fn == 'log_pixelre_sigma'),]   <- exp(draws[which(fn == 'log_pixelre_sigma'),])
    draws[which(fn == 'log_cre_sigma'),]      <- exp(draws[which(fn == 'log_cre_sigma'),])
    draws[which(fn == 'log_nidre_sigma'),]    <- exp(draws[which(fn == 'log_nidre_sigma'),])
    
    fn[fn == 'logtau']           <- 'tau'
    fn[fn == 'logkappa']         <- 'kappa'
    fn[fn == 'log_pixelre_sigma'] <- 'nugget_SD'
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
read_inla_prior_matern <- function(prior_string, mesh_int){
  prior_list <- eval(parse(text=prior_string[1]))
  
  spde_list <- build_spde_prior(prior_list, mesh_int, st_gp_int_zero=FALSE)
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

#' @title Apply constraints to fitted value draws
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

#' @title Apply sum-to-one to fitted value draws
#' 
#' @description Given draws of fixed effect values generated from a fitted TMB
#'   model and a vector of which columns are stackers with a sum-to-one in that 
#'   model, apply transformations on constrained fixed effects to reproduce how 
#'   they were incorporated in the model (constraining them sum-to-one). See
#'   https://en.wikipedia.org/wiki/Dirichlet_distribution#Gamma_distribution for
#'   details on how a transformation of (log) gamma distributed RVs can be 
#'   transformed into RVs with a dirichlet distribution. If a stacker column has 
#'   a label of 0 (unconstrained, the default) or for any non-stacker columns, 
#'   the untransformed draws will be returned.
#' 
#' @param stacker_col_id int vector, which columns have sum-to-one (1) or not (0)
#' @param FE_draws matrix, fitted draws of coefficients
#' 
#' @return matrix of transformed beta coefficients for fixed effect draws
#' 
apply_sum_to_one <- function(stacker_col_id, FE_draws){
  unnormalized_draws <- t(exp(FE_draws[stacker_col_id==1,]))
  normalized_draws <- unnormalized_draws/rowSums(unnormalized_draws)
  FE_draws[stacker_col_id==1,] <- t(normalized_draws)
  
  return(FE_draws)
}
