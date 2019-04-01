## Functions to work with covariates

## Load and crop covariates to the modeled area
load_and_crop_covariates <- function(fixed_effects, simple_polygon, agebin=1) {

  # get selected covs
  selected_covs <- strsplit(fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]

  # Hard-code central directories
  central_cov_dir <- paste0('<<<< FILEPATH REDACTED >>>>') # central folder for covariates
  central_tv_covs <- c('mss','msw','unrakedmss','unrakedmsw','sevwaste','sevstunt','matedu_yrs','unrakedmatedu_yrs','wocba','evi','lights_new','LST_day','total_pop','rates','malaria','fertility','fertility_infill','fertility_smooth','urban_rural', 'land_cover', 'LST_avg', 'gpcp_precip', 'aridity_cruts','malaria_pfpr')
  central_ntv_covs <- c('access','irrigation','LF','LF_vector','reservoirs','aridity','elevation','annual_precip','PET','dist_rivers_lakes','dist_rivers_only','lat','lon','latlon')

  # Load all temporally-varying covariates
  evi             <- brick(paste0(central_cov_dir, 'EVI_stack.tif'))
  lights_new      <- brick(paste0(central_cov_dir, 'NTL_stack.tif'))
  LST_day         <- brick(paste0(central_cov_dir, 'LST_day_stack.tif'))
  total_pop       <- brick(paste0(central_cov_dir, 'WorldPop_allStages_stack.tif'))
  urban_rural     <- brick(paste0(central_cov_dir, 'GHS_settlement_model_stack.tif'))
  land_cover      <- brick(paste0(central_cov_dir, 'landCover_stack.tif'))
  LST_avg         <- brick(paste0(central_cov_dir, 'LST_avg_stack.tif'))
  gpcp_precip     <- brick(paste0(central_cov_dir, 'GPCP_precip_stack.tif'))
  aridity_cruts   <- brick(paste0(central_cov_dir, 'cruts_ard_stack.tif'))

  # Load all temporally-nonvarying covariates

  # Human/Cultural Synoptic Rasters
  access          <- brick(paste0(central_cov_dir, 'synoptic_humCul_stack.tif'))$synoptic_humCul_stack.1
  irrigation      <- brick(paste0(central_cov_dir, 'synoptic_humCul_stack.tif'))$synoptic_humCul_stack.2
  LF              <- brick(paste0(central_cov_dir, 'synoptic_humCul_stack.tif'))$synoptic_humCul_stack.3
  LF_vector       <- brick(paste0(central_cov_dir, 'synoptic_humCul_stack.tif'))$synoptic_humCul_stack.4
  reservoirs      <- brick(paste0(central_cov_dir, 'synoptic_humCul_stack.tif'))$synoptic_humCul_stack.5

  # Environmental/Physical Synoptic Rasters
  aridity           <- brick(paste0(central_cov_dir, 'synoptic_envPhy_stack.tif'))$synoptic_envPhy_stack.1
  dist_rivers_lakes <- brick(paste0(central_cov_dir, 'synoptic_envPhy_stack.tif'))$synoptic_envPhy_stack.2
  dist_rivers_only  <- brick(paste0(central_cov_dir, 'synoptic_envPhy_stack.tif'))$synoptic_envPhy_stack.3
  elevation         <- brick(paste0(central_cov_dir, 'synoptic_envPhy_stack.tif'))$synoptic_envPhy_stack.4
  annual_precip     <- brick(paste0(central_cov_dir, 'synoptic_envPhy_stack.tif'))$synoptic_envPhy_stack.5
  PET               <- brick(paste0(central_cov_dir, 'synoptic_envPhy_stack.tif'))$synoptic_envPhy_stack.6

  lat     <- raster(paste0(central_cov_dir, 'lat.tif'))
  lon     <- raster(paste0(central_cov_dir, 'lon.tif'))
  latlon  <- lat*lon

  malaria         <- brick(paste0(central_cov_dir, 'malaria_infant_death_rate_stack.tif'))
  if('malaria' %in% selected_covs)  values(malaria)=log(as.matrix(malaria)+.01)
  fertility         <- brick(paste0(central_cov_dir, 'fertility_stack.tif'))
  fertility_smooth         <- brick(paste0(central_cov_dir, 'fertility_smooth_stack.tif'))
  fertility_infill         <- brick(paste0(central_cov_dir, 'fertility_infill_stack.tif'))

  mss         <- brick(paste0(central_cov_dir, 'mss_stack.tif'))
  msw         <- brick(paste0(central_cov_dir, 'msw_stack.tif'))
  unrakedmss  <- brick(paste0(central_cov_dir, 'unrakedmss_stack.tif')) #'unrakedmss_stack.tif'))
  unrakedmsw  <- brick(paste0(central_cov_dir, 'unrakedmsw_stack.tif')) #'unrakedmsw_stack.tif'))
  sevstunt    <- brick(paste0(central_cov_dir, 'ss_stack.tif')) #'sevstunt_stack.tif'))
  sevwaste    <- brick(paste0(central_cov_dir, 'sw_stack.tif')) #'sevwaste_stack.tif'))

  wocba       <- brick(paste0(central_cov_dir, 'WOCBA_stack.tif'))
  malaria_pfpr       <- brick(paste0(central_cov_dir, 'malaria_pfpr_stack.tif'))

  matedu_yrs         <- brick(paste0(central_cov_dir, 'matedu_yrs_stack.tif')) #'matedu_yrs_stack.tif'))
  unrakedmatedu_yrs  <- brick(paste0(central_cov_dir, 'unrakedmatedu_yrs_stack.tif')) #'unrakedmatedu_yrs_stack.tif'))


  u5m_dir <-paste0('<<<< FILEPATH REDACTED >>>>')
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_m0.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_m1_11.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_2q1.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_2q3.Rdata'))
  load(paste0(u5m_dir,'data/raw/covariates/national_mr_5q0.Rdata'))
  rates = get(paste0('rates_',agebin))
  names(rates)=paste0('rates.',1:4)
  if('rates' %in% selected_covs)  rates=extend(rates,extent(-180, 180, -90, 90),keepres=TRUE)

  # Add names to layers
  names(access)         <- "access"
  names(irrigation)     <- "irrigation"
  names(LF)             <- "LF"
  names(LF_vector)      <- "LF_vector"
  names(reservoirs)     <- "reservoirs"
  names(aridity)        <- "aridity"
  names(elevation)      <- "elevation"
  names(annual_precip)  <- "annual_precip"
  names(PET)            <- "PET"
  names(dist_rivers_only)  <- "dist_rivers_only"
  names(dist_rivers_lakes) <- "dist_rivers_lakes"
  names(lat)            <- "lat"
  names(lon)  <- "lon"
  names(latlon) <- "latlon"


  for(c in central_tv_covs){
    tmp=get(c)
    names(tmp)=rep(paste0(c,'.',1:4))
    assign(c,tmp)
  }

  # Construct list of covariates to GAM and use in model from fixed_effects parameter equation.
  num_covs <- length(selected_covs)
  lcovs <- list()
  for(i in 1:num_covs) {
    message(selected_covs[i])
    this_cov <- selected_covs[i]
    if(this_cov %in% central_ntv_covs) { # Add if it is from the temporally non-varying list.
      lcovs[[i]] <- get(this_cov)
    }
    if(this_cov %in% central_tv_covs) { # Add if it is from the temporally varying list.
      lcovs[[i]] <- get(this_cov)
    }
    names(lcovs)[i] <- this_cov
  }

  # Make sure covariate layers line up with raster we are modeling over
  for(l in 1:length(lcovs)) {
    message(names(lcovs)[l])
    print(lcovs[[l]])
    lcovs[[l]]  <- crop(lcovs[[l]], extent(simple_polygon))
    print(lcovs[[l]])
    lcovs[[l]]  <- setExtent(lcovs[[l]], simple_polygon) #, keepres=TRUE, snap=TRUE)
    print(lcovs[[l]])
    lcovs[[l]]  <- mask(lcovs[[l]], simple_polygon)
  }

  return(lcovs)

}

# save a csv of contributions of different variables in brt
    save_brt_contributions <- function(brt_mod_obj = trans_covs[[1]],
                                       a=age,r=reg,pa=pathaddin,returnx=FALSE){
      x=data.table(rbind(cbind(brt_mod_obj[[1]]$contributions,year=2000,age=a,reg=r),
                         cbind(brt_mod_obj[[2]]$contributions,year=2005,age=a,reg=r),
                         cbind(brt_mod_obj[[3]]$contributions,year=2010,age=a,reg=r),
                         cbind(brt_mod_obj[[4]]$contributions,year=2015,age=a,reg=r)))
       x$var=gsub(pattern='\\.[0-9]',replacement='',x=x$var) # remove .# in varnames
       fwrite(x,paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/',
                                      indicator, '/output/', run_date,'/',
                                      'brt_contribs_bin',pa,'.csv'))
       if(returnx) return(x)
    }



## Function to take input data and covariate layers and fit GAM models
#   Arguments:
#     df = data.table with variables "latitude", "longitude", "N", and specified indicator
#     lcovs = list of raster layers or bricks (for temporally varying). Currently assumes
#             four layers named like U5M periods.
#   Returns: Two-item list (bricks of time-varying and non-time-varying covariates)

gam_covs <- function(df, indicator_family, lcovs) {
  coords <- df[, c('longitude', 'latitude'), with=FALSE]
  coords$lat=as.numeric(coords$latitude)
  coords$long=as.numeric(coords$longitude)
  coords <- coords[, c('long', 'lat'), with=FALSE]

  if(indicator_family=="binomial") response <- cbind(died = df[, get(indicator)], lived = df[, N] - df[, get(indicator)])
  if(indicator_family=="gaussian") response <- cbind(outcome = df[, get(indicator)])

  extra_data <- data.frame(year = df$year)

  # fit gam
  # This should take a few minutes
  system.time(trans <- gamTrans(coords=coords,
                                response=response,
                                covs=lcovs,
                                family = indicator_family,
                                extra_terms = ~ year,
                                extra_data = extra_data,
                                bam = TRUE,
                                predict = TRUE,
                                condition = NULL,
                                condition_covs = NULL,
                                s_args = list(bs = 'ts', k = 3),
                                samfrac = 0.1,
                                use.chol = TRUE))

  ### IF trans$trans IS RETURNED AS A LIST, IT IS SPLIT INTO TEMPORAL AND NON-TEMPORAL COVARIATES
  if(class(trans$trans)=='list') {
    temporal=TRUE
  } else {
    temporal=FALSE
  }

  message("CLAMP AND SAVE")
  if(!temporal){
    # THEY ALL PASS
    # use chi-squared stats to determine covariate usefulness
    keep <- which(summary(trans$model)$chi.sq > 0.1)
    trans_ras <- trans$trans[[keep]]

    # clamp covariates
    # find most extreme vaaues of transformed covariates that were observed
    vals <- extract(trans_ras, coords[idx_fit, ])
    sry <- apply(vals, 2, range, na.rm = TRUE)

    # clamp the covariates to these values
    for (i in 1:nlayers(trans_ras)) {
      range <- sry[, colnames(sry) == names(trans_ras)[i]]
      trans_ras[[i]][trans_ras[[i]] < range[1]] <- range[1]
      trans_ras[[i]][trans_ras[[i]] > range[2]] <- range[2]
    }

    return(trans_ras)

  }

  # temporally varying covariates are present, save them all separately
  # non varying ones will be save in the same covs_transformed location as before
  if(temporal){

    # first clamp and save non temporally varying
    message('time invariant covariates')
    trans_ras=trans$trans$nT_vals_trans
    # clamp covariates
    # find most extreme vaaues of transformed covariates that were observed
    vals <- extract(trans_ras, coords)
    sry <- apply(vals, 2, range, na.rm = TRUE)


    # clamp the covariates to these values
    for (i in 1:nlayers(trans_ras)) {
      range <- sry[, colnames(sry) == names(trans_ras)[i]]
      trans_ras[[i]][trans_ras[[i]] < range[1]] <- range[1]
      trans_ras[[i]][trans_ras[[i]] > range[2]] <- range[2]
    }

    # If you only specify one non-varying term, it simply gets name "layer" in the GAM function. Rather than fixing it in there
    # this will check if that's the case and rename it here.
    all_rasters <- list()
    if(length(names(trans_ras))==1) {
      for(cov in c('access','irrigation','elevation','africa_lf')) {
        if(cov %in% names(lcovs)) {
          names(trans_ras) <- cov
          all_rasters[[cov]] <- trans_ras
        }
      }
    }

    # Now, this same process for the individual temporally varying covariates
    count <- length(all_rasters) + 1
    for(n in names(trans$trans$T_vals_trans)){
      message(n)
      message(count)

      trans_ras=trans$trans$T_vals_trans[[n]]

      # clamp covariates
      # find most extreme vaaues of transformed covariates that were observed
      vals <- extract(trans_ras, coords)
      sry <- apply(vals, 2, range, na.rm = TRUE)


      # clamp the covariates to these values
      for (i in 1:nlayers(trans_ras)) {
        range <- sry[, colnames(sry) == names(trans_ras)[i]]
        trans_ras[[i]][trans_ras[[i]] < range[1]] <- range[1]
        trans_ras[[i]][trans_ras[[i]] > range[2]] <- range[2]
      }

      all_rasters[[n]] <- trans_ras
      count <- count + 1

    }

    return(all_rasters)

  }

}



#################################################################################
### Takes Covariates and Data and returns Boosted Regression Trees Outputs
## Inputs:
    # df = data.table with variables "latitude", "longitude", "N", and specified indicator
    # lcovs = list of raster layers or bricks (for temporally varying).
    #         output of load_and_crop_covariates()
    # years: analysis years. defaults to 2000, 2005, 2010, 2015
    # weight: character of weight variable name in df, defaults to null
    # tc: tree.complexity for BRT, defaults to 4
    # lr: learning.rate for BRT, defaults to 0.005
    # bf: bagging.fraction for BRT, defaults to 0.75
    # return_only_raster: if TRUE, only returns raster of results, otherwise returns a list of BRT model object and prediction raster
## Outputs: see above
# Note: depends on gbm and dismo packages
#################################################################################

brt_covs     <- function(df,
                         indicator_family   = indicator_family,
                         lcovs              = cov_layers,
                         years              = c(2000,2005,2010,2015),
                         weight             = NULL,
                         tc                 = 4,
                         lr                 = 0.005,
                         bf                 = 0.75,
                         return_only_raster = TRUE) {

    require(gbm)
    require(dismo)

    # take lcovs, see which are bricks (time varying) and which arent
    ntv <- tv <- c()
    for(i in 1:length(lcovs))
      if(inherits(lcovs[[i]],"RasterLayer")) ntv=c(ntv,i) else tv=c(tv,i)


    # make sure the tv covariates have the same number of years as the years argument
    for(i in tv){
      y=dim(lcovs[[i]])[3]
      if(y==length(years)){
        message(sprintf('The time-varying covariate `%s` has %i years of data. Matching to argument `years`: %s. With layer 1 as earliest year. If this is an issue, please rearrange your input RasterBrick for time-varying covariates. \n\n',names(lcovs)[i], y,paste(years,collapse=', ')))
      } else stop(sprintf('The time-varying covariate `%s` has %i years of data. Which does not match the number of years in argument `years`: %s.',names(lcovs)[i], y,paste(years,collapse=', ')))
    }

    # run BRT by years
    message(sprintf('Running BRT on %i years of data independently. Result will be a RasterBrick with %i layers.',length(years),length(years)))

    x<-list()
    for(i in 1:length(years)){
        # subset only to years we need
        d <- df[df$year==years[i],]
        d <- na.omit(d)

        # keep only the variables of interest
        coords           <- d[, c('longitude', 'latitude'), with=FALSE]
        coords$latitude  <- as.numeric(coords$latitude)
        coords$longitude <- as.numeric(coords$longitude)

        # BRT function we use has no binomial, only pois with offset
        if(indicator_family %in% c('binomial','poisson'))  {
          indicator_family = 'poisson'
          offset =log(d$N)
          message('WARNING: For Poisson to work, need to round decimals in the response')
          d[,eval(indicator):=round(d[,eval(indicator),with=FALSE],0)]
        } else offset = NULL

        d <- d[,c(indicator,weight),with=FALSE]

        # extract the values of the covariates
        c <- brick(lcovs[ntv])
        for(j in tv)
          c <- addLayer(c,lcovs[[j]][[i]])
        d <- cbind(d,extract(c,coords))


        # learning brt
        set.seed(123)
      
        mod <- try(
                   gbm.step(data             = d,
                            gbm.y            = 1,
                            gbm.x            = names(c),
                            offset           = offset,
                            family           = indicator_family,
                            weights          = weight,
                            tree.complexity  = tc,
                            learning.rate    = lr,
                            bag.fraction     = bf),silent=TRUE)
        if(is.null(mod)){
          message('First BRT attempt failed. Lowering Learning Rate by 1/10')
          mod <- try(
                   gbm.step(data             = d,
                            gbm.y            = 1,
                            gbm.x            = names(c),
                            offset           = offset,
                            family           = indicator_family,
                            weights          = weight,
                            tree.complexity  = tc,
                            learning.rate    = lr*.1,
                            bag.fraction     = bf))
        }
        if(is.null(mod)){
          message('Second BRT attempt failed. Lowering Original Learning rate by 1/1000 AGAIN')
          mod <- try(
                   gbm.step(data             = d,
                            gbm.y            = 1,
                            gbm.x            = names(c),
                            offset           = offset,
                            family           = indicator_family,
                            weights          = weight,
                            tree.complexity  = tc,
                            learning.rate    = lr*.001,
                            bag.fraction     = bf))
        }
        if(is.null(mod)){
          message('Third BRT attempt failed. Slow learn plus low tree complexity')
          mod <- try(
                   gbm.step(data             = d,
                            gbm.y            = 1,
                            gbm.x            = names(c),
                            offset           = offset,
                            family           = indicator_family,
                            weights          = weight,
                            tree.complexity  = 2,
                            learning.rate    = lr*.001,
                            bag.fraction     = bf))
        }

        if(is.null(mod)) stop('ALL BRT ATTEMPTS FAILED')

        # save prediction rasters
        p <- predict(c,mod,n.trees=mod$gbm.call$best.trees,type='response')

        # save the outputs
        x[[paste0('m_',i)]]=mod
        x[[paste0('brt.',i)]]=p
        #x

    } # closes years parallel loop


  #extract objects and prediction surfaces as two separate lists and save them
  for(i in 1:length(x))
    assign(names(x)[i],x[[i]])

  m=(mget(paste0('m_',1:length(years))))
  p=(mget(paste0('brt.',1:length(years))))



  if(return_only_raster) {
    return(brick(p))
  } else {
    return(list(m,brick(p)))
  }


}

###A faster version of load_gbd_covariates()
#' @param gbd_fixed_effects a vector of strings where the values are covariate_short_names of GBD covariates
#' @param year_ids a numeric vector of year ids. E.g. 2000:2015. Usually fufilled by year_list
#' @param gaul_list doesn't actually do anything anymore. Just there for backward compatibility
#' @param template a raster layer of the buffered modelling area. usually it will be cov_layers[[1]][[1]]
#' @param use_subnationals logical. If true, the function will replace admin0 with a subnational units where possible.
#' 
load_gbd_covariates = function(gbd_fixed_effects, year_ids, gaul_list, template, use_subnationals = F){
  
  #load the analysis shapefile
  world_shape = readRDS('<<<< FILEPATH REDACTED >>>>/GBD2016_analysis_final.rds')
  
  #if we are not using subnationals, keep only national units; otherwise remove the parent units
  if(!use_subnationals){
    world_shape = world_shape[world_shape$level==3,]
  }else{
    world_shape = world_shape[!world_shape$loc_id %in% unique(world_shape$parent_id),]
  }
  
  world_shape = crop(world_shape, template)
  
  afras = rasterize(world_shape, template, 'GAUL_CODE', fun = 'last')
  
  #check to make sure all of the relevant gauls are in the file
  if(!all(world_shape$GAUL_CODE %in% unique(afras))){
    afras = raster::rasterize(world_shape[!world_shape$GAUL_CODE %in% unique(afras),], afras, 'GAUL_CODE', fun = 'first', update = T)
  }
  
  
  #loop over requested covariates
  fetch_gbd_cov = function(gbd_name, afras){
    gbd <- load_gbd_data(gbd_type = "covariate",
                         gbd_name = gbd_name,
                         gaul_list = unique(afras),
                         measure_id = 6,
                         age_group_id = '2 3 4 5',
                         metric_id = 3,
                         year_ids = year_ids,
                         return_by_age_sex = 'no',
                         collapse_age_sex = TRUE)
    
    #for each gaul code and year, update the values
    blank = brick(lapply(year_list, function(x) afras * NA))
    
    for(yyy in 1:length(year_ids)){
      for(ggg in unique(afras)){
        blank[[yyy]][which(raster::getValues(afras)==ggg)] = gbd[name == ggg & year == year_list[yyy], mean]
      }
    }
    
    names(blank) = rep(gbd_name, times = dim(blank)[3])
    
    return(blank)
  }
  
  all_gbd_layers = lapply(gbd_fixed_effects, function(x) fetch_gbd_cov(x, afras))
  names(all_gbd_layers) = gbd_fixed_effects
  
  return(all_gbd_layers)
  
}

load_mbg_covariates <- function(mbg_fixed_effects, simple_polygon) {

  # Hard-code MBG covariates that have best models available
  mbg_covs <- c('edu_0')

  # Load MBG covariates available
  edu_0    <- brick('<<<< FILEPATH REDACTED >>>>/education/edu_0/output/best/edu_0.tif')
  edu_mean <- brick('<<<< FILEPATH REDACTED >>>>/education/edu_mean/output/best/edu_mean.tif')

  # Construct list of covariates to GAM and use in model from fixed_effects parameter equation.
  selected_covs <- strsplit(mbg_fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
  num_covs <- length(selected_covs)
  lcovs <- list()
  for(i in 1:num_covs) {
    this_cov <- selected_covs[i]
    lcovs[[i]] <- get(this_cov)
    names(lcovs[[1]]) <- gsub('period_', paste0(this_cov,'.'), names(lcovs[[1]]))
    names(lcovs)[i] <- this_cov
  }

  # Make sure covariate layers line up with raster we are modeling over
  for(l in 1:length(lcovs)) {
    lcovs[[l]]  <- extend(lcovs[[l]],extent(-180, 180, -90, 90),keepres=TRUE)
    lcovs[[l]]  <- crop(lcovs[[l]], extent(simple_polygon))
    lcovs[[l]]  <- setExtent(lcovs[[l]], simple_polygon)
    lcovs[[l]]  <- mask(lcovs[[l]], simple_polygon)
  }

  return(lcovs)
}


load_and_crop_covariates_annual <- function(covs,              
                                            measures,             
                                            simple_polygon,
                                            start_year  = 1998,
                                            end_year    = 2017,
                                            interval_mo = 60,
                                            agebin=1) {

  # covariate directory
  covdir <- '<<<< FILEPATH REDACTED >>>>'

  # pull a vector of all available covariates
  all_covs <- list.dirs(path = covdir, full.names = TRUE, recursive = FALSE)
  all_covs <- gsub(paste0(covdir,'/'),'',all_covs)

  # duration
  if(interval_mo %in% c(12,24,60)){
    inyrs = TRUE
    dur   = paste0(interval_mo/12,'y')
    pers  = length(seq(start_year,end_year,interval_mo/12))
  } else {
    dur   = paste0(interval_mo,'m')
    stop('Monthly covariates not currently supported.')
  }
  message(paste0('Duration set to ',dur,'.'))
  message(paste0('Will search for covariates starting in ',start_year,' and ending in ',end_year,'.'))
  message(paste0('Interval years are:',paste(seq(start_year,end_year,interval_mo/12),collapse=', '),'.\n'))


  # CHECK: Covs and measures are of the same length
  if(class(covs)!='character')       stop('covs argument must be a character vector')
  if(class(measures)!='character')   stop('measuress argument must be a character vector')
  if(length(covs)!=length(measures)) stop('covs and measures vectors must be of same length')

  # CHECK: look through fixed effects and make sure we have covariates avaiilable
  mismatch <- covs[!covs %in% all_covs]
  if(length(mismatch)>0){
    stop(paste('You have selected some covariates in fixed_effects which do not exist:\n',paste0(mismatch,collapse=', ')))
  }

  # Pull covariates
  i=1
  covlist = list()
  message('Finding your covariates and  searching for preferred measures and durations...\n')
  for(c in covs){
    message(paste(c,'\n'))
    durtmp  = dur
    perstmp = pers
    # check for measure preference
    if(!dir.exists(paste0(covdir,c,'/',measures[i],'/')))
      stop(paste('measure',measures[i],'for',c,'does not exist.'))
    # check for duration preference
    if(!dir.exists(paste0(covdir,c,'/',measures[i],'/',dur,'/'))){
      if(dir.exists(paste0(covdir,c,'/',measures[i],'/','synoptic'))){
        message(paste(c,'is synoptic only.'))
        durtmp  = 'synoptic'
        perstmp = 1
      } else {
        stop(paste(dur,'duration for measure',measures[i],'for',c,'does not exist.'))
      }
    }

    # load data
    message('Loading in data...\n')
    covlist[[c]] = list()
    for(p in 1:perstmp){
      # deal with tv an ntv differently in naming
      if(durtmp!='synoptic') {
        yrtmp = paste0(durtmp,'_',start_year+(interval_mo/12)*(p-1),'_00_00')
      } else {
        yrtmp = 'synoptic'
      }
      f <- paste0(covdir,c,'/',measures[i],'/',durtmp,'/',c,'_',measures[i],'_',yrtmp,'.tif')
      if(!file.exists(f)) stop(paste('Searched for the following file and it does not exist:',f))
      covlist[[c]][[p]]=raster(f)

    }

    # collapse list to a rasterBrick
    covlist[[c]]        <- stack(covlist[[c]][1:perstmp])
    names(covlist[[c]]) <- rep(paste0(c,'.',1:perstmp))

    # convert to regular raster if synoptic
    if(durtmp == 'synoptic') {
      #covlist[[c]]  <- raster(  covlist[[c]]  )
      covlist[[c]] <- covlist[[c]][[1]]
      names(covlist[[c]]) <- c
    }

    i=i+1
  }

  # Make sure covariate layers line up with raster we are modeling over
  message('Cropping to simple_polygon')
  for(l in 1:length(covlist)) {
    message(names(covlist)[l])
    #print(covlist[[l]])
    covlist[[l]]  <- crop(covlist[[l]], extent(simple_polygon))
    #  print(covlist[[l]])
    covlist[[l]]  <- setExtent(covlist[[l]], simple_polygon) #, keepres=TRUE, snap=TRUE)
    #print(covlist[[l]])
    covlist[[l]]  <- mask(covlist[[l]], simple_polygon)
  }

  return(covlist)

}

## save_standard_covariate ################################################

#' Save an mbg output as a standard covariate
#'
#' @param cov_name Name for the covariate (how will it be called in the
#'                 fixed effects parameter?).  In most cases will be the same
#'                 as indicator.
#' @param ig Indicator groups
#' @param run_date Run date to pull from to make this covariate
#' @param ind Indicator (if different from cov name).  If `indicator` is NULL,
#'            (default) then `indicator` is set to `cov_name`. 
#' @param raked Do you want to use the raked version?
#' @param year_list Vector (numeric) of years, with length the same as the number
#'                  of layers in your modeled output .tif
#' @param measure Covariate measure (e.g. mean, median, etc.)
#' @return Creates a series of .tif files, one for each year, in the standard MBG
#'         covariate directory using appropriate naming convention.  Note that for
#'         consistency with other covariates, these are expanded to global rasters.
#'         Also saves a readme.txt file with details of the covariate creation.
#' @examples
#' save_standard_covariate(cov_name = "dpt3_cov_test",
#'                         ig = "vaccine",
#'                         run_date = "2017_10_27_12_55_26",
#'                         raked = T,
#'                         year_list = c(2000:2015),
#'                         measure = "mean")

save_standard_covariate <- function(cov_name,
                                    ig = indicator_group,
                                    run_date,                          
                                    ind = NULL, 
                                    raked = T, 
                                    year_list = c(2000:2016),
                                    measure = "mean") {

  # CONFIG ------------------------------------------------------------------
  if (is.null(ind)) ind <- cov_name

  # Set up directories
  std_cov_dir <- "<<<< FILEPATH REDACTED >>>>"
  sharedir <- paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", ind, "/output/", run_date, "/")

  # Load a global template
  template_file <- paste0(std_cov_dir, "worldpop/total/1y/worldpop_total_1y_2015_00_00.tif")
  template <- raster(template_file)

  # Read in raster files & format -------------------------------------------

  # Read in the covariate raster, trying different naming conventions
  if (raked == T) {
    raster_file <- paste0(sharedir, ind, "_", measure, "_raked_raster.tif")
  } else if (raked == F) {
    raster_file_1 <- paste0(sharedir, ind, "_", measure, "_raster.tif")
    raster_file_2 <- paste0(sharedir, ind, "_", measure, "_unraked_raster.tif")
    if (file.exists(raster_file_1)) {
      raster_file <- raster_file_1
    } else if (file.exists(raster_file_2)) {
      raster_file <- raster_file_2
    }
  }

  message("Loading raster file:")
  message(raster_file)
  r_brick <- brick(raster_file)

  if (length(year_list) != nlayers(r_brick)) {
    stop("Your year list does not match the number of layers in the raster brick.")
  }

  # Extend to global (using template)
  r_brick <- extend(r_brick, y = template, value = NA)

  # Write output files ------------------------------------------------------
  message("\n Writing covariate ", cov_name)

  ## Make dirs
  main_dir <- '<<<< FILEPATH REDACTED >>>>'
  out_dir <- paste0(main_dir, cov_name, '/', measure, '/1y/')
  dir.create(out_dir, recursive = T, showWarnings = F)

  ## Write raster layer for each year, raked and unraked
  for(i in 1:length(names(r_brick))) {
    subset_r_brick <- r_brick[[i]]
    year <- year_list[i]
    message(paste0("  Saving ", year))
    writeRaster(subset_r_brick,
                filename = paste0(out_dir, cov_name, "_", measure, '_1y_', year, "_00_00.tif"),
                format = "GTiff",
                overwrite = TRUE)
  }

  message('\nAll layers successfully saved in <<<< FILEPATH REDACTED >>>>')
  message(paste0('Covariate ', cov_name, ' with measure ', measure, ' ready to call in MBG config files in the fixed_effects parameter.'))
  
  # Save a readme file -----------------------------------------------------
  readme_file <- paste0(out_dir, "readme.txt")
  fileConn <- file(readme_file)
  writeLines(c("## MBG COVARIATE README #############################",
               "",
               paste0("Covariate generated by ", Sys.info()['user'], " on ", Sys.time(), "."),
               "",
               paste0("Covariate name: ", cov_name),
               "",
               "Details of MBG run used to produce this covariate:",
               paste0("  Indicator name: ", ind),
               paste0("  Indicator group: ", ig),
               paste0("  Run date: ", run_date),
               paste0("  Raked: ", raked), 
               "",
               "####################################################"), 
            fileConn)
  close(fileConn)
}