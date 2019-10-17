## Load and crop covariates to the modeled area
load_and_crop_covariates <- function(fixed_effects, simple_polygon, agebin=1) {

  # get selected covs
  selected_covs <- strsplit(fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]

  # Hard-code central directories
  central_cov_dir <- '<<<< FILEPATH REDACTED >>>>' # central folder per Lucas
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


  # some u5m additions
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

  # Add names to layersss
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
       fwrite(x,paste0('<<<< FILEPATH REDACTED >>>>'))
       if(returnx) return(x)
    }



## Function to take input data and covariate layers and fit GAM models
#   Arguments:
#     df = data.table with variables "latitude", "longitude", "N", and specified indicator
#     lcovs = list of raster layers or bricks (for temporally varying). Currently assumes
#             four layers named like U5M periods.
#   Returns: Two-item list (bricks of time-varying and non-time-varying covariates)
#   Notes: This function currently won't work if you have anything other than *four* distinct periods in your data and covariates.
#           This will be made more flexible.
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
    # find most extreme vaaues of transofmred covariates that were observed
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
    # find most extreme vaaues of transofmred covariates that were observed
    vals <- extract(trans_ras, coords)
    sry <- apply(vals, 2, range, na.rm = TRUE)


    # clamp the covariates to these values
    for (i in 1:nlayers(trans_ras)) {
      range <- sry[, colnames(sry) == names(trans_ras)[i]]
      trans_ras[[i]][trans_ras[[i]] < range[1]] <- range[1]
      trans_ras[[i]][trans_ras[[i]] > range[2]] <- range[2]
    }

    # If you only specify one non-varying term, it simply gets name "layer" in the GAM function. Rather than fixing it in there
    # I'm just going to check if that's the case and rename it here.
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
    #require(parallel)
    #require(foreach)
    #require(doMC)

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

    #registerDoMC(cores=length(years))
    #out <- foreach(i = 1:length(years),
    #        .packages=c('dismo', 'gbm', 'raster'),
    #        .export=c('df','years','indicator_family','weight','ntv','tv','lcovs')
    #        ) %dopar% {
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
        # TODO: throw a try-catch so if some years work it at least will return that, if it doesnt it will try different things (like changing the learning rate. )
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


###A faster version of load_gbd_covariates. I've left the old one for backwards capability (dccasey 8/23/2017)
#' @param covs A character vector listing GBD covariates/outputs to extract. For covariates, this
#'     should be indicated by covariate_name_short, while for outputs, this should be indicated by
#'     acause. Usually fulfilled by gbd_fixed_effects.
#' @param measures A character vector coresponding to 'covs' listing the type of GBD quantity for
#'     each item. Options are 'covariate' and 'outcome'. Usually fulfilled by gbd_fixed_effects_measures.
#' @param year_ids A numeric vector of year_ids. Usually fulfilled by year_list.
#' @param age_ids A string of age_ids. Usually fulfilled by gbd_fixed_effects_age.
#' @param template A raster layer of the buffered modelling area. usually it will be cov_layers[[1]][[1]].
#'     If NULL, a default template is loaded using load_and_crop_covariates_annual()
#' @param use_subnationals Logical. If true, the function will replace admin0 with a subnational units
#'     where possible. Use with caution because it's not been tested outside of Africa. It might not
#'     work for countries with multiple levels of subnational units (e.g. UK or India).
#' @param simple_polygon simple_polygon object used for the modeling region. made in load_simple_polygon
#' @param interval_mo number of months in a time unit. usually 12 to correspond 'annual' year_loadxs
load_gbd_covariates = function(covs, measures, year_ids, age_ids,
                               template, use_subnationals = F,
                               simple_polygon, interval_mo){

  # check to see if the template is class raster, if not it is most
  # likely NULL since we pass in cov_layers[[1]][[1]] by default and
  # that is only created if we loaded in geospatial covs. otherwise,
  # we load in a geospatial cov to use as template
  if(class(template) != 'RasterLayer'){
    message('Loading in raster template for GBD covs since template argument was not a RasterLayer')
    template <-  load_and_crop_covariates_annual(covs            = 'evi',
                                                 measures        = 'median',
                                                 simple_polygon  = simple_polygon,
                                                 start_year      = min(year_ids),
                                                 end_year        = max(year_ids),
                                                 interval_mo     = as.numeric(interval_mo))[[1]][[1]]
  }
  
  # Load the analysis shapefile
  ## since this is fixed and is not related to shapefile_version, we
  ## also fix the shapefile_version passed to load_gbd_data inside
  ## fetch_gbd_covs to be one that worked for this shapeset
  world_shape <- readRDS('<<<< FILEPATH REDACTED >>>>')
  shapefile_version <- '2018_08_28' ## to match the entries of world_shape

  # If we are not using subnationals, keep only national units; otherwise remove the parent units
  if (!use_subnationals) {
    world_shape <- world_shape[world_shape$level==3,]
  } else {
    world_shape <- world_shape[!world_shape$loc_id %in% unique(world_shape$parent_id),]
  }

  world_shape <- crop(world_shape, template)
  # we must skip using the link_table as it is not relevant to this shapefile
  afras <- rasterize_check_coverage(world_shape, template, 'GAUL_CODE', fun = 'last', link_table = NULL)
  afras <- crop_set_mask(afras, template)
  
  # Check to make sure all of the relevant gauls are in the file
  if (!all(world_shape$GAUL_CODE %in% unique(afras))) {
    afras <- rasterize_check_coverage(world_shape[!world_shape$GAUL_CODE %in% unique(afras),], afras, 'GAUL_CODE', fun = 'first', update = T)
    afras <- crop_set_mask(afras, template)
  }


  # Loop over requested covariates
  fetch_gbd_cov = function(name, measure, afras) {

    # Load country-level results
    message("  Loading: ", name)
    gbd <- load_gbd_data(gbd_type          = measure,
                         gbd_name          = name,
                         gaul_list         = unique(afras),
                         measure_id        = 5,
                         age_group_id      = age_ids,
                         metric_id         = 3,
                         year_ids          = year_ids,
                         return_by_age_sex = 'no',
                         shapefile_version = shapefile_version, ## should be most up-to-date modified GAUL
                                                                ## to match GBD2016_analysis_final.rds
                                                                ## world shapefile
                         collapse_age_sex  = TRUE)


    if (nrow(gbd) != nrow(unique(gbd[, list(name, year)]))) stop(paste0(name, "is not unique by location-year"))

    # For each gaul code and year, update the values
    blank = brick(lapply(year_list, function(x) afras * NA))

    for (yyy in 1:length(year_ids)) {
      for (ggg in unique(afras)) {
        blank[[yyy]][which(raster::getValues(afras) == ggg)] = gbd[name == ggg & year == year_list[yyy], mean]
      }
    }

    names(blank) <- rep(name, times = dim(blank)[3])

    return(blank)
  }

  all_gbd_layers <- lapply(1:length(covs), function(x) fetch_gbd_cov(covs[x], measures[x], afras))
  names(all_gbd_layers) <- covs

  return(all_gbd_layers)
}

load_mbg_covariates <- function(mbg_fixed_effects, simple_polygon) {

  # Hard-code MBG covariates that have best models available
  mbg_covs <- c('edu_0')

  # Load MBG covariates available
  edu_0    <- brick('<<<< FILEPATH REDACTED >>>>')
  edu_mean <- brick('<<<< FILEPATH REDACTED >>>>')

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


load_and_crop_covariates_annual <- function(covs,     # covs     = c('evi','lstday','irrigation')
                                            measures, # measures = c('median','median','mean')
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
    stop('THIS FUNCTION IS NOT CURRENTLY SET UP TO DO MONTHLY, CONTACT LBD CENTRAL CODE.')
  }
  message(paste0('Duration set to ',dur,'.'))
  message(paste0('Will search for covariates starting in ',start_year,' and ending in ',end_year,'.'))
  message(paste0('Interval years are:',paste(seq(start_year,end_year,interval_mo/12),collapse=', '),'.\n'))

  ## make vector of all periods
  all.pers <- seq(start_year,end_year,interval_mo/12)

  ## CHECK: Covs and measures are of the same length
  if(class(covs)!='character')       stop('covs argument must be a character vector')
  if(class(measures)!='character')   stop('measuress argument must be a character vector')
  if(length(covs)!=length(measures)) stop('covs and measures vectors must be of same length')

  ## CHECK: look through fixed effects and make sure we have covariates avaiilable
  mismatch <- covs[!covs %in% all_covs]
  if(length(mismatch)>0){
    stop(paste('You have selected some covariates in fixed_effects which do not exist:\n',paste0(mismatch,collapse=', ')))
  }

  ## Pull covariates
  i <- 1
  covlist = list()
  message('Finding your covariates and  searching for preferred measures and durations...\n')
  for(c in covs){
    message(paste(c,'\n'))
    durtmp  <- dur
    perstmp <- pers
    
    ## check for measure preference
    if(!dir.exists(paste0(covdir,c,'/',measures[i],'/')))
      stop(paste('measure',measures[i],'for',c,'does not exist.'))

    ## check for duration preference
    if(!dir.exists(paste0(covdir,c,'/',measures[i],'/',dur,'/'))){
      if(dir.exists(paste0(covdir,c,'/',measures[i],'/','synoptic'))){
        message(paste(c,'is synoptic only.'))
        durtmp <- 'synoptic'
        perstmp <- 1
      } else {
        stop(paste(dur,'duration for measure',measures[i],'for',c,'does not exist.'))
      }
    }

    ## check if all years are available for the given
    ## cov-measure-duration combo.
    ## if some times are missing, copy over from the closest
    ## available time giving preference to earlier over later years
    ## when presented with ties
    if( durtmp != 'synoptic' ){ ## if we have different rasters for different times
      ## get all available periods by checking available files in the raster library
      all.file.pers <- list.files(paste0(covdir,c,'/',measures[i],'/',durtmp,'/'))
      all.file.pers <- all.file.pers[grep(durtmp, all.file.pers)]
      ## these filenames are of the form:
      ## <covariate>_<measure>_<duration>_<period>.tif. We want to
      ## grab period by splitting off the file ending, and then
      ## splitting on underscores
      all.file.pers <- unique(unlist(lapply(strsplit(all.file.pers, split = '.', fixed = T), function(x){x[1]})))
      all.file.pers <- as.numeric(sort(unlist(lapply(strsplit(all.file.pers, split = '_', fixed = T), function(x){x[length(x) - 2]}))))

      ## check to see if we have all the ones we need/want
      if(length(setdiff(all.pers, all.file.pers)) > 0){ ## i.e. we're missing periods
        missing.pers <- setdiff(all.pers, all.file.pers)
        message("WARNING! You are trying to load a raster covariate but the following years are missing:")
        message(paste(missing.pers, collapse = ' ', sep = ''))
        message("WARNING! We will map adjacent nearby years to these missing periods to fill in your dataset")
      }else{
        missing.pers <- NULL
      }
    }
    
    ## load data
    message('Loading in data...\n')
    covlist[[c]] = list()
    for(p in 1:perstmp){
      # deal with tv an ntv differently in naming
      if(durtmp!='synoptic') {
        this.yr <- all.pers[p]

        ## setup filepath to correct year, unless it's mising, then we overwrite with a neighbor period
        yrtmp = paste0(durtmp,'_',start_year+(interval_mo/12)*(p-1),'_00_00')
        ## grab a non-missing period if the period is missing
        ## grab nearest existing neighbor - preference to older periods in case of tie
        if(!is.null(missing.pers)){
          if(this.yr %in% missing.pers){
            ## get distance to existing years
            per.dist <- abs(this.yr - all.file.pers)
            yrtmp = paste0(durtmp,'_',all.file.pers[which.min(per.dist)],'_00_00')
            message(sprintf("WARNING! We are substituting in %s data from period: %i to use as if it were for period: %i",
                            c, all.file.pers[which.min(per.dist)], this.yr))
          }
        }
        
      } else {
        yrtmp = 'synoptic'
      }
      f <- paste0(covdir,c,'/',measures[i],'/',durtmp,'/',c,'_',measures[i],'_',yrtmp,'.tif')
      if(!file.exists(f)) stop(paste('Searched for the following file and it does not exist:',f))
      covlist[[c]][[p]]=raster(f)

    }

    ## collapse list to a rasterBrick
    covlist[[c]]        <- stack(covlist[[c]][1:perstmp])
    names(covlist[[c]]) <- rep(paste0(c,'.',1:perstmp))

    ## convert to regular raster if synoptic
    if(durtmp == 'synoptic') {
      #covlist[[c]]  <- raster(  covlist[[c]]  )
      covlist[[c]] <- covlist[[c]][[1]]
      names(covlist[[c]]) <- c
    }

    i <- i+1 ## update measure index for next cov
  }

  # Make sure covariate layers line up with raster we are modeling over
  message('Cropping to simple_polygon')
  for(l in 1:length(covlist)) {
    message(names(covlist)[l])
    covlist[[l]]  <- crop(covlist[[l]], extent(simple_polygon))
    ## setting the extent to simple_polygon here can cause issues when trying to align the covs to simple_raster later...
    ## I'm 99.99% sure this is safe to remove - but it's been in our code for ages. I left it here in case we were wrong and it's actually needed - azimmer
    ## covlist[[l]]  <- setExtent(covlist[[l]], simple_polygon) #, keepres=TRUE, snap=TRUE)
    covlist[[l]]  <- raster::mask(covlist[[l]], simple_polygon)
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
                                    year_list = c(2000:2015),
                                    measure = "mean") {

  # CONFIG ------------------------------------------------------------------
  # runtime configuration for central cluster
  if (Sys.info()["sysname"] == "Linux") {
    j_root <- '<<<< FILEPATH REDACTED >>>>'
  } else {
    stop("Run this on the cluster, please.")
  }

  if (is.null(ind)) ind <- cov_name

  # Set up directories
  std_cov_dir <- '<<<< FILEPATH REDACTED >>>>'
  sharedir <- '<<<< FILEPATH REDACTED >>>>'

  # Load a global template
  template_file <- '<<<< FILEPATH REDACTED >>>>'
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
  out_dir <- '<<<< FILEPATH REDACTED >>>>'
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



#' Checks for covariate issues
#'
#' This function looks for issues with covariates that will cause issues down the line and tries to fix them. The two issues looked for here are uniform covariates and pixel coverage. Uniform Covariates: If one of the extracted variables does not vary over the data area. Pixel Coverage: If one of the covariates is missing in large parts of your model region, these areas will be NA in your results later on. This often happens if you are using a modeled covariate in a new or partially new region which was not previously modelled. This function takes in all parallel model objects that would need to change if a covariate was removed for the code below in the parallel mode to work. 
#'
#' @param cc object to change: cs_covs
#' @param afe object to change: all_fixed_effects
#' @param afeb object to change: all_fixed_effects_brt
#' @param acl object to change: all_cov_layers
#' @param tc object to change: the_covs
#' @param check_uniform boolean, do a check for uniform covariates, defaults to TRUE
#' @param check_pixelcount boolean, do a check for covariates with too few pixel (ie too low of geographic coverage in the region), defaults to TRUE
#' @param check_pixelcount_thresh for a pixelcount check, what proportion of the maximum observed pixel coverage is needed to keep the covariate in? must be between 0 and 1. defaults to 0.95
#'
check_for_cov_issues   <- function(cc                      = cs_covs,
                                   afe                     = all_fixed_effects,
                                   afeb                    = all_fixed_effects_brt,
                                   fe                      = fixed_effects,
                                   tc                      = the_covs,
                                   acl                     = all_cov_layers,
                                   check_uniform           = TRUE,
                                   check_pixelcount        = TRUE,
                                   check_pixelcount_thresh = 0.95){
  
  # check that threshold is between 0 and 1
  if(check_pixelcount_thresh < 0 | check_pixelcount_thresh > 1){
    stop('check_pixelcount_thresh must be between 0 and 1')
  }
  
  # make a few useful objects to track names and dropped covs
  covs <- as.character(cc$cs_df$name)
  covs <- covs[!grepl('gaul_code', covs)]
  dropcovs <- c()
  
  # run through extracted data and find any non-varying covariates
  if(check_uniform == TRUE){
    message('Checking for uniform covariates ... ')
    for(covname in covs){
      if (all( abs(na.omit(cc$covs[[covname]]) - mean(na.omit(cc$covs[[covname]]))) == 0 )){
        message(sprintf('WARNING: %s did not vary in your data and is being removed as a covariate', covname))
        dropcovs <- c(dropcovs,covname)
      }
    }
    if(length(dropcovs) == 0){
      message('  All clear')
    }
  }
  
  # check to see if any of the covariates have significant area missing, if so drop them 
  # this typically happens if you use a modelled surface as a covariate in a region that was missing one or more countries
  # count the number of pixels in each covariate, make sure they are within 95% non-NA of the max covered covariate
  if(check_pixelcount == TRUE){
    message('Checking for pixel counts ... ') 
    px_cnt <- c()
    for(covname in covs){
      dimtoget <- max(dim(acl[[covname]])[3]) # grab just the most recent (or synoptic) layer, assuming they are all the same
      px_cnt   <- c(px_cnt, length(cellIdx(acl[[covname]][[dimtoget]]))) 
    }
    threshdrops <- covs[px_cnt/max(px_cnt) < check_pixelcount_thresh]
    if(length(threshdrops) == 0){
      message('  All clear')
    } else {
      message(sprintf('WARNING: the following covariates had < %i per cent pixel coverage in this region and are being dropped:\n %s', 
                      round(check_pixelcount_thresh,2)*100, paste(threshdrops,collapse=', ')))
      dropcovs <- c(dropcovs, threshdrops)
    }
  }
  
  # make sure dropcovs is unique
  dropcovs <- unique(dropcovs)
  
  # if any non-varying covariates were detected, then remove them from any objects
  if(length(dropcovs) > 0){
    for(dc in dropcovs){
      
      # drop from the cs_covs object, which itself is a list
      cc$covs[, (dc) := NULL]
      cc$cs_df <- cc$cs_df[!as.character(cc$cs_df$name) %in% dc,]
      
      # drop from the covariate layers list of rasters
      acl <- replace(acl, dc, NULL)
    }
  } else {
    message('No non-varying covariates were detected in the data. Yay!')
  }

  # redo the all fixed effects pseudo-formula strings
  fe   <- paste0(format_covariates(fe)  [!format_covariates(fe)   %in% dropcovs], collapse = ' + ')
  afe  <- paste0(format_covariates(afe) [!format_covariates(afe)  %in% dropcovs], collapse = ' + ')
  afeb <- paste0(format_covariates(afeb)[!format_covariates(afeb) %in% dropcovs], collapse = ' + ')
  tc   <- tc[!tc %in% dropcovs]
  
  # return the objects as a named list that will later get assigned to the environment within the parallel_model.R script
  return(list( cs_covs               = cc,
               all_fixed_effects     = afe,
               all_fixed_effects_brt = afeb,
               fixed_effects         = fe,
               the_covs              = tc,
               all_cov_layers        = acl))
}



