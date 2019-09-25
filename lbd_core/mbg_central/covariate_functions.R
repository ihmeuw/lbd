## Load and crop covariates to the modeled area
load_and_crop_covariates <- function(fixed_effects, simple_polygon, agebin=1) {

  # get selected covs
  selected_covs <- strsplit(fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]

  # Hard-code central directories
  central_cov_dir <- "<<<< FILEPATH REDACTED >>>>"
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
load_gbd_covariates = function(covs, measures, year_ids, age_ids, template, use_subnationals = F, simple_polygon, interval_mo){

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
  world_shape <- readRDS('"<<<< FILEPATH REDACTED >>>>"')

  # If we are not using subnationals, keep only national units; otherwise remove the parent units
  if (!use_subnationals) {
    world_shape <- world_shape[world_shape$level==3,]
  } else {
    world_shape <- world_shape[!world_shape$loc_id %in% unique(world_shape$parent_id),]
  }

  world_shape <- crop(world_shape, template)
  afras <- rasterize(world_shape, template, 'GAUL_CODE', fun = 'last')

  # Check to make sure all of the relevant gauls are in the file
  if (!all(world_shape$GAUL_CODE %in% unique(afras))) {
    afras <- raster::rasterize(world_shape[!world_shape$GAUL_CODE %in% unique(afras),], afras, 'GAUL_CODE', fun = 'first', update = T)
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
  covdir <- "<<<< FILEPATH REDACTED >>>>"

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
    #print(covlist[[l]])
    covlist[[l]]  <- crop(covlist[[l]], extent(simple_polygon))
    #  print(covlist[[l]])
    covlist[[l]]  <- setExtent(covlist[[l]], simple_polygon) #, keepres=TRUE, snap=TRUE)
    #print(covlist[[l]])
    covlist[[l]]  <- raster::mask(covlist[[l]], simple_polygon)
  }

  return(covlist)

}
