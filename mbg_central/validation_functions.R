
##################################################################################################
##################################################################################################
## Overview:  Will search saved results and model objects to produce m draws of estimates
#             for each holdout area. Note this is written as a unit (ie for one reg, holdout),
#             an will likely be run in a loop or apply function
#
## Inputs:
#          indicator_group: you know
#          indicator: oh you know
#          run_date: you definitely know
#          reg: string, region abbreviation, or just 'africa'
#          addl_strat: ie c('age'=1), make age=0 for others
#          holdout: holdout number
#
## Outputs: A table with one row per area-holdout, with p as estimated from data
#           and m draws of p_hat for each area
##################################################################################################
##################################################################################################

aggregate_validation <- function(holdoutlist = stratum_qt,
                                 cell_draws_filename = '%s_cell_draws_eb_bin%i_%s_%i%s.RData',
                                 years = c(2000,2005,2010,2015),
                                 indicator_group,
                                 indicator,
                                 run_date,
                                 reg,
                                 holdout,
                                 vallevel = "",
                                 addl_strat = c('age'=0),
                                 sbh_wgt = NULL,
                                 returnaggonly=TRUE,
                                 iter=i){

  require(raster); require(data.table)

  ## Load data
  datdir <- sprintf('<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/',indicator_group,indicator,run_date)

  # cell draws
  cdfile <- sprintf(paste0(datdir,cell_draws_filename),indicator,addl_strat,reg,holdout,vallevel)
  load(cdfile) # loads cell pred

  # holdout data
  d <- data.frame(holdoutlist[[sprintf('region__%s___%s__%i',reg,names(addl_strat),addl_strat)]])
  if(holdout!=0) {
    d_oos      <- d[d$fold==holdout,]
  } else {
    d_oos      <- d
    d_oos$fold = 0
  }
  # Acct for weights
  d_oos$indic_orig <-   d_oos[,indicator]
  if(is.null(sbh_wgt)){
    d_oos$exposure       <- d_oos$N*d_oos$weight
    d_oos[,indicator]    <- d_oos[,indicator] *d_oos$weight
  } else { # either a variable set (char) or numeric auto down weight
    if(class(sbh_wgt)=='character'){
      d_oos$exposure     <- d_oos$N*d_oos$weight*d_oos[,sbh_wgt]
      d_oos[,indicator]  <- d_oos$died*d_oos$weight*d_oos[,sbh_wgt]
    }
    if(class(sbh_wgt)=='numeric'){
      d_oos$exposure      <- d_oos$N*d_oos$weight*sbh_wgt
      d_oos[,indicator]   <- d_oos[,indicator] *d_oos$weight*sbh_wgt
    }
  }


  draws=dim(cell_pred)[2]

  # load simple raster
  load(paste0('<<<< FILE PATH REDACTED >>>>/simple_raster',reg,'.RData'))

  # for each draw in cell_draws, pull the p
  r_list = list()
#  pb <- txtProgressBar(style = 3)
  for(i in 1:draws){
#    setTxtProgressBar(pb, i / draws)
    r_list[[length(r_list)+1]] <- insertRaster(simple_raster, matrix(cell_pred[,i],ncol = length(years)))
  }
#  close(pb)

  # extract probabilities from draw rasters
  #message(paste('Extracting estimated probabilities at data locations for',draws,'draws rasters.'))
  ycol = match(d_oos$year,years)
  #pb <- txtProgressBar(style = 3)
  for(i in 1:draws){
  #  setTxtProgressBar(pb, i / draws)require(package, lib.loc = package_lib, character.only=TRUE)
    # extract at locations
    t=raster::extract(r_list[[i]],cbind(d_oos$longitude,d_oos$latitude))
    # match year and put into df
    d_oos[,paste0(indicator,'hat_',i)]=      sapply(seq_along(ycol), function(x){t[x,ycol[x]]}) * d_oos$exposure
    d_oos[,paste0(indicator,'hat_full_',i)]= sapply(seq_along(ycol), function(x){t[x,ycol[x]]}) * d_oos$N

  }
  #close(pb)

  # get binomial draws for each point-draw, and estimate the coverage based on X%ci
  # message('Doing coverage using binomial draws. ') # NAs happen if points were in water
  x <- matrix(NA,ncol=draws,nrow=nrow(d_oos))
  for(i in 1:draws){
    x[,i] <- rbinom(nrow(d_oos),
                    size=round(d_oos$N,0),
                    prob=d_oos[,paste0(indicator,'hat_full_',i)]/d_oos$N)
  }
  for(c in c(50,80,95)){
    coverage = c/100
    li=apply(x,1,quantile,p=(1-coverage)/2,na.rm=T)
    ui=apply(x,1,quantile,p=coverage+(1-coverage)/2,na.rm=T)
    d_oos[,paste0('clusters_covered_',c)] = d_oos[,'indic_orig']>=li & d_oos[,'indic_orig'] <= ui
  }

  # collapse into areas
  d_oos[,names(addl_strat)]<-addl_strat
  d_oos$total_clusters=1
  #message(names(d_oos))
  res <- d_oos[c(names(addl_strat),'fold','ho_id','year',indicator,'exposure',paste0('clusters_covered_',c(50,80,95)),'total_clusters',paste0(indicator,'hat_',1:draws))]
  f   <- as.formula(paste('.~ho_id+year+fold',names(addl_strat),sep='+'))
  resfull <- copy(res)
  res   <- aggregate(f,data=res,FUN=sum)

  # transform back to probability
  res$p = res[,indicator]/res$exposure
  res$region = reg
  res$oos = res$fold != 0
  for(i in 1:draws){
    res[,paste0('phat_',i)]=res[,paste0(indicator,'hat_',i)]/res$exposure
    res[,paste0(indicator,'hat_',i)]=NULL
  }
  res[,indicator]=NULL

  message(sprintf('%i is finished in function as %s',iter,paste0(class(res),collapse=' ')))


  ## return results table
  if(returnaggonly==TRUE){
    return(data.table(res))
  } else {
    # transform back to probability
    resfull$p = resfull[,indicator]/resfull$exposure
    resfull$region = reg
    resfull$oos = resfull$fold != 0
    for(i in 1:draws){
      resfull[,paste0('phat_',i)]=resfull[,paste0(indicator,'hat_',i)]/resfull$exposure
      resfull[,paste0(indicator,'hat_',i)]=NULL
    }
    resfull[,indicator]=NULL
    return(list(agg=data.table(res),cluster=data.table(resfull)))
  }

}


aggregate_validation_dev <- function( holdoutlist = stratum_qt,
                                      cell_draws_filename = '%s_cell_draws_eb_bin%i_%s_%i%s.RData',
                                      years = c(2000,2005,2010,2015),
                                      indicator_group,
                                      indicator,
                                      run_date,
                                      reg,
                                      holdout,
                                      vallevel = "",
                                      addl_strat = c('age'=0),
                                      iter=i){

  require(raster); require(data.table)

  ## Load data
  datdir <- sprintf('<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/',indicator_group,indicator,run_date)

  # cell draws
  cdfile <- sprintf(paste0(datdir,cell_draws_filename),indicator,addl_strat,reg,holdout,vallevel)
  #message(paste0('Loading cell draws from ',cdfile))
  load(cdfile)

  # holdout data
  d          <- as.data.table(holdoutlist[[sprintf('region__%s',reg)]])
  if(holdout!=0) {
    d_oos      <- d[d$fold==holdout,]
  } else {
    d_oos      <- d
    d_oos$fold = 0
  }
  d_oos <- d_oos[, exposure := N * weight]
  d_oos <- d_oos[, (indicator) := get(indicator) * weight]
  draws=dim(cell_pred)[2]
  # load simple raster
  load(paste0('<<<< FILE PATH REDACTED >>>>/simple_raster',reg,'.RData'))

  # for each draw in cell_draws, pull the p
  #message(paste('Pulling estimated rates for',dim(cell_pred)[2],'draws into rasters.'))
  r_list = list()
  #  pb <- txtProgressBar(style = 3)
  for(i in 1:draws){
    #    setTxtProgressBar(pb, i / draws)
    r_list[[length(r_list)+1]] <- insertRaster(simple_raster, matrix(cell_pred[,i],ncol = length(years)))
  }
  #  close(pb)

  # extract probabilities from draw rasters
  #message(paste('Extracting estimated probabilities at data locations for',draws,'draws rasters.'))
  ycol = match(d_oos$year,years)
  #pb <- txtProgressBar(style = 3)
  for(i in 1:draws){
    #  setTxtProgressBar(pb, i / draws)require(package, lib.loc = package_lib, character.only=TRUE)
    # extract at locations
    t=raster::extract(r_list[[i]],cbind(d_oos$longitude,d_oos$latitude))
    # match year and put into df
    d_oos[,paste0(indicator,'hat_',i)]= sapply(seq_along(ycol), function(x){t[x,ycol[x]]}) * d_oos$exposure
  }
  #close(pb)

  # get binomial draws for each point-draw, and estimate the coverage based on X%ci (95% for now)
  #message('Doing coverage using binomial draws. ')
  # get binomial draws for each point-draw, and estimate the coverage based on X%ci (95% for now)
  #message('Doing coverage using binomial draws. ')
  x=matrix(rbinom(nrow(d_oos)*100,size=round(d_oos$exposure,0),prob=d_oos[,get(paste0(indicator,'hat_',i))]/d_oos$exposure),ncol=100)
  for(c in c(50,80,95)){
    coverage = c/100
    li=apply(x,1,quantile,p=(1-coverage)/2,na.rm=T)
    ui=apply(x,1,quantile,p=coverage+(1-coverage)/2,na.rm=T)
    d_oos <- d_oos[, li := li]
    d_oos <- d_oos[, ui := ui]
    d_oos[get(indicator) >= li & get(indicator) <= ui, covered := 1]
    d_oos[is.na(covered), covered := 0]
    setnames(d_oos, 'covered', paste0('clusters_covered_',c))
  }


  # collapse into areas
  d_oos[,names(addl_strat)]<-addl_strat
  d_oos$total_clusters=1
  #message(names(d_oos))
  res <- d_oos[,c(names(addl_strat),'fold','ho_id','year',indicator,'exposure',paste0('clusters_covered_',c(50,80,95)),'total_clusters',paste0(indicator,'hat_',1:draws)), with=F]
  f   <- as.formula(paste('.~ho_id+year+fold',names(addl_strat),sep='+'))
  res   <- aggregate(f,data=res,FUN=sum)

  # transform back to probability
  res$p = res[,indicator]/res$exposure
  res$region = reg
  res$oos = res$fold != 0
  for(i in 1:draws){
    res[,paste0('phat_',i)]=res[,paste0(indicator,'hat_',i)]/res$exposure
    res[,paste0(indicator,'hat_',i)]=NULL
  }
  res[,indicator]=NULL

  message(sprintf('%i is finished in function as %s',iter,paste0(class(res),collapse=' ')))
  ## return results table
  return(data.table(res))

}

## ################################
## is_oos_preds
##
## this function takes the data that went into our MBG framework and
## returns an array containing: location, observed value, in sample
## predictions, and out of sample predictions (when applicable), country, year, N
## #################################

is_oos_preds_testing <- function(rd = run_date,
                         all.data = df,
                         cell_draws_filename = '%s_cell_draws_eb_bin%i_%s_%i.RData', ## in sprintf notation
                         holdouts = 5, ## number of holdouts. if zero only does in sample
                         reg,
                         years = 2000:2015,
                         indic = indicator,
                         indic_group = indicator_group,
                         holdoutlist = NULL ## if null, only does in sample
                         ){

  ## place to look for things
  output.dir <- sprintf("<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/",
                        indic_group, indic, rd)

  ## Load data
  datdir <- sprintf('<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/',indic_group,indic,rd)

  ## holdout data
  if(!is.null(holdoutlist)){
    d <- data.frame(holdoutlist[[sprintf('region__%s', reg)]])
  }else{
    d <- df
  }

  ## load the simple raster for this region
  load(paste0('<<<< FILE PATH REDACTED >>>>/simple_raster',reg,'.RData'))

  ## setup the data in the data order to return
  if(holdouts == 0){
    return.df <- d
  }else{
    return.df <- d[d$fold == 1, ]
    for(i in 2:holdouts){
      return.df <- rbind(return.df, d[d$fold == i, ])
    }
  }
  return.df <- return.df[, c('longitude', 'latitude', 'year', 'country', indic, 'N')]
  return.df$OOS <- NA
  return.df$IS  <- NA

  ## ####################
  ## get out of sample ##
  ## ####################

  OOS <- NULL

  ## loop through the holdouts
  if(holdouts!=0) {
    for(hh in 1:holdouts){

      ## load the associated preds
      load(sprintf(paste0('<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/', cell_draws_filename),
                   indic_group, indic, rd, indic, 0, reg, hh))

      ## average across the draws
      mean.cell.pred <- rowMeans(cell_pred)

      ## turn into a raster
      temp.rast <- insertRaster(simple_raster, matrix(mean.cell.pred, ncol = length(years)))

      ## get the OOS part of the data
      d.oos <- d[d$fold == hh, ]

      ## extract the values at the OOS locations by year
      temp.oos <- numeric(nrow(d.oos))
      for(yr in years){
        yr.rows <- which(d.oos$year == yr)
        if(length(yr.rows) > 0){
          temp.oos[yr.rows] <- raster::extract(y = cbind(d.oos$longitude, d.oos$latitude)[yr.rows, ],
                                               x = temp.rast[[which(years == yr)]])
        }
      }

      ## add to OOS vec
      OOS <- c(OOS, temp.oos)

    }

    return.df$OOS <- OOS

  }

  ## ################
  ## get in sample ##
  ## ################

  ## regardless of whether holdouts==0 or not, we can do in sample extraction
  hh <- 0
  d.is <- d
  d.is$fold <- 0

  ## load the associated preds
  load(sprintf(paste0('<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/', cell_draws_filename),
               indic_group, indic, rd, indic, 0, reg, hh))

  ## average across the draws
  mean.cell.pred <- rowMeans(cell_pred)

  ## turn into a raster
  temp.rast <- insertRaster(simple_raster, matrix(mean.cell.pred, ncol = length(years)))

  ## extract the values at the OOS locations
  is <- numeric(nrow(d.is))
  for(yr in years){
    yr.rows <- which(d.is$year == yr)
    if(length(yr.rows) > 0){
      is[yr.rows] <- raster::extract(y = cbind(d.is$longitude, d.is$latitude)[yr.rows, ],
                                           x = temp.rast[[which(years == yr)]])
      }
  }

  return.df$IS <- is

  ## return data with IS and OOS columns
  return(return.df)

}

## ###################################################################
## get_is_oos_draws()
##
## this function makes a dataframe of all your data that went into the
## model and appends columns of draw values from cell_preds
##
## INPUT
##
##   ind_gp: indicator_group
##   ind:    indicator
##   rd:     run_date
##   ind_fm: indicator_family (binomial, gaussian)
##   model_domain: larger domain you modelled over (e.g. africa even
##     if you use subregions)
##   nperiod: number of periods/years in model
##   years: vector of years. should be of length nperiod
##   get.oos: should we also get OOS extracted values?
##   write.to.file: if true writes final df to standard output dir
##
## OUTPUT
##
##   a data frame with:
##     nrow = number data observations
##     columns for each draw and some identifying columns:
##       holdout_id, value, sample size, region
##
## #####################################################################
get_is_oos_draws <- function(ind_gp,
                             ind,
                             rd,
                             ind_fm = 'binomial',
                             model_domain = 'africa',
                             age = 0,
                             nperiod = 16,
                             yrs = 2000:2015,
                             write.to.file = FALSE,
                             get.oos = FALSE,
                             year_col = 'original_year'
                             ){

  ## ##################################################################
  ## load in data, get domain templates and add on necessary columns ##
  ## ##################################################################

  mod.dir <- sprintf('<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/', ind_gp, ind, rd)

  message('loading in model domain templates')
  ## load in some shape data for the model_domain
  gaul_list <- get_gaul_codes(model_domain)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                             buffer = 0.4)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

  ## load in regions that were used in modeling
  all.regions <- get_output_regions(mod.dir)

  ## load raw data from modelling domain
  message('loading in input data used in model')
  if(!get.oos){
    df <- fread(paste0(mod.dir, 'input_data.csv')) ## raw input data

    ## add on gaul codes
    indicator <- ind
    Regions <- get_output_regions(in_dir = paste0('<<<< FILE PATH REDACTED >>>>/', ind_gp, '/', ind, '/output/', rd))
    df <- add_gauls_regions(df = df,
                            simple_raster = simple_raster)

    ## these regions may be wrong if people used custom regions...
    df$region <- NULL

    ## make a matrix of gauls and regions to merge onto df
    region.gaul <- data.frame(matrix(ncol = 2, nrow = 0))
    colnames(region.gaul) <- c('region', 'gaul')
    for(rr in all.regions){
      rr.gauls <- get_gaul_codes(rr)
      region.gaul <- rbind(region.gaul,
                           data.frame(region = rep(rr, length(rr.gauls)),
                                      gaul = rr.gauls))
    }
    region.gaul$region <- as.character(region.gaul$region)

    ## merge on region
    df <- merge(df, region.gaul, by.x = "GAUL_CODE", by.y = 'gaul', all.x = TRUE)

  }else{
    df <- as.data.table(do.call(rbind, readRDS(paste0(mod.dir, 'stratum.rds')))) ## holdoutlist
  }

  ## rename year column for convenience
  setnames(df, year_col, "the_year_col")

  ## make year to period key
  yr.per <- cbind(yrs, 1:nperiod)
  colnames(yr.per) <- c('yr', 'per')

  ## load in admin1 and 2 relevant shapefiles
  message('loading in shapefiles to determine ad1 and ad2 membership')
  ## Load admin2 raster
  admin_level <- 2
  shapes <- readRDS("<<<< FILE PATH REDACTED >>>>/g_2015_2014_2_modified/g_2015_2014_2_modified.rds")
  
  ## Fix several african disputed territories (assign them to a
  ## country so we can include them in our estimates and they dont
  ## drop out)
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==61013]=133
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40760]=40765
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40762]=145
  shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==102  ]=6
  admin2_shapefile <- shapes

  locs <- SpatialPoints(cbind(df$longitude, df$latitude), proj4string = CRS(proj4string(admin2_shapefile)))
  adm.df <- sp::over(locs, admin2_shapefile)
  df$ad1 <- adm.df$ADM1_CODE
  df$ad2 <- adm.df$ADM2_CODE
  ## set up ad0 for pvtables and consistency
  df$ad0 <- df$country

  ## ######################
  ## get in sample draws ##
  ## ######################

  if(get.oos){ ## save a copy of the df for oos if we're doing oos
    df.oos <- df
  }

  ## now, for each region we need to load in draws, make a raster,
  ## extract values and insert them into the df

  message('~~~~~GETTING IS DRAWS~~~~~')
  for(rr in all.regions){

    message(sprintf('loading cell preds for: %s', rr))
    load(paste0(mod.dir, sprintf('%s_cell_draws_eb_bin%i_%s_0.RData', ind, age, rr)))

    n.pred <- ncol(cell_pred)

    ## load the simple raster template
    gaul_list <- get_gaul_codes(rr)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                               buffer = 0.4)
    subset_shape   <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
    raster_list    <- build_simple_raster_pop(subset_shape)
    template       <- raster_list[['simple_raster']]
    ## pop_raster     <- raster_list[['pop_raster']]

    ## get rows in region
    df.r  <- which(df$region == rr)
    loc.r <- cbind(df$longitude, df$latitude)[df.r, ]

    ## make a matrix to fill with draws
    draw.mat <- matrix(ncol = ncol(cell_pred), nrow = length(df.r))

    ## getting draws
    message('::extracting cell preds')
    ## make a matrix full of ids to easily pull the correct rows
    idx.r <- 1:nrow(cell_pred)
    idx.rast <- insertRaster(template,
                             matrix(idx.r,
                                    ncol = nperiod))

    ## get ids for data observations
    for(yy in 1:nperiod){

      ## extract data by year
      yr.r <- which(df[df.r, ]$the_year_col == yrs[yy])
      idx.yy <- raster::extract(idx.rast[[yy]], matrix(loc.r[yr.r, ], ncol = 2))

      ## enter draws into draw.mat
      draw.mat[yr.r, ] <- cell_pred[idx.yy, ]
    }

    if(which(all.regions == rr) == 1){
      domain.draws <- matrix(ncol = ncol(cell_pred), nrow = nrow(df))
      colnames(domain.draws) <- paste0("draw", 1:ncol(cell_pred))
    }
    domain.draws[df.r, ] <- draw.mat
  }

  df$fold <- 0
  df <- cbind(df, domain.draws) ## now we have fold 0, i.e. in sample draws

  ## Tack on draws of hyperparameters needed for predictive draws to
  ## make coverage later (i.e., tau for Gaussian models)
  if(ind_fm=='gaussian') {
    pull_region_taus <- function(rr) {
      load(paste0(mod.dir, 'inla_draws/inla_draws_', rr, '.RData'))
      hyper_name <- 'Log precision for the Gaussian observations -- in user scale'
      taus <- c()
      for(i in 1:length(draws)) taus[i] <- draws[[i]][['hyperpar']][[hyper_name]]
      taus <- data.table(tau=taus, index=paste0('tau_', 1:length(taus)), region=rr)
      taus <- dcast(taus, region ~ index, value.var = "tau")
      return(taus)
    }
    all_taus <- rbindlist(lapply(all.regions, pull_region_taus))
    df <- merge(df, all_taus, by='region')
  }

  ## ##########
  ## get OOS ##
  ## ##########

  if(get.oos){
    message('~~~~~GETTING OOS DRAWS~~~~~')
    ## now we take the df.oos and do the same thing but load cell
    ## preds by holdout

    domain.draws <- matrix(ncol = ncol(cell_pred), nrow = nrow(df.oos))
    ## we should already have a cell_preds loaded since IS happens first
    colnames(domain.draws) <- paste0("draw", 1:ncol(cell_pred))

    all.holds <- sort(unique(df.oos$fold))

    for(rr in all.regions){
      message(sprintf('loading cell preds for region: %s', rr))

      ## load the simple raster template
      gaul_list <- get_gaul_codes(rr)
      simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                                 buffer = 0.4)
      subset_shape   <- simple_polygon_list[[1]]
      simple_polygon <- simple_polygon_list[[2]]
      raster_list    <- build_simple_raster_pop(subset_shape)
      template       <- raster_list[['simple_raster']]

      df.oos.r  <- which(df.oos$region == rr)
      loc.r <- cbind(df.oos$longitude, df.oos$latitude)[df.oos.r, ]

      ## make a matrix to fill with draws for the region
      draw.mat <- matrix(ncol = ncol(cell_pred), nrow = length(df.oos.r))

      for(hh in all.holds){
        message(sprintf('::loading cell preds for fold: %i', hh))

        load(paste0(mod.dir, sprintf('%s_cell_draws_eb_bin%i_%s_%i.RData', ind, age, rr, hh)))

        if(hh == 1){
          idx.r <- 1:nrow(cell_pred)
          idx.rast <- insertRaster(template,
                               matrix(idx.r,
                                      ncol = nperiod))
        }

        ## get rows in holdout
        df.oos.rh  <- which(df.oos[df.oos.r, fold] == hh) ## index within region in our holdout
        loc.rh <- loc.r[df.oos.rh, ] ## locations for region and holdout
        df.oos.rh.idx <- df.oos.r[df.oos.rh] ## rows in df.oos for region and holdout

        ## get ids for data observations
        message('::::extracting cell preds across years')
        for(yy in 1:nperiod){

          ## extract data by year
          yr.rh <- which(df.oos[df.oos.rh.idx, ]$the_year_col == yrs[yy]) ## index within region and holdout for yy
          idx.yy <- raster::extract(idx.rast[[yy]], matrix(loc.rh[yr.rh, ], ncol = 2))

          ## enter draws into draw.mat
          draw.mat[df.oos.rh[yr.rh], ] <- cell_pred[idx.yy, ] ## draws in region rr,
                                                              ## and in holdout hh,
                                                              ## and in year yy
        }
      }
      domain.draws[df.oos.r, ] <- draw.mat ## draws for a region
    }

    df.oos <- cbind(df.oos, domain.draws)

    ## stick together is and oos

    df <- rbind(df, df.oos)
  }

  # Rename back to original year column name
  setnames(df, "the_year_col", year_col)
  
  if(write.to.file){
    write.csv(df, paste0(mod.dir, 'output_draws_data.csv'))
  }
  return(df)

}

## ####################################
## a function to plot raking factors ##
## ####################################
## this function plots a scatterplot of GBD estimates vs MBG estimates
plot.rfs <- function(ind.gp = indicator_group,
                     ind = indicator,
                     rd = run_date,
                     output.dir = si.fig.dir,
                     title = "Comparison to GBD 2016 in\n" ## region gets pasted on to this
){

  require(ggrepel) ## for labels

  regions <- get_output_regions(in_dir = paste0('<<<< FILE PATH REDACTED >>>>/',
                                              ind.gp, '/',
                                              ind, '/output/',
                                              rd))

  for(rr in regions){

    ## convert rrs to full names
    if(rr == 'essa') rr_name = "Eastern Sub-Saharan Africa"
    if(rr == 'wssa') rr_name = "Western Sub-Saharan Africa"
    if(rr == 'name') rr_name = "North Africa"
    if(rr == 'sssa') rr_name = "Southern Sub-Saharan Africa"
    if(rr == 'cssa') rr_name = 'Central Sub-Saharan Africa'

    in_dir  <- paste0('<<<< FILE PATH REDACTED >>>>/', ind.gp, '/', ind, '/output/', rd)
    default_rf_path <- paste0(in_dir, '/', ind, '_rf.csv')
    all_rfs <- fread(default_rf_path)
    gaul_list = get_gaul_codes(rr)
    rfs <- all_rfs[name %in% gaul_list, ]
    loc_names <- setDT(read.csv("<<<< FILE PATH REDACTED >>>>/gaul_to_loc_id.csv",stringsAsFactors = F))
    setnames(rfs, "name", "GAUL_CODE")
    rfs <- merge(rfs, loc_names, by="GAUL_CODE")
    rfs[, Year:= as.factor(year)]
    max_val = max(max(rfs[,.(rake_to_mean, geo_mean)],na.rm=T),na.rm= T)

    ## plot w/o country labels
    gg_rfs <- ggplot(data = rfs, aes(x = rake_to_mean, y = geo_mean)) +
    geom_point(aes(color = Year)) +
    ylab("MBG Mean") +
    xlab("GBD Mean") +
    theme_bw() +
    xlim(0, max_val) +
    ylim(0, max_val) +
    geom_abline(slope = 1) +
    ggtitle(paste0(title, rr_name))

    assign(sprintf('%s_rf', rr), gg_rfs)

    ## plot w/ country labels
    gg_rfs <- ggplot(data = rfs, aes(x = rake_to_mean, y = geo_mean)) +
    geom_point(aes(color = Year)) +
    geom_text_repel(aes(label = ihme_lc_id),
                    segment.color = 'grey80') +
    ylab("MBG Mean") +
    xlab("GBD Mean") +
    theme_bw() +
    xlim(0, max_val) +
    ylim(0, max_val) +
    geom_abline(slope = 1) +
    ggtitle(paste0(title, rr_name))

    assign(sprintf('%s_rf_labs', rr), gg_rfs)

  }

  ## stick them all together
  require(gridExtra)
  margin = theme(plot.margin = unit(rep(.5, 4), "cm"))
  all.rfs <- grid.arrange(cssa_rf + margin,
                          essa_rf + margin,
                          name_rf + margin,
                          sssa_rf + margin,
                          wssa_rf + margin,
                          ncol=2)
  ggsave(filename = sprintf('%s%s_all_rfs.png',
                            output.dir, ind),
         all.rfs, width = 12, height = 16)

  ## stick them all together
  require(gridExtra)
  margin = theme(plot.margin = unit(rep(.5, 4), "cm"))
  all.rfs <- grid.arrange(cssa_rf_labs + margin,
                          essa_rf_labs + margin,
                          name_rf_labs + margin,
                          sssa_rf_labs + margin,
                          wssa_rf_labs + margin,
                          ncol=2)
  ggsave(filename = sprintf('%s%s_all_rfs_labs.png',
                            output.dir, ind),
         all.rfs, width = 12, height = 16)
}

## #################################################################
## ~~~~~~~~~~~ make table of INLA model results ~~~~~~~~~~~~~~~~~ ##
## #################################################################
## takes in standard model run inputs outputs a table of fixed effect,
## spatio-temporal hyperparameter, and random effects parameter
## summaries.
## note: takes a little while since it has to recreate the
## SPDE INLA object since neither we nor INLA saved that object

model_fit_table_list <- function(regions, rd=run_date, holdout = 0,
                                 age = 0,
                                 ind= indicator,
                                 ind_gp = indicator_group,
                                 sharedir = sprintf('<<<< FILE PATH REDACTED >>>>/%s/%s',indicator_group,indicator)){
  ## load models
  require(INLA)
  message(sprintf('Pulling together results for %s models',rd))

  tlist=list()

  for(rr in regions){
    message(sprintf('::on region %s',rr))
    reg  <-  rr

    message("::::loading in pre-INLA objects to get spde")
    pathaddin  <-  paste0('_bin',age,'_',rr,'_',holdout)
    load(paste0('<<<< FILE PATH REDACTED >>>>/', ind_gp, '/',
                ind, '/model_image_history/', rd, pathaddin,
                '.RData'))

    modnames = c('gam','gbm','ridge','enet','lasso')

    full_raster_list <- cov_list

    for(mm in modnames){
      if(min(na.omit(values(full_raster_list[[mm]][[1]]))) < 0){
        message(sprintf("un-logiting: %s", mm))
        full_raster_list[[mm]] <- ilogit(full_raster_list[[mm]])
      }
    }

    ## for stacking, overwrite the columns matching the model_names so
    ## that we can trick inla into being our stacker
    df = df[,paste0(child_model_names) := lapply(child_model_names,
                                                 function(x) get(paste0(x,'_cv_pred')))]

    ## Create SPDE INLA stack
    input_data <- build_mbg_data_stack(df = df,
                                       fixed_effects = all_fixed_effects,
                                       mesh_s = mesh_s,
                                       mesh_t = mesh_t,
                                       use_ctry_res = use_inla_country_res,
                                       use_nugget = use_inla_nugget)

    spde <- input_data[[2]]
    ## this is what we neede!

    message('::::loading in INLA fit\n')
    f <-  sprintf('%s/output/%s/inla_model_fit_pre_preds_%s_holdout_%i.RDS',
                         sharedir,rd,reg, holdout)
    res_fit <- readRDS(f)

    ## now we extract what we need from the fit to get transformed spatial params
    res.field <- inla.spde2.result(res_fit, 'space', spde, do.transf=TRUE)

    ## nominal range at 0.025, 0.5, 0.975 quantiles
    range   <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.range.nominal[[1]])
    nom.var <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.variance.nominal[[1]])
    spat.hyps <- rbind(range, nom.var)
    rownames(spat.hyps) <- c('Nominal Range', 'Nominal Variance')

    ## other hyperparmas
    hyps <- summary(res_fit)$hyperpar[-(1:2), ] ## first two rows are
                                                ## theta1, theta2 which
                                                ## we have in range and
                                                ## nom.var

    colnames(spat.hyps) <- colnames(hyps)[3:5]
    ## fixed effects
    fixed <- summary(res_fit)$fixed[,1:6]

    ## combine them all and just keep three quantiles

    all.res <- rbind(fixed[, 3:5],
                     spat.hyps,
                     hyps[, 3:5])
    tlist[[rr]] <- all.res
  }
  return(tlist)
}

## ############################################################
## ~~~~~~~~~~~function to plot residual errors ~~~~~~~~~~~~~ ##
## ############################################################
## takes in a gaul_list over the entire modelling domain (e.g. africa)
## takes in subset shape for entire modelling domain
## takes in the data.frame that get_is_oos_draws() outputs
## save.dir determines output path
##
## saves plots to output dir and also returns them

plot_abs_errors <- function(gaul_list = gaul_list,
                            df, ## takes output from get_is_oos_draws()
                            sample = 'BOTH',## sample == "IS" or "OOS", or "BOTH"
                            subset_shape = subset_shape,
                            ind = indicator,
                            ind_gp = indicator_group,
                            rd = run_date,
                            save.dir,
                            year_col = 'original_year') {

  if(ind == 'wasting_mod_b') nice.name = "Wasting"
  if(ind == 'stunting_mod_b') nice.name = "Stunting"
  if(ind == 'underweight_mod_b') nice.name = "Underweight"

  ## setup the dataframe
  subset_shape <- subset_shape[subset_shape$GAUL_CODE %in% gaul_list, ]

  ## rename year col for convenience
  df <- copy(as.data.table(df))
  setnames(df, year_col, "the_year_col")

  ## calculate residual: count/N - pred
  ## calculate residual: count/N - pred
  phat <- base::rowMeans(draws.df[, grep('draw', colnames(df)), with = FALSE], na.rm = TRUE)
  phat[is.nan(phat)] <- NA
  df$phat <- phat
  df$pobs <- df[[ind]] / df[['N']]
  df$abs_error =  df$pobs - df$phat
  df <- subset(df, !is.na(abs_error))

  if(sample == 'IS')   to.do <- c(1, 0)
  if(sample == 'OOS')  to.do <- c(0, 1)
  if(sample == 'BOTH') to.do <- c(1, 1)


  full.df <- df
  if(to.do[1] == 1){ ## is

    df <- subset(full.df, fold == 0)

    if(length(df[, GAUL_CODE]) != 0) {
      this_shape.dt <- data.table(fortify(subset_shape))
      redwhiteblue <- c(scales::muted('blue'),
                         'white',
                         scales::muted('red'))
      ## plot gg
      gg.is <- ggplot(df, aes(longitude, latitude)) +
      geom_point(aes(color=abs_error,
                     size = N,
                     alpha = weight)) +
      coord_fixed() +
      geom_path(data=this_shape.dt, aes(x=long, y=lat, group=group),
                color='black', lwd=.1) +
      scale_color_gradientn(colours=redwhiteblue,
                            values=c(-1,0,1), limits=c(-1,1),
                            na.value = "#000000",
                            rescaler = function(x, ...) x,
                            oob = identity) +
      guides(color=guide_colorbar(title="Absolute\nerror",
                                  label=TRUE,
                                  ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      theme(panel.margin = unit(0,"lines"),
            plot.margin = unit(c(0,0,0,0),"lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~the_year_col) +
      ggtitle(paste0(nice.name, ' absolute error'))

      ggsave(filename = sprintf('%s%s_abs_error_plot_IS.png', save.dir, ind),
             plot = gg.is, width = 12, height = 12, units = 'in')
    }else{
      gg.is <- NULL
    }
  }

  if(to.do[2] == 1){ ## oos

    df <- subset(full.df, fold != 0)

    if(length(df[, GAUL_CODE]) != 0) {
      this_shape.dt <- data.table(fortify(subset_shape))
      redwhiteblue <- c(scales::muted('blue'),
                        'white',
                        scales::muted('red'))
      ## plot gg
      gg.oos <- ggplot(df, aes(longitude, latitude)) +
      geom_point(aes(color=abs_error,
                     size = N,
                     alpha = weight)) +
      coord_fixed() +
      geom_path(data=this_shape.dt, aes(x=long, y=lat, group=group),
                color='black', lwd=.1) +
      scale_color_gradientn(colours=redwhiteblue,
                            values=c(-1,0,1), limits=c(-1,1),
                            na.value = "#000000",
                            rescaler = function(x, ...) x,
                            oob = identity) +
      guides(color=guide_colorbar(title="Absolute\nerror",
                                  label=TRUE,
                                  ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      theme(panel.margin = unit(0,"lines"),
            plot.margin = unit(c(0,0,0,0),"lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~the_year_col) +
      ggtitle(paste0(nice.name, ' absolute error'))

      ggsave(filename = sprintf('%s%s_abs_error_plot_OOS.png', save.dir, ind),
             plot = gg.oos, width = 12, height = 12, units = 'in')

    }else{
      gg.oos <- NULL
    }
  }

  return(list(is = gg.is,
              oos = gg.oos))
}




###########
## PLOT OF PLOTS FUNCTION FOR TESTING AGGREGATED OOS RESULTS IN STACKING
plot_of_plots <- function(run_date,
                          reg,
                          age,
                          nfolds=5) {

  # load in neeed data
  folds <- c()
  for(i in 1:nfolds){
    if(file.exists(sprintf('%s/output/%s/fin_bin%i_%s_%i',sharedir,run_date,age,reg,i))){
      folds <- c(folds,i)
    }
  }
  message(paste('The Following Folds (of',nfolds,') are ready and will be used:'))
  message(paste(folds,collapse=', '))

  if(length(folds)>0){
    d <- data.table()
    message('\nLoading data:')
    for(i in folds){
      message(i)
      patt <- sprintf('%s_aggval_bin%i_%s_%i',indicator,age,reg,i)
      for(f in list.files(path   =sprintf('%s/output/%s/',sharedir,run_date),pattern=patt)){
        message(f)
         tmp       <- fread(sprintf('%s/output/%s/%s',sharedir,run_date,f))
         tmp$model <- gsub('.csv','',gsub(paste0(patt,'_'),'',f))
         d         <- rbind(d,tmp)
       }
    }

    # simplify d for now
    d <- d[,c('p','mean','error','model','ho_id','year','fold','age','exposure','clusters_covered_95', 'total_clusters'),with=F]

    # mean error
    me  <- aggregate(error~model,d,mean)
    names(me) <- c('model','mean_error')

    # RMSE
    rmse <- aggregate(error~model,d,function(x){sqrt(mean(x^2))})
    names(rmse) <- c('model','rmse')

    # correlation
    require(plyr)
    corr <- ddply(d,"model",function(x) cor(x$p,x$mean))
    names(corr) <- c('model','correlation')

    res <- data.frame(model=corr$model)
    for(pv in c('me','rmse','corr'))
      res <-merge(res,get(pv),by='model')



    # plot
    require(ggplot2); require(grid); require(gridExtra)
    get_legend<-function(myggplot){
      tmp <- ggplot_gtable(ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)
    }
    p1 <- ggplot(res,aes(x=correlation,y=rmse,colour=model,shape=model))+geom_point()+
                    theme_bw()+
                    theme(legend.position="none")
    p2 <- ggplot(res,aes(x=correlation,y=mean_error,colour=model,shape=model))+geom_point()+
                    theme_bw()+
                    theme(legend.position="none")
    p3 <- ggplot(res,aes(x=rmse,y=mean_error,colour=model,shape=model))+geom_point()+
                    theme_bw()
    legend <- get_legend(p3)
    p3 <- p3 + theme(legend.position="none")

    message(paste0('saving ',sprintf('%s/output/%s/plot_of_plots_bin%i_%s.pdf',sharedir,run_date,age,reg)))
    pdf(sprintf('%s/output/%s/plot_of_plots_bin%i_%s.pdf',sharedir,run_date,age,reg))
      grid.arrange(p1, p2, p3, legend, ncol = 2, top = sprintf('OOS Predictive Validity: Age Bin %i, %s. %i/%i folds analyzed',age,reg,length(folds),nfolds))
    dev.off()
 } else {
   message('NO DATA WRITTEN YET, SO NO PLOTS')
 }
 #return(plot)
}

####################################################################################################################
## INPUT: # data file matching <<<< FILE PATH REDACTED >>>>/child_growth_failure/wasting_mod_b/output/2017_06_29_01_24_15/output_draws_data.csv
## OUTPUT SUMMARIZED PREDICTIVE VALIDITY METRICS
##
####################################################################################################################

get_pv_table <- function(d,
                        indicator,
                        indicator_group,
                        rd,
                        aggregate_on, # column of spatial aggregation (admin1 or admin2?)
                        draws = 1000, # number of draws
                        coverage_probs = c(95), # coverage percentage
                        result_agg_over =  c('year','oos'), # final table aggregates over
                        weighted = TRUE, # weight PV metrics on SS?
                        family = 'binomial', # distribution
                        plot = TRUE,
                        plot_ci = FALSE,
                        plot_ci_level = 95,
                        ci_color = "grey", 
                        point_alpha = 1,
                        point_color = "black",
                        plot_title = indicator,
                        plot_ncol = 4,
                        save_csv = T,
                        out.dir) {

  require(boot)
  require(ggplot2)
  require(tidyr)

  d <- data.table(d)

  ## Get binomial predictions
  message('Simulating predictive draws')
  if (family == 'binomial') {
    x <- sapply(1:draws, function(i) rbinom(nrow(d), round(d$N), d[[paste0('draw', i)]]))
    d[, Y := get(indicator) * round(N)/N] # adjust observed number of cases for effect of rounding N
  }
  if (family == 'gaussian') {
    message('NICK: THIS SHOULD WORK IF YOU HAVE DRAWS OF PRECISION NAMED tau_1,...tau_1000 BUT THIS IS UNTESTED')
    x <- sapply(1:draws, function(i) rnorm(nrow(d), sqrt(d[[paste0('tau_', i)]]), d[[paste0('draw', i)]]))
    d[, Y := get(indicator)]
  }

  message(paste0('...NOTE: missing predictions for ', sum(is.na(x[,1])), ' of ', nrow(x), ' (', round(100*mean(is.na(x[,1]))), '%) data points.'))

  ## Get coverage for each data point
  message('Calculate coverage')
  for(c in coverage_probs){
    one_side <- (1 - c/100)/2
    ui <- t(apply(x, 1, quantile, c(one_side, 1 - one_side), na.rm = T))
    d[, paste0('clusters_covered_', c) := Y >= ui[,1] & Y <= ui[,2]]  ## NOTE THAT THIS DOES NOT ACCOUNT FOR ROUNDING N!
  }

  ## Collapse data and draws spatially
  message('Collapse data and draws spatially')
  d[, oos := (fold != 0)]
  d[, p := get(indicator) / N]
  d[, exposure := N * weight]

  by_vars <- unique(c('oos', 'year', result_agg_over, aggregate_on))
  collapse_vars <- c('p', paste0('clusters_covered_', coverage_probs), paste0('draw', 1:draws))

  res <- d[!is.na(draw1),
            c(list(total_clusters = .N, exposure = sum(exposure)),
              lapply(.SD, function(x) weighted.mean(x, exposure, na.rm=T))),
            keyby = by_vars, .SDcols = collapse_vars]
  res[, mean_draw := rowMeans(.SD), .SDcols = paste0('draw', 1:draws)]
  res[, error := p - mean_draw]
  res[, abs_error := abs(error)]

  ## Collapse to calculate predictive validity metrics
  message('Outputting predictive validity metrics for your models')
  weighted.rmse <- function(error, w) sqrt( sum(  (w/sum(w)) * ((error)^2) ) )
  if (weighted) res$weight <- res$exposure else res$weight <- 1
  res2 <- res[, c(lapply(.SD, function(x) weighted.mean(x, weight)),
                 rmse = weighted.rmse(error, weight),
                 median_SS = median(exposure),
                 cor = corr(cbind(p, mean_draw), weight)),
             by = result_agg_over,
             .SDcols = c("error", "abs_error", "mean_draw", "p", paste0("clusters_covered_", coverage_probs))]
  setnames(res2, c("error", "abs_error", "p", paste0("clusters_covered_", coverage_probs)),
                 c("me", "mae", "mean_p", paste0("coverage_", coverage_probs)))

  # Make plots
  if(plot==TRUE){
    message('Making plots of aggregated data and estimates')

    if (aggregate_on == "country") agg_title <- "Country"
    if (aggregate_on == "ad1") agg_title <- "Admin 1"
    if (aggregate_on == "ad2") agg_title <- "Admin 2"

    for (oosindic in unique(res$oos)) {

      if (!("region" %in% names(res))) res[, region := "all"] #allow looping over regions

      for (reg in unique(res$region)) {

        plot_file <- paste0(out.dir, indicator, '_validation_plot_', 
                            paste(c(aggregate_on, result_agg_over), collapse="_"), "_", 
                            ifelse(oosindic, "OOS", "IS"), 
                            ifelse(reg == "all", "", paste0("_", reg)),
                            '.png')
        message(paste('...saving plot here:', plot_file))
        png(plot_file, width = 12, height = 12, units = 'in', res = 350)

        fdata <- res[oos == oosindic,]
        if (reg != "all") fdata <- res[region == reg,]

        if (plot_ci) {
          fdata[, upper := apply(.SD, 1, quantile, p = 0.01*(plot_ci_level + (100 - plot_ci_level)/2), rm.na=TRUE), .SDcols = paste0('draw', 1:draws)]
          fdata[, lower := apply(.SD, 1, quantile, p = 0.01*((100 - plot_ci_level)/2), rm.na=TRUE), .SDcols = paste0('draw', 1:draws)]
          limits <- fdata[, range(c(p, mean_draw, lower, upper))]
        } else {
          limits <- fdata[, range(c(p, mean_draw))]
        }

        gg <- ggplot(fdata, aes(x = p, y = mean_draw, size = weight)) +
          geom_abline(intercept=0, slope=1, color = 'red') +
          geom_point(colour = point_color, alpha = point_alpha) +
          scale_size_area() +
          xlim(limits) +
          ylim(limits) +
          coord_equal() +
          theme_bw() +
          theme(strip.background = element_rect(fill="white")) +
          labs(x = 'Data Estimate', 
               y = 'Mean Prediction', 
               size = "Weight",
               title = paste0("Validation Plot for ", plot_title, " by ", agg_title),
               subtitle = paste0("OOS: ", oosindic, ifelse(is.null(reg), "", paste0(" | Region: ", reg))))

        if(plot_ci) {
          gg <- gg +  geom_errorbar(aes(ymin = lower, ymax = upper), colour = point_color, width = 0, size = .3, alpha = min(point_alpha, 0.2))
        }
        if (length(setdiff(result_agg_over, 'oos')) > 0) {
          gg <- gg + facet_wrap(as.formula(paste("~", paste(setdiff(result_agg_over, 'oos'), collapse = "+"))))
        }

        plot(gg)
        dev.off()

      }
    }

    if (plot == T & length(coverage_probs) > 1) {
      # Plot multiple coverage levels
      library(tidyr)
      library(ggplot2)
      str_match <- stringr::str_match

      # Assign a region if none present to allow looping over regions
      if (!("region" %in% names(res2))) res2[, region := "all"]
      for (reg in unique(res2$region)) {
        reg_addin <- ifelse(reg == "all", "", paste0("_", reg))
        res2_cov <- subset(res2, select = c("Year", "OOS", names(res2)[grepl("% Cov.", names(res2))]))
        res2_cov <- res2_cov %>% gather(coverage, observed_coverage, -Year, -OOS) %>% as.data.table
        res2_cov[, expected_coverage := as.numeric(str_match(coverage, "(.*)% Cov.")[,2])]
        res2_cov[, observed_coverage := observed_coverage * 100]

        reg_addin <- ifelse(is.null(reg), "", paste0("_", reg))
        filename <- paste0(out.dir, indicator, '_calibration_plot_', aggregate_on, reg_addin, '.png')
        png(filename, width = 12, height = 6, units = 'in', res = 350)

        gg <- ggplot(res2_cov, aes(x = expected_coverage, y = observed_coverage, group = Year, color = Year)) +
                geom_line(aes(x = expected_coverage, y = expected_coverage), color = "red", alpha = 0.2) +
                geom_point(aes(x = expected_coverage, y = expected_coverage), shape = 8, color = "red") +
                geom_point() +
                geom_line(alpha = 0.2) +
                theme_bw()+
                labs(x = "Expected coverage", y = "Observed coverage",
                     title = paste0("Validation Plot for ", plot_title, " by ", agg_title),
                     subtitle = ifelse(is.null(reg), "", paste0("Region: ", reg))) +
                facet_wrap(~OOS, ncol = plot_ncol) +
                scale_x_continuous(expand = c(0, 0), limits = c(0,100)) + 
                scale_y_continuous(expand = c(0, 0), limits = c(0,100)) + 
                coord_equal()   

        plot(gg)
        dev.off()
      }
    }
  }


  if (save_csv) {
    # Save final table if desired
    filename <- sprintf("<<<< FILE PATH REDACTED >>>>/%s/%s/output/%s/summary_metrics/%s_metrics%s.csv",
                         indicator_group, indicator, run_date, aggregate_on,
                         ifelse(("region" %in% result_agg_over), "_by_region", ""))
    message(paste0("Saving csv to ", filename, "..."))
    write.csv(res2, file = filename)
  }

  if (plot == T & length(coverage_probs) > 1) {
    message('Making coverage calibration plots')

    fdata <- res2[, unique(c("oos", result_agg_over, paste0("coverage_", coverage_probs))), with=F]
    fdata <- melt(fdata, id.vars = unique(c("oos", result_agg_over)), value.name = "observed_coverage", variable.name = "coverage")
    fdata[, expected_coverage := as.numeric(gsub("coverage_", "", coverage))]
    fdata[, observed_coverage := observed_coverage * 100]
    fdata$group <- apply(fdata[, setdiff(result_agg_over, 'oos'), with=F], 1, paste, collapse=" ")
    if (sum(!is.na(fdata$group)) == 0) fdata[, group := "All"]
    fdata[, oos := factor(oos, c(T, F), c("Out of Sample", "In Sample"))]

    png(paste0(out.dir, indicator, '_calibration_plot_', paste(c(aggregate_on, result_agg_over), collapse="_"), '.png'),
        width = 12, height = 12, units = 'in', res = 350)

    limits <- fdata[, range(c(observed_coverage, expected_coverage))]
    gg <- ggplot(fdata, aes(x = expected_coverage, y = observed_coverage, group = group, color = group)) +
      geom_abline(intercept = 0, slope = 1, color = "red") +
      geom_point() +
      geom_line(alpha = 0.2) +
      scale_color_discrete(name = "") +
      facet_wrap(~ oos, ncol = plot_ncol) +
      coord_equal() +
      xlim(limits) +
      ylim(limits) +
      theme_bw()+
      labs(x = "Expected coverage", y = "Observed coverage")

    plot(gg)
    dev.off()
  }

  # format
  setorderv(res2, result_agg_over)
  setnames(res2, result_agg_over, ifelse(result_agg_over == 'oos', 'OOS', gsub("(^.)", "\\U\\1", result_agg_over, perl=T)))
  setnames(res2, c("me", "mean_draw", "mean_p", "rmse", "median_SS", "cor", paste0("coverage_", coverage_probs)),
           c("Mean Err.", "Mean Pred.", "Mean Obs.", "RMSE", "Median SS", "Corr.", paste0(coverage_probs, "% Cov.")))

  return(res2)
}