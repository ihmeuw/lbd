
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
  datdir <- sprintf('/%s/',indicator_group,indicator,run_date)

  # cell draws
  cdfile <- sprintf(paste0(datdir,cell_draws_filename),indicator,addl_strat,reg,holdout,vallevel)
  #message(paste0('Loading cell draws from ',cdfile))
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
  load(paste0('<<<< FILEPATH REDACTED >>>>>/simple_raster',reg,'.RData'))

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
  datdir <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/',indicator_group,indicator,run_date)

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
  load(paste0('<<<< FILEPATH REDACTED >>>>>/simple_raster',reg,'.RData'))

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
    t=extract(r_list[[i]],cbind(d_oos$longitude,d_oos$latitude))
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

###################################################################################
## this function is meant to plot your raw data (or something as close
## to raw as possible) in a format that closely resembles your mbg
## model output for comparison
##
## Inputs -
##   data:    cleaned data.frame that will go into MBG run
##   *_col:   string name of relevant columns in data
##   smooth:  tuning for thin-plate spline smoothness
##   raster_tmp: raster template to match its grid
##
## Output - layered raster brick/stack. different layers for different time periods
###################################################################################

plot_5q0 <- function(data = df,
                     smooth     = 1,
                     age_col    = "age",
                     yr_col     = "year",
                     died_col   = "died",
                     N_col      = "N",
                     weight_col = "weight",
                     long_col   = "longitude",
                     lat_col    = "latitude",
                     dt_col     = "data_type", ## data type col
                     save_dir   = "~",
                     dim        = 10 ## height/width of pdf images
                     ){

  ## load some libraries!
  library(rgeos)
  library(raster)
  library(shapefiles)
  library(fields)
  library(akima)

  ## get the dataframe we need set up
  d <- as.data.frame(data)    ## for indexing carefulness
  d[[died_col]]  <- d[[died_col]] * d[[weight_col]]
  d[[N_col]]     <- d[[N_col]]    * d[[weight_col]]
  d <- d[, c(long_col, lat_col, yr_col, age_col, died_col, N_col, dt_col)]
  yrs <- sort(unique(d[[yr_col]]))

  if(FALSE){ ## for interpolating onto raster locations vvvvv
    ## load the population raster - we'll match the cells it uses
    ## must be on cluster!
    ## also, pull out info about the raster we'll use to predict and make the raster object
    pop <- brick('<<<< FILEPATH REDACTED >>>>>/pop_stack.tif')
    ## Convert raster to SpatialPointsDataFrame
    r.pts <- rasterToPoints(pop, spatial=TRUE)
    proj4string(r.pts)

    ## reproject sp object
    geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    r.pts <- spTransform(r.pts, CRS(geo.prj))
    proj4string(r.pts)

    ## Assign coordinates to @data slot, display first 6 rows of data.frame
    r.pts@data <- data.frame(r.pts@data, long=coordinates(r.pts)[,1],
                           lat=coordinates(r.pts)[,2])
    head(r.pts@data)
    ## make a matrix of locations where we want to project our smoothed estimate
    proj.ll <- r.pts@data[, 5:6]

  }

  ## loop through the years and make a raster for each year
  fq0_list_cbh <- list(NULL)
  fq0_list_sbh <- list(NULL)
  fq0_list <- list(NULL)


#  for(dt in c("cbh", "sbh")){
#    print(paste0("~~~ ON DATA TYPE: ", dt))
    for(y in yrs){

      print(paste0("On year: ", y))

      dy <- d[which(d[[yr_col]] == y), ]#& d[[dt_col]] == dt), ] ## get the subset of the data for the year and data type

      ## aggregate died_col by lat-long location and combine aggregates into one frame
      f.dyd <- formula(paste0(died_col, "~", long_col, "+", lat_col, "+", age_col))
      dya <- aggregate(f.dyd, dy, sum)

      f.dyn <- formula(paste0(N_col, "~", long_col, "+", lat_col, "+", age_col))
      dyn <- aggregate(f.dyn, dy, sum)
      dya <- merge(dya, dyn)

      ## merge on data_type
      dyu <- unique(dy[, c(long_col, lat_col, dt_col)])
      dya <- merge(dya, dyu)

      ## now, for each unique lat-long, calculate 5q0
      ll   <- unique(cbind(dya[[long_col]], dya[[lat_col]]))
      colnames(ll) <- c(long_col, lat_col)
      fq0    <- numeric(dim(ll)[1])
      for(i in 1:length(fq0)){

        ## get the part of the dataset at this location
        llr  <- which(dya[[long_col]] == ll[i, 1] & dya[[lat_col]] == ll[i, 2])
        dll  <- dya[llr, ]
        ages <- dll[[age_col]]

        ## get all the probabilities at this location
        if(1 %in% ages){
          r1 <- which(ages == 1)
          p1 <- sum(dll[[died_col]][r1]) / sum(dll[[N_col]][r1])
        }else{
          p1 <-  0
        }
        if(2 %in% ages){
          r2 <- which(ages == 2)
          p2 <- sum(dll[[died_col]][r2]) / sum(dll[[N_col]][r2])
        }else{
          p2 <-  0
        }
        if(3 %in% ages){
          r3 <- which(ages == 3)
          p3 <- sum(dll[[died_col]][r3]) / sum(dll[[N_col]][r3])
        }else{
          p3 <-  0
        }
        if(4 %in% ages){
          r4 <- which(ages == 4)
          p4 <- sum(dll[[died_col]][r4]) / sum(dll[[N_col]][r4])
        }else{
          p4 <-  0
        }

        ## calculate 5q0 at this location
        fq0[i] <- 1 - (1 - p1) ^ 1 *(1 - p2) ^ 11 * (1 - p3) ^ 24 * (1 - p4) ^ 24
    
      }
    
      fq0_list[[which(yrs %in% y)]] <- list('loc' = ll, 'fq0' = fq0)

    }
#  }

  ## add year names to list
  names(fq0_list_sbh) <- names(fq0_list_cbh) <- yrs
  names(fq0_list) <- yrs

  ## save pdfs for each year
  library(fields)
  ## load africa country shapefile
  as <- shapefile("<<<< FILEPATH REDACTED >>>>>africa_ad0.shp")

  ## set the colors
  breaks <- c(0, 25,
              26:50,
              51:200,
              1000)
  col.f1 <- colorRampPalette(c("#e58bba", "#f2e8b5"))
  col.f2 <- colorRampPalette(c("#f2e8b5", "#ed152e"))
  col   <- c("#74039E",
             col.f1(25),
             col.f2(150),
             "#ED152E")


  for(i in 1:4){
    pdf(paste0(save_dir, "/raw_fq0_quilt_all", names(fq0_list)[i], ".pdf"), height = dim, width = dim)
    par(mfrow = c(1, 1),
        mar   = c(5, 4, 4, 4))
    plot(as,  main = "",)
    quilt.plot(fq0_list[[i]][['loc']][, 1],
               fq0_list[[i]][['loc']][, 2],
               fq0_list[[i]][['fq0']] * 1000,
               main = "CBH",
               add = T,
               breaks = breaks, col = col,
               FUN = median)
    plot(as, add = T)

    dev.off()
  }



  if(FALSE){

    ## akima::interp() attempt

    save(list = ls(), file = "~/plot.RData")
    load("~/plot.RData")

    rows <- base:::sample(size = 2e4 ,1:length(fq0))
    ##  rows <- 1:length(fq0)
    library(akima)
    system.time(fit  <- interp(ll[rows, 1], ll[rows, 2], fq0[rows], ))

    pdf("~/interp.pdf", width=10, height=10)
    image (fit)
    contour(fit, add=TRUE)
    points (ll[rows, ], pch = 3)
    dev.off()
  }

  if(FALSE){
    ## Thin plate spline attempt

    ## now we make the thin-plate smoothed version
    dircos<- function(x1){
      coslat1 <- cos((x1[, 2] * pi)/180)
      sinlat1 <- sin((x1[, 2] * pi)/180)
      coslon1 <- cos((x1[, 1] * pi)/180)
      sinlon1 <- sin((x1[, 1] * pi)/180)
      cbind(coslon1*coslat1, sinlon1*coslat1, sinlat1)
    }

    tps.name  <- paste0("tps.fit.", y)
    system.time( ## this is slow
      tps.fit <- Tps(dircos(ll), fq0)
    )
    assign(tps.name, tps.fit)

    ## project where we want to get estimates
    ## preds <- predict(tps.fit, x=as.matrix(proj.ll))

    ## now we just need to store these back in the raster correctly (or make a new raster)
  }

  print("Done making images")
  message("Done making images")

}

#######################################################################################
## this one does a coarser summed estimate
## it sums exp in a box, sums deaths in a box and maps 1-(1-sum(death)/sum(exp))^60
## WARNING: still in development though all the code needed to produce maps is located internally
plot_5q0_summed <- function(data = df,
                            smooth     = 1,
                            age_col    = "age",
                            yr_col     = "year",
                            died_col   = "died",
                            N_col      = "N",
                            weight_col = "weight",
                            long_col   = "longitude",
                            lat_col    = "latitude",
                            dt_col     = "data_type", ## data type col
                            save_dir   = "~",
                            dim        = 10 ## height/width of pdf images
                     ){

  ## load some libraries!
  library(rgeos)
  library(raster)
  library(shapefiles)
  library(fields)
  library(akima)

  ## get the dataframe we need set up
  d <- as.data.frame(data)    ## for indexing carefulness
  d[[died_col]]  <- d[[died_col]] * d[[weight_col]]
  d[[N_col]]     <- d[[N_col]]    * d[[weight_col]]
  d <- d[, c(long_col, lat_col, yr_col, age_col, died_col, N_col, dt_col)]
  yrs <- sort(unique(d[[yr_col]]))

  if(FALSE){ ## for interpolating onto raster locations vvvvv
    ## load the population raster - we'll match the cells it uses
    ## must be on cluster!
    ## also, pull out info about the raster we'll use to predict and make the raster object
    pop <- brick('<<<<< FILEPATH REDACTED >>>>>/pop_stack.tif')
    ## Convert raster to SpatialPointsDataFrame
    r.pts <- rasterToPoints(pop, spatial=TRUE)
    proj4string(r.pts)

    ## reproject sp object
    geo.prj <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
    r.pts <- spTransform(r.pts, CRS(geo.prj))
    proj4string(r.pts)

    ## Assign coordinates to @data slot, display first 6 rows of data.frame
    r.pts@data <- data.frame(r.pts@data, long=coordinates(r.pts)[,1],
                           lat=coordinates(r.pts)[,2])
    head(r.pts@data)
    ## make a matrix of locations where we want to project our smoothed estimate
    proj.ll <- r.pts@data[, 5:6]

  }

  ## set the colors
  breaks <- c(0, 25,
              26:50,
              51:200,
              202)
  col.f1 <- colorRampPalette(c("#e58bba", "#f2e8b5"))
  col.f2 <- colorRampPalette(c("#f2e8b5", "#ed152e"))
  col   <- c("#74039E",
             col.f1(25),
             col.f2(150),
             "#ED152E")

  library(fields)
  ## load africa country shapefile
  as <- shapefile("<<<< FILEPATH REDACTED >>>>>/africa_ad0.shp")

  for(y in yrs){

    print(paste0("On year: ", y))

    pdf(paste0(save_dir, "/raw_fq0_new_all", y, ".pdf"), height = dim, width = dim)
    par(mfrow = c(1, 1),
        mar = c(5, 4, 4, 5))

      dy <- d[which(d[[yr_col]] == y), ]

      ## get the prob of dying in the age bin from the sums in the square
      ## first, make a grid that will be constant for y and dt
      gr.y.dt <- discretize.image(x = dy[, 1:2], m = 64, n = 64, grid = NULL, boundary.grid = FALSE)$grid

      for(a in 1:4){
        ## subset the data
        dya <- subset(dy, get(age_col) == a)

        ## now we use the hacked quilt.plot to sum on grid squares
        sum.deaths <- as.image(Z = dya$died, x = dya[, 1:2], nx = 64, ny = 64, na.rm = TRUE,
                               grid = gr.y.dt, FUN = sum)
        sum.exps   <- as.image(Z = dya$N, x = dya[, 1:2], nx = 64, ny = 64, na.rm = TRUE,
                               grid = gr.y.dt, FUN = sum)

        sum.ratio  <- sum.deaths$z / sum.exps$z
        not.na.z <- which(!is.na(sum.ratio))

        if(a == 1){
          dyp  <- cbind(not.na.z,
                        sum.ratio[not.na.z])
          colnames(dyp) <- c("n.na.z", paste0("s.r", a))
        }else{
          dy.t <-  cbind(not.na.z,
                         sum.ratio[not.na.z])
          colnames(dy.t) <- c("n.na.z", paste0("s.r", a))
          dyp <- merge(dyp, dy.t, all = T, ny = "n.na.z")
        }
      }

      ## fill in NAs with zeros
      dyp <- as.matrix(dyp)
      dyp[which(is.na(dyp))] <- 0

      ## make 5q0 estimates
      fq0 <- 1 - (1 - dyp[, 2]) ^ 1 *  (1 - dyp[, 3]) ^ 11 *  (1 - dyp[, 4]) ^ 24 *  (1 - dyp[, 5]) ^ 24

      ## mask values higher than 200 to 201 for plotting purposes
      ## the color scale stops at 200 anyways
      fq0[which(fq0 > 200)] <- 201

      ## set up quilt.plot
      plot.im <-  as.image(Z = dy$died, x = dy[, 1:2], nx = 64, ny = 64, na.rm = TRUE,
                           grid = gr.y.dt, FUN = sum)

      ## fill z with our numbers
      z <- plot.im$z
      z[dyp[, 1]] <- fq0
      plot.im$z <- z * 1000 ## convert to counts

      ## now we can plot it
      plot(as,  main = "All")
      image.plot(plot.im, nlevel = length(col), col = col, add = TRUE, breaks = breaks)
    plot(as, add = TRUE)

      dev.off()

    }

#    dev.off()
#  }

  ## add year names to list
  names(fq0_list_sbh) <- names(fq0_list_cbh) <- yrs

  ## save pdfs for each year
  library(fields)
  ## load africa country shapefile
  as <- shapefile("<<<< FILEPATH REDACTED >>>>>/africa_ad0.shp")

  ## set the colors
  breaks <- c(0, 25,
              26:50,
              51:200,
              1000)
  col.f1 <- colorRampPalette(c("#e58bba", "#f2e8b5"))
  col.f2 <- colorRampPalette(c("#f2e8b5", "#ed152e"))
  col   <- c("#74039E",
             col.f1(25),
             col.f2(150),
             "#ED152E")


  for(i in 1:4){
    pdf(paste0(save_dir, "/raw_fq0_quilt_cbh_sbh", names(fq0_list_cbh)[i], ".pdf"), height = dim, width = dim * 2)
    par(mfrow = c(1, 2),
        mar   = c(5, 4, 4, 10))
    plot(as,  main = "CBH",)
    quilt.plot(fq0_list_cbh[[i]][['loc']][, 1],
               fq0_list_cbh[[i]][['loc']][, 2],
               fq0_list_cbh[[i]][['fq0']] * 1000,
               main = "CBH",
               add = T,
               breaks = breaks, col = col,
               FUN = median)
    plot(as, add = T)

    plot(as,  main = "SBH",)
    quilt.plot(fq0_list_sbh[[i]][['loc']][, 1],
               fq0_list_sbh[[i]][['loc']][, 2],
               fq0_list_sbh[[i]][['fq0']] * 1000,
               main = "SBH",
               add = T,
               breaks = breaks, col = col,
               FUN = median)
    plot(as, add = T)

    dev.off()
  }



  if(FALSE){

    ## akima::interp() attempt

    save(list = ls(), file = "~/plot.RData")
    load("~/plot.RData")

    rows <- base:::sample(size = 2e4 ,1:length(fq0))
    ##  rows <- 1:length(fq0)
    library(akima)
    system.time(fit  <- interp(ll[rows, 1], ll[rows, 2], fq0[rows], ))

    pdf("~/interp.pdf", width=10, height=10)
    image (fit)
    contour(fit, add=TRUE)
    points (ll[rows, ], pch = 3)
    dev.off()
  }

  if(FALSE){
    ## Thin plate spline attempt

    ## now we make the thin-plate smoothed version
    dircos<- function(x1){
      coslat1 <- cos((x1[, 2] * pi)/180)
      sinlat1 <- sin((x1[, 2] * pi)/180)
      coslon1 <- cos((x1[, 1] * pi)/180)
      sinlon1 <- sin((x1[, 1] * pi)/180)
      cbind(coslon1*coslat1, sinlon1*coslat1, sinlat1)
    }

    tps.name  <- paste0("tps.fit.", y)
    system.time( ## this is slow
      tps.fit <- Tps(dircos(ll), fq0)
    )
    assign(tps.name, tps.fit)

    ## project where we want to get estimates
    ## preds <- predict(tps.fit, x=as.matrix(proj.ll))

    ## now we just need to store these back in the raster correctly (or make a new raster)
  }

  print("Done making images")
  message("Done making images")

}

## ################################
## is_oos_preds
##
## this function takes the data that went into our MBG framework and
## returns an array containing: location, observed value, in sample
## predictions, and out of sample predictions (when applicable), country, year, N
##
## INPUT:
##
## OUTPUT:
##
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
  output.dir <- sprintf("<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/",
                        indic_group, indic, rd)

  ## Load data
  datdir <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/',indic_group,indic,rd)

  ## holdout data
  if(!is.null(holdoutlist)){
    d <- data.frame(holdoutlist[[sprintf('region__%s', reg)]])
  }else{
    d <- df
  }

  ## load the simple raster for this region
  load(paste0('<<<< FILEPATH REDACTED >>>>>/simple_raster',reg,'.RData'))

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
      load(sprintf(paste0('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/', cell_draws_filename),
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
  load(sprintf(paste0('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/', cell_draws_filename),
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
                             get.oos = FALSE
                             ){

  ## ##################################################################
  ## load in data, get domain templates and add on necessary columns ##
  ## ##################################################################

  mod.dir <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/', ind_gp, ind, rd)

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
    Regions <- get_output_regions(in_dir = paste0('<<<< FILEPATH REDACTED >>>>>/', ind_gp, '/', ind, '/output/', rd))
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
    df <- as.data.table(do.call(rbind, readRDS(paste0(mod.dir, 'stratum_qt.RDs')))) ## holdoutlist
  }


  ## make year to period key
  yr.per <- cbind(yrs, 1:nperiod)
  colnames(yr.per) <- c('yr', 'per')

  ## load in admin1 and 2 relevant shapefiles
  message('loading in shapefiles to determine ad1 and ad2 membership')
  ## Load admin2 raster
  admin_level <- 2
  shapes <- readRDS("<<<< FILEPATH REDACTED >>>>>/g_2015_2014_2_modified.rds")
  
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
    
    ## get rows in region
    df.r  <- which(df$region == rr)
    loc.r <- cbind(df$longitude, df$latitude)[df.r, ]

    ## make a matrix to fill with draws
    draw.mat <- matrix(ncol = n.pred, nrow = length(df.r))

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
      yr.r <- which(df[df.r, ]$original_year == yrs[yy])
      idx.yy <- raster::extract(idx.rast[[yy]], matrix(loc.r[yr.r, ], ncol = 2))

      ## enter draws into draw.mat
      draw.mat[yr.r, ] <- cell_pred[idx.yy, 1:n.pred]
    }

    if(which(all.regions == rr) == 1){
      domain.draws <- matrix(ncol = n.pred, nrow = nrow(df))
      colnames(domain.draws) <- paste0("draw", 1:n.pred)
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
    
    n.pred <- ncol(cell_pred)
    if(rd == '2017_11_06_23_09_00' |
       rd == '2017_11_06_23_08_52' |
       rd == '2017_11_06_23_08_39'){n.pred <- 250} ## these have 1000 draws, the rest have 250


    domain.draws <- matrix(ncol = n.pred, nrow = nrow(df.oos))
    ## we should already have a cell_preds loaded since IS happens first
    colnames(domain.draws) <- paste0("draw", 1:n.pred)

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
      draw.mat <- matrix(ncol = n.pred, nrow = length(df.oos.r))

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
          yr.rh <- which(df.oos[df.oos.rh.idx, ]$original_year == yrs[yy]) ## index within region and holdout for yy
          idx.yy <- raster::extract(idx.rast[[yy]], matrix(loc.rh[yr.rh, ], ncol = 2))

          ## enter draws into draw.mat
          draw.mat[df.oos.rh[yr.rh], ] <- cell_pred[idx.yy, 1:n.pred] ## draws in region rr,
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

  regions <- get_output_regions(in_dir = paste0('<<<< FILEPATH REDACTED >>>>>',
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

    in_dir  <- paste0('<<<<< FILEPATH REDACTED >>>>>', ind.gp, '/', ind, '/output/', rd)
    default_rf_path <- paste0(in_dir, '/', ind, '_rf.csv')
    all_rfs <- fread(default_rf_path)
    gaul_list = get_gaul_codes(rr)
    rfs <- all_rfs[name %in% gaul_list, ]
    loc_names <- setDT(read.csv("<<<< FILEPATH REDACTED >>>>>/gaul_to_loc_id.csv",stringsAsFactors = F))
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
                                 sharedir = sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s',indicator_group,indicator)){
  ## load models
  require(INLA)
  message(sprintf('Pulling together results for %s models',rd))

  tlist=list()

  for(rr in regions){
    message(sprintf('::on region %s',rr))
    reg  <-  rr

    message("::::loading in pre-INLA objects to get spde")
    pathaddin  <-  paste0('_bin',age,'_',rr,'_',holdout)
    load(paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/',
                ind, '/model_image_history/', rd, pathaddin,
                '.RData'))

    modnames = c('gam','gbm','ridge','enet','lasso')

    full_raster_list <- cov_list

    ## hotfix!! inv.logit ## TODO ## TODO don't need this if either
    ##   1) resolve logit fits
    ##   2) have new runs where I didn't logit stackers
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
                            save.dir) {

  if(ind == 'wasting_mod_b') nice.name = "Wasting"
  if(ind == 'stunting_mod_b') nice.name = "Stunting"
  if(ind == 'underweight_mod_b') nice.name = "Underweight"

  ## setup the dataframe
  subset_shape <- subset_shape[subset_shape$GAUL_CODE %in% gaul_list, ]

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
      guides(color=guide_colorbar(title="Bias",
                                  label=TRUE,
                                  ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      theme(panel.margin = unit(0,"lines"),
            plot.margin = unit(c(0,0,0,0),"lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~original_year) +
      ggtitle(paste0(nice.name, 'bias'))

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
      guides(color=guide_colorbar(title="Bias",
                                  label=TRUE,
                                  ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      theme(panel.margin = unit(0,"lines"),
            plot.margin = unit(c(0,0,0,0),"lines"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~original_year) +
      ggtitle(paste0(nice.name, ' bias'))

      ggsave(filename = sprintf('%s%s_abs_error_plot_OOS.png', save.dir, ind),
             plot = gg.oos, width = 12, height = 12, units = 'in')

    }else{
      gg.oos <- NULL
    }
  }

  return(list(is = gg.is,
              oos = gg.oos))
}


####################################################################################################################
## INPUT: # data file matching 
## OUTPUT SUMMARIZED PREDICTIVE VALIDAITY METRICS
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
                        out.dir) {

  d <- data.table(d)

  # Acct for polygon weights
  d$indic_orig     <-  d[[indicator]]
  d$exposure       <-  d$N*d$weight
  d[[indicator]]   <-  d[[indicator]]*d$weight

  # get binomial predictions
  message('Simulating predictive draws')
  x <- matrix(NA,ncol=draws,nrow=nrow(d))
  pb <- txtProgressBar(min = 1, max = draws, initial = 1)
  if(family == 'binomial'){
    for(i in 1:draws){
      x[,i] <- rbinom(nrow(d),
                      size=round(d$N,0),
                      prob=d[[paste0('draw',i)]])
      setTxtProgressBar(pb,i)
    }
  }
  if(family == 'gaussian'){
    for(i in 1:draws){
      x[,i] <- rnorm(n=nrow(d),
                     sd=(d[[paste0('tau_',i)]])^(-1/2),
                     mean=d[[paste0('draw',i)]])
      setTxtProgressBar(pb,i)
    }
  }

  # show warning for those that got dropped
  message(paste0('NOTE: missing predictions for ', sum(is.na(x[,1])),' of ',nrow(x),' (',round(sum(is.na(x[,1]))/nrow(x)*100,0),'%) data points.'))

  # get coverages
  message('Getting Coverage')
  for(c in coverage_probs){
    message(paste0('For ',c,'% coverage.'))
    coverage <- c/100
    li       <- apply(x,1,quantile,p=(1-coverage)/2,na.rm=T)
    ui       <- apply(x,1,quantile,p=coverage+(1-coverage)/2,na.rm=T)
    d[,paste0('clusters_covered_',c)] = d[['indic_orig']]>=li & d[['indic_orig']] <= ui
  }

  d$oos <- d$fold != 0

  # collapse data into areas of spatial aggregation
  message('Collapse into areas of spatial aggregation')
  d$total_clusters <- 1
  res <- d[,c('oos','year',indicator,aggregate_on,'exposure',paste0('clusters_covered_',coverage_probs),'total_clusters',paste0('draw',1:draws)),with=FALSE]
  for(i in 1:draws) res[[paste0('draw',i)]] = res[[paste0('draw',i)]]*res$exposure # expected number for binomial
  for(c in coverage_probs) res[[paste0('clusters_covered_',c)]] = res[[paste0('clusters_covered_',c)]]*res$exposure
  f     <- as.formula(paste('.~oos+year',aggregate_on,sep='+'))
  res   <- data.table(aggregate(f,data=res,FUN=sum))

  # transform back to probability/coverage
  res$p   <- res[[indicator]]/res$exposure ## true data p
  for(c in coverage_probs) res[[paste0('clusters_covered_',c)]] = res[[paste0('clusters_covered_',c)]] / res$exposure

  for(i in 1:draws){ ## phat is our estimated p from mbg modeling
    res[[paste0('phat_',i)]]=res[[paste0('draw',i)]]/res$exposure
    res[[paste0('draw',i)]]=NULL
  }
  res[[indicator]]=NULL

  # estimates at aggregated level
  res$mean_phat   <- apply(res[,paste0('phat_',1:draws),with=FALSE],1,mean)
  res[, error := p - mean_phat]
  res$var_phat <- apply(res[,paste0('phat_',1:draws),with=FALSE],1,var)
  res$cov <- res$mean_phat / res$var_phat

  # choose to weight on SS or not.
  if(weighted==TRUE){
      res$weight = res$exposure
  } else {
      res$weight = 1
  }

  require(boot)
  # quick define a weighted RMSE
  weighted.rmse <- function(error, w) sqrt( sum(  (w/sum(w)) * ((error)^2) ) )

  message('Outputting predictive validity metrics for your models')
  res2<- data.table(cbind(res[,.(me   = weighted.mean(error,w=weight)),      by=result_agg_over],
                          mean_phat   = res[,.(mean_phat   = weighted.mean(mean_phat,w=weight)),  by=result_agg_over]$mean_phat,
                          mean_p      = res[,.(mean_p      = weighted.mean(p,   w=weight)), by=result_agg_over]$mean_p,
                          rmse        = res[,.(rmse        = weighted.rmse(error,w=weight)),      by=result_agg_over]$rmse,
 #                         cov         = res[,.(mean_cov    = weighted.mean(cov, w=weight)),  by=result_agg_over]$mean_cov, 
                          median_SS   = res[,.(avg_exp     = median(exposure)),                   by=result_agg_over]$avg_exp,
                          cor         = res[,.(cor         = corr(cbind(p,mean_phat), w=weight )),by=result_agg_over]$cor,
                          coverage_95 = res[,.(coverage_95 = weighted.mean(clusters_covered_95, w=weight)) ,  by=result_agg_over]$coverage_95))

  
  setorder(res2, oos, year)
  colnames(res2) <- c('Year', 'OOS', 'Mean Err.',
                      'Mean Pred.', 'Mean Obs.', 'RMSE',
                      #'CoV', 
                      'Median SS', 'Corr.', '95% Cov.')


  if(plot==TRUE){
    require(ggplot2)
    res$SS = res$weight

    for(oosindic in unique(res$oos)){
      message(paste0('Saving plot here: ', paste0(out.dir, indicator, '_validation_plot_', aggregate_on,'OOS_', oosindic, '.png')))
      png(paste0(out.dir, indicator, '_validation_plot_', aggregate_on,'OOS_', oosindic, '.png'), width = 12, height = 12, units = 'in', res = 350)

      tmp = subset(res,oos==oosindic)
      tmp$upper   <- apply(tmp[,paste0('phat_',1:draws),with=FALSE],1,quantile,p=0.01*(plot_ci_level+(100-plot_ci_level)/2),rm.na=TRUE)
      tmp$lower   <- apply(tmp[,paste0('phat_',1:draws),with=FALSE],1,quantile,p=0.01*((100-plot_ci_level)/2),rm.na=TRUE)
      tmp         <- na.omit(tmp[,c('p','mean_phat','SS','upper','lower','year'),with=FALSE])
      gg=ggplot(tmp,aes(x=p,y=mean_phat,size=SS))+
      geom_abline(intercept=0,slope=1,colour='red')+
      geom_point(colour='grey')+
      theme_bw()+
      xlab('Data Estimate') +
      ylab('Mean Prediction')+
      facet_wrap(~year,ncol=2)+
      theme(strip.background = element_rect(fill="white"))+
      ggtitle(paste0('OOS: ',oosindic))
      if(plot_ci) gg <- gg +  geom_errorbar(aes(ymin=lower, ymax=upper), colour="black", width=0, size=.3, alpha = 0.2)
      gg <- gg + geom_abline(intercept=0,slope=1,colour='red')
      plot(gg)
      dev.off()
    }

  }


  return(res2)
}
