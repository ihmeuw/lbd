## this file contains functions to make different types of holdouts for spatio-temporal data

#######################
#######################
### HOLDOUT OPTIONS ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######################
#######################


####################
## TOTALLY RANDOM ##
####################

## 1) As it sounds, totally random points across space and time


##############
## IN SPACE ##
##############

## 1) random in space (this is just totally random with one timepoint)
## 2) small 'tesselations' aggregated up until certain population is reached
##    a) with quadtree
##    b) NOT DONE: with weighted k-means (weights inversely proportional to sample size)
## 3) by admin2
## 4) countries
## 5) larger regions (e.g. East Africa, West Africa, ...)


##############
## IN  TIME ##
##############

## 1) random in time
## 2) proportional to data amount in the period
## 3) full years, randomly selected across the duration of data
## 4) build up different datasets as if we were moving chronologically in time
##    i.e. first set is only yr 1, second set is yrs 1 & 2, third set is yrs 1, 2, 3 ...


##############
## IN   S-T ##
##############

## all the combos of the space time sets!


#####################
#####################
### CODING STARTS ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####################
#####################

## ## load packages
## library(raster)
## library(seegSDM)
## library(seegMBG)
## library(rgdal)
## library(rgeos)
## library(foreach)
## library(doParallel)
## library(data.table)

####################
## TOTALLY RANDOM ##
####################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## makes folds completely at random across all data
##
## INPUTS:
## data: cleaned data.table/frame to be split into folds
## n_folds: number of folds
## strat_cols: vector of column string names to
##    stratify over when making folds. if NULL, holdout
##    sets are made across the entire data structure
##
## OUTPUTS: 2 item list
## 1st item: 1 by nrow(data) vector containing integers identifying folds
## 2nd item: matrix of stratification combinations used to make holdouts
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
rand_folds <- function(data,
                       n_folds=5,
                       strat_cols=NULL,
                       ...){

  ## make fold vectors stratifying by specified columns
  if(!is.null(strat_cols)){

    ## get different stratum
    all_strat <- get_strat_combos(data=data, strat_cols=strat_cols)

    ## make a vector to identify folds and assign completely at random by strata
    fold_vec <- rep(NA, nrow(data))

    for(strat in 1:nrow(all_strat)){

      strata <- as.data.frame(all_strat[strat, ])
      colnames(strata) <- colnames(all_strat)

      ## get rows in current strata
      strat_rows <- get_strat_rows(data=data,
                                   strata=strata)

      ## assign fold numbers uniformly (with rounding)
      fold_s  <- cut(seq(1, sum(strat_rows)),
                     breaks = n_folds, labels = 1:n_folds)
      fold_s <- as.numeric(as.character(fold_s))

      ## randomize the numbers within the strata
      fold_vec[which(strat_rows==1)]  <- sample(fold_s)
    }

    ## check to make sure we got all rows
    if(sum(is.na(fold_vec) > 0)){
      message("Warning! Check why some data rows didn't get assigned to folds")
    }

  }else{ ## make folds across all data
    fold_vec <- sample(cut(seq(1, nrow(data)),
                           breaks = n_folds, labels = 1:n_folds))
  }


  ## and return a list containing
  ## 1) vector of folds
  ## 2) matrix containing all stratum
  if(is.null(strat_cols)){
    return(list(folds   = fold_vec,
                stratum = NA))
  }else{
    return(list(folds   = fold_vec,
                stratum = all_strat))
  }
}

##############
## IN SPACE ##
##############

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## (1) random in space.
## this is the same as totally random with one time point
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## (2) small aggregates in space - QUADTREE
##
## INPUTS:
##
## OUTPUTS:
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
quadtree_folds <- function(xy, ## 2 col matrix of xy locs
                           ss, ## vector of sample size at each loc - if all 1s, it aggregates by # of points
                           mb, ## minimum allowed pts in a quadtree region
                           ts, ## target sample size in each region
                           n_folds, ## number of folds,
                           plot_fn  = NULL,  ## if true plots data and quad tree must be .png
                           plot_shp = NULL,  ## to add shapefile outlines to plot
                           save_qt  = TRUE,  ## if desired, can save quadtree regions to shapefiles
                           ...,
                           t_fold = 1, ## if multiple t_folds. to save shapefiles
                           stratum = 1 ## for multiple stratum to save shapefiles
#                           pathaddin = ""   ## file path addin to specifiy specifis of model run
#                           run_date = run_date
                           ){

  ## this function segregates data into spatial regions. it splits
  ## data down until the target sample size is reached or until a
  ## minimum allowable number of points are in the bin

  ## it then breaks the data into folds by ensuring that all data in
  ## any one spatial region stays together and such that the sample
  ## size sums in each fold are relatively similar

  ## points at the same place pose a problem because you can't split them
  ## so, first we make a subset of the data using only unique xy combos
  ## quadtree is run on the unique list and then we buid back up to
  ## the full dataset and make folds

  ## make data.table
  library(data.table)
  dt_xy <- data.table(cbind(xy, ss))

  ## round to 5 decimal points to ensure matching
  cols <- names(dt_xy)[1:2]
  dt_xy[, (cols) := round(.SD, 5), .SDcols=cols]

  ## griyo by lat and long combos
  dt_xy <- dt_xy[, loc_id := .GRP, by = c('long','lat')]

  ## unique xy locations by group and sum ss in groups
  un_xy <- dt_xy[, list(long=mean(long), lat=mean(lat), ss=sum(ss)), by=c('loc_id')]

  ## to allow for deeper recursions in quadtree
  options(expressions=500000)

  ## get quad-tree regions on unique locations and ss sum at those locations
  system.time(qt <- quadtree_ct(xy=as.matrix(un_xy[, .(long, lat)]),
                                ss=as.numeric(un_xy[, ss]),
                                target_ss=ts,
                                min_in_bin=1,
                                rand = T))

  ## in order to match these back to all the points, we collect locations and ids
  ## this can take a little while...
  ids <- id(qt)
  ids <- na.omit(ids) ## qt alg adds a few NA rows sometimes...

  ## now we match ids to our original dataset
  dt_xy[, row:=1:nrow(dt_xy) ]
  setkey(dt_xy, long, lat)
  dt_xy[, qt_id:= -1.0]
  system.time(
    for(i in 1:nrow(ids)){
      ## this next line is to fix weird dropping 0 issues resulting in non-matches
      loc <- data.table(long = round(ids[i, 2], 6), lat = round(ids[i, 3], 6))
      dt_xy[.(loc), qt_id:=ids[i, 1] ]
    }
  )

  ## now we stratify - selecting by qt_id and ensuring sum of ss in
  ## folds is close to even across folds
  if(length(unique(dt_xy[, qt_id])) < n_folds){
    message("Warning: Fewer quadtree sections than n_folds!! Things may break. Either increase data amount or decrease ts")
  }

  cts_in_qt <- dt_xy[, sum(ss), by = qt_id]
  fold_vec <- make_folds_by_poly(cts_in_polys = as.matrix(cts_in_qt),
                                 pt_poly_map  = as.numeric(dt_xy[,qt_id]),
                                 n_folds      = n_folds)

  ## now we 'unsort' dt_xy to make it match the orde of the original data
  dt_xy[, fold_vec := fold_vec]
  setkey(dt_xy, row)
  fold_vec <- dt_xy[, fold_vec]
  ho_id    <- dt_xy[, qt_id]

  ## plot if desired
  if(!is.null(plot_fn)){
    xylim <- cbind(x=c(min(xy[,1]), max(xy[,1])), y=c(min(xy[,2]), max(xy[,2])))
    png(paste0(plot_fn), width=4000, height=4000)
    title <- paste0("Quadtree with threshold of ", ts)
    if(!is.null(plot_shp)){
      plot(plot_shp, , xlab="x", ylab="y", main=title)
    }else{
      plot(xylim, type="n", xlab="x", ylab="y", main=title)
    }
    lines(qt, xylim, col="Gray")
    cols=hsv(fold_vec/n_folds, 1, 1)
    library(scales)
    points(dt_xy[,.(long, lat)], col=alpha(cols, alpha=0.5), pch=16, cex=2)
    dev.off()
  }

  ## save quadtree rectangles to shapefile if desired
  if(save_qt){
    output_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date)

    xylim <- cbind(x=c(min(un_xy[,2]), max(un_xy[,2])), y=c(min(un_xy[,3]), max(un_xy[,3])))
    polys <- cell(qt, xylim)
    polys_attr <- data.frame(id=unique(polys$id))
    library(shapefiles)
    polys_shapefile <- convert.to.shapefile(polys, polys_attr, "id", 5)
    write.shapefile(polys_shapefile, paste0(output_dir, '/', 'spat_holdout_stratum_', stratum, '_t_fold', t_fold), arcgis=TRUE)
  }

  return(cbind(fold_vec, ho_id))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## (2) small aggregates in space - WEIGHTED K-MEANS
## not yet (or ever?) implemented. quadtree looks pretty good
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## (3) admin2 in space
##
## INPUTS:
##
## OUTPUTS:
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
ad2_folds <- function(admin_shps,
                      shape_ident="gaul_code",
                      admin_raster,
                      ss=1,
                      xy,
                      n_folds,
                      mask_shape,
                      mask_raster,
                      ...){
  ## this function takes in admin2 (or any mutually exclusive and
  ## collectively exhaustive) shapefiles covering the domain of your
  ## data and splits your data into folds of approximately equal
  ## sample size using admin2 units to split the data

  ## admin_shps: file location of all pertinent shapefiles to use when folding (.grd)
  ## shape_ident: string identifying data col in shapefile used to uniquely identify polygons
  ## admin_raster: file location of associated raster for admin_shps (.shp)
  ## data: complete dataset that you want to fold
  ## strat_cols: vector of column string names to
  ##    stratify over when making folds. if NULL, holdout
  ##    sets are made across the entire data structure
  ## ss: vector of sample sizes for each row. if <1>, assumes all have ss=1
  ## xy: matrix of xy coordinates
  ## n_folds: number of folds to make
  ## mask_shape: shapefile file location for boundary of area to be folded
  ## mask_raster: raster to project mask_shape onto

  library(raster)

  ## make a mask for ther region we care about
  mask <- rasterize(shapefile(mask_shape), raster(mask_raster))*0

  ## get raster cells in mask
  cell_idx <- cellIdx(mask)

  ## load raster and shapefile for admin units
  rast      <- raster(admin_shps)
  rast_cell <- extract(rast, cell_idx)
  shp_full  <- shapefile(admin_raster)
  shp       <- shp_full@data[c('name', 'gaul_code')]
  ## plot(shp_full, col=1:length(shp_full))

  ## match raster polys to gaul
  rast_code <- match(rast_cell, shp$gaul_code)
  rast_code <- shp$gaul_code[rast_code]

  ## get number of datapoints in shapefiles
  shp_cts <- get_sample_counts(ss          = ss,
                               xy          = xy,
                               shapes      = shp_full,
                               shape_ident = shape_ident)
  pt_shp_map  <- shp_cts$pt_poly_map
  cts_in_shps <- shp_cts$ct_mat

  ## make the folds by using the polygons (ad2 units) as the

  ## sampling unit to divide the data into the folds while keeping a
  ## similar sample size in each fold
  fold_vec <- make_folds_by_poly(cts_in_polys = cts_in_shps,
                                 pt_poly_map  = pt_shp_map,
                                 n_folds      = n_folds)
  ho_id <- pt_shp_map



  return(cbind(fold_vec,
               ho_id))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## (4) countries in space
##
## INPUTS:
##
## OUTPUTS:
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
ct_folds <- function(xy,   ## xy location matrix
                     ct,   ## country vec
                     ss=1, ## sample size vec (or 1 if equi-ss)
                     n_folds=5,
                     ...){

  if(length(unique(yr) < n_folds)){
    message("Too many folds for too few countries! Expand your horizons")
    stop()
  }

  if(length(ss)==1) ss <- rep(1, nrow(xy))

  ## first we find the sample size in each of the countries
  library(data.table)
  dt <- data.table(long = xy[,1],
                   lat  = xy[,2],
                   ss   = ss,
                   ct   = ct)

  ## get sample size totals in each country
  cts_in_ct <- dt[, sum(ss), by=ct]

  ## make the folds
  fold_vec <- make_folds_by_poly(cts_in_polys = as.data.frame(cts_in_ct),
                                 pt_poly_map  = as.character(dt[,ct]),
                                 n_folds      = n_folds)


  return(fold_vec)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## (4) random folds in space
##
## INPUTS:
##
## OUTPUTS:
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
rand_s_folds <- function(xy,   ## xy location matrix
                         ss=1, ## sample size vec (or 1 if equi-ss)
                         n_folds=5,
                         ...){

  if(length(ss) < n_folds){
    message("Too many folds for too few countries! Expand your horizons")
    stop()
  }

  if(length(ss) != nrow(xy)){
    message("length of ss and nrow(xy) must match!")
    stop()
  }

  if(length(ss)==1) ss <- rep(1, nrow(xy))

  total.ct <- sum(ss)
  max.fold.ct <- ceiling(total.ct/n_folds)

  ## make a vector to store folds
  fold.vec <- numeric(nrow(xy))

  ## randomize the order of the rows
  rand.ord <- sample(1:nrow(xy))

  ## randomly decide if the first one will include the final poly
  flip <- sample(0:1, 1)
  if(flip==0){
    include.final <- rep(0:1, ceiling(n_folds/2))
  }else{
    include.final <- rep(1:0, ceiling(n_folds/2))
  }

  fold.sums <- rep(NA, n_folds)
  start.ind <- 1
  stop.ind  <- 1
  for(fold in 1:(n_folds-1)){

    ## check threshhold
    while(sum(ss[rand.ord[start.ind:stop.ind]]) < max.fold.ct){
      stop.ind <- stop.ind + 1
    }

    ## check if final row is included
    if(include.final[fold]==1){
      stop.ind <- stop.ind - 1
    }

    ## identify the selected rows as being in this fold
    fold.vec[rand.ord[start.ind:stop.ind]] <- fold

    ## record total in fold
    total.ct[fold] <- sum(ss[rand.ord[start.ind:stop.ind]])

    ## adjust indices
    start.ind <- stop.ind + 1
    stop.ind  <- start.ind + 1
  }

  ## the last fold is everything else
  stop.ind <- nrow(xy)
  fold.vec[rand.ord[start.ind:stop.ind]] <- n_folds
  total.ct[n_folds] <- sum(ss[rand.ord[start.ind:stop.ind]])

  ## print ss sums in folds
  message("The sum in each different fold is: \n")
  for(i in 1:n_folds){
    message(total.ct[i])
  }

  return(fold.vec)
}


##############
## IN  TIME ##
##############

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 1) random in time - already done with completely random
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 2) proportional to data amount in the period
## already done. this is equivalent to stratifying by time the
## stratified holdouts can then be recombined across time to yield
## proportional to data in the time period
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
proptime_folds <- function(yr, n_folds = NULL, ss=NULL, ...){

  ## so all we have to do is return a 'fold vector' with integers for unique years

  yrf  <- as.factor(yr)
  fold_vec <- as.numeric(yrf)

  fold_list <- list(NULL)
  for(i in sort(unique(fold_vec))){
    fold_list[[i]] <- which(fold_vec == i)
  }
  return(fold_list)

}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 3) full years, randomly selected across the duration of data
##
## INPUTS:
##
## OUTPUTS:
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
yr_folds <- function(xy,   ## xy location matrix
                     yr,   ## yr vec
                     ss=1, ## sample size vec (or 1 if equi-ss)
                     n_folds=5,
                     ...){

  if(length(unique(yr) < n_folds)){
    message("Too many folds for too few years! Try again in a decade")
    stop()
  }

  if(length(ss)==1) ss <- rep(1, nrow(xy))

  ## first we find the sample size in each of the countries
  library(data.table)
  dt <- data.table(long = xy[,1],
                   lat  = xy[,2],
                   ss   = ss,
                   yr   = yr)

  ## get sample size totals in each country
  cts_in_yr <- dt[, sum(ss), by=ct]

  ## make the folds
  fold_vec <- make_folds_by_poly(cts_in_polys = as.data.frame(cts_in_yr),
                                 pt_poly_map  = as.character(dt[,ct]),
                                 n_folds      = n_folds)

  fold_list <- list(NULL)
  for(i in sort(unique(fold_vec))){
    fold_list[[i]] <- which(fold_vec == i)
  }
  return(fold_list)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## 4) build up different datasets as if we were moving chronologically in time
##    e.g. 2000. 2000 & 2005. 2000, 2005 & 2010. ...
##
## INPUTS:
## yr: 1 by #datapt vector containing the year for each datapt
##
## OUTPUTS: list with as many entries as unique years
## each item in the list contains all data rows that should be in that fold
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
chrono_folds <- function(yr, ss=NULL, n_folds=NULL, ...){

  ## get unique years
  yrs <- sort(unique(yr))

  ## make the list
  chronos <- list(NULL)
  for(y in yrs){
    chronos[[which(yrs==y)]] <- which(yr <= y)
  }

  return(chronos)
}

###################
## IN SPACE-TIME ##
###################
## first we obey the temporal choice, then split in space
## currently, these are all combos of space holdout funcs and time holdout funcs.


##############################
## OVERALL WRAPPER FUNCTION ##
##############################
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##
## STEPS:
## 1st: subset data by strata if specified
## 2nd: fold by temporal strategy if chosen
## 3rd: fold by spatial  strategy if chosen
##
## INPUTS:
##
## data: full dataset to split into test/train
## n_folds: number of folds in test/train splits
## spat_strat: spatial  holdout strategy. one of: c('rand', 'poly', 'qt', 'ct')
## temp_strat: temporal holdout strategy. one of: c('rand', 'prop', 'yr', 'chrono')
## seed: RNG seed in case you'd like to be able to recreate folds
##
## OUTPUTS: 2 item list
## 1st item: 1 by nrow(data) vector containing integers identifying folds
## 2nd item: matrix of stratification combinations used to make holdouts
##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
make_folds <- function(data,
                       n_folds,
                       spat_strat = 'rand',
                       temp_strat = 'rand',
                       spte_strat =  NULL,
                       save.file = paste0('<<<< FILEPATH REDACTED >>>>>',
                                          indicator_group, '/',
                                          indicator, '/output/',
                                          run_date, '/stratum.rds'), ## place to save the final object from the function
                       ...,
                       strat_cols, ## needed for stratifying
                       seed
                       ){

  ##~#####################
  ## Prepare for battle ##
  ##~#####################

  ## for some reason, functions that work on data.frames don't all work on data.tables
  data <- as.data.frame(data)

  message('Identifying folds for validation')

  if(!is.null(spte_strat)){
    message("B/C you've chosen a space-time strategy, any input into either spat_strat or temp_strat strategies will be ignored")
    spat_strat=NA
    temp_strat=NA
  }

  if(temp_strat == "yr" |
     temp_strat == "chrono"){
    message("B/C you've chosen temp_strat==('yr'|'chrono'), I've set spat_strat to be NA")
    spat_strat=NA
  }

  message(paste("You've selected:",
                paste0("spat_strat = ", spat_strat),
                paste0("temp_strat = ", temp_strat),
                paste0("spte_strat = ", spte_strat),
                sep="\n"))

  ## check for unused arguments and give warning if any won't be used
  params <- list(...)
  optional_params <- c('ts', 'mb', 'plot_fn', 'plot_shp', 'shp_fn',
                       'admin_shps', 'shape_ident', 'admin_raster',
                       'ss_col', 'mask_shape', 'mask_raster',
                       'long_col', 'lat_col', 'yr_col', 'ts_vec')
  unused_params <- setdiff(names(params),optional_params)
  if(length(unused_params)){
    stop('You entered unused parameters! ', paste(unused_params,collapse = ', '))
  }

  ## set the seed if desired
  if(!(missing(seed))) set.seed(seed)


  ##~##############################
  ## FIRST, split data by strata ##
  ##~##############################
  print("Making stratum")

  if(!missing(strat_cols)){

    ## get different stratum
    all_strat <- get_strat_combos(data=data, strat_cols=strat_cols)

    ## make a list. each element is a different strata of data
    stratum <- list(NULL)
    for(s in 1:nrow(all_strat)){

      ## get all combos
      strata <- as.data.frame(all_strat[s, ])
      colnames(strata) <- colnames(all_strat)

      ## get rows in current strata
      strat_rows <- get_strat_rows(data=data,
                                   strata=strata)

      ## enter the data strata into the list
      stratum[[s]] <- data[strat_rows, ]
    }
  }else{
    stratum[[1]] <- data
  }

  m1 <- paste0("Length of data is: ", nrow(data))
  m2 <- paste0("Length of data in all strata is: ", sum(unlist(lapply(stratum, nrow))))
  print(m1)
  print(m2)
  print("These should be the same number. If not, maybe you have NAs to clean?")

  ##~#############################
  ## SECOND, make folds in time ##
  ##~#############################

  print("Making time folds")

  if(!is.na(temp_strat)){
    ## we have a dictionary of temporal methods to call from
    t_dict <- list('rand'  = NULL,
                   'prop'  = proptime_folds,
                   'yr'    = yr_folds,
                   'chrono'= chrono_folds
                   )

    ## match the choice to the function and run it
    t_fun   <- t_dict[[temp_strat]]
    if(!is.null(t_fun)){
      ## TODO Parallelize
      t_stratum <- lapply(1:length(stratum),
                          function(x){
                            yr <- stratum[[x]][[yr_col]]
                            ss <- stratum[[x]][[ss_col]]
                            ## get the fold indices
                            t_folds <- t_fun(yr=yr, ss=ss,
                                             n_folds=n_folds,
                                             ...)
                            ## make & return the subsetted datasets
                            lapply(1:length(t_folds),
                                   function(y){
                                     stratum[[x]][t_folds[[y]], ]
                                   })

                          })
    }else{
      t_stratum <- list(NULL)
      for(i in 1:length(stratum)){
        t_stratum[[i]]  <- stratum[[i]]
      }
    }
  }else{
    t_stratum <- list(NULL)
    for(i in 1:length(stratum)){
      t_stratum[[i]]  <- stratum[[i]]
    }
  }

  ## identify time folds in each stratum and unlist back to stratum
  stratum  <- NULL
  for(i in 1:length(t_stratum)){
    for(j in 1:length(t_stratum[[i]])){
      t_stratum[[i]][[j]] <- cbind(t_stratum[[i]][[j]],
                                   rep(j, nrow(t_stratum[[i]][[j]])))
      colnames(t_stratum[[i]][[j]])[ncol(t_stratum[[i]][[j]])] <- "t_fold"
    }
    stratum[[i]] <- do.call(rbind, t_stratum[[i]])
  }

  m1 <- paste0("Length of data is: ", nrow(data))
  m2 <- paste0("Length of data in all strata is: ", sum(unlist(lapply(stratum, nrow))))
  print(m1)
  print(m2)
  print("These should be the same number. If not, maybe you have NAs to clean?")


  ##~#############################
  ## THIRD, make folds in space ##
  ##~#############################

  print("Making space folds")

  if(!is.na(spat_strat)){
    ## we have a dictionary of spatial methods to call from
    s_dict <- list('rand' = rand_s_folds,
                   'poly' = ad2_folds,
                   'qt'   = quadtree_folds,
                   'ct'   = ct_folds
                   )

    ## match the choice to the function and run it on the double
    ## list of data
    s_fun <- s_dict[[spat_strat]]
    if(!is.null(s_fun)){
      ## TODO parallelize
      s_stratum <- lapply(1:length(stratum),
                          function(x){

                            message(paste0("Making folds for stratum: ", x))
                            ## first col is fold second col is holdout ID
                            s_fold_hoid <- matrix(ncol = 2, nrow = nrow(stratum[[x]]))

                            lapply(1:max(stratum[[x]][['t_fold']]),
                                   function(y){

                                     message(paste0("On temp fold: ", y))
                                     t_fold_r <- which(stratum[[x]][['t_fold']] == y)
                                     xy <- cbind(stratum[[x]][t_fold_r,
                                                              long_col],
                                                 stratum[[x]][t_fold_r,
                                                              lat_col])
                                     colnames(xy) <- c("long", "lat")
                                     ss <- stratum[[x]][t_fold_r, ss_col]

                                     ## IF 'ts_vec' has been input, as opposed to just 'ts', we need to reassign that here
                                     extra.args <- list(...)
                                     if(!is.null(extra.args[['ts_vec']])){
                                       extra.args[['ts']] <- extra.args[['ts_vec']][x]
                                       extra.args[['ts_vec']] <- NULL
                                     }
                                     all.args <- c(list('xy' = xy,
                                                        'ss' = ss,
                                                        'n_folds' = n_folds,
                                                        'stratum' = x,
                                                        't_folds' = y),
                                                   extra.args)

                                     ## get the fold indices
                                     s_folds <- do.call(s_fun, all.args)

                                     s_fold_hoid[t_fold_r, ]  <- s_folds

                                                      #...)

                                   })
                          })

      ## unpack fold indices back to stratum
      for(i in 1:length(stratum)){
        fold_hoid <- matrix(ncol = 2, nrow = nrow(stratum[[i]]))

        for(j in sort(unique(stratum[[i]][['t_fold']]))){
          t_fold_r <- which(stratum[[i]][['t_fold']] == j)
          fold_hoid[t_fold_r, ] <- s_stratum[[i]][[j]]
          colnames(fold_hoid) <- c("fold", "ho_id")
        }
        stratum[[i]] <- cbind(stratum[[i]], fold_hoid)
      }
    }
  }

  m1 <- paste0("Length of data is: ", nrow(data))
  m2 <- paste0("Length of data in all strata is: ", sum(unlist(lapply(stratum, nrow))))
  print(m1)
  print(m2)

  ##~######################################################
  ## FOURTH, if not the others, make folds in space-time ##
  ##~######################################################

  if(!is.null(spte_strat)){

    ## we have a dictionary of space-time methods to call from
    st_dict <- list()

    ## match the choice to the function and run it
    st_fun   <- st_dict[[spte_strat]]
    st_folds <- st_fun(xy=xy, yr=yr, ss=ss, n_folds=n_folds, ...)
  }

  ##~###########################################
  ## Lastly, if completely random is selected ##
  ##~###########################################



  ## ############################
  ## FINISH w/ post-processing ##
  ## ############################

  ## name the folds
  names(stratum) <- name_strata(all_strat = all_strat)

  ## save a copy of stratum
  saveRDS(file = save.file, object = stratum)

  ## return
  return(stratum)

}

#######################
## UTILITY FUNCTIONS ##
#######################

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## stratification functions
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

get_strat_combos <- function(data=data,
                             strat_cols=strat_cols,
                             ...){

  ## this function creates all unique combos that we need to stratify over when making folds
  ## e.g. if gender and age_bin vectors are specified this function
  ## will return a list with all unique combos of age_bin and gender

  ## get all unique items from  each column
  unique_list <- list(NULL)
  for (i in 1:length(strat_cols)){
    unique_list[[i]] <- sort(unique(data[, strat_cols[i]]))
  }

  ## make a dataframe of all combos and return it
  all_combos <- expand.grid(unique_list)
  colnames(all_combos) <- strat_cols
  return(all_combos)

}

get_strat_rows <- function(data=data,
                           strata,
                           ...){

  ## this function returns all the rows in a strata

  if(length(strata) < 1){
    message("Need to identify some strata!")
    stop()
  }

  ## loop through and intersect all rows we want
  good_rows <- data[, colnames(strata)[1] ] == strata[[1]]

  if(length(strata) > 1){
    for(c in 2:ncol(strata)){
      tmp_rows <- data[, colnames(strata)[c] ] == strata[[c]]
      good_rows <- good_rows * tmp_rows ## intersect them
    }
  }

  good_rows <- which(good_rows == 1)

  return(good_rows)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## shapefile countss
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
get_sample_counts <- function(ss = 1,
                              xy,
                              shapes,
                              shape_ident   = "gaul_code",
                              ...){
  ## this function takes in data and relevant shapefiles and returns
  ## how many of the datapoints falls into each of the shapefiles as
  ## well as the mapping vector between points and the shape_ident
  ## of the polygons

  ## data: full dataset
  ## ss: vector of sample sizes in each row. if <1> it assumes all have ss=1
  ## xy: matrix of xy coordinates
  ## shapes: all relevant loaded shapefiles (e.g loaded visa raster::shapefiles() )
  ## shape_ident: in a spatial polygons dataframe, there is data associated with each
  ##    entry. this variable identifies which of these data cols you'd like to use to
  ##    refer to the different polygons in the SPDF. if using admin2, leave as "gaul_code"
  ##    but, if you make your own set of shapes, you may want to select another col

  library(sp)

  ## grab relevant cols
  if(length(ss) == 1){
    data <- cbind(xy, rep(1, nrow(nrow))) ## samplesize
  }else{
    data <- cbind(xy, ss)
  }
  colnames(data) <- c("long", "lat", "ss")

  ## make sure all relevant cols are truly numeric
  for(i in 1:3){
    data[, i] <- as.numeric(as.character(data[, i]))
  }

  data <- as.data.frame(data)

  ## prepare the data to extract points into the shapefile
  coordinates(data) <- ~long+lat
  proj4string(data) <- proj4string(shapes)

  ## find which shapes the data falls into
  ## this takes the longest
  row_shapes <- over(data, shapes)
  row_shapes <- row_shapes[, shape_ident]

  ## WARNING! If some of your pts don't fit into the poly shapes
  ## they get randomly assigned to folds at the end
  if(sum(is.na(row_shapes)) > 0 ){
    message("Warning!! Some of your pts don't fall into any of the polygon shapes!!")
    message("They will get randomly assigned to folds")
    png("~/pts_outside_polys.png", width=1080, height=1080)
    plot(shapes)
    points(data[which(is.na(row_shapes)), ],col="red", pch=3)
    dev.off()
  }

  ## find the sample size in each shape
  all_shapes <- sort(unique(shapes@data[,shape_ident]))
  n_shapes   <- length(all_shapes)
  samp_size_shape <- cbind(all_shapes, rep(NA, n_shapes))
  for(i in 1:n_shapes){
    shape_id  <- samp_size_shape[i, 1]
    data_rows <- which(row_shapes==shape_id)
    samp_size_shape[i, 2] <- sum(data$ss[data_rows])
  }
  colnames(samp_size_shape) <- c(shape_ident, "count")

  return(list(ct_mat=samp_size_shape,
              pt_poly_map=row_shapes))
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## make approximate equal sample size folds by subsets
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
make_folds_by_poly <- function(cts_in_polys,
                               pt_poly_map,
                               n_folds,
                               ...){
  ## cts_in_polys:matrix containing two cols. 1st col is shape_ident
  ##    which identifies polygons in shapefile. 2nd col is count in shape
  ## pt_poly_map: vector as long as dataframe. each entry contains
  ##    the map from data to polygon shape_ident

  ## this function randomizes admin unit order and then sequentially
  ## adds to each fold until the sample size in the fold reaches
  ## 1/n_folds of the total count

  ## to make sure that the last one isn't much smaller, half of the
  ## folds take out the last poly that pushed them over the
  ## allowable count

  ## randomize the order

  if(n_folds > nrow(cts_in_polys)){
    message("You have too little data somewhere to split into ", n_folds, " folds")
    message("Check the sample size in each strata, and each strata after time holdouts")
    stop()
  }

  rand_ord <- sample(1:nrow(cts_in_polys))

  ## randomly decide if the first one will include the final poly
  flip <- sample(0:1, 1)
  if(flip==0){
    include_final <- rep(0:1, ceiling(n_folds/2))
  }else{
    include_final <- rep(1:0, ceiling(n_folds/2))
  }

  ## get the sample size threshhold in each poly
  total_ct <- sum(cts_in_polys[,2])
  max_fold_ct <- ceiling(total_ct/n_folds)

  ## add polys to folds
  fold_sums <- rep(NA, n_folds)
  start.ind <- 1
  stop.ind  <- 1
  for(fold in 1:(n_folds-1)){

    ## check threshhold
    while(sum(cts_in_polys[rand_ord[start.ind:stop.ind], 2]) < max_fold_ct){
      stop.ind <- stop.ind + 1
    }

    ## check if final poly is included
    if(include_final[fold]==1){
      stop.ind <- stop.ind - 1
    }

    ## store all the polys in the fold
    assign(paste0('polys_in_fold_', fold),
           cts_in_polys[rand_ord[start.ind:stop.ind], 1])

    ## record total in fold
    total_ct[fold] <- sum(cts_in_polys[rand_ord[start.ind:stop.ind], 2])

    ## adjust indices
    start.ind <- stop.ind + 1
    stop.ind  <- start.ind + 1
  }

  ## and the last fold is everything else
  stop.ind <- length(rand_ord)
  fold <- n_folds
  assign(paste0('polys_in_fold_', fold),
         cts_in_polys[rand_ord[start.ind:stop.ind], 1])
  total_ct[fold] <- sum(cts_in_polys[rand_ord[start.ind:stop.ind], 2])

  message("The sum in each different fold is: \n")
  for(i in 1:n_folds){
    message(total_ct[i])
  }

  ## now we have the polys in different folds. we make a
  ## vector for which data rows are in the folds
  fold_vec <- rep(NA, length(pt_poly_map))
  for(fold in 1:n_folds){
    fold_rows <- which(pt_poly_map %in% get( (paste0('polys_in_fold_', fold) )))
    fold_vec[fold_rows] <- fold
  }

  ## lastly, some of the points may not have fallen in the shapefiles
  ## for the moment they get randomly assigned to folds
  if(sum(is.na(pt_poly_map)) > 0){
    last_folds <- cut(seq(1, sum(is.na(pt_poly_map))),
                      breaks = n_folds, labels = 1:n_folds)
    last_folds <- sample(as.numeric(as.character(last_folds)))
    fold_vec[which(is.na(pt_poly_map))] <- last_folds
  }

  ## return a vector placing each data row in a fold
  return(fold_vec)

}


## make strata names for final list of holdout datasets
name_strata  <-  function(all_strat){
  ## takes in a matrix containing combinations of all strata
  ##      (e.g. the output from get_strat_combos() )
  ## returns a vector containing the names of the strata
  ## NAMING CONVENTION: stratcol1__strata1___stratcol2__strata2___...


  as <- all_strat
  strat_names <- paste(colnames(as)[1], as[, 1], sep = "__")
  if(ncol(as) > 1){
    for(i in 2:ncol(as)){
      strat_names <- paste(strat_names,
                           paste(colnames(as)[i], as[, i], sep = "__"),
                           sep = "___")
    }
  }

  return(strat_names)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## quadtree functions ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## adapted from:
## http://gis.stackexchange.com/questions/31236/how-can-i-generate-irregular-grid-containing-minimum-n-points

## quadtree by points
quadtree <- function(xy, k=1) {
  d = dim(xy)[2]
  quad <- function(xy, i, id=1) {
    if (nrow(xy) < k*d) {
      rv = list(id=id, value=xy)
      class(rv) <- "quadtree.leaf"
    }
    else {
      q0 <- (1 + runif(1,min=-1/2,max=1/2)/dim(xy)[1])/2 # Random quantile near the median
      x0 <- quantile(xy[,i], q0)
      j <- i %% d + 1 # (Works for octrees, too...)
      rv <- list(index=i, threshold=x0,
                 lower=quad(xy[xy[,i] <= x0, ], j, id*2),
                 upper=quad(xy[xy[,i] > x0, ], j, id*2+1))
      class(rv) <- "quadtree"
    }
    return(rv)
  }
  quad(xy, 1)
}

## my hacked version. quadtree by sum of value at each point.
## is ss=1 for all points you get back to original quadtree
quadtree_ct <- function(xy, ss, target_ss, min_in_bin=5, rand = T) {
  ## this function quadtrees by sample size
  ## rand introduces some randomness to the "median"
  print(paste0("Aiming for ts between: ", target_ss, " and ", 2 * target_ss))
  d = dim(xy)[2]
  quad <- function(xy, i, id=1, ss, target_ss) {
    if (sum(ss) < target_ss * 2 | length(xy)/2 <= min_in_bin){
      rv = list(id=id, value=xy)
      class(rv) <- "quadtree.leaf"
    }
    else {
      if(rand){
        q0 <- (1 + runif(1,min=-1/10,max=1/10)/dim(xy)[1])/2 # Random quantile near the median
      }else{
        q0 <- 1 / 2 ## no randomness, just the median
      }
      x0 <- quantile(xy[,i], q0)
      j <- i %% d + 1 # (Works for octrees, too...)
      rv <- list(index=i, threshold=x0,
                 lower=quad(xy[xy[,i] <= x0, ], j, id*2,  ss[xy[,i] <= x0], target_ss),
                 upper=quad(xy[xy[,i] > x0, ], j, id*2+1, ss[xy[,i] >  x0], target_ss))
      class(rv) <- "quadtree"
    }
    return(rv)
  }
  quad(xy=xy, i=1, id=1, ss=ss, target_ss=target_ss)
}

## plotting functions
points.quadtree <- function(q, alpha=0.1, ...) {
  points(q$lower, alpha, ...); points(q$upper, alpha, ...)
}

points.quadtree.leaf <- function(q, alpha=0.1, ...) {
  library(scales)
  points(q$value, col=alpha(q$id, alpha=alpha), ...)
}

lines.quadtree <- function(q, xylim, ...) {
  i <- q$index
  j <- 3 - q$index
  clip <- function(xylim.clip, i, upper) {
    if (upper) xylim.clip[1, i] <- max(q$threshold, xylim.clip[1,i]) else
      xylim.clip[2,i] <- min(q$threshold, xylim.clip[2,i])
    xylim.clip
  }
  if(q$threshold > xylim[1,i]) lines(q$lower, clip(xylim, i, FALSE), ...)
  if(q$threshold < xylim[2,i]) lines(q$upper, clip(xylim, i, TRUE), ...)
  xlim <- xylim[, j]
  xy <- cbind(c(q$threshold, q$threshold), xlim)
  lines(xy[, order(i:j)],  ...)
}
lines.quadtree.leaf <- function(q, xylim, ...) {} # Nothing to do at leaves!


## get cell boundaries to make shapefiles
cell <- function(q, xylim, ...) {
  if (class(q)=="quadtree") f <- cell.quadtree else f <- cell.quadtree.leaf
  f(q, xylim, ...)
}
cell.quadtree <- function(q, xylim, ...) {
  i <- q$index
  j <- 3 - q$index
  clip <- function(xylim_clip, i, upper) {
    if(upper){
      xylim_clip[1, i] <- max(q$threshold, xylim_clip[1,i])
    }else{
      xylim_clip[2,i] <- min(q$threshold, xylim_clip[2,i])
    }
    xylim_clip
  }
  d <- data.frame(id=NULL, x=NULL, y=NULL)
  if(q$threshold > xylim[1,i]) d <- cell(q$lower, clip(xylim, i, FALSE), ...)
  if(q$threshold < xylim[2,i]) d <- rbind(d, cell(q$upper, clip(xylim, i, TRUE), ...))
  d
}
cell.quadtree.leaf <- function(q, xylim) {
  data.frame(id = q$id,
             x = c(xylim[1,1], xylim[2,1], xylim[2,1], xylim[1,1], xylim[1,1]),
             y = c(xylim[1,2], xylim[1,2], xylim[2,2], xylim[2,2], xylim[1,2]))
}


## get ids for each point
id <- function(q, ...){
  d <- data.frame(id=NULL, x=NULL, y=NULL)
  if (class(q)=="quadtree") f <- id.quadtree else f <- id.quadtree.leaf
  f(q, xylim, ...)
}
id.quadtree <- function(q, ...) {
  rbind(id(q$lower), id(q$upper))
}
id.quadtree.leaf <- function(q, ...) {
  ## print(q$id) ## for debugging
  if(length(q$value)==0){
    data.frame(id = q$id,
               x  = NA,
               y  = NA)
  }else if(length(q$value)==2){
    data.frame(id = q$id,
               x  = q$value[1],
               y  = q$value[2])
  }else{
    data.frame(id = q$id,
               x  = q$value[, 1],
               y  = q$value[, 2])
  }
}
