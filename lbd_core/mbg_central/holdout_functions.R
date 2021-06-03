#' @title Create Random Folds
#' @description Randomly divide data into folds, optionally within strata
#'
#' @param data data to be divided into folds. each row is an observation
#' @param n_folds number of folds to divide data into
#' @param strat_cols vector with names of columns to stratify by
#' 
#' @return returns a list containing 1) vector of folds and 2) matrix containing all stratum
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

  if(is.null(strat_cols)){
    return(list(folds   = fold_vec,
                stratum = NA))
  }else{
    return(list(folds   = fold_vec,
                stratum = all_strat))
  }
}

#' @title Create Quadtree Folds
#' @description This function segregates data into spatial regions. It splits 
#' data down until the target sample size is reached or until a minimum allowable 
#' number of points are in the bin. It then breaks the data into folds by 
#' ensuring that all data in any one spatial region stays together and such that 
#' the sample size sums in each fold are relatively similar. Points at the same 
#' place pose a problem because you can't split them so, first we make a subset 
#' of the data using only unique xy combos. Quadtree is run on the unique list 
#' and then we buid back up to the full dataset and make folds.
#'
#' @param xy 2 col matrix of xy locs
#' @param ss vector of sample size at each loc - if all 1s, it aggregates by # of points
#' @param ts target sample size in each region
#' @param n_folds number of folds
#' @param plot_fn if not null plots data and quad tree must be name of file ending in .png
#' @param plot_shp shapefile to add shapefile outlines to plot
#' @param save_qt T/F save quadtree regions to shapefiles
#' @param t_folds which t_fold is this, used when saving shapefiles
#' @param stratum which stratum is this, used when saving shapefiles
#' 
#' @return Matrix where first column is the fold and the second column is the quadtree id
#' 
#' @examples
#' \dontrun{
#' df <- fread(file.path(fp_list$temp_root, '/geospatial/U5M_africa/data/clean/fully_processed.csv',
#'             stringsAsFactors = FALSE))
#' df$long <- as.numeric(as.character(gsub(",", "", df$long)))
#' df$lat  <- as.numeric(as.character(gsub(",", "", df$lat)))
#' df <- df[-which(df$lat > 90), ]
#' data <- df
#' plot_shp <- shapefile(file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/shapefiles/africa_ad2.shp'))
#' xy <- data[,c('long', 'lat')]
#' ss <- data$exposed
#' ts <- 1e6
#' plot_fn <- 'quadtree.png'
#' 
#' quadtree_folds(xy = xy,
#'                ss = ss,
#'                ts = ts,
#'                n_folds = 5,
#'                plot_shp = plot_shp,
#'                plot_fn = plot_fn
#'                )
#' }
quadtree_folds <- function(xy,
                           ss, 
                           ts, 
                           n_folds,
                           plot_fn  = NULL,
                           plot_shp = NULL, 
                           save_qt  = TRUE, 
                           ...,
                           t_folds = 1,
                           stratum = 1
                           ){

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
  ## TODO: make this faster somehow... Roy suggested a forwardfill option
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
      plot(plot_shp, xlab="x", ylab="y", main=title)
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
#    if(time_stamp==TRUE) output_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/', run_date)
#    if(time_stamp==FALSE) output_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/scratch')
#    dir.create(output_dir, showWarnings = FALSE)
    output_dir <- paste0(fp_list['mbg_root'], indicator_group, '/', indicator, '/output/', run_date)
    qt_shape_dir <- paste0(output_dir, '/holdout_shapefiles/')
    dir.create(qt_shape_dir)

    xylim <- cbind(x=c(min(un_xy[,2]), max(un_xy[,2])), y=c(min(un_xy[,3]), max(un_xy[,3])))
    polys <- cell(qt, xylim)
    polys_attr <- data.frame(id=unique(polys$id))
    library(shapefiles)
    polys_shapefile <- convert.to.shapefile(polys, polys_attr, "id", 5)
    write.shapefile(polys_shapefile, paste0(qt_shape_dir, 'spat_holdout_stratum_', stratum, '_t_fold', t_folds), arcgis=TRUE)
  }

  return(cbind(fold_vec, ho_id))
}

#' @title Create Folds by Admin2
#' @description This function takes in admin2 (or any mutually exclusive and
#' collectively exhaustive) shapefiles covering the domain of your data and 
#' splits your data into folds of approximately equal sample size using admin2 
#' units to split the data
#'
#' @param admin_raster file location of all pertinent shapefiles to use when folding (.grd)
#' @param shape_ident string identifying data col in shapefile used to uniquely identify polygons
#' @param admin_shps file location of associated raster for admin_raster (.shp)
#' @param ss vector of sample sizes for each row. if <1>, assumes all have ss=1
#' @param xy matrix of xy coordinates
#' @param n_folds number of folds to make
#' @param mask_shape shapefile file location for boundary of area to be folded
#' @param mask_raster raster to project mask_shape onto
#' 
#' @return Matrix where first column is the fold and the second column is the admin2 id
#' 
#' @examples 
#' \dontrun{
#' df <- read.csv(file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/fully_processed.csv'),
#'                stringsAsFactors = FALSE)
#' df$long <- as.numeric(as.character(gsub(",", "", df$long)))
#' df$lat  <- as.numeric(as.character(gsub(",", "", df$lat)))
#' df <- df[-which(df$lat>90), ]
#' shp_full <- shapefile(file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/shapefiles/africa_ad2.shp'))
#' folds <- ad2_folds(admin_raster=file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/shapefiles/ad2_raster.grd'),
#'                    shape_ident="gaul_code",
#'                    admin_shps=file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/shapefiles/africa_ad2.shp'),
#'                    data=df,
#'                    strat_cols=NULL,
#'                    ss=data$exposed,
#'                    xy=cbind(df$long, df$lat),
#'                    n_folds=5,
#'                    mask_shape=file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/shapefiles/africa_simple.shp'),
#'                    mask_raster=file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/shapefiles/ad0_raster'))
#' library(scales)
#' cols <- folds
#' cols[which(cols==1)] <- "cyan"
#' cols[which(cols==2)] <- "red"
#' cols[which(cols==3)] <- "blue"
#' cols[which(cols==4)] <- "green"
#' cols[which(cols==5)] <- "magenta"
#' png("~/check_folds_plot.png", width=1080, height=1080)
#' plot(shp_full)
#' points(df$long, df$lat, col=alpha(cols, alpha=0.01), pch=16)
#' dev.off()
#' }
ad2_folds <- function(admin_raster,
                      shape_ident="gaul_code",
                      admin_shps,
                      ss=1,
                      xy,
                      n_folds,
                      mask_shape,
                      mask_raster,
                      ...){ 

  library(raster)

  ## make a mask for ther region we care about
  mask <- rasterize_check_coverage(shapefile(mask_shape), 
                                   raster(mask_raster), 
                                   field=names(shapefile(mask_shape))[1])*0

  ## get raster cells in mask
  cell_idx <- cellIdx(mask)

  ## load raster and shapefile for admin units
  rast      <- raster(admin_raster)
  rast_cell <- raster::extract(rast, cell_idx)
  shp_full  <- shapefile(admin_shps)
  shp       <- shp_full@data[shape_ident]
  ## plot(shp_full, col=1:length(shp_full))

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

#' @title Create Folds by Country
#' @description This function takes in admin2 (or any mutually exclusive and
#' collectively exhaustive) shapefiles covering the domain of your data and 
#' splits your data into folds of approximately equal sample size using admin2 
#' units to split the data
#'
#' @param xy matrix of xy coordinates
#' @param ct country vector
#' @param ss vector of sample sizes for each row. if <1>, assumes all have ss=1
#' @param n_folds number of folds to make
#' 
#' @return Vector with fold id
#' 
#' @examples 
#' \dontrun{
#' df <- read.csv(file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/fully_processed.csv'),
#'                 stringsAsFactors = FALSE)
#' shp_full <- shapefile(file.path(fp_list$temp_root, 'geospatial/U5M_africa/data/clean/shapefiles/africa_ad2.shp'))
#' df$long <- as.numeric(as.character(gsub(",", "", df$long)))
#' df$lat  <- as.numeric(as.character(gsub(",", "", df$lat)))
#' df$country <- gsub(pattern='Guinea-Bissau\\"',replacement="Guinea-Bissau", df$country)
#' df <- df[-which(df$lat>90), ]
#' data <- df
#' xy <- data[,c('long', 'lat')]
#' ss <- data$exposed
#' ct <- data$country
#' folds <- ct_folds(xy, ct, ss)
#' plot(shp_full)
#' library(scales)
#' points(xy, col=alpha(folds, alpha=0.25), pch=".")
#' }
ct_folds <- function(xy,
                     ct,  
                     ss=1,
                     n_folds=5,
                     ...){

  if(length(unique(ct)) < n_folds){
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


#' @title Create Random Folds in Space
#' @description This function randomly assigns each lat-long point to a fold
#'
#' @param xy matrix of xy coordinates
#' @param ss vector of sample sizes for each row. if <1>, assumes all have ss=1
#' @param n_folds number of folds to make
#' 
#' @return vector with fold id
rand_s_folds <- function(xy,   ## xy location matrix
                         ss=1, ## sample size vec (or 1 if equi-ss)
                         n_folds=5,
                         ...){

  if(length(ss) < n_folds){
    message("Too many folds for too few countries! Expand your horizons")
    stop()
  }

  if(length(ss)==1) ss <- rep(1, nrow(xy))
  
  if(length(ss) != nrow(xy)){
    message("length of ss and nrow(xy) must match!")
    stop()
  }

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

#' @title Create Proportional Folds in Time
#' @description This function takes in which time period each observation is in
#' and returns a list of length equal to number of unique time periods
#'
#' @param yr vector of time period for each observation
#' 
#' @return list of folds where each is a vector of which observations correspond to that fold
proptime_folds <- function(yr, ...){

  ## so all we have to do is return a 'fold vector' with integers for unique years

  yrf  <- as.factor(yr)
  fold_vec <- as.numeric(yrf)

  fold_list <- list(NULL)
  for(i in sort(unique(fold_vec))){
    fold_list[[i]] <- which(fold_vec == i)
  }
  return(fold_list)

}

#' @title Create Folds Across Time Periods
#' @description This function creates folds chronologically. Starting with the 
#' first year not yet assigned a fold, subsequent years are combined until a 
#' minimum sample size and unique observations threshold is met. The minimum 
#' year where the cumlative totals are sufficiently high is found and then all 
#' available years less than or equal to that year are assigned that fold. 
#'
#' @param yr vector of time period for each observation
#' @param n_folds number of folds
#' @param ss vector of sample sizes for each observation 
#' @param ts target sample size in each fold
#' 
#' @return list of folds where each is a vector of which observations correspond 
#' to that fold
proptime_folds_combine <- function(yr, n_folds, ss, ts, ...) {

  # count up the sample size by year, and the number of unique observations by year
  totals <- data.table(year = sort(unique(yr)),
                       total_n = tapply(ss, yr, sum),
                       total_poly = tapply(ss, yr, length))

  # starting with the first year, combine years until a minimum sample size and unique observations threshold is met.
  totals[, t_fold := as.numeric(NA)]
  fold <- 0

  while(sum(is.na(totals$t_fold)) > 0) {
    # iterate the fold number
    fold <- fold + 1

    # recalculate cumulative totals for all years that haven't been assigned a fold yet
    totals[is.na(t_fold), c("cum_total_n", "cum_total_poly") := list(cumsum(total_n), cumsum(total_poly))]

    # find the minimum year where the cumlative totals are sufficiently high
    yy <- totals[cum_total_n > (n_folds * ts * 2) & cum_total_poly > n_folds,][1, year]

    # assign all years not already assigned and less than or equal to this year to a new fold
    if (!is.na(yy)) {
      totals[is.na(t_fold) & year <= yy, t_fold := fold]

    # if there is no unassigned year where we reach the threshold (ie, if we run out of time), combine remaining years into the previous fold
    } else {
      totals[is.na(t_fold), t_fold := fold - 1]
    }

    totals[, c("cum_total_n", "cum_total_poly") := NULL]
  }

  # create the folds based on this mapping
  fold_list <- lapply(1:max(totals$t_fold), function(fold) which(yr %in% totals[t_fold == fold, year]))
  return(fold_list)

}

#' @title Create Folds by Year
#' @description Deprecated? I can't seem to get this function to work.
#'
#' @param xy xy location matrix
#' @param yr vector of time period for each observation
#' @param ss vector of sample sizes for each observation
#' @param n_folds number of folds
#' @param ct ?? vector of countries for each observation
#' 
#' @return list of folds where each is a vector of which observations correspond 
#' to that fold
yr_folds <- function(xy,  
                     yr,   
                     ss=1, 
                     n_folds=5,
                     ct,
                     ...){

  if(length(unique(yr)) < n_folds){
    message("Too many folds for too few years! Try again in a decade")
    stop()
  }

  if(length(ss)==1) ss <- rep(1, nrow(xy))

  ## first we find the sample size in each of the countries
  library(data.table)
  dt <- data.table(long = xy[,1],
                   lat  = xy[,2],
                   ss   = ss,
                   yr   = yr,
                   ct   = ct)

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

#' @title Create Folds Chronologically in Time
#' @description Build up different datasets as if we were moving chronologically 
#' in time. E.g. 2000. 2000 & 2005. 2000 & 2005 & 2010. ...
#'
#' @param yr vector of time period for each observation
#' 
#' @return list with as many entries as unique years each item in the list 
#' contains all data rows that should be in that fold
chrono_folds <- function(yr,  
                         ...){

  ## get unique years
  yrs <- sort(unique(yr))

  ## make the list
  chronos <- list(NULL)
  for(y in yrs){
    chronos[[which(yrs==y)]] <- which(yr <= y)
  }

  return(chronos)
}

#' @title Make Folds
#' @description Wrapper function for creating folds in data to use as holdout in
#' validation. Steps: 1) subset data by strata if specified. 2) fold by space-time
#' strategy if chosen 3) fold by temporal strategy if chosen. 4) fold by spatial 
#' strategy if chosen.
#'
#' @param data full dataset to split into test/train
#' @param n_folds number of folds in test/train splits
#' @param spat_strat spatial  holdout strategy. one of: c('rand', 'poly', 'qt', 'ct')
#' @param temp_strat temporal holdout strategy. one of: c('rand', 'prop', 'prop_comb', 'yr', 'chrono')
#' @param spte_strat spatio-temporal holdout stratgey. Can only be 'nids'
#' @param indicator indicator
#' @param indicator_group group indicator belongs to
#' @param run_date run date
#' @param save.file place to save the final object from the function
#' @param strat_cols vector of columns to stratify by
#' @param n_folds number of folds
#' @param ct vector of countries for each observation
#' @param seed RNG seed in case you'd like to be able to recreate folds
#' 
#' @return 2 item list: 1) 1 by nrow(data) vector containing integers identifying folds
#' and 2) matrix of stratification combinations used to make holdouts
#' 
#' @examples 
#' \dontrun{
#' n <- 500
#' yr <- c(2000, 2005, 2010, 2015)
#' xy <- matrix(runif(n * length(yr) * 2), ncol = 2)
#' ss <- sample(1:100, size = n * length(yr), replace = TRUE)
#' age <- sample(1:4, size=n*length(yr), replace=TRUE)
#' yrs <- sample(yr, size = n * length(yr), replace = TRUE)
#' df <- as.data.frame(cbind(xy, ss, yrs, age))
#' colnames(df) <- c("longitude", "latitude", "ss", "year", "age")
#' ## Note: this does not run wihtout indicator, indicator_group, run_date, and
#' ## save.file specified
#' folds = make_folds(data = df, n_folds = 5, spat_strat = 'rand', temp_strat = 'prop',
#'                   long_col = 'longitude', lat_col = 'latitude', strat_cols = 'age')

#' }
make_folds <- function(data,
                       n_folds,
                       spat_strat = NULL,
                       temp_strat = NULL,
                       spte_strat = NULL,
                       indicator = use_global_if_missing("indicator"),
                       indicator_group = use_global_if_missing("indicator_group"),
                       run_date = use_global_if_missing("run_date"),
                       save.file = paste0(fp_list['mbg_root'],
                                          indicator_group, '/',
                                          indicator, '/output/',
                                          run_date, '/stratum.rds'),
                       ...,
    #                   long_col,   ## needed for any spat_strat
    #                   lat_col,    ## needed for any spat_strat
    #                   ss_col,     ## needed for most strats
    #                   yr_col,     ## needed for most temp_strats
    #                   ct_col,     ## needed for ct spat_strat
                       strat_cols, 
                       seed
                       ){

  ## setup for testing purposes
  ## n_folds = 5
  ## spat_strat = "qt"
  ## temp_strat = "prop"
  ## spte_strat = NULL
  ## strat_cols = "age_bin"
  ## yr_col = "yr"
  ## ss_col = "exposed"
  ## lat_col = "lat"
  ## long_col = "long"
  ## ts = 1000


  ## ## make a dataset
  ## n = 10000
  ## data = data.frame(lat  = runif(n, 10, 30),
  ##                   long = runif(n, 0, 40),
  ##                   yr   = sample(rep(c(1995, 2000, 2005, 2010), n / 4), n),
  ##                   age_bin = sample(1:4, n, replace = T),
  ##                   exposed = sample(25:35, n, replace = T))

  ##~#####################
  ## Prepare for battle ##
  ##~#####################

  ## for some reason, functions that work on data.frames don't all work on data.tables
  data <- as.data.frame(data)

  message('Identifying folds for validation')

  if (!is.null(spte_strat)) {
    message("B/C you've chosen a space-time strategy, any input into either spat_strat or temp_strat strategies will be ignored")
    spat_strat <- NULL
    temp_strat <- NULL
  }

  if (!is.null(temp_strat)) {
    if (temp_strat == "yr" | temp_strat == "chrono") {
      message("B/C you've chosen temp_strat==('yr'|'chrono'), I've set spat_strat to be NA")
      spat_strat <- NULL
    }
  }

  message(paste("You've selected:",
                paste0("spat_strat = ", spat_strat),
                paste0("temp_strat = ", temp_strat),
                paste0("spte_strat = ", spte_strat),
                sep="\n"))

  ## check for unused arguments and give warning if any won't be used
  params <- list(...)
  optional_params <- c('ts', 'plot_fn', 'plot_shp', 'shp_fn',
                       'admin_raster', 'shape_ident', 'admin_shps',
                       'ss_col', 'mask_shape', 'mask_raster',
                       'long_col', 'lat_col', 'yr_col', 'ts_vec')
  unused_params <- setdiff(names(params),optional_params)
  if(length(unused_params)){
    stop('You entered unused parameters! ', paste(unused_params,collapse = ', '))
  }

  #if additional parameters were passed to the function, assign them to the local environment
  if(length(params) > 0) {
    list2env(params, environment())
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

  if (!is.null(temp_strat)) {
    print("Making time folds")

    ## we have a dictionary of temporal methods to call from
    t_dict <- list('rand'      = NULL,
                   'prop'      = proptime_folds,
                   'prop_comb' = proptime_folds_combine,
                   'yr'        = yr_folds,
                   'chrono'    = chrono_folds)

    ## match the choice to the function and run it
    t_fun   <- t_dict[[temp_strat]]
    if (!is.null(t_fun)) {
      ## TODO Parallelize
      t_stratum <- lapply(1:length(stratum),
                          function(x){
                            yr <- stratum[[x]][[yr_col]]
                            ss <- stratum[[x]][[ss_col]]
                            ## get the fold indices
                            t_folds <- t_fun(yr=yr,
                                             ss=ss,
                                             n_folds = n_folds,
                                             ...)
                            ## make & return the subsetted datasets
                            lapply(1:length(t_folds),
                                   function(y){
                                     stratum[[x]][t_folds[[y]], ]
                                   })

                          })
    } else {
      t_stratum <- lapply(stratum, function(x) list(x))
    }

  } else {
    t_stratum <- lapply(stratum, function(x) list(x))
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

  if (!is.null(spat_strat)) {
    print("Making space folds")

    ## we have a dictionary of spatial methods to call from
    s_dict <- list('rand' = rand_s_folds,
                   'poly' = ad2_folds,
                   'qt'   = quadtree_folds,
                   'ct'   = ct_folds)

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
                                                        't_folds' = y,
                                                        "indicator" = indicator,
                                                        "indicator_group" = indicator_group,
                                                        "run_date" = run_date),
                                                   extra.args)

                                     ## get the fold indices
                                     s_folds <- do.call(s_fun, all.args)

                                     s_fold_hoid[t_fold_r, ]  <- s_folds

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

  if(!is.null(spte_strat)) {

    # get space-time folds
    if (spte_strat == "nids") {
      st_folds <- lapply(stratum, function(x) {
        nid_folds(nids = x$nid, yr = x[[yr_col]], ss = x[[ss_col]], n_folds)
      })

    } else {
      stop(paste(spte_strat, "is not a valid option for spte_strat"))
    }

    # copy fold IDs back to stratum
    for (i in 1:length(stratum)) {
      for (j in 1:length(st_folds[[i]])) {
        stratum[[i]][st_folds[[i]][[j]], "fold"] <- j
      }
    }
  }

  ##~###########################################
  ## Lastly, if completely random is selected ##
  ##~###########################################



  ## ############################
  ## FINISH w/ post-processing ##
  ## ############################

  ## recombine temporal stratification if needed (i.e. if )

  ## check it.
  ## table(data$fold,data$age_bin,dnn=list('fold','age_bin'))

  ## ## clean it up
  ## data=data[order(as.numeric(row.names(data))),]
  ## data$elig=NULL

  ## ## save to disk
  ## write.csv(data,
  ##           file = paste0('data/clean/mortality_combined.csv'),
  ##           row.names = FALSE)

  ## name the folds
  names(stratum) <- name_strata(all_strat = all_strat)

  ## save a copy of stratum
  saveRDS(file = save.file, object = stratum)

  ## return
  return(stratum)

}

#' @title Get All Combinations of Stratifying Columns
#' @description Creates a data frame of all combination of stratifying columns.
#' E.g. if gender and age_bin vectors are specified this function will return a 
#' list with all unique combos of age_bin and gender
#'
#' @param data full dataset to split into test/train
#' @param strat_cols vector of columns to stratify by
#' 
#' @return data frame containing one row for each combination of strata
get_strat_combos <- function(data=data,
                             strat_cols=strat_cols,
                             ...){

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

#' @title Get All Rows in a Stratum
#' @description Creates a vector with the rows in the data that correspond to 
#' a particular combination of strata
#'
#' @param data full dataset to split into test/train
#' @param strata datafame with one row containing the values for a partiuclar 
#' combination of strata (columns are the strata)
#' 
#' @return vector of rows in data that have that combination of strata
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



#' @title Get Sample Counts
#' @description This function takes in data and relevant shapefiles and returns
#' how many of the datapoints falls into each of the shapefiles as well as the 
#' mapping vector between points and the shape_ident of the polygons
#'
#' @param ss vector of sample sizes in each row. if <1> it assumes all have ss=1
#' @param xy matrix of xy coordinates
#' @param shapes all relevant loaded shapefiles (e.g loaded visa raster::shapefiles() )
#' @param shape_ident in a spatial polygons dataframe, there is data associated with each
#' entry. this variable identifies which of these data cols you'd like to use to
#' refer to the different polygons in the SPDF. if using admin2, leave as "gaul_code"
#' but, if you make your own set of shapes, you may want to select another col
#' 
#' @return list of number of observations that fall into each of the shapefiles 
#' and a vector that maps between the observations and shape_ident of the polygon
get_sample_counts <- function(ss = 1,
                              xy,
                              shapes,
                              shape_ident   = "gaul_code",
                              ...){
  library(sp)

  data <- cbind(xy, ss)
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
  ## TODO: map these to nearest polys somehow
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


#' @title Make Folds by Polygon
#' @description This function randomizes admin unit order and then sequentially
#' adds to each fold until the sample size in the fold reaches 1/n_folds of the 
#' total count. To make sure that the last one isn't much smaller, half of the
#' folds take out the last poly that pushed them over the allowable count
#'
#' @param cts_in_polys matrix containing two cols. 1st col is shape_ident which 
#' identifies polygons in shapefile. 2nd col is count in shape
#' @param pt_poly_map vector as long as dataframe. each entry contains the map 
#' from data to polygon shape_ident
#' @param n_folds number of folds
#' 
#' @return vector placing each data row in a fold
make_folds_by_poly <- function(cts_in_polys,
                               pt_poly_map,
                               n_folds,
                               ...){
  ## randomize the order

  if(n_folds > sum(as.numeric(cts_in_polys[, 2]) > 0)){
    message("You have too little data somewhere to split into ", n_folds, " folds")
    message("Check the sample size in each strata, and each strata after time holdouts")
    message("Still, we will assign the available data randomly to folds...")

    fold_vec <- sample(x = 1:n_folds, size =  sum(as.numeric(cts_in_polys[, 2]) > 0), replace = F)

    message("The sample size sum in each fold is: \n")
    pp.ind <- 1
    for(i in 1:n_folds){
      if(i %in% fold_vec){
        message(sprintf('Fold %i: %s', i, cts_in_polys[which(cts_in_polys[, 1] == pt_poly_map[pp.ind]), 2]))
        pp.ind <- pp.ind + 1
      }else{
        message('0')
      }
    }

    return(fold_vec)
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
  total_ct <- sum(as.numeric(cts_in_polys[,2]))
  max_fold_ct <- ceiling(total_ct/n_folds)

  ## add polys to folds
  fold_sums <- rep(NA, n_folds)
  start.ind <- 1
  stop.ind  <- 1
  for(fold in 1:(n_folds-1)){

    ## check threshhold
    while(sum(as.numeric(cts_in_polys[rand_ord[start.ind:stop.ind], 2])) < max_fold_ct &
          stop.ind < nrow(cts_in_polys)){
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
    total_ct[fold] <- sum(as.numeric(cts_in_polys[rand_ord[start.ind:stop.ind], 2]))

    ## adjust indices
    start.ind <- stop.ind + 1
    stop.ind  <- start.ind + 1
  }

  ## and the last fold is everything else
  stop.ind <- length(rand_ord)
  fold <- n_folds
  assign(paste0('polys_in_fold_', fold),
         cts_in_polys[rand_ord[start.ind:stop.ind], 1])
  total_ct[fold] <- sum(as.numeric(cts_in_polys[rand_ord[start.ind:stop.ind], 2]))

  message("The sample size sum in each fold is: \n")
  for(i in 1:n_folds){
    message(sprintf('Fold %i: %s', i, total_ct[i]))
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
  ## TODO: map these to nearest polys somehow
  if(sum(is.na(pt_poly_map)) > 0){
    last_folds <- cut(seq(1, sum(is.na(pt_poly_map))),
                      breaks = n_folds, labels = 1:n_folds)
    last_folds <- sample(as.numeric(as.character(last_folds)))
    fold_vec[which(is.na(pt_poly_map))] <- last_folds
  }

  ## return a vector placing each data row in a fold
  return(fold_vec)

}


#' @title Make Strata Names for Final List of Holdout Datasets
#' @description This function randomizes admin unit order and then sequentially
#' adds to each fold until the sample size in the fold reaches 1/n_folds of the 
#' total count. To make sure that the last one isn't much smaller, half of the
#' folds take out the last poly that pushed them over the allowable count
#'
#' @param all_strat matrix containing combinations of all strata (e.g. the 
#' output from get_strat_combos() )
#' 
#' @return vector containing the names of the strata, e.g.
#'  strat1name__strat1val1___strat2name__strat2val1___...,
#'  strat1name__strat1val5___strat2name__strat2val3___...
name_strata  <-  function(all_strat){


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

#' @title Quadtree
#' @description Not used?? Code from:
#'  http://gis.stackexchange.com/questions/31236/how-can-i-generate-irregular-grid-containing-minimum-n-points
#'
#' @param xy matrix of xy locations
#' @param k 
#' 
#' @return 
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

#' @title Quadtree by Sample Size
#' @description Based on code from:
#'  http://gis.stackexchange.com/questions/31236/how-can-i-generate-irregular-grid-containing-minimum-n-points
#'
#' @param xy matrix of xy locations
#' @param ss sample size
#' @param target_ss lower bound of target sample size in a bin
#' @param min_in_bin minimum number of points in a bin
#' @param rand should there be randomness to the "median"
#' 
#' @return Quadtree object
#' 
#' @examples 
#' \dontrun{
#' n <- 20          # Points per cluster
#' n_centers <- 10  # Number of cluster centers
#' sd <- 1/2        # Standard deviation of each cluster
#' set.seed(17)
#' centers <- matrix(runif(n_centers*2, min=c(-90, 30), max=c(-75, 40)), ncol=2, byrow=TRUE)
#' xy <- matrix(apply(centers, 1, function(x) rnorm(n*2, mean=x, sd=sd)), ncol=2, byrow=TRUE)
#' ss <- c(rep(1, n*n_centers/2), rep(5, n*n_centers/2))
#' n <- 1           # Points per cluster
#' n_centers <- 5   # Number of cluster centers
#' sd <- 1/2        # Standard deviation of each cluster
#' centers <- matrix(runif(n_centers*2, min=c(-90, 30), max=c(-75, 40)), ncol=2, byrow=TRUE)
#' xy <- rbind(xy,
#'            matrix(apply(centers, 1, function(x) rnorm(n*2, mean=x, sd=sd)), ncol=2, byrow=TRUE))
#' ss <- c(ss, rep(20, n*n_centers))
#' system.time(qt <- quadtree_ct(xy, ss, 20, min_in_bin=1))
#' xylim <- cbind(x=c(min(xy[,1]), max(xy[,1])), y=c(min(xy[,2]), max(xy[,2])))
#' png("qt_test.png")
#' plot(xylim, type="n", xlab="x", ylab="y", main="Quadtree w/ max SS = 20")
#' lines(qt, xylim, col="Gray")
#' points(qt, pch=16, cex=0.5, alpha=0.9)
#' points(xy, col=as.factor(ss), pch=16)
#' legend('bottomright', legend=c('1', '5', '20'), col=1:3, pch=16)
#' dev.off()
#' }
quadtree_ct <- function(xy, ss, target_ss, min_in_bin=5, rand = T) {
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

#' @title Plot Quadtree Points
#' @description Add quadtree points
#' 
#' @param q quadtree object
#' @param alpha transparency of points in plot
#' 
#' @return None
points.quadtree <- function(q, alpha = 0.1, ...) {
  points(q$lower, alpha, ...); points(q$upper, alpha, ...)
}

#' @title Plot Quadtree Leaf Points
#' @description Add quadtree leaf points
#' 
#' @param q quadtree object
#' @param alpha transparency of points in plot
#' 
#' @return None
points.quadtree.leaf <- function(q, alpha = 0.1, ...) {
  library(scales)
  points(q$value, col = alpha(q$id, alpha = alpha), ...)
}

#' @title Plot Lines for quadtree
#' @description Add lines for quadtree
#'
#' @param q quadtree object
#' @param xylim limits for xy
#' 
#' @return None
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

#' @title Plot Lines for quadtree leaves
#' @description Add lines for quadtree leaves
#'
#' @param q quadtree object
#' @param xylim limits for xy
#' 
#' @return None
lines.quadtree.leaf <- function(q, xylim, ...) {} # Nothing to do at leaves!


#' @title Get Cell Boundaries to Make Shapefiles
#' @description get cell boundaries to make shapefiles
#'
#' @param q quadtree or quadtree leaf object
#' @param xylim limits for xy
#' 
#' @return None
cell <- function(q, xylim, ...) {
  if (class(q)=="quadtree") f <- cell.quadtree else f <- cell.quadtree.leaf
  f(q, xylim, ...)
}

#' @title Get Cell Boundaries for Quadtree Object to Make Shapefiles
#' @description get cell boundaries for quadtree object to make shapefiles
#'
#' @param q quadtree object
#' @param xylim limits for xy
#' 
#' @return data frame of boundaries
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

#' @title Get Cell Boundaries for Quadtree Leaf Object to Make Shapefiles
#' @description get cell boundaries for quadtree leaf object to make shapefiles
#'
#' @param q quadtree leaf object
#' @param xylim limits for xy
#' 
#' @return dataframe of boundaries
cell.quadtree.leaf <- function(q, xylim) {
  data.frame(id = q$id,
             x = c(xylim[1,1], xylim[2,1], xylim[2,1], xylim[1,1], xylim[1,1]),
             y = c(xylim[1,2], xylim[1,2], xylim[2,2], xylim[2,2], xylim[1,2]))
}


#' @title Get IDs for each point
#' @description get ids for each point
#'
#' @param q quadtree or quadtree leaf object
#' 
#' @return dataframe of ids
id <- function(q, ...){
  d <- data.frame(id=NULL, x=NULL, y=NULL)
  if (class(q)=="quadtree") f <- id.quadtree else f <- id.quadtree.leaf
  f(q, xylim, ...)
}

#' @title Get IDs for each point
#' @description get ids for each point
#'
#' @param q quadtree object
#' 
#' @return dataframe of ids
id.quadtree <- function(q, ...) {
  rbind(id(q$lower), id(q$upper))
}

#' @title Get IDs for each point
#' @description get ids for each point
#'
#' @param q quadtree leaf object
#' 
#' @return dataframe of ids
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


#' @title Create Folds Based on NIDs
#' @description Randomly divide data into folds, based on NIDs (spatio-temporal 
#' holdout strategy)
#'
#' @param nids vector of NIDs for each data observation
#' @param yr vector of years of each observation
#' @param ss sample size of each data observation
#' @param n_folds number of oflds
#' 
#' @return returns a list of folds with each element containing which data points 
#' are in that fold
nid_folds <- function(nids, yr, ss, n_folds) {

  # calculate the total sample size by NID
  ss_by_nid <- tapply(ss, nids, sum)

  # randomly sort the NIDs
  ss_by_nid <- ss_by_nid[sample(1:length(ss_by_nid))]

  # calculate a running total sample size in the new sort order
  cumulative_ss <- Reduce(sum, ss_by_nid, accumulate = T)

  # identify five roughly equal folds based on the cumulative sample size
  target_fold_ss <- sum(ss) / n_folds
  brks <- sapply(1:5, function(x) which.min(abs(cumulative_ss - x * target_fold_ss)))
  folds <- cut(cumulative_ss, breaks = c(0, cumulative_ss[brks]), labels = F)

  # return in the proper format
  fold_list <- lapply(1:n_folds, function(x) {
    which(nids %in% names(ss_by_nid)[folds == x])
  })
  return(fold_list)

}