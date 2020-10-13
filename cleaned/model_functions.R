source('FILEPATH.R')
package_list <- c('seegSDM', 'gbm', 'rgdal', 'dismo', 'maptools')
load_R_packages(package_list)

## get.vals.yrs(): function for extracting year indexed covariate values and appending new column to input data
## returns: data.frame of data with extracted covariate values appended
## params:
## l.lim 	-> the oldest year of available raster data
## u.lim 	-> the most recent year of available raster data
## data	  -> data.frame w/ 'year' colm and 'longitude' 'latitude' colms in col 2:3
## brick 	-> raster brick with time indexed layers (lyr name contains year and month)
## beg   	-> first index of year in brick lyr name
## end   	-> last index of month in brick lyr name
## name  	-> name of extracted values column

get.vals.yrs <- function(l.lim, u.lim, data, brick, yrbeg, yrend, name){
  yrs <- as.numeric(levels(data.frame(table(data$year))[,1]))
  df.out <- data.frame()
  old.dat <- data[data$year<=l.lim,]
  if (nrow(old.dat) != 0){
    old.dat <- cbind(old.dat,
                     data.frame(extract(brick[[which(substr(names(brick), yrbeg, yrend) %in% l.lim)]],
                                        old.dat[,c('long', 'lat')])))
    names(old.dat) <- append(names(data), name)
    df.out <- rbind(df.out, old.dat)
  }
  new.dat <- data[data$year>=u.lim,]
  if (nrow(new.dat) != 0){
    new.dat <- cbind(new.dat,
                     data.frame(extract(brick[[which(substr(names(brick), yrbeg, yrend) %in% u.lim)]],
                                        new.dat[,c('long', 'lat')])))
    names(new.dat) <- append(names(data), name)
    df.out <- rbind(df.out, new.dat)
  }
  yrs <- yrs[yrs>l.lim & yrs<u.lim]
  if (length(yrs) == 0){
    return(df.out)
  } else {
    for (i in 1:length(yrs)){
      mid <- data[data$year==yrs[i],]
      mid <- cbind(mid,
                   data.frame(extract(brick[[which(substr(names(brick), yrbeg, yrend) %in% yrs[i])]],
                                      mid[,c('long', 'lat')])))
      names(mid) <- append(names(data), name)
      df.out <- rbind(df.out, mid)
    }
    return(df.out)
  }
}

## get.vals.yrs_mos(): function for extracting year_month indexed covariate values and appending new column to input data
## returns: data.frame of data with extracted covariate values appended
## params:
## l.lim 	-> the oldest year_mo combo of available raster data
## u.lim 	-> the most recent year_mo combo of available raster data
## data	  -> data.frame w/ 'year' colm and 'longitude' 'latitude' colms in col 2:3
## brick 	-> raster brick with time indexed layers (lyr name contains year and month)
## beg   	-> first index of year in brick lyr name
## end   	-> last index of month in brick lyr name
## name  	-> name of extracted values column

get.vals.yrs_mos <- function(l.lim, u.lim, data, brick, beg, end, name){
  yrs <- as.numeric(levels(data.frame(table(data$year))[,1]))
  mos <- as.numeric(levels(data.frame(table(data$month))[,1]))
  # combine yrs and mos to match indices on covariate names...
  # ...for entire range (make list of all year/month combinations)
  yrs_mos <- c()
  for (i in 1:(length(yrs))){
    for (j in 1:(length(mos))){
      if (mos[j] <= 9){
        yrs_mos <- append(yrs_mos, paste0(yrs[i],'_0',mos[j]))
      } else {
        yrs_mos <- append(yrs_mos, paste0(yrs[i],'_',mos[j]))
      }
    }
  }
  # ...for given data time-stamps (combine year/month from data)
  yrs_mos.data <- c()
  for (i in 1:(nrow(data))){
    if (data$month[i] <= 9){
      yrs_mos.data <- append(yrs_mos.data, paste0(data$year[i],'_0',data$month[i]))
    } else {
      yrs_mos.data <- append(yrs_mos.data, paste0(data$year[i],'_',data$month[i]))
    }
  }
  ## make output data.frame object
  df.out <- data.frame()
  # extract values...
  # ...prior to existing covariate date ranges
  old.dat <- data[yrs_mos.data < l.lim,]
  if (nrow(old.dat) > 0){
    yrs_mos.data_old <- paste0(substr(l.lim, 1, 4), "_", substr(yrs_mos.data[yrs_mos.data < l.lim],6,7))
    for (i in 1:length(yrs_mos.data_old)){
      if (substr(yrs_mos.data_old[i], 6, 7) < substr(l.lim, 6, 7)) {
        yrs_mos.data_old[i] <- gsub(substr(yrs_mos.data_old[i], 5, 7), substr(l.lim, 5, 7), yrs_mos.data_old[i])
      }
      old.ext <- cbind(old.dat[i,],
                       data.frame(extract(brick[[which(substr(names(brick), beg, end) %in% yrs_mos.data_old[i])]],
                                          old.dat[i, c('long', 'lat')])))
      names(old.ext) <- append(names(data), name)
      df.out <- rbind(df.out, old.ext)
    }
  }
  # ...following existing covariate date ranges
  new.dat <- data[yrs_mos.data >= u.lim,]
  if (nrow(new.dat) > 0){
    yrs_mos.data_new <- paste0(substr(u.lim, 1, 4), "_", substr(yrs_mos.data[yrs_mos.data >= u.lim],6,7))
    for (i in 1:length(yrs_mos.data_new)){
      if (substr(yrs_mos.data_new[i], 6, 7) > substr(u.lim, 6, 7)) {
        yrs_mos.data_new[i] <- gsub(substr(yrs_mos.data_new[i], 5, 7), substr(u.lim, 5, 7), yrs_mos.data_new[i])
      }
      new.ext <- cbind(new.dat[i,],
                       data.frame(extract(brick[[which(substr(names(brick), beg, end) %in% yrs_mos.data_new[i])]],
                                          new.dat[i, c('long', 'lat')])))
      names(new.ext) <- append(names(data), name)
      df.out <- rbind(df.out, new.ext)
    }
  }
  # ...from covariates corresponding to data time-stamps
  yrs_mos <- yrs_mos[yrs_mos>=l.lim & yrs_mos<u.lim]
  mid.idx <- which(yrs_mos.data %in% yrs_mos)
  yrs_mos.data <- yrs_mos.data[mid.idx]
  mid <- data[mid.idx,]
  if (nrow(mid) == 0){
    return(df.out)
  } else {
    for (i in 1:length(yrs_mos.data)){
      mid.ext <- cbind(mid[i,],
                       data.frame(extract(brick[[which(substr(names(brick), beg, end) %in% yrs_mos.data[i])]],
                                          mid[i, c('long', 'lat')])))
      names(mid.ext) <- append(names(data), name)
      df.out <- rbind(df.out, mid.ext)
    }
    return(df.out)
  }
}

## run_optimizerPy(): function for running brt hyperparameter optimization script in python
## returns: NULL, launches python script
## params:
## python           -> directory to pyton binary to run script with
## funcs.file_path  -> directory to brt hyperparameter optimization script
## funcs.file       -> name of brt hyperparameter optimization script
## bounds.file_path -> directory to search space boundary file
## bounds.file      -> name of search space boundary file
## data.file_path   -> directory to csv of to model environmental suitability of
## data.loc         -> name of csv of model environmental suitability of
## optimizer        -> string specifying choice of non-parametric model to optimize with: brt, gp, rf
##                      -- brt = Boosted Regression Tree, gp = Gaussian Process, rf = Random Forest
## learner          -> string specifying choice of BRT to optimize
## cv_folds         -> number of cross-validation folds to confirm parameter selection
## n_calls          -> number of iterations of optimization to perform
## jobnum           -> index of 
## col_start        -> index of first column of environmental data present in data.loc
##                      -- assumes col_start through final column in data.loc are intended for modelling, ie - contain environmental covariate data

run_optimizerPy <- function(python,
                            funcs.file_path,
                            funcs.file,
                            bounds.file_path,
                            bounds.file,
                            data.file_path,
                            data.loc,
                            optimizer,
                            learner,
                            cv_folds,
                            n_calls,
                            jobnum,
                            col_start)
{
  system(
    paste0(python, ' ',
           funcs.file_path, funcs.file, ' ',
           bounds.file_path, bounds.file, ' ',
           data.file_path, ' ',
           data.file_path, data.loc, ' ',
           optimizer, ' ',
           learner, ' ',
           cv_folds, ' ',
           n_calls, ' ',
           jobnum, ' ',
           col_start))
}

## calc_tbin_prob(): function for computing the probability of time bins as they occurr in the base dataset
## returns: tuple of year frequencies and month frequencies
## params:
## dat_orig -> data.frame containing columns: year_start, year_end, month_start, month_end 

calc_tbin_prob <- function(dat_orig){
  year <- c()
  month <- c()
  
  # sample from w/in time bins: mo_start, yr_start -> mo_end, yr_end
  for (i in 1:nrow(dat_orig)){
    yr_start <- dat_orig$year_start[i]
    yr_end <- dat_orig$year_end[i]
    mo_start <- dat_orig$month_start[i]
    if (is.na(mo_start)) mo_start <- 1
    mo_end <- dat_orig$month_end[i]
    if (is.na(mo_end)) mo_end <- 1
    
    # get no. of months occur in each year
    yr_list_freq <- c()
    for (t in yr_start:yr_end){
      if (t == yr_start){
        tt <- mo_start
        while (!tt%%12) tt <- tt + 1
      } else if (t == yr_end) {
        tt <- 1
        while (tt!=mo_end) tt <- tt + 1
      } else {
        tt <- 12
      }
      yr_list_freq <- append(yr_list_freq, tt)
    }
    
    # make list of years with necessary number of copies - sample likelihood proportional to what exists in data
    yr_list <- c()
    yrs <- yr_start:yr_end
    for (yr in 1:length(yrs)) yr_list <- append(yr_list, rep(yrs[yr], yr_list_freq[yr]))
    
    # make lists of months and years including sampled date ranges
    if (!is.na(yr_start) && !is.na(yr_end)){
      year[i] <- safe_sample(yr_list)
      if (!is.na(mo_start) && !is.na(mo_end)){
        if (yr_start == yr_end){
          month[i] <- safe_sample(mo_start:mo_end)
        } else if (year[i] == yr_start){
          month[i] <- safe_sample(mo_start:12)
        } else if (year[i] == yr_end) {
          month[i] <- safe_sample(1:mo_end)
        } else {
          month[i] <- safe_sample(1:12)
        }
      } else {
        month[i] <- NA
      }
    } else {
      year[i] <- NA
    }
    if (is.na(mo_start) || is.na(mo_end)){
      month[i] <- NA
    }
  }
  # fit distribution to month data, get probability of a given bin (representative of that month)
  mo_dens <- density(na.omit(month))
  mo_freqs <- c()
  for (mo in 1:12){
    a <- mo -1
    if (a < min(mo_dens$x)) a <- min(mo_dens$x)
    mo_freqs <- append(integrate.xy(mo_dens$x, mo_dens$y, a, mo), mo_freqs)
  }
  # again ~ years
  yr_dens <- density(na.omit(year))
  yr_range <- min(year):max(year)
  yr_freqs <- c()
  for (yr in yr_range){
    a <- yr -1
    if (a < min(yr_dens$x)) a <- min(yr_dens$x)
    yr_freqs <- append(integrate.xy(yr_dens$x, yr_dens$y, a, yr), yr_freqs)
  }
  yr_freqs <- abs(yr_freqs)
  return(list(yr_freqs=yr_freqs, mo_freqs=mo_freqs))
}

## point.data(): function for creating time-stamped point data
## returns: data.frame of point occurrences with estimated time-bin columns for missing values
## params:
## dat_orig  -> base data.frame to subset occurrences from with columsn names: shape_type, long, lat, PA
## yr_freqs  -> list from calc_tbin_prob() output
## mo_freqs  -> list from calc_tbin_prob() output

point.data <- function(dat_orig, yr_freqs, mo_freqs){
  occ <- cbind(dat_orig[, c('shape_type', 'long', 'lat')], cbind(year, month))
  occ <- occ[occ$shape_type == 'point',]
  occ <- occ[,-which(names(occ) %in% 'shape_type')]
  names(occ) <- c('long','lat','year','month')
  
  # impute time data based on frequency of those months in the dataset
  for (i in 1:nrow(occ)){
    if (is.na(occ$year[[i]])) occ$year[[i]] <- sample(yr_range, 1, prob=yr_freqs)
    if (is.na(occ$month[[i]])) occ$month[[i]] <- sample(1:12, 1, prob=mo_freqs)
  }
  occ <- cbind(PA = rep(1, nrow(occ)), occ)
  return(occ)
}

## bg.data(): function for creating time-stamped background data
## returns: data.frame of background points with estimated time-bin columns for missing values
## params:
## dat_orig  -> base data.frame to get number of rows from
## buffer    -> spatial shape to sample background points from
## bias  -> bias raster defining regions to avoid sampling
## yr_freqs  -> list from calc_tbin_prob() output
## mo_freqs  -> list from calc_tbin_prob() output

bg.data <- function(dat_orig, buffer, bias, yr_freqs, mo_freqs) {
  bg_smpsz <- nrow(dat_orig)
  if (is.na(bias)){
     bg <- data.frame(spsample(x=buffer, n=bg_smpsz, type='random')) # w/out pop. biasing
  } else {
    bg <- data.frame(bgSample(n = bg_smpsz,
                              raster = mask(bias, buffer, updatevalue=0),
                              replace=F, prob=T))
  }
  
  # impute time-stamps based on existing data
  bg_years <- c()
  for (i in 1:bg_smpsz){
    bg_years[[i]] <- sample(yr_range, 1, prob=yr_freqs)
  }
  bg <- cbind(bg, bg_years)
  bg_months <- c()
  for (i in 1:bg_smpsz){
    bg_months[[i]] <- sample(1:12, 1, prob=mo_freqs)
  }
  bg <- cbind(bg, data.frame(bg_months))
  names(bg) <- c('long', 'lat', 'year', 'month')
  
  # combine occurence and background data to extract covariate values for
  dat <- cbind(PA = rep(0, nrow(bg)), bg)
  return(dat)
}

## bg.data(): function for creating time-stamped point data from buffer regions
## returns: data.frame of occurrences with estimated time-bin columns for missing values
## params:
## dat_orig     -> base data.frame to get number of rows from
## buffer_data  -> data.frame of buffer data rows
## bias  -> bias raster defining regions to avoid sampling
## yr_freqs  -> list from calc_tbin_prob() output
## mo_freqs     -> list from calc_tbin_prob() output

buffer.data <- function(dat_orig, buffer_data, bias, yr_freqs, mo_freqs) {
  # draw from buffer regions - 1 draw / buffer
  smp_buffs <- data.frame()
  for (i in 1:nrow(buffer_data)) {
    # reference individual buffer
    buffer <- buffer_data[i,]
    # draw sample from w/in buffer
    pts <- data.frame(bgSample(n = 1,
                               raster = mask(bias, buffer, updatevalue=0),
                               replace=F, prob=T))
    names(pts) <- c('long', 'lat')
    smp_buffs <- rbind(smp_buffs, pts)
  }
  
  # subset out buffer data
  occ_buffs <- cbind(dat_orig, cbind(year, month))
  occ_buffs <- occ_buffs[occ_buffs$poly_type == 'buffer',]
  # match buffer data to already triaged buffer data from 'make_bg_buffers.R'
  occ_buffs <- occ_buffs[which(buffer_data$occ_id %in% occ_buffs$occ_id),]
  
  # impute month data based on frequency of those months in the dataset
  for (i in 1:nrow(occ_buffs)){
    if (is.na(occ_buffs$year[[i]])) occ_buffs$year[[i]] <- sample(yr_range, 1, prob=yr_freqs)
    if (is.na(occ_buffs$month[[i]])) occ_buffs$month[[i]] <- sample(1:12, 1, prob=mo_freqs)
  }
  
  buf.dat <- cbind(PA = rep(1, nrow(smp_buffs)), cbind(smp_buffs, occ_buffs[,c('year','month')]))
  return(buf.dat)
}

## bg.data(): function for creating time-stamped point data from spatial shape regions
## returns: data.frame of point occurrences with estimated time-bin columns for missing values
## params:
## dat_orig      -> base data.frame to get number of rows from
## polygon_data  -> data.frame of polygon data rows
## bias  -> bias raster defining regions to avoid sampling
## yr_freqs  -> list from calc_tbin_prob() output
## mo_freqs      -> list from calc_tbin_prob() output

polygon.data <- function(dat_orig, polygon_data, bias, yr_freqs, mo_freqs) {
  # subset out polygons
  occ_polys <- cbind(dat_orig, cbind(year, month))
  occ_polys <- occ_polys[occ_polys$shape_type == 'polygon',]
  occ_polys <- occ_polys[occ_polys$poly_type != 'buffer',]
  occ_polys <- occ_polys[occ_polys$poly_type != "",]
  occ_polys <- occ_polys[!is.na(occ_polys$poly_reference),]
  # match buffer data to already triaged buffer data from 'FILEPATH.R'
  occ_polys <- occ_polys[which(polygon_data$occ_id %in% occ_polys$occ_id),]
  
  # impute month data based on frequency of those months in the dataset
  for (i in 1:nrow(polygon_data)){
    if (is.na(occ_polys$year[[i]])) occ_polys$year[[i]] <- sample(yr_range, 1, prob=yr_freqs)
    if (is.na(occ_polys$month[[i]])) occ_polys$month[[i]] <- sample(1:12, 1, prob=mo_freqs)
  }
  
  # sample polygons for occurrence points
  smp.occ <- data.frame()
  for (i in 1:nrow(polygon_data)){
    # reference individual polygon
    polygon <- polygon_data[i,]
    # draw from selected polygon
    pts <- data.frame(bgSample(n = 1,
                               raster = mask(bias, polygon, updatevalue=0),
                               replace=F, prob=T))
    names(pts) <- c('long', 'lat')
    smp.occ <- rbind(smp.occ, pts)
  }
  
  smp.dat <- cbind(PA = rep(1, nrow(smp.occ)), cbind(smp.occ, occ_polys[1:nrow(polygon_data),c('year','month')]))
  return(smp.dat)
}

## dat.all(): function for taking time-stamped point and extracting temporally varying environmental data for each point
##             -- performs extraction from specified file system
## returns: data.frame of point occurrences with associated environmental covariate data
## params:
## smp.dat        -> data.frame of polygon sampled data with time stamps
## buf.dat        -> data.frame of buffer sampled data with time stamps
## pt.dat         -> data.frame of point data with time stamps
## bg.dat         -> data.frame of background point data with time stamps
## cov_name_list  -> data.frame with columns: cov_name, measure, measure_yrly, monthly, beg, end, l_lim, u_lim, units, static
## cov_dir        -> string specifying the location of environmental covariate data
## stack_out_dir  -> string specifying the location of environmental covariate data, independent of cov_dir

# extract covariate data from raster stack by long/lat and year/(month <- if exists)
dat.all <- function(smp.dat, buf.dat, pt.dat, bg.dat, cov_name_list, cov_dir, stack_out_dir){
  dat_all <- rbind(smp.dat, buf.dat, bg.dat, pt.dat)
  
  for (j in 1:nrow(cov_name_list)){
    # get requisite data from master list
    measure <- cov_name_list$measure[j]
    cov_name <- cov_name_list$cov_name[j]
    monthly <- as.logical(cov_name_list$monthly[j])
    beg <- cov_name_list$beg[j]
    end <- cov_name_list$end[j]
    u_lim <- toString(cov_name_list$u_lim[j])
    l_lim <- toString(cov_name_list$l_lim[j])
    static <- as.logical(cov_name_list$static[j])
    if (static) {
      print(cov_name)
      f_dir <- paste0(cov_name, '/', measure, '/synoptic/', cov_name, '_', measure, '_synoptic.tif')
      cov <- raster(paste0(cov_dir, '/', f_dir))
      names(cov) <- cov_name
      dat_all <- cbind(dat_all, extract(cov, dat_all[,c('long', 'lat')]))
      names(dat_all)[ncol(dat_all)] <- cov_name
    } else {
      print(cov_name)
      if (monthly){
        t_bin <- 'monthly'
      } else {
        t_bin <- 'yearly'
      }
      # read in stack
      cov.flist <- list.files(path=stack_out_dir, pattern=paste0(cov_name, '_', t_bin), full.names=T)
      cov.flist <- cov.flist[tools::file_ext(cov.flist)=='grd']
      cov_stack <- do.call(stack, lapply(cov.flist, brick))
      # extract and append covariate data as column to dat_all w/ colm name 'cov_name'
      if (monthly){
        dat_all <- get.vals.yrs_mos(l_lim, u_lim, dat_all, cov_stack, beg, end, cov_name)
      } else {
        dat_all <- get.vals.yrs(l_lim, u_lim, dat_all, cov_stack, beg, end, cov_name)
      }
    }
  }
  
  # omit missing values
  dat_all <- na.omit(dat_all)
  return(dat_all)
}

## run_brt_model(): function for taking time-stamped point and extracting temporally varying environmental data for each point
##             -- performs extraction from specified file system
## returns: brt model predictions and fit statistics. if final=T, outputs additional information
## params:
## data       -> data.frame of model input data with environmental data columns matching cov_names
## par        -> list of BRT hyperparameters (provided by optimizers.py)
## covs       -> environmental covariate raster brick/stack for specified (contemporary) values to make environmental suitability predictions
## cov_names  -> names given in dat.all() to columns of environmental covariate data and found in "data" parameter
## final      -> boolean; if T, function will output model object

run_brt_model <- function(data, par, covs, cov_names, final=T) {
  # define weighting function (bernoulli)
  wt.func <- function(PA) ifelse(PA == 1, 1, sum(PA) / sum(1 - PA))
  wt <- wt.func(data$PA)
  # run model
  model <- gbm(data$PA~.,
               distribution = 'bernoulli',
               data = data[,cov_names],
               n.trees = par$n.trees,
               interaction.depth = par$interaction.depth,
               shrinkage = par$shrinkage,
               n.minobsinnode = par$n.minobsinnode,
               weights = wt,
               bag.fraction = 0.8,
               train.fraction = 0.8,
               verbose = F)
  # get fit statistics
  gbm.pred <- predict.gbm(model,
                          dat_all[,cov_names],
                          n.trees = model$n.trees,
                          type = 'response')
  stats <- calcStats(data.frame(dat_all$PA, gbm.pred))
  # only predict if not train/test split to validate par
  if (final) {
    # effect curves
    effects <- lapply(1:length(cov_names),
                      function(i) plot(model, i, return.grid = TRUE))
    # get relative influence
    relinf <- summary(model, plotit = FALSE)
    # create prediction raster based on most comtemporary covariate values
    names(covs) <- cov_names
    pred.raster <- predict(covs, model, type = 'response', n.trees = model$n.trees)
    # get coordinates
    coords <- dat_all[,c('long', 'lat')]
    # format like seegSDM::runBRT() result
    model.list <- list(model = model,
                       effects = effects,
                       relinf = relinf,
                       pred = pred.raster,
                       coords = coords)
    return(list(stats=stats,
                model.list=model.list,
                pred.raster=pred.raster))
  } else {
    return(list(stats=stats))
  }
}
