# function file

## Make time stamp in standardized format.
make_time_stamp <- function(time_stamp) {
  
  run_date <- gsub("-","_",Sys.time())
  run_date <- gsub(":","_",run_date)
  run_date <- gsub(" ","_",run_date)
  
  if(time_stamp==FALSE) run_date <- 'scratch'
  
  return(run_date)
  
}

firstDay <- function (year, month) {
  # given a year and month, return a Date object of the first day of that month
  date_string <- paste(year,
                       month,
                       '01',
                       sep = '-')
  
  date <- as.Date (date_string)
  
  return (date)
  
}

lastDay <- function (year, month) {
  # given a year and month, return a Aate object of the last day of that month
  next_month <- ifelse(month == 12,
                       1,
                       month + 1)
  
  next_year <- ifelse(month == 12,
                      year + 1,
                      year)
  
  next_date_string <- paste(next_year,
                            next_month,
                            '01',
                            sep = '-')
  next_date <- as.Date(next_date_string)
  date <- next_date - 1
  
  return (date)
}

sentenceCase <- function (text) {
  # given a vector of text strings `text`, convert to sentence case
  
  # convert all to lower case
  text <- tolower(text)
  
  # split at spaces
  text_list <- strsplit(text,
                        ' ')
  
  text_list <- lapply(text_list,
                      function(x) {
                        x[1] <- paste(toupper(substring(x[1], 1, 1)),
                                      substring(x[1], 2),
                                      sep = "")
                        x <- paste(x, collapse = ' ')
                        return(x)
                      })
  
  text_vector <- unlist(text_list)
  
  return (text_vector)
  
}

# define functions
firstTwo <- function (text) {
  # given a vector of text strings `text` subset each to only the first two
  # words (bits separated by spaces) and return this as a vector.
  text <- as.character(text)
  text_list <- strsplit(text, ' ')
  text_list <- lapply(text_list, '[', 1:2)
  text_list <- lapply(text_list, paste, collapse = ' ')
  text_vector <- unlist(text_list)
  return (text_vector)
}

rasterizeSpecies <- function(species,
                             shape,
                             raster,
                             buffer = NULL,
                             folder = 'FILEPATH') {
  
  # first buffer, then rasterize IUCN range map for species
  shp <- shape[shape$BINOMIAL == species, ]
  
  if (!is.null(buffer)) {
    # convert buffer from kilometres to decimal degrees, assuming at the equator
    buffer <- buffer / 111.32
    
    # buffer by this amount
    shp <- gBuffer(shp,
                   width = buffer)
  }
  
  # rasterize the shapefile 
  tmp <- rasterize(shp,
                   raster,
                   field = 1,
                   background = 0,
                   fun = 'first')
  
  writeRaster(tmp,
              filename = paste0('FILEPATH',
                                folder,
                                '/',
                                gsub(' ', '_', species)),
              format = 'GTiff',
              overwrite = TRUE)
  
  rm(tmp)
  
  return (NULL)
}

tidySpecies <- function (filename, template) {
  # load a raster if it contains any of the species' range,
  # mask and resave it, else delete it
  tmp <- raster(filename)
  if (!is.na(maxValue(tmp)) && maxValue(tmp) == 1) {
    tmp <- mask(tmp,
                template)
    
    writeRaster(tmp,
                file = filename,
                overwrite = TRUE)
    
  } else {
    
    rm(tmp)
    
    file.remove(filename)
    
  }
  
  return (NULL)
  
}

subsamplePolys <- function (data, ...) {
  # given a presence-background dataset, with multiple rows for some of the
  # occurrence records, subset it so that there's only one randomly selected
  # point from each polygon and then take a bootstrap it using `subsample`.
  # Dots argument is passed to subsample.
  
  # index for background records (outbreak id = 0)
  bg_idx <- data$outbreak_id == 0
  
  # subset to get occurrence section only
  occ <- data[!bg_idx, ]
  
  # get the different outbreaks
  u <- unique(occ$outbreak_id)
  
  # loop through, picking an index for each based on the number available
  occ_idx <- sapply(u,
                    function (id, occ) {
                      idx <- which(occ$outbreak_id == id)
                      sample(idx, 1)
                    },
                    occ)
  
  # get the subsetted dataset
  dat <- rbind(occ[occ_idx, ],
               data[bg_idx, ])
  
  # randomly subsample the dataset
  ans <- subsample(dat,
                   n = nrow(dat),
                   ...)
  
  # remove the outbreak ID column
  ans <- ans[, -which(names(ans) == 'outbreak_id')]
  
  return (ans)
}

# change the polygon IDs of an SPDF so it can be rbinded to something else
makeUniform <- function (SPDF) {
  pref <- substitute(SPDF)  #just putting the file name in front.
  newSPDF <- spChFIDs(SPDF,
                      as.character(paste(pref,
                                         rownames(as(SPDF,
                                                     "data.frame")),
                                         sep = "_")))
  return (newSPDF)
}

summarizeStats <- function (path) {
  
  # load validation stats
  stats <- read.csv(paste0(path, 'FILEPATH.csv'),
                    row.names = 1)
  
  auc  <- c(as.character(round(mean(stats$auc,
                                    na.rm = TRUE),
                               2)),
            as.character(round(sd(stats$auc,
                                  na.rm = TRUE),
                               2)))
  
  # load relative influence stats
  relinf <- read.csv(paste0(path,
                            'FILEPATH.csv'),
                     stringsAsFactors = FALSE)
  
  ans <- c(auc_mean = auc[1],
           auc_sd = auc[2],
           cov1 = relinf[1, 1],
           relinf1 = as.character(round(relinf[1, 2],
                                        1)),
           cov2 = relinf[2, 1],
           relinf2 = as.character(round(relinf[2, 2],
                                        1)),
           cov3 = relinf[3, 1],
           relinf3 = as.character(round(relinf[3, 2],
                                        1)),
           cov4 = relinf[4, 1],
           relinf4 = as.character(round(relinf[4, 2],
                                        1)),
           cov5 = relinf[5, 1],
           relinf5 = as.character(round(relinf[5, 2],
                                        1)))
  
  return (ans)
  
}

thresholdRisk <- function (risk_raster,
                           occ,
                           proportion = 1) {
  
  # given a raster `risk_raster` giving risk on the (0,1] level and a 2d
  # dataframe `occ` with columns named 'lat' and 'long' giving the latitudes
  # and longitudes of known occurrence records,
  # find the threshold value so that `proportion` fraction of the records
  # fall in areas classified as 'at risk' and return the thresholded map
  
  # extract risk values for the occurrence data
  occ_risk <- extract(risk_raster[[1]],
                      occ[, c('long', 'lat')])
  
  # remove any missing data
  occ_risk <- na.omit(occ_risk)
  
  # get the relevant quantile
  thresh <- quantile(occ_risk,
                     1 - proportion,
                     na.rm = TRUE)
  
  # classify the raster
  at_risk_raster <- risk_raster > thresh
  
  # return this
  return (at_risk_raster)
  
}

# function for extracting year indexed covariate values and appending new column to input data
get.vals.yrs <- function(l.lim, 		# the oldest year index of available raster data
                         u.lim, 		# the most recent year index of available raster data
                         data,  		# data.frame w/ 'year' colm and longitude latitude colms named 'long' and 'lat'
                         brick, 		# raster brick with time indexed layers (lyr name contains year)
                         yrbeg, 		# first index of year in brick lyr name
                         yrend,  		# last index of year in brick lyr name
                         name   		# name of extracted values column
){
  yrs <- as.numeric(levels(data.frame(table(data$year))[,1]))
  df.out <- data.frame()
  old.dat <- data[data$year<=l.lim,]
  if (nrow(old.dat) != 0){
    old.dat <- cbind(old.dat,
                     data.frame(raster::extract(brick[[which(substr(names(brick), yrbeg, yrend) %in% l.lim)]],
                                        old.dat[,c('long', 'lat')])))
    names(old.dat) <- append(names(data), name)
    df.out <- rbind(df.out, old.dat)
  }
  new.dat <- data[data$year>=u.lim,]
  if (nrow(new.dat) != 0){
    new.dat <- cbind(new.dat,
                     data.frame(raster::extract(brick[[which(substr(names(brick), yrbeg, yrend) %in% u.lim)]],
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
                   data.frame(raster::extract(brick[[which(substr(names(brick), yrbeg, yrend) %in% yrs[i])]],
                                      mid[,c('long', 'lat')])))
      names(mid) <- append(names(data), name)
      df.out <- rbind(df.out, mid)
    }
    return(df.out)
  }
}

# function for extracting year_month indexed covariate values and appending new column to input data
get.vals.yrs_mos <- function(l.lim, 	# the oldest year_mo combo of available raster data
                             u.lim, 	# the most recent year_mo combo of available raster data
                             data,  	# data.frame w/ 'year' colm and 'longitude' 'latitude' colms in col 2:3
                             brick, 	# raster brick with time indexed layers (lyr name contains year and month)
                             beg,   	# first index of year in brick lyr name
                             end,   	# last index of month in brick lyr name
                             name   	# name of extracted values column
){
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
    month <- data$month[i]
    if (month <= 9) month <- paste0(0, data$month[i])
    yrs_mos.data <- append(yrs_mos.data, paste0(data$year[i],'_',month))
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
                       data.frame(raster::extract(brick[[which(substr(names(brick), beg, end) %in% yrs_mos.data_old[i])]],
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
                       data.frame(raster::extract(brick[[which(substr(names(brick), beg, end) %in% yrs_mos.data_new[i])]],
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
                       data.frame(raster::extract(brick[[which(substr(names(brick), beg, end) %in% yrs_mos.data[i])]],
                                          mid[i, c('long', 'lat')])))
      names(mid.ext) <- append(names(data), name)
      df.out <- rbind(df.out, mid.ext)
    }
    return(df.out)
  }
}

# function for submitting optimizer jobs
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

path_converter <- function(winpath, disc){
  if (disc == "ADDRESS"){
    path <- gsub("FILEAPTH", "FILEPATH", gsub("\\\\", "/", winpath))
  }
  return(path)
}

# copied from sfsmisc
integrate.xy <- function(x,fx, a,b, use.spline = TRUE, xtol = 2e-8)
{
  if(is.list(x)) {
    fx <- x$y; x <- x$x
    if(length(x) == 0)
      stop("list 'x' has no valid $x component")
  }
  if((n <- length(x)) != length(fx))
    stop("'fx' must have same length as 'x'")
  
  if(is.unsorted(x)) { i <- sort.list(x); x <- x[i]; fx <- fx[i] }
  if(any(i <- duplicated(x))) {
    n <- length(x <- x[!i])
    ## we might have to check that the same fx[] are duplicated
    ## otherwise either give an error or take the mean() of those...
    fx <- fx[!i]
  }
  if(any(diff(x) == 0))
    stop("bug in 'duplicated()' killed me: have still multiple x[]!")
  
  if(missing(a)) a <- x[1]
  else if(any(a < x[1])) stop("'a' must NOT be smaller than min(x)")
  if(missing(b)) b <- x[n]
  else if(any(b > x[n])) stop("'b' must NOT be larger  than max(x)")
  if(length(a) != 1 && length(b) != 1 && length(a) != length(b))
    stop("'a' and 'b' must have length 1 or same length !")
  else {
    k <- max(length(a),length(b))
    if(any(b < a))    stop("'b' must be elementwise >= 'a'")
  }
  
  if(use.spline) {
    xy <- spline(x,fx, n = max(1024, 3*n))
    if(xy$x[length(xy$x)] < x[n]) {
      if(TRUE) cat("working around spline(.) BUG")
      xy$x <- c(xy$x,  x[n])
      xy$y <- c(xy$y, fx[n])
    }
    x <- xy$x; fx <- xy$y
    n <- length(x)
  }
  
  ab <- unique(c(a,b))
  BB <- abs(outer(x,ab,"-")) < (xtol * max(b - a))
  if(any(j <- 0 == colSums(BB))) { # the j-th element(s) of ab are not in x[]
    y <- approx(x,fx, xout = ab[j])$y
    x <- c(ab[j],x)
    i <- sort.list(x)
    x <- x[i];  fx <- c(y,fx)[i];  n <- length(x)
  }

  dig0 <- floor(-log10(xtol)) 
  f.match <- function(x,table,dig) match(signif(x,dig), signif(table,dig))

  d <- dig0; while(anyNA(ai <- f.match(a,x, d))) d <- d - 1/8 ; ai <- rep_len(ai, k)
  d <- dig0; while(anyNA(bi <- f.match(b,x, d))) d <- d - 1/8 ; bi <- rep_len(bi, k)
  dfx <- fx[-c(1,n)] * diff(x,lag = 2)
  r <- numeric(k)
  for (i in 1:k) {
    a <- ai[i];  b <- bi[i]
    r[i] <- (x[a+1] - x[a])*fx[a] + (x[b] - x[b-1])*fx[b] +
      sum(dfx[seq(a, length = max(0,b-a-1))])
  }
  r/2
}

# sample while checking for length of 1
safe_sample <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } else {
    return(sample(x,1))
  }
}

get_top_months <- function(data, nmonths) {
  exact_dates <- data.frame(data[data$obs_date_high == data$obs_date_low, "obs_date_low"], stringsAsFactors = FALSE)
  names(exact_dates) <- "date"
  for (i in 1:nrow(exact_dates)) {
    exact_dates[i, "year"] <- as.character(format(as.Date(exact_dates$date[i]), "%Y"))
    exact_dates[i, "month"] <- as.character(format(as.Date(exact_dates$date[i]), "%m"))
  }
  
  #just get month and year
  exact_dates <- exact_dates[,-1]
  
  combination_counts <- melt(table(exact_dates))
  top_months <- head(combination_counts[order(-combination_counts$value), c("year", "month")], nmonths)
  return(top_months)
}

get_map_for_month <- function(month, year, out_dir) {
  source('FILEPATH.R')
  source('FILEPATH.R')
  source('FILEPATH.R')
  library(reshape2)
  
  configs <- read.csv(paste0(out_dir, '/FILEPATH.csv'), stringsAsFactors = FALSE)
  in_dir <- configs$in_dir
  njobs <- configs$njobs
  project <- configs$project
  smry_slots <- configs$smry_slots
  if (month <= 9) month <- paste0(0, month)

  #make pred stack
  if (!file.exists(paste0(out_dir, "/FILEPATH", year, "_", month, ".grd"))) {
    cov_pred.run <- paste0(in_dir, '/FILEPATH.R')
    qsub(paste0("make_preds_", year, "_", month),
         cov_pred.run,
         pass=c(out_dir, out_dir, 'TRUE', year, month),
         log=T,
         slots=2,
         submit=T,
         proj=project,
         user='ahcastle')
    check_loc_results(month, paste0(out_dir, '/FILEPATH'), prefix=paste0("pred_stack_", year, '_'), postfix=".grd")
  }
  
  #make monthly pred
  
  mo_pred.run <- paste0(in_dir, "/FILEPATH.R")
  pred_mo_dir <- paste0(out_dir, '/FILEPATH', year, '_', month)
  
  if (!file.exists(pred_mo_dir)) dir.create(pred_mo_dir)
  for (jobnum in 1:njobs) {
    if (!file.exists(paste0(pred_mo_dir, '/FILEPATH', jobnum, '.tif'))){
      qsub(paste0("mp_", year, '_', month, "_", jobnum),
           mo_pred.run,
           pass=c(out_dir, out_dir, jobnum, month, year),
           log=T,
           slots=2,
           # threads=smry_threads,
           # mem=smry_mem,
           submit=T,
           proj=project,
           user='USERNAME')
    }
  }
  check_loc_results(c(1:njobs), pred_mo_dir, prefix="preds_",postfix=".tif")
  
  #summary
  smry_file <- configs$smry_file
  sry.run <- paste0(in_dir, "/FILEPATH.R")
  qsub(paste0("sry", year, '_', month),
       sry.run,
       pass=c(out_dir, out_dir, year, 'TRUE', month),
       log=T,
       slots=smry_slots,
       # threads=smry_threads,
       # mem=smry_mem,
       submit=T,
       proj=project,
       user='USERNAME')
}

month_sample <- function(month_start, month_end, year_start, year_end) {
  make_vector <- function(string) {
    vector <- as.numeric(unlist(strsplit(string, split = ";")))
  }
  
  month_start <- make_vector(month_start)
  month_end <- make_vector(month_end)
  year_start <- make_vector(year_start)
  year_end <- make_vector(year_end)
  
  total_months <- data.frame(month = numeric(0), year = numeric(0))
  for (i in 1:length(month_start)) {
    year_advance <- year_end[i] - year_start[i]
    month_end_temp <- month_end[i] + 12*year_advance
    month_list <- month_start[i]:month_end_temp
    year_list <- c()
    
    for (j in 1:length(month_list)) {
      month <- month_list[j] - 12*floor(month_list[j]/12)
      year <- year_start[i] + floor(month_list[j]/12)
      if (month == 0) {
        month <- 12
        year <- year - 1
      }
      total_months <- rbind(total_months, c(month, year))
    }
  }
  names(total_months) <- c("month", "year")
  
  random_index <- sample(1:nrow(total_months), size = 1)
  sample_month <- total_months[random_index, "month"]
  sample_year <- total_months[random_index, "year"]
  return(c(sample_month, sample_year))
} 

threshold_buffer <- function(buffer, raster, threshold) {
  #threshold should be 0-100
  raster_final <- rasterize(buffer, raster, getCover = T)
  raster_final[raster_final < threshold] <- 0
  raster_final[raster_final >= threshold] <- 1
  return(raster_final)
}

row_in_df <- function(row, field_list, df_list) {
  row[is.na(row)] <- "NA"
  in_df <- NA
  i <- 1
  while ((is.na(in_df) | !in_df) & (i <= length(df_list))) {
    df <- get(df_list[i])
    df[is.na(df)] <- "NA"
    j <- 1
    while (j <= length(field_list)) {
      df <- df[which(row[,field_list[j]] == df[,field_list[j]]),]
      if (!nrow(df)) {
        in_df <- FALSE
      } else {
        if (j == length(field_list)) {
          in_df <- TRUE
          break
        }
      }
      j <- j + 1
    }
    i <- i + 1
    if (i > length(df_list) & is.na(in_df)) {
      in_df <- TRUE
      break
    }
  }
  return(in_df)
}

which_row_in_df <- function(row, field_list, df) {
  row[is.na(row)] <- "NA"
  in_df <- NA
  return_list <- c()
  df[is.na(df)] <- "NA"
  j <- 1
  while (is.na(in_df)) {
    index_list_new <- which(row[,field_list[j]] == df[,field_list[j]])
    if (j == 1) {
      index_list <- index_list_new
    } else {
      index_list <- index_list[index_list %in% index_list_new]
    }
    if (!length(index_list)) {
      in_df <- FALSE
    }
    j <- j + 1
    if (j > length(field_list) & is.na(in_df)) {
      in_df <- TRUE
      return_list <- index_list
    }
  }
  return(return_list)
}

rvf_date0 <- function(data_file, cutoff = NA) {
  dat_orig <- read.csv(data_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  year <- c()
  month <- c()
  count <- 0
  # sample from w/in time bins: mo_start, yr_start -> mo_end, yr_end
  for (i in 1:nrow(dat_orig)){
    yr_start <- as.character(dat_orig$year_start[i])
    yr_end <- as.character(dat_orig$year_end[i])
    mo_start <- as.character(dat_orig$month_start[i])
    mo_end <- as.character(dat_orig$month_end[i])
    count <- count + 1
    if (!is.na(yr_start) && !is.na(yr_end)){
      if (!is.na(mo_start) && !is.na(mo_end)){
        month_year <- month_sample(mo_start, mo_end, yr_start, yr_end)
        month[i] <- month_year[1]
        year[i] <- month_year[2]
      } else {
        month[i] <- NA
        year[i] <- safe_sample(yr_start:yr_end)
      }
    } else {
      month[i] <- NA
      year[i] <- NA
    }
  }
  month <- as.numeric(month)
  year <- as.numeric(year)
  mo_dens <- density(na.omit(month), from = 0, to = 12)
  mo_freqs <- c()
  for (mo in 1:12){
    freq <- integrate.xy(mo_dens$x, mo_dens$y, mo-1, mo)
    if (freq < 0) freq <- 0
    mo_freqs <- append(freq, mo_freqs)
  }
  # Again ~ years
  yr_range <- min(na.omit(dat_orig$year_start)):max(na.omit(dat_orig$year_end))
  yr_dens <- density(na.omit(year), from = yr_range[1] -1)
  yr_freqs <- c()
  for (yr in yr_range){
    freq <- integrate.xy(yr_dens$x, yr_dens$y, yr-1, yr)
    if (freq < 0) freq <- 0
    yr_freqs <- append(freq, yr_freqs)
    
  }
  
  for (i in 1:nrow(dat_orig)) {
    if (is.na(year[i]) & is.na(month[i])) {
      yr_end <- as.numeric(dat_orig$year_end[i])
      mo_end <- as.numeric(dat_orig$month_end[i])
      if (!is.na(yr_end)) {
        if (!is.na(mo_end)) {
          yr_index <- which(yr_range == yr_end)
          freqs_temp <- yr_freqs[1:yr_index]
          prop <- mo_end/12
          diff <- freqs_temp[yr_index] - prop*freqs_temp[yr_index]
          freqs_temp[yr_index] <- prop*freqs_temp[yr_index]
          freqs_temp <- freqs_temp + diff/length(freqs_temp)
          year[i] <- sample(yr_range[1:yr_index], 1, prob = freqs_temp)
          if (year[i] == yr_end) {
            month[i] <- sample(1:mo_end, 1, prob = mo_freqs[1:mo_end])
          } else {
            month[i] <- sample(1:12, 1, prob = mo_freqs)
          }
        } else {
          yr_index <- which(yr_range == yr_end)
          if (yr_index == 1) {
            year[i] <- yr_end
          } else {
            year[i] <- sample(yr_range[1:yr_index], 1, prob = yr_freqs[1:yr_index])
          }
          month[i] <- sample(1:12, 1, prob = mo_freqs)
        }
      } else {
        year[i] <- sample(yr_range, 1, prob = yr_freqs)
        month[i] <- sample(1:12, 1, prob = mo_freqs)
      }
    } else {
      if (is.na(year[i])) year[i] <- sample(yr_range, 1, prob = yr_freqs)
      if (is.na(month[i])) month[i] <- sample(1:12, 1, prob = mo_freqs)
    }
  }
  return(cbind(month, year))
}

rvf_date1 <- function(data_file, cutoff = NA) {
  dat_orig <- read.csv(data_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))

  year <- c()
  month <- c()
  for (i in 1:nrow(dat_orig)){
    yr_start <- as.character(dat_orig$year_start[i])
    yr_end <- as.character(dat_orig$year_end[i])
    mo_start <- as.character(dat_orig$month_start[i])
    mo_end <- as.character(dat_orig$month_end[i])
    if ((!is.na(dat_orig$diagnostic[i]) & !is.na(dat_orig$clinical[i]))&(dat_orig$diagnostic[i] == "PCR" || dat_orig$clinical[i] == "symptomatic")) {
      if (!is.na(yr_start) & !is.na(yr_end)){
        if (!is.na(mo_start) & !is.na(mo_end)){
          month_year <- month_sample(mo_start, mo_end, yr_start, yr_end)
          month[i] <- month_year[1]
          year[i] <- month_year[2]
        } else {
          month[i] <- NA
          year[i] <- safe_sample(yr_start:yr_end)
        }
      } else {
        year[i] <- NA
        if (!is.na(mo_start) & !is.na(mo_end)) {
          month[i] <- NA
        } else {
          month[i] <- safe_sample(mo_start:mo_end)
        }
      } 
    } else {
      month[i] <- NA
      year[i] <- NA
    }
  }
  
  month <- as.numeric(month)
  year <- as.numeric(year)
  mo_dens <- density(na.omit(month), from = 0, to = 12)
  mo_freqs <- c()
  for (mo in 1:12){
    freq <- integrate.xy(mo_dens$x, mo_dens$y, mo-1, mo)
    if (freq < 0) freq <- 0
    mo_freqs <- append(mo_freqs, freq)
  }
  # Again ~ years
  yr_range <- min(na.omit(dat_orig$year_start)):max(na.omit(dat_orig$year_end))
  yr_dens <- density(na.omit(year), from = yr_range[1] -1)
  yr_freqs <- c()
  for (yr in yr_range){
    freq <- integrate.xy(yr_dens$x, yr_dens$y, yr-1, yr)
    if (freq < 0) freq <- 0
    yr_freqs <- append(yr_freqs, freq)
    
  }
  
  for (i in 1:nrow(dat_orig)) {
    if (is.na(year[i]) & is.na(month[i])) {
      yr_end <- as.numeric(dat_orig$year_end[i])
      mo_end <- as.numeric(dat_orig$month_end[i])
      if (!is.na(yr_end)) {
        if (!is.na(mo_end)) {
          yr_index <- which(yr_range == yr_end)
          freqs_temp <- yr_freqs[1:yr_index]
          prop <- mo_end/12
          diff <- freqs_temp[yr_index] - prop*freqs_temp[yr_index]
          freqs_temp[yr_index] <- prop*freqs_temp[yr_index]
          if (yr_end != yr_range[1]) {
            year[i] <- sample(yr_range[1:yr_index], 1, prob = freqs_temp)
          } else {
            year[i] <- yr_end
          }
          if (year[i] == yr_end) {
            month[i] <- sample(1:mo_end, 1, prob = mo_freqs[1:mo_end])
          } else {
            month[i] <- sample(1:12, 1, prob = mo_freqs)
          }
        } else {
          yr_index <- which(yr_range == yr_end)
          if (yr_index == 1) {
            year[i] <- yr_end
          } else {
            year[i] <- sample(yr_range[1:yr_index], 1, prob = yr_freqs[1:yr_index])
          }
          month[i] <- sample(1:12, 1, prob = mo_freqs)
        }
      } else {
        year[i] <- sample(yr_range, 1, prob = yr_freqs)
        month[i] <- sample(1:12, 1, prob = mo_freqs)
      }
    } else {
      if (is.na(year[i])) year[i] <- sample(yr_range, 1, prob = yr_freqs)
      if (is.na(month[i])) month[i] <- sample(1:12, 1, prob = mo_freqs)
    }
  }
  return(cbind(month, year))
}

rvf_date2 <- function(dat_orig, cutoff = NA, yr_min, yr_max) {
  names(dat_orig)[which(names(dat_orig) == "year")] <- "year_lit"
  gbd_regions <- read.csv("FILEPATH.csv", stringsAsFactors = FALSE)
  
  year <- c()
  month <- c()
  yr_range <- yr_min:yr_max
  for (i in 1:nrow(dat_orig)){
    yr_start <- as.character(dat_orig$year_start[i])
    yr_end <- as.character(dat_orig$year_end[i])
    mo_start <- as.character(dat_orig$month_start[i])
    mo_end <- as.character(dat_orig$month_end[i])
    if ((!is.na(dat_orig$diagnostic[i]) & !is.na(dat_orig$clinical[i]))&(dat_orig$diagnostic[i] == "PCR" || dat_orig$clinical[i] == "symptomatic")) {
      if (!is.na(yr_start) & !is.na(yr_end)){
        if (!is.na(mo_start) & !is.na(mo_end)){
          month_year <- month_sample(mo_start, mo_end, yr_start, yr_end)
          month[i] <- month_year[1]
          year[i] <- month_year[2]
        } else {
          month[i] <- NA
          year[i] <- safe_sample(yr_start:yr_end)
        }
      } else {
        year[i] <- NA
        if (!is.na(mo_start) & !is.na(mo_end)) {
          if (mo_start <= mo_end) {
            month[i] <- safe_sample(mo_start:mo_end)
          } else {
            mos <- 1:12
            no_mos <- mo_end:mo_start
            yeah_mos <- mos[-which(mos %in% no_mos)]
            yeah_mos <- append(yeah_mos, c(mo_start, mo_end))
            month[i] <- safe_sample(yeah_mos)
          }
        } else {
          month[i] <- NA
        }
      } 
    } else {
      month[i] <- NA
      year[i] <- NA
    }
  }
  
  
  
  data.all <- cbind(dat_orig, cbind(month, year))
  row_list <- c()
  for (i in 1:length(month)) {
    if (!is.na(month[i]) & !is.na(year[i])) {
      row_list <- append(row_list, i)
    }
  }
  data.complete <- data.all[row_list,]
  countries <- unique(data.all$country)[!is.na(unique(data.all$country))]
  region_count <- 0
  super_count <- 0
  nothin_count <- 0
  nothin_countries <- c()
  for (i in 1:length(countries)) {
    country <- countries[i]
    data.use <- data.complete[data.complete$country == country,]
    if (nrow(data.use) < cutoff) {
      region <- gbd_regions[gbd_regions$ihme_loc_id == country, "region_id"]
      region_countries <- gbd_regions[gbd_regions$region_id == region, "ihme_loc_id"]
      data.use <- data.complete[data.complete$country %in% region_countries,]
      region_count <- region_count + 1
      if (nrow(data.use) < cutoff) {
        super_region <- gbd_regions[gbd_regions$ihme_loc_id == country, "super_region_id"]
        super_region_countries <- gbd_regions[gbd_regions$super_region_id == super_region, "ihme_loc_id"]
        data.use <- data.complete[data.all$country %in% super_region_countries,]
        super_count <- super_count + 1
        if (nrow(data.use) < cutoff) {
          data.use <- data.complete
          nothin_count <- nothin_count + 1
          nothin_countries <- append(nothin_countries, country)
        }
      }
    }
    mo_list <- sort(data.use$month)
    middle <- mo_list[ceiling(length(mo_list)/2)]
    shift <- 6 - middle
    shift_list <- mo_list + shift
    for (j in 1:length(shift_list)) {
      if (shift_list[j] > 12) shift_list[j] <- shift_list[j] -12
      if (shift_list[j] < 1) shift_list[j] <- shift_list[j] +12
    }
    
    mo_dens <- density(na.omit(shift_list), from = 0, to = 12)
    mo_freqs <- c()
    for (mo in 1:12){
      freq <- integrate.xy(mo_dens$x, mo_dens$y, mo-1, mo)
      if (freq < 0) freq <- 0
      mo_freqs <- append(mo_freqs, freq)
    }
    if (shift != 0) mo_freqs <- c(tail(mo_freqs, -shift), head(mo_freqs, shift))
    
    yr_dens <- density(na.omit(data.use$year), from = (yr_min - 1), to = yr_max)
    yr_freqs <- c()
    for (yr in yr_range){
      freq <- integrate.xy(yr_dens$x, yr_dens$y, yr-1, yr)
      if (freq < 0) freq <- 0
      yr_freqs <- append(yr_freqs, freq)
    }
    
    new_row_yr <- cbind(country, t(yr_freqs))
    if (i == 1) {
      year_frequencies <- data.frame(new_row_yr, stringsAsFactors = F)
    } else {
      year_frequencies <- rbind(year_frequencies, new_row_yr)
    }
    
    new_row_mo <- cbind(country, t(mo_freqs))
    if (i == 1) {
      month_frequencies <- data.frame(new_row_mo, stringsAsFactors = F)
    } else {
      month_frequencies <- rbind(month_frequencies, new_row_mo)
    }
  }
  
  for (i in 1:nrow(dat_orig)) {
    if (!is.na(data.all$country[i])) {
      mo_freqs <- as.numeric(month_frequencies[month_frequencies$country == data.all$country[i], 2:13])
      yr_freqs <- as.numeric(year_frequencies[year_frequencies$country == data.all$country[i], -1])
      if (is.na(year[i]) & is.na(month[i])) {
        print(i)
        yr_end <- as.numeric(dat_orig$year_end[i])
        mo_end <- as.numeric(dat_orig$month_end[i])
        if (!is.na(yr_end)) {
          if (!is.na(mo_end)) {
            yr_index <- which(yr_range == yr_end)
            freqs_temp <- yr_freqs[1:yr_index]
            prop <- mo_end/12
            diff <- freqs_temp[yr_index] - prop*freqs_temp[yr_index]
            freqs_temp[yr_index] <- prop*freqs_temp[yr_index]
            if (yr_end != yr_range[1]) {
              year[i] <- sample(yr_range[1:yr_index], 1, prob = freqs_temp)
            } else {
              year[i] <- yr_end
            }
            if (year[i] == yr_end) {
              month[i] <- sample(1:mo_end, 1, prob = mo_freqs[1:mo_end])
            } else {
              month[i] <- sample(1:12, 1, prob = mo_freqs)
            }
          } else {
            yr_index <- which(yr_range == yr_end)
            if (yr_index == 1) {
              year[i] <- yr_end
            } else {
              year[i] <- sample(yr_range[1:yr_index], 1, prob = yr_freqs[1:yr_index])
            }
            month[i] <- sample(1:12, 1, prob = mo_freqs)
          }
        } else {
          year[i] <- sample(yr_range, 1, prob = yr_freqs)
          month[i] <- sample(1:12, 1, prob = mo_freqs)
        }
      } else {
        if (is.na(year[i])) year[i] <- sample(yr_range, 1, prob = yr_freqs)
        if (is.na(month[i])) month[i] <- sample(1:12, 1, prob = mo_freqs)
      }
    } else {
      year[i] <- NA
      month[i] <- NA
    }
  }
  return(cbind(month, year))
}

rvf_date <- function(data_file, method, cutoff, yr_min = NULL, yr_max = NULL) {
  if (method == 0) {
    dates <- rvf_date0(data_file)
  } else if (method == 1) {
    dates <- rvf_date1(data_file)
  } else if (method == 2) {
    dates <- rvf_date2(data_file, cutoff, yr_min, yr_max)
  }
  return(dates)
}

get_freqs <- function(dat_orig, yr_min, yr_max) {

  year <- c()
  month <- c()
  for (i in 1:nrow(dat_orig)){
    yr_start <- as.character(dat_orig$year_start[i])
    yr_end <- as.character(dat_orig$year_end[i])
    mo_start <- as.character(dat_orig$month_start[i])
    mo_end <- as.character(dat_orig$month_end[i])
    if ((!is.na(dat_orig$diagnostic[i]) & !is.na(dat_orig$clinical[i]))&(dat_orig$diagnostic[i] == "PCR" || dat_orig$clinical[i] == "symptomatic")) {
      if (!is.na(yr_start) & !is.na(yr_end)){
        if (!is.na(mo_start) & !is.na(mo_end)){
          month_year <- month_sample(mo_start, mo_end, yr_start, yr_end)
          month[i] <- month_year[1]
          year[i] <- month_year[2]
        } else {
          month[i] <- NA
          year[i] <- safe_sample(yr_start:yr_end)
        }
      } else {
        year[i] <- NA
        if (!is.na(mo_start) & !is.na(mo_end)) {
          if (mo_start <= mo_end) {
            month[i] <- safe_sample(mo_start:mo_end)
          } else {
            mos <- 1:12
            no_mos <- mo_end:mo_start
            yeah_mos <- mos[-which(mos %in% no_mos)]
            yeah_mos <- append(yeah_mos, c(mo_start, mo_end))
            month[i] <- safe_sample(yeah_mos)
          }
        } else {
          month[i] <- NA
        }
      } 
    } else {
      month[i] <- NA
      year[i] <- NA
    }
  }
  
  month <- as.numeric(month)
  year <- as.numeric(year)
  mo_dens <- density(na.omit(month), from = 0, to = 12)
  mo_freqs <- c()
  for (mo in 1:12){
    freq <- integrate.xy(mo_dens$x, mo_dens$y, mo-1, mo)
    if (freq < 0) freq <- 0
    mo_freqs <- append(mo_freqs, freq)
  }
  # Again ~ years
  yr_range <- yr_min:yr_max
  yr_dens <- density(na.omit(year), from = yr_range[1] -1)
  yr_freqs <- c()
  for (yr in yr_range){
    freq <- integrate.xy(yr_dens$x, yr_dens$y, yr-1, yr)
    if (freq < 0) freq <- 0
    yr_freqs <- append(yr_freqs, freq)
    
  }
  return(append(mo_freqs, yr_freqs))
}

library(PresenceAbsence)

auc2 <- function (DATA,
                  st.dev = TRUE,
                  which.model = 1,
                  na.rm = FALSE) {
  if (is.logical(st.dev) == FALSE) {
    stop ("'st.dev' must be of logical type")
  }
  if (is.logical(na.rm) == FALSE) {
    stop ("'na.rm' must be of logical type")
  }
  if (sum(is.na(DATA)) > 0) {
    if (na.rm == TRUE) {
      NA.rows <- apply(is.na(DATA), 1, sum)
      warning (length(NA.rows[NA.rows > 0]), " rows ignored due to NA values")
      DATA <- DATA[NA.rows == 0, ]
    } else {
      return (NA)
    }
  }
  if (length(which.model) != 1) {
    stop ("this function will only work for a single model, 'which.model' must be of length one")
  }
  if (which.model < 1 || round(which.model) != which.model) {
    stop ("'which.model' must be a positive integer")
  }
  if (which.model + 2 > ncol(DATA)) {
    stop ("'which.model' must not be greater than number of models in DATA")
  }
  DATA <- DATA[, c(1, 2, which.model + 2)]
  DATA[DATA[, 2] > 0, 2] <- 1
  OBS <- DATA[, 2]
  PRED <- DATA[, 3]
  if (length(OBS[OBS == 1]) == 0 || length(OBS[OBS == 1]) == 
      nrow(DATA)) {
    if (st.dev == FALSE) {
      return (NaN)
    } else {
      return (data.frame(AUC = NaN, AUC.sd = NaN))
    }
  }
  rm(DATA)
  PRED.0 <- PRED[OBS == 0]
  PRED.1 <- PRED[OBS == 1]
  N <- length(PRED)
  n0 <- as.double(length(PRED.0))
  n1 <- as.double(length(PRED.1))
  R <- rank(PRED, ties.method = "average")
  R0 <- R[OBS == 0]
  R1 <- R[OBS == 1]
  U <- n0 * n1 + (n0 * (n0 + 1)) / 2 - sum(R0)
  AUC <- U / (n0 * n1)

  rm(PRED)
  rm(OBS)
  if (st.dev == FALSE) {
    return (AUC = AUC)
  } else {
    RR0 <- rank(PRED.0, ties.method = "average")
    RR1 <- rank(PRED.1, ties.method = "average")
    pless.0 <- (R0 - RR0) / n1
    pless.1 <- (R1 - RR1) / n0
    var.0 <- var(pless.0)
    var.1 <- var(pless.1)
    var.AUC <- (var.0 / n0) + (var.1 / n1)
    st.dev.AUC <- var.AUC ^ 0.5
    return (data.frame(AUC = AUC, AUC.sd = st.dev.AUC))
  }
}

rmse <- function(truth, prediction)
  # root mean squared error of prediction from true probability
{
  sqrt(mean(abs(prediction - truth) ^ 2))
}

devBern <- function (truth, prediction)
  # predictive deviance from a bernoulli distribution
{
  -2 * sum(dbinom(truth, 1, prediction, log = TRUE))
}

calcStats <- function(df) {
  
  # if any elements of df are NAs, return NAs
  if (any(is.na(df))) {
    
    results <- c(deviance = NA,
                 rmse = NA,
                 kappa = NA,
                 auc = NA,
                 sens = NA,
                 spec = NA,
                 pcc = NA,
                 kappa_sd = NA,
                 auc_sd = NA,
                 sens_sd = NA,
                 spec_sd = NA,
                 pcc_sd = NA,
                 thresh = NA)
    
  } else {
    
    # add an id column (needed for PresenceAbsence functions)
    df <- data.frame(id = 1:nrow(df), df)
    
    # ~~~~~~~~~~~
    # copntinuous probability metrics
    
    # bernoulli deviance
    dev <- devBern(df[, 2], df[, 3])
    
    # root mean squared error
    rmse <- rmse(df[, 2], df[, 3])
    
    # auc (using my safe version - see above)
    auc <- auc2(df, st.dev = TRUE)
    
    # ~~~~~~~~~~~  
    # discrete classification metrics
    
    # calculate the 'optimum' threshold - one at which sensitivity == specificity
    opt <- optimal.thresholds(df, threshold = 101, which.model = 1, 
                              opt.methods = 3)
    
    # create confusiuon matrix at this threshold
    confusion <- cmx(df, threshold = opt[1, 2])
    
    # kappa (using threshold)
    kappa <- Kappa(confusion, st.dev = TRUE)
    
    # sensitivity and specificity using threshold
    sens <- sensitivity(confusion, st.dev = TRUE)
    spec <- specificity(confusion, st.dev = TRUE)
    
    # proportion correctly classified using threshold
    pcc <- pcc(confusion, st.dev = TRUE)
    
    # create results vector
    results <- c(deviance = dev,
                 rmse = rmse,
                 kappa = kappa[, 1],
                 auc = auc[, 1],
                 sens = sens[, 1],
                 spec = spec[, 1],
                 pcc = pcc[, 1],
                 kappa_sd = kappa[, 2],
                 auc_sd = auc[, 2],
                 sens_sd = sens[, 2],
                 spec_sd = spec[, 2],
                 pcc_sd = pcc[, 2],
                 thresh = opt[1, 2])
    
    
  }
  
  # and return it
  return (results)
}
