
## Various Useful Functions

## Make time stamp in standardized format.
make_time_stamp <- function(time_stamp) {

  run_date <- gsub("-","_",Sys.time())
  run_date <- gsub(":","_",run_date)
  run_date <- gsub(" ","_",run_date)

  if(time_stamp==FALSE) run_date <- 'scratch'

  return(run_date)
}

# given a year and month, return a Date object of the first day of that month
firstDay <- function (year, month) {
  date_string <- paste(year, month, '01', sep = '-')
  date <- as.Date (date_string)
  return (date)
}

# given a year and month, return a Date object of the last day of that month
lastDay <- function (year, month) {
  next_month <- ifelse(month == 12, 1, month + 1)
  next_year <- ifelse(month == 12, year + 1, year)
  next_date_string <- paste(next_year, next_month, '01', sep = '-')
  next_date <- as.Date(next_date_string)
  date <- next_date - 1
  return (date)
}

# given a vector of text strings `text`, convert to sentence case
sentenceCase <- function (text) {

  # convert all to lower case
  text <- tolower(text)
  # split at spaces
  text_list <- strsplit(text, ' ')

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

# given a vector of text strings `text` subset each to only the first two
# words (bits separated by spaces) and return this as a vector.
firstTwo <- function (text) {
  text <- as.character(text)
  text_list <- strsplit(text, ' ')
  text_list <- lapply(text_list, '[', 1:2)
  text_list <- lapply(text_list, paste, collapse = ' ')
  text_vector <- unlist(text_list)
  return (text_vector)
}

# first buffer, then rasterize IUCN range map for species
rasterizeSpecies <- function(species,
                             shape,
                             raster,
                             buffer = NULL,
                             folder = 'FILEPATH/') {

  shp <- shape[shape$BINOMIAL == species, ]

  if (!is.null(buffer)) {
    # convert buffer from kilometres to decimal degrees, assuming at the equator
    buffer <- buffer / 111.32
    # buffer by this amount
    shp <- gBuffer(shp, width = buffer)
  }

  # rasterize the shapefile
  tmp <- rasterize(shp, raster, field = 1, background = 0, fun = 'first')
  writeRaster(tmp,
              filename = paste0('~/FILEPATH/',
                                folder,
                                '/',
                                gsub(' ', '_', species)),
              format = 'GTiff',
              overwrite = TRUE)
  rm(tmp)
  return (NULL)
}

# load a raster if it contains any of the species' range,
# mask and resave it, else delete it
tidySpecies <- function (filename, template) {
  tmp <- raster(filename)
  
  if (!is.na(maxValue(tmp)) && maxValue(tmp) == 1) {
    tmp <- mask(tmp, template)
    writeRaster(tmp, file = filename, overwrite = TRUE)
  } else {
    rm(tmp)
    file.remove(filename)
  }
  return (NULL)
}

# given a presence-background dataset, with multiple rows for some of the
# occurrence records, subset it so that there's only one randomly selected
# point from each polygon and then take a bootstrap it using `subsample`.
# Dots argument is passed to subsample.
subsamplePolys <- function (data, ...) {

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
  pref <- substitute(SPDF)  
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

# given a raster `risk_raster` giving risk on the (0,1] level and a 2d
# dataframe `occ` with columns named 'lat' and 'long' giving the latitudes
# and longitudes of known occurrence records,
# find the threshold value so that `proportion` fraction of the records
# fall in areas classified as 'at risk' and return the thresholded map
thresholdRisk <- function (risk_raster,
                           occ,
                           proportion = 1) {

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

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

path_converter <- function(winpath, disc){
  if (disc == "ADDRESS"){
    path <- gsub("ADDRESS:", "/ADDRESS", gsub("\\\\", "/", winpath))
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
    fx <- fx[!i]
  }
  if(any(diff(x) == 0))
    stop("bug in 'duplicated()' still have multiple x[]!")

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
      if(TRUE) cat("working around spline(.)")
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


  dig0 <- floor(-log10(xtol)) #
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

threshold_buffer <- function(buffer, raster, threshold) {
  #threshold should be 0-100
  raster_final <- rasterize(buffer, raster, getCover = T)
  raster_final[raster_final < threshold] <- 0
  raster_final[raster_final >= threshold] <- 1
  return(raster_final)
}

get_fit_stats <- function(dat_all, pred_raster, mo = "00") {
  monthly <- FALSE
  if (mo != "00") monthly <- TRUE

  dat0.pts <- dat_all[dat_all$PA==0, c('long', 'lat')]
  dat0.preds <- raster::extract(pred_raster, dat0.pts)
  
  dat1.pts <- dat_all[dat_all$PA==1, c('long', 'lat')]
  dat1.preds <- raster::extract(pred_raster, dat1.pts)
  
  pos_mean_preds <- na.omit(dat1.preds)
  neg_mean_preds <- na.omit(dat0.preds)
  
  sensitivity <- c()
  fpr <- c()
  dif = 0.01
  for (thresh in seq(0.00, 1, dif)){
    tp <- length(pos_mean_preds[which(pos_mean_preds >= thresh)])
    fn <- length(pos_mean_preds[which(pos_mean_preds < thresh)])
    tn <- length(neg_mean_preds[which(neg_mean_preds <= thresh)])
    fp <- length(neg_mean_preds[which(neg_mean_preds > thresh)])
    sensitivity <- append(sensitivity, tp / (tp + fn))
    fpr <- append(fpr, 1 - (tn / (fp + tn)))
  }
  
  distance <- c()
  for (i in 1:length(fpr)){
    distance <- append(distance, dist(rbind(c(fpr[i], sensitivity[i]), c(0, 1))))
  }
  opt_thresh <- max(seq(0.00, 1, dif)[which(distance == min(distance))])
  auc <- integrate.xy(fpr, sensitivity, 0, 1)
  thresh <- opt_thresh
  tp <- length(pos_mean_preds[which(pos_mean_preds >= thresh)])
  fn <- length(pos_mean_preds[which(pos_mean_preds < thresh)])
  tn <- length(neg_mean_preds[which(neg_mean_preds <= thresh)])
  fp <- length(neg_mean_preds[which(neg_mean_preds > thresh)])
  rmse <- sqrt(mean(append((1-pos_mean_preds), (neg_mean_preds - 0))^2))
  error_rate <- (fp+fn)/(tp+tn+fp+fn)
  
  fit_stats <- data.frame(opt_thresh, auc, rmse, error_rate)
  return(fit_stats)
}
