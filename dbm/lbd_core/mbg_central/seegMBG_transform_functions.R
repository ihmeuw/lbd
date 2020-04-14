# functions to transform covariates

#' @name pcaTrans
#' @rdname pcaTrans
#'
#' @title PCA transform a \code{Raster*} object
#'
#' @description Given a matrix of point coordinates and a Raster* object of
#'  covariates, return a \code{Raster*} object of the same size with the layers
#'  giving principal components from a PCA rotation based on the data values
#'  at the coordinates.
#'
#' @param coords a two-column matrix or dataframe giving the location about
#' which to carry out the PCA analysis
#'
#' @param covs a \code{Raster*} object containing covariates to rotate
#'
#' @return a \code{Raster*} object with the same extent, resolutiona and number
#'  of layers as \code{covs} but with each layer giving the location on a
#'  different principal component axis.
#'
#' @family transform
#'
#' @export
#' @import raster
#'
pcaTrans <- function(coords, covs) {
  
  
  # get covariate data at point locations
  vals <- extract(covs, coords)
  vals <- na.omit(vals)
  
  # do pca analysis
  pca <- prcomp(vals, retx = FALSE, center = TRUE, scale. = TRUE)
  
  # find non-missing cells
  cell_idx <- cellIdx(covs)
  
  # extract covariate values
  vals <- raster::extract(covs, cell_idx)
  
  # convert to a data.frame
  vals <- data.frame(vals)
  
  # get PCA predictions
  vals_trans <- predict(pca, vals)
  
  # set new raster values
  trans_ras <- insertRaster(raster = covs,
                            new_vals = vals_trans,
                            idx = cell_idx)
  return (trans_ras)
  
}



#' @name gamTrans
#' @rdname gamTrans
#'
#' @title Carry out a model-based covariate transformation using a GAM
#'
#' @description Define an optimal set of univariate covariate
#'  transformations of a set of model covariates by fitting a generalised
#'  additive model with univariate smoothers to data, and then using the
#'  smoothers to spline-transform the covariates.
#'  This makes use of the\code{type = 'terms'} argument in
#'  \code{\link{predict.gam}}.
#'  This function also makes use of
#'
#' @param coords a two-column matrix of coordinates of records
#'
#' @param response an object acting as thge response object in the GAM
#'  model (e.g. a vector of counts, or a matrix for binomial data)
#'
#' @param covs a \code{Raster*} object giving the spatial covariates
#'  for the main part of the model
#'
#' @param family the distribution family for the gam
#'
#' @param condition an optional vector of 1s and 0s of the same length as
#'  the number of records in \code{coords} and \code{response} and stating
#'  whether the record should also be modelled using covariates in
#'  \code{condition_covs} (1 if so and 0 if not). This enables the construction
#'  of slightly more complex models, such as those with an explicitly modelled
#'  observation process. This is achieved by passing \code{condition} to the
#'  \code{by} argument in \code{mgcv::s} when fitting smooths for
#'  the condition covariates, as well as adding the condition as an intercept.
#'
#' @param condition_covs an optional \code{Raster*} object giving the spatial covariates
#'  for the conditional part of the model
#'
#' @param extra_terms an optional formula object (of the form \code{~ s(x, k = 2)}
#'  or similar which can be concatenated onto the model formula)
#'  specifying further model components (in \code{extra_data}) not provided in
#'  the spatial covariates.
#'
#' @param extra_data an optional dataframe giving the covariates referred to in
#'  \code{extra_terms}
#'
#' @param bam whether to fit the model using \code{mgcv::bam} (the default),
#'  otherwise \code{mgcv::gam} is used instead
#'
#' @param s_args a named list of additional arguments to pass to the smoother on
#'   each covariate. For example, this may include the smoother type (\code{bs})
#'   or the basis dimension (\code{k}). See \code{\link[mgcv]{s}} for the list
#'   of available arguments.
#'
#' @param predict whether to transform the rasters after fitting the model.
#'  If set to \code{FALSE} this can enable model tweaking before the final
#'  transformations are applied, without the computational cost of prediction
#'
#' @param \dots other arguments to be passed to \code{mgcv::bam} or
#'  \code{mgcv::gam}
#'
#' @return a three-element named list containing:
#'  \itemize{
#'    \item{model}{the fitted \code{bam} or \code{gam} model object}
#'    \item{trans}{if \code{predict = TRUE} a \code{Raster*} object of the
#'     same extent, resolution and number of layers as \code{covs}, but with
#'     the values of each layer having been optimally spline-transformed.
#'     Otherwise \code{NULL}}
#'    \item{trans_cond}{if \code{predict = TRUE} and \code{condition} is not
#'     \code{NULL} a \code{Raster*} object of the same extent, resolution and
#'     number of layers as \code{condition_covs}, but with the values of each layer
#'     having been optimally spline-transformed. Otherwise \code{NULL}}
#'  }
#'
#' @export
#' @import mgcv
#' @import raster
#'
gamTrans <- function(coords,
                     response,
                     covs,
                     family = gaussian,
                     condition = NULL,
                     condition_covs = NULL,
                     extra_terms = NULL,
                     extra_data = NULL,
                     bam = TRUE,
                     s_args = list(),
                     predict = TRUE,
                     ...) {
  
  # test to see if covs are passed as a rasterbrick or a list
  # an RBrick indicates there are no temporal covariates
  # a list indicates either a mix, or only temporal
  # within the list RLayers are non-temporal, RBricks are temporally varying
  if(inherits(covs,'list')) {
    temporal = TRUE
  } else {
    temporal = FALSE
  } 
  
  # run a test to see that we have years in the extra_data if we have temporal covs
  if(temporal & is.null(extra_data))
    stop('You have temporally-varying covariates, but not temporally varying data. Please include in extra data argument.')
  
  # test we have the same number of periods in data as in all temporal covs
  if(temporal & !is.null(extra_data)){
    for(rb in covs) {
      if(class(rb)=="RasterBrick"&dim(rb)[3]!=length(unique(extra_data$year)))
        stop('One of your temporally-varying bricks does not have the correct number of layers. Assumes year is the temporal variable in your extra_data.')
    }
  }
  
  
  # whether there's a conditional bit
  cond <- !is.null(condition)
  
  stopifnot(inherits(extra_terms, 'formula'))
  
  # check inputs
  stopifnot(inherits(covs, 'Raster')|inherits(covs, 'list'))
  if (cond)
    stopifnot(inherits(condition_covs, 'Raster'))
  
  # add 'cond_' onto the conditional covariate names to prevent naming conflicts
  # with the disease model
  if (cond)
    names(condition_covs) <- paste0('cond_', names(condition_covs))
  
  # get covariate names
  cov_names <- names(covs)
  if (cond)
    cond_names <- names(condition_covs)
  
  
  # ~~~~~~~~~~~~~
  # build formula
  
  cov_terms_string <- paste(sprintf('s(%s, %s)',
                                    cov_names,
                                    parseArgsS(s_args)),
                            collapse = ' + ')
  
  # cov_terms <- reformulate(cov_terms_string)
  # f <- response ~ 1
  # f <- f + cov_terms
  
  # if required, add extra terms
  if (!is.null(extra_terms))
    cov_terms_string=paste(cov_terms_string,'+',as.character(extra_terms)[2])
  
  
  f<-formula(paste('response~1+',cov_terms_string))
  
  
  # if required, add conditional terms
  if (cond) {
    cond_terms_string <- paste(sprintf('s(%s, %s, by = condition)',
                                       cond_names,
                                       parseArgsS(s_args)),
                               collapse = ' + ')
    
    cond_terms <- reformulate(cond_terms_string)
    
    f <- f + cond_terms + ~ condition_intercept
  }
  
  # assign any objects in the arguments of l into this environment
  # so they can be accessed by gam/bam
  if (length(s_args) > 0) {
    for (i in 1:length(s_args))
      assign(names(s_args)[i], s_args[[i]])
  }
  
  # ~~~~~~~~~~~~~
  # get training data
  
  # extraction for temporally varying and non-temporally varing covariates is slightly different
  message('Extracting covariates at data locations')
  if(temporal){
    # split apart temporally varying and non temporally varying covariates
    nT_covs=T_covs=list()
    for(c in names(covs)){
      if(class(covs[[c]])=="RasterLayer")
        nT_covs[[c]]=covs[[c]]
      if(class(covs[[c]])=="RasterBrick")
        T_covs[[c]]=covs[[c]]
    }
    time.varying.cov.names=names(T_covs)
    # initiate 
    if (!is.null(extra_data)){
      data <- extra_data
    } else {
      stop('you have temporal covariates with no year extra data terms.')
    }
    
    
    # begin data with nontemporally varying because that is easy
    if(length(nT_covs)>0){
      nT_covs=brick(nT_covs)
      crs(nT_covs)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
      data <- data.frame(cbind(data,extract(nT_covs, coords)))
    }
    
    
    # Note, we expect rasterbricks to go from earliest to latest period order, it will match with watch we get in extra_data to do this
    periods=sort(unique(data$year))
    tnames=c()
    for(n in names(T_covs)){
      names(T_covs[[n]])=paste0(n,periods)
      tnames=c(tnames,names(T_covs[[n]]))
      data <- data.frame(cbind(data,extract(T_covs[[n]], coords)))
    }
    
    # match years of time varying covariates and data points, then clean up the data frame
    # assumes a 4 year time period (years, in most cases 2000,2005,2010,2015)
    for(n in tnames)
      data[,n][data$year!= as.numeric(substr(n,nchar(n)-3,nchar(n)))]=NA
    for(n in names(T_covs))
      data[,n]=rowSums(data[,paste0(n,periods)],na.rm=T)
    
    data=data[,!names(data) %in% tnames]
    
    # save all covs in a brick, will need them later
    covs<-brick(covs)
  }
  
  
  
  # extract covariates if have no temporally varying
  if(!temporal) {
    
    data <- data.frame(extract(covs, coords))
    
    if (!is.null(extra_data)){
      data <- extra_data
    } 
  } 
  
  
  # optionally combine this with the conditional and extra data
  if (cond)
    data <- cbind(data,
                  data.frame(extract(condition_covs,
                                     coords)),
                  condition_intercept = condition,
                  condition)
  
  
  
  
  # find any missing values and remove corresponding rows
  rem_idx <- badRows(data)
  data <- data[!rem_idx, ]
  
  if (is.vector(response)) {
    response <- response[!rem_idx]
  } else {
    response <- response[!rem_idx, ]
  }
  
  
  
  message('Fitting GAM')
  # fit the model
  if (bam) {
    m <- mgcv::bam(f, data = data, family = family)
  } else {
    m <- mgcv::gam(f, data = data, family = family)
  }
  
  
  
  # ~~~~~~~~~~
  # optionally apply transformations
  
  if (predict) {
    
    # find index for non-missing cells
    cell_idx <- cellIdx(covs)
    
    
    
    # transform the main covariates
    
    # extract covariate values
    message('Extracting Covariates at all cells (takes a few minutes)')
    vals <- raster::extract(covs, cell_idx)
    
    
    # convert to a data.frame
    vals <- data.frame(vals)
    
    # add on the extra data
    vals <- cbind(vals, extra_data[rep(1, nrow(vals)), , drop = FALSE])
    
    # optionally add on dummy data for the condition intercept and variables
    if (cond) {
      cond_data <- data.frame(matrix(0,
                                     nrow = 1,
                                     ncol = length(cond_names)))
      names(cond_data) <- cond_names
      vals <- cbind(vals,
                    condition = 0,
                    condition_intercept = 0,
                    cond_data)
    }
    
    
    
    # get the transformations of these values
    message('Grab Transformed Values')
    
    if( temporal ) {
      
      # need to account for time-varying variables here first, predict with 
      # them in roating order
      # loop through the years of covariates to transform them. 
      valslist = list()
      for(p in 1:length(periods)){ # 1:4
        tmp=vals
        for(n in names(covs)){  # loop through all covariate names
          if(length(grep(paste0('.',p),n))>0){    # find ones that have the period number (ie 1-4, not in 2000-2015 format) in them
            names(tmp)[names(tmp)==n]=substr(n,1,nchar(n)-2)
          }
          if(length(grep(paste(paste0('.',1:4),collapse='|'),n))>0){           # drop all variables not in the time bin
            tmp[,n]=NULL
          }
        }
        valslist[[length(valslist)+1]]=data.frame(tmp)
        
      }
      # so now we have a list of suitable data frames with which to predict out transformations
      
      # thow a check for prediction frame names matching those in the estimation data used
      
      # predict them out  
      vals_trans=list()
      for(p in 1:length(periods)){ # 1:4
        
        # quick fix, BAD BAD, Get lucas to make sure we have no NA values
        for(n in names(T_covs)){
          if(length(valslist[[p]][,n][is.na(valslist[[p]][,n])])!=0) {
            message(paste0(n, ': ', length(valslist[[p]][,n][is.na(valslist[[p]][,n])]), ' rows missing!'))
            message('...filling in with nearby average. Explore why you have NAs in this covariate.')
            valslist[[p]][,n][is.na(valslist[[p]][,n])]=mean(valslist[[p]][,n],na.rm=T) # do the ones around it..
          }
        }
        for(n in names(nT_covs)){
          if(length(valslist[[p]][,n][is.na(valslist[[p]][,n])])!=0) {
            message(paste0(n, ': ', length(valslist[[p]][,n][is.na(valslist[[p]][,n])]), ' rows missing!'))
            valslist[[p]][,n][is.na(valslist[[p]][,n])]=mean(valslist[[p]][,n],na.rm=T) # do the ones around it..
            message('...filling in with nearby average. Explore why you have NAs in this covariate.')
          }
        }
        
        
        message(paste0('Predicting Transform for Time Period ',p))
        
        # vals_trans[[p]] <- predict(m,
        #                            newdata = valslist[[p]],
        #                            type = 'terms',exclude='s(fertility)')
        vals_trans[[p]] <- predict(m,
                                   newdata = valslist[[p]],
                                   type = 'terms')
        
        # format the names
        colnames(vals_trans[[p]]) <- gsub(')', '', colnames(vals_trans[[p]]))
        colnames(vals_trans[[p]]) <- gsub('^s\\(', '', colnames(vals_trans[[p]]))
      }
      message('Prediction Complete')
      
      # now explode back out the temporally varying covariates and prep them to be returned in the same format they came in
      # in other words a list with RLs and RBs.. 
      message('Raster Brick non temporally varying covariates')
      nT_vals_trans=as.matrix(vals_trans[[1]][,colnames(vals_trans[[1]]) %in% names(nT_covs)]) # first remove all temporals
      nT_vals_trans=insertRaster(raster = covs,
                                 new_vals = nT_vals_trans,
                                 idx = cell_idx)
      
      
      # save temporals in their own lists, to be bricked
      message('Raster Brick and List temporally varying covariates')
      T_vals_trans=list()
      for(n in names(T_covs)){
        print(n)
        temp=data.frame(id=1:nrow(vals_trans[[1]]))
        for(p in 1:length(periods)){
          temp[,paste0(n,'.',p)]=vals_trans[[p]][,n]
        }
        temp$id=NULL
        assign(n,temp)
        T_vals_trans[[n]]=insertRaster(raster = covs,
                                       new_vals = as.matrix(get(n)),
                                       idx = cell_idx)
      }
      
      
      # set new raster values
      trans_ras <- list(nT_vals_trans=nT_vals_trans,T_vals_trans=T_vals_trans)
      
      
      
      
      
      
      
      
      
    } else { # if not temporal
      vals_trans <- predict(m,
                            newdata = vals,
                            type = 'terms')
      # find any condition terms
      cond_terms_idx <- grep('):condition$', colnames(vals_trans))
      
      # remove the condition ones and the condition index
      if (length(cond_terms_idx) > 0) {
        vals_trans <- vals_trans[, -cond_terms_idx]
      }
      vals_trans <- vals_trans[, !(colnames(vals_trans) %in% c('condition', 'condition_intercept'))]
      
      # format the names
      colnames(vals_trans) <- gsub(')', '', colnames(vals_trans))
      colnames(vals_trans) <- gsub('^s\\(', '', colnames(vals_trans))
      
      
      # keep only the names that are in the raster
      # gets rid of year and age bin
      vals_trans <- vals_trans[, colnames(vals_trans) %in% cov_names]
      
      # set new raster values
      trans_ras <- insertRaster(raster = covs,
                                new_vals = vals_trans,
                                idx = cell_idx)
    }
    
    
    # optionally apply the transformation to the conditional terms
    
    if (cond) {
      
      # remove the previous variables
      rm(vals, vals_trans)
      
      # extract covariate values
      vals <- raster::extract(condition_covs, cell_idx)
      
      # convert to a data.frame
      vals <- data.frame(vals)
      
      # add on the extra data
      vals <- cbind(vals, extra_data[rep(1, nrow(vals)), ])
      
      # add on dummy data for the main variables
      main_data <- data.frame(matrix(0,
                                     nrow = 1,
                                     ncol = length(cov_names)))
      
      names(main_data) <- cov_names
      
      vals <- cbind(vals,
                    condition = 1,
                    condition_intercept = 0,
                    main_data)
      
      # get the transformations of these values
      vals_trans <- predict(m,
                            newdata = vals,
                            type = 'terms')
      
      # keep only the condition ones
      vals_trans <- vals_trans[, grep('):condition$', colnames(vals_trans))]
      
      # format the names
      colnames(vals_trans) <- gsub('):condition$', '', colnames(vals_trans))
      colnames(vals_trans) <- gsub('^s\\(', '', colnames(vals_trans))
      
      # keep only the names that are in the raster
      vals_trans <- vals_trans[, colnames(vals_trans) %in% cond_names]
      
      # remove the `cond_` bit
      colnames(vals_trans) <- gsub('^cond_', '', colnames(vals_trans))
      
      # set new raster values
      trans_cond_ras <- insertRaster(raster = condition_covs,
                                     new_vals = vals_trans,
                                     idx = cell_idx)
      
    } else { # if not cond
      trans_cond_ras <- NULL
    }
    
  } else {    # if not poredicting, set these to NULL
    
    trans_ras <- trans_cond_ras <- NULL
    
  }
  
  # return the three components
  return (list(model = m,
               trans = trans_ras,
               trans_cond = trans_cond_ras))
  
}


addQuotes <- function (x) {
  # If x is a character string, add (escaped) quotation marks
  if (is.character(x)) {
    x <- sprintf('\"%s\"', x)
  }
  return (x)
}

parseArgsS <- function(l) {
  # parse a list of additional arguments to smoothers in gamTrans
  stopifnot(is.list(l))
  l_string <- paste(names(l),
                    lapply(l, addQuotes),
                    sep = ' = ',
                    collapse = ', ')
  return (l_string)
}
