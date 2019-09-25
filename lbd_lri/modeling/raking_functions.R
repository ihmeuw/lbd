#' a function to rake MBG results linearly. Where gbd = sum(mbg) * scalar
#'
#' @param cell_pred matrix. draws from an mbg model in the cell pred form
#' @param rake_targets data.table with at least 3 columns (specified by lyv) with gaul code, year and val to rake towards
#' @param lyv character vector of length 3. Denotes the gaul_code column, year column and value column-- in that order
#' @param year_list numeric vector. Vector of years implied by the cell pred.
#' @param simple_raster raster. Should be the rasterized version of a shapefile or otherwise denoting the gaul/admin codes specified in rake targets
#' @param weight_brick rasterbrick. Raster brick where the values are the spatial weightings (usually population).
#' @param rake_method character. One of 'linear' or 'logit'.
#' @param fractional_pixel logical. Should pixels be assigned fractionally?
#' @usage a = calculate_raking_factors(cell_pred, aaa[[1]], year_list = 2000:2015, simple_raster = simple_raster, weight_brick = pop_raster_annual, rake_method = 'linear')
#' @usage b = calculate_raking_factors(cell_pred, aaa[[1]], year_list = 2000:2015, simple_raster = simple_raster, weight_brick = pop_raster_annual, rake_method = 'logit')
#' @import data.table
calculate_raking_factors = function(cell_pred, rake_targets, lyv = c('name','year', 'mean'), year_list, 
                                    simple_raster, weight_brick, rake_method = c('linear', 'logit'), fractional_pixel = F, MaxJump = 10, MaxIter = 80, FunTol = 1e-5, cores = 1){
  #check to make sure rake targets is a data table
  setDT(rake_targets)
  
  #check to make sure rake targets has valid columns
  if(!all(lyv %in% names(rake_targets))){
    stop('rake_targets does not contain all the columns specified by lyv')
  }
  
  #scoping
  rake_targets = copy(rake_targets[, lyv, with = F])
  setnames(rake_targets, lyv, c('loc', 'year', 'target'))
  
  
  #check to make sure all country years are represented
  sr_locs = unique(simple_raster[])
  sr_locs = sr_locs[!is.na(sr_locs)]
  cys = setDT(expand.grid(loc = sr_locs, year = year_list))
  start_row = nrow(cys)
  cys = merge(cys, rake_targets, by= c('loc','year'), all.x= T)
  
  if(start_row != nrow(cys)){
    stop('The location years implied by simple_raster and year_list are not matched in rake_targets')
  }
  
  #check to make sure weight_brick has the same number of years as year list
  if(dim(weight_brick)[3] != length(year_list)){
    stop('year_list implies a different number of time steps than the weight brick')
  }
  
  #check to make sure simple raster, weight_brick, and cell pred all have the proper dimensions
  stopifnot(dim(simple_raster)[1:2] == dim(weight_brick)[1:2])
  
  #check to make sure cell pred is an accurate subset of simple raster and weight brick
  if(!dim(cell_pred)[1] / length(cellIdx(simple_raster)) == dim(weight_brick)[3]){
    stop('cell_pred and simple_raster dimensions are not aligned')
  }
  
  #match rake_method
  rake_method = match.arg(rake_method)
  
  ##assuming checks pass, calculate a data table from which raking factors can be calculated
  sri = cellIdx(simple_raster)
  ngoodpixels = length(sri)
  
  #figure out which pixels over time have good values
  sri_years = unlist(lapply(seq(year_list)-1, function(x) sri + (x * ncell(simple_raster))))
  nyears = length(year_list)
  
  #first format the weight_brick
  raker = data.table(pixel_means = rowMeans(cell_pred),
                     pixel_xyt_id = sri_years,
                     pixel_xy_id  = sri,
                     loc = rep.int(simple_raster[sri], nyears),
                     year = as.vector(outer(Y = year_list, X = rep.int(1, ngoodpixels))),
                     weight = weight_brick[][sri_years])
  raker[,cell_pred_id := .I]
  
  #specify by_vars
  byvars = c('loc', 'year')
  
  #remove NA weight values from raker
  pre = dim(raker)[1]
  raker = raker[!is.na(weight) & !is.na(pixel_means), ]
  post = dim(raker)[1]
  
  if(pre != post){
    warning(paste(pre - post, 'pixels (over the whole cube) were removed because of NA weight_brick or pixel values'))
  }
  
  if(rake_method=='linear'){
    #collapse by country year to get the linear adjustment
    rak = raker[, list(px_mean = sum(pixel_means * weight, na.rm = T), sumweight = sum(weight, na.rm = T)), by = byvars]
    rak[, start_point := px_mean/sumweight]
    #merge on the target
    rak = merge(rak, rake_targets, all.x = T)
    rak[, raking_factor := target/start_point]
    
  } else if(rake_method == 'logit'){
    
    #for each country year, find the logit raking factor
    #redefine cys
    cys = unique(raker[,.(loc, year)])
    rf = lapply(seq(nrow(cys)), function(x) {
      #year and location
      theloc = cys[x,loc]
      theyear = cys[x,year]
      ret = NewFindK(gbdval    = rake_targets[loc == theloc & year == theyear,.(target)],
                     pixelval  = cell_pred[raker[loc == theloc & year == theyear, cell_pred_id],], #pass the cell pred rows that corrospond to this country year
                     weightval = raker[loc == theloc & year == theyear, weight],
                     MaxJump   = MaxJump,
                     MaxIter   = MaxIter,
                     FunTol    = FunTol)
      
      return(ret)
      
    })
    
    cys[,raking_factor := unlist(rf)]
    cys = merge(cys, rake_targets, all.x = T)
    
    #calculate the starting point post hoc for standardized reporting
    rak = raker[, list(px_mean = sum(pixel_means * weight, na.rm = T), sumweight = sum(weight, na.rm = T)), by = byvars]
    rak[, start_point := px_mean/sumweight]
    
    rak = merge(rak, cys, all.x = T, by = c('loc','year'))
  }
  
  #raking factors at the cys level
  rak = rak[, .(loc, year, start_point, target, raking_factor)]
  
  if(fractional_pixel){
    px_rak = merge(raker, rak, by = c('loc','year'), all.x = T)
  }else{
    px_rak = NA
  }
  
  return(list(rak, px_rak))
  
}


NewFindK <- function(gbdval, pixelval, weightval, MaxIter = 40, MaxJump = 10, FunTol = 1e-5){
  
  NumIter <- ceiling(-log2(FunTol / MaxJump))
  
  if(NumIter > MaxIter){
    stop(paste("Maximum number of iterations is less than projected iterations required:", NumIter / MaxIter))
  }
  
  CurrentError <- EvalDiff(gbdval, pixelval, weightval)
  if (CurrentError > 0){
    Range <- c(0, MaxJump)
  } else {
    Range <- c(-MaxJump, 0)
  }
  
  a <- Range[1]
  b <- Range[2]
  F_a <- EvalDiff(gbdval, NewPixelVal(a, pixelval), weightval)
  F_b <- EvalDiff(gbdval, NewPixelVal(b, pixelval), weightval)
  
  if (F_a * F_b > 0){
    stop("Your estimates are WAY too far away from GBD")
  } else {
    i <- 1
    c <- (a + b) / 2
    F_c <- EvalDiff(gbdval, NewPixelVal(c, pixelval), weightval)
    Success <- (abs(F_c) <= FunTol) 
    while (!Success & i < NumIter){
      if (sign(F_c) == sign(F_a)){
        a <- c
        F_a <- F_c
      } else {
        b <- c
        F_b <- F_c
      }
      c <- (a + b) / 2
      F_c <- EvalDiff(gbdval, NewPixelVal(c, pixelval), weightval)
      Success <- (abs(F_c) <= FunTol)
      i <- i + 1
    }
    if (Success){
      return(c)
    } else {
      return(sprintf("Need more iterations, current output: K = %g, F(K) = %g",c, F_c))
    }
  }
}

NewPixelVal <- function(K, vals){
  return(ilogit(logit(vals)+K))
}

#vals: matrix
#weightval: vector
NewEst <- function(vals, weightval){
  
  vals = vals * weightval
  vals = apply(vals, 2, sum)
  vals = vals / sum(weightval)
  
  return(mean(vals))
}

EvalDiff <- function(gbdval, vals, weightval){
  return(gbdval - NewEst(vals, weightval))
}

logit <- function(x) {
  log(x/(1-x))
}
ilogit <- function(x) {
  exp(x)/(1+exp(x))
}

#' A function to apply raking factors to a cell pred such that rowMeans(cell_pred) summed by country/year equal GBD
#' 
#' @param cell_pred Matrix. Cell_pred object
#' @param simple_raster Raster. Should be the rasterized version of a shapefile or otherwise denoting the gaul/admin codes specified in rake_dt
#' @param rake_dt data.table Data.table output of calculate_raking_factors. Or at the very least, a three column data.table with columns for loc, year and raking_factor.
#' @param logit_rake logical. Determines whether logit raking or linear scaling occurs
#' @param fractional_pixel logical. Determines whether the raking factor should be applied via the fractional pixel approach
#' @param force_simp_ras_dt_match logical. Determines whether the function should break if it detects a difference between number of country years models and number of country years
#'                                         supplied to rake.
apply_raking_factors = function(cell_pred, simple_raster, rake_dt, logit_rake = F, fractional_pixel = T, force_simp_ras_dt_match = T){
  
  cpdim = dim(cell_pred)
  
  #check to make sure simple raster and cell_pred play nice
  nyears = dim(cell_pred)[1] %% length(cellIdx(simple_raster)) != 0
  if(nyears){
    stop('good cells in simple raster is not a multiple of cell pred rows')
  }
  nyears = dim(cell_pred)[1] / length(cellIdx(simple_raster))
  thelocs = unique(as.vector(simple_raster))
  thelocs = thelocs[!is.na(thelocs)]
  
  if(force_simp_ras_dt_match){
    if(fractional_pixel == F & nrow(rake_dt) != (length(thelocs) * nyears)){
      stop('Combination of inputs suggest a disconnect between rake_dt and cell pred')
    }
  }
  if(!fractional_pixel){
    #create comparison to cell pred to fill in the raking factors
    dt = rbindlist(lapply(unique(rake_dt[,year]), function(x) data.table(cell_xy_id = cellIdx(simple_raster),
                                                                         loc =simple_raster[][cellIdx(simple_raster)],
                                                                         year = x)))
    
    dt[, id:= .I]
    
    #merge on the raking factors
    rake_dt = merge(dt, rake_dt, all.x = T, by = c('loc','year'))
    
    #enfore prespecified order
    setorder(rake_dt, id)
  }else{
    #find missing pixels
    mis_cells = setdiff(1:nrow(cell_pred), rake_dt[,cell_pred_id])
    
    if(length(mis_cells)>0){
      miss = data.table(cell_pred_id = mis_cells, raking_factor = NA, area_fraction = 0)
      rake_dt = rbind(rake_dt[,.(cell_pred_id, raking_factor, area_fraction)], miss)
    }
    setorder(rake_dt, cell_pred_id)
    
    #This no longer works if fractional pixel raking is being implemented
    stopifnot(nrow(rake_dt) == nrow(cell_pred))
    
  }
  
  if(!logit_rake){
    cell_pred = cell_pred * rake_dt[,raking_factor]
  }else{
    cell_pred = invlogit(logit(cell_pred) + rake_dt[,raking_factor])
  }
  
  return(cell_pred)
  
}

#overwrite cellIdx function to reduce outward dependancies
cellIdx =function(x) which(!is.na(getValues(x[[1]])))