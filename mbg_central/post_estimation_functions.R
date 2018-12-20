#################################################################################
### Pull GBD estimates from database and return as a data table
## Inputs:
##    gbd_type = "covariate" or "output", depending which database you need to pull from
##    gbd_name = GBD cov_name_short if "covariate" and GBD cause_id if "output"
##    gaul_list = list of GAUL codes you want to pull
##    measure_id = if "output", which measure you want to pull (defaults to incidence)
##    age_group_id = if "output", which age group you want to pull (defaults to Under 5)
##    metric_id = if "output", which metric you want to pull (defaults to rate)
## Outputs:
##    Returns 3-column data.table where "name" = GAUL code, "year" = year, and
##      "mean" = value. Consistent with expected input for calc_raking_factors.
#################################################################################
load_gbd_data     <- function(gbd_type,
                              gbd_name,
                              gaul_list,
                              measure_id = 6,
                              age_group_id = 1,
                              metric_id = 3,
                              year_ids = c(2000,2005,2010,2015),
                              return_by_age_sex = "no",
                              collapse_age_sex = FALSE,
                              gbd_round_id = 4){

  gaul_to_loc_id <- fread("<<<< FILEPATH REDACTED >>>>>/gaul_to_loc_id.csv")
  names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"

  if(gbd_type=="covariate") {
    
    library('dbplyr', lib.loc = '<<<< FILEPATH REDACTED >>>>>/mbg_pkgs_conda')
    metadata <- get_covariate_metadata()
    
    source('<<<< FILEPATH REDACTED >>>>>/get_covariate_estimates.R')
    gbd_estimates <- get_covariate_estimates(covariate_id = metadata[covariate_name_short == gbd_name, covariate_id], gbd_round_id = gbd_round_id)

    # Collapse to all age, all sex if needed
    if(collapse_age_sex==TRUE & length(unique(gbd_estimates[, age_group_id]))!=1) {
      gbd_estimates <- get_covariate_estimates(covariate_name_short = gbd_name, age_group_id = age_group_id)
      source('<<<< FILEPATH REDACTED >>>>>/get_population.R')
      gbd_pops <- get_population(age_group_id=paste(unique(gbd_estimates[, age_group_id]), collapse = ' '),
                                 location_id=paste(unique(gbd_estimates[, location_id]), collapse = ' '),
                                 year_id=paste(unique(gbd_estimates[, year_id]), collapse = ' '),
                                 sex_id=paste(unique(gbd_estimates[, sex_id]), collapse = ' '))
      gbd_estimates <- merge(gbd_estimates, gbd_pops, by=c('location_id','sex_id','age_group_id','year_id'))
      gbd_estimates <- gbd_estimates[, lapply(.SD, weighted.mean, w=population, na.rm=TRUE), by=c('location_id', 'year_id'), .SDcols='mean_value' ]
    }

    gaul_to_loc_id <- fread('<<<< FILEPATH REDACTED >>>>>/gaul_to_loc_id.csv')
    names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"
    gbd_estimates <- gbd_estimates[year_id %in% year_ids,]
    gbd_estimates <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
    gbd_estimates <- gbd_estimates[GAUL_CODE %in% gaul_list,]
    names(gbd_estimates)[names(gbd_estimates)=="GAUL_CODE"] <- "name"
    names(gbd_estimates)[names(gbd_estimates)=="year_id"] <- "year"
    names(gbd_estimates)[names(gbd_estimates)=="mean_value"] <- "mean"
    if(return_by_age_sex=='no') gbd_estimates <- gbd_estimates[, c('name', 'year', 'mean'), with = FALSE]
    if(return_by_age_sex=='yes') gbd_estimates <- gbd_estimates[, c('name', 'year', 'mean', 'sex_id', 'age_group_id'), with = FALSE]

    return(gbd_estimates)

  }

  if(gbd_type=="output") {

    loc_ids <- gaul_to_loc_id[GAUL_CODE %in% gaul_list]
    loc_ids <- loc_ids$location_id
    loc_ids <- toString(loc_ids)
    loc_ids <- gsub(",","",loc_ids)

    source('<<<< FILEPATH REDACTED >>>>>/get_outputs.R')
    gbd_estimates <- get_outputs(topic = "cause",
                                 version = "best",
                                 gbd_round_id = gbd_round_id,
                                 cause_id = gbd_name,
                                 measure_id = measure_id,
                                 metric_id = metric_id,
                                 age_group_id = age_group_id,
                                 location_id = loc_ids,
                                 year_id = year_ids)
    all_data <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
    names(all_data)[names(all_data)=="GAUL_CODE"] <- "name"
    names(all_data)[names(all_data)=="year_id"] <- "year"
    names(all_data)[names(all_data)=="val"] <- "mean"
    all_data <- as.data.table(all_data)
    all_data <- all_data[, c('name', 'year', 'mean'), with = FALSE]

    return(all_data)

  }

}




#################################################################################
### Pull population data and return as a raster
## Inputs:
    # simple_raster:
## Outputs:
## Notes:
    # Population data in a temporary location, Lucas to find permanent home for it
    # Later will need to be explicit about years (for now everyone is doing 2000-2015, 4 yr intervals)
#################################################################################
get_population_data     <- function(simple_raster){

  # load population raster
  # RB 27SEPT: Note this is a temporary location, and is only Africa so updates will be necessary
  pop<- brick(sprintf('<<<< FILEPATH REDACTED >>>>>/pop_stack.tif',root))

  # make sure population is cropped and extented to same as simple_raster
  # this is important, otherwise IDX wont work.
  if(!is.null(simple_raster)){
    pop <- mask(crop(pop,simple_raster),simple_raster)
    extent(pop)=extent(simple_raster)
  }
  return(pop)
}




#################################################################################
### load admin raster
## Inputs:
    # admin_level: 0,1, or 2 are accepted.
    # disag_fact: factor to increase resolution (needed whith small districts), 50 makes it 100m
    # simple_raster: keeps resolution and extent
## Outputs: returns raster layer with gaul codes for values
## Notes:
    # Will need to update location once lucas gives these a permanent home
#################################################################################
load_admin_raster  <- function(admin_level, simple_raster, disag_fact=NULL){

    if(!admin_level %in% c(0,1,2)) stop("admin_level must be either 0, 1, or 2")

    # load admin raster
    #adm <- raster(sprintf('%stemp/geospatial/U5M_africa/data/clean/shapefiles/ad%i_raster.grd',root,admin_level))
    if(!is.null(disag_fact)){
      sr = disaggregate(simple_raster,fact=disag_fact)
    } else {
      sr = simple_raster
    }

  # UPDATED: master gaul admin shapefiles
    if(admin_level %in% c(0,1)) shapes <- shapefile(paste0("<<<< FILEPATH REDACTED >>>>>", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, ".shp"))
    if(admin_level == 2) shapes <- shapefile(paste0("<<<< FILEPATH REDACTED >>>>>", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, "_modified.shp"))

  # fix several african disputed territories (assign them to a country so we can include them in our estimates and they dont drop out)
    #  61013     Ilemi triangle   to Kenya (gaul 133)
    #  40760   Hala'ib triangle   to Egypt (gaul 40765)
    #  40762    Ma'tan al-Sarra   to Libya (gaul 145)
    #  102              Abyei     to Sudan (gaul 6)
    shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==61013]=133
    shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40760]=40765
    shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==40762]=145
    shapes[[paste0('ADM', admin_level,'_CODE')]][shapes[[paste0('ADM', admin_level,'_CODE')]]==102  ]=6

  # crop
    cropped_shapes <- crop(shapes, extent(sr), snap="out")

    ## Fix rasterize
    initial_raster <- rasterize(cropped_shapes, sr, field = paste0('ADM', admin_level,'_CODE'))
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])!=0) {
      rasterized_shape <- merge(rasterize(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),], sr, field = paste0('ADM', admin_level,'_CODE')), initial_raster)
    }
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])==0) {
      rasterized_shape <- initial_raster
    }
    #rasterized_shape <- rasterize(cropped_shapes, sr, field=paste0('ADM', admin_level,'_CODE'))
    masked_shapes <- mask(x=rasterized_shape, mask=sr)

    return(masked_shapes)

}




#################################################################################
### Use admin and population rasters and returns a weights raster
## Inputs:
    # admin_level: 0,1, or 2 are accepted. Level to which weights aggregate to 1
    # simple_raster: this is the raster template, should be the one that cell_preds was made with.
    # countries: should typically be template, by just a list of countries you want included in the output
## Outputs: returns a list with the following:
    # pop_wts: a cell_idx by years matrix of population weights, each number representing a cell
    # admin_level: 0,1,or 2 which admin level this represents
    # adm_code: which admin code (gaul) each cell is coded to
## Notes:
    # Later will need to be explicit about years (for now everyone is doing 2000-2015, 4 yr intervals)
#################################################################################
make_population_weights  <- function(admin_level,
                                     simple_raster,
                                     pop_raster,
                                     gaul_list) {

  # load admin raster, crop down to simple_raster
  adm <- load_admin_raster(admin_level, simple_raster) 
  adm <- crop(adm,simple_raster)
  adm <- mask(adm,simple_raster)

  if(admin_level==0) adm[adm==1013965]=227

  # load population data
  pop <- pop_raster

  # do using cell index (vectorized)
  cell_idx <- seegSDM:::notMissingIdx(simple_raster)
  adm_cell <- extract(adm, cell_idx)

  pop_layer_names <- names(pop_raster)

  make_weights_vector <- function(x) {

    ad2_code <- extract(adm$layer, cell_idx)
    ad2_code[is.na(ad2_code)] <- 0
    pop_cell <- extract(pop[[x]], cell_idx)
    pop_cell[is.na(pop_cell)] <- 0

    message(paste0(length(unique(ad2_code)), " unique admin codes in this level."))

    pop_totals_adm <- aggregate(pop_cell~ad2_code, FUN=sum)
    pop_totals_adm$pop_cell[pop_totals_adm$pop_cell==0] <- 0.0001

    ## replicate totals for all cells in area
    pop_totals_adm_cell <- as.vector(pop_totals_adm)[match(adm_cell,
                                                           pop_totals_adm$ad2_code),]
    pop_totals_adm_cell$ad2_code <- NULL

    ## get  population weights for each cell
    pop_wt_adm <- pop_cell / pop_totals_adm_cell
    pop_wt_adm[is.na(pop_wt_adm)] <- 0
    pop_wt_adm_vector <- pop_wt_adm$pop_cell

    wt_sum_ad2 <- tapply(pop_wt_adm_vector, ad2_code, sum)
    stopifnot(all.equal(wt_sum_ad2, round(wt_sum_ad2)))

    names(pop_wt_adm_vector) <- x
    return(pop_wt_adm_vector)

  }

  pop_weights_by_layer <- lapply(pop_layer_names, make_weights_vector)
  pop_wt_matrix <- do.call(cbind, pop_weights_by_layer)
  rownames(pop_wt_matrix) <- cell_idx
  colnames(pop_wt_matrix) <- pop_layer_names

  # return
  return(list('pop_wt'     =as.matrix(pop_wt_matrix),
              'admin_level'=admin_level,
              'adm_code'   =adm_cell))

}



#################################################################################
### Wrapper for condSim from seegMBG
## Inputs:
    # pop_wts_object: list object that is output of make_population_weights()
    # cell_pred: Cells by Draws matrix which is output from predict_mbg()
    # years: vector of years in your analysis (should match population as well)
    # fun: Mean by default. Other functions include gini
## Outputs: returns a matrixs of (number of admin units * number of years in analysis)
##          rows by #draws in cell_pred columns
#################################################################################
make_condSim  <- function(pop_wts_object,
                          cell_pred,
                          admin_level,
                          years     = c(2000,2005,2010,2015),
                          gaul_list = gaul_list,
                          fun       = 'mean',
                          summarize = FALSE){
  # ease
  p=pop_wts_object

  message(sprintf('Aggregating results by draw to the admin %i level. If you want a different aggregation level, re-run make_population_weights with a different value for the admin_level argument. \n\n',p$admin_level))

  # make pop_wt long
  pop_wt=as.vector(p$pop_wt)

  # throw exception if data are non-conformable
  if(length(pop_wt)!=dim(cell_pred)[1])
    stop('Population Weights not same number of cells as cell predictions.')


  # get admin gauls (or names) at each cell. Repeat for number of years
  yrs = length(pop_wt)/length(p$adm_code)
  message(sprintf('There are %i years of data here. Modeler please be sure your population years and data years are aligned. Current input is %s \n\n',yrs,paste(years,collapse=' ')))

  p$adm_code[p$adm_code == 'NA'] <- NA
  periods <- rep(years, each = length(p$adm_code))
  adm_code_all <- paste(p$adm_code, periods, sep = '_')

  # get NA mask
  good_cells <- which(!is.na(rep(p$adm_code,yrs)))
  good_cells <- good_cells[which(!is.na(cell_pred[good_cells,1]))]
  good_cells <- good_cells[which(!is.na(pop_wt[good_cells]))]

  # rescale pop weights if some got dropped
  require(data.table)
  dt = data.table(pw=pop_wt[good_cells],ad=adm_code_all[good_cells])
  dt$pw[dt$pw==0]=0.00000001 # so no areas have zero pop which leads to buggy behavior in 2 or 3 districts
  dt[,tot:=sum(pw),by=ad]
  dt[,pw:=pw/tot]
  pw = dt$pw
  ad = dt$ad

  # actual condSim run
  if(fun=='mean'){
    cond_sim_adm <- condSim(vals    = cell_pred[good_cells,],
                            weights = pw, #pop_wt[good_cells],
                            group   = ad) #adm_code_all[good_cells])
  } else {
    cond_sim_adm <- condSim(vals    = cell_pred[good_cells,],
                            weights = pw, #pop_wt[good_cells],
                            group   = ad, #adm_code_all[good_cells],
                            fun     = fun)
  }

  # summarize if requested
  if(summarize==TRUE){
    message(sprintf('Returning Summary (mean) made from %i draws. \n\n',dim(cell_pred)[2]))
    cond_sim_adm=(apply(cond_sim_adm   , 1, mean))
  } else{
    message(sprintf('Returning %i draws. \n\n',dim(cell_pred)[2]))
  }

  # Make sure we only return results for admin units in the list of adm0 gaul codes requested. 
  # Get list of admin* codes within selected admin0 codes
  adm_gaul_list <- get_gaul_codes_subnat(gaul_list, admin_level)
  # Convert named vector to dataframe
  if(summarize) {
    tmp <- cbind(read.table(text = names(cond_sim_adm)), val = cond_sim_adm)
    # Make new column of just gaul codes
    tmp$gaul_code <- gsub("_.*", "", tmp$V1)
    # Subset
    tmp <- tmp[tmp$gaul_code %in% adm_gaul_list,]
    # Convert back to dataframe
    clean_cond_sim_adm <- tmp$val
    names(clean_cond_sim_adm) <- tmp$V1

    return(clean_cond_sim_adm)

  } else { return(cond_sim_adm) }


}



#################################################################################
### Split country/years (in format 'country_year') out of rownames of a condSim df
## Inputs:
    # condSim_object: data matrix that make_condSim() spits out
## Outputs: returns matrix with name and year in two columns of same length
#################################################################################
split_geo_names <- function(condSim_object){
  splits <- strsplit(gsub("C\xf4te d'Ivoire","Cote dIvoire" ,rownames(condSim_object)), split = '_')
  ctry   <- sapply(splits, '[', 1)
  year   <- as.numeric(sapply(splits, '[', 2))
  ret    <- cbind('name' = ctry,
                  'year' = as.numeric(year))
  rownames(ret) <- NULL
  return (ret)
}


#################################################################################
### Takes Geo condSim data and merges estimates to others (typically GBD) and get raking factors
## Inputs:
    # agg_geo_est: named vector from make_condSim(...,summarize=TRUE)
    # rake_to:     sometime to rake to, typically from load_gbd_data()
## Outputs: data table with merged mean estimates and raking factors
#################################################################################
calc_raking_factors <- function(agg_geo_est = cond_sim_adm0,
                                gaul_list   = gaul_list,
                                rake_to     = gbd) {

  message('WARNING: function will not work as expected if agg_geo_est and rake_to are not aggregated at the same admin level.')

  if(sum(colnames(rake_to) %in% c('name','year','mean')) != 3)
      stop('rake_to should be a data table with column names `name`, `year`, and `mean`')

  if(dim(t(as.matrix(agg_geo_est)))[1]!=1 | is.null(names(agg_geo_est)))
    stop('agg_geo_est must be a named vector returned from make_condSim() with summarize option equal to TRUE')



  # transpose and rename agg_geo_est
  agg_geo_est <- data.table(split_geo_names(as.matrix(agg_geo_est)),agg_geo_est)
  agg_geo_est$year <- as.numeric(agg_geo_est$year)
  agg_geo_est$name <- as.numeric(agg_geo_est$name)

  # merge
  merged <- merge(agg_geo_est,rake_to, by=c('name','year'),all.x=T)
  names(merged)[names(merged)=='agg_geo_est'] <- 'geo_mean'
  names(merged)[names(merged)=='mean'] <- 'rake_to_mean'

  # calculate raking factors
  merged[,raking_factor:=rake_to_mean/geo_mean,]

  return(merged)
}

#################################################################################
### Rakes estimates to gold standard estimates (usually GBD)
## Inputs:
    # raking_factors: output from calc_raking_factors
    # pop_wts_object: output from make_population_weights
    # cell_pred: Cells by Draws matrix which is output from predict_mbg()
## Outputs: returns cell_preds raked
#################################################################################
rake_predictions <- function(raking_factors,
                             pop_wts_object,
                             cell_pred,
                             logit_rake=FALSE){

  ## Replicate adm codes by years
  adm_codes_all = pop_wts_object$adm_code
  adm_codes_all[adm_codes_all == 'NA'] <- NA
  periods <- rep(unique(raking_factors$year), each = length(adm_codes_all))
  adm_codes_all <- paste(adm_codes_all, periods, sep = '_')

  # merge adm codes and raking factor
  z=data.frame(ID=adm_codes_all,order=1:length(adm_codes_all))
  z$ID<-as.character(z$ID)
  raking_factors$ID = paste(raking_factors$name,raking_factors$year,sep='_')
  factors_all=merge(z,raking_factors[,c('ID','raking_factor'),with=FALSE],by="ID",all.x=T)
  factors_all=factors_all[order(factors_all$order), ]

  # make a results frame
  raked<-matrix(NA,nrow=nrow(cell_pred),ncol=ncol(cell_pred))

  if(dim(factors_all)[1]!=dim(raked)[1])
    stop('Factors_all and Raked dimensions do not match')

  # multiply through the raking factor
  message(sprintf('raking %i draws',dim(raked)[2]))
  pb <- txtProgressBar(min=0,max=dim(raked)[2],initial=0)
  if(logit_rake==FALSE) {
    for(i in 1:dim(raked)[2]){
          setTxtProgressBar(pb,i)
          raked[,i] <- cell_pred[,i]*factors_all$raking_factor
    }
  }
  if(logit_rake==TRUE) {
    for(i in 1:dim(raked)[2]){
      setTxtProgressBar(pb,i)
      raked[,i] <- ilogit(logit(cell_pred[,i]) + factors_all$raking_factor)
    }
  }
  close(pb)

  return(raked)

}


#################################################################################
### CI Range estimate given draws
#################################################################################
cirange = function(x){
    z=quantile(x,probs=c(.025,.975),na.rm=T)
    return(z[2]-z[1])
}



#################################################################################
### Takes in raked or raw draw-level estimates and makes stat summary rasters
## Inputs:
    # draw_level_cell_pred: Cells by Draws matrix which is output from predict_mbg() or from rake_predictions()
    # mask: Should be the simple_raster
    # return_as_raster: If TRUE returns as raster, else as table
    # summary_stat: ie mean, cirange, quantile, sd
## Outputs: Summary table or raster of the cell_pred table put in
#################################################################################
make_cell_pred_summary    <- function(draw_level_cell_pred,
                                      mask                 = simple_raster,
                                      return_as_raster     = TRUE,
                                      summary_stat         = 'mean',
                                      ...){

    # custom summary functions here
    lower <- function(x) {
      # Simply get and return a percentile
      output <- quantile(x, 2.5 / 100, na.rm = T)
      return(output)
    }
    
    upper <- function(x) {
      # Simply get and return a percentile
      output <- quantile(x, 97.5 / 100, na.rm = T)
      return(output)
    }
  
    qcd <- function(x, quantile_low = 5, quantile_high = 95) {

      # "quantile  coefficient of dispersion"
      low <- quantile_low/100
      high <- quantile_high/100
      quantiles <- quantile(x, c(low, high), na.rm = T)

      q_low <- as.numeric(quantiles[1])
      q_high <- as.numeric(quantiles[2])

      output <- (q_high - q_low) / (q_high + q_low)

      return(output)
    }

    iqr <- function(x, quantile_low = 25, quantile_high = 75) {

      # "interquantile" range - can specify
      low <- quantile_low / 100
      high <- quantile_high / 100
      quantiles <- quantile(x, c(low, high), na.rm = T)

      q_low <- as.numeric(quantiles[1])
      q_high <- as.numeric(quantiles[2])

      output <- (q_high - q_low)

      return(output)

    }

    iqr_traditional <- function(x, quantile_low = 25, quantile_high = 75) {

      # "interquantile" range - can specify
      low <- quantile_low / 100
      high <- quantile_high / 100
      quantiles <- quantile(x, c(low, high), na.rm = T)

      q_low <- as.numeric(quantiles[1])
      q_high <- as.numeric(quantiles[2])

      output <- (q_high - q_low)

      return(output)

    }

    get_percentile <- function(x, percentile) {

      # Simply get and return a percentile
      percentile <- as.numeric(percentile)
      output <- quantile(x, percentile / 100, na.rm = T)
      return(output)

    }

    p_below <- function(x, value, equal_to = F) {

      # probability <  (or <= if equal_to = T) target value
      value <- as.numeric(value)
      if (equal_to == T) output <- sum(x <= value)
      if (equal_to == F) output <- sum(x < value)

      output <- output / length(x)
      return(output)

    }

    p_above <- function(x, value, equal_to = F) {

     # probability > (or >= if equal_to = T) target value
      value <- as.numeric(value)
      if (equal_to == T) output <- sum(x >= value)
      if (equal_to == F) output <- sum(x > value)

      output <- output / length(x)
      return(output)

    }

     # make summary
     summ <- apply(draw_level_cell_pred, 1, summary_stat, ...)

     # put it in a raster
     if(return_as_raster){
       yrs = dim(draw_level_cell_pred)[1]/length(cellIdx(mask))
       message(sprintf('Making a RasterBrick with %i layers',yrs))
       summ <- insertRaster(mask,  matrix(summ,  ncol = yrs))
     }


     return(summ)

}


## make_admin_pred_summary ###################################################

#' Take an admin pred object and an sp_hierarchy list and generate a sorted,
#' cleaned table for a given set of summary statistics
#'
#' @param admin_pred the admin draws object (e.g. admin_1, etc)
#' @param sp_hierarchy_list spatial hierarchy object from the admin draws RData files
#' @param summary_stat summary statistic to apply
#' @param ... any other arguments to pass to the `summary_stats` functions
#'            (note that currently will be passed to all functions)
#' @return data table of summary stats and admin info
#' @examples
#' # load("indicator_raked_admin_draws_eb_bin0.RData") # Replace with your admin draws
#' make_admin_pred_summary(admin_2, 
#'                         sp_hierarchy_list, 
#'                         c("mean", "cirange", "upper", "lower"))

make_admin_pred_summary    <- function(admin_pred,
                                       sp_hierarchy_list, 
                                       summary_stats = 'mean',
                                       ...){
  
  ### Get set up  
  str_match <- stringr::str_match
  
  # Split up your input data
  ad_code <- subset(admin_pred, select = grep("ADM[0-9]_CODE", names(admin_pred)))
  year <- subset(admin_pred, select = "year")
  pop <- subset(admin_pred, select = "pop")
  draws <- subset(admin_pred, select = grep("V[0-9]*", names(admin_pred)))
  
  # Get the admin level
  ad_code_name <- names(ad_code)
  ad_level <- as.numeric(str_match(ad_code_name,"ADM([0-9])_CODE")[,2])
  
  # Get all admin levels
  all_admin_levels <- as.numeric(str_match(names(sp_hierarchy_list), "ADM([0-9])_CODE")[,2])
  all_admin_levels <- unique(all_admin_levels)[!is.na(unique(all_admin_levels))]
  
  ### Make summary
  summ_list <- lapply(summary_stats, function(ss) {
    apply(draws, 1, ss, ...)})
  
  if (length(summ_list) > 1) {
    output_df <- as.data.table(cbind(year, ad_code, do.call(cbind, summ_list)))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  } else if (length(summ_list) == 1) {
    output_df <- as.data.table(cbind(year, ad_code, summ_list[[1]]))
    names(output_df)[grep("V[0-9]*", names(output_df))] <- summary_stats
  }
  
  ### Add on identifiers
  
  # Drop smaller admin levels
  drop_levels <- all_admin_levels[all_admin_levels > ad_level]
  
  if (length(drop_levels) > 0) {
    for (d in drop_levels) {
      drop_cols <- names(sp_hierarchy_list)[grepl(paste0("ADM",d), names(sp_hierarchy_list))]
      sp_hierarchy_list <- subset(sp_hierarchy_list, select = !(names(sp_hierarchy_list) %in% drop_cols))
    }
  }
  sp_hierarchy_list <- unique(sp_hierarchy_list)
  
  output_df <- merge(sp_hierarchy_list, output_df, all.y = T, all.x = F)
  
  # Clean up & sort
  order_vars <- c(paste0("ADM", 0:ad_level, "_NAME"), "year")
  setorderv(output_df, order_vars)
  
  # Get col order
  ad_col_order <- sort(names(output_df)[grep("AD*_*", names(output_df))])
  cols_to_order <- c(ad_col_order, "region", "year", summary_stats)
  if (length(cols_to_order) < length(names(output_df))) {
    other_cols <- names(output_df)[!(names(output_df) %in% cols_to_order)]
    cols_to_order <- c(cols_to_order, other_cols)
  }
  
  setcolorder(output_df, cols_to_order)
  
  return(output_df)
  
}

## summarize_admins ################################################

#' Function to summarize admin_pred objects
#'
#' This is a wrapper for `make_admin_pred_summary()`
#'
#' @param ind indicator
#' @param ig indicator_group
#' @param summstats Summary statistics (functions) to compute.
#'                  Order will be the order they are written in csv
#'                  This is passed to `make_admin_pred_summary()`
#' @param raked Raked (T), unraked (F), or both (`c(T,F)`)?
#' @param ad_levels Admin levels to summarize (0, 1, and 2 by default)
#' @return Writes csv files to `sharedir/pred_derivatives/admin_summaries/`
#' @examples
#' summarize_admins(summstats = c("mean", "lower", "upper", "cirange"), 
#'                  ad_levels = c(0,1,2), 
#'                  raked     = c(T,F))

summarize_admins <- function(ind = indicator,
                             ig = indicator_group,
                             summstats = c("mean", "lower", "upper", "cirange"),
                             raked = c(T,F),
                             ad_levels = c(0,1,2)) {

  sharedir       <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s',ig,ind)
  input_dir <- paste0(sharedir, "/output/", run_date, "/")
  output_dir <- paste0(input_dir, "/pred_derivatives/admin_summaries/")
  dir.create(output_dir, recursive = T, showWarnings = F)

  # Convert raked to character
  rr <- character()
  if (T %in% raked) rr <- c(rr, "raked")
  if (F %in% raked) rr <- c(rr, "unraked")

  # Summarize and save admin preds
  for (rake in rr) {
    load(paste0(input_dir, ind, "_", rake, "_admin_draws_eb_bin0_0.RData"))
    for (ad in ad_levels) {
      message(paste0("Summarizing ", ind, ": admin ", ad, " (", rake, ")"))
      ad_summary_table <- make_admin_pred_summary(admin_pred = get(paste0("admin_", ad)),
                                                  sp_hierarchy_list,
                                                  summary_stats = summstats)
      fwrite(ad_summary_table, 
             file = paste0(output_dir, ind, "_admin_", ad, "_", rake, "_summary.csv"))
    }
  }
}

#################################################################################
### Saves post-estimation output files in the proper directory
## Inputs:
## Outputs:
#################################################################################
save_post_est   <- function(x,
                            filetype,
                            filename,
                            indic = indicator){

   output_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indic, '/output/', run_date)


   dir.create(output_dir, showWarnings = FALSE)

   filetype = tolower(filetype)
   if(!filetype %in% c('rdata','raster','csv'))
      stop('filetype argument has to be either rdata or raster or csv')


   if(filetype=='raster')
       writeRaster(
          x,
          file = paste0(output_dir, '/', indic,'_',filename),
          format='GTiff',
          overwrite = TRUE
        )

   if(filetype=='rdata')
       save(
          x,
          file = paste0(output_dir, '/', indic,'_',filename,'.RData'),
          compress = TRUE
        )

   if(filetype=='csv')
       write.csv(
          x,
          file = paste0(output_dir, '/', indic,'_',filename,'.csv')
        )


}



#################################################################################
### Pull cell preds
## Inputs:
## Outputs:
#################################################################################

load_cell_preds <- function(indicator_group,
                            indicator,
                            rd = run_date,
                            region,
                            agebin,
                            u5m=FALSE,
                            other='',
                            ageasindic=TRUE){


if(u5m){
   if(ageasindic==FALSE){
     load(paste0('<<<< FILEPATH REDACTED >>>>>',indicator_group,'/',indicator,'/output/',rd,'/',
               indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'NA.RData')) # the 0 are no holdout
   } else {
     load(paste0('<<<< FILEPATH REDACTED >>>>>',indicator_group,'/',indicator,'_age',agebin,'/output/',rd,'/',
               indicator,'_age',agebin,'_cell_draws_eb_bin',agebin,'_',region,'_0',other,'.RData'))
   }
} else {
   load(paste0('<<<< FILEPATH REDACTED >>>>>',indicator_group,'/',indicator,'/output/',rd,'/',
               indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0.RData')) # the 0 are no holdouts
}
   return(cell_pred)
}



load_cell_preds_stack <- function(indicator_group,
                            indicator,
                            rd = run_date,
                            region,
                            agebin){

   load(paste0('<<<< FILEPATH REDACTED >>>>>',indicator_group,'/',indicator,'/output/',rd,'/',
               indicator,'_cell_draws_eb_bin',agebin,'_',region,'_0_stacked_results.RData')) # the 0 are no holdouts

   return(cell_pred)
}


fit_stats <- function(is_data,
                      oos_data,
                      indicator,
                      indicator_group,
                      run_date,
                      pathaddin) {

  ### Fit statistics
  model_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date)
  image_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/model_image_history')

  ## Append in- and out-of-sample coords
  is_data <- is_data[, c('latitude','longitude','period','N',indicator), with=FALSE]
  is_data <- is_data[, oos := 0]
  oos_data <- oos_data[, c('latitude','longitude','period', 'N',indicator), with=FALSE]
  oos_data <- oos_data[, oos := 1]
  all_data <- rbind(is_data, oos_data)
  if(indicator_family=='binomial') all_data <- all_data[, obs := get(indicator) / N]
  if(indicator_family!='binomial') all_data <- all_data[, obs := get(indicator)]

  ## Extract predictions at input data and bind to input datatable
  message("Loading cell_preds and inputs...")
  load(paste0(image_dir,'/', run_date, pathaddin, '.RData'))
  mean_preds <- brick(paste0(model_dir, '/', indicator, '_prediction_eb', pathaddin))
  load(paste0(model_dir, '/', indicator, '_cell_draws_eb', pathaddin, '.RData'))

  cell_pred.dt <- as.data.table(cell_pred)
  cols <- names(cell_pred.dt)
  message("Making upper and lower credible interval rasters...")
  cell_pred.dt <- cell_pred.dt[ , upper := apply(.SD, 1, quantile, c(.975), na.rm=TRUE), .SDcols=cols]
  cell_pred.dt <- cell_pred.dt[ , lower := apply(.SD, 1, quantile, c(.025), na.rm=TRUE), .SDcols=cols]
  upper <- cell_pred.dt[, upper]
  lower <- cell_pred.dt[, lower]

  upper_raster <- insertRaster(simple_raster,
                               matrix(upper,
                                      ncol = length(unique(all_data$period))))
  lower_raster <- insertRaster(simple_raster,
                               matrix(lower,
                                      ncol = length(unique(all_data$period))))

  message("Extracting upper and lower values at all coordinates...")
  # extract cluster covariates
  all_data$longitude<-as.numeric(all_data$longitude)
  all_data$latitude <-as.numeric(all_data$latitude )

  # add dummy column for time-varying covariates
  all_data$pred <- NA
  all_data$upper <- NA
  all_data$lower <- NA

  # loop through time varying predictions to insert them
  for (period in sort(unique(all_data$period))) {

    # find matching rows
    message(paste0('Extracting predictions for period ',period))
    idx_tv<- which(all_data$period == period)

    # get period prediction raster
    craster  <- mean_preds[[period]]
    crs(craster)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

    uraster  <- upper_raster[[period]]
    crs(uraster)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

    lraster  <- lower_raster[[period]]
    crs(lraster)="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

    # extract values
    extr <- extract(craster, all_data[idx_tv, c('longitude', 'latitude'),with=F])
    upper_extr <- extract(uraster, all_data[idx_tv, c('longitude', 'latitude'),with=F])
    lower_extr <- extract(lraster, all_data[idx_tv, c('longitude', 'latitude'),with=F])

    # add into df
    all_data[['pred']][idx_tv]=extr
    all_data[['upper']][idx_tv]=upper_extr
    all_data[['lower']][idx_tv]=lower_extr

  }

  ## Make fit statistics
  ## Add to df to combine later however we want
  message("Calculating all fit statistics at each coordinate...")

  # Drop missing predictions
  message(paste0('Predictions missing for ', length(all_data[is.na(pred), pred]), '/', length(all_data[ , pred]), ' observations'))
  message('If you have missings, this either due to missing covariate values (bad) or points that were in the buffer zone (fine).')
  message('Dropping missing rows...')
  df_no_nas <- all_data[!is.na(pred), ]

  # Coverage
  df_no_nas <- df_no_nas[obs >= lower & obs <= upper, covered := 1]
  df_no_nas <- df_no_nas[obs < lower | obs > upper, covered := 0]

  # Error
  df_no_nas <- df_no_nas[, error := obs - pred]

  # Make statistics over entire IS dataset
  mean_error <- mean(df_no_nas[oos==0, error])
  mean_absolute_error <- mean(df_no_nas[oos==0, abs(error)])
  rmse <- sqrt(mean(df_no_nas[oos==0, error]^2))
  coverage <- mean(df_no_nas[oos==0, covered])

  # Make statistics over entire OOS dataset
  oos_mean_error <- mean(df_no_nas[oos==1, error])
  oos_mean_absolute_error <- mean(df_no_nas[oos==1, abs(error)])
  oos_rmse <- sqrt(mean(df_no_nas[oos==1, error]^2))
  oos_coverage <- mean(df_no_nas[oos==1, covered])

  # Average over all observations
  message("IN-SAMPLE:")
  message(paste0('                       Average error:    ', round(mean_error, 3)))
  message(paste0('                       Average MAE:      ', round(mean_absolute_error, 3)))
  message(paste0('                       Average RMSE:     ', round(rmse, 3)))
  message(paste0('                       Average coverage: ', round(coverage, 3)))

  message("OUT-OF-SAMPLE:")
  message(paste0('                       Average error:    ', round(oos_mean_error, 3)))
  message(paste0('                       Average MAE:      ', round(oos_mean_absolute_error, 3)))
  message(paste0('                       Average RMSE:     ', round(oos_rmse, 3)))
  message(paste0('                       Average coverage: ', round(oos_coverage, 3)))

  fit_folder <- paste0('<<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/fit_stats')
  message(paste0('Saving fit statistics in ', fit_folder))
  dir.create(fit_folder, showWarnings = FALSE)
  message(pathaddin)
  write.csv(df_no_nas, file = paste0(fit_folder, '/fit_stats_', pathaddin, '.csv'))

}

fit_stats_ho_id <- function(draws,
                            draw_prefix,
                            observed) {

  message('Your draws must contain columns oos and ho_id.')
  message('Summarizing draws called phat_*...')
  cols <- grep(draw_prefix, names(draws), value=TRUE)
  draws <- draws[, mean := .(mean = rowMeans(.SD)), by = ho_id, .SDcols=cols]
  draws <- draws[ , upper := apply(.SD, 1, quantile, c(.975), na.rm=TRUE), .SDcols=cols]
  draws <- draws[ , lower := apply(.SD, 1, quantile, c(.025), na.rm=TRUE), .SDcols=cols]

  message('Calculating fit statistics at each ho_id...')
  # Coverage
  draws <- draws[get(observed) >= lower & get(observed) <= upper, covered := 1]
  draws <- draws[get(observed) < lower | get(observed) > upper, covered := 0]

  # Error
  draws <- draws[, error := get(observed) - mean]

  # Make statistics over entire OOS dataset
  oos_mean_error <- mean(draws[oos==1, error])
  oos_mean_absolute_error <- mean(draws[oos==1, abs(error)])
  oos_rmse <- sqrt(mean(draws[oos==1, error]^2))
  oos_coverage <- mean(draws[oos==1, covered])

  message("OUT-OF-SAMPLE:")
  message(paste0('                       Average error:    ', round(oos_mean_error, 3)))
  message(paste0('                       Average MAE:      ', round(oos_mean_absolute_error, 3)))
  message(paste0('                       Average RMSE:     ', round(oos_rmse, 3)))
  message(paste0('                       Average coverage: ', round(oos_coverage, 3)))

  message('Returning draws with fit stat columns added.')
  return(draws)

}

compile_results_table <- function(ind_gp,
                                  ind,
                                  rd,
                                  measure = 'mortality',
                                  baseline_year = 2000,
                                  year_to_map = 2015,
                                  goal_threshold = 0.001,
                                  metric = 'sdgprob',
                                  out.dir ## location to save csvs and rasters
                                  ) {

  message("loading and setting up")
  results_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', ind_gp, '/', ind, '/output/', rd, '/table_', baseline_year)
  dir.create(paste0(results_dir, '/summary_tables'))
  all_files <- list.files(results_dir, pattern = measure, full.names = TRUE)
  if(metric=='percent') all_files <- all_files[grepl('_percent2010', all_files)]
  if(metric=='sdgprob') all_files <- all_files[!grepl('_percent2010', all_files)]
  all_admins <- rbindlist(lapply(all_files, fread))

  message("combining admins")
  all_admins <- all_admins[is.na(mean), mean := goal_threshold]
  all_admins <- all_admins[is.na(upper), upper := goal_threshold]
  all_admins <- all_admins[is.na(lower), lower := goal_threshold]
  setnames(all_admins, 'admin2', 'ADM2_CODE')
  setnames(all_admins, 'admin1', 'ADM1_CODE')
  setnames(all_admins, 'admin0', 'ADM0_CODE')
  gaul_to_loc <- fread('<<<< FILEPATH REDACTED >>>>>/gaul_to_loc_id.csv')
  all_gauls <- unique(all_admins[, ADM0_CODE])
  for(admin_level in c(0,1,2)) {
    if(admin_level %in% c(0,1)) admin_names <- read.dbf(paste0("<<<< FILEPATH REDACTED >>>>>/", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, ".dbf"))
    if(admin_level == 2) admin_names <- read.dbf(paste0("<<<< FILEPATH REDACTED >>>>>", admin_level, "/g2015_2014_", admin_level, "/g2015_2014_", admin_level, "_modified.dbf"))
    admin_names <- as.data.table(admin_names)
    admin_names <- admin_names[, c(paste0('ADM', admin_level, '_CODE'), paste0('ADM', admin_level, '_NAME')), with=FALSE]
    all_admins <- merge(all_admins, admin_names, by=paste0('ADM', admin_level, '_CODE'), all.x=TRUE)
  }

  message("saving csvs")
  write.csv(all_admins,
            paste0(results_dir, '/summary_tables/', ind, '_', measure, '_', metric, '_summary_table.csv'))

  ## Save copy of all_admins (by level) for Lucas
  for(adm in c(0,1,2)) {
    if(adm == 0) {
      this_admin_data <- all_admins[is.na(ADM1_CODE) & is.na(ADM2_CODE), ]
      this_admin_data <- this_admin_data[year == year_to_map, ]
    }
    if(adm == 1) {
      this_admin_data <- all_admins[!is.na(ADM1_CODE) & is.na(ADM2_CODE), ]
      this_admin_data <- this_admin_data[year == year_to_map, ]
    }
    if(adm == 2) {
      this_admin_data <- all_admins[!is.na(ADM1_CODE) & !is.na(ADM2_CODE), ]
      this_admin_data <- this_admin_data[year == year_to_map, ]
    }
    write.csv(this_admin_data, paste0(out.dir, ind, '_probMeetSDG_', year_to_map, 'ad', adm, '_summary_table.csv'))
  }

  message("making africa rasters")
  make_africa_raster <- function(gaul, pixel_year) {

    # Convert raster to SpatialPointsDataFrame
    message(paste0('rasterizing ', gaul, '...'))
    pixel_probs <- fread(paste0(results_dir, '/pixels/', measure, '_', gaul, '_probs_', pixel_year, '.csv'))
    simple_raster <- raster(paste0(results_dir, '/simple/', gaul, '.tif'))
    probs_raster <- insertRaster(simple_raster, matrix(pixel_probs[, p_goal], ncol = 1))
    return(probs_raster)

  }

  africa_raster <- lapply(unique(all_admins[, ADM0_CODE]), make_africa_raster, year_to_map)
  final_africa_raster = do.call(raster::merge, africa_raster)
  writeRaster(final_africa_raster,
              file = paste0(out.dir, ind, '_probMeetSDG_', year_to_map),
              format='GTiff',
              overwrite = TRUE)

}


make_africa_raster <- function(gaul) {

  # Convert raster to SpatialPointsDataFrame
  message(paste0('rasterizing ', gaul, '...'))
  pixel_probs <- fread(paste0(results_dir, '/pixels/', measure, '_', gaul, '.csv'))
  pixel_probs <- pixel_probs[year==2015, ]
  simple_raster <- raster(paste0(results_dir, '/simple/', gaul, '.tif'))
  probs_raster <- insertRaster(simple_raster, matrix(pixel_probs[, mean], ncol = 1))
  return(probs_raster)

}

#################################################################################
### Generates aggregated estimates in case of logit raking
## Modified from calc_raking_factors()
#################################################################################
get_aggs <- function(agg_geo_est = cond_sim_adm0,
                     gaul_list   = gaul_list,
                     rake_to     = gbd) {

  message('WARNING: function will not work as expected if agg_geo_est and rake_to are not aggregated at the same admin level.')

  if(sum(colnames(rake_to) %in% c('name','year','mean')) != 3)
      stop('rake_to should be a data table with column names `name`, `year`, and `mean`')

  if(dim(t(as.matrix(agg_geo_est)))[1]!=1 | is.null(names(agg_geo_est)))
    stop('agg_geo_est must be a named vector returned from make_condSim() with summarize option equal to TRUE')

  # transpose and rename agg_geo_est
  agg_geo_est <- data.table(split_geo_names(as.matrix(agg_geo_est)),agg_geo_est)
  agg_geo_est$year <- as.numeric(agg_geo_est$year)
  agg_geo_est$name <- as.numeric(agg_geo_est$name)

  # merge
  merged <- merge(agg_geo_est,rake_to, by=c('name','year'),all.x=T)
  names(merged)[names(merged)=='agg_geo_est'] <- 'geo_mean'
  names(merged)[names(merged)=='mean'] <- 'rake_to_mean'

  return(merged)

}

#' extract_from_brick
#'
#' @description Extract values from a raster or rasterbrick to lat/lon points by relevant year of data.
#'
#' @author Rebecca Stubbs
#'
#' @param data A data.table or data.frame (that will be converted to a data.table) with
#'             columns for latitude, longitude, and year.
#' @param yearcol A string name of a column that defines the year of the data within the data object
#' @param varname A string name of what you want the column of extracted values to be
#' @param raster_years A vector that describes the years that serve as bands in the rasterbrick, in order.
#' @param rb A rasterbrick with the bands that correspond to the raster_years.
#' @return Returns the data object, sorted based on year, with a new column with the extracted values (name
#' as specified based on the varname parameter).

extract_from_brick<-function(data,
                             yearcol="original_year",
                             rb,
                             varname="extracted_values",
                             raster_years=2000:2015){

  # Making sure you've passed the right arguments to the function:
  if(sum(c("longitude","latitude",yearcol) %in% names(data))!=3){
    stop("You need to have the columns longitude, latitude, and a specified column that describes year (as the yearcol parameter) in the data object passed to this function.")
  }
  if(varname=="extracted_values"){warning("The variable added to the data.table will be named 'extracted_values' unless you supply a different name to the varname parameter.")}

  if(min(data[[yearcol]])<min(raster_years)){stop("You have years of raw input data that extend beyond the min or max of the raster years supplied to the function. Consider setting years outside those bounds to the min/max of the raster years.")}

  # Ordering the data object by year
  data<-data[order(data[[yearcol]])]

  # Renaming the rasterbrick layers to be descriptive of year
  names(rb)<-paste0("year_",as.character(raster_years))

  # Define a vector to store the results
  extracted_values<-c()

  for (yr in raster_years){
    message(paste0("Extracting values for ",yr))
    # Select only the data points from that year
    singleyr<-data[data[[yearcol]]==yr,]

    # Convert to spatial object
    coordinates(singleyr) <- ~longitude + latitude

    # setting the proj4string (the projection information) of this layer to be the
    # same as the results raster
    proj4string(singleyr)<-proj4string(rb)

    # Extract raster values to points
    values<-raster::extract(rb[[paste0("year_",as.character(yr))]], singleyr)
    extracted_values<-c(extracted_values,values)
  }

  data[[varname]]<-extracted_values
  return(data)
}

#' get_weighted_correlation_raw_estimated
#' @description Returns the weighted correlation between the raw data points input into a model run,
#' and the unraked estimates.
#' @author Rebecca Stubbs
#'
#' @param indicator_group The string of the indicator group.
#' @param indicator The string name of the indicator.
#' @param rundate The string name of the rundate.
#'
#' @return Returns the output from a call to the weights::wtd.cors() function
#' between the raw data points and the estimates.

get_weighted_correlation_raw_estimated<-function(indicator_group,
                                                 indicator,
                                                 rundate){

  source(paste0("<<<< FILEPATH REDACTED >>>>>",Sys.getenv("LOGNAME"),"/mbg/mbg_central/","prep_functions.R"))
  load_libs(c('data.table','raster','sp','weights'))

  input_folder<-paste0(root,"<<<< FILEPATH REDACTED >>>>>")

  message("Reading in the raw data")
  # Read in CSV, reduce table size by eliminating irrelevant columns
  raw_data<-fread(paste0(input_folder,indicator,".csv"))
  raw_data[,raw:=raw_data[[indicator]]/N]
  raw_data<-raw_data[,list(year=original_year,latitude,longitude,raw,N,weight)]

  message("Loading in the mean raster (Unraked)")
  # Load in Raster of Values
  results_tif<-raster::stack(paste0("<<<< FILEPATH REDACTED >>>>>",indicator_group,"/",indicator,"/output/",rundate,"/",indicator,"_mean_raster.tif"))
  names(results_tif)<-as.character(seq(2000,2015))

  raw_estimate_comparison<-list()

  message("Extracting values to points")
  for(yr in seq(2000,2015)){
    print(yr)
    raw_singleyr<-raw_data[year==yr,]

    # Convert to spatial object
    coordinates(raw_singleyr) <- ~longitude + latitude

    # setting the proj4string (the projection information) of this layer to be the
    # same as the results raster
    proj4string(raw_singleyr)<-proj4string(results_tif)

    # Extract raster values to points
    values<-raster::extract(results_tif[[paste0("X",as.character(yr))]], raw_singleyr)

    # Going back into table-space, and adding the data.table to a list to rbind later
    raw_singleyr<-as.data.table(raw_singleyr)
    raw_singleyr[,estimate:=values]
    raw_estimate_comparison[[as.character(yr)]]<-raw_singleyr
  }

  message("Calculating correlation coefficient using the weights::wtd.cors function")
  # Combining together the raw estimates into 1 data.table
  raw_estimate_comparison<-rbindlist(raw_estimate_comparison)
  cor_results<-weights::wtd.cors(x=raw_estimate_comparison$raw,
                                 y=raw_estimate_comparison$estimate,
                                 weight=raw_estimate_comparison$N * raw_estimate_comparison$weight)

  return(cor_results)
}


prep_postest <- function(indicator,
                         indicator_group,
                         run_date,
                         save_objs) {

  # Save a list of objects in a standard location for parallel scripts to pull from
  main_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/')
  temp_dir <- paste0(main_dir, "temp_post_est/")
  temp_file <- paste0(temp_dir, "post_est_temp_objs.RData")
  dir.create(temp_dir, showWarnings = F)
  save(list = save_objs, file = temp_file)
}

post_load_combine_save <- function(regions    = strata,
                                   summstats  = c('mean','cirange','upper','lower'),
                                   raked      = c('raked','unraked'),
                                   rf_table   = TRUE,
                                   run_summ   = TRUE){

  rake_addin <- character()
  if ("unraked" %in% raked) {
    lookup_dir <- paste0(sharedir, "/output/", run_date, "/")
    ur <- length(grep(paste0(indicator, ".*unraked.*raster.tif"), list.files(lookup_dir)))
    if (ur > 0) rake_addin <- c(rake_addin, unraked = "_unraked")
    if (ur == 0) rake_addin <- c(rake_addin, unraked = "")
  }

  if ("raked" %in% raked) {
    rake_addin <- c(rake_addin, raked = "_raked")
  }    

  # loop through and combine all rasters
  message("\nCombining rasters...")
  for(rake in rake_addin){
    message(names(rake_addin)[which(rake_addin == rake)])
    rr <- rake
    for(ss in summstats){
      message(paste0('  ',ss))
      rlist <- list()
      for(reg in regions){
        message(paste0('    ',reg))
        rlist[[reg]] <-
          brick(sprintf('%s/output/%s/%s_%s%s_%s_raster.tif',sharedir,run_date,indicator,reg,rake,ss))
      }
      if(ss=='cirange') ssname = 'range' else ssname = ss # naming convention
      save_post_est(do.call(raster::merge,unname(rlist)),'raster',
                            paste0(ssname,rr,'_raster'))
    }
  }

  # do rf also
  if(rf_table){
    message('RF table')
    rflist <- list()
    for(reg in regions){
      rflist[[reg]] <-
        read.csv(sprintf('%s/output/%s/%s_%s_rf.csv',sharedir,run_date,indicator,reg))
    }
    save_post_est(do.call(rbind.fill,rflist),'csv','rf')
  }

  # make a run summary graph
  if(run_summ) {
    graph_run_summary(run_date = run_date,
                      indicator_group = indicator_group, 
                      indicator = indicator)
  }

}

clean_after_postest <- function(indicator,
                                indicator_group,
                                run_date,
                                strata,
                                delete_region_rasters = F) {

  # Delete intermediate files that we no longer need
  main_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/')
  temp_dir <- paste0(main_dir,'temp_post_est/')
  
  # Deleting - be careful!
  unlink(temp_dir, recursive = T) 

  # If desired, insert code here to delete other temporary objects (eg. region-specific rasters)
  grep_string <- paste0(indicator, "_(", paste(strata, collapse = "|"), ").*_raster.tif")
  region_rasters <- grep(grep_string,list.files(main_dir), value=T) %>%
                       paste0(main_dir, .)
  if(delete_region_rasters == T) unlink(region_rasters)
}

submit_aggregation_script <- function(indicator, indicator_group, run_date, raked, 
                                      pop_measure, overwrite, ages, holdouts, regions, 
                                      repo, log_dir, geos_node = F, slots = 8) {

  # Takes vectors of regions, holdouts, and ages and submits qsubs for each of these

  dir.create(log_dir)
  dir.create(paste0(log_dir, '/errors'))
  dir.create(paste0(log_dir, '/output'))

  shell.file <- ifelse(geos_node,
                       "/mbg_central/r_shell_geos.sh",
                       "/mbg_central/r_shell.sh")

  proj.flag <- ifelse(geos_node,
                      "-P proj_geo_nodes -l gn=TRUE",
                      "-P proj_geospatial")

  qsubs_to_make <- expand.grid(regions, holdouts, ages, raked)

  for (i in 1:nrow(qsubs_to_make)) {
    region <- qsubs_to_make[i, 1]
    holdout <- qsubs_to_make[i, 2]
    age <- qsubs_to_make[i, 3]
    rake <- qsubs_to_make[i, 4]

    qsub <- paste0('qsub -e ', log_dir, '/errors -o ', log_dir, '/output',
                   ' -cwd -pe multi_slot ', slots, ' ', proj.flag, 
                   ' -N ',indicator, '_', region, '_aggregate ', 
                   repo, shell.file, ' ', repo, '<<<< FILEPATH REDACTED >>>>>/aggregate_results.R ', 
                   indicator, ' ',         # arg 4
                   indicator_group, ' ',   # arg 5
                   run_date, ' ',          # arg 6
                   rake, ' ',              # arg 7
                   pop_measure, ' ',       # arg 8
                   overwrite, ' ',         # arg 9
                   age, ' ',               # arg 10
                   holdout, ' ',           # arg 11
                   region)                 # arg 12

   system(qsub)
 }

 return(NULL)

}


