
##################################################################################################
##################################################################################################
## Overview:  Will search saved results and model objects to produce m draws of estimates
#             for each holdout area. Note this is written as a unit (ie for one reg, holdout),
#             an will likely be run in a loop or apply function
#
## Inputs:
#          indicator_group: indicator group
#          indicator: indicator
#          run_date: run date
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
  datdir <- '<<<< FILEPATH REDACTED >>>>'

  # cell draws
  cdfile <- '<<<< FILEPATH REDACTED >>>>'
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
  load(paste0('<<<< FILEPATH REDACTED >>>>/simple_raster',reg,'.RData'))

  # for each draw in cell_draws, pull the p
  r_list = list()
  for(i in 1:draws){
    r_list[[length(r_list)+1]] <- insertRaster(simple_raster, matrix(cell_pred[,i],ncol = length(years)))
  }

  # extract probabilities from draw rasters
  ycol = match(d_oos$year,years)
  for(i in 1:draws){
    # extract at locations
    t=raster::extract(r_list[[i]],cbind(d_oos$longitude,d_oos$latitude))
    # match year and put into df
    d_oos[,paste0(indicator,'hat_',i)]=      sapply(seq_along(ycol), function(x){t[x,ycol[x]]}) * d_oos$exposure
    d_oos[,paste0(indicator,'hat_full_',i)]= sapply(seq_along(ycol), function(x){t[x,ycol[x]]}) * d_oos$N

  }

  # get binomial draws for each point-draw, and estimate the coverage based on X%ci
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


## ################################
## is_oos_preds_testing
##
## this function takes the data that went into our MBG framework and
## returns an array containing: location, observed value, in sample
## predictions, and out of sample predictions (when applicable), country, year, N
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
  output.dir <- '<<<< FILEPATH REDACTED >>>>'

  ## Load data
  datdir <- '<<<< FILEPATH REDACTED >>>>'

  ## holdout data
  if(!is.null(holdoutlist)){
    d <- data.frame(holdoutlist[[sprintf('region__%s', reg)]])
  }else{
    d <- df
  }

  ## load the simple raster for this region
  load(paste0('<<<< FILEPATH REDACTED >>>>/simple_raster',reg,'.RData'))

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
      load(paste0('<<<< FILEPATH REDACTED >>>>', cell_draws_filename))

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
  load(sprintf(paste0('<<<< FILEPATH REDACTED >>>>', cell_draws_filename))

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
##   shapefile_version: String specifying shapefile version to pull
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
                             ind_fm = "binomial",
                             age = 0,
                             nperiod = 16, # FIXME: superfluous, can just use yrs
                             yrs = 2000:2015,
                             write.to.file = FALSE,
                             get.oos = FALSE,
                             year_col = "original_year",
                             shapefile_version = "current",
                             ...) {
  
  
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
  ##     if you use subregions) [DEPRECATED]
  ##   nperiod: number of periods/years in model
  ##   years: vector of years. should be of length nperiod
  ##   get.oos: should we also get OOS extracted values?
  ##   write.to.file: if true writes final df to standard output dir
  ##   shapefile_version: String specifying shapefile version to pull
  ##
  ## OUTPUT
  ##
  ##   a data frame with:
  ##     nrow = number data observations
  ##     columns for each draw and some identifying columns:
  ##       holdout_id, value, sample size, region
  ##
  ## #####################################################################
  
  ###############
  ## Load data ##
  ###############
  
  message("Load input data used in model")
  
  # load regions that were used in modeling
  mod.dir <- '<<<< FILEPATH REDACTED >>>>'
  all.regions <- get_output_regions(mod.dir)
  
  # load raw data
  if (!get.oos) {
    df <- fread(paste0(mod.dir, "input_data.csv")) ## raw input data
    df <- merge_with_ihme_loc(df, shapefile_version = shapefile_version)
  } else {
    df <- lapply(paste0(mod.dir, "stratum_", all.regions, ".rds"), readRDS)
    for (i in 1:length(df)) df[[i]] <- data.table(df[[i]][[1]])
    df <- rbindlist(df)
  }
  
  # rename year column for convenience
  setnames(df, year_col, "the_year_col")
  
  ###################
  ## Assign admins ##
  ###################
  
  message("Identify ad1 and ad2 membership")
  
  # load admin2 shapefile (also contains admin1)
  admin2_shapefile <- rgdal::readOGR(get_admin_shapefile(admin_level = 2, version = shapefile_version))
  for (v in grep("CODE", names(admin2_shapefile@data))) admin2_shapefile@data[[v]] <- as.numeric(as.character(admin2_shapefile@data[[v]]))
  
  # identify the admin2 (and by extension, the admin1) each point belongs to
  locs <- sp::SpatialPoints(cbind(df$longitude, df$latitude), proj4string = CRS(proj4string(admin2_shapefile)))
  adm.df <- sp::over(locs, admin2_shapefile)
  
  # for those that don't fall inside a polygon, assign them to the nearest polygon (this tends to happen on coastlines)
  # do this by country to speed things up and to make sure the point ends up at least in the correct admin0
  for (ctry in unique(df[is.na(adm.df[, 1]), GAUL_CODE])) {
    ii <- which(is.na(adm.df[, 1]) & df$GAUL_CODE == ctry)
    temp_shape <- admin2_shapefile[admin2_shapefile@data$ADM0_CODE == ctry, ]
    distmat <- gDistance(locs[ii], temp_shape, byid = T)
    jj <- apply(distmat, 2, which.min)
    adm.df[ii, ] <- temp_shape@data[jj, ]
    rm(ii, jj, temp_shape, distmat)
  }
  
  # copy admin information into df
  df$ad1 <- adm.df$ADM1_CODE
  df$ad2 <- adm.df$ADM2_CODE
  df$ad0 <- df$GAUL_CODE # for consistency...
  
  ###############
  ## Get draws ##
  ###############
  
  message("Get draws")
  
  # loop over regions
  df_all <- rbindlist(lapply(all.regions, function(rr) {
    message(paste("...Region:", rr))
    
    # load the simple raster template
    message("......load simple raster template")
    gaul_list <- get_adm0_codes(rr, shapefile_version = shapefile_version)
    simple_polygon_list <- load_simple_polygon(
      gaul_list = gaul_list, buffer = 0.4, subset_only = TRUE,
      shapefile_version = shapefile_version
    )
    subset_shape <- simple_polygon_list[[1]]
    raster_list <- build_simple_raster_pop(subset_shape)
    template <- raster_list[["simple_raster"]]
    
    # subset data to this region
    df.r <- df[region == rr, ]
    loc.r <- as.matrix(df.r[, list(longitude, latitude)])
    
    # 'nudge' coordinates to make them fall inside the nearest cell non-NA cell in template if they
    # are in an NA cell (this can happen near coastlines)
    check <- raster::extract(template, loc.r)
    miss <- which(is.na(check))
    for (ii in miss) {
      dist <- replace(raster::distanceFromPoints(template, loc.r[ii, ]), is.na(template), NA)
      loc.r[ii, ] <- raster::xyFromCell(template, raster::which.min(dist)[1])
    }
    
    # make a raster of cell_pred row ids
    id_raster <- insertRaster(template, matrix(1:(length(cellIdx(template)) * nperiod), ncol = nperiod))
    
    # get cell_pred row ids for each data point
    for (yy in 1:nperiod) {
      this_year <- which(df.r$the_year_col == yrs[yy])
      if (length(this_year) == 0) next
      df.r[this_year, cell_pred_id := raster::extract(id_raster[[yy]], loc.r[this_year, ])]
    }
    
    # if out of sample metrics are requested, duplicate the data to create separate rows for in and out of sample
    if (get.oos) {
      df.r <- rbind(df.r, cbind(df.r[, -"fold", with = F], fold = 0), use.names = T)
    } else {
      df.r[, fold := 0]
    }
    
    # loop over holdouts
    for (this_fold in sort(unique(df.r$fold))) {
      
      # load cell pred objects
      message(paste("......load cell preds for holdout", this_fold))
      load('<<<< FILEPATH REDACTED >>>>')
      
      # extract draws
      df.r[fold == this_fold, paste0("draw", 1:ncol(cell_pred)) := as.data.table(cell_pred[cell_pred_id, ])]
      
    }
    
    # return combined draws and draws
    return(df.r)
  }))
  
  # rename year column back to original
  setnames(df_all, "the_year_col", year_col)
  
  # save combined data and draws to file and return
  if (write.to.file) write.csv(df_all, paste0(mod.dir, "output_draws_data.csv"), row.names = F)
  return(df_all)
}

## ####################################
## a function to plot raking factors ##
## ####################################
## this function plots a scatterplot of GBD estimates vs MBG estimates
plot.rfs <- function(ind.gp = indicator_group,
                     ind = indicator,
                     rd = run_date,
                     output.dir = si.fig.dir,
                     shapefile_version = 'current', 
                     title = "Comparison to GBD 2016 in\n" ## region gets pasted on to this
){

  require(ggrepel) ## for labels

  regions <- get_output_regions(in_dir = '<<<< FILEPATH REDACTED >>>>')

  for(rr in regions){

    ## convert rrs to full names
    if(rr == 'essa') rr_name = "Eastern Sub-Saharan Africa"
    if(rr == 'wssa') rr_name = "Western Sub-Saharan Africa"
    if(rr == 'name') rr_name = "North Africa"
    if(rr == 'sssa') rr_name = "Southern Sub-Saharan Africa"
    if(rr == 'cssa') rr_name = 'Central Sub-Saharan Africa'

    in_dir  <- '<<<< FILEPATH REDACTED >>>>'
    default_rf_path <- paste0(in_dir, '/', ind, '_rf.csv')
    all_rfs <- fread(default_rf_path)
    gaul_list = get_adm0_codes(rr, shapefile_version = shapefile_version)
    rfs <- all_rfs[name %in% gaul_list, ]
    loc_names <- setDT(get_location_code_mapping(shapefile_version = shapefile_version))
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
  ggsave(filename = '<<<< FILEPATH REDACTED >>>>',
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
                                 sharedir = '<<<< FILEPATH REDACTED >>>>'){
  ## load models
  require(INLA)
  message(sprintf('Pulling together results for %s models',rd))

  tlist=list()

  for(rr in regions){
    message(sprintf('::on region %s',rr))
    reg  <-  rr

    message("::::loading in pre-INLA objects to get spde")
    pathaddin  <-  paste0('_bin',age,'_',rr,'_',holdout)
    load(paste0('<<<< FILEPATH REDACTED >>>>', rd, pathaddin,
                '.RData'))

    modnames = c('gam','gbm','ridge','enet','lasso')

    full_raster_list <- cov_list

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
    f <-  '<<<< FILEPATH REDACTED >>>>'
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
    if(file.exists('<<<< FILEPATH REDACTED >>>>')){
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
      for(f in list.files(path   ='<<<< FILEPATH REDACTED >>>>',pattern=patt)){
        message(f)
        tmp       <- fread('<<<< FILEPATH REDACTED >>>>')
        tmp$model <- gsub('.csv','',gsub(paste0(patt,'_'),'',f))
        d         <- rbind(d,tmp)
      }
    }

    # simplify d for now
    d <- d[,c('p','mean','error','model','ho_id','year','fold','age','exposure','clusters_covered_95', 'total_clusters'),with=F]

    # coverage
    #d[,cov:=clusters_covered_95/total_clusters]

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
      pdf(NULL) # Workaround for bug in ggplot_gtable causing empty Rplots.pdf to be created
      tmp <- ggplot_gtable(ggplot_build(myggplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      graphics.off()
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

    message(paste0('saving ', '<<<< FILEPATH REDACTED >>>>'))
    pdf(s'<<<< FILEPATH REDACTED >>>>')
    grid.arrange(p1, p2, p3, legend, ncol = 2, top = sprintf('OOS Predictive Validity: Age Bin %i, %s. %i/%i folds analyzed',age,reg,length(folds),nfolds))
    dev.off()
  } else {
    message('NO DATA WRITTEN YET, SO NO PLOTS')
  }
  #return(plot)
}


####################################################################################################################
## INPUT: # data file matching <<<< FILEPATH REDACTED >>>>/output_draws_data.csv
## OUTPUT SUMMARIZED PREDICTIVE VALIDITY METRICS
##
####################################################################################################################

## get_pv_table ################################################

#' @title Get predictive validity table
#' @description Generate plots and tables of metrics of predictive validity
#'
#' @param d data.table containing data and IS/OOS draws from `run_is_oos()`
#' @param indicator indicator
#' @param indicator_group indicator_group
#' @param rd run date
#' @param aggregate_on column of spatial aggregation ("country", "ad1", or "ad2" usually) - can pass vector
#' @param draws number of draws (numeric)
#' @param coverage_probs probability at which you want to assess coverage of predictive intervals
#'                       can be a single number or a vector; if a vector can produce calibration plots
#'                       assessing x\% coverage at various cutoffs
#' @param result_agg_over what should the final table aggregate over (vector)? For instance
#'                        c("year", "oos") are the defaults.  If you include "region" here
#'                        the function will produce a separate set of validation plots,
#'                        one for each region
#' @param weighted do you want to weight PV metrics on sample size? (Boolean)
#' @param family distribution to use ("binomial" or "gaussian")
#' @param plot produce plots? (Boolean)
#' @param plot_by produce separate plots for an element of `result_agg_over`?
#'                for instance, `plot_by = region` with `result_agg_over = c("year", "oos", "region"`
#'                will produce a separate plot for each region/oos combo
#' @param plot_by_title title for captions for "plot_by" object, i.e. "Region" (defaults to plot_by)
#' @param plot_ci plot CIs (Boolean)
#' @param plot_ci_level what ci level to plot (numeric, i.e. "95" is the 95\% UI)
#' @param ci_color what color do you want the CI lines to be?
#' @param point_alpha how transparent do you want the points themselves to be (numeric, [0-1])
#' @param point_color what color do you want the points to be in your validation plots?
#' @param plot_title main title for your plots (defaults to "indicator")
#' @param plot_ncol how many columns do you want your plots to have?
#' @param save_csv save a csv table of your predictive validity metrics?
#' @param out.dir where do you want these results to be saved?
#'
#' @return data.table of PV metrics; plots and tables of PV metrics depending on the combinations of options above
#'
#' @examples
#' 
#' \dontrun{
#' #' # A PV table by modeling regions with plots (one for each modeling region) and semi-transparent points
#' # Also includes calibration plots since there are multiple coverage levels
#' 
#' # Get in and out of sample draws
#' run_in_oos <- get_is_oos_draws(
#'   ind_gp = indicator_group,
#'   ind = indicator,
#'   rd = run_date,
#'   ind_fm = "binomial",
#'   model_domain = Regions,
#'   age = 0,
#'   nperiod = length(year_list),
#'   yrs = year_list,
#'   get.oos = as.logical(makeholdouts),
#'   write.to.file = TRUE,
#'   year_col = "year"
#' )
#' 
#' # Load the IS/OOS draws created above
#' draws.df <- fread('<<<< FILEPATH REDACTED >>>>')
#' 
#' # run the PV table function
#' pvtable.reg <- get_pv_table(
#'   d = draws.df,
#'   indicator_group = indicator_group,
#'   rd = run_date,
#'   indicator = indicator,
#'   result_agg_over = c("year", "oos", "region"),
#'   coverage_probs = seq(from = 5, to = 95, by = 5),
#'   aggregate_on = "ad2",
#'   draws = as.numeric(samples),
#'   out.dir = out_dir,
#'   plot = TRUE,
#'   plot_by = "region",
#'   plot_by_title = "Region",
#'   plot_ci = TRUE,
#'   point_alpha = 0.1,
#'   point_color = "black",
#'   ci_color = "gray",
#'   plot_title = plot_title
#' )
#' }
#' 
#' @export
get_pv_table <- function(d,
                         indicator,
                         indicator_group,
                         rd,
                         aggregate_on,
                         draws = 1000,
                         coverage_probs = c(95),
                         result_agg_over = c("year", "oos"),
                         weighted = TRUE,
                         family = "binomial",
                         plot = TRUE,
                         plot_by = NULL,
                         plot_by_title = NULL,
                         plot_ci = FALSE,
                         plot_ci_level = 95,
                         ci_color = "grey",
                         point_alpha = 1,
                         point_color = "black",
                         plot_title = indicator,
                         plot_ncol = 4,
                         save_csv = T,
                         out.dir) {
  
  str_match <- stringr::str_match
  
  d <- data.table(d)
  
  if (!is.null(plot_by)) {
    if (!(plot_by %in% result_agg_over)) {
      stop("If you specify `plot_by`, you must also include that item in `result_agg_over`")
    }
  }
  
  ## Get binomial predictions
  message("Simulating predictive draws")
  if (family == "binomial") {
    x <- sapply(1:draws, function(i) rbinom(nrow(d), round(d$N), d[[paste0("draw", i)]]))
    d[, Y := get(indicator) * round(N) / N] # adjust observed number of cases for effect of rounding N
  }
  if (family == "gaussian") {
    x <- sapply(1:draws, function(i) rnorm(nrow(d), sqrt(d[[paste0("tau_", i)]]), d[[paste0("draw", i)]]))
    d[, Y := get(indicator)]
  }
  
  message(paste0("...NOTE: missing predictions for ", sum(is.na(x[, 1])), " of ", nrow(x), " (", round(100 * mean(is.na(x[, 1]))), "%) data points."))
  
  ## Get coverage for each data point
  message("Calculate coverage")
  one_side <- (1 - coverage_probs / 100) / 2
  ui <- t(apply(x, 1, quantile, c(one_side, 1 - one_side), na.rm = T))
  for (c in coverage_probs) {
    c_id <- which(c == coverage_probs)
    d[, paste0("clusters_covered_", c) := Y %between% list(ui[, c_id], ui[, c_id + length(coverage_probs)])]
  }
  
  ## Collapse data and draws spatially
  message("Collapse data and draws spatially")
  d[, oos := (fold != 0)]
  d[, p := get(indicator) / N]
  d[, exposure := N * weight]
  
  # Create list of pv metrics for all levels of aggregation in `aggregate_on`
  message("Outputting predictive validity metrics for your models at each level of aggregation")
  pv_table_list <- lapply(aggregate_on, function(agg_on) {
    message(paste0("  ", agg_on))
    by_vars <- unique(c("oos", "year", result_agg_over, agg_on))
    collapse_vars <- c("p", paste0("clusters_covered_", coverage_probs), paste0("draw", 1:draws))
    
    res <- d[!is.na(draw1),
             c(
               list(total_clusters = .N, exposure = sum(exposure)),
               lapply(.SD, function(x) weighted.mean(x, exposure, na.rm = T))
             ),
             keyby = by_vars, .SDcols = collapse_vars
             ]
    res[, mean_draw := rowMeans(.SD), .SDcols = paste0("draw", 1:draws)]
    res[, error := p - mean_draw]
    res[, abs_error := abs(error)]
    
    ## Collapse to calculate predictive validity metrics
    if (weighted) res$weight <- res$exposure else res$weight <- 1
    res2 <- res[, c(lapply(.SD, function(x) weighted.mean(x, weight)),
                    rmse = rmse(error, weight),
                    median_SS = median(exposure),
                    cor = boot::corr(cbind(p, mean_draw), weight)
    ),
    by = result_agg_over,
    .SDcols = c("error", "abs_error", "mean_draw", "p", paste0("clusters_covered_", coverage_probs))
    ]
    setnames(
      res2, c("error", "abs_error", "p", paste0("clusters_covered_", coverage_probs)),
      c("me", "mae", "mean_p", paste0("coverage_", coverage_probs))
    )
    
    return(list(res = res, res2 = res2))
  })
  
  names(pv_table_list) <- aggregate_on
  
  # Make plots
  if (plot == TRUE) {
    message("Making plots of aggregated data and estimates")
    
    # Get unique levels of `plot_by` and set up a plot_by title if needed
    if (!is.null(plot_by)) {
      plot_by_levels <- unique(pv_table_list[[1]][["res"]][, get(plot_by)])
    } else {
      plot_by_levels <- NULL
    }
    
    if (is.null(plot_by_title)) plot_by_title <- plot_by
    
    # Create a table of things to plot
    plot_table <- CJ(
      aggregate_on = aggregate_on,
      oos = unique(pv_table_list[[1]][["res"]]$oos),
      plot_by_value = if (is.null(plot_by_levels)) NA else plot_by_levels
    )
    
    message("...saving plots here: ", out.dir)
    # Loop over plots
    for (i in 1:nrow(plot_table)) {
      
      # Grab items from plot table
      agg_on <- plot_table[i, aggregate_on]
      oosindic <- plot_table[i, oos]
      pb_val <- plot_table[i, plot_by_value]
      
      # Set up titles
      if (agg_on == "country") agg_title <- "Country"
      if (agg_on == "ad1") agg_title <- "Admin 1"
      if (agg_on == "ad2") agg_title <- "Admin 2"
      
      res <- pv_table_list[[agg_on]][["res"]]
      res2 <- pv_table_list[[agg_on]][["res2"]]
      
      # Make a validation plot -----------------------------------------------------
      
      # Set up filename and file
      plot_filename <- paste0(
        indicator, "_validation_plot_",
        paste(c(
          as.character(agg_on),
          setdiff(result_agg_over, "oos")
        ),
        collapse = "_"
        ), "_",
        ifelse(oosindic, "OOS", "IS"),
        ifelse(is.na(pb_val), "", paste0("_", pb_val)),
        ".png"
      )
      message(paste("    ", plot_filename))
      png(paste0(out.dir, plot_filename), width = 12, height = 12, units = "in", res = 350)
      
      # Subset data
      fdata <- res[oos == oosindic, ]
      if (!is.na(pb_val)) {
        setnames(fdata, plot_by, "plot_by_column") # convenience
        fdata <- fdata[plot_by_column == pb_val, ]
      }
      
      # Set up CI bar limits; range as defaults
      if (plot_ci) {
        fdata[, upper := apply(.SD, 1, quantile, p = 0.01 * (plot_ci_level + (100 - plot_ci_level) / 2), rm.na = TRUE), .SDcols = paste0("draw", 1:draws)]
        fdata[, lower := apply(.SD, 1, quantile, p = 0.01 * ((100 - plot_ci_level) / 2), rm.na = TRUE), .SDcols = paste0("draw", 1:draws)]
        limits <- fdata[, range(c(p, mean_draw, lower, upper))]
      } else {
        limits <- fdata[, range(c(p, mean_draw))]
      }
      
      # The plot code itself
      gg <- ggplot(fdata, aes(x = p, y = mean_draw, size = weight)) +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        geom_point(colour = point_color, alpha = point_alpha) +
        scale_size_area() +
        ggplot2::xlim(limits) +
        ggplot2::ylim(limits) +
        coord_equal() +
        theme_bw() +
        theme(strip.background = element_rect(fill = "white")) +
        labs(
          x = "Data Estimate",
          y = "Mean Prediction",
          size = "Weight",
          title = paste0("Validation Plot for ", plot_title, " by ", agg_title),
          subtitle = paste0("OOS: ", oosindic, ifelse(is.na(pb_val), "", paste0(" | ", plot_by_title, ": ", pb_val)))
        )
      
      if (plot_ci) {
        gg <- gg + geom_errorbar(aes(ymin = lower, ymax = upper), colour = point_color, width = 0, size = .3, alpha = min(point_alpha, 0.2))
      }
      if (length(setdiff(result_agg_over, "oos")) > 0) {
        gg <- gg + facet_wrap(as.formula(paste("~", paste(setdiff(result_agg_over, c("oos", plot_by)), collapse = "+"))),
                              ncol = plot_ncol
        )
      }
      
      plot(gg)
      dev.off()
      
      # Make a calibration plot -----------------------
      if (plot == T & length(coverage_probs) > 1) {
        
        # Set up a subset of the data for plotting
        fdata <- res2[, unique(c("oos", result_agg_over, paste0("coverage_", coverage_probs))), with = F]
        fdata <- fdata[oos == oosindic]
        if (!is.na(pb_val)) {
          setnames(fdata, plot_by, "plot_by_column") # convenience
          fdata <- fdata[plot_by_column == pb_val, ]
        }
        fdata <- melt(fdata,
                      id.vars = names(fdata)[!grepl("coverage", names(fdata))],
                      value.name = "observed_coverage",
                      variable.name = "coverage"
        )
        fdata[, expected_coverage := as.numeric(gsub("coverage_", "", coverage))]
        fdata[, observed_coverage := observed_coverage * 100]
        fdata$group <- apply(fdata[, result_agg_over[!(result_agg_over %in% c("oos", "plot_by_column", plot_by))], with = F], 1, paste, collapse = " ")
        if (sum(!is.na(fdata$group)) == 0) fdata[, group := "All"]
        
        # Create filename
        
        cplot_filename <- paste0(
          indicator, "_calibration_plot_",
          paste(c(
            as.character(agg_on),
            setdiff(result_agg_over, "oos")
          ),
          collapse = "_"
          ), "_",
          ifelse(oosindic, "OOS", "IS"),
          ifelse(is.na(pb_val), "", paste0("_", pb_val)),
          ".png"
        )
        message(paste("    ", cplot_filename))
        png(paste0(out.dir, cplot_filename), width = 12, height = 12, units = "in", res = 350)
        
        # Set limits and plot
        limits <- fdata[, range(c(observed_coverage, expected_coverage))]
        gg <- ggplot(fdata, aes(x = expected_coverage, y = observed_coverage, group = group, color = group)) +
          geom_abline(intercept = 0, slope = 1, color = "red") +
          geom_point() +
          geom_line(alpha = 0.2) +
          scale_color_discrete(name = "") +
          coord_equal() +
          xlim(limits) +
          ylim(limits) +
          theme_bw() +
          labs(x = "Expected coverage", y = "Observed coverage")
        
        plot(gg)
        dev.off()
      } # End calibration plot if statement
    } # End plot table loop
  } # End if (plot==T) loop
  
  # Format, save (if desired), and return `res2` objects
  output_list <- lapply(aggregate_on, function(agg_on) {
    res2 <- copy(pv_table_list[[agg_on]][["res2"]])
    setorderv(res2, result_agg_over)
    setnames(res2, result_agg_over, ifelse(result_agg_over == "oos", "OOS", gsub("(^.)", "\\U\\1", result_agg_over, perl = T)))
    setnames(
      res2, c("me", "mae", "mean_draw", "mean_p", "rmse", "median_SS", "cor", paste0("coverage_", coverage_probs)),
      c("Mean Err.", "Mean Abs. Err.", "Mean Pred.", "Mean Obs.", "RMSE", "Median SS", "Corr.", paste0(coverage_probs, "% Cov."))
    )
    
    # Save final tables if desired
    if (save_csv) {
      a_pathaddin <- ifelse(length(setdiff(result_agg_over, "oos")) > 0,
                            paste0("_by_", paste0(setdiff(result_agg_over, "oos"), collapse = "_")),
                            ""
      )
      filename <- paste0(out.dir, agg_on, "_metrics", a_pathaddin, ".csv")
      
      message(paste0("Saving csv to ", filename, "..."))
      write.csv(res2, file = filename)
    }
    return(res2)
  })
  
  names(output_list) <- aggregate_on
  
  return(output_list)
}


## get_admin_pv function ###################################################

#' @title get_admin_pv
#' @description This function collapses input data to admin 0/admin 1/admin 2 using the `input_aggregate_admin` function
#' See `input_aggregate_admin` function documentation for input data requirements 
#' 
#' To use this function, your input dataset must have a column that is the sum of the sample weights
#'
#' This function pulls the input data from the model directory and aggregated admin summaries created in the `aggregate_results.R` script
#' 
#' This function returns average RMSE, Bias, Mean Absolute Error, and Standard Error for admin 0 /admin 1 / admin 2
#'
#' @param indicator indicator name used in file structure for mbg
#' @param indicator_group indicator group
#' @param strata Regions specified from your model
#' @param run_date  model run date
#' @param input_data If specified, provides the preloaded input data so the function does not look in your model directory
#' @param indicator_family If specified as Gaussian, this makes sure to not divide the prevalence by N which is required for binomial indicators
#' @param nic_col This is the name of the unique survey variable, usually labeled "nid" but other teams might use different terminology
#' @param samp_col This is the name of the column that contains the sum of the sample weights for the collapsed points.
#' @param shapefile_version String indicating version of shapefile to pull
#' @param input_file If your team does not use the standard `'<<<< FILEPATH REDACTED >>>>'` input data filepath, supply the filepath to your data here
#' @param admin0_file If your team does not use the standard `'<<<< FILEPATH REDACTED >>>>'` directory to save unraked aggregated admin summaries, supply the correct filepath here
#' @param admin1_file If your team does not use the standard `'<<<< FILEPATH REDACTED >>>>'` directory to save unraked aggregated admin summaries, supply the correct filepath here
#' @param admin2_file If your team does not use the standard `'<<<< FILEPATH REDACTED >>>>'` directory to save unraked aggregated admin summaries, supply the correct filepath here
#' 
#' @return a table with bias, rmse, mae, and se for each admin level
#' @export


get_admin_pv <- function(indicator, 
                         indicator_group, 
                         run_date, 
                         input_file = NULL, 
                         shapefile_version = 'current',
                         samp_col = 'sum_of_sample_weights',
                         strata = strata,
                         indicator_family = 'binomial',
                         nid_col = 'nid',
                         admin0_file = NULL,
                         admin1_file = NULL,
                         admin2_file = NULL
){
  message('Pulling in input data...')
  if(is.null(input_file)){
    if(!file.exists(paste0('<<<< FILEPATH REDACTED >>>>', indicator, '.csv'))) stop('Indicator input data does not exist in /share - please specify a file path with the input_file arg')
    input.dt <- fread(paste0('<<<< FILEPATH REDACTED >>>>', indicator, '.csv'))
  } else{
    input.dt <- fread(input_file)
  }
  
  if(!('point' %in% colnames(input.dt))) stop('Your input dataset needs a binary var called point')
  if(!(samp_col %in% colnames(input.dt))) stop('Your input dataset needs a column that is the sum of sample weights')
  
  message("Collapsing data to admin levels...")
  admin_data <- input_aggregate_admin(indicator = indicator, indicator_group, run_date = run_date, sample_column = samp_col, 
                                      input_data = input.dt, regions = strata, indicator_family = ind_fam, svy_id = nid_col, shapefile_version = shapefile_version)
  
  ad0_data <- admin_data$ad0
  ad1_data <- admin_data$ad1 
  ad2_data <- admin_data$ad2 
  
  #Read in admin level results
  message("Pulling in admin results")
  if(is.null(admin0_file)){
    ad0_mbg <- fread(paste0('<<<< FILEPATH REDACTED >>>>',
                            indicator, '_admin_0_unraked_summary.csv'))
  } else{
    ad0_mbg <- fread(admin0_file)
  }
  if(is.null(admin1_file)){
    ad1_mbg <- fread(paste0('<<<< FILEPATH REDACTED >>>>',
                            indicator, '_admin_1_unraked_summary.csv'))
  } else{
    ad1_mbg <- fread(admin1_file)
  }
  if(is.null(admin2_file)){
    ad2_mbg <- fread(paste0('<<<< FILEPATH REDACTED >>>>',
                            indicator, '_admin_2_unraked_summary.csv'))
  } else{
    ad2_mbg <- fread(admin2_file)
  }
  
  #Merge nid and mbg results
  
  ad0 <- merge(ad0_data, ad0_mbg, by = c('ADM0_NAME', 'ADM0_CODE', 'year'))
  ad1 <- merge(ad1_data, ad1_mbg, by = c('ADM0_NAME', 'ADM0_CODE', 'ADM1_NAME', 'ADM1_CODE', 'year'))
  ad2 <- merge(ad2_data, ad2_mbg, by = c('ADM0_NAME', 'ADM0_CODE', 'ADM1_NAME', 'ADM1_CODE','ADM2_NAME', 'ADM2_CODE', 'year'))
  
  
  # Biasa: paNID-mean(pa,yNID,1mbg, ... ,pa,yNID,ndrawsmbg)   
  # MAEa: |paNID-mean(pa,yNID,1mbg, ... ,pa,yNID,ndrawsmbg)| 
  # SEa: (paNID-mean(pa,yNID,1mbg, ... ,pa,yNID,ndrawsmbg))2
  #and then we average over admin(-years) and take the square root to get RMSE across admin(-years)
  message("Calculate predictive validity...")
  
  ad0[, bias := outcome - mean]
  ad0[, mae := abs(bias)]
  ad0[, se := bias^2]
  
  ad1[, bias := outcome - mean]
  ad1[, mae := abs(bias)]
  ad1[, se := bias^2]
  
  ad2[, bias := outcome - mean]
  ad2[, mae := abs(bias)]
  ad2[, se := bias^2]
  
  pv_table <- data.table(ad_level = c(0,1,2),
                         rmse = c(sqrt(mean(ad0$se, na.rm = T)), sqrt(mean(ad1$se, na.rm = T)), sqrt(mean(ad2$se, na.rm = T))),
                         bias = c(weighted.mean(ad0$bias, ad0$N, na.rm = T), weighted.mean(ad1$bias, ad1$N, na.rm = T), weighted.mean(ad2$bias, ad2$N, na.rm = T)), 
                         mae = c(weighted.mean(ad0$mae, ad0$N, na.rm = T), weighted.mean(ad1$mae, ad1$N, na.rm = T), weighted.mean(ad2$mae, ad2$N, na.rm = T)), 
                         se = c(weighted.mean(ad0$se, ad0$N, na.rm = T), weighted.mean(ad1$se, ad1$N, na.rm = T), weighted.mean(ad2$se, ad2$N, na.rm = T)))
  
  return(pv_table)
} 