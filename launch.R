## ##########################################################################
## ##########################################################################
## This script is a generic launch script used to launch any/all CGF
## models it should be called from within a launch<indicator>.R script
## which is used to set the indicator and specify the config.
## ##########################################################################
## ##########################################################################

## Load libraries and miscellaneous MBG project functions.
repo <- '<<<< FILEPATH REDACTED >>>>>/lbd_core/'
setwd(repo)
user_repo <- '<<<< FILEPATH REDACTED >>>>>'
root <- ifelse(Sys.info()[1]=="Windows", "<<<< FILEPATH REDACTED >>>>>", "<<<< FILEPATH REDACTED >>>>>")
package_lib <- ifelse(grepl("geos", Sys.info()[4]),
                      paste0(root,'temp/geospatial/geos_packages'),
                      paste0(root,'temp/geospatial/packages'))

package_list <- c('foreign', 'data.table','raster','rgdal', 'INLA',
                  'seegSDM','seegMBG','plyr','dplyr', 'boot',
                  'ggplot2', 'rgeos', 'tictoc')

if(Sys.info()[1] == 'Windows'){
  stop('STOP! you will overwrite these packages if you run from windows\n
        STOP! also, lots of this functions wont work so get on the cluster!')
} else {
  for(package in package_list)
    library(package, lib.loc = package_lib, character.only=TRUE)
  for(funk in list.files(path=user_repo, recursive=TRUE,pattern='functions')){
    if(length(grep('central',funk))!=0 | length(grep(indicator_group,funk))!=0){
        message(paste0('loading ',funk))
        source(funk)
    }
  for(funk in list.files(recursive=TRUE,pattern='functions')){
    if(length(grep('central',funk))!=0 | length(grep(indicator_group,funk))!=0){
        message(paste0('loading ',funk))
        source(funk)
    }
  }
}


## Read config file and save all parameters in memory
config <- load_config(repo = user_repo,
                      indicator_group = indicator_group,
                      indicator = indicator)

## Create run date in correct format unless it's in the config as 'NULL'
run_date <- make_time_stamp(time_stamp)

## Create directory structure for this model run
create_dirs(indicator_group = indicator_group,
            indicator = indicator)

## Load simple polygon template to model over
gaul_list <- get_gaul_codes('africa')
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                           buffer = 0.4)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]

## Load input data and make sure it is cropped to modeling area
df <- load_input_data(indicator = indicator,
                      simple = simple_polygon,
                      removeyemen = TRUE,
                      use_share = TRUE)

## Add GAUL_CODE and region to df given the lat/longs. We need this to
## stratify in holdout function.
df <- add_gauls_regions(df = df,
                        simple_raster = simple_raster)

## Run function to create holdouts (returns list of data.frames with
## an additional "folds" column)
table(df$region, df$year)
table(df$region, df$original_year)

## should we run holdouts??
if(n_ho_folds > 1){
  
  stratum_qt <- stratum_ho <- readRDS(paste0('<<<<< FILEPATH REDACTED >>>>>',
                                             indicator_group, '/',
                                             indicator, '/output/',
                                             run_date, "/",
                                             "stratum_qt.RDs"))

  ## resave them for later processing
  saveRDS(file = paste0('<<<<< FILEPATH REDACTED >>>>>', indicator_group,
                        '/', indicator, '/output/', run_date, "/",
                        "stratum_qt.RDs"),
          object = stratum_qt)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~  Parallel MBG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~ Submit job by strata/holdout  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## This is where we begin to do things by region

parallel_script <- 'parallel_model_full'

## last few things to prep
slots <- as.numeric(slots)
sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>>/mbg/%s/%s',indicator_group,indicator)

## loop through and launch models by strata and or region
if(n_ho_folds > 1){
  all_holdouts <- 0:n_ho_folds
}else{
  all_holdouts <- 0
}
loopvars <- NULL
strata <- Regions
for(r in strata){
  for(holdout in all_holdouts) {

    ## make qsub and save object for parallel script
    qsub <- make_qsub_share(reg = r,
                            saveimage = TRUE,
                            test = as.numeric(test),
                            holdout = holdout,
                            log_location = 'sharedir',
                            geo_nodes = as.logical(use_geos_nodes),
                            memory = as.numeric(mem),
                            additional_job_name = eval(parse(text = jn)), ## from config
                            cores = as.numeric(slots))
    
    ## subimt qsub
    system(qsub)
    
    ## keep track of loops iters submitted
    loopvars <- rbind(loopvars, c(r,holdout))
    
  } ## holdouts
}  ## strata

## ~~~~~~~~~~~~~~~~~~~~~~~~ Post-Estimation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## check to make sure models are done before continuing
waitformodelstofinish()

## Save strata for Shiny to use in producing aggregated fit statistics
dir.create(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/',
                  indicator, '/output/', run_date, '/fit_stats'))
save(strata, file = paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group,
                           '/', indicator, '/output/', run_date,
                           '/fit_stats/strata.RData'))

## Post estimation to be done by strata (in my case, this is just region)

## Set strata as character vector of each strata (in my case, just
## stratifying by region whereas U5M stratifies by region/age)
strata <- Regions

gbd <- as.data.table(read.csv(paste0(repo,"child_growth_failure/data_prep/", indicator, "_rake_2016.csv"),
                              stringsAsFactors = FALSE))
new_gbd <- interpolate_gbd(gbd)

for(reg in strata){
  message(reg)

  load(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/',
              indicator, '/output/', run_date, '/', indicator,
              '_cell_draws_eb_bin0_', reg, '_0.RData'))

  ## get aggregated estimates for all admin0. Aggregate to level you rake to
  simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg),
                                             buffer = 0.4,
                                             subset_only = FALSE)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

  ## Pull 2000-2015 annual population brick using new covariates function
  pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',
                                                       measures = 'total',
                                                       simple_polygon = simple_polygon,
                                                       start_year  = min(eval(parse(text = year_list))),
                                                       end_year    = max(eval(parse(text = year_list))),
                                                       interval_mo = 12,
                                                       agebin=1)
  pop_raster_annual  <- pop_raster_annual[[1]]
  pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
  pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
  pop_raster_annual  <- mask(pop_raster_annual, simple_raster)

  ## Create population weights using the annual brick and feed custom year argument to aggregation function
  pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                          simple_raster = simple_raster,
                                          pop_raster    = pop_raster_annual,
                                          gaul_list     = get_gaul_codes(reg))

  ## find national level estimates
  cond_sim_raw_adm0 <- make_condSim(pop_wts_object = pop_wts_adm0,
                                    gaul_list      = get_gaul_codes(reg),
                                    admin_level    = 0,
                                    cell_pred      = cell_pred,
                                    summarize      = TRUE,
                                    years          = eval(parse(text = year_list)))

  ## Get raking factors (old method)
  rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
                              gaul_list   = gaul_list,
                              rake_to     = gbd)

  ## interpolate raking factors between 5 year bitsb
  for(cc in intersect(unique(rf$name), unique(gbd$name))){
    cc.r <- which(rf$name == cc)

    rf$raking_factor[cc.r[2:5]]   <- approx(x=c(2000, 2005),
                                            y=rf$raking_factor[c(cc.r[1], cc.r[6])],
                                            xout=2001:2004)$y
    rf$raking_factor[cc.r[7:10]]  <- approx(x=c(2005, 2010),
                                            y=rf$raking_factor[c(cc.r[6], cc.r[11])],
                                            xout=2006:2009)$y
    rf$raking_factor[cc.r[12:15]] <- approx(x=c(2010, 2015),
                                            y=rf$raking_factor[c(cc.r[11], cc.r[16])],
                                            xout=2011:2014)$y
    rf$raking_factor[cc.r[17]]    <- 2 * rf$raking_factor[cc.r[16]] - rf$raking_factor[cc.r[15]]
  }

  ## rake cell preds
  raked_cell_pred <- rake_predictions(raking_factors = rf,
                                      pop_wts_object = pop_wts_adm0,
                                      cell_pred      = cell_pred,
                                      logit_rake     = FALSE)

  ## summarize unraked predictions for each cell
  mean_raster <- make_cell_pred_summary( draw_level_cell_pred = cell_pred,
                                        mask                 = simple_raster,
                                        return_as_raster     = TRUE,
                                        summary_stat         = 'median')

  ## make mean raked raster
  raked_mean_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                                              mask                 = simple_raster,
                                              return_as_raster     = TRUE,
                                              summary_stat         = 'median')

  ## make range raked rasterx
  raked_range_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                                               mask                 = simple_raster,
                                               return_as_raster     = TRUE,
                                               summary_stat         = 'cirange')

  raked_lower_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                                               mask                 = simple_raster,
                                               return_as_raster     = TRUE,
                                               summary_stat         = 'lower')

  raked_upper_raster <- make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                                               mask                 = simple_raster,
                                               return_as_raster     = TRUE,
                                               summary_stat         = 'upper')

  ## assign things with region name
  assign(sprintf('%s_rf',reg),rf)
  assign(sprintf('%s_mean_raster',reg),mean_raster)
  assign(sprintf('%s_raked_range_raster',reg),raked_range_raster)
  assign(sprintf('%s_raked_lower_raster',reg),raked_lower_raster)
  assign(sprintf('%s_raked_upper_raster',reg),raked_upper_raster)
  assign(sprintf('%s_raked_mean_raster',reg),raked_mean_raster)

  ## also save raked cell preds
  saveRDS(file = paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group,
                        '/', indicator, '/output/', run_date, "/",
                        reg, "_raked_cell_pred.RDs"),
          object = raked_cell_pred)

  ## remove things
  rm(mean_raster)
  rm(raked_mean_raster)
  rm(raked_lower_raster)
  rm(raked_upper_raster)
  rm(raked_cell_pred)
  rm(cell_pred)
  rm(pop_wts_adm0)
  rm(range_raked_raster)

}


## combine regions raster
rf <- do.call(rbind.fill, list(name_rf,
                               essa_rf,
                               sssa_rf,
                               wssa_rf,
                               cssa_rf))
m = do.call(raster::merge,list(name_mean_raster,
                               essa_mean_raster,
                               sssa_mean_raster,
                               wssa_mean_raster,
                               cssa_mean_raster))
m_raked = do.call(raster::merge,list(name_raked_mean_raster,
                                     essa_raked_mean_raster,
                                     sssa_raked_mean_raster,
                                     wssa_raked_mean_raster,
                                     cssa_raked_mean_raster))
l_raked = do.call(raster::merge,list(name_raked_lower_raster,
                                     essa_raked_lower_raster,
                                     sssa_raked_lower_raster,
                                     wssa_raked_lower_raster,
                                     cssa_raked_lower_raster))
u_raked = do.call(raster::merge,list(name_raked_upper_raster,
                                     essa_raked_upper_raster,
                                     sssa_raked_upper_raster,
                                     wssa_raked_upper_raster,
                                     cssa_raked_upper_raster))
r_raked = do.call(raster::merge,list(name_raked_range_raster,
                                     essa_raked_range_raster,
                                     sssa_raked_range_raster,
                                     wssa_raked_range_raster,
                                     cssa_raked_range_raster))

## save!
save_post_est(rf,'csv','rf')
save_post_est(m,'raster','mean_raster')
save_post_est(m_raked,'raster','mean_raked_raster')
save_post_est(r_raked,'raster','range_raked_raster')
save_post_est(l_raked,'raster','lower_raked_raster')
save_post_est(u_raked,'raster','upper_raked_raster')
## save_post_est(r,'raster','range_raster')


pdf(paste0(sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/',indicator_group, indicator, run_date), indicator, '_mean_plot.pdf'), width = 16, height = 16)
for(i in c(1, 6, 11, 16)){
  raster::plot(m_raked[[i]])
}
dev.off()


## #######################
## setup some locations ##
## #######################
commondir <- sprintf('<<<< FILEPATH REDACTED >>>>>/common_inputs')
mod.dir <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/',
                   indicator_group, indicator, run_date)
si.fig.dir <- paste0(mod.dir, 'si_figs/')
dir.create(si.fig.dir)



## #########################################
## launch aggregation to adm0, adm1, adm2 ##
## #########################################
if(as.numeric(n_ho_folds) > 0){
  holdouts <- 0:as.numeric(n_ho_folds)
}else{
  holdouts <- 0
}

holdouts <- 0 ## i just want to do this for holdout 0

submit_aggregation_script(indicator       = indicator,
                          indicator_group = indicator_group, 
                          run_date        = run_date, 
                          raked           = c(TRUE, FALSE), 
                          pop_measure     = pop_measure, 
                          overwrite       = T, 
                          ages            = 0, # Note: can take vector of ages 
                          holdouts        = holdouts,
                          regions         = strata, 
                          corerepo            = repo, 
                          log_dir         = paste0(sharedir, "/output/", run_date, "/"), 
                          geo_node       = TRUE, 
                          slots           = 8)

waitforaggregation(rd = run_date, indic = indicator, ig = indicator_group,
                   ages     = 0,
                   regions  = strata,
                   holdouts = holdouts,
                   raked    = c(T, F))

combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages     = 0, 
                    regions  = strata,
                    holdouts = holdouts,
                    raked    = c(T, F))

summarize_admins(summstats = c("mean", "upper", "lower", "cirange"), 
                 ad_levels = c(0,1,2), 
                 raked     = c(T, F))


## #######################
## make summary metrics ##
## #######################
run_is_oos <- get_is_oos_draws(ind_gp = indicator_group,
                               ind = indicator,
                               rd = run_date,
                               ind_fm = 'binomial',
                               model_domain = 'africa',
                               age = 0,
                               nperiod = length(eval(parse(text = year_list))),
                               yrs = eval(parse(text = year_list)),
                               get.oos = TRUE,
                               write.to.file = TRUE)
## for admin0
draws.df <- fread(sprintf("<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/output_draws_data.csv",
                          indicator_group, indicator, run_date))
#draws.df$ad0 <- draws.df$country

samples <- length(grep('draw', colnames(draws.df)))

## ad0
ad0.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            aggregate_on='ad0',
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            plot = T, 
                            out.dir = si.fig.dir)
## make coef of var approx
ad0.pvtable$CoV <-  ad0.pvtable$'RMSE' / ad0.pvtable$'Mean\ Pred.' * 100
cols <- names(ad0.pvtable)[-(1:2)]
ad0.pvtable[,(cols) := round(.SD,3), .SDcols=cols]
ad0.pvtable[, c(4, 5)] <- NULL
ad0.pvtable <- ad0.pvtable[, c(1, 2, 5, 3, 4, 8, 6, 7), with = FALSE]
write.csv(ad0.pvtable,
          file = sprintf("%s/ad0_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)

## ad1
ad1.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            aggregate_on='ad1',
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            plot = T, 
                            out.dir = si.fig.dir)
## make coef of var approx
ad1.pvtable$CoV <-  ad1.pvtable$'RMSE' / ad1.pvtable$'Mean\ Pred.' * 100
cols <- names(ad1.pvtable)[-(1:2)]
ad1.pvtable[,(cols) := round(.SD,3), .SDcols=cols]
ad1.pvtable[, c(4, 5)] <- NULL
ad1.pvtable <- ad1.pvtable[, c(1, 2, 5, 3, 4, 8, 6, 7), with = FALSE]
write.csv(ad1.pvtable,
          file = sprintf("%s/ad1_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)

## ad2
ad2.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            aggregate_on='ad2',
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            plot = T, 
                            out.dir = si.fig.dir)
## make coef of var approx
ad2.pvtable$CoV <-  ad2.pvtable$'RMSE' / ad2.pvtable$'Mean\ Pred.' * 100
cols <- names(ad2.pvtable)[-(1:2)]
ad2.pvtable[,(cols) := round(.SD,3), .SDcols=cols]
ad2.pvtable[, c(4, 5)] <- NULL
ad2.pvtable <- ad2.pvtable[, c(1, 2, 5, 3, 4, 8, 6, 7), with = FALSE]
write.csv(ad2.pvtable,
          file = sprintf("%s/ad2_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)

## ho_id (quadtree)
hoid.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            aggregate_on='ho_id',
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            plot = T, 
                            out.dir = si.fig.dir)
## make coef of var approx
hoid.pvtable$CoV <-  hoid.pvtable$'RMSE' / hoid.pvtable$'Mean\ Pred.' * 100
cols <- names(hoid.pvtable)[-(1:2)]
hoid.pvtable[,(cols) := round(.SD,3), .SDcols=cols]
hoid.pvtable[, c(4, 5)] <- NULL
hoid.pvtable <- hoid.pvtable[, c(1, 2, 5, 3, 4, 8, 6, 7), with = FALSE]
write.csv(hoid.pvtable,
          file = sprintf("%s/hoid_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)
