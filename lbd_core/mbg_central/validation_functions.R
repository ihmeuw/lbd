#' @title Get IS and OOS Draws
#' @description This function makes a dataframe of all your data that went into the
#' model and appends columns of draw values from cell_preds
#'
#' @param ind_gp indicator_group
#' @param ind indicator
#' @param rd run date
#' @param ind_fm indicator_family (binomial, gaussian)
#' @param age age group
#' @param nperiod number of periods/years in model
#' @param years vector of years. should be of length nperiod
#' @param write.to.file if true writes final df to standard output dir
#' @param get.oos should we also get OOS extracted values?
#' @param year_col string, name of year column in data frame
#' @param shapefile_version string specifying shapefile version to pull
#' 
#' @return data frame with prediction values at each observation location
#' 
get_is_oos_draws <- function(ind_gp,
                             ind,
                             rd,
                             ind_fm = "binomial",
                             age = 0,
                             nperiod = 16,
                             yrs = 2000:2015,
                             write.to.file = FALSE,
                             get.oos = FALSE,
                             year_col = "original_year",
                             shapefile_version = "current",
                             ...) {

  ###############
  ## Load data ##
  ###############

  message("Load input data used in model")

  # load regions that were used in modeling
  mod.dir <- file.path(fp_list['mbg_root'], ind_gp, ind, 'output', rd)
  all.regions <- get_output_regions(mod.dir)

  # load raw data
  if (!get.oos) {
    df <- fread(file.path(mod.dir, "input_data.csv")) ## raw input data
    df <- merge_with_ihme_loc(df, shapefile_version = shapefile_version)
  } else {
    df <- rbindlist(lapply(readRDS(file.path(mod.dir, "stratum.rds")), data.table))
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
      year_pts <- loc.r[this_year, ]
      #When there is a single point, transpose so that
      #raster::extract processes it correctly
      if (class(year_pts) == 'numeric') {
        year_pts <- t(year_pts)
      }
      df.r[this_year, cell_pred_id := raster::extract(id_raster[[yy]], year_pts)]
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
      load(file.path(mod.dir, sprintf("%s_cell_draws_eb_bin%i_%s_%i.RData", ind, age, rr, this_fold)))

      # extract draws
      df.r[fold == this_fold, paste0("draw", 1:ncol(cell_pred)) := as.data.table(cell_pred[cell_pred_id, ])]
    }

    # return combined draws and draws
    return(df.r)
  }))

  # rename year column back to original
  setnames(df_all, "the_year_col", year_col)

  # save combined data and draws to file and return
  if (write.to.file) write.csv(df_all, file.path(mod.dir, "output_draws_data.csv"), row.names = F)
  return(df_all)
}


#' @title Create scatterplot of Raking Factors
#' @description This function plots a scatterplot of GBD estimates vs MBG estimates
#'
#' @param ind_gp indicator_group
#' @param ind indicator
#' @param rd run date
#' @param output.dir output directory to save files
#' @param shapefile_version string specifying shapefile version to pull
#' @param title title for plots
#' 
#' @return N/A Saves scatterplots in output.dir
#' 

plot.rfs <- function(ind.gp = indicator_group,
                     ind = indicator,
                     rd = run_date,
                     output.dir = si.fig.dir,
                     shapefile_version = 'current',
                     title = "Comparison to GBD 2016 in\n" ## region gets pasted on to this
){

  require(ggrepel) ## for labels

  regions <- get_output_regions(in_dir = paste0(fp_list['mbg_root'],
                                                ind.gp, '/',
                                                ind, '/output/',
                                                rd))

  for(rr in regions){

    ## convert rrs to full names
    if(rr == 'essa') rr_name = "Eastern Sub-Saharan Africa"
    if(rr == 'wssa') rr_name = "Western Sub-Saharan Africa"
    if(rr == 'name') rr_name = "North Africa"
    if(rr == 'sssa') rr_name = "Southern Sub-Saharan Africa"
    if(rr == 'cssa') rr_name = 'Central Sub-Saharan Africa'

    in_dir  <- paste0(fp_list['mbg_root'], ind.gp, '/', ind, '/output/', rd)
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
  ggsave(filename = sprintf('%s%s_all_rfs_labs.png',
                            output.dir, ind),
         all.rfs, width = 12, height = 16)
}

#' @title Make table of INLA model results
#' @description This function takes in standard model run inputs and outputs a table 
#' of fixed effect,spatio-temporal hyperparameter, and random effects parameter
#' summaries. Note: takes a little while since it has to recreate the
#' SPDE INLA object.
#'
#' @param regions regions
#' @param rd run date
#' @param holdout holdout group
#' @param age age group
#' @param ind indicator
#' @param ind_gp indicator group
#' @param sharedir filepath were model runs are saved
#' 
#' @return list of tables (each entry is a region) containing the model summaries
#' 

model_fit_table_list <- function(regions, 
                                 rd=run_date, 
                                 holdout = 0,
                                 age = 0,
                                 ind= indicator,
                                 ind_gp = indicator_group,
                                 sharedir = file.path(fp_list['mbg_root'], indicator_group, indicator)){
  ## load models
  require(INLA)
  message(sprintf('Pulling together results for %s models',rd))

  tlist=list()

  for(rr in regions){
    message(sprintf('::on region %s',rr))
    reg  <-  rr

    message("::::loading in pre-INLA objects to get spde")
    pathaddin  <-  paste0('_bin',age,'_',rr,'_',holdout)
    load(paste0(fp_list['mbg_root'], ind_gp, '/',
                ind, '/model_image_history/', rd, pathaddin,
                '.RData'))

    modnames = c('gam','gbm','ridge','enet','lasso')

    full_raster_list <- cov_list

    ## hotfix!! inv.logit ## TODO ## TODO don't need this if either
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
                                       use_ctry_res = use_country_res,
                                       use_nugget = use_inla_nugget)

    spde <- input_data[[2]]
    ## this is what we neede!

    message('::::loading in INLA fit\n')
    f <-  sprintf('%s/output/%s/inla_model_fit_pre_preds_%s_holdout_%i.RDS',
                  sharedir,rd,reg, holdout)
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

#' @title Plot Absolute Residual Errors
#' @description This function plots residual errors (i.e. observed value - predicted value)
#'
#' @param gaul_list list of GAUL codesfor the entire modeling domain
#' @param df data frame that is output from `get_is_oos_draws()`
#' @param sample string indicating whether to whether to do IS, OOS, of BOTH
#' @param subset_shape spatial polygon of entirme modeling domain
#' @param ind indicator
#' @param ind_gp indicator group
#' @param rd run date
#' @param save.dir filepath to save plots
#' @param year_col string indicating year column in df
#' 
#' @return Plots showing the difference between observed and predicted values
#' 

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
#' draws.df <- fread(file.path(
#'   fp_list['mbg_root'], indicator_group, indicator, 'output',
#'   run_date, 'output_draws_data.csv'
#' ))
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
#' @import boot
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
    message("NICK: THIS SHOULD WORK IF YOU HAVE DRAWS OF PRECISION NAMED tau_1,...tau_1000 BUT THIS IS UNTESTED")
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
    weighted.rmse <- function(error, w) {
      sqrt(sum((w / sum(w)) * ((error) ^ 2)))
    }
    if (weighted) res$weight <- res$exposure else res$weight <- 1
    res2 <-
      res[, c(lapply(.SD, function(x) weighted.mean(x, weight)),
                    rmse = weighted.rmse(error, weight),
                    median_SS = median(exposure),
                    cor = boot::corr(cbind(p, mean_draw), weight)),
          by = result_agg_over,
          .SDcols = c("error", "abs_error", "mean_draw", "p", paste0("clusters_covered_", coverage_probs))]

    setnames(res2, c("error", "abs_error", "p", paste0("clusters_covered_", coverage_probs)),
             c("me", "mae", "mean_p", paste0("coverage_", coverage_probs)))

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
      png(file.path(out.dir, plot_filename), width = 12, height = 12, units = "in", res = 350)

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
        png(file.path(out.dir, cplot_filename), width = 12, height = 12, units = "in", res = 350)

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
      filename <- file.path(out.dir, sprintf("%s_metrics%s.csv", agg_on, a_pathaddin))
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
#' @author Lauren Woyczynski, \email{lpw3@uw.edu}
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
#' @param input_file If your team does not use the standard `<mbg_root>/input_data/{indicator}.csv` input data filepath, supply the filepath to your data here
#' @param admin0_file If your team does not use the standard `<mbg_root>',indicator_group, '/', indicator, '/output/',run_date, '/pred_derivatives/admin_summaries/',indicator, '_admin_0_unraked_summary.csv'` directory to save unraked aggregated admin summaries, supply the correct filepath here
#' @param admin1_file If your team does not use the standard `<mbg_root>,indicator_group, '/', indicator, '/output/',run_date, '/pred_derivatives/admin_summaries/',indicator, '_admin_1_unraked_summary.csv'` directory to save unraked aggregated admin summaries, supply the correct filepath here
#' @param admin2_file If your team does not use the standard `<mbg_root>,indicator_group, '/', indicator, '/output/',run_date, '/pred_derivatives/admin_summaries/',indicator, '_admin_2_unraked_summary.csv'` directory to save unraked aggregated admin summaries, supply the correct filepath here
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
    if(!file.exists(paste0(fp_list['mbg_root'], '/input_data/', indicator, '.csv'))) stop('Indicator input data does not exist in /share - please specify a file path with the input_file arg')
    input.dt <- fread(paste0(fp_list['mbg_root'], '/input_data/', indicator, '.csv'))
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
    ad0_mbg <- fread(paste0(fp_list['mbg_root'],
                            indicator_group, '/',
                            indicator, '/output/',
                            run_date, '/pred_derivatives/admin_summaries/',
                            indicator, '_admin_0_unraked_summary.csv'))
  } else{
    ad0_mbg <- fread(admin0_file)
  }
  if(is.null(admin1_file)){
    ad1_mbg <- fread(paste0(fp_list['mbg_root'],
                            indicator_group, '/',
                            indicator, '/output/',
                            run_date, '/pred_derivatives/admin_summaries/',
                            indicator, '_admin_1_unraked_summary.csv'))
  } else{
    ad1_mbg <- fread(admin1_file)
  }
  if(is.null(admin2_file)){
    ad2_mbg <- fread(paste0(fp_list['mbg_root'],
                            indicator_group, '/',
                            indicator, '/output/',
                            run_date, '/pred_derivatives/admin_summaries/',
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


#' @title Plot Hyperparameters
#' @description This function creates a plot of hyperparameters for a model run. The plots show the prior and posterior distributions by region for each hyperparameter present in the model.
#'
#' @param indicator indicator name
#' @param indicator_group indicator group
#' @param run_date  model run date
#' @param age age as specified in loopvars of modelling run
#' @param holdout holdout number as specified in loopvars of modelling run
#' @param save_file default NULL, if NULL creates an `INLA_hyperparameters.pdf` file in the run date folder. Else takes a file path of where to save the plot to. 
#' @param regions default Regions in global env. A vector of region names used in the model run.
#'
#' @return Null
plot_hyperparameters <- function(indicator, 
                                 indicator_group,
                                 run_date, 
                                 age, 
                                 holdout, 
                                 save_file = NULL,
                                 regions = Regions) {
  
  # get regions
  run_dir <- paste0(fp_list$mbg_root, indicator_group, '/', indicator, '/output/', run_date, "/")
  config <- fread(paste0(run_dir, 'config.csv'))
  
  # extract prior and posterior distributions from INLA model objects
  message("Load models & extract priors and posteriors for hyper-parameters")
  dist <- rbindlist(lapply(regions, function(r) {
    
    # load model
    message(paste0('...', r))
    load(paste0(run_dir, indicator, "_model_eb_bin", age, "_", r, "_", holdout, ".RData"))
    
    # extract hyper-priors from INLA (based on plot.inla() code)
    all.hyper <- INLA:::inla.all.hyper.postprocess(res_fit$all.hyper)
    hyper <- res_fit$marginals.hyperpar
    id <- strsplit(sapply(hyper, attr, 'hyperid'), split = '\\|')
    prior <- rbindlist(lapply(names(id), function(x) {
      print(x)
      if (grepl("Theta. for", x)) range <- c(-5, 5)
      # if (grepl("Range for", x)) range <- c(1, 1000)
      # if (grepl("Stdev for", x)) range <- c(1, 1000)
      if (grepl("GroupRho for", x)) range <- c(-0.999, 0.999)
      if (grepl("Group PACF. for", x)) range <- c(-0.999, 0.999)
      if (grepl("Precision for", x)) range <- c(1, 1000)
      p <- INLA:::inla.get.prior.xy(section = tolower(id[[x]][2]), hyperid = id[[x]][1], all.hyper = all.hyper, range = range, intern = F)
      if (grepl("Precision for", x)) {
        p <- inla.tmarginal(function(x) sqrt(1/x), p, method = 'linear')
        x <- gsub("Precision for", "SD for", x)
      }
      data.table(region = r, type = 'prior', name = x, x = p$x, y = p$y)
    }))
    
    # extract corresponding posteriors from INLA
    post <- rbindlist(lapply(names(hyper), function(x) {
      p <- hyper[[x]]
      if (grepl("Precision for", x)) {
        p <- inla.tmarginal(function(x) sqrt(1/x), p, method = 'linear')
        x <- gsub("Precision for", "SD for", x)
      }
      # if (x == "Range for space") {
      #   x <- "Nominal range"
      #   p[, 'x'] <- sqrt(8) / exp(p[, 'x'])
      # }
      # if (x == "Stdev for space") {
      #   x <- "Nominal variance"
      #   t1 <- hyper[["Theta1 for space"]]
      # }
      data.table(region = r, type = 'posterior', name = x, x = p[, 'x'], y = p[, 'y'])
    }))
    
    # combine
    all <- rbind(prior, post)
    all[, name := factor(name, unique(name))]
    return(all)
  }))
  
  # make plots
  message("Plotting hyper-parameters")
  
  if (is.null(save_file)) save_file <- paste0(run_dir, "/inla_hyperparameters.pdf")
  pdf(save_file, width = 14, height = 8)
  gg <- ggplot(dist[y > 1e-8,], aes(x = x, y = y, color = region, linetype = type)) +
    facet_wrap(~ name, scales = "free") +
    geom_line() +
    labs(x = '', y = '', title = "Hyper-parameter prior and posterior distributions") +
    theme_bw()
  print(gg)
  dev.off()
  
  return(dist)
}