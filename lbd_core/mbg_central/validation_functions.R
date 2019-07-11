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
##   model_domain: larger domain you modelled over 
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
                             ind_fm = 'binomial',
                             model_domain = 'africa', 
                             age = 0,
                             nperiod = 16,
                             yrs = 2000:2015,
                             write.to.file = FALSE,
                             get.oos = FALSE,
                             year_col = 'original_year',
                             shapefile_version = 'current'){

  ###############
  ## Load data ##
  ###############

  message("Load input data used in model")

  # load regions that were used in modeling
  mod.dir <- "<<<< FILEPATH REDACTED >>>>"
  all.regions <- get_output_regions(mod.dir)

  # load raw data
  if(!get.oos) {
    df <- fread(paste0(mod.dir, 'input_data.csv')) ## raw input data
    df <- merge_with_ihme_loc(df, shapefile_version = shapefile_version)
  } else {
    df <- rbindlist(lapply(readRDS(paste0(mod.dir, 'stratum.rds')), data.table))
  }

  # rename year column for convenience
  setnames(df, year_col, "the_year_col")

  ###################
  ## Assign admins ##
  ###################

  message('Identify ad1 and ad2 membership')

  # load admin2 shapefile (also contains admin1)
  admin2_shapefile <- rgdal::readOGR(get_admin_shapefile(admin_level=2, version = shapefile_version))
  for (v in grep("CODE", names(admin2_shapefile@data))) admin2_shapefile@data[[v]] <- as.numeric(as.character(admin2_shapefile@data[[v]]))

  # identify the admin2 (and by extension, the admin1) each point belongs to
  locs <- SpatialPoints(cbind(df$longitude, df$latitude), proj4string = CRS(proj4string(admin2_shapefile)))
  adm.df <- sp::over(locs, admin2_shapefile)

  # for those that don't fall inside a polygon, assign them to the nearest polygon (this tends to happen on coastlines)
  # do this by country to speed things up and to make sure the point ends up at least in the correct admin0
  for (ctry in unique(df[is.na(adm.df[,1]), GAUL_CODE])) {
    ii <- which(is.na(adm.df[,1]) & df$GAUL_CODE == ctry)
    temp_shape <- admin2_shapefile[admin2_shapefile@data$ADM0_CODE == ctry,]
    distmat <- gDistance(locs[ii], temp_shape, byid = T)
    jj <- apply(distmat, 2, which.min)
    adm.df[ii,] <- temp_shape@data[jj,]
    rm(ii, jj, temp_shape, distmat)
  }

  # copy admin information into df
  df$ad1 <- adm.df$ADM1_CODE
  df$ad2 <- adm.df$ADM2_CODE
  df$ad0 <- df$GAUL_CODE # for consistency

  ###############
  ## Get draws ##
  ###############

  message('Get draws')

  # loop over regions
  df_all <- rbindlist(lapply(all.regions, function(rr) {

    message(paste('...Region:', rr))

    # load the simple raster template
    message('......load simple raster template')
    gaul_list <- get_adm0_codes(rr, shapefile_version = shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 0.4, subset_only = TRUE,
                                               shapefile_version = shapefile_version)
    subset_shape <- simple_polygon_list[[1]]
    raster_list <- build_simple_raster_pop(subset_shape)
    template <- raster_list[['simple_raster']]

    # subset data to this region
    df.r <- df[region == rr,]
    loc.r <- as.matrix(df.r[, list(longitude, latitude)])

    # 'nudge' coordinates to make them fall inside the nearest cell non-NA cell in template if they
    # are in an NA cell (this can happen near coastlines)
    check <- raster::extract(template, loc.r)
    miss <- which(is.na(check))
    for (ii in miss) {
      dist <- replace(distanceFromPoints(template, loc.r[ii,]), is.na(template), NA)
      loc.r[ii,] <- xyFromCell(template, which.min(dist)[1])
    }

    # make a raster of cell_pred row ids
    id_raster <- insertRaster(template, matrix(1:(length(cellIdx(template)) * nperiod), ncol = nperiod))

    # get cell_pred row ids for each data point
    for (yy in 1:nperiod) {
      this_year <- which(df.r$the_year_col == yrs[yy])
      if (length(this_year) == 0) next
      df.r[this_year, cell_pred_id := raster::extract(id_raster[[yy]], loc.r[this_year,])]
    }

    # if out of sample metrics are requested, duplicate the data to create separate rows for in and out of sample
    if (get.oos) {
      df.r <- rbind(df.r, cbind(df.r[, -"fold", with=F], fold = 0), use.names = T)
    } else {
      df.r[, fold := 0]
    }

    # loop over holdouts
    for (this_fold in sort(unique(df.r$fold))) {

      # load cell pred objects
      message(paste('......load cell preds for holdout', this_fold))
      load(paste0(mod.dir, sprintf('%s_cell_draws_eb_bin%i_%s_%i.RData', ind, age, rr, this_fold)))

      # extract draws
      df.r[fold == this_fold, paste0("draw", 1:ncol(cell_pred)) := as.data.table(cell_pred[cell_pred_id,])]

    }

    # return combined draws and draws
    return(df.r)

  }))

  # rename year column back to original
  setnames(df_all, "the_year_col", year_col)

  # save combined data and draws to file and return
  if (write.to.file) write.csv(df_all, paste0(mod.dir, 'output_draws_data.csv'), row.names = F)
  return(df_all)

}

####################################################################################################################
## INPUT: # data file matching <<<< FILEPATH REDACTED >>>>/output_draws_data.csv
## OUTPUT SUMMARIZED PREDICTIVE VALIDITY METRICS
##
####################################################################################################################

## get_pv_table ################################################

#' Generate plots and tables of metrics of predictive validity
#'
#' @param d data.table containing data and IS/OOS draws from `run_is_oos()`
#' @param indicator indicator
#' @param indicator_group indicator_group
#' @param rd run date
#' @param aggregate_on column of spatial aggregation ("country", "ad1", or "ad2" usually) - can pass vector
#' @param draws number of draws (numeric)
#' @param coverage_probs probability at which you want to assess coverage of predictive intervals
#'                       can be a single number or a vector; if a vector can produce calibration plots
#'                       assessing x% coverage at various cutoffs
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
#' @param plot_ci_level what ci level to plot (numeric, i.e. "95" is the 95% UI)
#' @param ci_color what color do you want the CI lines to be?
#' @param point_alpha how transparent do you want the points themselves to be (numeric, [0-1])
#' @param point_color what color do you want the points to be in your validation plots?
#' @param plot_title main title for your plots (defaults to "indicator")
#' @param plot_ncol how many columns do you want your plots to have?
#' @param save_csv save a csv table of your predictive validity metrics?
#' @param out.dir where do you want these results to be saved?
#'
#' @return data.table of PV metrics
#' @export plots and tables of PV metrics depending on the combinations of options above
#'
#' @examples
#'
#' # A PV table by modeling regions with plots (one for each modeling region) and semi-transparent points
#' # Also includes calibration plots since there are multiple coverage levels
#'
#' # Get in and out of sample draws
#' run_in_oos <- get_is_oos_draws(ind_gp = indicator_group,
#'                                ind = indicator,
#'                                rd = run_date,
#'                                ind_fm = 'binomial',
#'                                model_domain = Regions,
#'                                age = 0,
#'                                nperiod = length(year_list),
#'                                yrs = year_list,
#'                                get.oos = as.logical(makeholdouts),
#'                                write.to.file = TRUE,
#'                                year_col = "year")
#'
#' # Load the IS/OOS draws created above
#' draws.df <- fread(sprintf("<<<< FILEPATH REDACTED >>>>/output_draws_data.csv"))
#'
#' # run the PV table function
#' pvtable.reg <- get_pv_table(d = draws.df,
#'                             indicator_group = indicator_group,
#'                             rd = run_date,
#'                             indicator=indicator,
#'                             result_agg_over = c("year", "oos", "region"),
#'                             coverage_probs = seq(from=5,to=95,by=5),
#'                             aggregate_on='ad2',
#'                             draws = as.numeric(samples),
#'                             out.dir = out_dir,
#'                             plot = TRUE,
#'                             plot_by = "region",
#'                             plot_by_title = "Region",
#'                             plot_ci = TRUE,
#'                             point_alpha = 0.1,
#'                             point_color = "black",
#'                             ci_color = "gray",
#'                             plot_title = plot_title)

get_pv_table <- function(d,
                        indicator,
                        indicator_group,
                        rd,
                        aggregate_on,
                        draws = 1000,
                        coverage_probs = c(95),
                        result_agg_over =  c('year','oos'),
                        weighted = TRUE,
                        family = 'binomial',
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

  require(boot)
  require(ggplot2)
  require(tidyr)
  str_match <- stringr::str_match

  d <- data.table(d)

  if (!is.null(plot_by)) {
    if (!(plot_by %in% result_agg_over)) {
      stop("If you specify `plot_by`, you must also include that item in `result_agg_over`")
    }
  }

  ## Get binomial predictions
  message('Simulating predictive draws')
  if (family == 'binomial') {
    x <- sapply(1:draws, function(i) rbinom(nrow(d), round(d$N), d[[paste0('draw', i)]]))
    d[, Y := get(indicator) * round(N)/N] # adjust observed number of cases for effect of rounding N
  }
  if (family == 'gaussian') {
    x <- sapply(1:draws, function(i) rnorm(nrow(d), sqrt(d[[paste0('tau_', i)]]), d[[paste0('draw', i)]]))
    d[, Y := get(indicator)]
  }

  message(paste0('...NOTE: missing predictions for ', sum(is.na(x[,1])), ' of ', nrow(x), ' (', round(100*mean(is.na(x[,1]))), '%) data points.'))

  ## Get coverage for each data point
  message('Calculate coverage')
  one_side <- (1 - coverage_probs/100)/2
  ui <- t(apply(x, 1, quantile, c(one_side, 1 - one_side), na.rm = T))
  for (c in coverage_probs) {
    c_id <- which(c == coverage_probs)
    d[, paste0('clusters_covered_', c) := Y %between% list(ui[, c_id], ui[, c_id + length(coverage_probs)])]
  }

  ## Collapse data and draws spatially
  message('Collapse data and draws spatially')
  d[, oos := (fold != 0)]
  d[, p := get(indicator) / N]
  d[, exposure := N * weight]

  # Create list of pv metrics for all levels of aggregation in `aggregate_on`
  message('Outputting predictive validity metrics for your models at each level of aggregation')
  pv_table_list <- lapply(aggregate_on, function(agg_on) {
    message(paste0("  ", agg_on))
    by_vars <- unique(c('oos', 'year', result_agg_over, agg_on))
    collapse_vars <- c('p', paste0('clusters_covered_', coverage_probs), paste0('draw', 1:draws))

    res <- d[!is.na(draw1),
             c(list(total_clusters = .N, exposure = sum(exposure)),
               lapply(.SD, function(x) weighted.mean(x, exposure, na.rm=T))),
             keyby = by_vars, .SDcols = collapse_vars]
    res[, mean_draw := rowMeans(.SD), .SDcols = paste0('draw', 1:draws)]
    res[, error := p - mean_draw]
    res[, abs_error := abs(error)]

    ## Collapse to calculate predictive validity metrics
    weighted.rmse <- function(error, w) sqrt( sum(  (w/sum(w)) * ((error)^2) ) )
    if (weighted) res$weight <- res$exposure else res$weight <- 1
    res2 <- res[, c(lapply(.SD, function(x) weighted.mean(x, weight)),
                    rmse = weighted.rmse(error, weight),
                    median_SS = median(exposure),
                    cor = corr(cbind(p, mean_draw), weight)),
                by = result_agg_over,
                .SDcols = c("error", "abs_error", "mean_draw", "p", paste0("clusters_covered_", coverage_probs))]
    setnames(res2, c("error", "abs_error", "p", paste0("clusters_covered_", coverage_probs)),
             c("me", "mae", "mean_p", paste0("coverage_", coverage_probs)))

    return(list(res = res, res2 = res2))
  })

  names(pv_table_list) <- aggregate_on

  # Make plots
  if(plot==TRUE){
    message('Making plots of aggregated data and estimates')

    # Get unique levels of `plot_by` and set up a plot_by title if needed
    if (!is.null(plot_by)) {
      plot_by_levels <- unique(pv_table_list[[1]][["res"]][, get(plot_by)])
    }  else {
      plot_by_levels <- NULL
    }

    if (is.null(plot_by_title)) plot_by_title <- plot_by

    # Create a table of things to plot
    plot_table <- CJ(aggregate_on = aggregate_on,
                     oos = unique(pv_table_list[[1]][["res"]]$oos),
                     plot_by_value = if (is.null(plot_by_levels)) NA else plot_by_levels)

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
      plot_filename <- paste0(indicator, '_validation_plot_',
                              paste(c(as.character(agg_on),
                                      setdiff(result_agg_over, "oos")),
                                    collapse="_"), "_",
                              ifelse(oosindic, "OOS", "IS"),
                              ifelse(is.na(pb_val), "", paste0("_", pb_val)),
                              '.png')
      message(paste('    ', plot_filename))
      png(paste0(out.dir, plot_filename), width = 12, height = 12, units = 'in', res = 350)

      # Subset data
      fdata <- res[oos == oosindic,]
      if (!is.na(pb_val)) {
        setnames(fdata, plot_by, "plot_by_column") # convenience
        fdata <- fdata[plot_by_column == pb_val,]
      }

      # Set up CI bar limits; range as defaults
      if (plot_ci) {
        fdata[, upper := apply(.SD, 1, quantile, p = 0.01*(plot_ci_level + (100 - plot_ci_level)/2), rm.na=TRUE), .SDcols = paste0('draw', 1:draws)]
        fdata[, lower := apply(.SD, 1, quantile, p = 0.01*((100 - plot_ci_level)/2), rm.na=TRUE), .SDcols = paste0('draw', 1:draws)]
        limits <- fdata[, range(c(p, mean_draw, lower, upper))]
      } else {
        limits <- fdata[, range(c(p, mean_draw))]
      }

      # The plot code itself
      gg <- ggplot(fdata, aes(x = p, y = mean_draw, size = weight)) +
        geom_abline(intercept=0, slope=1, color = 'red') +
        geom_point(colour = point_color, alpha = point_alpha) +
        scale_size_area() +
        xlim(limits) +
        ylim(limits) +
        coord_equal() +
        theme_bw() +
        theme(strip.background = element_rect(fill="white")) +
        labs(x = 'Data Estimate',
             y = 'Mean Prediction',
             size = "Weight",
             title = paste0("Validation Plot for ", plot_title, " by ", agg_title),
             subtitle = paste0("OOS: ", oosindic, ifelse(is.na(pb_val), "", paste0(" | ", plot_by_title, ": ", pb_val))))

      if(plot_ci) {
        gg <- gg +  geom_errorbar(aes(ymin = lower, ymax = upper), colour = point_color, width = 0, size = .3, alpha = min(point_alpha, 0.2))
      }
      if (length(setdiff(result_agg_over, 'oos')) > 0) {
        gg <- gg + facet_wrap(as.formula(paste("~", paste(setdiff(result_agg_over, c('oos', plot_by)), collapse = "+"))),
                              ncol = plot_ncol)
      }

      plot(gg)
      dev.off()

      # Make a calibration plot -----------------------
      if (plot == T & length(coverage_probs) > 1) {

        # Set up a subset of the data for plotting
        fdata <- res2[, unique(c("oos", result_agg_over, paste0("coverage_", coverage_probs))), with=F]
        fdata <- fdata[oos == oosindic]
        if (!is.na(pb_val)) {
          setnames(fdata, plot_by, "plot_by_column") # convenience
          fdata <- fdata[plot_by_column == pb_val,]
        }
        fdata <- melt(fdata, id.vars = names(fdata)[!grepl("coverage", names(fdata))],
                      value.name = "observed_coverage",
                      variable.name = "coverage")
        fdata[, expected_coverage := as.numeric(gsub("coverage_", "", coverage))]
        fdata[, observed_coverage := observed_coverage * 100]
        fdata$group <- apply(fdata[, result_agg_over[!(result_agg_over %in% c('oos', 'plot_by_column', plot_by))], with=F], 1, paste, collapse=" ")
        if (sum(!is.na(fdata$group)) == 0) fdata[, group := "All"]

        # Create filename

        cplot_filename <- paste0(indicator, '_calibration_plot_',
                                 paste(c(as.character(agg_on),
                                         setdiff(result_agg_over, "oos")),
                                       collapse="_"), "_",
                                 ifelse(oosindic, "OOS", "IS"),
                                 ifelse(is.na(pb_val), "", paste0("_", pb_val)),
                                 '.png')
        message(paste('    ', cplot_filename))
        png(paste0(out.dir, cplot_filename), width = 12, height = 12, units = 'in', res = 350)

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
          theme_bw()+
          labs(x = "Expected coverage", y = "Observed coverage")

        plot(gg)
        dev.off()

      } # End calibration plot if statement
    }# End plot table loop
  } # End if (plot==T) loop

  # Format, save (if desired), and return `res2` objects
  output_list <- lapply(aggregate_on, function(agg_on) {
    res2 <- copy(pv_table_list[[agg_on]][["res2"]])
    setorderv(res2, result_agg_over)
    setnames(res2, result_agg_over, ifelse(result_agg_over == 'oos', 'OOS', gsub("(^.)", "\\U\\1", result_agg_over, perl=T)))
    setnames(res2, c("me", "mae", "mean_draw", "mean_p", "rmse", "median_SS", "cor", paste0("coverage_", coverage_probs)),
             c("Mean Err.", "Mean Abs. Err.", "Mean Pred.", "Mean Obs.", "RMSE", "Median SS", "Corr.", paste0(coverage_probs, "% Cov.")))

    # Save final tables if desired
    if (save_csv) {
      a_pathaddin <- ifelse(length(setdiff(result_agg_over, 'oos')) > 0,
                            paste0("_by_", paste0(setdiff(result_agg_over, 'oos'), collapse = "_")),
                            "")
      filename <- paste0(out.dir, agg_on, "_metrics", a_pathaddin, ".csv")

      message(paste0("Saving csv to ", filename, "..."))
      write.csv(res2, file = filename)
    }
    return(res2)
  })

  names(output_list) <- aggregate_on

  return(output_list)
}
