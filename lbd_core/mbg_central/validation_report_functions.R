## Function to run each country in a region
report_each_gaul <- function(gaul, quilt_plot_list, target_plot_list, shapefile_version = 'current') {

  message(paste0('REPORT FOR ', gaul, ':'))

  ##################################################################################
  ## Face validity checks
  ##################################################################################
  ##  - Raking factors
  ##  - DHS admin1 comparison
  ##  - Semivariograms
  ##  - Time series, GBD vs. MBG
  ##################################################################################
  ## Plot 2000, 2005, 2010, 2015 predictions, raked and unraked
  pred_maps <- extract_raked_and_unraked(gaul,
                                         results,
                                         results_raked,
                                         master_list,
                                         admin0 = admin0)

  ## Plot raking factors (all + by country)
  country_rf_plot <- val_raking_factors(gaul_list = gaul,
                                        rfs = all_rfs)

  ## DHS admin1 validation (all + by country)

  ## (???) Comparison to some threshold (full Bayesian draw-level)

  ## Semivariograms of raw data and predictions at raw data locations (by country)
  vario_and_violin_plots <- plot_varios_and_violins(indicator = indicator,
                                                    indicator_group = indicator_group,
                                                    run_date = run_date,
                                                    this_gaul = gaul)

  ## GBD and MBG admin0 estimates by year with credible intervals
  country_time_series_plots <- summarize_admin2(gaul = gaul,
                                                indicator = indicator,
                                                indicator_group = indicator_group,
                                                run_date = run_date,
                                                nperiod = total_periods,
                                                master_list = master_list)

  if(indicator == 'edu_mean') {

    ## Load premade .csv of DHS admin1 aggregates for edu_mean and an SPDF with all the relevant location_code/shapefile entries.
    load('/share/geospatial/mbg/validation_polygons/dhs_admin1.RData')
    dhs_admin1_data <- fread('/share/geospatial/mbg/validation_data/edu_mean_dhs_admin1.csv')

    ## Make plot for each NID in this country
    loc_names <- get_location_code_mapping(shapefile_version = shapefile_version)
    this_iso3 <- loc_names[GAUL_CODE %in% gaul, ihme_lc_id]
    dhs_admin1_data <- dhs_admin1_data[iso3 == this_iso3 & year >= 2000, ]
    iso3_nid_list <- unique(dhs_admin1_data[, nid])

    ## Plot comparison
    plot_each_nid <- function(this_nid) {
      unraked_dhs_plot <- compare_admin_estimates(nid_list = this_nid,
                                                  results_raster = results,
                                                  compare_source_title = 'DHS Admin1',
                                                  raked_title = 'unraked',
                                                  outcome_title = 'Mean years',
                                                  compare_data = dhs_admin1_data,
                                                  compare_spdf = poly_shapes_all,
                                                  master_list = master_list)

      return(unraked_dhs_plot)
    }
    unraked_dhs_plots <- lapply(iso3_nid_list, plot_each_nid)

    plot_each_nid <- function(this_nid) {
      raked_dhs_plot <- compare_admin_estimates(nid_list = this_nid,
                                                results_raster = results_raked,
                                                compare_source_title = 'DHS Admin1',
                                                raked_title = 'raked',
                                                outcome_title = 'Mean years',
                                                compare_data = dhs_admin1_data,
                                                compare_spdf = poly_shapes_all,
                                                master_list = master_list)
      return(raked_dhs_plot)
    }
    raked_dhs_plots <- lapply(iso3_nid_list, plot_each_nid)

  }

  ##################################################################################
  ## Statistical validity
  ##################################################################################
  ##  - Covariate relative importance
  ##  - In-sample and out-of-sample error (quiltplots, tables)
  ##################################################################################
  message('VALIDATION STEP 2: Running statistical validity functions...')

  ## Covariate importance

  ## Quiltplots (bias, coverage, RMSE), also summarize into tabulated tables

  ##################################################################################
  ## Make reports
  ##################################################################################
  ##  - All gauls
  ##  - By country
  ##################################################################################
  message('Making final reports for all countries and each country individually...')

  # plot(vario_and_violin_plots[[2]], main="PACF")
  # grid.echo()
  # pacf <- grid.grab()
  # pacf$vp <- viewport(layout.pos.row = 11:14, layout.pos.col = 9:16)

  gLegend<-function(a.gplot){
    pdf(NULL) # Workaround for bug in ggplot_gtable causing empty Rplots.pdf to be created
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    graphics.off()
    return(legend)
  }
  
  # grab your legends using the predefined functions, then state their grid location
  map.legend <- gLegend(pred_maps[['raked']])
  map.legend$vp <- viewport(layout.pos.row = 1:10, layout.pos.col = 17:18)

  # Initialize plot with master title
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(28, 18)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  # Plot all pred maps
  grid.text("Unraked", vp = vplayout(1,1:16), gp=gpar(fontsize=20, fontface='bold'))
  # print(pred_maps[['unraked']][[1]] + theme(legend.position="none"), vp = vplayout(2:5, 1:4))
  # print(pred_maps[['unraked']][[2]] + theme(legend.position="none"), vp = vplayout(2:5, 5:8))
  # print(pred_maps[['unraked']][[3]] + theme(legend.position="none"), vp = vplayout(2:5, 9:12))
  # print(pred_maps[['unraked']][[4]] + theme(legend.position="none"), vp = vplayout(2:5, 13:16))
  print(pred_maps[['unraked']] + theme(legend.position="none"), vp = vplayout(2:7, 1:16))
  grid.text("Raked", vp = vplayout(8,1:16), gp=gpar(fontsize=20, fontface='bold'))
  # print(pred_maps[['raked']][[1]] + theme(legend.position="none"), vp = vplayout(7:10, 1:4))
  # print(pred_maps[['raked']][[2]] + theme(legend.position="none"), vp = vplayout(7:10, 5:8))
  # print(pred_maps[['raked']][[3]] + theme(legend.position="none"), vp = vplayout(7:10, 9:12))
  # print(pred_maps[['raked']][[4]] + theme(legend.position="none"), vp = vplayout(7:10, 13:16))
  print(pred_maps[['raked']] + theme(legend.position="none"), vp = vplayout(9:14, 1:16))
  grid.draw(map.legend)
  # Plot vario
  #print(vario_and_violin_plots[[1]] + theme(legend.position="none"), vp = vplayout(11:14, 1:16))
  #grid.draw(pacf)
  # Plot raking factors and violins
  print(country_rf_plot + theme(legend.position="none"), vp = vplayout(15:21, 1:8))
  print(vario_and_violin_plots[[3]] + theme(legend.position="none"), vp = vplayout(15:21, 9:16))
  # Admin plots
  print(country_time_series_plots, vp = vplayout(22:28, 1:16))

  ## Make new page with threshold plot and covariate importance
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(28, 18)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  this_target_plot <- target_plot_list[[as.character(gaul)]]
  this_target_plot$vp <- viewport(layout.pos.row = 1:16, layout.pos.col = 1:18)
  grid.draw(this_target_plot)

  ## Make new page with quilt plots, if they exist
  if(!is.null(quilt_plot_list[[as.character(gaul)]])) {
    quilt.legend <- gLegend(quilt_plot_list[[as.character(gaul)]])
    quilt.legend$vp <- viewport(layout.pos.row = 1:15, layout.pos.col = 17:18)
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(28, 18)))
    vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
    # Plot all quilt maps for years that exist
    grid.text("In-Sample", vp = vplayout(1,1:16), gp=gpar(fontsize=20, fontface='bold'))
    print(quilt_plot_list[[as.character(gaul)]] + theme(legend.position="none"), vp = vplayout(2:16, 1:15))
    grid.draw(quilt.legend)
  }

  ## Plot DHS admin1 comparison for education
  if(indicator == 'edu_mean') {
    if(length(unraked_dhs_plots)!=0) {
      for(i in 1:length(unraked_dhs_plots)) {

        grid.newpage()
        pushViewport(viewport(layout = grid.layout(28, 18)))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        print(unraked_dhs_plots[[i]], vp = vplayout(1:13, 1:16))
        print(raked_dhs_plots[[i]], vp = vplayout(15:27, 1:16))

      }
    }
  }

  return(NULL)

}

extract_raked_and_unraked <- function(gaul,
                                      results,
                                      results_raked,
                                      master_list,
                                      admin0,
                                      load_simple_raster = FALSE,
                                      shapefile_version = 'current') {

  if(load_simple_raster==TRUE) {
    gaul_list <- gaul
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                               buffer = 0.4,
                                               subset_only = T,
                                               shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
  }
  if(load_simple_raster==FALSE) {
    simple_raster <- master_list[[paste0('list_', gaul, '.simple_', gaul)]]
  }
  results_subset  <- crop(results, extent(simple_raster))
  results_subset  <- setExtent(results_subset, simple_raster)
  results_subset  <- mask(results_subset, simple_raster)

  raked_results_subset  <- crop(results_raked, extent(simple_raster))
  raked_results_subset  <- setExtent(raked_results_subset, simple_raster)
  raked_results_subset  <- mask(raked_results_subset, simple_raster)

  # Convert raster to SpatialPointsDataFrame
  unraked_preds.sp <- rasterToPoints(results_subset, spatial=TRUE)
  projection <- proj4string(unraked_preds.sp)

  # reproject sp object
  unraked_preds.sp <- spTransform(unraked_preds.sp, CRS(projection))
  unraked_preds.sp@data <- data.frame(unraked_preds.sp@data, long=coordinates(unraked_preds.sp)[,1],lat=coordinates(unraked_preds.sp)[,2])
  unraked_preds.dt <- data.table(unraked_preds.sp@data)
  names(unraked_preds.dt)[names(unraked_preds.dt) == "lat"] = "latitude"
  names(unraked_preds.dt)[names(unraked_preds.dt) == "long"] = "longitude"

  color_list <- c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')
  results_limits=c(min(c(minValue((results_subset)), minValue((raked_results_subset)))),
                   max(c(maxValue((results_subset)), maxValue((raked_results_subset)))))

  admin0 <- admin0[admin0$ADM0_CODE %in% gaul, ]
  admin0.dt <- data.table(fortify(admin0))

  ## Reshape for the gg
  unraked_preds.dt = melt(unraked_preds.dt, id.vars = c("longitude", "latitude"), measure.vars = grep(indicator, names(unraked_preds.dt), value = T))
  unraked_preds.dt <- unraked_preds.dt[, year := as.character(variable)]
  unraked_preds.dt <- unraked_preds.dt[, year := as.numeric(gsub(paste0(indicator, '_mean_raster.'), '', year))]
  unraked_preds.dt <- unraked_preds.dt[year %in% c(1,6,11,16), ]
  unraked_preds.dt <- unraked_preds.dt[, year := year + 2000 - 1]

  results_gg <- ggplot(unraked_preds.dt,aes(longitude,latitude)) +
    geom_raster(aes(fill=value)) +
    coord_fixed() +
    theme_minimal() +
    geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='black', lwd=.1) +
    scale_fill_gradientn(colours=(color_list), limits=results_limits, na.value = "#000000") +
    guides(fill=guide_colorbar(title="Mean outcome", label=TRUE, ticks=FALSE)) +
    scale_x_continuous("", breaks=NULL) +
    scale_y_continuous("", breaks=NULL) +
    theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) +
    facet_wrap(~year, ncol = 4)

  # Convert raster to SpatialPointsDataFrame
  raked_preds.sp <- rasterToPoints(raked_results_subset, spatial=TRUE)
  projection <- proj4string(raked_preds.sp)

  # reproject sp object
  raked_preds.sp <- spTransform(raked_preds.sp, CRS(projection))
  raked_preds.sp@data <- data.frame(raked_preds.sp@data, long=coordinates(raked_preds.sp)[,1],lat=coordinates(raked_preds.sp)[,2])
  raked_preds.dt <- data.table(raked_preds.sp@data)
  names(raked_preds.dt)[names(raked_preds.dt) == "lat"] = "latitude"
  names(raked_preds.dt)[names(raked_preds.dt) == "long"] = "longitude"

  ## Reshape for the gg
  raked_preds.dt = melt(raked_preds.dt, id.vars = c("longitude", "latitude"), measure.vars = grep(indicator, names(raked_preds.dt), value = T))
  raked_preds.dt <- raked_preds.dt[, year := as.character(variable)]
  raked_preds.dt <- raked_preds.dt[, year := as.numeric(gsub(paste0(indicator, '_mean_raked_raster.'), '', year))]
  raked_preds.dt <- raked_preds.dt[year %in% c(1,6,11,16), ]
  raked_preds.dt <- raked_preds.dt[, year := year + 2000 - 1]

  raked_results_gg <- ggplot(raked_preds.dt,aes(longitude,latitude)) +
    geom_raster(aes(fill=value)) +
    coord_fixed() +
    theme_minimal() +
    geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='black', lwd=.1) +
    scale_fill_gradientn(colours=(color_list), limits=results_limits, na.value = "#000000") +
    guides(fill=guide_colorbar(title="Mean outcome", label=TRUE, ticks=FALSE)) +
    scale_x_continuous("", breaks=NULL) +
    scale_y_continuous("", breaks=NULL) +
    theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) +
    facet_wrap(~year, ncol = 4)

  all_plots <- list(results_gg, raked_results_gg)
  names(all_plots) <- c('unraked','raked')
  return(all_plots)

}

val_raking_factors <- function(rfs,
                               shapefile_version = 'current', 
                               gaul_list) {

  rfs <- rfs[name %in% gaul_list, ]
  loc_names <- get_location_code_mapping(shapefile_version = shapefile_version)
  setnames(rfs, "name", "GAUL_CODE")
  rfs <- merge(rfs, loc_names, by="GAUL_CODE")
  gg_rfs <- ggplot(data = rfs,
                   aes(x = rake_to_mean,
                       y = geo_mean,
                       color = year)) +
    geom_text(label = rfs$year) +
    geom_line(aes(x = geo_mean, y = geo_mean), col='black', size=1) +
    ylab("MBG Mean") +
    xlab("GBD Mean") +
    theme_minimal() +
    theme(legend.position="none")
  return(gg_rfs)

}

parse_region_name <- function(i, regions) {

  region_start <- grep('bin', regions[[1]]) + 1
  this_reg_name <- regions[[i]][region_start]
  ## Handle custom regions that may have an '_' in them.
  for(add in (region_start+1):length(regions[[i]])) {
    if(!grepl('RData', regions[[i]][add])) this_reg_name <- paste0(this_reg_name, '_', regions[[i]][add])
  }
  return(this_reg_name)

}

## updated and moved to misc_functions.R Sept 2018
## get_output_regions <- function(in_dir) {

##   regions <- list.files(in_dir, pattern = paste0(indicator, '_cell_draws_eb_'))
##   regions <- strsplit(regions, '_')
##   regions <- unlist(lapply(1:length(regions), parse_region_name, regions = regions))
##   regions <- unique(regions)
##   return(regions)

## }


## Pull draws by country and save template rasters by country in memory
pull_country_draws <- function(reg, periods, raked = "", in_dir,
                               pop_measure, start_year, end_year,
                               admin2_shapes, all_region_pops,
                               shapefile_version = 'current') {

  ## Load simple_raster and pop_raster into memory for this GAUL_CODE
    gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                               buffer = 0.4,
                                               subset_only = TRUE,
                                               shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]

  ## Load draws into memory for this GAUL_CODE.
    message(paste0('Loading draws for ', reg,'...'))
    load(paste0(in_dir, '/', indicator, raked, '_cell_draws_eb_bin0_',reg,'_0.RData'))
    draws <- as.data.table(cell_pred)
    rm(cell_pred)

  ## Get dt or draws for this region indexed by GAUL_CODE.
    cell_idx <- seegSDM:::notMissingIdx(simple_raster)
    coords <- xyFromCell(simple_raster, seegSDM:::notMissingIdx(simple_raster))
    template <- raster::extract(simple_raster,coords)
    draws_by_gaul <- draws[, GAUL_CODE := rep(template,periods)]

  ## Pull out draws/templates for each GAUL_CODE in this region.
    message('Pulling draws by country for ', reg, '...')
    make_country_list <- function(gaul, reg_simple_raster) {
      message(gaul)
      # Save draws
      country_draws <- draws_by_gaul[GAUL_CODE==gaul,]
      country_draws <- country_draws[, GAUL_CODE := NULL]
      # Save a template raster and population raster
      country_simple_raster <- reg_simple_raster
      country_simple_raster[country_simple_raster!=gaul] <- NA
      # gaul_list <- gaul
      # country_shape  <- subset_shape[subset_shape@data$GAUL_CODE %in% gaul_list, ]
      # raster_list    <- build_simple_raster_pop(country_shape)
      # simple_raster  <- raster_list[['simple_raster']]
      # pop_raster     <- raster_list[['pop_raster']]
      # Get population brick for all periods
      pop_raster_annual  <- all_region_pops
      pop_raster_annual  <- crop(pop_raster_annual, extent(country_simple_raster))
      pop_raster_annual  <- setExtent(pop_raster_annual, country_simple_raster)
      pop_raster_annual  <- mask(pop_raster_annual, country_simple_raster)
      # Get raster of admin2 codes
      admin_level <- 2
      shapes <- admin2_shapes
      cropped_shapes <- crop(shapes, extent(country_simple_raster), snap="out")
        ## Fix rasterize
        initial_raster <- rasterize_check_coverage(cropped_shapes, country_simple_raster, field = paste0('ADM', admin_level,'_CODE'))
        if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])!=0) {
          rasterized_shape <- merge(rasterize_check_coverage(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),], country_simple_raster, field = paste0('ADM', admin_level,'_CODE')), initial_raster)
        }
        if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])==0) {
          rasterized_shape <- initial_raster
        }
        masked_shapes <- mask(x=rasterized_shape, mask=country_simple_raster)
      # Add items to list
      return_list <- list(country_draws,
                          pop_raster_annual,
                          country_simple_raster,
                          masked_shapes)
      names(return_list) <- c(paste0('draws_', gaul), paste0('pops_', gaul), paste0('simple_', gaul), paste0('admin2_', gaul))
      return(return_list)
    }
    reg_gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
    reg_gaul_list <- reg_gaul_list[!(reg_gaul_list %in% 40762)]
    return_list <- lapply(reg_gaul_list, make_country_list, 
                          reg_simple_raster = simple_raster)
    names(return_list) <- paste0('list_', reg_gaul_list)

  ## Return list of draws, pops, templates
    return(return_list)

}


summarize_admin2 <- function(gaul, indicator, indicator_group, run_date, nperiod, master_list) {

  message(paste0('Working on ', gaul, '...'))
  #survey_year <- 2015
  #nperiod <- total_periods
  #this_period <- survey_year - 2000 + 1

  ## Load draws and simple_raster by country. Create pixel_index so we don't mess up sorting.
    cell_pred.dt <- as.data.table(master_list[[paste0('list_', gaul, '.draws_', gaul)]])
    cell_pred.dt[cell_pred.dt < 0] <- 0
    names <- names(cell_pred.dt)
    cols <- names(cell_pred.dt)
    period_index <- c()
    for(i in 1:nperiod) {
      period_index <- c(period_index, rep(paste0("period_", i), length(cell_pred.dt$V1)/nperiod))
    }
    pixel_id <- c(rep(1:(length(cell_pred.dt$V1)/nperiod),
                      nperiod))
    cell_pred.dt <- cbind(cell_pred.dt, period_index, pixel_id)

  ## Summarize and subset
    #cell_pred.dt <- cell_pred.dt[period_index == paste0('period_', this_period), ]

  ## Get pops
    country_pops <- master_list[[paste0('list_', gaul, '.pops_', gaul)]]
    country_pops <- crop(country_pops, extent(master_list[[paste0('list_', gaul, '.simple_', gaul)]]))
    country_pops <- setExtent(country_pops, master_list[[paste0('list_', gaul, '.simple_', gaul)]])
    country_pops <- mask(country_pops, master_list[[paste0('list_', gaul, '.simple_', gaul)]])

  ## Get admin2 codes
    country_admin2 <- master_list[[paste0('list_', gaul, '.admin2_', gaul)]]
    country_admin2 <- crop(country_admin2, extent(master_list[[paste0('list_', gaul, '.simple_', gaul)]]))
    country_admin2 <- setExtent(country_admin2, master_list[[paste0('list_', gaul, '.simple_', gaul)]])
    country_admin2 <- mask(country_admin2, master_list[[paste0('list_', gaul, '.simple_', gaul)]])

  ## Get ids for all cells
    cell_idx <- seegSDM:::notMissingIdx(master_list[[paste0('list_', gaul, '.simple_', gaul)]])

  ## Make full datatable of draws with admin2 and pop info
    geo.dt <- cell_pred.dt
    #admin2_codes <- rep(extract(country_admin2, cell_idx), nperiod)
    pull_period_pops <- function(period) {
      geo.subset <- geo.dt[period_index == paste0('period_', period), ]
      geo.subset <- geo.subset[, pops := raster::extract(country_pops[[period]], cell_idx)]
      geo.subset <- geo.subset[, admin2 := raster::extract(country_admin2, cell_idx)]
      return(geo.subset)
    }
    geo.dt <- rbindlist(lapply(1:nperiod, pull_period_pops))
    geo.dt <- geo.dt[is.na(pops), pops := 0] # weighted.mean doesn't like NA weights

  ## Make national and admin2 population-weighted means for each draw
    natl_mean_draws <- geo.dt[, lapply(.SD, weighted.mean, w=pops, na.rm=TRUE), by=c('period_index'), .SDcols=grep("^V", names(geo.dt)) ]
    admin2_mean_draws <- geo.dt[, lapply(.SD, weighted.mean, w=pops, na.rm=TRUE), by=c('period_index', 'admin2'), .SDcols=grep("^V", names(geo.dt)) ]
  ## Make national and admin2 summaries across draws
    natl_mean_draws <- natl_mean_draws[, lower := apply(.SD, 1, quantile, c(.025), na.rm=T), .SDcols=grep("^V", names(natl_mean_draws))]
    natl_mean_draws <- natl_mean_draws[, mean := apply(.SD, 1, mean), .SDcols=grep("^V", names(natl_mean_draws))]
    natl_mean_draws <- natl_mean_draws[, upper := apply(.SD, 1, quantile, c(.975), na.rm=T), .SDcols=grep("^V", names(natl_mean_draws))]
    admin2_mean_draws <- admin2_mean_draws[, lower := apply(.SD, 1, quantile, c(.025), na.rm=T), .SDcols=grep("^V", names(admin2_mean_draws))]
    admin2_mean_draws <- admin2_mean_draws[, mean := apply(.SD, 1, mean), .SDcols=grep("^V", names(admin2_mean_draws))]
    admin2_mean_draws <- admin2_mean_draws[, upper := apply(.SD, 1, quantile, c(.975), na.rm=T), .SDcols=grep("^V", names(admin2_mean_draws))]
  ## Make plots
    rfs <- fread(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_rf.csv'))
    natl_mean_draws <- natl_mean_draws[, year := as.numeric(gsub('period_','',period_index)) + 2000 - 1]
    natl_mean_draws <- natl_mean_draws[, c('mean','upper','lower','year'), with=FALSE]
    admin2_mean_draws <- admin2_mean_draws[, year := as.numeric(gsub('period_','',period_index)) + 2000 - 1]
    admin2_mean_draws <- admin2_mean_draws[, c('mean','upper','lower','year','admin2'), with=FALSE]
    high_admin2 <- admin2_mean_draws[order(-mean)]
    high_admin2 <- high_admin2[year == 2015, .SD[1:1]]
    high_admin2 <- high_admin2[, admin2]
    low_admin2 <- admin2_mean_draws[order(mean)]
    low_admin2 <- low_admin2[year == 2015, .SD[1:1]]
    low_admin2 <- low_admin2[, admin2]
    admin2_mean_draws <- admin2_mean_draws[admin2 %in% high_admin2, Category := 'MBG highest admin2']
    admin2_mean_draws <- admin2_mean_draws[admin2 %in% low_admin2, Category := 'MBG lowest admin2']
    admin2_mean_draws <- admin2_mean_draws[!is.na(Category), ]
    natl_mean_draws <- natl_mean_draws[, Category := 'MBG admin0']
    setnames(rfs, 'rake_to_mean', 'mean')
    rfs <- rfs[name==gaul, c('mean','year'), with = FALSE]
    rfs <- rfs[, Category := 'GBD admin0']
    all_data <- rbind(natl_mean_draws, admin2_mean_draws, rfs, fill=TRUE)

    gaul_gg <- ggplot() +
      ## Admin2s
      geom_line(data=all_data,
                aes(x = year,
                    y = mean,
                    group = Category,
                    color = Category,
                    size = Category)) +
      geom_ribbon(data=all_data,
                  aes(x = year,
                      ymin=lower,
                      ymax=upper,
                      group = Category,
                      fill = Category),
                  alpha=0.2) +
      scale_colour_manual(
        values = c('black', '#2ca25f', '#2b8cbe', '#e34a33')
      ) +
      scale_fill_manual(
        values = c('black', '#2ca25f', '#2b8cbe', '#e34a33')
      ) +
      scale_size_manual(
        values = c(2,2,1,1)
      ) +
      theme_minimal() +
      ylab('Mean\noutcome') +
      xlab('Year')

    return(gaul_gg)

}

plot_varios_and_violins <- function(indicator,
                                    indicator_group,
                                    run_date,
                                    this_gaul,
                                    shapefile_version = 'current') {

  in_dir  <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date)
  default_rf_path <- paste0(in_dir, '/', indicator, '_rf.csv')
  all_rfs <- fread(default_rf_path)

  ## Define path to .tif of results and raked results, load in raster bricks.
  default_raked_results_path <- brick(paste0(in_dir, '/', indicator, '_mean_raked_raster.tif'))
  results_raked <- brick(default_raked_results_path)
  default_results_path <- brick(paste0(in_dir, '/', indicator, '_mean_raster.tif'))
  results <- brick(default_results_path)
  total_periods <- length(names(results))
  regions <- get_output_regions(in_dir)
  for(reg in regions) {
    if(this_gaul %in% get_adm0_codes(reg, shapefile_version = shapefile_version)) this_region <- reg
  }
  input_data <- list.files(in_dir, pattern = paste0("input_data[a-zA-Z0-9_]*", this_region), full.names = T) %>% fread

  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  gauls <- this_gaul
  this_country <- gaul_to_loc_id[GAUL_CODE %in% gauls, ihme_lc_id]

  if(indicator != 'edu_mean') input_data <- input_data[, rate := get(indicator) / N]
  if(indicator == 'edu_mean') input_data <- input_data[, rate := get(indicator)]
  country_input <- input_data[country %in% this_country, ]
  country_input <- subset(country_input, year >= 1998)
  names(country_input)[names(country_input) == "year"] = "original_year"
  country_input <- country_input[original_year >= 1998 & original_year < 2003, year := 2000]
  country_input <- country_input[original_year >= 2003 & original_year < 2008, year := 2005]
  country_input <- country_input[original_year >= 2008 & original_year < 2013, year := 2010]
  country_input <- country_input[original_year >= 2013 & original_year < 2018, year := 2015]
  country_input_dt <- country_input
  coordinates(country_input) = ~longitude+latitude

  ## Calculate variogram for data
  v_all = variogram(log(rate+0.001)~1, country_input)
  v_all$year_category <- "All years"
  ## Do five-year period variograms
  vario_list <- list()
  for(this_year in unique(country_input_dt[, year])) {
    year_subset <- country_input_dt[year == this_year, ]
    coordinates(year_subset) = ~longitude+latitude
    v = variogram(log(rate+0.001)~1, year_subset)
    v$year_category <- as.character(paste0(this_year, ' window'))
    #message(paste0('vario_plot_', this_year))
    if('data.frame' %in% class(v) == FALSE) v <- NULL
    assign(paste0('vario_plot_', this_year), v)
    vario_list[[as.character(this_year)]] <- v
  }
  #vario_years <- ls()[(grep("vario_plot_", ls()))]
  #vario_years <- rbindlist(lapply(vario_years, get))
  vario_years <- rbindlist(vario_list)
  vario_years <- rbind(v_all, vario_years)
  ## Make variogram
  gg_country_vario <- ggplot(data=vario_years, aes(x = dist, y = gamma, col = year_category)) +
    geom_point() +
    geom_line() +
    ylab('Semivariance') +
    xlab('Distance') +
    theme_minimal() +
    guides(col=guide_legend(title="Year category"))

  ## Make PACF
  country_input_dt_acf <- country_input_dt[, list(new_rate = weighted.mean(x=rate,w=weight*N)), by = c('original_year')]
  country_input_dt_acf <- country_input_dt_acf[order(original_year)]
  #country_pacf <- pacf(country_input_dt_acf[, new_rate], main = "PACF")
  country_pacf <- 'placeholder' ## Super annoying to grab plot object right now

  ## Calculate variogram for predictions
  ## Extract preds for each year in country data
  default_raked_results_path <- paste0(in_dir, '/', indicator, '_mean_raked_raster.tif')
  results_raked <- brick(default_raked_results_path)
  default_results_path <- paste0(in_dir, '/', indicator, '_mean_raster.tif')
  results <- brick(default_results_path)
  extract_year_preds <- function(country_year) {
    period <- country_year - 2000 + 1
    input_data_year <- country_input_dt[year == country_year, ]
    preds_year <- results[[period]]
    preds_at_points <- raster::extract(preds_year, input_data_year[, c('longitude','latitude'), with = F])
    input_data_year <- input_data_year[, pred := preds_at_points]
    return(input_data_year)
  }
  country_input_wpreds <- rbindlist(lapply(unique(country_input_dt[,year]), extract_year_preds))
  country_input_wpreds <- country_input_wpreds[!is.na(pred),]
  country_input_wpreds <- country_input_wpreds[, data := rate]
  violin_vars <- c('data', 'pred')
  country_input_wpreds <- country_input_wpreds[, c('latitude','longitude','weight', violin_vars ,'year'), with=F]
  country_input_wpreds <- melt(country_input_wpreds, id.vars = c('longitude','latitude','weight','year'), measure.vars = violin_vars)
  country_input_wpreds <- country_input_wpreds[, country := this_country]
  if(file.exists(paste0(results_dir,'/', indicator, '_rf.csv'))) {
    rfs <- fread(paste0(results_dir,'/', indicator, '_rf.csv'))
    country_input_wpreds <- country_input_wpreds[, name := as.integer(this_gaul)]
    country_input_wpreds <- merge(country_input_wpreds, rfs, by=c('name','year'))
  }
  violin_countries <- this_country
  loc_names <- get_location_code_mapping(shapefile_version = shapefile_version)
  convert_to_iso3 <- function(x) {
    if(nchar(x) > 3) x <- loc_names[loc_name==x, ihme_lc_id]
    return(x)
  }
  violin <- ggplot(data=country_input_wpreds) +
    geom_violin(aes(x = variable, y = value, fill = country, weight = weight))
  if(file.exists(paste0(results_dir,'/', indicator, '_rf.csv'))) {
    violin <- violin + geom_hline(aes(yintercept=rake_to_mean))
  }
  violin <- violin +
    theme(axis.text.x = element_text(face="bold", size=10, angle=90)) +
    theme(legend.position="none") +
    xlab("") +
    ylab("Mean outcome") +
    facet_wrap(~country + year)

  all_plots <- list(gg_country_vario, country_pacf, violin)

  ## Return ggs for violin, vario, and pacf model object to plot() later.
  return(all_plots)

}

raking_factors_map <- function(indicator,
                               indicator_group,
                               run_date,
                               gaul,
                               shapefile_version, 
                               year_plot_list) {

  require(sp)
  library(rgeos)
  library(maptools)
  library(ggplot2)

  # extract admin0
  master_shape <- readOGR(dsn='/snfs1/temp/learl/data_coverage', layer='africa_ad0')
  names(master_shape)[names(master_shape)=="ADM0_CODE"] <- "GAUL_CODE"
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  master_shape <- merge(master_shape, gaul_to_loc_id, by="GAUL_CODE")
  subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% gaul, ]
  admin0.dt <- data.table(fortify(subset_shape))

  rf <- fread(paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/', indicator, '_rf.csv'))
  names(rf)[names(rf)=="name"] <- "GAUL_CODE"
  rf$raking_factor[is.infinite(rf$raking_factor)] <- -1
  rf <- rf[GAUL_CODE %in% gaul, ]

  # Loop over periods
  make_period_matrix <- function(period) {

    rf_shape <- merge(subset_shape,rf[rf$year==period,], by="GAUL_CODE")
    rf_shape@data$id = rownames(rf_shape@data)
    rf.pts <- fortify(rf_shape, region="id")
    rf.df = join(rf.pts, rf_shape@data, by="id")
    rf.df$raking_factor[is.na(rf.df$raking_factor)] <- 0
    rf.df$raking_factor <- as.numeric(rf.df$raking_factor)
    rf.df$year <- period
    return(rf.df)

  }
  all_rf_df <- rbindlist(lapply(year_plot_list, make_period_matrix))

  # Plot
  color_list <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
  redblue <- c('#313695','#ffffff','#a50026')
  admin2.gg <- ggplot(all_rf_df,aes(long,lat,group=group)) +
    geom_polygon(aes(fill=raking_factor)) +
    geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='black', size=.5) +
    scale_fill_gradientn(colours=redblue, values=c(0.5,1,2), limits=c(0.5,2), na.value = "#000000", rescaler = function(x, ...) x, oob = identity) +
    guides(fill=guide_colorbar(label=TRUE, ticks=FALSE)) +
    scale_x_continuous("", breaks=NULL) +
    scale_y_continuous("", breaks=NULL) +
    guides(fill=guide_legend(title="Raking\nfactor")) +
    coord_equal() +
    theme_minimal() +
    facet_wrap(~year, ncol=4)
  return(admin2.gg)

}


make_beta_plot <- function(indicator_group,
                           indicator,
                           run_date,
                           region) {

  # fit_dir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/fit_stats')
  # load(paste0(fit_dir,'/strata.RData'))
  # Regions=strata
  process_betas <- function(reg) {

    load(paste0('/share/geospatial/mbg/',indicator_group, '/', indicator, '/output/', run_date,'/', indicator, '_model_eb_bin0_', reg, '_0.RData'))
    model_data <- as.data.table(res_fit$summary.fixed)
    names(model_data)[names(model_data)=='0.025quant'] <- 'lower'
    names(model_data)[names(model_data)=='0.975quant'] <- 'upper'
    model_data[, cov := row.names(res_fit$summary.fixed)]
    model_data <- model_data[cov!='int']
    model_data[, region := reg]
    return(model_data)

  }
  model_data <- lapply(region, process_betas)
  model_data <- do.call(rbind.fill, model_data)
  model_data <- as.data.table(model_data)
  model_data <- model_data[!grepl('gaul_', cov), ] # Remove country fixed effects from this plot
  selected_gg <- ggplot(model_data, aes(x=cov, y=mean)) +
    geom_point() +
    geom_hline(yintercept=0, color='red') +
    geom_errorbar(aes(ymax = upper, ymin = lower), width=0.25) +
    theme(axis.text.x = element_text(angle = 75, hjust = 1, size=15)) +
    ggtitle('Selected model') +
    theme(strip.text.x = element_text(size = 15))

  return(selected_gg)

}


## ################################
## is_oos_preds
##
## this function takes the data that went into our MBG framework and
## returns an array containing: location, observed value, in sample
## predictions, and out of sample predictions (when applicable), country, year, N
##
## INPUT:
##
## OUTPUT:
##
## #################################

is_oos_preds <- function(rd = run_date,
                         all.data = df,
                         cell_draws_filename = '%s_cell_draws_eb_bin%i_%s_%i.RData', ## in sprintf notation
                         holdouts = 5, ## number of holdouts. if zero only does in sample
                         reg,
                         years = 2000:2015,
                         indic = indicator,
                         indic_group = indicator_group,
                         holdoutlist = NULL, ## if null, only does in sample
                         fun = "mean", ## function to sapply over rows
                         shapefile_version = 'current', 
                         ... ## additional arguments can be passed to `fun`
){

  ## Make data.table
  df <- as.data.table(df)

  ## place to look for things
  output.dir <- sprintf("/share/geospatial/mbg/%s/%s/output/%s/",
                        indic_group, indic, rd)

  ## Load data
  datdir <- sprintf('/share/geospatial/mbg/%s/%s/output/%s/',indic_group,indic,rd)

  ## holdout data
  if(!is.null(holdoutlist)){
    d <- data.frame(holdoutlist[[sprintf('region__%s', reg)]])
  }else{
    d <- df
  }

  ## load the simple raster for this region
  message("Loading simple raster...")
  if(file.exists(paste0('/share/geospatial/shapefiles/simple_raster',reg,'.RData'))) {
    load(paste0('/share/geospatial/shapefiles/simple_raster',reg,'.RData'))
  } else {
    simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(regions,
                                                                          shapefile_version = shapefile_version),
                                               buffer = 0.4,
                                               subset_only = TRUE,
                                               shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
  }

  ## setup the data in the data order to return
  if(holdouts == 0){
    return.df <- d
  }else{
    return.df <- d[d$fold == 1, ]
    for(i in 2:holdouts){
      return.df <- rbind(return.df, d[d$fold == i, ])
    }
  }
  return.df <- return.df[, c('longitude', 'latitude', 'year', 'country', 'source',
                              indic, 'N','weight'), with=FALSE]

  return.df$OOS <- NA
  return.df$IS  <- NA

  ## ####################
  ## get out of sample ##
  ## ####################
  OOS <- NULL

  ## loop through the holdouts
  if(holdouts!=0) {
    message("Starting out-of-sample... ")
    for(hh in 1:holdouts){

      ## load the associated preds
      message(paste0("Loading cell preds for holdout ", hh, "..."))
      load(sprintf(paste0('/share/geospatial/mbg/%s/%s/output/%s/', cell_draws_filename),
                   indic_group, indic, rd, indic, 0, reg, hh))

      ## obtain desired summary of cell pred
      message("Summarizing cell preds...")
      ## if mean, will use rowMeans() for speed; otherwise, apply
      if (fun == "mean") {
        ## average across the draws
        summary.cell.pred <- rowMeans(cell_pred)
      } else {
        summary.cell.pred <- apply(cell_pred, 1, get(fun), ...)
      }

      ## turn into a raster
      temp.rast <- insertRaster(simple_raster, matrix(summary.cell.pred, ncol = length(years)))

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
  message("Starting in-sample...")

  ## regardless of whether holdouts==0 or not, we can do in sample extraction
  hh <- 0
  d.is <- d
  d.is$fold <- 0

  ## load the associated preds
  message("Loading cell preds...")
  load(sprintf(paste0('/share/geospatial/mbg/%s/%s/output/%s/', cell_draws_filename),
               indic_group, indic, rd, indic, 0, reg, hh))

  ## obtain desired summary of cell pred
  message("Summarizing cell preds...")
  ## if mean, will use rowMeans() for speed; otherwise, apply
  if (fun == "mean") {
    ## average across the draws
    summary.cell.pred <- rowMeans(cell_pred)
  } else {
    summary.cell.pred <- apply(cell_pred, 1, get(fun), ...)
  }

  ## turn into a raster
  temp.rast <- insertRaster(simple_raster, matrix(summary.cell.pred, ncol = length(years)))

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


load_region_input_data <- function(indicator_group, indicator, run_date, reg) {

  model_path <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date)
  df <- list.files(model_path, pattern = paste0("input_data[a-zA-Z0-9_]*", reg), full.names = T) %>% fread
  return(df)

}

plot_quilt <- function(gaul_list, df, sample, subset_shape, shapefile_version = 'current') { ## sample == "IS" or "OOS"

  df <- df[!is.na(get(sample)), ]
  subset_shape <- subset_shape[subset_shape$GAUL_CODE %in% gaul_list, ]
  df <- as.data.table(df)
  loc_names <- get_location_code_mapping(shapefile_version = shapefile_version)
  setnames(df, "country", "ihme_lc_id")
  df <- merge(df, loc_names, by="ihme_lc_id")
  df <- df[GAUL_CODE %in% gaul_list, ]
  if(length(df[, GAUL_CODE]) != 0) {
    make_quilt_matrix <- function(this_year) {

      ## Run quilt function from fields library, suppress inline plot
      df_year <- df[year == this_year, ]

      if(nrow(unique(df_year[,c("longitude", "latitude")])) > 1) {
        # checks to see if there are at least 2 points, otherwise quilt.plot fails
        quilt_x <- df_year[, longitude]
        quilt_y <- df_year[, latitude]
        quilt_z <- df_year[, get(sample)]
        pdf("NULL")
        quilt <- fields::quilt.plot(quilt_x, quilt_y, quilt_z, main="Absolute error")
        dev.off()

        ## Pull spatial info from quilt object
        long_matrix <- melt(quilt$z) # Make matrix of values from quilt.plot long df
        z_raster <- rasterFromXYZ(long_matrix) # Convert first two columns as lon-lat and third as value
        proj4string(z_raster) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
        z_raster <- setExtent(z_raster, extent(min(quilt$x), max(quilt$x), min(quilt$y), max(quilt$y)))

        # Convert raster to SpatialPointsDataFrame
        z_raster.sp <- rasterToPoints(z_raster, spatial=TRUE)
        projection <- proj4string(z_raster.sp)

        # reproject sp object
        z_raster.sp <- spTransform(z_raster.sp, CRS(projection))
        z_raster.sp@data <- data.frame(z_raster.sp@data, long=coordinates(z_raster.sp)[,1],lat=coordinates(z_raster.sp)[,2])
        z_raster.dt <- data.table(z_raster.sp@data)
        names(z_raster.dt)[names(z_raster.dt) == "lat"] = "latitude"
        names(z_raster.dt)[names(z_raster.dt) == "long"] = "longitude"
        z_raster.dt <- z_raster.dt[, year := this_year]
      } else {
        z_raster.dt <- data.table(
          value = numeric(),
          longitude = numeric(),
          latitude = numeric(),
          year = numeric())
      }
      return(z_raster.dt)
    }
    quilt_list <- rbindlist(lapply(sort(unique(df[, year])), make_quilt_matrix))
    color_list <- c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')
    this_shape.dt <- data.table(fortify(subset_shape))
    redblue <- c('#313695','#ffffff','#a50026')
    # plot gg
    quilt_gg <- ggplot(quilt_list,aes(longitude,latitude)) +
      geom_raster(aes(fill=value)) +
      coord_fixed() +
      theme_minimal() +
      geom_path(data=this_shape.dt, aes(x=long, y=lat, group=group), color='black', lwd=.1) +
      #scale_fill_gradientn(colours=(color_list), limits=c(min(df[, get(sample)]), max(df[, get(sample)])), na.value = "#000000") +
      scale_fill_gradientn(colours=redblue, values=c(min(df[, get(sample)]),0,max(df[, get(sample)])), limits=c(min(df[, get(sample)]),max(df[, get(sample)])), na.value = "#000000", rescaler = function(x, ...) x, oob = identity) +
      guides(fill=guide_colorbar(title="Absolute\nerror", label=TRUE, ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) +
      facet_wrap(~year)
    return(quilt_gg)
  }
  if(length(df[, GAUL_CODE]) == 0) {
    return(NULL)
  }
}


compare_admin_estimates <- function(nid_list, results_raster,
                                    compare_source_title, raked_title,
                                    outcome_title, compare_data,
                                    compare_spdf, master_list,
                                    shapefile_version, 
                                    by_source = FALSE) {

  ## Define function to extract preds/pops over comparison polygons, and calculate pop-weighted mean outcome over polygons.
  compile_polygons <- function(this_nid, estimates_raster) {
    message(this_nid)
    nid_results <- compare_data[nid == this_nid, ]

    gaul <- gaul_convert(nid_results[, iso3][1], shapefile_version = shapefile_version)
    message(paste0('gaul', gaul))
    country_pops <- master_list[[paste0('list_', gaul, '.pops_', gaul)]]
    country_pops <- crop(country_pops, extent(master_list[[paste0('list_', gaul, '.simple_', gaul)]]))
    country_pops <- setExtent(country_pops, master_list[[paste0('list_', gaul, '.simple_', gaul)]])
    country_pops <- mask(country_pops, master_list[[paste0('list_', gaul, '.simple_', gaul)]])
    message('crop')
    country_estimates <- crop(estimates_raster, extent(master_list[[paste0('list_', gaul, '.simple_', gaul)]]))
    country_estimates <- setExtent(country_estimates, master_list[[paste0('list_', gaul, '.simple_', gaul)]])
    country_estimates <- mask(country_estimates, master_list[[paste0('list_', gaul, '.simple_', gaul)]])

    all_data <- merge(compare_spdf, nid_results, by=c('location_code','shapefile'))
    all_data <- all_data[!is.na(all_data@data$outcome),]
    all_data$geo_mean <- 0
    for(shape_x in unique(all_data$location_code)) {
      message(shape_x)
      test_poly <- all_data[all_data$location_code == shape_x, ]
      period <- test_poly$year[1] - 2000 + 1
      preds <- extract(country_estimates[[period]], test_poly)
      pops <- extract(country_pops[[period]], test_poly)
      all_data$geo_mean[all_data$location_code == shape_x] <- weighted.mean(preds[[1]], pops[[1]], na.rm=T)
    }
    this_data <- as.data.table(all_data)
    return(this_data)
  }

  ## Pull unraked estimates
  all_data <- rbindlist(lapply(unique(compare_data[year > 2000 & nid %in% nid_list, nid]), compile_polygons,
                               estimates_raster = results_raster))

  ## Calculate super simple bias estimate
  ## Fit simple lm
  all_data <- all_data[!is.na(outcome) & !is.na(geo_mean), ]
  bias_model <- lm(outcome ~ geo_mean, data = all_data)
  ## Intercept
  intercept <- summary(bias_model)$coefficients[1, 1]
  intercept_p <- summary(bias_model)$coefficients[1, 4]
  if(intercept_p <= 0.0005) intercept_p <- 0.0005
  ## Slope
  slope <- summary(bias_model)$coefficients[2, 1]
  slope_p <- summary(bias_model)$coefficients[2, 4]
  if(slope_p <= 0.0005) slope_p <- 0.0005
  summary(bias_model)
  ## Predict fitted values for line
  all_data$lm_pred <- predict(bias_model)

  ## Make title for plot depending on if comparing a single source (grab year/iso3) or many sources
  if(length(unique(all_data[, nid]))==1) add_title <- paste0('\n', all_data[, iso3][1], ' ', all_data[, year][1], ' ', compare_source_title, '\n')
  if(length(unique(all_data[, nid]))!=1) add_title <- '\n'
  main_gg_title <- paste0(compare_source_title, ' vs. aggregated ', raked_title, ' MBG estimates', add_title, 'Intercept: ', round(intercept, 2), ' (p <= ', round(intercept_p, 4),')\nSlope: ', round(slope, 2), ' (p <= ', round(slope_p, 4),')')

  ## Plot gg for unraked
  comparison_gg <- ggplot() +
    geom_line(data = all_data,
              aes(x = lm_pred,
                  y = geo_mean)) +
    geom_point(data=all_data,
               aes(x=outcome,
                   y=geo_mean)) +
    geom_abline(slope=1, intercept=0, colour = "red") +
    xlab(paste0(outcome_title, ', ', compare_source_title)) +
    ylab(paste0(outcome_title, ', MBG estimates')) +
    ggtitle(main_gg_title) +
    theme_minimal()

  return(comparison_gg)

}

## Functions for comparison to SDG target

get_p_target <- function(draws, target_type, target) {

  ## Calculates probability of being at/over a given target

  # draws = row of n preds over which to apply the function
  # target_type = ">",">=","<","<="
  # target = SDG, etc - e.g. "0.8" for 80% vaccine coverage

  output <- sum(do.call(get(target_type), list(draws, target)))
  output <- output / length(draws)
  return(output)

}

plot_target <- function(gaul_list, df, sample, subset_shape, target, target_type, region=NULL, shapfile_version = 'current') {
  ## sample == "IS" or "OOS"
  ## region: NULL if a single GAUL, otherwise passes region name for plot
  # based on plot_quilt() by Nick G

  package_lib <- paste0('/home/j/temp/geospatial/packages')
  .libPaths(package_lib)
  cut2 <- Hmisc::cut2 # custom load one function from HMisc (avoid name conflicts)
  loc_names <- get_location_code_mapping(shapefile_version = shapefile_version)

  # get location name - either for single gaul or for region
  if (!is.null(region)) {
    loc_name <- paste0(": ", region)
  } else if (length(gaul_list) == 1) {
    loc_name <- paste0(": ", loc_names[GAUL_CODE == gaul_list, ihme_lc_id])
  } else {
    loc_name <- ""
  }

  message(paste0("\nRunning plot_target for ", substr(loc_name, 3, nchar(loc_name))))

  subset_shape <- subset_shape[subset_shape$GAUL_CODE %in% gaul_list, ]
  df <- as.data.table(df)
  setnames(df, "country", "ihme_lc_id")
  df <- merge(df, loc_names, by="ihme_lc_id")
  df <- df[GAUL_CODE %in% gaul_list, ]

  # Check for NAs & drop
  if (nrow(df[is.na(get(sample)),]) > 0) {
    warning(paste0("You have ", nrow(df[is.na(get(sample)),]), " rows in your data where ", sample, " is NA.",
                   " \nThis may reflect edge points from other countries, but check if large number. Dropping..."))
    df <- df[!is.na(get(sample)),]
  }

  # Create year bins
  df[year >= 1998 & year < 2003, bin_year := 2000]
  df[year >= 2003 & year < 2008, bin_year := 2005]
  df[year >= 2008 & year < 2013, bin_year := 2010]
  df[year >= 2013 & year < 2018, bin_year := 2015]

  # Determine whether the true data hit the target
  df[, outcome := get(indicator) / N]
  df[, hit_target := as.numeric(do.call(get(target_type), list(outcome, target)))]

  # Create an exposure variable (to encapsulate N and weight)
  df[, exposure := N*weight]

  # create graphs
  if(length(df[, GAUL_CODE]) != 0) {

     years <- sort(unique(df$bin_year))

    # by year
    make_year_table <- function(this_year) {
      df_year <- df[bin_year == this_year,]

      # Create bins by predicted probability for that year
      df_year$bins <- cut2(df_year[, get(sample)], cuts = seq(0,1,0.1))
      levels(df_year$bins) <- c("[0.0, 0.1)",
                                "[0.1, 0.2)",
                                "[0.2, 0.3)",
                                "[0.3, 0.4)",
                                "[0.4, 0.5)",
                                "[0.5, 0.6)",
                                "[0.6, 0.7)",
                                "[0.7, 0.8)",
                                "[0.8, 0.9)",
                                "[0.9, 1.0]")

      # calculate summary measures (weighted means) for each bin
      df_year_summ <- df_year[, list(wmean_pred = weighted.mean(get(sample), w = exposure),
                                     wmean_data = weighted.mean(hit_target, w = exposure),
                                     N_total = sum(N * weight)),
                                     by = bins]

     # now simulate using our predicted probability and the number of observations
     # try to capture uncertainty in sample sizes

     message(paste0(">> Simulating for year ", this_year, "..."))
     n_draws <- 500

     # simulate draws using N for each cluster and predicted probability
     # each cell represents the observed probability for a simulation
     # given N for the cluster and assuming that the predicted value is true

     draw_matrix <- mapply(function(draws, N, prob) rbinom(draws, N, prob) / N,
                             draws = n_draws,
                             N = df_year$N,
                             prob = df_year[,get(sample)]) %>% t

     # Function to obtain a weighted group mean for each draw
      group_wt_mean <- function(i) {
        probs <- draw_matrix[, i]
        df_wt <- subset(df_year, select = c("bins", "exposure"))
        df_wt[, prob := probs]
        df_summ <- df_wt[, list(wmean_p_hat = weighted.mean(prob, w = exposure)),
                                by = bins]
        setnames(df_summ, "wmean_p_hat", paste0("wmean_p_hat_", i))
      }

      # Create summary table (by draws) of weighted means by bin
      df_draw_summ <- lapply(1:ncol(draw_matrix), group_wt_mean) %>%
                      Reduce(merge, .)

      draw_cols <- names(df_draw_summ)[grepl("wmean_p_hat", names(df_draw_summ))]

      # generate 5th, 50th, and 95th %iles
      df_draw_summ <- cbind(subset(df_draw_summ, select = "bins"),
                            apply(subset(df_draw_summ, select = draw_cols),
                                  1,
                                  quantile,
                                  probs = c(0.05, 0.5, 0.95)) %>%
                            t %>%
                            as.data.table)

      df_summ <- merge(df_year_summ, df_draw_summ, by = "bins")
      df_summ[, year := this_year]

    }

    df_summ <- lapply(years, make_year_table) %>%
               rbindlist(.)

    setnames(df_summ, c("5%", "50%", "95%"), c("p5", "p50", "p95"))


    message(">> Building plots...")

    # make calibration plot
    plot_list <- lapply(c(2000, 2005, 2010, 2015), function(this_year) {
      df_plot <- df_summ[year == this_year, ]
      df_plot[, my_fill := "Predicted"]

      if(nrow(df_plot) > 0) {
        p_scatter <- ggplot(df_plot, aes(x = bins)) +
          geom_boxplot(stat = "identity",
                     aes(lower = p5, ymin = p5, middle = wmean_pred, upper = p95, ymax = p95, fill = my_fill),
                     colour = "darkgray") +
          geom_point(aes(y = wmean_data, shape = 'Observed'), colour = "black") +
          scale_x_discrete(drop = F) +
          theme_classic() +
          labs(x = "Prediction bins", y = "P(meets target)") +
          theme(axis.text.x=element_text(angle=45, hjust=1)) +
          scale_shape_manual('', values = c("Observed" = 19)) +
          scale_fill_manual('', values = c("Predicted" = "lightgray")) +
          ylim(0, 1)
      } else {
        p_scatter <- ggplot(df_plot) + geom_blank() + theme_classic()
      }
      return(p_scatter)
      })

    # make histograms
    hist_list <- lapply(c(2000, 2005, 2010, 2015),
      function(this_year){
        if(nrow(df_summ[year == this_year,]) > 0) {
          p_hist <- ggplot(df_summ[year == this_year,], aes(x = bins)) +
                 geom_bar(stat = "identity", aes(y = N_total)) +
                 scale_x_discrete(drop = F) +
                 theme_classic() +
                 labs(y = "N", x = "") +
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()) +
                 scale_shape_discrete(name = "Observed")

        } else {
          p_hist <- ggplot(df_summ[year == this_year,]) + geom_blank() + theme_classic()
        }
      })

    # Pull a legend
    g_legend<-function(a.gplot){
      pdf(NULL) # Workaround for bug in ggplot_gtable causing empty Rplots.pdf to be created
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      if (length(leg) > 0)  {
        legend <- tmp$grobs[[leg]]
      } else {
        legend <- NULL
      }
      graphics.off()
      return(legend)
    }

    my_legends <- lapply(plot_list, g_legend) # get the first one that works
    my_legends <- my_legends[vapply(my_legends, Negate(is.null), NA)]
    my_legend <- my_legends[[1]]

    plot_list <- lapply(plot_list, function(a_plot) {
      a_plot <- a_plot + theme(legend.position = "none")
      })

    lay <- rbind(c(5,5,6,6),
                 c(1,1,2,2),
                 c(1,1,2,2),
                 c(NA,NA,NA,NA),
                 c(7,7,8,8),
                 c(3,3,4,4),
                 c(3,3,4,4))

    # generate a plot for each year
    bin_years <- c(2000, 2005, 2010, 2015)
    year_plots <- lapply(1:length(bin_years),
      function(i) {
        year_plot <- arrangeGrob(grobs = list(plot_list[[i]], hist_list[[i]]),
                                 layout_matrix = rbind(2,1),
                                 heights = c(0.2, 0.8),
                                 top = textGrob(label = as.character(bin_years[i]),
                                                gp = gpar(fontsize = 24)))
        return(year_plot)
        })

    # generate 4-up of plots (both plot & hist)
    lay <- rbind(c(NA, NA),
                 c(1,2),
                 c(NA,NA),
                 c(3,4))

    all_plots <- arrangeGrob(grobs = year_plots, layout_matrix = lay, heights = c(0.1, 1, 0.2, 1))

    # add the final legend
    final_graph <- grid.arrange(all_plots, my_legend,
                                layout_matrix = rbind(c(1,NA,2)),
                                widths = c(7,0.2,0.8),
                                top = textGrob(label = paste0("Calibration", loc_name,
                                                              " (target ", as.character(target_type), " ",
                                                              as.character(target), ")"),
                                               gp = gpar(fontsize = 24)))
    return(final_graph)
   }
 }


check_for_holdouts <- function(dir) {
  cell_draw_files <- list.files(dir, pattern = "_cell_draws_")
  n_holdouts <- gsub('.*_([0-9]+).RData', '\\1', cell_draw_files) %>%
                  as.numeric %>%
                  max
}

## Pull draws by country and save template rasters by country in memory
pull_country_draws_all_admins <- function(reg, periods, raked = "",
                                          in_dir, pop_measure,
                                          start_year, end_year,
                                          admin2_shapes,
                                          admin1_shapes,
                                          all_region_pops,
                                          subtract = FALSE,
                                          subtract_cell_pred = '',
                                          shapefile_version = 'current') {

  ## Load simple_raster and pop_raster into memory for this GAUL_CODE
  gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list,
                                             buffer = 0.4,
                                             subset_only = TRUE,
                                             shapefile_version = shapefile_version)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]

  ## Load draws into memory for this GAUL_CODE.
  message(paste0('Loading draws for ', reg,'...'))
  load(paste0(in_dir, '/', indicator, raked, '_cell_draws_eb_bin0_',reg,'_0.RData'))
  if(raked=='_raked') {
    cell_pred <- raked_cell_pred
    rm(raked_cell_pred)
  }
  if(subtract == TRUE) {
    og_cell_pred <- cell_pred
    rm(cell_pred)
    load(subtract_cell_pred)
      if(raked=='_raked') {
        cell_pred <- raked_cell_pred
        rm(raked_cell_pred)
      }
    cell_pred <- og_cell_pred - cell_pred
  }
  draws <- as.data.table(cell_pred)
  rm(cell_pred)

  ## Get dt or draws for this region indexed by GAUL_CODE.
  cell_idx <- seegSDM:::notMissingIdx(simple_raster)
  coords <- xyFromCell(simple_raster, seegSDM:::notMissingIdx(simple_raster))
  template <- raster::extract(simple_raster,coords)
  draws_by_gaul <- draws[, GAUL_CODE := rep(template,periods)]

  ## Pull out draws/templates for each GAUL_CODE in this region.
  message('Pulling draws by country for ', reg, '...')
  make_country_list <- function(gaul) {
    message(gaul)
    # Save draws
    country_draws <- draws_by_gaul[GAUL_CODE==gaul,]
    country_draws <- country_draws[, GAUL_CODE := NULL]
    # Save a template raster and population raster
    gaul_list <- gaul
    country_shape  <- subset_shape[subset_shape@data$GAUL_CODE %in% gaul_list, ]
    raster_list    <- build_simple_raster_pop(country_shape)
    simple_raster  <- raster_list[['simple_raster']]
    pop_raster     <- raster_list[['pop_raster']]
    # Get population brick for all periods
    pop_raster_annual  <- all_region_pops
    pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
    pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
    pop_raster_annual  <- mask(pop_raster_annual, simple_raster)
    # Get raster of admin2 codes
    admin_level <- 2
    shapes <- admin2_shapes
    cropped_shapes <- crop(shapes, extent(simple_raster), snap="out")
    ## Fix rasterize
    initial_raster <- rasterize_check_coverage(cropped_shapes, simple_raster, field = paste0('ADM', admin_level,'_CODE'))
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])!=0) {
      rasterized_shape <- merge(rasterize_check_coverage(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),], simple_raster, field = paste0('ADM', admin_level,'_CODE')), initial_raster)
    }
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])==0) {
      rasterized_shape <- initial_raster
    }
    admin2_masked_shapes <- mask(x=rasterized_shape, mask=simple_raster)
    # Get raster of admin1 codes
    admin_level <- 1
    shapes <- admin1_shapes
    cropped_shapes <- crop(shapes, extent(simple_raster), snap="out")
    ## Fix rasterize
    initial_raster <- rasterize_check_coverage(cropped_shapes, simple_raster, field = paste0('ADM', admin_level,'_CODE'))
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])!=0) {
      rasterized_shape <- merge(rasterize_check_coverage(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),], simple_raster, field = paste0('ADM', admin_level,'_CODE')), initial_raster)
    }
    if(length(cropped_shapes[!cropped_shapes[[paste0('ADM', admin_level,'_CODE')]]%in%unique(initial_raster),])==0) {
      rasterized_shape <- initial_raster
    }
    admin1_masked_shapes <- mask(x=rasterized_shape, mask=simple_raster)
    # Add items to list
    return_list <- list(country_draws,
                        pop_raster_annual,
                        simple_raster,
                        admin2_masked_shapes,
                        admin1_masked_shapes)
    names(return_list) <- c(paste0('draws_', gaul), paste0('pops_', gaul), paste0('simple_', gaul), paste0('admin2_', gaul), paste0('admin1_', gaul))
    return(return_list)
  }
  reg_gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)
  reg_gaul_list <- reg_gaul_list[!(reg_gaul_list %in% 40762)]
  return_list <- lapply(reg_gaul_list, make_country_list)
  names(return_list) <- paste0('list_', reg_gaul_list)

  ## Return list of draws, pops, templates
  return(return_list)

}
