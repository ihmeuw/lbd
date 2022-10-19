get_is_oos_draws_inla <- function(ind_gp,
                                  ind,
                                  rd,
                                  age = 0,
                                  yrs = 2000:2015,
                                  write.to.file = FALSE,
                                  get.oos = FALSE,
                                  year_col = "original_year",
                                  shapefile_version = "current",
                                  holdouts,
                                  region,
                                  ...) {
  
  
  ## ###################################################################
  ## get_is_oos_draws_inla()
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
  mod.dir <- <<<< FILEPATH REDACTED >>>>
  all.regions = region
  
  outputdir <- <<<< FILEPATH REDACTED >>>> # new outputdir
  
  df_temp <- data.table()
  
  if (length(holdouts > 1) & holdouts != 0) {
    for (i in holdouts) {
      message(paste0("Processing draws for holdout ", i))
      pathaddin <- paste0('_bin', 0, '_', all.regions, '_', i)
      
      df <- fread(<<<< FILEPATH REDACTED >>>>)
      load(<<<< FILEPATH REDACTED >>>>)
      
      if (i == 0) {
        df$fold <- 0
        df$indicator_original <- df[, get(indicator)]
      }
      
      pred_draws <- inv.logit(pred_draws)
      pred_draws <- as.data.table(pred_draws)
      colnames(pred_draws) <- c(paste0("draw", 1:length(pred_draws)))
      
      temp <- as.data.table(cbind(df, pred_draws))
      
      df_temp <- rbindlist(list(df_temp, temp[fold == i,]), use.names = TRUE, fill = TRUE, idcol = FALSE)
    }
    
    df <- copy(df_temp)
    
    df[, (indicator) := indicator_original]
  } else {
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
      if (exists("pass_link_table")) {
        raster_list <- build_simple_raster_pop(subset_shape, link_table = pass_link_table)
      } else {
        raster_list <- build_simple_raster_pop(subset_shape)
      }

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
      # id_raster <- insertRaster(template, matrix(1:(length(cellIdx(template)) * nperiod), ncol = nperiod))
      id_raster <- insertRaster(template, matrix(1:(length(cellIdx(template)))))

      # get cell_pred row ids for each data point
      # for (yy in 1:nperiod) {
      #   this_year <- which(df.r$the_year_col == yrs[yy])
      #   if (length(this_year) == 0) next
      #   df.r[this_year, cell_pred_id := raster::extract(id_raster[[yy]], loc.r[this_year, ])]
      # }
      df.r[, cell_pred_id := raster::extract(id_raster, loc.r)]
      
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
        
        for (y in yrs) {
          message(<<<< FILEPATH REDACTED >>>>)
          load(<<<< FILEPATH REDACTED >>>>)
          df.r[fold == this_fold & the_year_col == y, paste0("draw", 1:ncol(cell_pred)) := as.data.table(cell_pred[cell_pred_id, ])]
        }
      }

      # return combined draws and draws
      return(df.r)
    }))
    
    df <- df.r
  }
    
  df <- as.data.table(df)
  
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
  
  # setnames(df, "ADM0_CODE", "GAUL_CODE")
  if ("ADM0_CODE" %in% colnames(df)) {
    setnames(df, "ADM0_CODE", "GAUL_CODE")
  }
  
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
  
  df_all <- copy(df)
  
  # rename year column back to original
  setnames(df_all, "the_year_col", year_col)
  
  # save combined data and draws to file and return
  if (write.to.file) write.csv(df_all, <<<< FILEPATH REDACTED >>>>, row.names = F)
  return(df_all)
}
