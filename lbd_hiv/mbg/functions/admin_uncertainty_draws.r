####################################################################################################
## Make table of uncertainty in admin draws
####################################################################################################

make_admin_change_uncertainty_tables <-
  function(run_date,
           indicator,
           indicator_group,
           regions,
           core_repo = paste0("/lbd_core/"),
           outdir,
           raked = T,
           years = c(2000, 2017),
           age   = 0, 
           sex   = 0) {

    # Read in year list from config
    share_dir <- paste("<<<< FILEPATH REDACTED >>>>")
    this_config <- fread(paste0(share_dir, '/config.csv'))
    year_list <- this_config[V1 == 'year_list', V2]
    if (is.character(year_list)) year_list <- eval(parse(text = year_list))

    # Load in admin draws
    if (raked == T) {
      load(paste0("<<<< FILEPATH REDACTED >>>>"))
    } else {
      load(paste0("<<<< FILEPATH REDACTED >>>>"))
    }
    lvl <- c(0, 1, 2)

    # Generate csv files that show the upper and lower for admin changes
    lapply(lvl, function(x) {

      # Grab correct admin estimate
      admin_pred <- get(paste0('admin_', x))
      admin_code <- paste0('ADM', x, '_CODE')

      # Grab admin codes
      ad_codes <-
        admin_pred %>% filter(year == max(year_list)) %>% pull(admin_code)

      # Grab draws form first year
      draws_initial <-
        admin_pred %>%
        filter(year == min(years)) %>%
        dplyr::select(grep("V[0-9]*", names(admin_pred)))

      # Grab draws from last year
      draws_final <-
        admin_pred %>%
        filter(year == max(years)) %>%
        dplyr::select(grep("V[0-9]*", names(admin_pred)))

      # Calculate difference in percentage space (mean, upper, lower)
      draw_dif <- (as.matrix(draws_final) - as.matrix(draws_initial)) * 100
      draw_dif <- as.data.table(draw_dif)
      draw_upper <-
        apply(draw_dif, 1, function(x)
          quantile(x, 0.975, na.rm = T))
      draw_mean <- apply(draw_dif, 1, mean, na.rm = T)
      # Combine upper lower and mean
      draw_summary <-
        apply(draw_dif, 1, function(x)
          quantile(x, 0.025, na.rm = T)) %>%
        data.table(ad_codes,
                   lower = .,
                   upper = draw_upper,
                   mean = draw_mean)
      # add a column indicating what % of draws are positive (an increase in prevalence)
      perc_pos_draws <- list()
      for (i in 1:nrow(draw_dif)) {
        c <- draw_dif[i,]
        c <- unlist(c)
        pos <- c[which(c > 0)]
        perc <- length(pos) / length(c)
        perc_pos_draws[[i]] <- perc
      }
      perc_pos_draws <- unlist(perc_pos_draws)
      draw_summary <- cbind(draw_summary, perc_pos_draws)
      # Label each admin_difference with country name
      if (x == 0) {
        admin_dif <-
          sp_hierarchy_list %>%
          dplyr::select(ADM0_CODE, ADM0_NAME) %>% unique %>%
          mutate(ADM0_CODE = as.numeric(as.character(ADM0_CODE))) %>%
          right_join(draw_summary, by = c("ADM0_CODE" = "ad_codes")) %>%
          arrange(mean) %>%
          dplyr::rename(lower_diff = lower,
                        upper_diff = upper,
                        mean_diff = mean)
        
      } else if (x == 1) {
        admin_dif <-
          sp_hierarchy_list %>%
          dplyr::select(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME) %>% unique %>%
          mutate(ADM1_CODE = as.numeric(as.character(ADM1_CODE))) %>%
          right_join(draw_summary, by = c("ADM1_CODE" = "ad_codes")) %>%
          arrange(mean) %>%
          dplyr::rename(lower_diff = lower,
                        upper_diff = upper,
                        mean_diff = mean)
      } else {
        admin_dif <-
          sp_hierarchy_list %>%
          dplyr::select(ADM0_CODE,
                        ADM0_NAME,
                        ADM1_CODE,
                        ADM1_NAME,
                        ADM2_CODE,
                        ADM2_NAME) %>%
          unique() %>%
          mutate(ADM2_CODE = as.numeric(as.character(ADM2_CODE))) %>%
          right_join(draw_summary, by = c("ADM2_CODE" = "ad_codes")) %>%
          arrange(mean) %>%
          dplyr::rename(lower_diff = lower,
                        upper_diff = upper,
                        mean_diff = mean)
      }
      write.csv(admin_dif, paste0(outdir, "/pred_derivatives/admin_summaries/admin", x, "_diff_uncertainty_", min(years), "_", max(years), "_bin", age, ".csv"))
    })
  }


####################################################################################################
## Create raster of pixel level differences
####################################################################################################

make_pixel_difference_uncertainty_rasters <-
  function(run_date,
           indicator,
           indicator_group,
           regions,
           summary_stat,
           core_repo = paste0("/lbd_core/"),
           outdir,
           age = 0,
           sex = 0) {
    
    # Read in year list from config
    share_dir <- paste("<<<< FILEPATH REDACTED >>>>")
    this_config <- fread(paste0(share_dir, '/config.csv'))
    year_list <- this_config[V1 == 'year_list', V2]
    if (is.character(year_list)) year_list <- eval(parse(text = year_list))

    # Create folder in outdir for these pixel difference rasters
    outdir_raster <- paste0(outdir, "/pixel_difference_rasters/")
    dir.create(outdir_raster)

    # Over all regions create difference raster --------------------------------------------------------------
    lapply(regions, function(reg) {
      
      holdout <- 0
      pathaddin <- paste0('_bin', age, '_', reg, '_', holdout)

      # Load cell pred
      load(paste0(share_dir, indicator, "_raked_cell_draws_eb", pathaddin, ".RData"))

      # Extract first/last year to make different cell pred
      first_year_rows <- nrow(raked_cell_pred) / length(year_list)
      first_year      <- raked_cell_pred[1:first_year_rows,]
      last_year       <- raked_cell_pred[(nrow(raked_cell_pred) - first_year_rows + 1):nrow(raked_cell_pred), ]
      dif_pred        <- (last_year - first_year) * 100
      rm(first_year, last_year, raked_cell_pred)
      perc_pos_draws <- list()
      for(i in 1:NROW(dif_pred)){
        c <- dif_pred[i, ]
        c <- unlist(c)
        pos <- c[which(c>0)]
        perc <- length(pos)/length(c)
        perc_pos_draws[[i]]<- perc
      }
      perc_pos_draws <- unlist(perc_pos_draws)
      # Load simple raster
      simple_polygon_list <-
        load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version),
                            shapefile_version = shapefile_version,
                            buffer = 0.4,
                            subset_only = FALSE)
      subset_shape   <- simple_polygon_list[[1]]
      raster_list    <- build_simple_raster_pop(subset_shape)
      simple_raster  <- raster_list[['simple_raster']]
      rm(simple_polygon_list, subset_shape, raster_list)

      # Apply summary statistic
      for (stat in summary_stat) {
        summ <- apply(dif_pred, 1, stat)
        summ <- insertRaster(simple_raster,  matrix(summ,  ncol = 1))
        writeRaster(summ, file = paste0(outdir_raster, '/', indicator, '_', reg, '_bin', age, 
                                        "_raked_2000_2017_dif_", stat, "_raster.tif"),
                    format = 'GTiff',overwrite = TRUE)}
      perc_pos_draws <- insertRaster(simple_raster,  matrix(perc_pos_draws,  ncol = 1))
      writeRaster(perc_pos_draws, file = paste0(outdir_raster, '/', indicator, '_', reg, '_bin', age,
                                      "_raked_2000_2017_dif_perc_pos_draws_raster.tif"),
                  format = 'GTiff',overwrite = TRUE)
    })

    # Combine all raster maps to create over all of Africa
    lapply(summary_stat, function(stat) {
      reg_list <- lapply(regions, function(reg) {
        brick(paste0(outdir_raster, '/', indicator, '_', reg, '_bin', age,  "_raked_2000_2017_dif_", stat, "_raster.tif"))
      })
    rlist <- do.call(raster::merge, unname(reg_list))
    
    

    # Write raster of combined regions
    writeRaster(rlist,
                file = paste0(outdir_raster,'/',indicator,'_',"africa", '_bin', age, "_raked_2000_2017_dif_",stat,"_raster.tif"),
                format = 'GTiff',
                overwrite = TRUE)
    })
    
    perc_draws <- lapply(regions, function(reg) {
      brick(paste0(outdir_raster, '/', indicator, '_', reg, '_bin', age, "_raked_2000_2017_dif_perc_pos_draws_raster.tif"))
    })
    perc_d <- do.call(raster::merge, unname(perc_draws))
    writeRaster(perc_d,
                file = paste0(outdir_raster,'/',indicator,'_',"africa", '_bin', age, "_raked_2000_2017_dif_perc_pos_draws_raster.tif"),
                format = 'GTiff',
                overwrite = TRUE)
  }


