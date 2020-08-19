## plot_stackers() ########################################################

#' Plot maps of stackers, mean raster of covariate, and data informing model
#'
#' @param reg region
#' @param ig indicator group
#' @param ind indicator
#' @param rd run date
#' @param ss vector of stacker names
#' @param yl year list (in vector form, e.g. `c(2000:2015)`)
#' @param zmin minimum value for color scheme (if NULL will calculate from data)
#' @param zmax maximum value for color scheme (if NULL will calculate from data)
#' @param sh_dir `/share` directory, including run date
#' @param highisbad should high values be colored in red ("bad")? Logical.
#' @param o_dir output directory
#' @param individual_countries should individual countries be graphed as well? Logical.
#' @return Writes a series of image fileswith maps of each stacker, mean covariate raster, and
#'         a map of input data for each year-region combination (and year-country if individual_countries = T)
#'         in prespecified folder structure
#' @examples
#' mclapply(Regions, function(r) {
#'   message(paste0("Making stacker maps for region: ", r))
#'   plot_stackers(reg = r, highisbad = F, individual_countries = T)
#' }, mc.cores = 5)
#'
#' # plot stackers & mean outcome
plot_stackers <- function(reg,
                          ig = indicator_group,
                          ind = indicator,
                          rd = run_date,
                          ss = stackers,
                          yl = year_list,
                          predict_years = year_list,
                          zmin = NULL, zmax = NULL,
                          sh_dir = sharedir,
                          highisbad = F,
                          o_dir = out_dir,
                          individual_countries = T,
                          shapefile_version = "current",
                          holdout = 0) {
  gaul_list <- get_adm0_codes(reg, shapefile_version = shapefile_version)

  if (!(file.exists(<<<< FILEPATH REDACTED >>>>)))) {
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = use_premade, shapefile_version = shapefile_version)
    subset_shape <- simple_polygon_list[[1]]
    simple_polygon <- simple_polygon_list[[2]]
  } else {
    load(paste0(<<<< FILEPATH REDACTED >>>>))
  }

  master_shape <- subset_shape

  # Set up output dir
  o_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  dir.create(o_dir, recursive = T, showWarnings = F)

  # Load stacker objects
  load(paste0(<<<< FILEPATH REDACTED >>>>))
  # subset cov_list to stackers
  stacker_list <- cov_list[which(names(cov_list) %in% ss)]

  stacker_list_backup <- stacker_list

  ## make a pathaddin that gets used widely
  pathaddin <- paste0('bin', 0, '_', reg, '_', holdout)
  
  # load mean unraked raster for estimates. Try region rasters first; whole raster & crop if not.
  if (parallel_pred_agg) {
    result_brick <- brick(paste0(<<<< FILEPATH REDACTED >>>>))
    # addLayer(s, r/2, r*2)
    for (i in 2:length(predict_years)) {
      result_brick <- addLayer(result_brick, raster(paste0(<<<< FILEPATH REDACTED >>>>)))
    }
  } else {
    if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
      result_brick <- brick(paste0(<<<< FILEPATH REDACTED >>>>))
    } else if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
      result_brick <- brick(paste0(<<<< FILEPATH REDACTED >>>>))
    } else if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
      result_brick <- brick(paste0(<<<< FILEPATH REDACTED >>>>))
    } else if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
      result_brick <- brick(paste0(<<<< FILEPATH REDACTED >>>>))
      result_brick <- crop(result_brick, master_shape)
      result_brick <- mask(result_brick, master_shape)
    } else if (file.exists(paste0(<<<< FILEPATH REDACTED >>>>))) {
      result_brick <- brick(paste0(<<<< FILEPATH REDACTED >>>>))
      result_brick <- crop(result_brick, master_shape)
      result_brick <- mask(result_brick, master_shape)
    } else {
      stop("Could not find unraked raster .tif.")
    }
  }
  
  prediction_start_index <- which(yl == min(predict_years))

  stacker_length <- length(stacker_list)

  # add mean raster to stacker list
  stacker_list[[ind]] <- result_brick

  # load input data
  input_df <- read.csv(paste0(<<<< FILEPATH REDACTED >>>>),
    stringsAsFactors = F
  ) %>%
    as.data.table()

  input_df[, outcome := get(indicator) / N]

  # KK- define Focal 3 color scale
  f3_cols <- c("#810f7c", "#9ebcda", "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#800026")
  v_cut <- c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.40, 1)
  v <- v_cut[2:length(v_cut)]

  # cats <- cut(input_df$lf_prev, v_cut) %>%
  cats <- cut(input_df[, get(ind)/N], v_cut) %>%
    unique() %>%
    levels()
  col_table <- data.table(colors = f3_cols, categories = cats)
  input_df$cat <- cut(input_df[, get(ind)/N], v_cut)
  input_df[get(ind) == 0, cat := "(0,0.01]"]

  # Cap N at 90% for plot interpretability
  input_df[, cap_N := ifelse(N >= quantile(N, 0.9), quantile(N, 0.9), N)]

  full_df <- input_df

  # Define function to save maps
  save_maps <- function(input_df, stacker_list, master_shape, result_brick, zmin, zmax, yl, predict_years, ind, ig, sh_dir, highisbad, o_dir, ctry = NULL, adm0 = NULL, shapefile_version) {
    if (!is.null(ctry)) {
      input_df <- subset(input_df, country == ctry)
      master_shape <- subset(master_shape, ADM0_CODE == adm0)
      stacker_list <- lapply(stacker_list, function(x) {
        x <- suppressMessages(crop(x, extent(master_shape)))
        x <- suppressMessages(mask(x, master_shape))
        extent(x) <- extent(master_shape)
        return(x)
      })
    }

    master_shape_df <- suppressMessages(fortify(master_shape))

    # define mins / maxes from global values across all years in stacker_list & input data
    if (is.null(zmin)) {
      zmin <- min(
        sapply(stacker_list, function(x) min(minValue(x))),
        min(input_df$outcome)
      )
    }
    if (is.null(zmax)) {
      zmax <- max(
        sapply(stacker_list, function(x) max(maxValue(x))),
        max(input_df$outcome)
      )
    }

    # Check # years correct
    if (nlayers(result_brick) != length(predict_years)) stop("Number of mean raster brick layers does not equal length of year list")

    # rearrange
    for (i in 1:length(predict_years)) {
      message(paste0("   year ", predict_years[i], "..."))

      make_gg_map <- function(rbrick, title, i, color_table = col_table) {
        r_df <- rbrick[[i]] %>%
          as("SpatialPixelsDataFrame") %>%
          as.data.frame()
        names(r_df) <- c("vals", "x", "y")
        r_df$cut_vals <- cut(r_df$vals, v_cut)

        col_vec <- color_table$colors
        names(col_vec) <- color_table$categories

        gg_result <- ggplot() +
          geom_raster(data = r_df, aes(x = x, y = y, fill = cut_vals), show.legend = T) +
          coord_equal() +
          theme_void() +
          labs(title = title) +
          scale_fill_manual(values = col_vec, labels = color_table$categories, name = "ICT Prevalence") +
          geom_path(
            data = master_shape_df,
            aes(x = long, y = lat, group = group)
          )

        return(gg_result)
      }

      gg_stackers_results <- lapply(1:length(stacker_list), function(n) {
        the_title <- paste0(names(stacker_list)[n], ": ", predict_years[i])
        the_rbrick <- stacker_list[[n]]
        if (n == (stacker_length + 1)) {
          return(make_gg_map(the_rbrick, the_title, i))
        } else {
          return(make_gg_map(the_rbrick, the_title, i))
        }
      })

      # Make data plot
      col_vec <- col_table$colors
      names(col_vec) <- col_table$categories

      gg_data <- ggplot() +
        geom_point(data = subset(input_df, year == predict_years[i]), aes(x = longitude, y = latitude, size = cap_N, alpha = weight, color = cat)) +
        coord_equal() +
        theme_void() +
        scale_color_manual(values = col_vec, labels = col_table$categories, name = "ICT Prevalence") +
        geom_path(
          data = master_shape_df,
          aes(x = long, y = lat, group = group)
        ) +
        scale_size_continuous(limits = c(NA, max(full_df$cap_N))) +
        labs(title = paste0("data: ", predict_years[i]))

      gg_stackers_results[[length(gg_stackers_results) + 1]] <- gg_data

      # Use first legend only
      the_legend <- g_legend(gg_stackers_results[[1]])
      gg_stackers_results <- lapply(gg_stackers_results, function(x) return(x + theme(legend.position = "none")))

      if (is.null(ctry)) {
        reg_dir <- paste0(o_dir, reg, "/")
        dir.create(reg_dir, recursive = T, showWarnings = F)
        fn <- paste0(<<<< FILEPATH REDACTED >>>>)
      } else if (!is.null(ctry)) {
        ctry_dir <- paste0(o_dir, ctry, "/")
        dir.create(ctry_dir, recursive = T, showWarnings = F)
        fn <- paste0(<<<< FILEPATH REDACTED >>>>)
      }

      png(
        filename = fn,
        width = 16, height = 9, units = "in",
        res = 200, pointsize = 10,
        type = "cairo-png"
      )

      multiplot(
        plotlist = gg_stackers_results,
        cols = ceiling(length(gg_stackers_results) / 2),
        legend = the_legend
      )

      dev.off()
    }
  }

  # Save map for entire region
  save_maps(input_df, stacker_list, master_shape, result_brick, zmin, zmax, yl, predict_years, ind, ig, sh_dir, highisbad, o_dir, ctry = NULL, adm0 = NULL, shapefile_version = shapefile_version)

  # Save maps for individual countries
  if (individual_countries == T) {
    gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
    gaul_list <- data.table(GAUL_CODE = get_adm0_codes(reg, shapefile_version = shapefile_version))
    gaul_list <- merge(gaul_list, gaul_to_loc_id, by = "GAUL_CODE")
    for (c in unique(gaul_list$ADM_CODE)) {
      if (!(get_adm0_codes(gaul_list[ADM_CODE == c, "ihme_lc_id"], shapefile_version = shapefile_version) %in% master_shape$ADM0_CODE)) {
        message(paste0("No shapes in master_shape corresponding to admin 0 code for ", c, " - skipping..."))
      } else {
        message(paste0("Saving stacker maps for country: ", gaul_list[ADM_CODE == c, "ihme_lc_id"]))
        save_maps(input_df, stacker_list, master_shape, result_brick, zmin, zmax, yl, predict_years, ind, ig, sh_dir, highisbad, o_dir, ctry = gaul_list[ADM_CODE == c, ihme_lc_id], adm0 = c, shapefile_version = shapefile_version)
      }
    }
  }
}
