## subnational_ts_plots ################################################

# Function to create subnational time-series plots with accompanying maps
#
# This function is meant to create three sets of PDFs that show
# time series plots of administrative-level estimates of your
# modeled quantity of interest:
#
# 1) Plots of country-level estimates (single PDF), optionally one page per region
# 2) Plots of admin1-level estimates (single PDF), one page per country
# 3) Plots of admin2-level estimates (one PDF per country), one page per admin1
#
# As currently written this function is quite repetitive. This allows for each of the
# maps to be tweaked more easily, but in the future it would be better to
# encapsulate some of the repeated bits in flexible functions that can be
# called for each of the categories of maps created above. This would
# also potentially allow for the use of non-standard shapefiles, which isn't
# well supported by the function at the moement.
#
# All data tables referenced below (and accompanying shapefiles) must have the
# standard ADM0_CODE, ADM0_NAME, etc. naming conventions and have those codes
# and names included as columns / shapefile attributes. In addition,
# the administrative estimates must have "mean", "upper", and "lower" columns
# to allow uncertainty intervals to be plotted.
#
# @author Jon Mosser, \email{jmosser@uw.edu}
#
# @param ad0_df admin-0 level data.frame or data.table from the standard aggregation code
# @param ad1_df admin-1 level data.frame or data.table from the standard aggregation code
# @param ad2_df admin-2 level data.frame or data.table from the standard aggregation code.
# @param ad0_shape admin0 SPDF object that corresponds to `ad0_df` above.
#                  If none provided, loads from GAUL as default
# @param ad1_shape admin1 SPDF object that corresponds to `ad1_df` above.
#                  If none provided, loads from GAUL as default
# @param ad2_shape admin2 SPDF object that corresponds to `ad2_df` above.
#                  If none provided, loads from GAUL as default
# @param ind_title title for your indicator, e.g. "DPT3 Coverage".
#                  For use in plots in multiple places.  String.
# @param out_dir output directory. String.
# @param out_filename_format output filename format, in `sprintf()` format
#                            with `%s` representing the place where the
#                            geographic identifier will go. String.
# @param val_range range of values for your indicator. Default: [0,1]
#                  Vector of length 2 (`c(min, max)`)
# @param highisbad if T, higher values will be more red. Boolean.
# @param ad0_map_regions if you want to have the national-level time series
#                        broken down by modeling region, you can pass a vector
#                        of region names here.  These will be processed by
#                        `get_gaul_codes()` so need to be identified in that
#                        function.
# @param ad0_map_region_titles corresponding to ad0_map_regions, you can pass
#                              titles for each of the map regions here (so that
#                              `essa` is identified as "Eastern Sub-Saharan Africa",
#                              for instance).  Vector of same length of ad0_map_regions
# @param ad1_map_countries If specified, which countries do you want to plot at the admin1/admin2 level.
#                          If null, the default is plotting all countries in the specified ad0_map_regions
# @param plot_levels which administrative levels do you want to plot? Vector containing
#                    any combination of "ad0", ad1", and/or "ad2".
# @param title_plot_size if you want, you can specify a font size to be used for ind_title in text_element.
#                        If none provided, default value is assigned
# @param title_grob_size if you want, you can specify a font size to be used for ind_title in text_grob.
#                        If none provided, default value is assigned
# @param multiple_runs should the function plot multiple model runs/ multiple indicators on the same plot? Boolean
# @param plot_data should the function plot the collapsed admin input data? Boolean
# @param ad0_data Collapsed admin0 input data, obtained using input_aggregate_admin function found below
# @param ad1_data Collapsed admin1 input data, obtained using input_aggregate_admin function found below
# @param ad2_data Collapsed admin2 input data, obtained using input_aggregate_admin function found below
# @param verbose should the function update progress along the way? Boolean.
# @param shapefile_version specifies the shapefile version to be used. `current` or date `yyyy_mm_dd`
# @param tol tolerance to be used when simplifying polygons. default `0.001` creates orphaned holes in GADM shapefile, can use `0.0001` to fix but the function slows down and outputs are larger.
# @return
# Saves a series of PDFs as described above in `out_dir`

subnational_ts_plots <- function(ad0_df,
                                 ad1_df,
                                 ad2_df,
                                 ad0_shape = NULL,
                                 ad1_shape = NULL,
                                 ad2_shape = NULL,
                                 ind_title,
                                 out_dir,
                                 out_filename_format = "subnational_ts_plots_%s.pdf",
                                 val_range = c(0, 1),
                                 highisbad = T,
                                 ad0_map_regions = NULL,
                                 ad0_map_region_titles = NULL,
                                 ad1_map_countries = NULL,
                                 plot_levels = c("ad0", "ad1", "ad2"),
                                 multiple_runs = F,
                                 plot_data = F,
                                 ad0_data = NULL,
                                 ad1_data = NULL,
                                 ad2_data = NULL,
                                 title_plot_size = NULL,
                                 title_grob_size = NULL,
                                 verbose = F,
                                 mda_data,
                                 plot_stacker_trends = T,
                                 modeling_shapefile_version = NULL,
                                 tol = 0.0001,
                                 ind_gp = indicator_group,
                                 ind = indicator,
                                 rd = run_date) {

  ################################################################################
  # 0. SETUP #####################################################################
  # Load packages, define themes, and set up objects needed for plotting #########
  ################################################################################

  # Load libraries and functions -------------------------------------------------

  library(ggrepel)
  library(gridExtra)
  library(ggplot2)

  if (is.null(modeling_shapefile_version)) {
    plot_stacker_trends <- plot_stacker_trends
    configurations <- fread(<<<< FILEPATH REDACTED >>>>)
    modeling_shapefile_version <- configurations[V1 == "modeling_shapefile_version", V2]
    stacked_fixed_effects <- configurations[V1 == "stacked_fixed_effects", V2]
    stacked_fixed_effects <- trim(strsplit(stacked_fixed_effects, "\\+")[[1]])
    message(paste0("Using modeling_shapefile_version specified in config: ", modeling_shapefile_version))
  }

  # Check to see if plot_data = T that admin data is attached
  if (plot_data == T) {
    if ("ad0" %in% plot_levels & is.null(ad0_data)) {
      stop(paste(
        "To plot data, you must include aggregated ad0_data.",
        "See input_aggregate_admin function"
      ))
    }
    if ("ad1" %in% plot_levels & is.null(ad1_data)) {
      stop(paste(
        "To plot data, you must include aggregated ad1_data.",
        "See input_aggregate_admin function"
      ))
    }
    if ("ad2" %in% plot_levels & is.null(ad2_data)) {
      stop(paste(
        "To plot data, you must include aggregated ad2_data.",
        "See input_aggregate_admin function"
      ))
    }
  }

  if (plot_data == T & multiple_runs == T) {
    message(paste(
      "You are choosing to plot multiple indicators/model runs with data.\n",
      "This is usually done if plotting multiple model runs from the same indicator"
    ))
  }

  # If plotting multiple model runs, make sure there is a column title 'run'
  if (multiple_runs == T) {
    if (!"run" %in% names(ad0_df)) stop("Must have a column title 'run' that differentiates between model runs if plotting multiple indicators/model runs. See add_run_label function")
    if (!"run" %in% names(ad1_df)) stop("Must have a column title 'run' that differentiates between model runs if plotting multiple indicators/model runs. See add_run_label function")
    if (!"run" %in% names(ad2_df)) stop("Must have a column title 'run' that differentiates between model runs if plotting multiple indicators/model runs. See add_run_label function")
  }

  plot_stacker_trends <- plot_stacker_trends
  configurations <- fread(<<<< FILEPATH REDACTED >>>>)
  stackers <- configurations[V1 == "stacked_fixed_effects", V2]
  stackers <- trim(strsplit(stackers, "\\+")[[1]])

  # Create directories -----------------------------------------------------------

  dir.create(out_dir, showWarnings = F)

  # Prepare data and shapefiles --------------------------------------------------

  # Merge on ihme_lc_ids
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = modeling_shapefile_version)
  gaul_to_loc_id <- subset(gaul_to_loc_id, select = c("GAUL_CODE", "ihme_lc_id"))
  setnames(gaul_to_loc_id, "GAUL_CODE", "ADM0_CODE")

  if (plot_data == T) {
    ad0_data$ihme_lc_id <- NULL
    ad1_data$ihme_lc_id <- NULL
    ad2_data$ihme_lc_id <- NULL

    ad0_data <- merge(copy(as.data.table(ad0_data)), gaul_to_loc_id, all.x = T, all.y = F, by = "ADM0_CODE")
    ad1_data <- merge(copy(as.data.table(ad1_data)), gaul_to_loc_id, all.x = T, all.y = F, by = "ADM0_CODE")
    ad2_data <- merge(copy(as.data.table(ad2_data)), gaul_to_loc_id, all.x = T, all.y = F, by = "ADM0_CODE")
  }

  # Set up default values for title font size if not specified
  if (is.null(title_plot_size)) title_plot_size <- 20
  if (is.null(title_grob_size)) title_grob_size <- 30

  if (verbose == T) message("Simplifying shapes")

  # Wrapper for gSimplify() to allow it to work with SPDFs
  simplify_spdf <- function(the_spdf, ...) {
    simple_spdf <- gSimplify(the_spdf, topologyPreserve = T, ...)
    return(SpatialPolygonsDataFrame(simple_spdf, the_spdf@data))
  }

  # Subset and simplify shapefiles for memory & speed
  ad0_codes <- unique(ad0_df$ADM0_CODE)
  # If any region of Africa is being used, use full map of Africa so it does not look disjointed
  if (any(ad0_map_regions %in% c("cssa", "wssa", "essa", "sssa", "name"))) ad0_codes <- union(ad0_codes, get_adm0_codes("Africa", shapefile_version = modeling_shapefile_version))
  subset_codes <- ad0_codes

  # Use simplified GAUL shapefiles as default if none provided, if provided simplify the shapefile for speed
  if (is.null(ad0_shape)) {
    ad0_shape_simple <- readOGR(get_admin_shapefile(admin_level = 0, raking = F, version = modeling_shapefile_version)) %>%
      subset(ADM0_CODE %in% subset_codes) %>%
      simplify_spdf(tol = tol)
  } else {
    ad0_shape_simple <- subset(ad0_shape, ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  }

  if (is.null(ad1_shape)) {
    ad1_shape_simple <- readOGR(get_admin_shapefile(admin_level = 1, raking = F, version = modeling_shapefile_version)) %>%
      subset(ADM0_CODE %in% subset_codes) %>%
      simplify_spdf(tol = tol)
  } else {
    ad1_shape_simple <- subset(ad1_shape, ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  }

  if (is.null(ad2_shape)) {
    ad2_shape_simple <- readOGR(get_admin_shapefile(admin_level = 2, raking = F, version = modeling_shapefile_version)) %>%
      subset(ADM0_CODE %in% subset_codes) %>%
      simplify_spdf(tol = tol)
  } else {
    ad2_shape_simple <- subset(ad2_shape, ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  }

  # Make sure that the admin codes are converted to numeric if they are factors
  if (is.factor(ad0_shape_simple@data$ADM0_CODE)) ad0_shape_simple@data$ADM0_CODE <- as.numeric(levels(ad0_shape_simple@data$ADM0_CODE))[ad0_shape_simple@data$ADM0_CODE]
  if (is.factor(ad1_shape_simple@data$ADM1_CODE)) ad1_shape_simple@data$ADM1_CODE <- as.numeric(levels(ad1_shape_simple@data$ADM1_CODE))[ad1_shape_simple@data$ADM1_CODE]
  if (is.factor(ad2_shape_simple@data$ADM2_CODE)) ad2_shape_simple@data$ADM2_CODE <- as.numeric(levels(ad2_shape_simple@data$ADM2_CODE))[ad2_shape_simple@data$ADM2_CODE]

  # renaming the ADMX_NAME columns using the names in ad2_shape, because this functions matches
  # everything on names and things get dropped if the names don't line up perfectly
  # This is an alternative to rewriting the entire function to match on coes
  admin_shp_data_adm0 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm0 <- unique(admin_shp_data_adm0[, c("ADM0_CODE", "ADM0_NAME")])
  admin_shp_data_adm0$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm0$ADM0_CODE))
  admin_shp_data_adm0$ADM0_NAME <- as.character(admin_shp_data_adm0$ADM0_NAME)
  gaul_to_loc_id_adm0 <- merge(gaul_to_loc_id, admin_shp_data_adm0, by = "ADM0_CODE")

  admin_shp_data_adm1 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm1 <- unique(admin_shp_data_adm1[, c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME")])
  admin_shp_data_adm1$ADM1_CODE <- as.integer(as.character(admin_shp_data_adm1$ADM1_CODE))
  admin_shp_data_adm1$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm1$ADM0_CODE))
  admin_shp_data_adm1$ADM1_NAME <- as.character(admin_shp_data_adm1$ADM1_NAME)
  admin_shp_data_adm1$ADM0_NAME <- as.character(admin_shp_data_adm1$ADM0_NAME)
  gaul_to_loc_id_adm1 <- merge(gaul_to_loc_id, admin_shp_data_adm1, by = "ADM0_CODE")

  admin_shp_data_adm2 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm2 <- unique(admin_shp_data_adm2[, c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "ADM2_CODE", "ADM2_NAME")])
  admin_shp_data_adm2$ADM2_CODE <- as.integer(as.character(admin_shp_data_adm2$ADM2_CODE))
  admin_shp_data_adm2$ADM1_CODE <- as.integer(as.character(admin_shp_data_adm2$ADM1_CODE))
  admin_shp_data_adm2$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm2$ADM0_CODE))
  admin_shp_data_adm2$ADM2_NAME <- as.character(admin_shp_data_adm2$ADM2_NAME)
  admin_shp_data_adm2$ADM1_NAME <- as.character(admin_shp_data_adm2$ADM1_NAME)
  admin_shp_data_adm2$ADM0_NAME <- as.character(admin_shp_data_adm2$ADM0_NAME)
  gaul_to_loc_id_adm2 <- merge(gaul_to_loc_id, admin_shp_data_adm2, by = "ADM0_CODE")

  ad0_df[, c("ADM0_NAME")] <- NULL
  ad1_df[, c("ADM0_NAME", "ADM1_NAME", "ADM0_CODE")] <- NULL
  ad2_df[, c("ADM0_NAME", "ADM1_NAME", "ADM2_NAME", "ADM0_CODE", "ADM1_CODE")] <- NULL

  ad0_df <- merge(copy(as.data.table(ad0_df)), gaul_to_loc_id_adm0, all.x = T, all.y = F, by = "ADM0_CODE")
  ad1_df <- merge(copy(as.data.table(ad1_df)), gaul_to_loc_id_adm1, all.x = T, all.y = F, by = "ADM1_CODE")
  ad2_df <- merge(copy(as.data.table(ad2_df)), gaul_to_loc_id_adm2, all.x = T, all.y = F, by = "ADM2_CODE")

  # Add simple world shapefile if region is specified as Africa
  if ("ad0" %in% plot_levels & "Africa" %in% ad0_map_regions) world_shape_simple <- readRDS(<<<< FILEPATH REDACTED >>>>)

  # Define some utility functions -----------------------------------------------------------------------------------------------------
  # Function to fix aspect ratios
  aspect_ratio_plot <- function(plt, aspect_ratio = 1, expand = 0.05) {
    xlims <- ggplot_build(plt)$layout$panel_scales$x[[1]]$range$range
    ylims <- ggplot_build(plt)$layout$panel_scales$y[[1]]$range$range

    # Try a different way if null - varies by ggplot build
    if (is.null(xlims)) xlims <- ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range
    if (is.null(ylims)) ylims <- ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range

    xmid <- mean(xlims)
    xrange <- xlims[2] - xlims[1]
    ymid <- mean(ylims)
    yrange <- ylims[2] - ylims[1]

    if (xrange / yrange < aspect_ratio) {
      new_yrange <- yrange * (1 + expand)
      new_xrange <- new_yrange * aspect_ratio
    }

    if (xrange / yrange > aspect_ratio) {
      new_xrange <- xrange * (1 + expand)
      new_yrange <- new_xrange / aspect_ratio
    }

    new_xlim <- c(xmid - 0.5 * new_xrange, xmid + 0.5 * new_xrange)
    new_ylim <- c(ymid - 0.5 * new_yrange, ymid + 0.5 * new_yrange)

    return(plt + expand_limits(
      x = new_xlim,
      y = new_ylim
    ))
  }

  # Title generate makes the title Grob that goes in the right hand corner of the pdf
  title_generate <- function(ind_title,
                               year_list,
                               admin = 0,
                               ad0_reg_title = "",
                               ctry = "",
                               ad1_df_ad1 = "") {
    arrangeGrob(textGrob("", gp = gpar(fontsize = 30)),
      textGrob(str_wrap(ind_title, 18),
        gp = gpar(fontsize = title_grob_size, fontface = "bold")
      ),
      textGrob("", gp = gpar(fontsize = 10)),
      textGrob(if (admin == 0) ifelse(ad0_reg_title == "All countries", NULL, ad0_reg_title) else ctry,
        gp = gpar(fontsize = 20)
      ),
      textGrob("", gp = gpar(fontsize = 10)),
      textGrob(if (admin == 0) "National estimates" else if (admin == 1) "By First-level Administrative Unit" else paste0("Admin 1: ", unique(ad1_df_ad1$ADM1_NAME)),
        gp = gpar(fontsize = 15)
      ),
      textGrob("", gp = gpar(fontsize = 10)),
      textGrob(paste0(min(year_list), " - ", max(year_list)),
        gp = gpar(fontsize = 18)
      ),
      ncol = 1,
      heights = c(30, 25, 10, 25, 10, 15, 10, 20)
    )
  }

  # Plot overlay defines where each object (ggplot and title) are placed on the pdf.
  # This depends on if data is being plotted, and if multiple runs are used
  plot_overlay <- function(plot_data, multiple_runs) {
    if (plot_data == T) {
      if (multiple_runs == T) {
        lay <- rbind(
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 4, 4, 4, NA),
          c(1, 1, 1, 1, 1, 4, 4, 4, NA),
          c(1, 1, 1, 1, 1, 4, 4, 4, 5)
        )
      } else {
        lay <- rbind(
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 2, 2, 3, 3),
          c(1, 1, 1, 1, 1, 4, 4, 4, 4),
          c(1, 1, 1, 1, 1, 4, 4, 4, 4),
          c(1, 1, 1, 1, 1, 4, 4, 4, 4)
        )
      }
    } else {
      lay <- rbind(
        c(1, 1, 1, 1, 2, 2, 3, 3),
        c(1, 1, 1, 1, 2, 2, 3, 3),
        c(1, 1, 1, 1, 4, 4, 4, 4),
        c(1, 1, 1, 1, 4, 4, 4, 4),
        c(1, 1, 1, 1, 4, 4, 4, 4),
        c(1, 1, 1, 1, NA, NA, NA, 5)
      )
    }
    return(lay)
  }

  # Custom empty theme
  theme_empty <- theme_classic() +
    theme(
      axis.line = element_blank(), axis.text.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )

  # Add IHME logo for plotting later
  ihme_logo <- png::readPNG(<<<< FILEPATH REDACTED >>>>)
  ihme_grob <- rasterGrob(ihme_logo)

  # Add slightly different version for plotting multiple runs
  ihme_grob_multiple <- rasterGrob(ihme_logo, y = 0.2)

  # Only list countries that are within the specified regions
  if (!is.null(ad0_map_regions)) {
    ad0_code_list <- lapply(ad0_map_regions, get_adm0_codes, shapefile_version = modeling_shapefile_version)
    names(ad0_code_list) <- ad0_map_region_titles
  } else {
    ad0_code_list <- list(unique(ad0_df$ADM0_CODE))
    names(ad0_code_list) <- "All countries"
  }

  # Only list countries that are within the given regions
  countries <- sort(unique(ad0_df[ADM0_CODE %in% unlist(ad0_code_list)]$ADM0_NAME))

  # Get years
  year_list <- unique(ad0_df$year)

  ################################################################################
  # I. NATIONAL-LEVEL PLOTS FOR ALL COUNTRIES ####################################
  ################################################################################

  if ("ad0" %in% plot_levels) {
    pdf_filename <- paste0(out_dir, sprintf(out_filename_format, "all_countries"))

    if (verbose == T) {
      message("\n############################################################")
      message("############################################################")
      message(paste0("Plotting national-level estimates by region"))
      message(paste0("  Writing output to ", pdf_filename))
    }

    pdf(
      file = pdf_filename,
      height = 10,
      width = 18
    )

    for (i in 1:length(ad0_code_list)) {

      # Set up names & titles for this region ---------------------------------------
      ad0_reg_codes <- ad0_code_list[[i]]
      ad0_reg_title <- names(ad0_code_list)[i]

      ad0_df_reg <- subset(ad0_df, ADM0_CODE %in% ad0_reg_codes)
      if (plot_data == T) ad0_data_reg <- subset(ad0_data, ADM0_CODE %in% ad0_reg_codes)

      if (verbose == T) message(paste0("  ", ad0_reg_title))

      # First, time series plot of all admin0s in the region
      n_ad0s <- length(unique(ad0_df_reg$ADM0_CODE))

      # Try to wrap names of countries as needed
      wrap_width <- 18
      if (n_ad0s > 0 & n_ad0s <= 16) wrap_width <- 18
      if (n_ad0s > 16 & n_ad0s <= 25) wrap_width <- 12
      if (n_ad0s > 25) wrap_width <- 9

      # Set up plot names; include national-level estimates; ensure ordering correct.
      # Use location ID as label modeling over more than 25 countries
      if (n_ad0s < 25) {
        ad0_df_reg[, plot_name := stringr::str_wrap(ADM0_NAME, width = wrap_width)]
      } else {
        ad0_df_reg[, plot_name := stringr::str_wrap(ihme_lc_id, width = wrap_width)]
      }

      # Admin 0 time series plot --------------------------------------------------------------
      if (multiple_runs == F) gg_ad0_ts <- time_series(ad0_df_reg, admin = 0, val_range, title_plot_size, ind_title, plot_stkrs = plot_stacker_trends, stkrs = stackers)
      if (multiple_runs == T) gg_ad0_ts <- time_series_multiple(ad0_df_reg, admin = 0, val_range, title_plot_size, ind_title)

      # A locator map --------------------------------------------------------------
      gg_locator_map <- location_map_draw(subset(ad0_shape_simple, ADM0_CODE %in% ad0_reg_codes), ad0_shape_simple)
      if (ad0_reg_title == "Africa") gg_locator_map <- location_map_draw(subset(world_shape_simple, ADM0_CODE %in% ad0_reg_codes), world_shape_simple)

      # A labeled map by countries -------------------------------------------------
      ad0_shape_reg <- subset(ad0_shape_simple, ADM0_CODE %in% ad0_reg_codes)

      centroids <- gCentroid(ad0_shape_reg, byid = T, id = ad0_shape_reg$ADM0_CODE)
      centroids <- as.data.frame(centroids) %>%
        cbind(rownames(.), .) %>%
        as.data.table(.) %>%
        setnames(., names(.), c("ADM0_CODE", "x", "y"))
      centroids$ADM0_CODE <- as.numeric(levels(centroids$ADM0_CODE))[centroids$ADM0_CODE]
      centroids <- merge(centroids, as.data.table(ad0_shape_reg), by = "ADM0_CODE")

      plot_df <-
        ad0_df_reg %>%
        filter(year == max(year_list)) %>%
        dplyr::select(ADM0_CODE, ihme_lc_id, mean)

      # If plotting multiple runs, make sure there arent duplicated ADM0 codes
      if (multiple_runs == T) {
        plot_df <- plot_df %>%
          distinct() %>%
          data.table()
      }

      centroids <- merge(centroids, subset(plot_df, select = c("ADM0_CODE", "ihme_lc_id")))

      # Merge and fortify outside of ggplot to speed things up
      if (multiple_runs == F) ad0_shape_reg <- suppressWarnings(merge(ad0_shape_reg, plot_df))

      ad0_shape_reg@data$id <- rownames(ad0_shape_reg@data)
      ad0_shape_reg_df <- fortify(ad0_shape_reg, region = "id")
      ad0_shape_reg_df <- merge(ad0_shape_reg_df, ad0_shape_reg@data, by = "id")
      ad0_shape_reg_df <- as.data.table(ad0_shape_reg_df)

      # Merge and fortify outside of ggplot to speed things up
      ad0_shape_reg <- suppressWarnings(merge(ad0_shape_reg, plot_df))

      if (highisbad) distiller_direction <- -1
      if (!highisbad) distiller_direction <- 1

      # Prepare last year maps
      if (multiple_runs == F) gg_lastyear_map <- last_year_map(ad0_shape_reg_df, centroids, admin = 0, val_range, distiller_direction, max(year_list))
      if (multiple_runs == T) gg_lastyear_map <- last_year_map_multiple(ad0_shape_reg_df, centroids, admin = 0, max(year_list))

      # Correct the aspect ratio if not plotting data
      gg_lastyear_map <- aspect_ratio_plot(gg_lastyear_map, 4 / 3, 0.25)

      # Set up the final layout ---------------------------------------------------------------------------

      # Create a title grob
      title_grob <- title_generate(ind_title, year_list, title_grob_size, admin = 0, ad0_reg_title = ad0_reg_title)

      # Create the overall plot by arranging the grobs
      lay <- plot_overlay(plot_data, multiple_runs)
      if (plot_data == T & multiple_runs == F) master_plot <- arrangeGrob(gg_ad0_ts, gg_locator_map, title_grob, gg_lastyear_map, layout_matrix = lay)
      if (plot_data == T & multiple_runs == T) master_plot <- arrangeGrob(gg_ad0_ts, gg_locator_map, title_grob, gg_lastyear_map, ihme_grob_multiple, layout_matrix = lay)
      if (plot_data == F) {
        master_plot <- arrangeGrob(gg_ad0_ts, gg_locator_map, title_grob, gg_lastyear_map, ihme_grob,
          layout_matrix = lay, heights = c(1, 1, 1, 1, 1, 0.6)
        )
      }
      grid.draw(master_plot)
      if (i != length(ad0_code_list)) plot.new()
    } # END `for (i in 1:length(ad0_code_list))`
    dev.off() # For PDF
  } # END `if` wrapper for ad0

  ################################################################################
  # II. ADMIN-1-LEVEL PLOTS FOR ALL COUNTRIES ####################################
  ################################################################################

  if ("ad1" %in% plot_levels) {
    pdf_filename <- <<<< FILEPATH REDACTED >>>>

    if (verbose == T) {
      message("\n############################################################")
      message("############################################################")
      message(paste0("Plotting Admin-1-level estimates by country\n"))
      message(paste0("  Writing output to ", pdf_filename, "\n"))
    }

    pdf(
      file = pdf_filename,
      height = 10,
      width = 18
    )

    # Loop over countries
    # If specified, only map over listed countries
    if (!is.null(ad1_map_countries)) {
      message(paste0(
        "You have chosen to map only specified countries:\n",
        paste(ad1_map_countries, collapse = " ")
      ))
      countries <- ad1_map_countries
    }


    for (ctry in sort(countries)) {
      if (verbose == T) {
        message(paste0("  --> ", ctry, "..."))
      }

      # Get the ad0 code for this country
      ad0_code <- unique(ad0_df[ADM0_NAME == ctry]$ADM0_CODE)

      # Create subsets of ad1 and ad2 dfs just for this country (for convenience)
      ad0_df_ctry <- subset(ad0_df, ADM0_CODE == ad0_code)
      ad1_df_ctry <- subset(ad1_df, ADM0_CODE == ad0_code)
      ad2_df_ctry <- subset(ad2_df, ADM0_CODE == ad0_code)

      # Create a lookup table of ad1s
      ad1_table <- unique(subset(ad1_df_ctry, select = c("ADM1_NAME", "ADM1_CODE")))

      # Subset maps
      ad0_national <- subset(ad0_shape_simple, ADM0_CODE == ad0_code)
      ad1_national <- subset(ad1_shape_simple, ADM0_CODE == ad0_code)
      ad2_national <- subset(ad2_shape_simple, ADM0_CODE == ad0_code)

      # Wrap facet labels & make sure ADMIN 1 comes first
      n_ad1s <- length(unique(ad1_table$ADM1_CODE))

      if (n_ad1s > 0 & n_ad1s <= 16) wrap_width <- 18
      if (n_ad1s > 16 & n_ad1s <= 25) wrap_width <- 12
      if (n_ad1s > 25) wrap_width <- 9

      # Set up plot names; include national-level estimates; ensure ordering correct
      ad1_df_ctry[, plot_name := stringr::str_wrap(ADM1_NAME, width = wrap_width)]
      ad0_df_ctry[, plot_name := "NATIONAL"]
      ad1_df_plot <- rbind(ad1_df_ctry, ad0_df_ctry, fill = T)
      ad1_lvls <- c("NATIONAL", setdiff(unique(ad1_df_plot$plot_name), "NATIONAL"))
      ad1_df_plot$plot_name <- factor(ad1_df_plot$plot_name, levels = ad1_lvls)

      # Admin 1 time series plot ----------------------------------------------------

      if (multiple_runs == F) gg_ad1ts <- time_series(ad1_df_plot, admin = 1, val_range, title_plot_size, ind_title, plot_stkrs = plot_stacker_trends, stkrs = stackers)
      if (multiple_runs == T) gg_ad1ts <- time_series_multiple(ad1_df_plot, admin = 1, val_range, title_plot_size, ind_title)

      # A locator map --------------------------------------------------------------

      gg_locator_map <- location_map_draw(subset(ad0_shape_simple, ADM0_CODE == ad0_code), ad0_shape_simple)

      # A labeled map by admin1 ---------------------------------------------------

      centroids <- gCentroid(ad1_national, byid = T, id = ad1_national$ADM1_CODE)
      centroids <- as.data.frame(centroids) %>%
        cbind(rownames(.), .) %>%
        as.data.table(.) %>%
        setnames(., names(.), c("ADM1_CODE", "x", "y"))
      centroids$ADM1_CODE <- as.numeric(levels(centroids$ADM1_CODE))[centroids$ADM1_CODE]
      centroids <- merge(centroids, as.data.table(ad1_national), by = "ADM1_CODE")
      plot_df <-
        ad1_df_ctry %>%
        filter(year == max(year_list)) %>%
        dplyr::select(ADM1_CODE, mean)

      # Merge and fortify outside of ggplot to speed things up
      if (multiple_runs == F) ad1_national <- suppressWarnings(merge(ad1_national, plot_df))
      if (multiple_runs == T) {
        plot_df <- plot_df %>%
          distinct() %>%
          data.table()
      }

      ad1_national@data$id <- rownames(ad1_national@data)
      ad1_national_df <- fortify(ad1_national, region = "id")
      ad1_national_df <- merge(ad1_national_df, ad1_national@data, by = "id")
      ad1_national_df <- as.data.table(ad1_national_df)

      if (highisbad) distiller_direction <- -1
      if (!highisbad) distiller_direction <- 1

      # Prepare for plotting
      if (multiple_runs == F) gg_ad1_map <- last_year_map(ad1_national_df, centroids, admin = 1, val_range, distiller_direction, max(year_list))
      if (multiple_runs == T) gg_ad1_map <- last_year_map_multiple(ad1_national_df, centroids, admin = 1, max(year_list))

      # Correct the aspect ratio if not plotting data
      gg_ad1_map <- aspect_ratio_plot(gg_ad1_map, 4 / 3, 0.25)

      # Set up the final layout ----------------------------------------------------

      # Create a title grob
      title_grob <- title_generate(ind_title, year_list, title_grob_size, admin = 1, ctry = ctry)

      # Create the overall plot by arranging the grobs
      lay <- plot_overlay(plot_data, multiple_runs)
      if (plot_data == T & multiple_runs == F) master_plot <- arrangeGrob(gg_ad1ts, gg_locator_map, title_grob, gg_ad1_map, layout_matrix = lay)
      if (plot_data == T & multiple_runs == T) master_plot <- arrangeGrob(gg_ad1ts, gg_locator_map, title_grob, gg_ad1_map, ihme_grob_multiple, layout_matrix = lay)
      if (plot_data == F) {
        master_plot <- arrangeGrob(gg_ad1ts, gg_locator_map, title_grob, gg_ad1_map, ihme_grob,
          layout_matrix = lay, heights = c(1, 1, 1, 1, 1, 0.6)
        )
      }
      grid.draw(master_plot)

      if (ctry != sort(countries)[length(sort(countries))]) plot.new()
    } # END ` for (ctry in sort(countries))`

    dev.off() # For PDF
  } # END `if` wrapper for ad1

  ################################################################################
  # III. ADMIN-2-LEVEL PLOTS BY COUNTRY ##########################################
  ################################################################################

  if ("ad2" %in% plot_levels) {
    if (verbose == T) {
      message("\n############################################################")
      message("############################################################")
      message(paste0("Plotting Admin-2-level estimates by country\n"))
    }
    # If specified, only map over listed countries
    if (!is.null(ad1_map_countries)) {
      message(paste0(
        "You have chosen to map only specified countries:\n",
        paste(ad1_map_countries, collapse = " ")
      ))
      countries <- ad1_map_countries
    }

    # make graph color and shapes mapping uniform across all graphs
    if (ind_gp == "lf") {
      data_collect_method_vals <- sort(unique(ad2_data$data_collect_method))
    } else if (ind_gp == "oncho") {
      data_collect_method_vals <- sort(unique(ad2_data$source))
    }
    diagnostic_vals <- sort(unique(ad2_data$diagnostic))

    # Loop over countries and create a PDF for each country
    for (ctry in sort(countries)) {

      # Create subsets of ad1 and ad2 dfs just for this country (for convenience)
      ad0_df_ctry <- subset(ad0_df, ADM0_NAME == ctry)
      ad1_df_ctry <- subset(ad1_df, ADM0_NAME == ctry)
      ad2_df_ctry <- subset(ad2_df, ADM0_NAME == ctry)

      # ad0_df_ctry$ihme_lc_id.y <- NULL
      # setnames(ad0_df_ctry, "ihme_lc_id.x", "ihme_lc_id")
      
      # Create a lookup table of ad1s
      ad1_table <- unique(subset(ad1_df_ctry, select = c("ADM1_NAME", "ADM1_CODE")))

      # Get the ad0 code for this country
      ad0_code <- unique(ad0_df_ctry$ADM0_CODE)

      # Subset maps
      ad0_national <- subset(ad0_shape_simple, ADM0_CODE == ad0_code)
      ad1_national <- subset(ad1_shape_simple, ADM0_CODE == ad0_code)
      ad2_national <- subset(ad2_shape_simple, ADM0_CODE == ad0_code)

      if (plot_data == T) {
        ad0_data_ctry <- subset(ad0_data, ADM0_NAME == ctry)
        ad1_data_ctry <- subset(ad1_data, ADM0_NAME == ctry)
        ad2_data_ctry <- subset(ad2_data, ADM0_NAME == ctry)
      }

      # Set up PDF filename
      pdf_filename <- paste0(out_dir, sprintf(out_filename_format, paste0(unique(ad0_df_ctry$ihme_lc_id), "_by_admin_2")))

      pdf(
        file = pdf_filename,
        height = 10,
        width = 18
      )

      if (verbose == T) {
        message(paste0("  --> ", ctry, "..."))
        message(paste0("      Writing file to ", pdf_filename))
      }

      # Create a lookup table of ad1s
      ad1_table <- unique(subset(ad1_df_ctry, select = c("ADM1_NAME", "ADM1_CODE")))

      for (ad1_code in ad1_table$ADM1_CODE) {
        # Get ad2 and ad1 tables just for this admin1 code
        # For this loop, ad2_df_ad1 means that it's a data frame of ad2s within the ad1
        # As opposed to ad1_df_ctry (ad1s within the country), etc.

        ad2_df_ad1 <- subset(ad2_df_ctry, ADM1_CODE == ad1_code)
        ad1_df_ad1 <- subset(ad1_df_ctry, ADM1_CODE == ad1_code)
        ad1_df_ad1[, plot_name := "ADMIN 1"]

        if (plot_data == T) ad2_data_ad1 <- subset(ad2_data_ctry, ADM1_CODE == ad1_code)

        if (verbose == T) message(paste0("     --> ", unique(ad1_df_ad1$ADM1_NAME)))

        # Create a single df for plotting

        # Wrap facet labels & make sure ADMIN 1 comes first
        n_ad2s <- length(unique(ad2_df_ad1$ADM2_CODE))

        if (n_ad2s > 0 & n_ad2s <= 16) wrap_width <- 18
        if (n_ad2s > 16 & n_ad2s <= 25) wrap_width <- 12
        if (n_ad2s > 25) wrap_width <- 9

        ad2_df_ad1[, plot_name := stringr::str_wrap(ADM2_NAME, width = wrap_width)]
        plot_df <- copy(ad2_df_ad1)
        ad2_lvls <- c(unique(ad2_df_ad1$plot_name))
        plot_df$plot_name <- factor(plot_df$plot_name, levels = ad2_lvls)

        if (plot_data == T) {
          ad2_data_ad1[, plot_name := stringr::str_wrap(ADM2_NAME, width = wrap_width)]
          ad2_data_ad1$plot_name <- factor(ad2_data_ad1$plot_name, levels = ad2_lvls)
        }

        # subset MDA implementation data
        mda_line <- mda_data[ADM2_CODE %in% ad2_df_ad1$ADM2_CODE, ]
        mda_line[, plot_name := stringr::str_wrap(ADM2_NAME, width = wrap_width)]
        mda_line$plot_name <- factor(mda_line$plot_name, levels = ad2_lvls)

        if (sum(mda_line$mda_start, na.rm = TRUE) == 0) {
          mda_line <- NULL
        }

        # Admin 2 time-series plot -------------------------------------------------

        if (plot_data == F & multiple_runs == F) gg_ad2ts <- time_series(plot_df, admin = 2, val_range, title_plot_size, ind_title, mda = mda_line, plot_stkrs = plot_stacker_trends, stkrs = stackers)
        if (plot_data == T & multiple_runs == F) gg_ad2ts <- time_series_data(plot_df, ad2_data_ad1, admin = 2, val_range, title_plot_size, ind_title, mda = mda_line, collect_method_vals = data_collect_method_vals, diag_vals = diagnostic_vals, plot_stkrs = plot_stacker_trends, stkrs = stackers)
        if (plot_data == F & multiple_runs == T) gg_ad2ts <- time_series_multiple(plot_df, admin = 2, val_range, title_plot_size, ind_title, mda = mda_line)
        if (plot_data == T & multiple_runs == T) gg_ad2ts <- time_series_multiple_data(plot_df, ad2_data_ad1, admin = 2, val_range, title_plot_size, ind_title, mda = mda_line, collect_method_vals = data_collect_method_vals, diag_vals = diagnostic_vals)

        # A locator map ------------------------------------------------------------

        gg_locator_map <- ggplot() +
          geom_polygon_quiet(
            data = subset(ad1_national, ADM1_CODE == ad1_code),
            aes(x = long, y = lat, group = group),
            fill = "red"
          ) +
          geom_path_quiet(
            data = ad1_national,
            aes(x = long, y = lat, group = group),
            size = 0.2
          ) +
          geom_path_quiet(
            data = ad0_national,
            aes(x = long, y = lat, group = group),
            size = 0.5
          ) +
          coord_equal() +
          theme_empty


        # A labeled map of the admin1 unit by admin2 -----------------------------------------
        ad2_ad1 <- subset(ad2_national, ADM1_CODE == ad1_code)

        centroids <- gCentroid(ad2_ad1, byid = T, id = ad2_ad1$ADM2_CODE)
        centroids <- as.data.frame(centroids) %>%
          cbind(rownames(.), .) %>%
          as.data.table(.) %>%
          setnames(., names(.), c("ADM2_CODE", "x", "y"))
        centroids$ADM2_CODE <- as.numeric(levels(centroids$ADM2_CODE))[centroids$ADM2_CODE]
        centroids <- merge(centroids, as.data.table(ad2_ad1), by = "ADM2_CODE")

        plot_df <-
          ad2_df_ad1 %>%
          filter(year == max(year_list)) %>%
          dplyr::select(ADM2_CODE, mean)

        # Merge and fortify outside of ggplot to speed things up
        if (multiple_runs == T) {
          plot_df <- plot_df %>%
            distinct() %>%
            data.table()
        }
        if (multiple_runs == F) ad2_ad1 <- suppressWarnings(merge(ad2_ad1, plot_df, duplicateGeoms = T))

        ad2_ad1@data$id <- rownames(ad2_ad1@data)
        ad2_ad1_df <- fortify(ad2_ad1, region = "id")
        ad2_ad1_df <- merge(ad2_ad1_df, ad2_ad1@data, by = "id")
        ad2_ad1_df <- as.data.table(ad2_ad1_df)

        if (highisbad) distiller_direction <- -1
        if (!highisbad) distiller_direction <- 1

        # Prepare for plotting
        if (multiple_runs == F) gg_ad2_map <- last_year_map(ad2_ad1_df, centroids, admin = 2, val_range, distiller_direction, max(year_list))
        if (multiple_runs == T) gg_ad2_map <- last_year_map_multiple(ad2_ad1_df, centroids, admin = 2, max(year_list))

        if (plot_data == F) gg_ad2_map <- aspect_ratio_plot(gg_ad2_map, 4 / 3, 0.25)

        # Set up the final layout ------------------------------------------------------------------
        # Create a title grob
        title_grob <- title_generate(ind_title, year_list, title_grob_size, admin = 2, ad1_df_ad1 = ad1_df_ad1, ctry = ctry)

        # Create the overall plot by arranging the grobs
        lay <- plot_overlay(plot_data, multiple_runs)
        if (plot_data == T & multiple_runs == F) master_plot <- arrangeGrob(gg_ad2ts, gg_locator_map, title_grob, gg_ad2_map, layout_matrix = lay)
        if (plot_data == T & multiple_runs == T) master_plot <- arrangeGrob(gg_ad2ts, gg_locator_map, title_grob, gg_ad2_map, ihme_grob_multiple, layout_matrix = lay)
        if (plot_data == F) {
          master_plot <- arrangeGrob(gg_ad2ts, gg_locator_map, title_grob, gg_ad2_map, ihme_grob,
            layout_matrix = lay, heights = c(1, 1, 1, 1, 1, 0.6)
          )
        }

        grid.draw(master_plot)

        if (ad1_code != ad1_table$ADM1_CODE[length(ad1_table$ADM1_CODE)]) plot.new()
      } # END `for (ad1_code in ad1_table$ADM1_CODE)`

      dev.off() # For PDF (by country)
    } # END ` for (ctry in sort(countries))`
  } # END `if` wrapper for ad2
}
