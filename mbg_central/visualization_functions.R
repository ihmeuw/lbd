########################################################################
# Useful functions for visualizing MBG outputs #########################
########################################################################

## subnational_ts_plots ################################################

#' Function to create subnational time-series plots with accompanying maps
#'
#' This function is meant to create three sets of PDFs that show
#' time series plots of administrative-level estimates of your
#' modeled quantity of interest:
#'
#' 1) Plots of country-level estimates (single PDF), optionally one page per region
#' 2) Plots of admin1-level estimates (single PDF), one page per country
#' 3) Plots of admin2-level estimates (one PDF per country), one page per admin1
#'
#' As currently written this function is quite repetitive. This allows for each of the
#' maps to be tweaked more easily, but in the future it would be better to
#' encapsulate some of the repeated bits in flexible functions that can be
#' called for each of the categories of maps created above. This would
#' also potentially allow for the use of non-standard shapefiles, which isn't
#' well supported by the function at the moement.
#'
#' All data tables referenced below (and accompanying shapefiles) must have the
#' standard ADM0_CODE, ADM0_NAME, etc. naming conventions and have those codes
#' and names included as columns / shapefile attributes. In addition,
#' the administrative estimates must have "mean", "upper", and "lower" columns
#' to allow uncertainty intervals to be plotted.
#'
#' @param ad0_df admin-0 level data.frame or data.table from the standard aggregation code
#' @param ad1_df admin-1 level data.frame or data.table from the standard aggregation code
#' @param ad2_df admin-2 level data.frame or data.table from the standard aggregation code.
#' @param ad0_shape admin0 SPDF object that corresponds to `ad0_df` above.
#'                  If none provided, loads from GAUL as default
#' @param ad1_shape admin1 SPDF object that corresponds to `ad1_df` above.
#'                  If none provided, loads from GAUL as default
#' @param ad2_shape admin2 SPDF object that corresponds to `ad2_df` above.
#'                  If none provided, loads from GAUL as default
#' @param ind_title title for your indicator, e.g. "DPT3 Coverage".
#'                  For use in plots in multiple places.  String.
#' @param out_dir output directory. String.
#' @param out_filename_format output filename format, in `sprintf()` format
#'                            with `%s` representing the place where the
#'                            geographic identifier will go. String.
#' @param val_range range of values for your indicator. Default: [0,1]
#'                  Vector of length 2 (`c(min, max)`)
#' @param highisbad if T, higher values will be more red. Boolean.
#' @param ad0_map_regions if you want to have the national-level time series
#'                        broken down by modeling region, you can pass a vector
#'                        of region names here.  These will be processed by
#'                        `get_gaul_codes()` so need to be identified in that
#'                        function.
#' @param ad0_map_region_titles corresponding to ad0_map_regions, you can pass
#'                              titles for each of the map regions here (so that
#'                              `essa` is identified as "Eastern Sub-Saharan Africa",
#'                              for instance).  Vector of same length of ad0_map_regions
#' @param ad1_map_countries If specified, which countries do you want to plot at the admin1/admin2 level. 
#'                          If null, the default is plotting all countries in the specified ad0_map_regions
#' @param plot_levels which administrative levels do you want to plot? Vector containing
#'                    any combination of "ad0", ad1", and/or "ad2".
#' @param title_plot_size if you want, you can specify a font size to be used for ind_title in text_element.
#'                        If none provided, default value is assigned
#' @param title_grob_size if you want, you can specify a font size to be used for ind_title in text_grob.
#'                        If none provided, default value is assigned
#' @param multiple_runs should the function plot multiple model runs/ multiple indicators on the same plot? Boolean
#' @param plot_data should the function plot the collapsed admin input data? Boolean
#' @param ad0_data Collapsed admin0 input data, obtained using input_aggregate_admin function found below
#' @param ad1_data Collapsed admin1 input data, obtained using input_aggregate_admin function found below
#' @param ad2_data Collapsed admin2 input data, obtained using input_aggregate_admin function found below
#' @param verbose should the function update progress along the way? Boolean.
#' @param shapefile_version specifies the shapefile version to be used. `current` or date `yyyy_mm_dd`
#' @param tol tolerance to be used when simplifying polygons. default `0.001` creates orphaned holes in GADM shapefile, can use `0.0001` to fix but the function slows down and outputs are larger.
#'
#' @return
#' Saves a series of PDFs as described above in `out_dir`
#'
#' @examples
#'
#'
#' # Set up directories and files  ########################################
#'
#' # Set run_date, indicator, indicator_group, out_dir per your preferences
#'
#' share_dir <- '<<<< FILEPATH REDACTED >>>>'
#' in_dir <- '<<<< FILEPATH REDACTED >>>>'
#'
#' in_file_ad0 <- paste0(in_dir, indicator, "_admin_0_raked_summary.csv")
#' in_file_ad1 <- paste0(in_dir, indicator, "_admin_1_raked_summary.csv")
#' in_file_ad2 <- paste0(in_dir, indicator, "_admin_2_raked_summary.csv")
#'
#' # Prepare inputs #######################################################
#'
#' ad0_df <- fread(in_file_ad0)
#' ad1_df <- fread(in_file_ad1)
#' ad2_df <- fread(in_file_ad2)
#'
#' # Drop Ma'tan al-Sarra if present
#' ad0_df <- subset(ad0_df, ADM0_CODE != 40762)
#' ad1_df <- subset(ad1_df, ADM0_CODE != 40762)
#' ad2_df <- subset(ad2_df, ADM0_CODE != 40762)
#'
#' # Load GAUL shapefiles #################################################
#'
#' # Pre-loading them here and using the "africa" background map instead of default GAUL
#' # since it's a bit prettier - cuts off far-outlying islands, etc.)
#'
#' ad0_shape <- readRDS('<<<< FILEPATH REDACTED >>>>')
#' ad1_shape <- readOGR(get_admin_shapefile(admin_level=1))
#' ad2_shape <- readOGR(get_admin_shapefile(admin_level=2))
#'
#' # Run the plotting code ################################################
#'
#' subnational_ts_plots(ad0_df = ad0_df,
#'                      ad1_df = ad1_df,
#'                      ad2_df = ad2_df,
#'                      ad0_shape = ad0_shape,
#'                      ad1_shape = ad1_shape,
#'                      ad2_shape = ad2_shape,
#'                      ind_title = "My Indicator",
#'                      out_dir = out_dir,
#'                      highisbad = F,
#'                      val_range = c(0,1),
#'                      ad0_map_regions = c("cssa", "essa", "name", "sssa", "wssa"),
#'                      ad0_map_region_titles = c("Central Sub-Saharan Africa",
#'                                                "Eastern Sub-Saharan Africa",
#'                                                "Northern Africa",
#'                                                "Southern Sub-Saharan Africa",
#'                                                "Western Sub-Saharan Africa"),
#'                      verbose = T)
#'
#' ## Running plotting code for multiple model runs or multiple indicators###################################
#'
#'  Set up directories and files  ########################################
#'
#' # Set run_dates (for multiple model runs), indicator (or indicators), indicator_group, out_dir per your preferences
#'
#' run_dates <- c("run_date_1", "run_date_2,... "run_date_x)
#' run_label <- c("indicator 1", "indicator 2",... "indicator x")
#' indicator <- c("indicator 1", "indicator 2", ... "indicator x")
#' share_dir <- '<<<< FILEPATH REDACTED >>>>'
#' in_dir <- '<<<< FILEPATH REDACTED >>>>'
#'
#' in_file_ad0 <- paste0(in_dir, indicator, "_admin_0_raked_summary.csv")
#' in_file_ad1 <- paste0(in_dir, indicator, "_admin_1_raked_summary.csv")
#' in_file_ad2 <- paste0(in_dir, indicator, "_admin_2_raked_summary.csv")
#'
#' # Prepare inputs #######################################################
#'
#' Read in all mode admin aggregations, adding run label to each with add_run_label function
#' ad0_df <- lapply(in_file_ad0, fread) %>% add_run_label(run_label)
#' ad1_df <- lapply(in_file_ad1, fread) %>% add_run_label(run_label)
#' ad2_df <- lapply(in_file_ad2, fread) %>% add_run_label(run_label)
#'
#' # Drop Ma'tan al-Sarra if present
#' ad0_df <- subset(ad0_df, ADM0_CODE != 40762)
#' ad1_df <- subset(ad1_df, ADM0_CODE != 40762)
#' ad2_df <- subset(ad2_df, ADM0_CODE != 40762)
#'
#' Only load shapefiles if you do not intend to use the GAUL shapefiles; the function has been written
#' so that if you are using GAUL you do not have to simplify the shape each time. This will make the function much quicker
#'
#' # Run the plotting code ################################################
#'
#' subnational_ts_plots(ad0_df = ad0_df,
#'                      ad1_df = ad1_df,
#'                      ad2_df = ad2_df,
#'                      ind_title = "My Indicator",
#'                      out_dir = out_dir,
#'                      highisbad = F,
#'                      val_range = c(0,1),
#'                      ad0_map_regions = c("cssa", "essa", "name", "sssa", "wssa"),
#'                      ad0_map_region_titles = c("Central Sub-Saharan Africa",
#'                                                "Eastern Sub-Saharan Africa",
#'                                                "Northern Africa",
#'                                                "Southern Sub-Saharan Africa",
#'                                                "Western Sub-Saharan Africa"),
#'                      multiple_runs = T,
#'                      verbose = T)
#'
#'
#' # to see example of running plots with data, see input_aggregate_admin function. See add_run_label function for questions containing how to specify "run" column
#'

subnational_ts_plots <- function(ad0_df,
                                 ad1_df,
                                 ad2_df,
                                 ad0_shape = NULL,
                                 ad1_shape = NULL,
                                 ad2_shape = NULL,
                                 ind_title,
                                 out_dir,
                                 out_filename_format = "subnational_ts_plots_%s.pdf",
                                 val_range = c(0,1),
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
                                 shapefile_version = 'current',
                                 tol = 0.001) {
  
  ################################################################################
  # 0. SETUP #####################################################################
  # Load packages, define themes, and set up objects needed for plotting #########
  ################################################################################
  
  # Load libraries and functions -------------------------------------------------
  
  library(ggrepel)
  library(gridExtra)
  library(ggplot2)
  library(stringr)
  
  # Check to see if plot_data = T that admin data is attached
  if (plot_data == T) {
    if ("ad0" %in% plot_levels & is.null(ad0_data)) stop(paste("To plot data, you must include aggregated ad0_data.",
                                                               "See input_aggregate_admin function"))
    if ("ad1" %in% plot_levels & is.null(ad1_data)) stop(paste("To plot data, you must include aggregated ad1_data.",
                                                               "See input_aggregate_admin function"))
    if ("ad2" %in% plot_levels & is.null(ad2_data)) stop(paste("To plot data, you must include aggregated ad2_data.",
                                                               "See input_aggregate_admin function"))
  }
  
  if (plot_data == T & multiple_runs == T) message(paste("You are choosing to plot multiple indicators/model runs with data.\n",
                                                         "This is usually done if plotting multiple model runs from the same indicator"))
  
  # If plotting multiple model runs, make sure there is a column title 'run'
  if (multiple_runs == T){
    if (!"run" %in% names(ad0_df)) stop("Must have a column title 'run' that differentiates between model runs if plotting multiple indicators/model runs. See add_run_label function")
    if (!"run" %in% names(ad1_df)) stop("Must have a column title 'run' that differentiates between model runs if plotting multiple indicators/model runs. See add_run_label function")
    if (!"run" %in% names(ad2_df)) stop("Must have a column title 'run' that differentiates between model runs if plotting multiple indicators/model runs. See add_run_label function")
  }
  # Create directories -----------------------------------------------------------
  
  dir.create(out_dir, showWarnings = F)
  
  # Prepare data and shapefiles --------------------------------------------------
  
  subset_codes <- get_adm0_codes(ad0_map_regions, shapefile_version = shapefile_version)
  
  # Wrapper for gSimplify() to allow it to work with SPDFs
  simplify_spdf <- function(the_spdf, ...) {
    simple_spdf <- gSimplify(the_spdf, topologyPreserve = T, ...)
    return(SpatialPolygonsDataFrame(simple_spdf, the_spdf@data))
  }
  
  # Use simplified GAUL shapefiles as default if none provided, if provided simplify the shapefile for speed
  if (is.null(ad0_shape)){
    ad0_shape_simple <- readOGR(get_admin_shapefile(admin_level = 0, raking = F, version = shapefile_version)) %>% subset(ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  } else {
    ad0_shape_simple <- subset(ad0_shape, ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  }
  
  if (is.null(ad1_shape)){
    ad1_shape_simple <- readOGR(get_admin_shapefile(admin_level = 1, raking = F, version = shapefile_version)) %>% subset(ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  } else {
    ad1_shape_simple <- subset(ad1_shape, ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  }
  
  if (is.null(ad2_shape)){
    ad2_shape_simple <- readOGR(get_admin_shapefile(admin_level = 2, raking = F, version = shapefile_version)) %>% subset(ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  } else {
    ad2_shape_simple <- subset(ad2_shape, ADM0_CODE %in% subset_codes) %>% simplify_spdf(tol = tol)
  }
  
  # Merge on ihme_lc_ids
  gaul_to_loc_id <- get_location_code_mapping(shapefile_version = shapefile_version)
  gaul_to_loc_id <- subset(gaul_to_loc_id, select = c("GAUL_CODE", "ihme_lc_id"))
  setnames(gaul_to_loc_id,"GAUL_CODE", "ADM0_CODE")
  
  #renaming the ADMX_NAME columns using the names in ad2_shape, because this functions matches
  #everything on names and things get dropped if the names don't line up perfectly
  #This is an alternative to rewriting the entire function to match on coes
  admin_shp_data_adm0 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm0 <- unique(admin_shp_data_adm0[,c("ADM0_CODE", "ADM0_NAME")])
  admin_shp_data_adm0$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm0$ADM0_CODE))
  admin_shp_data_adm0$ADM0_NAME <- as.character(admin_shp_data_adm0$ADM0_NAME)
  gaul_to_loc_id_adm0 <- merge(gaul_to_loc_id, admin_shp_data_adm0, by = "ADM0_CODE")
  
  admin_shp_data_adm1 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm1 <- unique(admin_shp_data_adm1[,c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME")])
  admin_shp_data_adm1$ADM1_CODE <- as.integer(as.character(admin_shp_data_adm1$ADM1_CODE))
  admin_shp_data_adm1$ADM0_CODE <- as.integer(as.character(admin_shp_data_adm1$ADM0_CODE))
  admin_shp_data_adm1$ADM1_NAME <- as.character(admin_shp_data_adm1$ADM1_NAME)
  admin_shp_data_adm1$ADM0_NAME <- as.character(admin_shp_data_adm1$ADM0_NAME)
  gaul_to_loc_id_adm1 <- merge(gaul_to_loc_id, admin_shp_data_adm1, by = "ADM0_CODE")
  
  admin_shp_data_adm2 <- as.data.table(ad2_shape_simple)
  admin_shp_data_adm2 <- unique(admin_shp_data_adm2[,c("ADM0_CODE", "ADM0_NAME", "ADM1_CODE", "ADM1_NAME", "ADM2_CODE", "ADM2_NAME")])
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
  
  ad0_df <- merge(copy(as.data.table(ad0_df)), gaul_to_loc_id_adm0, all.x=T, all.y=F, by = 'ADM0_CODE')
  ad1_df <- merge(copy(as.data.table(ad1_df)), gaul_to_loc_id_adm1, all.x=T, all.y=F, by = 'ADM1_CODE')
  ad2_df <- merge(copy(as.data.table(ad2_df)), gaul_to_loc_id_adm2, all.x=T, all.y=F, by = 'ADM2_CODE')
  
  if (plot_data == T){
    ad0_data <- merge(copy(as.data.table(ad0_data)), gaul_to_loc_id, all.x=T, all.y=F, by = 'ADM0_CODE')
    ad1_data <- merge(copy(as.data.table(ad1_data)), gaul_to_loc_id, all.x=T, all.y=F, by = 'ADM0_CODE')
    ad2_data <- merge(copy(as.data.table(ad2_data)), gaul_to_loc_id, all.x=T, all.y=F, by = 'ADM0_CODE')
  }
  
  # Set up default values for title font size if not specified
  if (is.null(title_plot_size)) title_plot_size <- 20
  if (is.null(title_grob_size)) title_grob_size <- 30
  
  if (verbose == T) message("Simplifying shapes")
  
  # Subset and simplify shapefiles for memory & speed
  ad0_codes <- unique(ad0_df$ADM0_CODE)
  # If any region of Africa is being used, use full map of Africa so it does not look disjointed
  if (any(ad0_map_regions %in% c("cssa", "wssa", "essa", "sssa", "name"))) ad0_codes <- union(ad0_codes, get_adm0_codes("Africa", shapefile_version = shapefile_version))
  subset_codes <- ad0_codes
  
  # Make sure that the admin codes are converted to numeric if they are factors
  if (is.factor(ad0_shape_simple@data$ADM0_CODE)) ad0_shape_simple@data$ADM0_CODE <- as.numeric(levels(ad0_shape_simple@data$ADM0_CODE))[ad0_shape_simple@data$ADM0_CODE]
  if (is.factor(ad1_shape_simple@data$ADM1_CODE)) ad1_shape_simple@data$ADM1_CODE <- as.numeric(levels(ad1_shape_simple@data$ADM1_CODE))[ad1_shape_simple@data$ADM1_CODE]
  if (is.factor(ad2_shape_simple@data$ADM2_CODE)) ad2_shape_simple@data$ADM2_CODE <- as.numeric(levels(ad2_shape_simple@data$ADM2_CODE))[ad2_shape_simple@data$ADM2_CODE]
  
  # Add simple world shapefile if region is specified as Africa
  if ("ad0" %in% plot_levels & "africa" %in% tolower(ad0_map_regions)) world_shape_simple <- readRDS('<<<< FILEPATH REDACTED >>>>/world_shape_simple.rds')

  # Define some utility functions -----------------------------------------------------------------------------------------------------
  # Function to fix aspect ratios
  aspect_ratio_plot <- function(plt, aspect_ratio = 1, expand=0.05){
    
    xlims <- ggplot_build(plt)$layout$panel_scales$x[[1]]$range$range
    ylims <- ggplot_build(plt)$layout$panel_scales$y[[1]]$range$range
    
    # Try a different way if null - varies by ggplot build
    if (is.null(xlims)) xlims <- ggplot_build(plt)$layout$panel_scales_x[[1]]$range$range
    if (is.null(ylims)) ylims <- ggplot_build(plt)$layout$panel_scales_y[[1]]$range$range
    
    xmid <- mean(xlims); xrange <- xlims[2]-xlims[1]
    ymid <- mean(ylims); yrange <- ylims[2]-ylims[1]
    
    if (xrange/yrange < aspect_ratio) {
      new_yrange <- yrange*(1+expand)
      new_xrange <- new_yrange*aspect_ratio
    }
    
    if (xrange/yrange > aspect_ratio) {
      new_xrange = xrange*(1+expand)
      new_yrange <- new_xrange/aspect_ratio
    }
    
    new_xlim <- c(xmid - 0.5*new_xrange, xmid + 0.5*new_xrange)
    new_ylim <- c(ymid - 0.5*new_yrange, ymid + 0.5*new_yrange)
    
    return(plt+expand_limits(x=new_xlim,
                             y=new_ylim))
  }
  
  # Title generate makes the title Grob that goes in the right hand corner of the pdf
  title_generate <- function(ind_title,
                             year_list,
                             admin = 0,
                             ad0_reg_title = "",
                             ctry = "",
                             ad1_df_ad1 = ""){
    arrangeGrob(textGrob("", gp = gpar(fontsize = 30)),
                textGrob(str_wrap(ind_title, 18),
                         gp = gpar(fontsize = title_grob_size, fontface = "bold")),
                textGrob("", gp = gpar(fontsize = 10)),
                textGrob(if (admin == 0) ifelse(ad0_reg_title == "All countries", NULL, ad0_reg_title) else ctry,
                         gp = gpar(fontsize = 20)),
                textGrob("", gp = gpar(fontsize = 10)),
                textGrob(if (admin == 0) "National estimates" else if (admin == 1) "By First-level Administrative Unit" else paste0("Admin 1: ", unique(ad1_df_ad1$ADM1_NAME)),
                         gp = gpar(fontsize = 15)),
                textGrob("", gp = gpar(fontsize = 10)),
                textGrob(paste0(min(year_list), " - ", max(year_list)),
                         gp = gpar(fontsize = 18)),
                ncol = 1,
                heights = c(30, 25, 10, 25, 10, 15, 10, 20))
  }
  
  # Plot overlay defines where each object (ggplot and title) are placed on the pdf.
  # This depends on if data is being plotted, and if multiple runs are used
  plot_overlay <- function(plot_data, multiple_runs){
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
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  # Add IHME logo for plotting later
  ihme_logo <- png::readPNG('<<<< FILEPATH REDACTED >>>>/ihme_logo.png')
  ihme_grob <- rasterGrob(ihme_logo)

  # Add slightly different version for plotting multiple runs
  ihme_grob_multiple <- rasterGrob(ihme_logo, y = 0.2)
  
  # Only list countries that are within the specified regions
  if (!is.null(ad0_map_regions)) {
    ad0_code_list <- lapply(ad0_map_regions, get_adm0_codes, shapefile_version = shapefile_version)
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
    
    pdf(file = pdf_filename,
        height = 10,
        width = 18)
    
    for (i in 1:length(ad0_code_list)) {
      
      # Set up names & titles for this region ---------------------------------------
      ad0_reg_codes <- ad0_code_list[[i]]
      ad0_reg_title <- names(ad0_code_list)[i]
      
      ad0_df_reg <- subset(ad0_df, ADM0_CODE %in% ad0_reg_codes)
      
      # Make sure you are only modeling the countries in specified region, this will cause problems if doing map of Africa
      if (tolower(ad0_map_regions[i]) != "africa") ad0_df_reg <- subset(ad0_df_reg, region == ad0_map_regions[i])
      
      if (plot_data == T) ad0_data_reg <- subset(ad0_data, ADM0_CODE %in% ad0_reg_codes)
      
      if (verbose == T) message(paste0("  ", ad0_reg_title))
      
      # First, time series plot of all admin0s in the region
      n_ad0s <- length(unique(ad0_df_reg$ADM0_CODE))
      
      # Try to wrap names of countries as needed
      if (n_ad0s >  0 & n_ad0s <= 16) wrap_width <- 18
      if (n_ad0s > 16 & n_ad0s <= 25) wrap_width <- 12
      if (n_ad0s > 25)                wrap_width <- 9
      
      # Set up plot names; include national-level estimates; ensure ordering correct.
      # Use location ID as label modeling over more than 25 countries
      if (n_ad0s < 25) {
        ad0_df_reg[, plot_name := stringr::str_wrap(ADM0_NAME, width = wrap_width)]
      } else {
        ad0_df_reg[, plot_name := stringr::str_wrap(ihme_lc_id, width = wrap_width)]
      }
      
      # When plotting data, aligning plot names depending on number of admin0 in region
      if (plot_data == T){
        if (n_ad0s < 25) {
          ad0_data_reg[, plot_name := stringr::str_wrap(ADM0_NAME, width = wrap_width)]
        } else {
          ad0_data_reg[, plot_name := stringr::str_wrap(ihme_lc_id, width = wrap_width)]
        }
      }
      
      # Admin 0 time series plot --------------------------------------------------------------
      if (plot_data == F & multiple_runs == F) gg_ad0_ts <- time_series(ad0_df_reg, admin = 0, val_range, title_plot_size, ind_title)
      if (plot_data == T & multiple_runs == F) gg_ad0_ts <- time_series_data(ad0_df_reg, ad0_data_reg, admin = 0, val_range, title_plot_size, ind_title)
      if (plot_data == F & multiple_runs == T) gg_ad0_ts <- time_series_multiple(ad0_df_reg, admin = 0, val_range, title_plot_size, ind_title)
      if (plot_data == T & multiple_runs == T) gg_ad0_ts <- time_series_multiple_data(ad0_df_reg, ad0_data_reg, admin = 0, val_range, title_plot_size, ind_title)
      
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
        ad0_df_reg %>% filter(year == max(year_list)) %>% dplyr::select(ADM0_CODE, ihme_lc_id, mean)
      
      # If plotting multiple runs, make sure there arent duplicated ADM0 codes
      if (multiple_runs == T) plot_df <- plot_df %>% distinct() %>% data.table()
      
      centroids <- merge(centroids, subset(plot_df, select = c("ADM0_CODE", "ihme_lc_id")))
      
      # Merge and fortify outside of ggplot to speed things up
      if (multiple_runs == F) ad0_shape_reg <- suppressWarnings(merge(ad0_shape_reg, plot_df))
      
      ad0_shape_reg@data$id = rownames(ad0_shape_reg@data)
      ad0_shape_reg_df = fortify(ad0_shape_reg, region="id")
      ad0_shape_reg_df <- merge(ad0_shape_reg_df, ad0_shape_reg@data, by = "id")
      ad0_shape_reg_df <- as.data.table(ad0_shape_reg_df)
      
      if(highisbad) distiller_direction <- -1
      if(!highisbad) distiller_direction <- 1
      
      # Prepare last year maps
      if (multiple_runs == F) gg_lastyear_map <- last_year_map(ad0_shape_reg_df, centroids, admin = 0, val_range, distiller_direction, max(year_list))
      if (multiple_runs == T) gg_lastyear_map <- last_year_map_multiple(ad0_shape_reg_df, centroids, admin = 0, max(year_list))
      
      # Correct the aspect ratio if not plotting data
      if (plot_data == F) gg_lastyear_map <- aspect_ratio_plot(gg_lastyear_map, 4/3, 0.25)
      
      # Set up the final layout ---------------------------------------------------------------------------
      
      # Create a title grob
      title_grob <- title_generate(ind_title, year_list, title_grob_size, admin = 0, ad0_reg_title = ad0_reg_title)
      
      # Create the overall plot by arranging the grobs
      lay <- plot_overlay(plot_data, multiple_runs)
      if (plot_data == T & multiple_runs == F) master_plot <- arrangeGrob(gg_ad0_ts, gg_locator_map, title_grob, gg_lastyear_map, layout_matrix = lay)
      if (plot_data == T & multiple_runs == T) master_plot <- arrangeGrob(gg_ad0_ts, gg_locator_map, title_grob, gg_lastyear_map, ihme_grob_multiple, layout_matrix = lay)
      if (plot_data == F) master_plot <- arrangeGrob(gg_ad0_ts, gg_locator_map, title_grob, gg_lastyear_map, ihme_grob,
                                                     layout_matrix = lay, heights = c(1,1,1,1,1,0.6))
      grid.draw(master_plot)
      if (i != length(ad0_code_list)) plot.new()
      
      
    } # END `for (i in 1:length(ad0_code_list))`
    dev.off() # For PDF
  } # END `if` wrapper for ad0
  
  ################################################################################
  # II. ADMIN-1-LEVEL PLOTS FOR ALL COUNTRIES ####################################
  ################################################################################
  
  if ("ad1" %in% plot_levels) {
    
    pdf_filename <- paste0(out_dir, sprintf(out_filename_format, "all_countries_by_admin1"))
    
    if (verbose == T) {
      message("\n############################################################")
      message("############################################################")
      message(paste0("Plotting Admin-1-level estimates by country\n"))
      message(paste0("  Writing output to ", pdf_filename, "\n"))
      
    }
    
    pdf(file = pdf_filename,
        height = 10,
        width = 18)
    
    # Loop over countries
    # If specified, only map over listed countries
    if (!is.null(ad1_map_countries)){
      message(paste0("You have chosen to map only specified countries:\n",
                     paste(ad1_map_countries, collapse = " ")))
      countries <- ad1_map_countries
    }
    
    
    for (ctry in sort(countries)) {
      
      if(verbose == T) {
        message(paste0("  --> ", ctry, "..."))
      }
      
      # Get the ad0 code for this country
      ad0_code <- unique(ad0_df[ADM0_NAME == ctry]$ADM0_CODE)
      
      # Create subsets of ad1 and ad2 dfs just for this country (for convenience)
      ad0_df_ctry <- subset(ad0_df, ADM0_CODE == ad0_code)
      ad1_df_ctry <- subset(ad1_df, ADM0_CODE == ad0_code)
      ad2_df_ctry <- subset(ad2_df, ADM0_CODE == ad0_code)
      
      
      if (plot_data == T) {
        ad0_data_ctry <- subset(ad0_data, ADM0_CODE == ad0_code)
        ad1_data_ctry <- subset(ad1_data, ADM0_CODE == ad0_code)
        ad2_data_ctry <- subset(ad2_data, ADM0_CODE == ad0_code)
      }
      
      # Create a lookup table of ad1s
      ad1_table <- unique(subset(ad1_df_ctry, select = c("ADM1_NAME", "ADM1_CODE")))
      
      # Subset maps
      ad0_national <- subset(ad0_shape_simple, ADM0_CODE == ad0_code)
      ad1_national <- subset(ad1_shape_simple, ADM0_CODE == ad0_code)
      ad2_national <- subset(ad2_shape_simple, ADM0_CODE == ad0_code)
      
      # Wrap facet labels & make sure ADMIN 1 comes first
      n_ad1s <- length(unique(ad1_table$ADM1_CODE))
      
      if (n_ad1s >  0 & n_ad1s <= 16) wrap_width <- 18
      if (n_ad1s > 16 & n_ad1s <= 25) wrap_width <- 12
      if (n_ad1s > 25)                wrap_width <- 9
      
      # Set up plot names; include national-level estimates; ensure ordering correct
      ad1_df_ctry[, plot_name := stringr::str_wrap(ADM1_NAME, width = wrap_width)]
      ad0_df_ctry[, plot_name := "NATIONAL"]
      ad1_df_plot <- rbind(ad1_df_ctry, ad0_df_ctry, fill=T)
      ad1_lvls <- c("NATIONAL", setdiff(unique(ad1_df_plot$plot_name), "NATIONAL"))
      ad1_df_plot$plot_name <- factor(ad1_df_plot$plot_name, levels = ad1_lvls)
      
      if (plot_data == T) ad1_data_ctry[, plot_name := as.factor(stringr::str_wrap(ADM1_NAME, width = wrap_width))]
      
      # Admin 1 time series plot ----------------------------------------------------

      if (plot_data == F & multiple_runs == F) gg_ad1ts <- time_series(ad1_df_plot, admin = 1, val_range, title_plot_size, ind_title)
      if (plot_data == T & multiple_runs == F) gg_ad1ts <- time_series_data(ad1_df_plot, ad1_data_ctry, admin = 1, val_range, title_plot_size, ind_title)
      if (plot_data == F & multiple_runs == T) gg_ad1ts <- time_series_multiple(ad1_df_plot, admin = 1, val_range, title_plot_size, ind_title)
      if (plot_data == T & multiple_runs == T) gg_ad1ts <- time_series_multiple_data(ad1_df_plot, ad1_data_ctry, admin = 1, val_range, title_plot_size, ind_title)

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
        ad1_df_ctry %>% filter(year == max(year_list)) %>% dplyr::select(ADM1_CODE, mean)
      
      # Merge and fortify outside of ggplot to speed things up
      if (multiple_runs == F) ad1_national <- suppressWarnings(merge(ad1_national, plot_df))
      if (multiple_runs == T) plot_df <- plot_df %>% distinct() %>% data.table()
      
      ad1_national@data$id = rownames(ad1_national@data)
      ad1_national_df = fortify(ad1_national, region="id")
      ad1_national_df <- merge(ad1_national_df, ad1_national@data, by = "id")
      ad1_national_df <- as.data.table(ad1_national_df)
      
      if(highisbad) distiller_direction <- -1
      if(!highisbad) distiller_direction <- 1
      
      # Prepare for plotting
      if (multiple_runs == F) gg_ad1_map <- last_year_map(ad1_national_df, centroids, admin = 1, val_range, distiller_direction, max(year_list))
      if (multiple_runs == T) gg_ad1_map <- last_year_map_multiple(ad1_national_df, centroids, admin = 1, max(year_list))
      
      # Correct the aspect ratio if not plotting data
      if (plot_data == F) gg_ad1_map <- aspect_ratio_plot(gg_ad1_map, 4/3, 0.25)
      
      # Set up the final layout ----------------------------------------------------
      
      # Create a title grob
      title_grob <- title_generate(ind_title, year_list, title_grob_size, admin = 1, ctry = ctry)
      
      # Create the overall plot by arranging the grobs
      lay <- plot_overlay(plot_data, multiple_runs)
      if (plot_data == T & multiple_runs == F) master_plot <- arrangeGrob(gg_ad1ts, gg_locator_map, title_grob, gg_ad1_map, layout_matrix = lay)
      if (plot_data == T & multiple_runs == T) master_plot <- arrangeGrob(gg_ad1ts, gg_locator_map, title_grob, gg_ad1_map, ihme_grob_multiple, layout_matrix = lay)
      if (plot_data == F) master_plot <- arrangeGrob(gg_ad1ts, gg_locator_map, title_grob, gg_ad1_map, ihme_grob,
                                                     layout_matrix = lay, heights = c(1,1,1,1,1,0.6))
      grid.draw(master_plot)
      
      if (ctry != sort(countries)[length(sort(countries))]) plot.new()
      
    } # END ` for (ctry in sort(countries))`
    
    dev.off() # For PDF
    
  } #END `if` wrapper for ad1
  
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
    if (!is.null(ad1_map_countries)){
      message(paste0("You have chosen to map only specified countries:\n",
                     paste(ad1_map_countries, collapse = " ")))
      countries <- ad1_map_countries
    }
    
    # Loop over countries and create a PDF for each country
    for (ctry in sort(countries)) {
      
      # Create subsets of ad1 and ad2 dfs just for this country (for convenience)
      ad0_df_ctry <- subset(ad0_df, ADM0_NAME == ctry)
      ad1_df_ctry <- subset(ad1_df, ADM0_NAME == ctry)
      ad2_df_ctry <- subset(ad2_df, ADM0_NAME == ctry)
      
      # Create a lookup table of ad1s
      ad1_table <- unique(subset(ad1_df_ctry, select = c("ADM1_NAME", "ADM1_CODE")))
      
      # Get the ad0 code for this country
      ad0_code <- unique(ad0_df_ctry$ADM0_CODE)
      
      # Subset maps
      ad0_national <- subset(ad0_shape_simple, ADM0_CODE == ad0_code)
      ad1_national <- subset(ad1_shape_simple, ADM0_CODE == ad0_code)
      ad2_national <- subset(ad2_shape_simple, ADM0_CODE == ad0_code)
      
      if(plot_data == T) {
        ad0_data_ctry <- subset(ad0_data, ADM0_NAME == ctry)
        ad1_data_ctry <- subset(ad1_data, ADM0_NAME == ctry)
        ad2_data_ctry <- subset(ad2_data, ADM0_NAME == ctry)
      }
      
      # Set up PDF filename
      pdf_filename <- paste0(out_dir, sprintf(out_filename_format, paste0(unique(ad0_df_ctry$ihme_lc_id), "_by_admin_2")))

      pdf(file = pdf_filename,
          height = 10,
          width = 18)
      
      if(verbose == T) {
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
        
        if(plot_data == T) ad2_data_ad1 <- subset(ad2_data_ctry, ADM1_CODE == ad1_code)
        
        if (verbose == T) message(paste0("     --> ", unique(ad1_df_ad1$ADM1_NAME)))
        
        # Create a single df for plotting
        
        # Wrap facet labels & make sure ADMIN 1 comes first
        n_ad2s <- length(unique(ad2_df_ad1$ADM2_CODE))
        
        if (n_ad2s >  0 & n_ad2s <= 16) wrap_width <- 18
        if (n_ad2s > 16 & n_ad2s <= 25) wrap_width <- 12
        if (n_ad2s > 25)                wrap_width <- 9
        
        ad2_df_ad1[, plot_name := stringr::str_wrap(ADM2_NAME, width = wrap_width)]
        plot_df <- rbind(ad1_df_ad1, ad2_df_ad1, fill=T)
        ad2_lvls <- c("ADMIN 1", unique(ad2_df_ad1$plot_name))
        plot_df$plot_name <- factor(plot_df$plot_name, levels = ad2_lvls)
        
        if(plot_data == T) {
          ad2_data_ad1[, plot_name := stringr::str_wrap(ADM2_NAME, width = wrap_width)]
          ad2_data_ad1$plot_name <- factor(ad2_data_ad1$plot_name, levels = ad2_lvls)
        }
        
        # Admin 2 time-series plot -------------------------------------------------

        if (plot_data == F & multiple_runs == F) gg_ad2ts <- time_series(plot_df, admin = 2, val_range, title_plot_size, ind_title)
        if (plot_data == T & multiple_runs == F) gg_ad2ts <- time_series_data(plot_df, ad2_data_ad1, admin = 2, val_range, title_plot_size, ind_title)
        if (plot_data == F & multiple_runs == T) gg_ad2ts <- time_series_multiple(plot_df, admin = 2, val_range, title_plot_size, ind_title)
        if (plot_data == T & multiple_runs == T) gg_ad2ts <- time_series_multiple_data(plot_df, ad2_data_ad1, admin = 2, val_range, title_plot_size, ind_title)

        # A locator map ------------------------------------------------------------
        
        gg_locator_map <- ggplot() +
          geom_polygon_quiet(data = subset(ad1_national, ADM1_CODE == ad1_code),
                             aes(x = long, y = lat, group = group),
                             fill = "red") +
          geom_path_quiet(data = ad1_national,
                          aes(x=long, y=lat, group=group),
                          size = 0.2) +
          geom_path_quiet(data = ad0_national,
                          aes(x=long, y = lat, group=group),
                          size = 0.5) +
          coord_equal() +
          theme_empty
        
        
        # A labeled map of the admin1 unit by admin2 -----------------------------------------
        ad2_ad1 <- subset(ad2_national, ADM1_CODE == ad1_code)
        
        if(nrow(ad2_ad1) == 0){
          message("no polygon for ADM1_CODE", ad1_code, "skipping")
          next
        }
        
        centroids <- gCentroid(ad2_ad1, byid = T, id = ad2_ad1$ADM2_CODE)
        centroids <- as.data.frame(centroids) %>%
          cbind(rownames(.), .) %>%
          as.data.table(.) %>%
          setnames(., names(.), c("ADM2_CODE", "x", "y"))
        centroids$ADM2_CODE <- as.numeric(levels(centroids$ADM2_CODE))[centroids$ADM2_CODE]
        centroids <- merge(centroids, as.data.table(ad2_ad1), by = "ADM2_CODE")
        
        plot_df <- 
          ad2_df_ad1 %>% filter(year == max(year_list)) %>% dplyr::select(ADM2_CODE, mean)
        
        # Merge and fortify outside of ggplot to speed things up
        if (multiple_runs == T) plot_df <- plot_df %>% distinct() %>% data.table()
        if (multiple_runs == F) ad2_ad1 <- suppressWarnings(merge(ad2_ad1, plot_df, duplicateGeoms = T))
        
        ad2_ad1@data$id = rownames(ad2_ad1@data)
        ad2_ad1_df = fortify(ad2_ad1, region="id")
        ad2_ad1_df <- merge(ad2_ad1_df, ad2_ad1@data, by = "id")
        ad2_ad1_df <- as.data.table(ad2_ad1_df)
        
        if(highisbad) distiller_direction <- -1
        if(!highisbad) distiller_direction <- 1
        
        # Prepare for plotting
        if (multiple_runs == F) gg_ad2_map <- last_year_map(ad2_ad1_df, centroids, admin = 2, val_range, distiller_direction, max(year_list))
        if (multiple_runs == T) gg_ad2_map <- last_year_map_multiple(ad2_ad1_df, centroids, admin = 2, max(year_list))
        
        if (plot_data == F) gg_ad2_map <- aspect_ratio_plot(gg_ad2_map, 4/3, 0.25)
        
        # Set up the final layout ------------------------------------------------------------------
        # Create a title grob
        title_grob <- title_generate(ind_title, year_list, title_grob_size, admin = 2, ad1_df_ad1 = ad1_df_ad1, ctry = ctry)
        
        # Create the overall plot by arranging the grobs
        lay <- plot_overlay(plot_data, multiple_runs)
        if (plot_data == T & multiple_runs == F) master_plot <- arrangeGrob(gg_ad2ts, gg_locator_map, title_grob, gg_ad2_map, layout_matrix = lay)
        if (plot_data == T & multiple_runs == T) master_plot <- arrangeGrob(gg_ad2ts, gg_locator_map, title_grob, gg_ad2_map, ihme_grob_multiple, layout_matrix = lay)
        if (plot_data == F) master_plot <- arrangeGrob(gg_ad2ts, gg_locator_map, title_grob, gg_ad2_map, ihme_grob,
                                                       layout_matrix = lay, heights = c(1,1,1,1,1,0.6))
        
        grid.draw(master_plot)
        
        if (ad1_code != ad1_table$ADM1_CODE[length(ad1_table$ADM1_CODE)]) plot.new()
        
      } #END `for (ad1_code in ad1_table$ADM1_CODE)`
      
      dev.off() # For PDF (by country)
      
    } # END ` for (ctry in sort(countries))`
  } #END `if` wrapper for ad2
}

############################################################################################################################
## Utility functions used in visualization function ########################################################################
############################################################################################################################

## time_series plots  #############################

#' Functions to create time series plots at the specified admin level, faceting by the number of admin levels in that geography
#'
#' This function is called in the subnational_ts_function and is a way to generalize making time
#' series for each admin unit. These functions are all very similar but come in four flavors:
#'
#'
#' 1) time_series:  This plots the model predictions along with the uncertainty, faceted by the admin unit
#' 2) time_series_data: This plots the model predictions and uncertainty but adds the data aggregated to that admin level
#' (using input_aggregated_admin function), and it colors the points by the data source, as well as differentiating between
#' survey size and point and polygon with size and shape aesthetics respectively
#' 3) time_series_multiple: This plots multiple model predictions (either different model runs or multiple indicators, specified by model_runs = T in subnational_ts_plots)
#' on each time trend at the specified admin unit (0/1/2)
#' 4) time_series_multiple_data: This plots multiple model predictions and also includes the aggregated input data to that admin unit.
#' I only recommend using this if comparing different model runs of the same indicator, as plotting the data and multiple indicators makes the plots too busy
#'
#' @param model_data admin level data table from the standard aggregation code (typically ad0_df/ad1_df/ad2_df specified to certain geographic area)
#' @param input_data input data collapsed to admin level (typically ad0_data/ad1_data/ad2_data specified to certain geographic area)
#' @param admin specified admin level (0/1/2)
#' @param val_range Range of values for y axis
#' @param title_plot_size Size of plot title
#'
#' @return
#' This function returns a ggplot object
#'
#' @examples
#'
#' ad0_df_reg is the ad0_df specified to specific region, say for example 'cssa' (Central Sub-Saharan Africa)
#' if (plot_data == F & multiple_runs == F) ad0_time_series <- time_series(ad0_df_reg, admin = 0, val_range = c(0,1))

time_series = function(model_data,
                       admin = 0,
                       val_range = c(0,1),
                       title_plot_size = 30,
                       ind_title = ""){
  
  # Color palette
  carto_discrete <- rep(c("#7F3C8D","#11A579","#F2B701","#E73F74",
                          "#3969AC","#80BA5A","#E68310","#008695",
                          "#CF1C90","#f97b72","#4b4b8f","#A5AA99"), 3)
  gg_admin <-
    ggplot(data = model_data, aes(x = year, y = mean, ymin = lower, ymax = upper)) +
    geom_ribbon(alpha = 0.3) +
    geom_line() +
    theme_bw(base_size = 16) +
    coord_cartesian(ylim = val_range) +
    facet_wrap(~ plot_name) +
    theme(
      strip.background = element_blank(),
      plot.caption = element_text(hjust = 0.5),
      plot.title = element_text(size = title_plot_size, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(y = ind_title,
         x = "Year") +
         {if (admin == 0) labs(title = paste0(ind_title, " by Country"))} +
         {if (admin == 1) labs(title = paste0(ind_title, " by First-level Administrative Unit"),
                               caption = paste0("Time series depict first-level administrative units except for NATIONAL,\n",
                                                "which shows the time series for the entire country"))} +
                                                {if (admin == 2) labs(title = paste0(ind_title, " by Second-level Administrative Unit"),
                                                                      caption = paste0("Plots shown are for second-level administrative units except for ADMIN 1,\n",
                                                                                       "which shows the time series for the entire first-level administrative unit"))}
  return(gg_admin)
}

# Plot time series trends with aggregated input_data
time_series_data = function(model_data,
                            input_data,
                            admin = 0,
                            val_range = c(0,1),
                            title_plot_size = 30,
                            ind_title = ""){
  
  # Color palette
  carto_discrete <- rep(c("#7F3C8D","#11A579","#F2B701","#E73F74",
                          "#3969AC","#80BA5A","#E68310","#008695",
                          "#CF1C90","#f97b72","#4b4b8f","#A5AA99"), 3)
  gg_admin <-
    ggplot() +
    geom_ribbon(data = model_data, aes(x = year, ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line(data = model_data, aes(x = year, y = mean)) +
    theme_bw(base_size = 16) +
    geom_point(data = input_data,
               aes(x = year, y = outcome, size = N, shape = as.factor(point), fill = as.factor(source)), alpha = 0.7) +
    scale_fill_manual("Survey", values = carto_discrete, drop = F) +
    scale_shape_manual("Type", breaks = c("0", "1", "2", "3"), values = c(22, 21, 12, 10),
                       label = c("polygon", "point", "Subnationally \n Representative", "Subnationally \n Representative"), drop = F) +
    facet_wrap( ~ plot_name) +
    coord_cartesian(ylim = val_range) +
    scale_size(range = c(1,7)) +
    theme(
      strip.background = element_blank(),
      plot.caption = element_text(hjust = 0.5),
      plot.title = element_text(size = title_plot_size, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 7),
      legend.justification = "top"
    ) +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 5, shape = 21)),
           shape = guide_legend(order = 2, override.aes = list(size = 5))) +
    labs(y = ind_title,
         x = "Year") +
         {if (admin == 0) labs(title = paste0(ind_title, " by Country"))} +
         {if (admin == 1) labs(title = paste0(ind_title, " by First-level Administrative Unit"),
                               caption = paste0("Time series depict first-level administrative units except for NATIONAL,\n",
                                                "which shows the time series for the entire country"))} +
                                                {if (admin == 2) labs(title = paste0(ind_title, " by Second-level Administrative Unit"),
                                                                      caption = paste0("Plots shown are for second-level administrative units except for ADMIN 1,\n",
                                                                                       "which shows the time series for the entire first-level administrative unit"))}
  return(gg_admin)
}

# Plot time series plot with multiple indicators or multiple model runs, with run column specified
time_series_multiple <- function(model_data,
                                 admin = 0,
                                 val_range = c(0,1),
                                 title_plot_size = 30,
                                 ind_title = ""){
  
  # Color palette
  carto_discrete <- rep(c("#7F3C8D","#11A579","#F2B701","#E73F74",
                          "#3969AC","#80BA5A","#E68310","#008695",
                          "#CF1C90","#f97b72","#4b4b8f","#A5AA99"), 3)
  gg_admin <-
    ggplot(data = model_data, aes(x = year, y = mean, ymin = lower, ymax = upper)) +
    geom_ribbon(aes(fill = run), alpha = 0.1) +
    geom_line(aes(color = run, linetype = run), alpha = 0.9) +
    theme_bw(base_size = 16) +
    scale_color_manual("", values = carto_discrete) +
    scale_fill_manual("", values = carto_discrete) +
    scale_linetype_discrete("") +
    coord_cartesian(ylim = val_range) +
    facet_wrap(~plot_name) +
    theme(strip.background = element_blank(),
          plot.caption = element_text(hjust = 0.5),
          plot.title = element_text(size = title_plot_size, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.justification = "top") +
    labs(y = ind_title,
         x = "Year") +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 1),
           linetype = guide_legend(order = 1)) +
           {if (admin == 0) labs(title = paste0(ind_title, " by Country"))} +
           {if (admin == 1) labs(title = paste0(ind_title, " by First-level Administrative Unit"),
                                 caption = paste0("Time series depict first-level administrative units except for NATIONAL,\n",
                                                  "which shows the time series for the entire country"))} +
                                                  {if (admin == 2) labs(title = paste0(ind_title, " by Second-level Administrative Unit"),
                                                                        caption = paste0("Plots shown are for second-level administrative units except for ADMIN 1,\n",
                                                                                         "which shows the time series for the entire first-level administrative unit"))}
  return(gg_admin)
}

# Plot time series with multiple model runs, including data
time_series_multiple_data <- function(model_data,
                                      input_data,
                                      admin = 0,
                                      val_range = c(0,1),
                                      title_plot_size = 30,
                                      ind_title){
  
  # Color palette
  carto_discrete <- rep(c("#7F3C8D","#11A579","#F2B701","#E73F74",
                      "#3969AC","#80BA5A","#E68310","#008695",
                      "#CF1C90","#f97b72","#4b4b8f","#A5AA99"), 3)
                      
  gg_admin <-
    ggplot() +
    {if (length(unique(model_data$run)) <= 3) geom_ribbon(data = model_data,
                                                          aes(x = year, ymin = lower, ymax = upper, fill = run), alpha = 0.1)} +
                                                          {if (length(unique(model_data$run)) > 3) geom_ribbon(data = model_data,
                                                                                                               aes(x = year, ymin = lower, ymax = upper, fill = run), alpha = 0.03)} +
    geom_line(data = model_data, aes(x = year, y = mean, color = run, linetype = run), alpha = 0.8) +
    geom_point(data = input_data,
               aes(x = year, y = outcome, size = N, shape = point), fill = "grey", alpha = 0.7) +
    theme_bw(base_size = 16) +
    scale_color_manual("", values = carto_discrete) +
    scale_fill_manual("", values = carto_discrete) +
    scale_shape_manual("Type", breaks = c("0", "1", "2", "3"), values = c(22, 21, 12, 10),
                       label = c("polygon", "point", "Subnationally \n Representative", "Subnationally \n Representative"), drop = F) +
    scale_linetype_manual("", values = rep(c("solid", "dotdash", "dashed"), 5)) +
    coord_cartesian(ylim = val_range) +
    facet_wrap(~plot_name) +
    theme(strip.background = element_blank(),
          plot.caption = element_text(hjust = 0.5),
          plot.title = element_text(size = title_plot_size, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 7),
          legend.justification = "top") +
    labs(y = ind_title,
         x = "Year") +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 1),
           linetype = guide_legend(order = 1),
           shape = guide_legend(order = 2, override.aes = list(size = 5, fill = "grey"))) +
           {if (admin == 0) labs(title = paste0(ind_title, " by Country"))} +
           {if (admin == 1) labs(title = paste0(ind_title, " by First-level Administrative Unit"),
                                 caption = paste0("Time series depict first-level administrative units except for NATIONAL,\n",
                                                  "which shows the time series for the entire country"))} +
                                                  {if (admin == 2) labs(title = paste0(ind_title, " by Second-level Administrative Unit"),
                                                                        caption = paste0("Plots shown are for second-level administrative units except for ADMIN 1,\n",
                                                                                         "which shows the time series for the entire first-level administrative unit"))}
  return(gg_admin)
}


## location_map_draw  ###############################################################################################

#' Functions to create a locator map that highlights the admin level in question in red within a map background
#'
#' This function is called in the subnational_ts_function and is a way to generalize locator maps
#' that are useful to contextualize where the admin unit lies within a reigon or country.
#'
#' @param admin_shape SpatialPolygonsDataFrame specified to the specific admin unit, colored in red
#' @param surround_shape SpatialPolygonsDataFrame of surrounding geography to place admin unit.
#' If highlighted admin_shape is a country, then the surround_shape is typically the region it falls in
#'
#' @return
#' This function returns a ggplot object
#'
#' @examples
#'
#' ad0_shape_simple is the simplified ad0 shapefile for Africa, and ad0_reg_codes defines the
#' admin 0 codes for a specified region in Africa
#' admin_shape <- subset(ad0_shape_simple, ADM0_CODE %in% ad0_reg_codes)
#' surround_shape <- ad0_shape_simple
#' location_map <- location_map_draw(admin_shape, surround_shape)

location_map_draw <- function(admin_shape, surround_shape) {
  gg_location <-
    ggplot() +
    geom_polygon_quiet(data = admin_shape,
                       aes(x = long, y = lat, group = group),
                       fill = "red") +
    geom_path_quiet(data = surround_shape,
                    aes(x = long, y = lat, group = group),
                    size = 0.05,
                    alpha = 0.5) +
    geom_path_quiet(data = surround_shape,
                    aes(x = long, y = lat, group = group),
                    size = 0.3) +
    coord_equal() +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  #gg_location <- aspect_ratio_plot(gg_location, 4/3, 0.25)
  return(gg_location)
}

## last_year_map  ###############################################################################################

#' Functions to create a labeled map of the admin levels included in the region/country/admin1 level. The labels included are the
#' specified admin names unless over 25 geographies are included, in which case they are labeled by ihme loc id. This is done so larger regions
#' like Africa do not print out the entire country name, as it would take up too much space.
#'
#' This map comes in two flavors:
#' 1) last_year_map: This map colors the region/country/admin1 level by your indicator with the scale specified by val_range (default 0-1)
#' 2) last_year_map_multiple : When multiple model indicators/runs are included, the map is no longer colored by indicator and instead
#' all grey.
#'
#' This function is called in the subnational_ts_function and is pretty specific to this function. However, it could also serve as a
#' template for generalizing a labeled map function that also colors by a specific indicator.
#'
#' @param shape SpatialPolygonsDataFrame specified to the region/country/admin1 level that your model was run over
#' @param centroids Data table that species the x and y coordinates of the centroids for the shape of the map
#'                  and has the attached labels to display using ggrepel.
#' @param admin The specified admin level called in subnational_ts_function
#' @param val_range The specified value range for mean value of your indicator
#' @param distiller_direction Either 1 or -1, if high is bad, distiller direction is -1 with the red tones at high values of indicator
#' @param year This specifies the last year of your time series, which is used for creating a legend for the function
#'
#'
#' @return
#' This function returns a ggplot object
#'
#' @examples
#'
#' If ad0_shape_reg_df is a fortified data frame over a specific region, and centroids are the centroids for the defined
#' countries in that region labeled by country name.
#' gg_lastyear_map <- last_year_map(ad0_shape_reg_df, centroids, admin = 0)


last_year_map <- function(shape, centroids, admin = 0, val_range = c(0,1), distiller_direction = 1, year = 2015) {
  # Requirements and custom theme
  require(ggrepel)
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  
  gg <- ggplot() +
    geom_polygon_quiet(data = shape,
                       aes(x=long, y=lat,
                           group=group, fill = mean)) +
    geom_path_quiet(data = shape,
                    aes(x=long, y=lat, group=group)) +
    coord_equal() +
    theme_empty +
    labs(fill = paste0("Mean\nestimate","\n(", year, ")"))
  
  if (admin == 0 & nrow(centroids) > 25) {
    gg <- gg + geom_label_repel(data = centroids,
                                aes(x=x,
                                    y=y,
                                    label = ihme_lc_id),
                                point.padding = unit(0.01, "lines"),
                                box.padding = unit(1.5, "lines"),
                                min.segment.length = unit(0, "lines"),
                                segment.alpha = 0.5) +
      scale_fill_distiller(palette = "RdYlBu",
                           limits = val_range,
                           direction = distiller_direction)
  } else {
    label <- paste0("ADM", admin, "_NAME")
    gg = gg + geom_label_repel(data = centroids,
                               aes(x=x,
                                   y=y,
                                   label =  get(label)),
                               point.padding = unit(0.01, "lines"),
                               box.padding = unit(1.5, "lines"),
                               min.segment.length = unit(0, "lines"),
                               segment.alpha = 0.5) +
      scale_fill_distiller(palette = "RdYlBu",
                           limits = val_range,
                           direction = distiller_direction)
  }
  return(gg)
}

# Last year maps if multiple runs are provided, this changes the fill to grey
last_year_map_multiple <- function(shape, centroids, admin = 0, year = 2015) {
  # Requirements and custom theme
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  gg <- ggplot() +
    geom_polygon_quiet(data = shape,
                       aes(x=long, y=lat,
                           group=group), fill = "grey") +
    geom_path_quiet(data = shape,
                    aes(x=long, y=lat, group=group)) +
    coord_equal() +
    theme_empty
  
  if (admin == 0 & nrow(centroids) > 25) {
    gg <- gg + geom_label_repel(data = centroids,
                                aes(x=x,
                                    y=y,
                                    label = ihme_lc_id),
                                point.padding = unit(0.01, "lines"),
                                box.padding = unit(1.5, "lines"),
                                min.segment.length = unit(0, "lines"),
                                segment.alpha = 0.5)
  } else {
    label <- paste0("ADM", admin, "_NAME")
    gg = gg + geom_label_repel(data = centroids,
                               aes(x=x,
                                   y=y,
                                   label =  get(label)),
                               point.padding = unit(0.01, "lines"),
                               box.padding = unit(1.5, "lines"),
                               min.segment.length = unit(0, "lines"),
                               segment.alpha = 0.5)
  }
  return(gg)
}

# utility functions ###########################################################################

# Quiet version of geom polygon and path to supress messages
geom_polygon_quiet <- function(...) {suppressMessages(ggplot2::geom_polygon(...))}
geom_path_quiet     <- function(...) {suppressMessages(ggplot2::geom_path(...))}


#' This function is used to add a "run" column specifying a model run or model indicators to be used in multiple runs, when
#' plotting model predictions of multiple runs. The input is a list of model prediction of a specific admin unit, and the output is
#' a data table that binds the predictions together with a column titled "run" that specified what indicator/model run the model
#' prediction is from, dependent on the run label.
#'
#' @param prediction_list List where each component is a model prediction at a specified admin level
#' @param run_label This labels each model run/different indicator so that we can eventually facet by the run
#'                  column when creating a ggplot2 object
#'
#'
#' @return
#' This function returns a single data table that binds all of prediction list components together, adding the appropriate "run" column
#'
#' @examples
#' # Example of code used for adding multiple indicators to a subnational_ts_plot
#' run_dates <- c("2018_05_30_23_01_59", "2018_03_23_16_36_23")
#' run_label <- c("STI symp", "Condom use")
#' indicators <- c("condom_last_time", "sti_symptoms")
#'
#' share_dir <- '<<<< FILEPATH REDACTED >>>>'
#' in_dir <- '<<<< FILEPATH REDACTED >>>>'
#' out_dir <- '<<<< FILEPATH REDACTED >>>>'
#'
#' in_file_ad0 <- paste0(in_dir, indicators, "_admin_0_unraked_summary.csv")
#' in_file_ad1 <- paste0(in_dir, indicators, "_admin_1_unraked_summary.csv")
#' in_file_ad2 <- paste0(in_dir, indicators, "_admin_2_unraked_summary.csv")
#'
#' ad0_df <- lapply(in_file_ad0, fread) %>% add_run_label(run_label)
#' ad1_df <- lapply(in_file_ad1, fread) %>% add_run_label(run_label)
#' ad2_df <- lapply(in_file_ad2, fread) %>% add_run_label(run_label)\
#' These inputs are now ready to be passed into subnational_ts_plots with multiple_runs = T

add_run_label <- function(prediction_list, run_label = NULL){
  for (i in 1:length(prediction_list)){
    prediction_list[[i]][["run"]] <- run_label[i]
  }
  x <- rbindlist(prediction_list)
  x[, run := factor(run, run_label)]
  return(x)
}



## plot_ordered_admins_colors ################################################

#' "Traffic light" plots and maps to help visualize categories of administrative estimates
#'
#' This is currently written to work for three broad categories (in the original case,
#' vaccine coverage 0-60%, 60-90%, and 90-100%) which are colored red, yellow, and blue.
#' Borderline estimates (i.e. with CIs that overlap 60% or 90% in the original example) are
#' colored orange and green.  If CIs extend into all 3 categories, then the estimate
#' is colored gray.
#'
#' Currently, this is relatively inflexible - ideally would be extended to include
#' other category/color schemes and work with admin1s.  Feel free to submit a PR!
#'
#' @param P Param description
#' @param P Param description
#' @return
#' @examples
#'

plot_ordered_admins_colors <- function(df,
                                       ctry,
                                       cutoffs = c(0, 0.6, 0.9, 1.0),
                                       error_bars = "yaxis",
                                       ad2_shp = NULL,
                                       ad0_shp = NULL,
                                       shapefile_version = 'current') {
  
  ################################################################################
  # 0. SETUP #####################################################################
  ################################################################################
  
  # Check inputs
  if (length(cutoffs) != 4) stop("Currently only three-bin estimates supported - length(cutoffs) should be 4")
  
  # Use GAUL shapefiles as defaults
  if (is.null(ad0_shp)) ad0_shp <- rgdal::readOGR(get_admin_shapefile(admin_level=0, version = shapefile_version))
  if (is.null(ad2_shp)) ad2_shp <- rgdal::readOGR(get_admin_shapefile(admin_level=2, version = shapefile_version))
  
  # Subset df to the country of interest and sort
  if (ctry != "all") {
    df <- copy(subset(df, ihme_lc_id == ctry))
    country_title <- unique(df$ADM0_NAME)
  } else {
    country_title <- "Africa (all countries)"
  }
  
  df[, sort_mean := "mean"]
  
  xaxis_title <- "Second-level administrative unit (ordered by mean MBG estimate)"
  sort_table <- subset(df[order(mean)], select = c("mean", "ADM2_CODE"))
  sort_table[, plot_idx := .I]
  plot_df <- merge(df, subset(sort_table, select = c("ADM2_CODE", "plot_idx")), by = "ADM2_CODE")
  
  # Assign color categories
  plot_df[upper < 0.6,                                 category := "Coverage significantly < 60 %"]
  plot_df[lower < 0.6 & upper >= 0.6 & upper <= 0.9,   category := "Coverage near 60% with overlapping UIs"]
  plot_df[lower > 0.6 & upper < 0.9,                   category := "Coverage significantly above 60% and below 90%"]
  plot_df[lower >= 0.6 & lower <= 0.9 & upper > 0.9,   category := "Coverage near 90% with overlapping UIs"]
  plot_df[lower > 0.9,                                 category := "Coverage significantly > 90%"]
  plot_df[is.na(category),                             category := "Coverage UIs overlap all categories"]
  
  cat_levels <- c("Coverage significantly < 60 %",
                  "Coverage near 60% with overlapping UIs",
                  "Coverage significantly above 60% and below 90%",
                  "Coverage near 90% with overlapping UIs",
                  "Coverage significantly > 90%",
                  "Coverage UIs overlap all categories")
  
  plot_df$category <- factor(plot_df$category, levels = cat_levels)
  
  color_scheme <- c("#FF0000", "#FF9900", "#FFD700", "#006400", "#00008B", "#A9A9A9")
  names(color_scheme) <- cat_levels
  
  gg_p <- ggplot(data = plot_df,
                 aes(x = plot_idx,
                     y = mean,
                     color = category))
  
  if (error_bars == "yaxis" | error_bars == "both") {
    gg_p <- gg_p + geom_errorbar(data = plot_df,
                                 aes(x = plot_idx, ymin = lower, ymax = upper), width = 0, alpha = 0.4)
  }
  
  gg_p <- gg_p +
    geom_abline(slope = 0, intercept = 0.6, color = "black", linetype = "dashed", alpha = 0.6) +
    geom_abline(slope = 0, intercept = 0.9, color = "black", linetype = "dashed", alpha = 0.6) +
    geom_point(aes(size = pop_last), alpha = 0.7) +
    theme_classic() +
    labs(y = paste0("Mean DPT3 Coverage: 2016"),
         x = xaxis_title,
         title = country_title,
         shape = "Source",
         color = "Source",
         size = "Population: children < 5 years (2016)") +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks = seq(0,1,0.1), labels = scales::percent) +
    scale_colour_manual(name = "Category",values = color_scheme, drop = F) +
    scale_size_area(labels = scales::comma) +
    theme(legend.position="right") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.caption = element_text(hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(linetype = 0), order = 1),
           size = guide_legend(order = 2))
  
  # Make a map
  
  # Subset the shape
  ctry_code <- unique(plot_df$ADM0_CODE)
  plot_shape <- subset(ad2_shp, ADM0_CODE == ctry_code)
  plot_shape <- merge(plot_shape, subset(plot_df, select = c("ADM2_CODE", "category")))
  background_shape <- subset(ad0_shp, ADM0_CODE == ctry_code)
  
  plot_shape@data$id = rownames(plot_shape@data)
  plot_shape_df = fortify(plot_shape, region="id")
  plot_shape_df <- merge(plot_shape_df, plot_shape@data, by = "id")
  plot_shape_df <- as.data.table(plot_shape_df)
  
  # A custom, mostly blank theme to use for mapping
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  # Make plot
  gg_map <- ggplot() +
    geom_polygon(data = background_shape,
                 aes(x = long,
                     y = lat,
                     group = group),
                 fill = "white") +
    geom_polygon(data = plot_shape_df,
                 aes(x = long,
                     y = lat,
                     group = group,
                     fill = category)) +
    geom_path(data = plot_shape,
              aes(x = long,
                  y = lat,
                  group = group),
              size = 0.2,
              color = "black") +
    geom_path(data = background_shape,
              aes(x = long,
                  y = lat,
                  group = group),
              size = 0.6,
              color = "black") +
    theme_empty +
    scale_size_area(max_size = max_pt_size, labels = comma) +
    coord_equal(ratio = 1) +
    scale_fill_manual(name = "Category",values = color_scheme, drop = F) +
    guides(fill = F)
  
  return(list(scatter = gg_p,
              map = gg_map))
  
}

format_plot_obj <- function(plot_obj) {
  
  g_legend<-function(a.gplot){
    pdf(NULL) # Workaround for bug in ggplot_gtable causing empty Rplots.pdf to be created
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    graphics.off()
    return(legend)
  }
  
  legend <- g_legend(plot_obj$scatter)
  scatter <- plot_obj$scatter + guides(color = F, size = F)
  map <- plot_obj$map
  
  lay <- rbind(c(1,1,1,1,2,2),
               c(1,1,1,1,2,2),
               c(1,1,1,1,3,3),
               c(1,1,1,1,3,3))
  
  return(grid.arrange(scatter, legend, map, layout_matrix = lay))
  
}


## plot_ordered_admins_colors ################################################

#' "Traffic light" plots and maps to help visualize categories of administrative estimates
#'
#' This is currently written to work for three broad categories (in the original case,
#' vaccine coverage 0-60%, 60-90%, and 90-100%) which are colored red, yellow, and blue.
#' Borderline estimates (i.e. with CIs that overlap 60% or 90% in the original example) are
#' colored orange and green.  If CIs extend into all 3 categories, then the estimate
#' is colored gray.
#'
#' Currently, this is relatively inflexible - ideally would be extended to include
#' other category/color schemes and work with admin1s.  Feel free to submit a PR!
#'
#' @param df data.frame or data.table from the aggregation code.  Needs to have
#'           columns for `mean`, `upper`, `lower`, `ADM2_CODE`, `ADM0_CODE`, `year`, and also a
#'           also a `pop` column which can be obtained as in the example (from merging on the
#'           pops from the admin draws object)
#' @param ad0_name the `ADM0_NAME` that you're interested in plotting. Can loop over this
#'                 in a wrapper function to be able to produce plots for many countries
#'                 sequentially
#' @param plot_year which year in df would you like to plot?
#' @param indicator title for your indicator, e.g. "DPT3 Coverage"
#' @param pop_title title for your population, e.g. "Population: children < 5 years"
#' @param cutoffs vector of cutoffs for your plot. For instance `c(0, 0.6, 0.9, 1)`.
#'                Must be 4 items (3 bins) for now
#' @param error_bars would you like error bars? Either "yaxis" or "none"
#' @param ad2_shp admin2 shape for plotting
#' @param ad0_shp admin0 shape for national borders
#' @param verbose more updates from the function (Boolean)
#' @param shapefile_version string indicating version of shapefile to pull
#'
#' @return a grob object for the full plot. Note: need to use `grid.draw()` to plot grobs!
#'
#' @examples
#'
#' # Set up directories and files  ########################################
#' # Set `run_date`, `indicator`, `indicator_group`, `out_dir` as you wish
#'
#' share_dir <- '<<<< FILEPATH REDACTED >>>>'
#' in_dir <- '<<<< FILEPATH REDACTED >>>>'
#' in_file_ad2 <- '<<<< FILEPATH REDACTED >>>>'
#' ad2_df <- fread(in_file_ad2)
#'
#' # Merge on populations
#' load(paste0(share_dir, indicator, "_raked_admin_draws_eb_bin0_0.RData"))
#' pops <- subset(admin_2, year %in% unique(ad2_df$year), select = c("ADM2_CODE", "pop", "year"))
#' ad2_df <- merge(ad2_df, pops, by = c("ADM2_CODE", "year"), all.x= T, all.y = F)
#' rm(admin_2)
#'
#' # Use GAUL shapefiles as defaults
#' ad0_shp <- readRDS("<<<< FILEPATH REDACTED >>>>")
#' ad2_shp <- readOGR(get_admin_shapefile(admin_level=2))
#'
#' a_grob <- plot_ordered_admins_colors(df = ad2_df,
#'                                      ad0_name = "Senegal",
#'                                      plot_year = 2016,
#'                                      indicator_title = "DPT3 Coverage",
#'                                      pop_title = "Population: Children < 5 years",
#'                                      cutoffs = c(0,0.6, 0.9, 1),
#'                                      error_bars = "yaxis",
#'                                      ad2_shp = ad2_shp,
#'                                      ad0_shp = ad0_shp,
#'                                      verbose = T)
#'
#' grid.draw(a_grob)
#'

plot_ordered_admins_colors <- function(df,
                                       ad0_name,
                                       plot_year,
                                       indicator_title,
                                       pop_title,
                                       cutoffs,
                                       error_bars = "yaxis",
                                       ad2_shp = NULL,
                                       ad0_shp = NULL,
                                       verbose = F,
                                       shapefile_version = 'current'){
  
  ################################################################################
  # 0. SETUP #####################################################################
  ################################################################################
  
  if (verbose) message(paste0("Working on ", ad0_name, "..."))
  
  # Replace some functions to be quieter
  geom_polygon_quiet <- function(...) {suppressMessages(ggplot2::geom_polygon(...))}
  geom_path_quiet     <- function(...) {suppressMessages(ggplot2::geom_path(...))}
  
  # Check inputs
  if (length(cutoffs) != 4) stop("Currently only three-bin estimates supported - length(cutoffs) should be 4")
  
  # Use GAUL shapefiles as defaults
  if (is.null(ad0_shp)) ad0_shp <- rgdal::readOGR(get_admin_shapefile(admin_level=0, version = shapefile_version))
  if (is.null(ad2_shp)) ad2_shp <- rgdal::readOGR(get_admin_shapefile(admin_level=2, version = shapefile_version))
  
  # Subset df to the country of interest & year and sort
  df <- copy(as.data.table(df))
  df <- subset(df, year == plot_year)
  if (ad0_name != "all") {
    df <- copy(subset(df, ADM0_NAME == ad0_name))
    country_title <- unique(df$ADM0_NAME)
  } else {
    country_title <- "Africa (all countries)"
  }
  
  xaxis_title <- "Second-level administrative unit (ordered by mean MBG estimate)"
  sort_table <- subset(df[order(mean)], select = c("mean", "ADM2_CODE"))
  sort_table[, plot_idx := .I]
  plot_df <- merge(df, subset(sort_table, select = c("ADM2_CODE", "plot_idx")), by = "ADM2_CODE")
  
  # Assign color categories
  c_2 <- cutoffs[2]
  c_3 <- cutoffs[3]
  clab_2 <- paste0(c_2 * 100, "%")
  clab_3 <- paste0(c_3 * 100, "%")
  
  cat_levels <- c(paste0("Significantly < ", clab_2),
                  paste0("Near ", clab_2, " with overlapping UIs"),
                  paste0("Significantly above ", clab_2, " and below ", clab_3),
                  paste0("Near ", clab_3, " with overlapping UIs"),
                  paste0("Significantly > ", clab_3),
                  "UIs overlap all categories")
  
  plot_df[upper < c_2,                                 category := cat_levels[1]]
  plot_df[lower < c_2 & upper >= c_2 & upper <= c_3,   category := cat_levels[2]]
  plot_df[lower > c_2 & upper < c_3,                   category := cat_levels[3]]
  plot_df[lower >= c_2 & lower <= c_3 & upper > c_3,   category := cat_levels[4]]
  plot_df[lower > c_3,                                 category := cat_levels[5]]
  plot_df[is.na(category),                             category := cat_levels[6]]
  
  plot_df$category <- factor(plot_df$category, levels = cat_levels)
  
  # Set up a color scheme
  color_scheme <- c("#FF0000", "#FF9900", "#FFD700", "#006400", "#00008B", "#A9A9A9")
  names(color_scheme) <- cat_levels
  
  ################################################################################
  # I. Scatter plot #################################################################
  ################################################################################
  
  gg_p <- ggplot(data = plot_df,
                 aes(x = plot_idx,
                     y = mean,
                     color = category))
  
  if (error_bars == "yaxis" | error_bars == "both") {
    gg_p <- gg_p + geom_errorbar(data = plot_df,
                                 aes(x = plot_idx, ymin = lower, ymax = upper), width = 0, alpha = 0.4)
  }
  
  gg_p <- gg_p +
    geom_abline(slope = 0, intercept = c_2, color = "black", linetype = "dashed", alpha = 0.6) +
    geom_abline(slope = 0, intercept = c_3, color = "black", linetype = "dashed", alpha = 0.6) +
    geom_point(aes(size = pop), alpha = 0.7) +
    theme_classic(base_size = 16) +
    labs(y = paste0("Mean ", indicator_title, ": ", plot_year),
         x = xaxis_title,
         title = country_title,
         shape = "Source",
         color = "Source",
         size = pop_title) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,1), breaks = seq(0,1,0.1), labels = scales::percent) +
    scale_colour_manual(name = "Category",values = color_scheme, drop = F) +
    scale_size_area(labels = scales::comma) +
    theme(legend.position="right") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          plot.caption = element_text(hjust = 0.5)) +
    guides(color = guide_legend(override.aes = list(linetype = 0), order = 1),
           size = guide_legend(order = 2))
  
  ################################################################################
  # Create a map #################################################################
  ################################################################################
  
  # Subset the shape
  ctry_code <- unique(plot_df$ADM0_CODE)
  plot_shape <- subset(ad2_shp, ADM0_CODE == ctry_code)
  plot_shape <- merge(plot_shape, subset(plot_df, select = c("ADM2_CODE", "category")))
  background_shape <- subset(ad0_shp, ADM0_CODE == ctry_code)
  
  plot_shape@data$id = rownames(plot_shape@data)
  plot_shape_df = fortify(plot_shape, region="id")
  plot_shape_df <- merge(plot_shape_df, plot_shape@data, by = "id")
  plot_shape_df <- as.data.table(plot_shape_df)
  
  # A custom, mostly blank theme to use for mapping
  theme_empty <- theme_classic() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 30, face = "bold"),
          strip.background = element_blank(),
          strip.text.x = element_blank())
  
  # Make plot
  gg_map <- ggplot() +
    geom_polygon_quiet(data = background_shape,
                       aes(x = long,
                           y = lat,
                           group = group),
                       fill = "white") +
    geom_polygon_quiet(data = plot_shape_df,
                       aes(x = long,
                           y = lat,
                           group = group,
                           fill = category)) +
    geom_path_quiet(data = plot_shape,
                    aes(x = long,
                        y = lat,
                        group = group),
                    size = 0.2,
                    color = "black") +
    geom_path_quiet(data = background_shape,
                    aes(x = long,
                        y = lat,
                        group = group),
                    size = 0.6,
                    color = "black") +
    theme_empty +
    scale_size_area(max_size = max_pt_size, labels = comma) +
    coord_equal(ratio = 1) +
    scale_fill_manual(name = "Category",values = color_scheme, drop = F) +
    guides(fill = F)
  
  ################################################################################
  # Put it all together ##########################################################
  ################################################################################
  
  # Function to reorganize the plot
  
  format_plot_obj <- function(map, scatter) {
    
    g_legend<-function(a.gplot){
      pdf(NULL) # Workaround for bug in ggplot_gtable causing empty Rplots.pdf to be created
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      graphics.off()
      return(legend)
    }
    
    legend <- g_legend(scatter)
    scatter <- scatter + guides(color = F, size = F)
    
    lay <- rbind(c(1,1,1,1,2,2),
                 c(1,1,1,1,2,2),
                 c(1,1,1,1,3,3),
                 c(1,1,1,1,3,3))
    
    return(arrangeGrob(scatter, legend, map, layout_matrix = lay))
    
  }
  
  output_grob <- format_plot_obj(gg_map, gg_p)
  return(output_grob)
}

## input_aggregate_admin function ###################################################

#' This function correctly assigns and then collapsed input data to the admin0/ admin1/ admin2 level
#' to be used in the visualization function and add input data to plots
#'
#' This function pulls the input data from the model directory, however it requires that you include
#' a column that retains the sum of the sample weights when collapsing to lat/long or shapefiles in your collapse
#' code in order to correctly recollapse to the specified admin. For this reason, it is important to include this
#' information in the collapse code before using this function. You can name the column whatever you want, though the
#' default name is sum_of_sample_weights. This is an argument to the function, sample_column. The input data must also have
#' a column labeled "point" that is an indicator for whether the data is point data (1) or polygon data (0). 
#' 
#' You may also specify the input
#' data yourself by passing the loaded input data as an argument to the function, do this if your input data is not specified
#' in the model directory, though this should be the case for most using the mbg pipeline. 
#'
#' This function uses the sf package, which will currently only work when running on the lbd_singularity image or in
#' Rstudio on the cluster.
#'
#' @param indicator indicator name used in file structure for mbg
#' @param indicator_group indicator group
#' @param regions Regions specified from your model
#' @param run_date  model run date
#' @param input_data If specified, provides the preloaded input data so the function does not look in your model directory
#' @param indicator_family If specified as Gaussian, this makes sure to not divide the prevalence by N which is required for binomial indicators
#' @param svy_id This is the unique survey id, usually labeled "nid" but other teams might use different terminology
#' @param sample_column This is the name of the column that contains the sum of the sample weights for the collapsed points.
#' @param subnational_nids If specified, these are the nids of surveys that are not nationally representative, and appear on the data plots with a different shape
#' @param shapefile_version String indicating version of shapefile to pull
#' @param use_kish boolean - If TRUE, admin aggregated N's are calculated using Kish approx
#'
#' @return a list containing three data tables that correspond to collapsed admin levels, named ad0/ad1/ad2. This list should then
#' be assigned to parameters ad0_data/ad1_data/ad2_data for the visualization function if plot_data == T. The admin data will have columns
#' named outcome
#'
#' @examples
#'
#' # Set up directories and files  ########################################
#' # Set `run_date`, `indicator`, `indicator_group`, `sample_column` as you wish
#'
#' admin_data <- input_aggregate_admin(indicator = indicator, indicator_group, regions = c("cssa", "wssa", "essa", "sssa"), indicator_family = "binomial", sample_column = "sum_of_sample_weights")
#' ad0_data <- admin_data$ad0
#' ad1_data <- admin_data$ad1
#' ad2_data <- admin_data$ad2
#'
#' Running the plotting code with plot_data = T ##########################
#' #' subnational_ts_plots(ad0_df = ad0_df,
#'                      ad1_df = ad1_df,
#'                      ad2_df = ad2_df,
#'                      ad0_shape = ad0_shape,
#'                      ad1_shape = ad1_shape,
#'                      ad2_shape = ad2_shape,
#'                      ind_title = "My Indicator",
#'                      out_dir = out_dir,
#'                      highisbad = F,
#'                      val_range = c(0,1),
#'                      ad0_map_regions = c("cssa", "essa", "name", "sssa", "wssa"),
#'                      ad0_map_region_titles = c("Central Sub-Saharan Africa",
#'                                                "Eastern Sub-Saharan Africa",
#'                                                "Northern Africa",
#'                                                "Southern Sub-Saharan Africa",
#'                                                "Western Sub-Saharan Africa"),
#'                      verbose = T,
#'                      plot_data = T,
#'                      ad0_data = ad0_data,
#'                      ad1_data = ad1_data,
#'                      ad2_data = ad2_data)
#'
#'

input_aggregate_admin <- function(indicator,
                                  indicator_group,
                                  regions = c("cssa", "wssa", "essa", "sssa", "name"),
                                  run_date = NULL,
                                  input_data = NULL,
                                  indicator_family = "binomial",
                                  svy_id = "nid",
                                  sample_column = "sum_of_sample_weights",
                                  subnational_nids = NULL,
                                  shapefile_version = 'current',
                                  use_kish = FALSE) {
  
  # Make sure this is run on singularity container in Rstudio or lbd singularity to use sf package
  if(!is_lbd_singularity() & !is_rstudio(check_singularity = TRUE)){
    stop("must be run in lbd singularity or Rstudio singularity container")
  }
  
  require(dplyr)
  require(data.table)
  require(sf)
  require(ggplot2)
  require(stringr)
  
  # If input_data not passed to the function, load in input data from model, or from mbg if the appropriate column name is not specified.
  if (is.null(input_data)){
    mod_dir <- '<<<< FILEPATH REDACTED >>>>'
    message(paste0("Input data was not provided, so loading from model directory: ", mod_dir,
                   ". \n If this does not exist, consider passing input data as an argument to this function"))
    input_data <- fread(paste0(mod_dir, 'input_data.csv'))
    
    if (!sample_column %in% names(input_data)) {
      message("Sample weights column not found in model directory, pulling from input data on <<<< FILEPATH REDACTED >>>>")
      input_data <- fread(paste0("<<<< FILEPATH REDACTED >>>>", indicator, ".csv"))
      if (!sample_column %in% names(input_data)) stop("Sample weights column not found on J drive, make sure to add a sum of sample weights column)                                         in collapse code and specify columns name as sample_column argument in function")
    }
  } else {
    if (!sample_column %in% names(input_data)) {
      stop("Sample weights column not found in the provided input data. Make sure that this column is included and specified by the sample_column
           argument in the function call.")
    }
    }
  
  if (!"source" %in% names(input_data)) stop("Need to specify the source of the data under a source column in the input data")
  if (!"point"  %in% names(input_data)) stop("You need a column in your input data called 'points' that classifies point data as 1 and polygon as 0")
  
  # Load in shapefile
  admin_shp <- sf::st_read(get_admin_shapefile(admin_level = 2, version = shapefile_version), quiet=T)
  
  # Add country name (not just 3 letter abbreviation) to input data
  # Subset input data to given region defined by the regions argument 
  gaul_codes <- get_adm0_codes(regions, shapefile_version = shapefile_version) 
  
  gaul_to_loc_id <- 
    get_location_code_mapping(shapefile_version = shapefile_version) %>% 
    dplyr::select(GAUL_CODE, loc_name, ihme_lc_id) %>% 
    dplyr::rename(location_name = loc_name) %>%
    filter(GAUL_CODE %in% gaul_codes)
  
  #edit gaul_to_loc_id to use location names from admin shapefile, because the subnational viz function merges on names
  admin_shp_data <- as.data.table(admin_shp)
  admin_shp_data <- unique(admin_shp_data[,c("ADM0_CODE", "ADM0_NAME")])
  gaul_to_loc_id <- gaul_to_loc_id[, c("GAUL_CODE", "ihme_lc_id")]
  gaul_to_loc_id <- merge(gaul_to_loc_id, admin_shp_data, by.x = "GAUL_CODE", by.y = "ADM0_CODE")
  setnames(gaul_to_loc_id, "ADM0_NAME", "location_name")
  
  # Join input data to location names to match shapefile. Some of those countries need to be manually changed to fit shapefile names
  input_data <-
    input_data %>%
    left_join(gaul_to_loc_id, by = c("country" = "ihme_lc_id")) %>%
    rowwise() %>% 
    ungroup() %>% data.table() %>%
    setnames(c(svy_id, sample_column), c("svy_id", "sample_column"))
  
  missing <-
    input_data %>%
    filter(is.na(location_name)) %>% nrow()
  
  message(paste0(round(100* missing / nrow(input_data), 2), " % of data is from outside of specified regions: ", paste(regions, collapse = " ")))
  
  input_data <-
    input_data %>% filter(!is.na(location_name))
  
  # Subset shapefile to include countries we have data on
  countries <- input_data$location_name %>% unique()
  
  admin_shp <-
    admin_shp %>%
    mutate(ADM0_NAME = as.character(ADM0_NAME)) %>%
    filter(ADM0_NAME %in% countries)
  
  # Assign input data to correct admin0/admin1/dmin2 in one step
  message("Assigning lat/longs to correct admin level")
  input_admin <-
    input_data %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = st_crs(admin_shp)) %>%
    st_join(admin_shp) %>%
    st_set_geometry(NULL) %>%
    setnames(eval(indicator), "prev")
  
  missing <-
    input_admin %>%
    filter(is.na(ADM0_CODE), is.na(ADM1_CODE), is.na(ADM2_CODE)) %>% nrow()
  
  message(paste0(round(100* missing / nrow(input_admin), 3), " % of data could not be matched to an admin level"))
  
  input_admin <-
    input_admin %>%
    filter(!is.na(ADM0_CODE), !is.na(ADM1_CODE), !is.na(ADM2_CODE)) %>% data.table() #remove those that are not assigned to an admin0/admin2
  
  # Make sure they are assigned to the correct country, some are right over the border and so are pushed into a different area.
  # These are excluded (usually not substantial) becuase they then do not nest well with admin1 and admin2 estimates
  wrong_admin0 <- nrow(input_admin[location_name != ADM0_NAME])
  message(paste0(round(wrong_admin0/ nrow(input_admin), 3), " % of input data is matched to a different country.\nThese are usually located on the border and will be dropped for the visualization."))
  
  input_admin <- input_admin[location_name == ADM0_NAME]
  
  # If binomial make sure it is prevalence space
  if (indicator_family == "binomial") input_admin[, prev := prev / N]
  
  # Collapse to admin 0 level
  message("collapsing to admin 0 level")
  input_admin0 <-
    input_admin %>%
    group_by(svy_id, source, point, ADM0_NAME, ADM0_CODE) %>%
    dplyr::summarise(
      year = floor(median(year, na.rm = T)),
      outcome = weighted.mean(prev, sample_column),
      N = ifelse(use_kish, 
                 sum(as.numeric(sample_column))^2 / sum(as.numeric(sample_column^2 / N)), 
                 sum(N * weight))
    ) %>%
    ungroup() %>%
    data.table()
  
  # Collapse to admin 1 level
  message("collapsing to admin 1 level")
  input_admin1 <-
    input_admin %>%
    group_by(svy_id, source, point, ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE) %>%
    dplyr::summarise(
      year = floor(median(year, na.rm = T)),
      outcome = weighted.mean(prev, sample_column),
      N = ifelse(use_kish, 
                 sum(as.numeric(sample_column))^2 / sum(as.numeric(sample_column^2 / N)), 
                 sum(N * weight))
    ) %>%
    ungroup() %>%
    data.table()
  
  # Collapse to admin 2 level
  message("collapsing to admin 2 level")
  input_admin2 <-
    input_admin %>%
    group_by(svy_id, source, point, ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE) %>%
    dplyr::summarise(
      year = floor(median(year, na.rm = T)),
      outcome = weighted.mean(prev, sample_column),
      N = ifelse(use_kish, 
                 sum(as.numeric(sample_column))^2 / sum(as.numeric(sample_column^2 / N)), 
                 sum(N * weight))
    ) %>%
    ungroup() %>%
    data.table()
  
  # Change source/polygon to factor and shorten survey names to fit on plot legend
  input_clean <- function(df) {
    df %>%
      rowwise() %>% 
      mutate(source = ifelse(source == "MACRO_DHS", "DHS", source)) %>%
      mutate(source = ifelse(source == "MACRO_AIS", "AIS", source)) %>%
      mutate(source = ifelse(source == "UNICEF_MICS", "MICS", source)) %>%
      mutate(source = ifelse(source == "COUNTRY_SPECIFIC", "CS", source)) %>%
      mutate(source = ifelse(source == "WB_CWIQ", "CWIQ", source)) %>%
      mutate(source = ifelse(source == "WB_CWIQ", "CWIQ", source)) %>%
      mutate(source = ifelse(source == "WB_LSMS", "LSMS", source)) %>%
      mutate(source = ifelse(source == "WB_LSMS_ISA", "ISA", source)) %>%
      mutate(source = ifelse(source == "WB_PRIORITY_SURVEY", "PRI_S", source)) %>%
      mutate(source = ifelse(source == "ARAB_LEAGUE_PAPFAM", "PAPFAM", source)) %>%
      mutate(source = ifelse(source == "JHSPH_PERFORMANCE_MONITORING_ACCOUNTABILITY_SURVEY_PMA2020", "PMA", source)) %>% 
      mutate(source = ifelse(nchar(source) > 6, str_trunc(source, 6, ellipsis = ""), source)) %>% #truncate source if it is too long and not specified above
      ungroup() %>% 
      mutate(source = as.factor(source)) %>%
      mutate(point = as.factor(point)) %>%
      data.table()
  }
  
  input_admin0 <- input_clean(input_admin0)
  input_admin1 <- input_clean(input_admin1)
  input_admin2 <- input_clean(input_admin2)
  
  # If subnational NID's are included, assign them values 2 and 3 for point and polygon data that are subnationally representative, respectively
  subnational_nid_subset <- function(df) {
    df %>% 
      mutate(point = as.numeric(levels(point))[point]) %>% 
      rowwise() %>% 
      mutate(point = ifelse(svy_id %in% subnational_nids, point + 2, point)) %>% 
      ungroup() %>% 
      mutate(point = as.factor(point)) %>% 
      data.table()
  }
  
  if (!is.null(subnational_nids)){
    input_admin0 <- subnational_nid_subset(input_admin0)
    input_admin1 <- subnational_nid_subset(input_admin1)
    input_admin2 <- subnational_nid_subset(input_admin2)
  }
  
  # Final list contains each aggregated admin
  input_admins <- list(ad0 = input_admin0, ad1 = input_admin1, ad2 = input_admin2)
  return(input_admins)
  }



