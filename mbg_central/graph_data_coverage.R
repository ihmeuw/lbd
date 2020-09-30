## ###########################################################################
##
## GRAPH DATA COVERAGE
## Purpose: Functions to graph data coverage and completeness for the LBD team
##
## SECTIONS
##   1 - Main function (graph_data_coverage_values)
##   2 - Data prep and validation functions
##   3 - Shapefile prep and validation functions
##   4 - Graphics functions: scatter plots
##   5 - Graphics functions: data coverage maps
##   6 - Graphics functions: Shiny prep
##
## ###########################################################################


## ###########################################################################
## SECTION 1 - MAIN FUNCTION
## ###########################################################################

## graph_data_coverage_values() ----------------------------------------------
#'
#' @title Graph Data Coverage Values
#'
#' @description Produces standard diagnostic maps used to assess input data
#'   coverage across Local Burden of Disease teams.
#'
#' ** ## PARAMETERS ## **
#' ** Basic arguments **
#' @param df The input data.frame or data.table used to create the plot.
#'   Must contain at least the following columns, which will be coerced to the
#'   following data types:
#'   * 'nid' (integer):  NID of survey
#'   * 'country' (char): ISO3 code associated with the survey country
#'   * 'source' (char):  Name of the data source category
#'   * 'year' (integer): Year, either when the data was collected or when an
#'                         event occurred.
#'   * 'N' (numeric):    The sample size at the given location.
#'   * <VAR> (numeric):  The value of the outcome of interest at the survey site.
#'                         Often expressed as a rate
#'   * 'cluster_id' (character): Index assigning unique unique observations from
#'                              the input data. Typically used for debugging.
#'   * 'latitude' (numeric):  Latitude associated with the observation, if available
#'   * 'longitude' (numeric): Longitude associated with the observation, if available
#'   * 'shapefile' (char):    Shapefile associated with the observation. NOTE:
#'                              each observation must have either a valid
#'                              latitude and longitude or a shapefile and location
#'                              code, or else it will be dropped
#'   * 'location_code' (integer): Location code in the given polygon that is
#'                              associated with the observation.
#'   Rows containing NA or empty values in the first seven fields will be
#'   dropped. All rows must have either 'latitude' and 'longitude' fields filled,
#'   or they will be dropped. For more information about how input data is
#'   validated, see the function `dcp_validate_input_data()`.
#' @param var Name of the field in `df` that will be plotted on the maps.
#' @param title Title stub for the map. Full title will take the form
#'   "<title>: <region name>".
#' @param legend_title Title of the legend (color scale) on the map.
#' @param year_var Name of the field in `df` used to indicate which year an
#'   observation falls within. Years not between `year_min` and `year_max` will
#'   be dropped.
#' @param year_min The earliest year of data to include in plots. The standard
#'   for LBD is currently 2000.
#' @param year_max The latest year of data to include in plots. The standard
#'   for LBD is currently 2017.
#' @param region Dictates which area of the world 
#'   will be shown on the plots.
#'
#' ** Arguments controlling code execution and saving **
#' @param cores The number of CPUs available to the function. 
#' @param indicator The indicator being estimated for this data coverage plot.
#'   This argument controls where the data coverage plots will be saved
#' @param extra_file_tag (default `''`) If set to something other than an empty
#'   string, adds additional text to the filepath that would end in
#'   `.../<indicator>/` by default.
#' @param out_dir (default `NULL`) Set a custom save directory. Overrides all
#'   other parameters controlling filepaths above.
#' @param core_repo The directory
#'   used to read in other MBG functions.
#' @param log_dir (default `NULL`) Path to a directory where a log file will
#'   optionally be saved. Useful mainly for debugging.
#'
#' ** Arguments related to speedups and custom outputs **
#' @param fast_shapefiles (default `TRUE`) Load geometries in RDS format rather
#'   than interacting with shapefiles directly. This function should almost
#'   always be set to TRUE, and the RDS files prepared beforehand using the
#'   `synchronize_shapefile_directories()` function in
#'   `"shapefile_functions.R"`.
#' @param new_data_plots (default `FALSE`) Toggles whether to create New Data
#'   plots. Leaving this as the default, `FALSE`, speeds up function execution.
#' @param since_date (default `NULL`) If New Data plots are created, specifies
#'   the date beyond which all data is considered "new". Must be either `NULL`
#'   (in which case the default is the last date when this function was run)
#'   or a character of form `"YYYY-MM-DD"` (including leading zeroes for months
#'   and days).
#' @param annual_period_maps (default `FALSE`) Toggles whether to period maps
#'   for each year between the `year_min` and `year_max` rather than the
#'   standard five-year binned maps. If set to `TRUE`, ONLY annual maps will
#'   be created, and the standard four-period coverage plots will not be
#'   produced. Useful for team-specific data validation.
#' @param save_period_maps (default `TRUE`) Toggles whether or not to save
#'   maps for individual time periods in addition to the standard four-period
#'   data coverage plot. If `annual_period_maps` is `TRUE`, this argument must
#'   also be `TRUE`.
#' @param prep_shiny (default `FALSE`) Toggles whether or not to save data
#'   objects for the interactive Shiny visualization tool.
#' @param return_maps (default `TRUE`) Toggles whether or not to return ggplot
#'   objects for all period-specific maps.
#' @param debug (default `FALSE`) If set to `TRUE`, the function will fail if
#'   any rows of data are dropped due to mismatched types or missing data.
#'
#' ** Arguments controlling map and scatter symbology **
#' @param color_scheme (default `"classic"`)
#' @param color_scheme_scatter (default `"brewer"`)
#' @param high_is_bad (default `TRUE`)
#' @param cap (default `90`)
#' @param cap_type (default `"percentile"`) The type of numeric cap used to set
#'   the legend maximum. Must be one of `'percentile'`, `'absolute'`,
#'   or `'none'`.
#' @param legend_min (default `NA`) Absolute lower bound for the map legend.
#'   Defaults to the lowest observation in the dataset.
#' @param legend_max (default `NA`) Absolute upper bound for the map legend. If
#'   something other than `NA`, overrides the `cap` and `cap_type` arguments.
#' @param endemic_gauls (default `NULL`) This argument can be set as a character
#'   vector containing ISO3 codes. If not `NULL`, all countries not in this
#'   list will be greyed out on the final plot. This argument overrides the
#'   `stage_3_gray` argument.
#' @param stage_3_gray (default `TRUE`) If `TRUE`, greys out all stage 3
#'   countries on the map. Data from these countries will still be displayed.
#' @param simplify_polys (default `TRUE`) If `TRUE`, simplifies input polygons
#'   to speed up map saving. While this step requires an extra initial processing
#'   step, it typically reduces execution time overall.
#' @param tolerance (default `0.03`) If `simplify_polys` is `TRUE`, this argument
#'   controls the degree to which polygons will be simplified.
#' @param base_font_size (default `18`) The default text size used in the plots.
#' @param map_point_size (default `0.8`) The default size of point data on the
#'   map.
#' @param poly_line_width (default `0.2`) The default width of white lines
#'   surrounding polygons on the map. If you are visualizing very high-resolution
#'   polygon data, it might be a good idea to set this argument to `0`.
#' @param remove_rank (default `TRUE`) Testing this parameter - keep as default.
#'
#' @return If `return_maps` is `TRUE`, returns ggplot objects for maps made
#'   for each time period. Otherwise, returns `NULL`.
#'
#' @family Data Coverage Plot functions
#'
graph_data_coverage_values <- function(df,
                                       var,
                                       title,
                                       legend_title,
                                       year_var,
                                       year_min,
                                       year_max,
                                       region,
                                       
                                       cores,
                                       indicator,
                                       extra_file_tag = '',
                                       out_dir = NULL,
                                       core_repo = '<<<< FILEPATH REDACTED >>>>',
                                       log_dir = NULL,
                                       
                                       fast_shapefiles = TRUE,
                                       new_data_plots = FALSE,
                                       since_date = NULL,
                                       annual_period_maps = FALSE,
                                       save_period_maps = TRUE,
                                       prep_shiny = FALSE,
                                       return_maps = TRUE,
                                       debug = FALSE,
                                       
                                       color_scheme = "classic",
                                       color_scheme_scatter = "brewer",
                                       high_is_bad = TRUE,
                                       cap = 90,
                                       cap_type = "percentile",
                                       legend_min = NA,
                                       legend_max = NA,
                                       endemic_gauls = NULL,
                                       stage_3_gray = TRUE,
                                       simplify_polys = TRUE,
                                       tolerance = 0.03,
                                       base_font_size = 18,
                                       map_point_size = 0.8,
                                       poly_line_width = 0.2,
                                       
                                       remove_rank = TRUE
) {
  
  ## I. Load packages & miscellaneous setup ##################################
  message("\n#################################################################")
  message(paste0("** Generating ", indicator, " graphs for ", region, " **\n"))
  
  ## Load packages
  message("Loading packages...")
  # Source package imports function
  source(paste0(core_repo, '/setup.R'))
  # Load all external packages
  package_list <- c('data.table', 'ggplot2', 'parallel', 'doParallel', 'grid',
                    'gridExtra', 'stringr', 'RColorBrewer', 'rgdal', 'sp',
                    'raster', 'magrittr', 'dplyr', 'RMySQL', 'rgeos', 'tidyr')
  load_R_packages(package_list)
  # Load other MBG-specific functions
  source(paste0(core_repo, '/gbd_functions.R'))
  source(paste0(core_repo, '/prep_functions.R'))
  source(paste0(core_repo, '/shapefile_functions.R'))
  
  ## Manually keep order_var space_rank for now
  order_var <- 'space_rank'
  
  ## Define root
  '<<<< FILEPATH REDACTED >>>>'
  
  ## Set up ggplot base font size
  theme_set(theme_minimal(base_size = base_font_size))
  
  #pull in list of stage 3 countries to grey out and remove from scatterplots
  stage3 <- load_gaul_lookup_table()
  stage3 <- stage3[Stage=='3',]
  stage3[, iso3 := toupper(iso3)]
  
  ## II. Prep data ###########################################################
  
  ## 1. Validate input arguments and dataframe -------------------------------
  in_args_all <- as.list(mget(names(formals()),
                              sys.frame(sys.nframe()))
  )
  dcp_validate_input_args(in_args_all)
  
  df_prepped <- dcp_validate_input_data(df,
                                        var      = var,
                                        year_var = year_var,
                                        debug    = debug
  )
  
  ## 2. Pull in country table and merge onto input data ----------------------
  message("Pulling GBD region list...")
  # Pull in list of regions from GBD
  country_table <-suppressMessages(suppressWarnings(
    data.table(get_location_hierarchy(41))[, .(ihme_loc_id,
                                               level,
                                               location_name,
                                               region_name)]
  ))
  # For now, subset to countries only
  #   Note: could change this for subnationals if desired!
  #   May be helpful if looking at a country for which we have tons of data (e.g. India)
  country_table <- country_table[level == 3]
  country_table[, level := NULL]
  setnames(country_table, "ihme_loc_id", "country")
  
  
  ## 2. Load data & split off relevant data sets -----------------------------
  
  # 2. Load data & split off relevant data sets ------------------------------
  message("Preparing data sets...")
  # Create df_poly, df_point, and df_summary
  with_gbd_locs <- dcp_merge_with_gbd_locations(df_prepped, country_table, region)
  df         <- with_gbd_locs[["df"]]                # Full data set
  df_poly    <- with_gbd_locs[["df_poly"]]           # Polygon data only
  df_point   <- with_gbd_locs[["df_point"]]          # Point data only
  df_summary <- with_gbd_locs[["df_summary"]]        # One-row-per-NID summary
  rm(with_gbd_locs)
  
  # If new_data_plots is TRUE, check for new data
  # Save most recent data summary to a CSV
  message("Comparing against data from previous coverage plots...")
  df_summary <- dcp_find_new_data(df_summary, indicator, since_date)
  
  # Grab a country list
  list_output <- get_country_list(df, region)
  reg_title    <- list_output[["reg_title"]]    # Formatted region title
  region_list  <- list_output[["region_list"]]  # Subheadings for regions
  country_list <- list_output[["country_list"]] # ISO3 codes for countries
  rm(list_output)
  
  # Make a summary data set to use for graphing, subset to region
  df_graph <- dcp_make_df_graph(df_summary, country_list)
  
  # Similarly subset df_point and df_poly to just the region of interest
  df_graph_point <- df_point[country %in% country_list]
  df_graph_poly <- df_poly[country %in% country_list]
  
  
  # III. Load polygons #######################################################
  
  # 1. Load polygons for the right (map) side of the figure ------------------
  
  # The polygon map works by pulling shapefiles for each polygon in the
  #   data set, and then coloring those by the lastest year. Polygons are
  #   overlaid chronologically, so the most recent of two overlapping
  #   polygons will be visible.
  
  # Set up point / polygon graphing data sets
  #   will create plots for polygons (latest year of data = color)
  #   and points (color = year of data) separately
  
  # IF annual plots are being made, each individual year will be plotted
  if (annual_period_maps){
    df_graph_poly[, plot_year := year]
    df_graph_point[, plot_year := year]
  } else {
    # Otherwise, bin years to the four standard age bins
    df_graph_poly  <- dcp_bin_years(df_graph_poly)
    df_graph_point <- dcp_bin_years(df_graph_point)
  }
  
  # Only do this bit if there are polygons!
  if(nrow(df_graph_poly) > 0) {
    
    message("Loading individual shapefiles for mapping...")
    
    # Check for missing shapefiles & report if present
    shapefile_list <- dcp_check_for_missing_shapefiles(
      df_graph_poly,
      fast_shapefiles = fast_shapefiles
    )
    shapefiles <- shapefile_list[["shapefiles"]]
    not_real_shapefiles <- shapefile_list[["not_real_shapefiles"]]
    
    if (length(not_real_shapefiles) > 0) {
      warning(paste0("  DROPPING ", length(not_real_shapefiles),
                     " MISSING SHAPES (in shapefile column but not shapefile dir)!"))
      if (!is.null(log_dir)) {
        warning(paste0("  Writing list of missing shapefiles to ", '<<<< FILEPATH REDACTED >>>>',
                       "missing_shapes.csv")
        )
        not_real_shapefiles %>%
          as.data.table %>%
          setnames(., ".", "Missing shapefiles") %>%
          write.csv(., file = '<<<< FILEPATH REDACTED >>>>')
      } else {
        warning(paste0("  Missing shapefiles: ", paste(not_real_shapefiles,
                                                       collapse = ", "))
        )
      }
    }
    
    rm(shapefile_list) # clean up
    
    # Only pull shapefiles if there are any 'real' shapefiles left
    if (length(shapefiles) > 0){
      ## Pull shapes in parallel ---------------------------------------------
      
      # Make a table of shapefiles & associated location codes
      df_shape_loc <- unique(df_graph_poly[, c("shapefile", "location_code")])
      # Pull all polygons in parallel
      poly_list <- pull_polys_in_parallel(shape_loc_list = df_shape_loc,
                                          shapefile_col = "shapefile",
                                          location_code_col = "location_code",
                                          cores = cores,
                                          fast_shapefiles = fast_shapefiles)
      message("Done pulling polys.")
      
      # Find if any broken shapes and, if so, report
      poly_shapes_all <- poly_list[["poly_shapes_all"]]
      broken_shapes <- poly_list[["broken_shapes"]]
      
      if(length(broken_shapes) > 0) {
        message("Warning: the following shapes encountered errors :")
        print(broken_shapes)
      }
      
      if(simplify_polys == T) {
        poly_shapes_all <- simplify_spdf(poly_shapes_all, tol = tolerance)
      }
      rm(poly_list)
    }
  } else {
    
    message("There is no polygon data, so individual shapefiles were not pulled.")
    
    # If no polygons, set products of the above to NULL
    poly_shapes_all <- NULL
    not_real_shapefiles <- NULL
    broken_shapes <- NULL
  }
  
  # 2. Load master shapefile -------------------------------------------------
  
  # This takes a while, so let's do it once only
  if (!("master_shape_all" %in% ls())) {
    
    if (fast_shapefiles == T) {
      message("\nFast-loading master shapefile... ")
      assign("master_shape_all",
             readRDS('<<<< FILEPATH REDACTED >>>>'),
             envir = globalenv())
    } else {
      message("\nOpening master shapefile... (good time to go make a cup of coffee)")
      assign("master_shape_all",
             readOGR(dsn = '<<<< FILEPATH REDACTED >>>>',
                     layer = '<<<< LAYER NAME REDACTED >>>>'),
             envir = globalenv())
    }
  } else {
    message("'master_shape_all' is already in the environment.")
  }
  # Rename ADM0_CODE field in master shapefile
  names(master_shape_all)[names(master_shape_all) == "ADM0_CODE"] <- "GAUL_CODE"
  
  
  # IV. Make scatters and plots ##############################################
  
  # 0. Make sure that all other devices outside of the null device are turned
  #  off, then sink all random graphical output to NULL; function specifies png
  #  where graphs actually desired
  pdf(NULL)   # Failsafe to ensure that dev.off() does not throw an error
  graphics_level <- 2
  while (graphics_level > 1){
    # Keep turning off graphics devices until the null device (1) is reached
    graphics_level <- dev.off()
  }
  pdf(NULL)
  
  # 1. Make table and scatter plots for the left side of the graph -------------
  # Only execute this step if annual maps are not being made, as the annual
  #   map option does not produce the 4-period map!
  if (!annual_period_maps){
    # 1a. Make table for the left side of graph---------------------------------
    table_data <- dcp_make_table_new(df_summary, country_list, year_min, year_max)
    polys_total <- sum(table_data$Polygons)
    points_total <- sum(table_data$Points)
    n_total <- sum(table_data$N)
    # 1b. Make plot for left side of graph -------------------------------------
    message("Making scatter plot...")
    data_scatter_list <- make_data_scatterplots(df_graph   = df_graph,
                                                df_summary = df_summary,
                                                title      = title,
                                                reg_title  = reg_title,
                                                year_min   = year_min,
                                                year_max   = year_max,
                                                table_data = table_data,
                                                base_font_size = base_font_size,
                                                region_name = region,
                                                color_scheme_scatter = color_scheme_scatter,
                                                stage3 = stage3,
                                                stage_3_gray = stage_3_gray,
                                                new_data_plots = new_data_plots)
    g_data     <- data_scatter_list[["g_data"]]                    # scatter for all data
    g_data_new <- data_scatter_list[["g_data_new"]]                # scatter for just new data
    g_data_legend  <- data_scatter_list[["g_data_legend"]]         # legend for g_data
    g_data_new_legend  <- data_scatter_list[["g_data_new_legend"]] # legend for g_data_new
    rm(data_scatter_list)
  }
  
  # 2a: Define map basics --------------------------------------------------
  
  # Create a background map for the plot (just the country polygons)
  message("Constructing background map...")
  background_map_list <- make_background_map(region, endemic_gauls, simplify_polys,
                                             tolerance, fast_shapefiles,
                                             master_shape_all, stage3, stage_3_gray)
  background_outline <- background_map_list[[1]]
  background_map <- background_map_list[[2]]
  background_map_not_endemic <- background_map_list[[3]]
  background_extent <- background_map_list[[4]]
  rm(background_map_list)
  message("  Finished constructing background map.")
  
  color_list <- get_color_list(color_scheme)
  if (high_is_bad==TRUE) color_list <- rev(color_list)
  
  # 2b: Make some graphs ---------------------------------------------------
  if (annual_period_maps){
    # IF annual plots are being made, each year will have its own period map
    map_these_periods <- year_min:year_max
  } else {
    # Otherwise, plot the standard 4 maps
    map_these_periods <- c(2000,2005,2010,2015)
  }
  
  # Create a list to store period objects
  period_map_storage <- vector('list', length=length(map_these_periods))
  
  for(period in map_these_periods) {
    message(paste0("Making period map for ",period,"..."))
    g_map <- make_a_period_map(period, df, region, poly_shapes_all,
                               background_map, background_outline,
                               background_map_not_endemic, background_extent,
                               df_graph_poly, df_graph_point, not_real_shapefiles,
                               color_list, legend_title, log_dir, base_font_size,
                               map_point_size, cap, cap_type, legend_min, legend_max,
                               poly_line_width, annual_period_maps)
    period_map_storage[[as.character(period)]] <- g_map
  }
  
  
  # V. Save plots ###########################################################
  
  if (is.null(out_dir)) out_dir <- '<<<< FILEPATH REDACTED >>>>'
  #dir.create(out_dir, showWarnings = FALSE, recursive = T)
  
  # Make the main plots - 4-up
  # SKIP THIS STEP IF ANNUAL PERIOD MAPS ARE BEING MADE, as the binned period
  #  maps that go into the main plot have not been created.
  if (!annual_period_maps){
    out_file <- '<<<< FILEPATH REDACTED >>>>'
    message(paste0('Saving ', out_file))
    unlink(out_file) # Delete any old versions

    png(filename=out_file,
        units = "in",
        width = 24.33,
        height = 12,
        pointsize = base_font_size,
        res = 1000
    )

    dcp_make_4up_map(g_datamap = g_data,
                     g_data_legend = g_data_legend,
                     map_list = list(period_map_storage[['2000']],
                                     period_map_storage[['2005']],
                                     period_map_storage[['2010']],
                                     period_map_storage[['2015']]),
                     n_countries = length(unique(df_graph$country)),
                     reg_title = reg_title,
                     title = title,
                     base_font_size = base_font_size,
                     n_total = n_total,
                     polys_total = polys_total,
                     points_total = points_total
    )
    dev.off()
    
    if (new_data_plots){
      # Repeat the main plot, but for new data
      out_file <- '<<<< FILEPATH REDACTED >>>>'
      message(paste0('Saving ', out_file))
      unlink(out_file) # Delete any old versions
      
      png(filename=out_file,
          units = "in",
          width = 24.33,
          height = 12,
          pointsize = base_font_size,
          res = 450
      )
      dcp_make_4up_map(g_datamap = g_data_new,
                       g_data_legend = g_data_new_legend,
                       map_list = list(period_map_storage[['2000']],
                                       period_map_storage[['2005']],
                                       period_map_storage[['2010']],
                                       period_map_storage[['2015']]),
                       n_countries = length(unique(df_graph$country)),
                       reg_title = reg_title,
                       title = title,
                       base_font_size = base_font_size,
                       n_total = n_total,
                       polys_total = polys_total,
                       points_total = points_total
      )
      dev.off()
    }
  }
  
  # Make the year bin plots - one for each of four years
  # Save these plots only if anunual period plots are being made or if they have
  #  been specified in the input arguments!
  if (annual_period_maps | save_period_maps){
    for (period in map_these_periods) {
      
      out_file <- paste0(out_dir, 'map_', period, '_', region, '.pdf')
      message(paste0('Saving ', out_file))
      unlink(out_file) # Delete any old versions
      
      # png(filename=out_file,
      #     units = "in",
      #     width = 17,
      #     height = 10,
      #     pointsize = 24,
      #     res = 450)
      pdf(out_file)
      
      print(period_map_storage[[as.character(period)]])
      
      dev.off()
    }
  }
  
  
  # V. Finish up -------------------------------------------------------------
  
  if (prep_shiny) {
    # Save some shiny outputs
    prep_data_coverage_shiny(df, df_graph_poly, poly_shapes_all, var, indicator)
  }
  
  ## Return individual maps if requested
  if(return_maps==TRUE) return(period_map_storage)
  
}




## ###########################################################################
## SECTION 2 - DATA PREP AND VALIDATION FUNCTIONS
## ###########################################################################


#' Validate input arguments (excluding input data.frame)
#'
#' @description: Make sure that all input arguments are valid before beginning
#'   processing.
#'
#' @param input_list The input arguments to graph_data_coverage_values(), stored
#'   as a list.
#'
#' @return Returns NULL (stops execution if any arguments are invalid)
#'
dcp_validate_input_args <- function(input_list){
  message("Validating input arguments...")
  # Evaluate all arguments so you don't get the variables containing the arguments
  input_list <- lapply(input_list, eval)
  # Ensure that all TRUE/FALSE arguments are actually booleans
  boolean_args <- c("high_is_bad","return_maps","fast_shapefiles",
                    "simplify_polys","remove_rank","prep_shiny","new_data_plots",
                    "stage_3_gray","annual_period_maps","save_period_maps")
  for (b_a in boolean_args){
    if ("name" %in% class(input_list[[b_a]])) input_list[[b_a]] <- eval(input_list[[b_a]])
    if (!("logical" %in% class(input_list[[b_a]]))){
      stop(paste0("'",b_a,"' must be either TRUE or FALSE"))
    }
  }
  # If 'annual_period_maps' is true, then 'save_period_maps' must be on
  if (input_list[['annual_period_maps']] & !(input_list[['save_period_maps']])){
    stop(paste0("If 'annual_period_maps' is TRUE, then 'save_period_maps' must",
                " also be TRUE."))
  }
  # Check that all required character fields are actually characters
  char_args <- c("var","title","year_var","region","indicator","color_scheme",
                 "color_scheme_scatter","cap_type","core_repo")
  for (c_a in char_args){
    if(!("character" %in% class(input_list[[c_a]]))){
      stop(paste0("'",c_a,"' must be a of type 'character'."))
    }
  }
  # Check that all required numeric fields are actually numeric
  num_args <- c("year_min","year_max","cores","tolerance","base_font_size",
                "map_point_size","cap")
  for (n_a in num_args){
    if (!('numeric' %in% class(input_list[[n_a]]))){
      stop(paste0("'",n_a,"' must be a of type 'numeric'."))
    }
  }
  # Check that the core_repo path actually exists
  core_repo <- gsub('/$','',input_list[['core_repo']])
  if (!file.exists(core_repo)) stop("Check that the 'core_repo' path is correct.")
  # Check that the 'cap_type' is in a list of valid options
  if (!(input_list[['cap_type']] %in% c('percentile','absolute','none'))){
    stop("'cap_type' must be one of 'percentile','absolute', or 'none'.")
  }
  # Check that if out_dir is not null, then it is a character vector
  out_dir <- input_list[['out_dir']]
  if (!is.null(out_dir)){
    if (!('character' %in% class(out_dir))){
      stop("'out_dir' must be either a character vector or NULL.")
    }
  }
  # Check that the year ranges are somewhat reasonable and year_min < year_max
  year_min <- input_list[['year_min']]
  year_max <- input_list[['year_max']]
  if ((year_min<1900) | (year_min>2018) | (year_max<1900) | (year_max>2018)){
    stop("Check that your year_min and year_max are in a reasonable range...")
  }
  if (year_min > year_max){
    stop("Check that year_min is less than or equal to year_max.")
  }
  # Check that the region is valid
  # This check must be kept up to date with the get_country_list() function!
  valid_regions <- c('africa','africa_no_yem','south_asia',
                     'south_asia_ind_collaborators','se_asia','latin_america',
                     'south_america','south_america_mex','central_america',
                     'central_america_no_mex','eastern_europe','middle_east',
                     'stage2')
  if (!(input_list[['region']] %in% valid_regions)){
    stop(paste0("The mapping 'region' must be one of the following:\n -",
                paste(valid_regions, collapse='\n -')))
  }
  # Check that the color scheme names are valid
  valid_color_schemes <- c('classic','darker_middle','red_blue','carto_red_blue',
                           rownames(brewer.pal.info[brewer.pal.info$category == "seq",])
  )
  if (!(input_list[['color_scheme']] %in% valid_color_schemes)){
    stop(paste0("'color_scheme' must be one of the following:\n -",
                paste(valid_color_schemes, collapse='\n -')))
  }
  valid_scatter_color_schemes <- c('brewer','binary','carto1','carto2','carto3')
  if (!(input_list[['color_scheme_scatter']] %in% valid_scatter_color_schemes)){
    stop(paste0("'color_scheme_scatter' must be one of the following:\n -",
                paste(valid_scatter_color_schemes, collapse='\n -')))
  }
  # Check that 'since_date' is either NULL or in the proper date format
  if (!is.null(input_list[['since_date']])){
    if (!(grepl('[[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2}',
                input_list[['since_date']]))){
      stop(paste0("'since_date' must either be NULL or a character vector with",
                  " the format YYYY-MM-DD"))
    }
  }
  
  # The input arguments have been successfully validated
  message("All non-data input arguments are valid.")
  return(NULL)
}


#' Validate input data.frane
#'
#' Description: Checks that the input data contains the minimum columns needed
#'   to create the data coverage plots, and that those columns are coerced to the
#'   correct data types before continuing. This function also changes the mapped
#'   variable column to a column named 'outcome' in order to simplify future
#'   functions.
#'
#' @param df The input data.frame. Must contain at least the following columns,
#'   which will be coerced to the following data types:
#'   * 'nid' (integer):  This field was called 'svy_id' in the original code
#'   * 'country' (char): ISO3 code associated with the survey country
#'   * 'source' (char):  Name of the data source category
#'   * 'year' (integer): Year, either when the data was collected or when an
#'                         event occurred.
#'   * 'N' (numeric):    The sample size at the given location.
#'   * <VAR> (numeric):  The value of the outcome of interest at the survey site.
#'                         Often expressed as a rate
#'   * 'cluster_id' (integer): Index assigning unique unique observations from
#'                              the input data. Typically used for debugging.
#'   * 'latitude' (numeric):  Latitude associated with the observation, if available
#'   * 'longitude' (numeric): Longitude associated with the observation, if available
#'   * 'shapefile' (char):    Shapefile associated with the observation.
#'   * 'location_code' (integer): Location code in the given polygon that is
#'                              associated with the observation.
#'   Rows containing NA or empty values in the first seven fields will be
#'   dropped. All rows must have either 'latitude' and 'longitude' fields filled,
#'   or they will be dropped.
#' @param var The outcome that is being mapped. Must correspond to a field
#'   in the data.frame
#' @param debug Stops execution of the function rather than dropping any rows
#'   in the input data that do not meet the criteria above. Default FALSE.
#'
#' @return Returns a data.table with the following changes made:
#'   * Converted to a data.table
#'   * Columns coerced to their corresponding type
#'   * Mapped variable field (defined by param 'var') renamed "outcome"
#'   * All extra columns dropped
#'   * Adds a new field, 'pointpoly', indicating whether latitude and longitude
#'       data is available
#'   * Rows not meeting the required completeness criteria dropped
#'
dcp_validate_input_data <- function(df,
                                    var      = 'outcome',
                                    year_var = 'year',
                                    debug    = FALSE
){
  message("Validating input data...")
  # Check that the input data type is a data.frame
  if(!("data.frame" %in% class(df))) stop(paste0("The input data type must be a",
                                                 " data.frame or data.table."))
  # Coerces input data to data.table
  df <- as.data.table(df)
  # Check that all required columns are in the data
  keep_cols <- c('nid', 'country', 'source', 'N', 'cluster_id', 'latitude',
                 'longitude', 'shapefile', 'location_code', year_var, var)
  missing_cols <- keep_cols[!(keep_cols %in% names(df))]
  if (length(missing_cols) > 0) stop(paste0("The input data is missing the ",
                                            "following fields: ",
                                            paste(missing_cols, collapse=', '))
  )
  # Rename the var column to 'outcome'
  if (var != "outcome") df[,outcome := get(var)]
  # Rename the year_var column to 'year'
  if (year_var != "year") df[, year := get(year_var)]
  
  # Subset to only necessary columns
  keep_cols[length(keep_cols)] <- 'outcome'
  keep_cols[length(keep_cols) - 1] <- 'year'
  df <- df[,keep_cols,with=FALSE]
  
  # Helper function for coercing data fields
  coerce_field <- function(field, field_name, coerce_type){
    start_nas <- sum(is.na(field))
    # Convert the field to the specified data type,
    field <- suppressWarnings(eval(parse(text=paste0("as.",coerce_type,'(field)'))))
    # Remove empty character strings, replacing with NA
    if (coerce_type == "character"){
      field[field==""] <- NA
    }
    end_nas <- sum(is.na(field))
    if (end_nas > 0) message(paste0("  ",end_nas," NAs in field '",field_name,
                                    "' (", end_nas - start_nas,
                                    " introduced during coercion to ",
                                    coerce_type,")")
    )
    return(field)
  }
  # Coerce the following fields to integer vectors
  for (col in c('nid','year','location_code')){
    df[,temp := coerce_field(df[,get(col)], col, 'integer')]
    df[,(col) := NULL]
    setnames(df,'temp',col)
  }
  # Coerce the following fields to numeric vectors
  for (col in c('N','outcome','latitude','longitude')){
    df[,temp := coerce_field(df[,get(col)], col, 'numeric')]
    df[,(col) := NULL]
    setnames(df,'temp',col)
  }
  # Coerce the following fields to character vectors
  for (col in c('country','source','cluster_id','shapefile')){
    df[,temp := coerce_field(df[,get(col)], col, 'character')]
    df[,(col) := NULL]
    setnames(df,'temp',col)
  }
  # Helper function defining behavior for dropping rows based on a true-false
  #  series
  informatively_drop <- function(in_df, keep_condition, drop_reason, debug=FALSE){
    start_nrow <- nrow(in_df)
    problem_rows <- unique(in_df[!keep_condition, cluster_id])
    sub_df <- in_df[keep_condition,]
    dropped_nrow <- start_nrow - nrow(sub_df)
    if ( dropped_nrow > 0 ){
      message(paste0("  ",dropped_nrow," rows out of ",start_nrow, " total dropped",
                     " from the input data due to ",drop_reason,"."))
      message(paste0("  Cluster_ids dropped: ",paste(problem_rows, collapse=', ')))
      if ( debug ){
        stop("Validation stopped due to dropped rows ('debug' is on).")
      }
    }
    return(sub_df)
  }
  # Drop any rows containing NA values in any field that does NOT define
  #  geography
  geog_cols <- c('latitude','longitude','shapefile','location_code')
  survey_cols <- keep_cols[!(keep_cols %in% geog_cols)]
  not_missing_survey_cols <- apply(!is.na(df[,survey_cols,with=F]), 1, FUN=all)
  df <- informatively_drop(df,
                           keep_condition = not_missing_survey_cols,
                           drop_reason = paste0("missing identifiers\n   (not ",
                                                "including missing geography ",
                                                "information)"),
                           debug = debug
  )
  # Drop any rows that do not have the required geography data (either 'latitude'
  #  and 'longitude' or 'shapefile' and 'location_code')
  has_latlong <- !is.na(df$latitude) & !is.na(df$longitude)
  has_poly    <- !is.na(df$shapefile) & !is.na(df$location_code)
  has_geo_data <- has_latlong | has_poly
  df <- informatively_drop(df,
                           keep_condition = has_geo_data,
                           drop_reason = "missing geographic data",
                           debug = debug
  )
  final_nrow <- nrow(df)
  # The input data has been cleaned and validated.
  message("  After input validation and cleaning, ",final_nrow," rows of data ",
          "remain.")
  # Add a new field, 'pointpoly', indicating whether the data has latitude and
  #  longitude data available
  df[, pointpoly := "Point"]
  df[, is_point := !is.na(latitude) & !is.na(longitude)]
  df[ !(is_point), pointpoly := "Polygon"]
  df[, is_point := NULL]
  # Formatting for the field "source" - truncate to the first 15 characters
  df[, source := substr(source, 1, 15)]
  return(df)
}


#' Prep input data.table for map-making
#'
#' @description: Merge on country identifiers and get the input data into a
#'   standardized format.
#'
#' @param df Input data.table, after it has been cleaned
#' @param country_table Table of countries (admin0s) from the GBD shared DB
#' @param region Name of the region to be modeled
#'
#' @return Returns a list of four data.tables:
#'  * 'df': The full prepped data, merged with GBD locations
#'  * 'df_summary': The prepped data, collapsed to 1 row per country-year-NID
#'  * 'df_point': The prepped data, subset to point data only
#'  * 'df_poly': The prepped data, subset to poly data only
#'
dcp_merge_with_gbd_locations <- function(df, country_table, region) {
  message("Combining dataframe with GBD country data...")
  # 1. Merge on country identifiers ------------------------------------------
  
  df <- merge(df, country_table, by = "country", all = T)
  
  if (region %in% c("africa","africa_no_yem")) {
    # Rename "North Africa and Middle East" to just "North Africa"
    df <- df[region_name == "North Africa and Middle East",
             region_name := "North Africa"]
  } else if (region == "middle_east") {
    # Rename "North Africa and Middle East" to just "Middle East"
    df <- df[region_name == "North Africa and Middle East",
             region_name := "Middle East"]
  }
  
  # Drop non-matched countries & notify user
  num_na_rows <- nrow(df[is.na(country)])
  if (num_na_rows > 0){
    message(paste0("  Dropping ", nrow(df[is.na(country)]),
                   " rows without matches in GBD country table.")
    )
    message(paste0("  Countries affected: ",
                   paste(unique(df[is.na(country),location_name]), collapse=', '))
    )
    message(paste0("  Cluster_ids affected: ",
                   paste(unique(df[is.na(country),cluster_id]), collapse=', '))
    )
  }
  df <- df[!is.na(country),]
  
  # Truncate long country names
  df <- df[location_name == "Democratic Republic of the Congo", location_name := "DRC"]
  df <- df[location_name == "Central African Republic", location_name := "CAR"]
  df <- df[location_name == "Sao Tome and Principe", location_name := "STP"]
  df <- df[location_name == "United Arab Emirates", location_name := "UAE"]
  df <- df[location_name == "Equatorial Guinea", location_name := "Eq. Guinea"]
  df <- df[location_name == "Saint Vincent and the Grenadines",
           location_name := "St. Vin. & Grenadines"]
  
  
  # 2. Generate subsets of data for further analysis -------------------------
  
  # Split off a point & polygon data set for use later
  df_point <- df[pointpoly == "Point"]
  df_poly  <- df[pointpoly == "Polygon"]
  
  # Sum over the N of the group - includes rows for countries with no data
  df_summary <- df[, .(n = sum(N), count = .N), by = .(source, country, year,
                                                       pointpoly, location_name,
                                                       region_name, nid)]
  return(list("df" = df,
              "df_point" = df_point,
              "df_poly" = df_poly,
              "df_summary" = df_summary))
}


#' Find new data for new data coverage plots
#'
#' @description Take a df_summary object and check to see if there's new data; mark if so
#'
#' @param df_summary The summary prepped data.table object output by the
#'   dcp_merge_with_gbd_locations() function
#' @param indicator The indicator associated with this data, used to associate
#'   with a filename in the folder <<<< FILEPATH REDACTED >>>>
#' @param since_date The date used to compare against old plots, in the format
#'   produced by as.character(Sys.Date()) - ie. YYYY-MM-DD. Defaults to the
#'   last time data was added to the summary table
#'
#' @return The df_summary with a new column, 'new_data', indicating whether
#'   each row was added in the time since 'since_date'
#'
dcp_find_new_data <- function(df_summary, indicator, since_date) {
  
  # Look for an existing data summary table
  summary_dir <- '<<<< FILEPATH REDACTED >>>>'
  summary_file <- '<<<< FILEPATH REDACTED >>>>'
  
  df_summary$date <- as.character(Sys.Date())
  
  if (file.exists(summary_file)) {
    # read in old file
    df_summary_old <- fread(summary_file, stringsAsFactors = F)
    # grab svy ids & dates from old file
    # Rename svy_id to nid, if needed
    needs_rename <- ("svy_id" %in% names(df_summary_old)) & !("nid" %in% names(df_summary_old))
    if (needs_rename) df_summary_old[, nid := svy_id]
    # Only merge on the old df_summary if both merge columns exist
    if (("nid" %in% names(df_summary_old)) & ("date" %in% names(df_summary_old))){
      # Ensure that the nid_dates_old is unique by NID
      nid_dates_old <- unique(df_summary_old[, c("nid", "date")])
      nid_dates_old <- nid_dates_old[!is.na(nid)]
      nid_dates_old <- nid_dates_old[!duplicated(nid),]
      
      # Ensure that the nid and date fields are the correct data types
      nid_dates_old[, date_old := as.character(date)]
      suppressWarnings(nid_dates_old[, temp := as.integer(nid)])
      nid_dates_old[, nid := NULL]
      setnames(nid_dates_old, 'temp', 'nid')
      nid_dates_old <- nid_dates_old[, c("nid", "date_old"), with=F]
      
      # merge
      df_summary <- merge(df_summary,
                          nid_dates_old,
                          by = c('nid'),
                          all.x = T)
      
      # replace if an older date exists
      df_summary[(date != date_old) & !is.na(date_old), date := date_old]
      df_summary[, date_old := NULL]
    }
  }
  
  #replace dates for country-rows with no data with "na"
  df_summary[is.na(nid), date := NA]
  df_summary <- df_summary[order(location_name)]
  
  unlink(summary_file)
  write.csv(df_summary, file = summary_file, row.names=FALSE)
  
  # Mark which data is new with a variable "new_data" (for graphing)
  if (is.null(since_date)) {
    df_summary[date == max(df_summary$date, na.rm = T), new_data := 1]
    df_summary[!is.na(date) & is.na(new_data), new_data := 0]
  } else {
    df_summary[(as.Date(date) > as.Date(since_date)) & !is.na(date), new_data := 1]
    df_summary[!is.na(date) & is.na(new_data), new_data := 0]
  }
  
  return(df_summary)
  
}


#' Get list of countries in modelling region
#'
#' @description Pull a country list based on your data & region of choice
#'
#' @param df The dataframe output by the dcp_merge_with_gbd_locations() function
#' @param region Name of the region to be modeled
#'
#' @return returns a list with the following 3 objects -
#'  * 'reg_title': A character string of the region
#'  * 'region_list': A vector of the modelling regions associated with the region passed in
#'  * 'country_list': A vector of country ISO3 codes within the region
#'
get_country_list <- function(df, region) {
  
  # Note: need to figure out how to deal with Oceania and Central Asia
  
  if (region == 'africa') {
    reg_title <- "Africa"
    region_list <- c("North Africa",
                     "Central Sub-Saharan Africa",
                     "Eastern Sub-Saharan Africa",
                     "Western Sub-Saharan Africa",
                     "Southern Sub-Saharan Africa")
    country_list <- unique(df[region_name %in% region_list]$country)
    country_list <- unique(c(country_list, "YEM"))
    
    # This is tough because of NAME including middle east - need to manually remove and add some
    remove_countries <- c("AFG", "ARE", "IRN", "IRQ", "JOR", "OMN", "PSE", "SAU", "SYR", "TUR","KWT","LBN","QAT","BHR","CPV")
    
    country_list <- country_list[!(country_list %in% remove_countries)]
  }
  
  if (region == 'africa_no_yem') {
    reg_title <- "Africa"
    region_list <- c("North Africa",
                     "Central Sub-Saharan Africa",
                     "Eastern Sub-Saharan Africa",
                     "Western Sub-Saharan Africa",
                     "Southern Sub-Saharan Africa")
    country_list <- unique(df[region_name %in% region_list]$country)
    # This is tough because of NAME including middle east - need to manually remove and add some
    remove_countries <- c("AFG","ARE","IRN","IRQ","JOR","OMN","PSE","SAU", 
                          "SYR","TUR","KWT","LBN","QAT","BHR","CPV","YEM")
    country_list <- country_list[!(country_list %in% remove_countries)]
  }
  
  if (region %in% c('south_asia','south_asia_ind_collaborators')) {
    reg_title <- "South Asia"
    region_list <- c("South Asia")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    # add Sri Lanka
    country_list <- unique(c(country_list, "LKA"))
  }
  
  if (region == 'se_asia') {
    reg_title <- "Southeast Asia"
    region_list <- c("East Asia",
                     "Southeast Asia")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    # move Sri Lanka & the Maldives to south_asia
    remove_countries <- c("LKA", "MDV")
    country_list <- country_list[!(country_list %in% remove_countries)]
    
    # add PNG
    country_list <- unique(c(country_list, "PNG", "MNG"))
  }
  
  if (region == 'latin_america') {
    reg_title <- "Latin America and Caribbean"
    region_list <- c("Andean Latin America",
                     "Caribbean",
                     "Central Latin America",
                     "Tropical Latin America")
    country_list <- unique(df[region_name %in% region_list]$country)
    country_list <- unique(c(country_list, "CUB"))
  }
  
  if (region == 'south_america') {
    reg_title <- "South America"
    region_list <- c("Andean Latin America",
                     "Tropical Latin America")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    # add Venezuela & others
    add_countries <- remove_countries <- c("VEN", "COL", "GUY", "SUR")
    country_list <- unique(c(country_list, add_countries))
  }
  
  if (region == 'south_america_mex') {
    reg_title <- "South America"
    region_list <- c("Andean Latin America",
                     "Tropical Latin America")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    # add Venezuela & others
    add_countries <- remove_countries <- c("VEN", "COL", "GUY", "SUR", "MEX")
    country_list <- unique(c(country_list, add_countries))
  }
  
  if (region == 'central_america') {
    reg_title <- "Central Latin America and Carribean"
    region_list <- c("Central Latin America",
                     "Caribbean")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    # remove Venezuela & others
    remove_countries <- c("VEN", "COL", "GUY", "SUR")
    country_list <- country_list[!(country_list %in% remove_countries)]
  }
  
  if (region == 'central_america_no_mex') {
    reg_title <- "Central Latin America and Carribean"
    region_list <- c("Central Latin America",
                     "Caribbean")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    # remove Venezuela & others
    remove_countries <- c("VEN", "COL", "GUY", "SUR", "MEX")
    country_list <- country_list[!(country_list %in% remove_countries)]
  }
  
  if (region == 'eastern_europe') {
    reg_title <- "Eastern Europe"
    region_list <- c("Eastern Europe")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    remove_countries <- c("RUS")
    country_list <- country_list[!(country_list %in% remove_countries)]
  }
  
  if (region == 'middle_east') {
    reg_title <- "Middle East and Central Asia"
    region_list <- c("Middle East")
    country_list <- unique(df[region_name %in% region_list]$country)
    
    remove_countries <- c("EGY","SDN","TUN","DZA","LBY","MAR")
    country_list <- country_list[!(country_list %in% remove_countries)]
    
    # Add in several Central Asia countries:
    # Uzbekistan, Turkmenistan, Tajikistan, Kyrgyzstan, Lebanon
    country_list <- unique(c(country_list, "UZB", "TKM", "TJK", "KGZ", "LBN"))
    
  }
  
  if (region == 'stage2'){
    reg_title <- 'Low and Middle-Income Countries'
    # Load stage metadata and subset to stage 2
    lookup_table <- load_gaul_lookup_table()
    stage2 <- lookup_table[Stage != '3']
    # Get all unique regions and countries in stage 2
    region_list <- unique(stage2[,reg_name])
    country_list <- toupper( unique(stage2[,iso3]) )
  }
  
  return(list("reg_title" = reg_title,
              "region_list" = region_list,
              "country_list" = country_list))
}


#' Make df_graph
#'
#' @description Make the df_graph object - subsetted to countries in country_list & formatted
#'
#' @param df_summary The summary prepped data.table object output by the
#'   dcp_merge_with_gbd_locations() function
#' @param country_list Vector of ISO3 codes from get_country_list() function
#'
#' @return The df_summary df subsetted to the countries in country list and
#'   formatted to handle na values in geospatial variables and new_data
#'
dcp_make_df_graph <- function(df_summary, country_list) {
  df_graph <- df_summary[country %in% country_list]
  
  # Fix NAs so that they don't graph
  message(paste0("\nFound ", nrow(df_graph[is.na(pointpoly)]), " rows not designated by point or poly."))
  message("Typically, this indicates that there are countries with no data - check this assumption if a large number")
  message("Fixing to ensure no NAs in legend...")
  df_graph[is.na(pointpoly), n := 0]
  df_graph[is.na(pointpoly), source := unique(df_graph$source[!is.na(df_graph$source)])[1]]
  df_graph[is.na(pointpoly), pointpoly := "Point"]
  
  df_graph$location_name <- factor(df_graph$location_name,
                                   levels = rev(sort(unique(df_graph$location_name))))
  
  df_graph[new_data == 0, new_data_lab := "No"]
  df_graph[new_data == 1, new_data_lab := "Yes"]
  df_graph[is.na(new_data), new_data_lab := "No"]
  
  return(df_graph)
}


#' Bin Years
#'
#' @description Take a data table and bin the year variable
#'
#' @param df_to_bin The subsetted and formatted data.table object output by the
#'   dcp_make_df_graph() function split into points or polygons
#'
#' @return A data.table with binned plot_year variable
#'
dcp_bin_years <- function(df_to_bin) {
  df_to_bin <- df_to_bin[, survey := paste0(source,'_',country,'_',year)]
  df_to_bin <- subset(df_to_bin, year > 1999)
  df_to_bin <- df_to_bin[year > 1999 & year <= 2002, plot_year := 2000]
  df_to_bin <- df_to_bin[year > 2002 & year <= 2007, plot_year := 2005]
  df_to_bin <- df_to_bin[year > 2007 & year <= 2012, plot_year := 2010]
  df_to_bin <- df_to_bin[year > 2012, plot_year := 2015]
  
  return(df_to_bin)
}




## ###########################################################################
## SECTION 3 - SHAPEFILE PREP AND VALIDATION FUNCTIONS
## ###########################################################################

#' Simplify SpatialPolygonsDataFrame
#'
#' @description Reduce the number of vertices in a polygon,
#'   wrapper function for gSimplify.
#'
#' @param spdf A SpatialPolygonsDataFrame to be simplified
#' @param tol The tolerance which determines the degree of simplification
#'
#' @return A simplified SpatialPolygonsDataFrame
#'
simplify_spdf <- function(spdf, tol = tolerance) {
  df_spdf <- data.frame(spdf)
  spdf <- gSimplify(spdf, tol = tol, topologyPreserve = T)
  spdf <- SpatialPolygonsDataFrame(spdf, df_spdf)
  return(spdf)
}


#' Pull polygons for a single shapefile
#'
#' @description Pulls the polygons from a master shapefile based on input data.table.
#'
#' @param shape_loc A data.table with columns 'shapefile' and 'location_code', where
#'   'location_code' corresponds to the GAUL_CODE column of corresponding 'shapefile'
#' @param fast_shapefiles Boolean. If true, reads in shapefile from RDS, else reads 
#'   from shapefile directory
#'
#' @return A SpatialPolygonsDataFrame with polygons specified in shape_loc
#'
pull_polys <- function(shape_loc, fast_shapefiles) {
  # Function to pull in a shapefile and assign each polygon the latest year
  # that any data was collected within it (determines fill color)
  
  # This is computationally intensive if many shapefiles
  
  shape <- unique(shape_loc$shapefile)
  loc_codes <- unique(shape_loc$location_code)
  
  message(paste0("  Working on ", shape))
  
  # Read in the master shapefile
  if (fast_shapefiles == T) {
    master_shape <- fast_load_shapefile(shape)
  } else  {
    master_shape <- readOGR(dsn = '<<<< FILEPATH REDACTED >>>>',
                            layer = shape)
  }
  
  message(paste0("    CRS: ",crs(master_shape)))
  
  names(master_shape)[names(master_shape)=="GAUL_Code"] <- "GAUL_CODE"
  
  # Subsetting will break if NAs in row index (GAUL_CODE)
  master_shape <- master_shape[!is.na(master_shape$GAUL_CODE),]
  
  # Subset to the relevant data & shapefile bits (by gaul code)
  subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% loc_codes, ]
  
  # Remove the data that we don't need & make a standard output format
  return_shape <- subset_shape[names(subset_shape) %in% c("GAUL_CODE")]
  
  # Remove any duplicates
  return_shape <- return_shape[duplicated(return_shape$GAUL_CODE) == F, ]
  
  # Convert GAUL_CODE to numeric
  return_shape@data$GAUL_CODE <- as.numeric(as.character(return_shape@data$GAUL_CODE))
  
  return(return_shape)
  
}


#' Pull polygons in parallel
#'
#' @description Pulls the polygons from a master shapefile based on input data.table.
#'   Calls `pull_polys()` in parallel and combines shapefiles into a single SPDF.
#'
#' @param shape_loc_list A data.table with columns 'shapefile' and 'location_code', where
#'   'location_code' corresponds to the GAUL_CODE column of corresponding 'shapefile'
#' @param shapefile_col name of shapefile column in `shape_loc_list`
#' @param location_code_col name of location_code column in `shape_loc_list`
#' @param cores number of cores to use in parallelization
#'   (see `graph_data_coverage_values()` for more information)
#' @param fast_shapefiles Boolean. If true, reads in shapefile from RDS,
#'   else reads from shapefile directory
#'
#'
#' @return A list with 2 objects -
#'   * 'poly_shapes_all': A combined SPDF with all polygons combined
#'   * 'broken_shapes': A list of shapefiles that did not load correctly
#'
pull_polys_in_parallel <- function(shape_loc_list,
                                   shapefile_col = "shapefile",
                                   location_code_col = "location_code",
                                   cores,
                                   fast_shapefiles = fast_shapefiles) {
  
  
  # Format table & rename to standard format
  shape_loc_list <- shape_loc_list %>%
    as.data.table %>%
    unique %>%
    setnames(., c(shapefile_col, location_code_col),
             c("shapefile", "location_code"))
  
  # Check for NAs
  na_rows <- shape_loc_list[is.na(shapefile) | is.na(location_code)]
  if (nrow(na_rows) > 0) {
    message("WARNING: You have NAs in your input data.  Check original data source & fix, as these will be dropped!")
    message("Affected rows:")
    print(na_rows)
    warning(paste0("NAs found in input data (", nrow(na_rows), " rows affected). These are being dropped! Check your input data & original source!"))
  }
  
  # Drop NAs
  shape_loc_list <- shape_loc_list[!is.na(shapefile) & !is.na(location_code)]
  
  # Get shapefiles
  shapefiles <- unique(shape_loc_list$shapefile)
  
  # Set up cluster
  cores <- min(c(cores, length(shapefiles))) # just get enough cores for shapefiles
  
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  
  if (fast_shapefiles == T) {
    message("Fast-pulling polys...")
    poly_shape_list <- lapply(1:length(shapefiles), function(i) {
      shape <- shapefiles[i]
      shape_loc <- shape_loc_list[shapefile == shape]
      pull_polys(shape_loc, fast_shapefiles = T)})
  } else {
    message("Pulling polys in parallel. This may take a while ...")
    # Distribute packages
    clusterCall(cl, function(x) {
      .libPaths('<<<< FILEPATH REDACTED >>>>')
      library(rgdal)
      library(data.table)
    })
    
    poly_shape_list <- foreach (i = 1:length(shapefiles),
                                .export=c("j_root", "pull_polys", "fast_load_shapefile"),
                                .errorhandling = 'pass') %dopar% {
                                  shape <- shapefiles[i]
                                  shape_loc <- shape_loc_list[shapefile == shape]
                                  pull_polys(shape_loc, fast_shapefiles = FALSE)
                                }
    
    stopCluster(cl)
  }
  
  # Find broken shapefiles
  find_broken_shapes <- function(shape) {
    if ('error' %in% class(shape)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  names(poly_shape_list) <- shapefiles
  
  broken_shape_names <- names(poly_shape_list)[unlist(lapply(poly_shape_list, find_broken_shapes))]
  working_shape_names <- names(poly_shape_list)[!(names(poly_shape_list) %in% broken_shape_names)]
  
  broken_shape_list <- lapply(broken_shape_names, function(x) poly_shape_list[[x]])
  names(broken_shape_list) <- broken_shape_names
  working_shape_list <- lapply(working_shape_names, function(x) poly_shape_list[[x]])
  names(working_shape_list) <- working_shape_names
  
  # Remake IDs to be unique across all of poly_shape_list (needed for rbind below)
  makeID_shiny <- function(shape) {
    shapefile <- poly_shape_list[[shape]]
    shapefile <- spChFIDs(shapefile, paste0(shape, "_", row.names(shapefile@data)))
    return(shapefile)
  }
  
  working_shape_list <- lapply(names(working_shape_list), makeID_shiny)
  
  ## check to make sure that all the projections are the same
  proj.strings <- unlist(lapply(working_shape_list, proj4string))
  if(length(unique(proj.strings)) > 1){
    ## take the most common one and assign it to the others
    tab.strings <- table(proj.strings)
    common.proj <- names(tab.strings[which(tab.strings == max(tab.strings))])
    diff.strings.idx <- which(proj.strings != common.proj)
    for(ss in diff.strings.idx){
      message(sprintf("changing the projection on shapefile: %s to match the majority of the other proj strings", shapefiles[ss]))
      message(sprintf("--- it will be changed from: %s" , proj4string(working_shape_list[[ss]])))
      message(sprintf("--- to: %s", common.proj))
      message(sprintf("if those look like substantially different projections you should go back and fix the shapefile"))
      proj4string(working_shape_list[[ss]]) <- common.proj
    }
  }
  
  ## Combine all into one big SPDF
  poly_shapes_all <- do.call(rbind, working_shape_list)
  rm(working_shape_list)
  
  # Pull out the shapefile name and assign this as a new field
  poly_shapes_all@data$id <- rownames(poly_shapes_all@data)
  poly_shapes_all$shapefile <- gsub("_[^_]+$", "", poly_shapes_all$id)
  
  poly_shapes_all$location_code <- as.numeric(poly_shapes_all$GAUL_CODE)
  
  return(list("poly_shapes_all" = poly_shapes_all, "broken_shapes" = broken_shape_list))
  
}

#' Check for missing shapefiles
#'
#' @description Checks that all shapefiles listed in the df_graph_poly data.table exist
#'   in the location they are being pulled from (determined by `fast_shapefiles`).
#'
#' @param df_graph_poly A data.table that only has data matched to polygons
#' @param fast_shapefiles Boolean. If true, reads in shapefile from RDS,
#'   else reads from shapefile directory using readOGR
#'
#' @return A list with 2 objects -
#'   * 'shapefiles': A list of shapefiles that exist in directory
#'   * 'not_real_shapefiles': A list of shapefiles that do not exist in directory
#'
dcp_check_for_missing_shapefiles <- function(df_graph_poly, fast_shapefiles) {
  # Pull a list of shapefiles for the polygon data for the region in question
  shapefiles <- unique(df_graph_poly$shapefile) %>% as.character %>% tolower
  
  ## Check that all shapefile entries are real shapefiles and report bad entries
  if (fast_shapefiles){
    real_shapefiles <- gsub('.rds','',
                            tolower(list.files('<<<< FILEPATH REDACTED >>>>',
                                               pattern = '.rds'))
    )
  } else {
    real_shapefiles <- gsub('.shp','',
                            tolower(list.files('<<<< FILEPATH REDACTED >>>>',
                                               pattern = '.shp'))
    )
  }
  not_real_shapefiles <- shapefiles[!(shapefiles %in% real_shapefiles)]
  
  shapefiles <- shapefiles[(shapefiles %in% real_shapefiles)]
  
  return(list("shapefiles" = shapefiles,
              "not_real_shapefiles" = not_real_shapefiles))
}




## ###########################################################################
## SECTION 4 - GRAPHICS FUNCTIONS: SCATTER PLOTS
## ###########################################################################

#' Make Table
#'
#' @description Makes the table to the right of the data scatterplot with the counts of
#'   points and polygons in each country
#'
#' @param df_summary The df_summary data.table outputted by `dcp_merge_with_gbd_locations()`
#' @param country_list list of ISO3s outputted by `get_country_list()`
#' @param year_min earliest year, subsets the data to the minimum year
#' @param year_max latest year, subsets the data to the maximum year
#'
#' @return A data.table with the number of cumulative points and polygons in the
#'   year range by country
#'
dcp_make_table_new <- function(df_summary, country_list, year_min, year_max) {
  td <- as.data.table(df_summary)
  td <- td[(year >= year_min) &
             (year <= year_max) &
             (country %in% country_list),]
  td <- td[, c("country", "pointpoly", "n", "count"), with = F]
  td <- td[,Count := sum(count), by=.(country, pointpoly)]
  td <- td[,N := sum(n), by=.(country, pointpoly)]
  td <- distinct(td, country, pointpoly, N, Count)
  
  td[is.na(pointpoly), N := 0]
  td[is.na(pointpoly), Count := 0]
  
  td_n <- td[,.(country, N)]
  td_n <- td[, N := sum(N), by=country]
  td_n <- distinct(td, country, N)
  
  td_count <- td %>%
    .[, N := NULL] %>%
    spread(pointpoly, Count, fill = 0) %>%
    as.data.table %>%
    merge(., td_n, by = "country") %>%
    setnames("country", "Country")
  
  # Catch if either no point or no poly data
  if (!("Point" %in% names(td_count))) {
    td_count[, Point := 0]
  }
  if (!("Polygon" %in% names(td_count))) {
    td_count[, Polygon := 0]
  }
  
  setnames(td_count, c("Point", "Polygon"), c("Points", "Polygons"))
  
  # Add in countries with no data
  no_data_countries <- country_list[!(country_list %in% unique(td_count$Country))]
  ## Catch the cases where all countries are filled
  if (length(no_data_countries) > 0){
    no_data_td <- data.table(Country = no_data_countries,
                             Points = 0,
                             Polygons = 0,
                             N = 0)
    
    td_count <- rbind(td_count, no_data_td)
  }
  
  if ("<NA>" %in% names(td)) td[, "<NA>" := NULL]
  
  setcolorder(td_count, c("Country", "Points", "Polygons", "N"))
  
  return(td_count)
}


#' Extract Legend
#'
#' @description Extracts the legend from a 5-yr bin graph for final plot
#'
#' @param a.plot a ggplot or grob object
#'
#' @return A grob legend object
#'
gLegend<-function(a.plot){
  
  if ("ggplot" %in% class(a.plot)) {
    tmp <- ggplot_gtable(ggplot_build(a.plot))
  } else if ("grob" %in% class(a.plot)) {
    tmp <- .gplot
  }
  
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


#' Add color scheme to scatterplot
#'
#' @description Defines the color scheme used in the data scatterplot
#'
#' @param g_datamap a ggplot data scatter built in `build_g_data_scatter()`
#' @param color_scheme_scatter a character string specifying different color schemes
#'   options: `brewer`, `binary`, `carto1`, `carto2`, `carto3`
#'
#' @return g_datamap with a color scale applied
#'
add_color_scheme_scatter <- function(g_datamap, color_scheme_scatter) {
  if(color_scheme_scatter == "brewer") {
    g_datamap <- g_datamap + scale_color_brewer(palette = "Paired")
  } else if (color_scheme_scatter == "binary") {
    g_datamap <- g_datamap + scale_colour_manual(values = c("No" = "gray", "Yes" = "firebrick"))
  } else if (color_scheme_scatter == "carto1") {
    g_datamap <- g_datamap + scale_colour_manual(values = c("#5F4690","#1D6996","#38A6A5","#0F8554",
                                                            "#73AF48","#EDAD08","#E17C05","#CC503E",
                                                            "#94346E","#6F4070","#994E95","#666666"))
  } else if (color_scheme_scatter == "carto2") {
    g_datamap <- g_datamap + scale_colour_manual(values = c("#7F3C8D","#11A579","#3969AC","#F2B701",
                                                            "#E73F74","#80BA5A","#E68310","#008695",
                                                            "#CF1C90","#f97b72","#4b4b8f","#A5AA99"))
  } else if (color_scheme_scatter == "carto3") {
    g_datamap <- g_datamap + scale_colour_manual(values = c("#E58606","#5D69B1","#52BCA3","#99C945",
                                                            "#CC61B0","#24796C","#DAA51B","#2F8AC4",
                                                            "#764E9F","#ED645A","#CC3A8E","#A5AA99"))
  }
}


#' Build data scatterplot
#'
#' @description Builds the data scatterplot and underlying table on the left hand side of coverage plots
#'
#' @param df_graph data.table outputted by `dcp_make_df_graph()`
#' @param alpha_val number between 0 and 1 controlling transparency of points
#' @param by_color name of column in df_graph defining color groupings in plot
#' @param color_label title of color portion of the scatterplot legend
#' @param color_scheme_scatter a character string specifying different color schemes
#'   options: `brewer`, `binary`, `carto1`, `carto2`, `carto3`
#' @param size_lab title of the size portion of the scatterplot legend
#' @param title title of the data coverage plot
#' @param reg_title name of the region being mapped
#' @param region_name column in df_graph that has the region names for the plot
#' @param table_data data.table of point and polygon counts by country
#'   from `dcp_make_table_new()`
#' @param base_font_size font size. All font sizes in the plots are scaled off this
#' @param year_min minimum year, controls range of years shown in scatterplot
#' @param year_max maximum year, controls range of years shown in scatterplot
#' @param stage3 data.table of stage 3 countries to be removed
#' @param stage_3_gray boolean. If true, removes stage 3 countries from scatterplot
#'
#' @return returns a list of 2 objects -
#'   * 'g_datamap': scatterplot ggplot object
#'   * 'g_table': underlying table ggplot object
#'
build_g_data_scatter <- function(df_graph,
                                 alpha_val,
                                 by_color,
                                 color_label,
                                 color_scheme_scatter,
                                 size_lab,
                                 title,
                                 reg_title,
                                 region_name,
                                 table_data,
                                 base_font_size,
                                 year_min,
                                 year_max,
                                 stage3,
                                 stage_3_gray) {
  
  # Set up table data
  td <- copy(table_data)
  td[, N := NULL] # Don't display N column
  setnames(td, "Country", "country")
  td <- gather(td, type, value, -country) %>% as.data.table
  
  # Add on regions
  reg_table <- subset(df_graph, select = c("country", "location_name", "region_name")) %>%
    unique
  td <- merge(td, reg_table, by="country")
  
  if(stage_3_gray){
    td <- td[!(country %in% stage3$iso3),]
    df_graph <- df_graph[!(country %in% stage3$iso3),]
  }
  
  base_font_size <- round(base_font_size * 0.8)
  
  panel_spacing <- 2
  if (region_name %in% c("africa","africa_no_yem")) panel_spacing <- 1
  if (region_name == 'stage2') panel_spacing <- .5
  
  
  # Make the table part of the plot
  g_table <- ggplot(data = td,
                    aes(x = type,
                        y = location_name,
                        label = value)) +
    geom_text(size = base_font_size*(5/14)*(0.8)) +
    theme_minimal(base_size = base_font_size) + #
    scale_x_discrete(position = "top") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          strip.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = rel(1)),
          panel.spacing.y = unit(panel_spacing, "lines")) +
    facet_wrap(~region_name, scales="free_y", ncol = 1)
  
  # Builds a data scatterplot
  g_datamap <- ggplot(data = df_graph,
                      aes(x = year,
                          y = location_name,
                          group = get(by_color))) +
    geom_vline(xintercept = c(2000,2005,2010,2015)) +
    geom_point(aes(size = capped_n,
                   shape = pointpoly,
                   color = get(by_color)),
               alpha = alpha_val) +
    scale_size_area(limits = c(0, max(df_graph$capped_n[!is.na(df_graph$capped_n)])), max_size = 4) +
    theme_minimal(base_size = base_font_size) +
    xlim(year_min, year_max) +
    labs(shape = "Data Type",
         color = color_label,
         size = size_lab) +
    guides(color = guide_legend(order = 1, override.aes = list(shape = 15, size = 5)),
           shape = guide_legend(order = 2, override.aes = list(size = 4)),
           size = guide_legend(order = 3)) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          strip.text.x = element_text(size = rel(1.3),
                                      vjust = 0.5,
                                      margin = margin(b = 20)),
          panel.spacing.y = unit(panel_spacing, "lines"),
          plot.margin = margin(l = 12)) +
    facet_wrap(~region_name, scales="free_y", ncol = 1)
  
  g_datamap <- add_color_scheme_scatter(g_datamap, color_scheme_scatter)
  
  return(list(g_datamap, g_table))
}


#' Data Scatter Height Fix
#'
#' @description Fixes the height of the rows in the data scatterplot and combines
#'   underlying table and data scatter into a grob object
#'
#' @param g_data_scatter_list list outputted by `build_g_data_scatter()`
#' @param df_graph data.table outputted by `dcp_make_df_graph()`
#' @param stage3 data.table of stage 3 countries to be removed
#' @param stage_3_gray boolean. If true, removes stage 3 countries from scatterplot
#'
#' @return returns a grob object with combined table and scatterplot
#'
fix_g_data_scatter <- function(g_data_scatter_list, df_graph, stage3, stage_3_gray) {
  
  # This function takes the data map objects and makes it
  # so that all of the heights are evenly distributed
  # returns a grob object instead of a ggplot object
  if(stage_3_gray) {
    df_graph <- df_graph[!(country %in% stage3$iso3),]
  }
  
  height_fix <- function(g_data_scatter) {
    g_data_scatter <- ggplotGrob(g_data_scatter)
    
    # Figure out # rows per region
    regions_table <- unique(df_graph[, c("country", "region_name")])
    regions_table <- as.data.frame(table(regions_table$region_name))
    
    names(regions_table) <- c("region_name", "n")
    
    # This is a bit weird
    # identify the indices of the relevant heights from the grob layout
    idx <- g_data_scatter$layout$b[grepl("panel", g_data_scatter$layout$name)]
    
    # reset heights (relative heights, scaled to # rows)
    g_data_scatter$heights[idx] <- unit(regions_table$n, "null")
    
    return(g_data_scatter)
  }
  
  # Fix heights
  g_data_fixed_list <- lapply(g_data_scatter_list, height_fix)
  
  g_data_combined <- cbind(g_data_fixed_list[[1]], g_data_fixed_list[[2]], size = "last")
  # Get index of the 2nd set of panels
  idx <- g_data_combined$layout$l[grepl("panel", g_data_combined$layout$name)] %>% unique
  g_data_combined$widths[idx[1]] <- unit(4, "null")
  g_data_combined$widths[idx[2]] <- unit(1.2, "null")
  
  return(g_data_combined)
  
}


#' Build data scatterplot wrapper
#'
#' @description Wrapper for build_g_data_scatter and fix_g_data_scatter to create
#'   complete data scatterplot grob objects for all and new data
#'
#' @param df_graph data.table outputted by `dcp_make_df_graph()`
#' @param df_summary data.table outputted by `dcp_make_df_graph()`
#' @param title title of the data coverage plot
#' @param reg_title name of the region being mapped
#' @param year_min minimum year, controls range of years shown in scatterplot
#' @param year_max maximum year, controls range of years shown in scatterplot
#' @param table_data data.table of point and polygon counts by country
#'   from `dcp_make_table_new()`
#' @param base_font_size font size. All font sizes in the plots are scaled off this
#' @param region_name column in df_graph that has the region names for the plot
#' @param color_scheme_scatter a character string specifying different color schemes
#'   options: `brewer`, `binary`, `carto1`, `carto2`, `carto3`
#' @param stage3 data.table of stage 3 countries to be removed
#' @param stage_3_gray boolean. If true, removes stage 3 countries from scatterplot
#' @param new_data_plots boolean. If true, makes data scatterplots with data marked new from
#'   `dcp_find_new_data()`
#'
#' @return returns a list of 4 objects -
#'   * 'g_data_legend': grob object, legend of g_data
#'   * 'g_data_new_legend': grob object, legend of g_data_new
#'   * 'g_data': scatterplot grob object with all data
#'   * 'g_data_new': scatterplot grob object with new data
#'
make_data_scatterplots <- function(df_graph,
                                   df_summary,
                                   title,
                                   reg_title,
                                   year_min,
                                   year_max,
                                   table_data,
                                   base_font_size,
                                   region_name,
                                   color_scheme_scatter,
                                   stage3,
                                   stage_3_gray,
                                   new_data_plots = FALSE) {
  
  # Fix title if present
  if(title != "") {
    title <- paste0(title, " ")
  }
  
  # Set up heading for size (proportional to 'N')
  size_lab <- 'Sample Size'
  
  # Ensure that year_min, year_max are numeric
  year_min <- as.numeric(year_min)
  year_max <- as.numeric(year_max)
  
  # Cap shape size at percentile
  n_cap_pctile <- 0.9
  df_graph <- df_graph[, capped_n := n]
  df_graph <- df_graph[capped_n >= quantile(df_summary$n[!is.na(df_summary$n)], probs = n_cap_pctile),
                       capped_n := quantile(df_summary$n[!is.na(df_summary$n)], probs = n_cap_pctile)]
  
  g_data_list <- build_g_data_scatter(df_graph = df_graph,
                                      alpha_val = 0.9,
                                      by_color = "source",
                                      color_label = "Data Source",
                                      color_scheme_scatter = color_scheme_scatter,
                                      size_lab = size_lab,
                                      title = title,
                                      reg_title = reg_title,
                                      region_name = region_name,
                                      table_data = table_data,
                                      base_font_size = base_font_size,
                                      year_min = year_min,
                                      year_max = year_max,
                                      stage3 = stage3,
                                      stage_3_gray = stage_3_gray)
  # Pull legends and then remove
  g_data_legend <- gLegend(g_data_list[[1]])
  g_data_list[[1]] <- g_data_list[[1]] + theme(legend.position="none")
  g_data <- fix_g_data_scatter(g_data_list, df_graph, stage3, stage_3_gray)
  message("  Successfully created scatter for main plot.")
  
  # Optionally make new data scatterplots
  if (new_data_plots){
    message("  You have chosen to make New Data Plots. Making New Data scatters...")
    g_data_new_list <- build_g_data_scatter(df_graph = df_graph,
                                            alpha_val = 1,
                                            by_color = "new_data_lab",
                                            color_label = "New Data",
                                            color_scheme_scatter = "binary",
                                            size_lab = size_lab,
                                            title = title,
                                            reg_title = reg_title,
                                            region_name = region_name,
                                            table_data = table_data,
                                            base_font_size = base_font_size,
                                            year_min = year_min,
                                            year_max = year_max,
                                            stage3 = stage3,
                                            stage_3_gray = stage_3_gray)
    g_data_new_legend <- gLegend(g_data_new_list[[1]])
    g_data_new_list[[1]] <- g_data_new_list[[1]] + theme(legend.position="none")
    g_data_new <- fix_g_data_scatter(g_data_new_list, df_graph, stage3, stage_3_gray)
  } else {
    message("  You have chosen not to make New Data Plots. Skipping new data scatters.")
    g_data_new_list <- NA
    g_data_new_legend <- NA
    g_data_new <- NA
  }
  
  
  message("Done making scatterplots.")
  return(list("g_data_legend" = g_data_legend,
              "g_data_new_legend" = g_data_new_legend,
              "g_data" = g_data,
              "g_data_new" = g_data_new))
  
}



## ###########################################################################
## SECTION 5 - GRAPHICS FUNCTIONS: DATA COVERAGE MAPS
## ###########################################################################


#' Build background map
#'
#' @description Builds a background map of a region with just the country borders. Grays out
#'   stage 3 countries if stage_3_gray is true
#'
#' @param region name of region to be mapped
#' @param endemic_gauls list of gaul codes. if non-null, the maps will be gray except
#'   for these countries
#' @param simplify_polys boolean. if true, simplifies shapefiles by reducing
#'   the number of vertices
#' @param tolerance numeric value that influences the degree of simplification of polygons
#'   (see Douglas-Peuker algorithm)
#' @param master_shape_all a combined SPDF with all polygons necessary to make maps
#' @param stage3 data.table of stage 3 countries to be removed
#' @param stage_3_gray boolean. If true, removes stage 3 countries from scatterplot
#'
#' @return returns a list of 4 objects -
#'   * 'background_outline_gg': gg geom_path, outline of region
#'   * 'background_map_gg': gg geom_poly, background map of region with country borders
#'   * 'background_map_gg_not_endemic': gg geom_poly, background map
#'        with countries marked non-endemic or stage 3 grayed out
#'   * 'extent(background_map)': extent of polygon defining region
#'
make_background_map <- function(region,
                                endemic_gauls,
                                simplify_polys,
                                tolerance,
                                fast_shapefiles,
                                master_shape_all,
                                stage3,
                                stage_3_gray) {
  
  # Create a background map for the plot (just the country polygons)
  if (region == "africa") {
    gaul_list <- get_gaul_codes('africa_dcp')
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == "africa_no_yem") {
    gaul_list <- get_gaul_codes('africa_dcp-yem')
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'south_asia') {
    gaul_list <- get_gaul_codes('south_asia_dcp')
    gaul_list <- unique(c(gaul_list, 154, 231)) # add Sri Lanka & Maldives
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'south_asia_ind_collaborators') {
    gaul_list <- get_gaul_codes('south_asia')
    gaul_list <- unique(c(gaul_list, 154, 231)) # add Sri Lanka & Maldives & PNG
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
    message("  This region requires some extra shapefile processing...")
    background_map[background_map@data$STATUS == "Sovereignty unsettled",
                   "GAUL_CODE"] <- 115
    # Keep the GAUL_CODEs for later
    bg_map_gaul_order <- as.data.frame(unique(background_map@data$GAUL_CODE))
    names(bg_map_gaul_order) <- 'GAUL_CODE'
    # Dissolve
    background_map <- gUnaryUnion(background_map, id = background_map@data$GAUL_CODE)
    # Merge the GAUL_CODEs back on and convert to spatial polygons data frame
    row.names(background_map) <- as.character(1:length(background_map))
    background_map <- SpatialPolygonsDataFrame(background_map, bg_map_gaul_order)
  }
  if (region == 'se_asia') {
    gaul_list <- get_gaul_codes('se_asia_dcp')
    gaul_list <- gaul_list[!(gaul_list %in% c(154, 231))] # remove Sri Lanka, Maldives
    add_countries <- c("PNG", "TWN")
    add_gauls <- get_gaul_codes(add_countries)
    gaul_list <- unique(c(gaul_list, add_gauls))
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'latin_america') {
    gaul_list <- get_gaul_codes('latin_america_dcp')
    gaul_list <- unique(c(gaul_list, 86, 63)) # Add French Guiana and Cuba
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'central_america') {
    gaul_list <- get_gaul_codes('central_america')
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'central_america_no_mex') {
    gaul_list <- get_gaul_codes('central_america')
    remove_countries <- c("MEX")
    remove_gauls <- get_gaul_codes(remove_countries)
    gaul_list <- gaul_list[!(gaul_list %in% remove_gauls)]
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'south_america') {
    gaul_list <- get_gaul_codes('south_america')
    gaul_list <- c(gaul_list, 86) # Add French Guiana
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'south_america_mex') {
    gaul_list <- get_gaul_codes('south_america')
    gaul_list <- c(gaul_list, 86) # Add French Guiana
    add_countries <- c("MEX")
    add_gauls <- get_gaul_codes(add_countries)
    gaul_list <- unique(c(gaul_list, add_gauls))
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'middle_east') {
    gaul_list <- get_gaul_codes('middle_east_dcp')
    add_countries <- c("UZB", "TKM", "TJK", "KGZ", "TUR", "ISR")
    add_gauls <- get_gaul_codes(add_countries)
    gaul_list <- unique(c(gaul_list, add_gauls))
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'eastern_europe') {
    gaul_list <- get_gaul_codes('eastern_europe')
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == "stage2") {
    background_map <- copy(master_shape_all)
  }
  
  if(simplify_polys == T) background_map <- simplify_spdf(background_map, tol = tolerance)
  
  background_map@data$id <- rownames(background_map@data)
  background_map_df <- suppressMessages(fortify(background_map))
  background_map_df <- merge(background_map_df, background_map@data[, c("id", "GAUL_CODE")], by = "id")
  background_map_df <- as.data.table(background_map_df)
  
  background_outline_gg <- suppressMessages(
    geom_path(data = background_map,
              aes(x = long,
                  y = lat,
                  group = group)))
  
  background_map_gg <- suppressMessages(
    geom_polygon(data = background_map,
                 aes(x = long,
                     y = lat,
                     group = group),
                 fill = "white"))
  
  if (!is.null(endemic_gauls)) {
    background_map_df_not_endemic <- background_map_df[!(GAUL_CODE %in% endemic_gauls), ]
    background_map_gg_not_endemic <- geom_polygon(data = background_map_df_not_endemic,
                                                  aes(x = long,
                                                      y = lat,
                                                      group = group),
                                                  fill = "lightgray")
  } else if (stage_3_gray) {
    background_map_df_not_endemic <- background_map_df[(GAUL_CODE %in% stage3$GAUL_CODE), ]
    background_map_gg_not_endemic <- geom_polygon(data = background_map_df_not_endemic,
                                                  aes(x = long,
                                                      y = lat,
                                                      group = group),
                                                  fill = "lightgray")
  } else {
    background_map_gg_not_endemic <- NULL
  }
  
  return(list(background_outline_gg, background_map_gg, background_map_gg_not_endemic, extent(background_map)))
}


#' get color list
#'
#' @description Gets hex values for colors to be used in maps based on `color_scheme`
#'
#' @param color_scheme character string of desired map coloring.
#'   options: `classic`, `darker_middle`, `red_blue`, `carto_red_blue`
#'
#' @return returns a list of color hex values
#'
get_color_list <- function(color_scheme) {
  # Set up the color scale
  if (color_scheme == "classic") {
    color_list <- c('#a50026', '#d73027', '#f46d43', '#fdae61',
                    '#fee090', '#ffffbf', '#e0f3f8', '#abd9e9',
                    '#74add1', '#4575b4', '#313695')
    
  } else if (color_scheme == "darker_middle") {
    color_list <- c("#A50026", "#960633", "#880D41", "#79144F",
                    "#6B1B5D", "#5C216B", "#4E2879", "#3F2F87",
                    "#313695")
    
  } else if (color_scheme == "red_blue") {
    color_list <- c("#A50026", "#B22E3C", "#C05D52", "#CD8B69",
                    "#DBBA7F", "#E9E996", "#C4C595", "#9FA195",
                    "#7A7D95", "#555995", "#313695")
    
  } else if (color_scheme == "carto_red_blue") {
    color_list <- c("#008080", "#70a494", "#b4c8a8", "#f6edbd",
                    "#edbb8a", "#de8a5a", "#ca562c")
    
  } else if (color_scheme %in% rownames(brewer.pal.info[brewer.pal.info$category == "seq",])) {
    color_list <- rev(brewer.pal(9, color_scheme)[-1])
    
  } else {
    stop(paste(color_scheme, "is not a recognized color scheme"))
  }
  
  return(color_list)
}


#' Make a data coverage map
#'
#' @description Takes a background map and adds polygon and point data to it
#'
#' @param period the year or set of years to make a map for
#' @param df data.table output from `dcp_merge_with_gbd_locations()`
#' @param region name of region to be mapped
#' @param poly_shapes_all a combined SPDF with all needed polygons,
#'   output by `pull_polys_in_parallel()`
#' @param background_map gg geom_path, outline of region. output of `make_background_map()`
#' @param background_outline gg geom_poly, background map of region with country borders.
#'   output of `make_background_map()`
#' @param background_map_not_epidemic gg geom_poly, background map
#'   with countries marked non-endemic or stage 3 grayed out. output of `make_background_map()`
#' @param background_extent extent of polygon defining region.
#'   output of `make_background_map()`
#' @param df_graph_poly data.table output of `dcp_make_df_graph()`, subsetted to countries
#'   in region
#' @param df_graph_point data.table output of `dcp_make_df_graph()`, subsetted to countries
#'   in region
#' @param not_real_shapefiles list of shapefiles not found in loading directory,
#'   output of `dcp_check_for_missing_shapefiles()`
#' @param color_list list of hex colors for map color scale, output of `get_color_list()`
#' @param legend_title Title of the legend (color scale) on the map
#' @param log_dir Path to a directory where a log file will
#'   optionally be saved. Useful mainly for debugging
#' @param base_font_size font size. All font sizes in the plots are scaled off this
#' @param map_point_size Size of points on map
#' @param cap numeric value, see `cap_type`
#' @param cap_type The type of numeric cap used to set
#'   the legend maximum. Must be one of `'percentile'`, `'absolute'`,
#'   or `'none'`.
#' @param legend_min Absolute lower bound for the map legend.
#'   Defaults to the lowest observation in the dataset.
#' @param legend_max Absolute upper bound for the map legend. If
#'   something other than `NA`, overrides the `cap` and `cap_type` arguments.
#' @param poly_line_width The default width of white lines
#'   surrounding polygons on the map.
#' @param annual_period_maps boolean. If true, map title is the year, else a year bin
#'
#' @return returns a complete data coverage map for one period
#'
make_a_period_map <- function(period,
                              df,
                              region,
                              poly_shapes_all,
                              background_map,
                              background_outline,
                              background_map_not_endemic,
                              background_extent,
                              df_graph_poly,
                              df_graph_point,
                              not_real_shapefiles,
                              color_list,
                              legend_title,
                              log_dir,
                              base_font_size,
                              map_point_size,
                              cap,
                              cap_type,
                              legend_min,
                              legend_max,
                              poly_line_width,
                              annual_period_maps) {
  
  # A blank theme
  theme_empty <- theme_classic(base_size = base_font_size) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  ## Subset data to this period
  df_period_polys <- df_graph_poly[plot_year == period, ]
  df_period_points <- df_graph_point[plot_year == period, ]
  
  ## Drop if shapefile not found
  df_period_polys <- df_period_polys[!(shapefile %in% not_real_shapefiles)]
  
  ## Initialize map
  g_map <- ggplot() + background_map + background_map_not_endemic
  
  # Merge polygons for this period to master spdf and add each one to map
  merge_poly_data <- function(this_survey, poly_shapes_all) {
    survey_subset <- df_graph_poly[nid == this_survey, ]
    survey_subset[,outcome_weighting := outcome * N]
    survey_subset <- survey_subset[,lapply(.SD,sum),
                                   .SDcols = c('N','outcome_weighting'),
                                   by=c('location_code','shapefile','nid','source')]
    survey_subset[, outcome := outcome_weighting / N]
    survey_subset[, outcome_weighting := NULL]
    survey_shapes <- try(sp::merge(poly_shapes_all, survey_subset,
                                   by=c('location_code','shapefile'),
                                   duplicateGeoms = TRUE))
    if (class(survey_shapes) == "try-error"){
      message(paste(this_survey, "is not merging properly. Please ensure unique location codes in geography codebook."))
      return(NULL)
    }
    survey_shapes <- survey_shapes[!is.na(survey_shapes@data$outcome),]
    return(survey_shapes)
  }
  
  # Function for printing informative warning messages about erroneous input data
  geog_warnings <- function(bad_data_df, geog_col, msg){
    setnames(bad_data_df, geog_col, 'geog_col')
    # List translating between geography field names and reporting names
    geog_report_names <- list(cluster_id = "Cluster IDs",
                              location_code = 'GAUL Codes')
    # Get all unique NIDs in the bad data
    for (bad_nid in unique(bad_data_df[,nid])){
      nid_df <- bad_data_df[nid==bad_nid,]
      source_name <- paste(unique(nid_df[,source]), collapse=', ')
      geogs <- paste(unique(nid_df[,geog_col]), collapse=', ')
      message(paste0("    ",msg,
                     "\n      Source: ", source_name,
                     "\n      NID:    ", bad_nid,
                     "\n      ",geog_report_names[[geog_col]],": ",geogs))
    }
  }
  
  # error checking - only do this if there is polygon data
  # otherwise, will leave the map as is - empty.
  
  # set up list of bad polys
  bad_polys <- c()
  poly_outside_list <- list()
  
  if (nrow(df_period_polys) > 0) {
    
    all_period_polys <- lapply(unique(df_period_polys[, nid]), function(this_survey) {
      merge_poly_data(this_survey, poly_shapes_all)
    })
    assign(paste0('poly_data_', period), all_period_polys)
    
    for(i in 1:length(all_period_polys)) {
      survey_spdf <- all_period_polys[[i]]
      
      if(length(survey_spdf) == 0) {
        bad_polys <- c(bad_polys, unique(df_period_polys[, survey])[i])
      }
      else {
        poly_dt <- suppressMessages(fortify(survey_spdf) %>% as.data.table)
        poly_dt <- try(merge(poly_dt, survey_spdf@data, by = "id"))
        if ("try-error" %in% class(poly_dt)){
          message("printing head of survey_spdf@data to see survey details")
          head(survey_spdf@data)
          next
        }
        
        # new table for polys that are outside the extent of the background
        poly_outside <-poly_dt[lat < background_extent@ymin | lat > background_extent@ymax |
                                 long < background_extent@xmin | long > background_extent@xmax,]
        
        poly_drop <-poly_dt[lat < background_extent@ymin - 0.5 |
                              lat > background_extent@ymax + 0.5 |
                              long < background_extent@xmin - 0.5 |
                              long > background_extent@xmax + 0.5, ]
        
        ## if a given polygon falls outside of the extend of the background region, throw a warning
        if(nrow(poly_outside)>0){
          geog_warnings(bad_data_df = poly_outside,
                        geog_col = 'location_code',
                        msg = "There are polygons outside of the background region range.")
          poly_outside_list[[length(poly_outside_list) + 1]] <- poly_outside
        }
        
        if(nrow(poly_drop)>0) {
          geog_warnings(bad_data_df = poly_drop,
                        geog_col = 'location_code',
                        msg = paste0("The following polygons are VERY far away",
                                     " from background polygon and will be dropped.\n",
                                     "    Check the location of these polygons!"))
          for (grp in unique(poly_drop$group)) {
            poly_dt <- subset(poly_dt, group != grp)
          }
        }
        
        g_map <- g_map +
          geom_polygon(data = poly_dt,
                       aes(x = long,
                           y = lat,
                           group = group,
                           fill = outcome),
                       color = "white",
                       size = poly_line_width)
        
      }
    }
    
    # Drop if more than 0.5 degrees away from extent of background map
    if (length(poly_outside_list) > 0) {
      poly_outside_list <- rbindlist(poly_outside_list)
    }
    
    
    
  } else {
    
    assign(paste0('poly_data_', period), NULL)
    
  }
  
  # Report out bad polys
  if (length(bad_polys) > 0) {
    warning("BAD SPECIFIC POLYGON MERGES - polygons have been dropped")
    if (is.null(log_dir)) {
      warning(paste0("List of affected surveys: ", paste(bad_polys, collapse = ", ")))
    } else if (!is.null(log_dir)) {
      as.data.table(bad_polys) %>%
        write.csv(., file = paste0(log_dir, "bad_polys.csv"))
    }
  }
  
  ## determine if there are points that fall outside the background extent
  points_outside <- df_period_points[latitude < background_extent@ymin | latitude > background_extent@ymax |
                                       longitude < background_extent@xmin | longitude > background_extent@xmax,]
  
  point_drop <- df_period_points[latitude < background_extent@ymin - 0.5 |
                                   latitude > background_extent@ymax + 0.5 |
                                   longitude < background_extent@xmin - 0.5 |
                                   longitude > background_extent@xmax + 0.5, ]
  
  ## if there are points outside of the extent, throw a warning and print the survey names + number of points outside
  if(nrow(points_outside)>0){
    geog_warnings(bad_data_df = points_outside,
                  geog_col = 'cluster_id',
                  msg = "There are points outside of the background region range.")
  }
  if(nrow(point_drop)>0) {
    geog_warnings(bad_data_df = point_drop,
                  geog_col = 'cluster_id',
                  msg = paste0("The following surveys have points that are VERY",
                               " far away from background polygon and will be dropped.\n",
                               "    Check the location of these points!")
    )
  }
  
  df_period_points <- subset(df_period_points,
                             !(latitude < background_extent@ymin - 0.5 | latitude > background_extent@ymax + 0.5 |
                                 longitude < background_extent@xmin - 0.5 | longitude > background_extent@xmax + 0.5))
  
  
  ## Set up color scale
  range <- range(df$outcome, na.rm=T)
  ## If specified, set manual legend range values
  range[2] <- ifelse(!(is.na(legend_max)), legend_max, range[2])
  ## Set legend max based on percentiles, if specified
  ## This will only take effect if the percentile cap is greater than the minimum
  cap_greater <- (quantile(df$outcome, cap/100, na.rm=T) > range[1])
  if (cap_type == "percentile" & is.na(legend_max)){
    if(cap_greater == TRUE){
      range[2] <- quantile(df$outcome, cap/100, na.rm=T)
    } else {
      # Print a warning message
      message(paste0("WARNING: The ",cap,"% upper legend bound specified is not greater ",
                     "than the legend minimum. The legend upper bound will remain ",
                     "at ",range[2],"."))
    }
  }
  ## Set legend based on absolute value, if specified
  ## This will only take effect if the absolute cap is greater than the minimum
  cap_greater <- (cap > range[1])
  if (cap_type == "absolute" & is.na(legend_max)){
    if(cap_greater){
      range[2] <- cap
    } else {
      # Print a warning message
      message(paste0("WARNING: The ",cap," absolute legend bound specified is ",
                     "not greater than the legend minimum. The legend upper ",
                     "bound will remain at ",range[2],"."))
    }
  }
  
  ## If the legend bounds are the same, add a very small number to the upper bound
  if(diff(range) <= 0){
    message(paste0("WARNING: The legend showed no variation in data. Increasing ",
                   "the upper bound slightly to avoid plotting erorrs."))
    range[2] <- range[1] + 1E-5
  }
  
  digits <- 10^(floor(log10(diff(range))) - 1)
  limits <- c(digits*floor(range[1]/digits), digits*ceiling(range[2]/digits))
  brks <- round(seq(limits[1], limits[2], length.out = 5), -1*log10(digits))
  brks[c(1, 5)] <- limits
  labs <- format(brks)
  if (cap_type != "none" | !(is.na(legend_max))) labs[5] <- paste0(labs[5])
  
  ## Finally, set up title text based on whether or not annual period maps are
  ##  being made
  if (annual_period_maps){
    title_text <- as.character(period)
  } else {
    if(period == 2000){
      title_text <- paste0(period, "-", period+2)
    }else{
      title_text <- paste0(period-2, "-", period+2)
    }
  }
  
  ## Add points and finishing touches
  g_map <- g_map +
    background_outline +
    geom_point(data = df_period_points,
               aes(x = longitude,
                   y = latitude,
                   fill = outcome),
               alpha = 1,
               size = map_point_size,
               color = "black",
               shape = 21,
               stroke = 0.1) +
    coord_equal() +
    theme_empty +
    scale_fill_gradientn(colours=color_list, limits=limits, breaks=brks, labels=labs,
                         na.value = color_list[length(color_list)], name = paste0(legend_title, "\n")) +
    guides(fill = guide_colorbar(barheight = 15, nbin=1000)) +
    labs(fill = "Outcome",
         title = title_text) +
    theme(plot.title = element_text(size = rel(1.5),
                                    face = "bold"))
  return(g_map)
}


#' Combine map objects
#'
#' @description Combines all elements into the final data coverage plot
#'
#' @param g_datamap data scatterplot, grob output of `make_data_scatterplots()`
#' @param g_data_legend data scatterplot legend, grob output of `make_data_scatterplots()`
#' @param map_list list of 4 period maps, outputted by `make_a_period_map()`
#' @param n_countries numeric, number of countries being mapped
#' @param reg_title title of the region being mapped
#' @param title title of data coverage plot
#' @param base_font_size font size. All font sizes in the plots are scaled off this
#' @param n_total numeric, total number of rows in dataset
#' @param polys_total numeric, total number of polygons
#' @param points_total numeric, total number of points
#'
#' @return returns a complete data coverage plot
#'
dcp_make_4up_map <- function(g_datamap,
                             g_data_legend,
                             map_list,
                             n_countries,
                             reg_title,
                             title,
                             base_font_size,
                             n_total,
                             polys_total,
                             points_total) {
  
  #g_datamap = object for left side of map (lets you add new data map
  
  # grab your legends using the predefined functions, then state their grid location
  p.legend <- gLegend(map_list[[1]])
  p.legend$vp <- grid::viewport(layout.pos.row = 1:12, layout.pos.col = 9)
  
  # Note g_data_legend grabbed above while processing g_data
  g_data_legend$vp <- grid::viewport(layout.pos.row = 2:11, layout.pos.col = 4)
  
  # Add a title
  title_grob <- textGrob(paste0(title, ":\n ", reg_title),
                         gp=gpar(fontsize=base_font_size*1.5))
  title_grob$vp <- grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3)

  # Add notes at bottom
  note_grob <- textGrob(paste0("N: ", formatC(n_total, format="d", big.mark=","), "\n",
                               "Points: ", formatC(points_total, format="d", big.mark=","), "\n",
                               "Polygons: ", formatC(polys_total, format="d", big.mark=",")),
                        gp = gpar(fontsize = base_font_size))
  note_grob$vp <- grid::viewport(layout.pos.row = 11:12, layout.pos.col = 9)
  
  # Set up based on number of countries
  if(n_countries > 24) g_datamap$vp <-  grid::viewport(layout.pos.row = 2:11, layout.pos.col = 1:3)
  if(n_countries > 12 & n_countries <= 24) g_datamap$vp <-  grid::viewport(layout.pos.row = 3:10, layout.pos.col = 1:3)
  if(n_countries <= 12) g_datamap$vp <-  grid::viewport(layout.pos.row = 4:9, layout.pos.col = 1:3)
  
  # Initialize plot with master title
  grid.newpage()
  pushViewport(grid::viewport(layout = grid.layout(12, 9)))
  vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y, clip = "off")
  # Plot all data coverage maps
  #print(tbl, as.table=TRUE)
  grid.draw(g_datamap) # Note now a grob, so grid.draw
  print(map_list[[1]] + theme(legend.position="none"), vp = vplayout(1:6, 5:6))
  print(map_list[[2]] + theme(legend.position="none"), vp = vplayout(1:6, 7:8))
  print(map_list[[3]] + theme(legend.position="none"), vp = vplayout(7:12, 5:6))
  print(map_list[[4]] + theme(legend.position="none"), vp = vplayout(7:12, 7:8))
  grid.draw(p.legend)
  grid.draw(g_data_legend)
  grid.draw(title_grob)
  grid.draw(note_grob)
  
}


## ###########################################################################
## SECTION 2 - SHINY PREP FUNCTIONS
## ###########################################################################


#' Save data for shiny
#'
#' @description Save the requisite objects for the data coverage shiny
#'
#' @param df data.table output from `dcp_merge_with_gbd_locations()`
#' @param df_graph_poly data.table output of `dcp_make_df_graph()`, subsetted to countries
#'   in region
#' @param poly_shapes_all a combined SPDF with all needed polygons,
#'   output by `pull_polys_in_parallel()`
#' @param var Name of the field in `df` that will be plotted on the maps.
#' @param indicator The indicator being estimated for this data coverage plot.
#'   This argument controls where the data coverage plots will be saved
#'
prep_data_coverage_shiny <- function(df,
                                     df_graph_poly,
                                     poly_shapes_all,
                                     var,
                                     indicator) {
  
  ## Save polygons with data for Data Coverage Shiny
  makeID_save <- function(shape) {
    shapefile <- all_polygon_data[[shape]]
    if(length(shapefile)!=0) {
      shapefile <- spChFIDs(shapefile, paste0(shape, "_", row.names(shapefile@data)))
      return(shapefile)
    }
    else {
      message(paste0("Not saving for geo_data Shiny: ", unique(df_graph_poly[, survey])[shape]))
    }
  }
  
  # Check if any polygons
  if (nrow(df_graph_poly) > 0) {
    
    ## Merge polygons for this period to master spdf and add each one to map
    merge_poly_data <- function(this_survey) {
      survey_subset <- df_graph_poly[survey == this_survey, ]
      survey_shapes <- try(merge(poly_shapes_all,
                                 survey_subset,
                                 by=c('location_code','shapefile'),
                                 duplicateGeoms = TRUE))
      if (class(survey_shapes) == "try-error"){
        message(paste(this_survey, "is not merging properly. Please ensure unique location codes in geography codebook."))
        return(NULL)
      }
      survey_shapes <- survey_shapes[!is.na(survey_shapes@data$outcome),]
      return(survey_shapes)
    }
    
    all_polygon_data <- lapply(unique(df_graph_poly[, survey]), merge_poly_data)
    all_polygon_data <- lapply(1:length(all_polygon_data), makeID_save)
    all_polygon_data <- all_polygon_data[!sapply(all_polygon_data, is.null)]
    all_polygon_data <- do.call(rbind, all_polygon_data)
    all_polygon_data$polygon_survey <- paste0(all_polygon_data$source, '_', all_polygon_data$shapefile, '_', all_polygon_data$year)
    just_polygons <- all_polygon_data[c('GAUL_CODE','shapefile')]
    just_polygons$gaul_shape <- paste0(just_polygons$GAUL_CODE, just_polygons$shapefile)
    just_polygons <- just_polygons[which(!duplicated(just_polygons$gaul_shape)), ]
  } else {
    just_polygons <- list(NULL)
  }
  
  save(list = 'just_polygons', file = '<<<< FILEPATH REDACTED >>>>')
  setnames(df, var, indicator)
  write.csv(df, '<<<< FILEPATH REDACTED >>>>')
  
}
