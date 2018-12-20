# Functions to graph data coverage

########################################################################################
##### MAIN FUNCTION ####################################################################
########################################################################################

graph_data_coverage_values <- function(df,
                                       var,
                                       title,
                                       year_min,
                                       year_max,
                                       year_var,
                                       region,
                                       cores,
                                       indicator,
                                       high_is_bad = T,
                                       return_maps = T,
                                       legend_title,
                                       since_date = NULL,
                                       sum_by = NULL,
                                       out_dir = NULL,
                                       extra_file_tag = '',
                                       save_on_share = FALSE,
                                       endemic_gauls = NULL,
                                       log_dir = NULL,
                                       fast_shapefiles = TRUE,
                                       simplify_polys = TRUE,
                                       tolerance = 0.03,
                                       remove_rank = T,
                                       base_font_size = 18,
                                       prep_shiny = T,
                                       color_scheme = "classic",
                                       color_scheme_scatter = "brewer",
                                       cap = 90,
                                       cap_type = "percentile") {

  ## I. Set up & ensure correct packages ##############################################
  message("\n###################################################################")
  message(paste0("\nGenerating ", indicator, " graphs for ", region, "..."))

  ## Sink all random graphical output to NULL; will specify png where graphs actually desired
  pdf(NULL)

  ## Manually keep order_var space_rank for now
  order_var <- 'space_rank'

  ## Define root
  assign("<<<< FILEPATH REDACTED >>>>>", "<<<< FILEPATH REDACTED >>>>>", envir = .GlobalEnv)

  if (Sys.info()["sysname"] != "Linux") stop("This function should be run on the cluster only.")

  # Load packages
  message("\nLoading packages...")
  package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                          paste0(j_root,'<<<< FILEPATH REDACTED >>>>>'),
                          paste0(j_root,'<<<< FILEPATH REDACTED >>>>>'))
  package_list <- c('data.table', 'ggplot2', 'parallel', 'doParallel', 'gridExtra',
                    'stringr', 'RColorBrewer', 'rgdal', 'sp', 'raster', 'magrittr',
                    'dplyr', 'RMySQL', 'rgeos', 'tidyr')

  for(package in package_list) {
    library(package, lib.loc = package_lib, character.only=TRUE)
  }

  # Define repo
  repo <- paste0('<<<< FILEPATH REDACTED >>>>>', Sys.info()["user"],'/mbg/')
  if (dir.exists(repo) == F) repo <- "<<<< FILEPATH REDACTED >>>>>"

  # Set up for GBD functions & gaul_convert
  source(paste0(repo, 'mbg_central/gbd_functions.R'))
  source(paste0(repo, 'mbg_central/prep_functions.R'))

  # Set up ggplot base font size
  theme_set(theme_minimal(base_size = base_font_size))

  ## II. Prep data #####################################################################

  # 1. Pull in country table -----------------------------------------------------------
  message("Pulling GBD region list...")
  # Pull in list of regions from GBD
  country_table <-suppressMessages(suppressWarnings(
                  data.table(get_location_hierarchy(41))[, .(ihme_loc_id,
                                                              level,
                                                              location_name,
                                                              location_name_short,
                                                              region_name)]
                  ))
  # For now, subset to countries only
  #   Note: could change this for subnationals if desired!
  #   May be helpful if looking at a country for which we have tons of data (e.g. India)

  country_table <- country_table[level == 3]
  setnames(country_table, "ihme_loc_id", "country")

  # 2. Load data & split off relevant data sets ---------------------------------------
  message("Preparing data sets...")
  # Create df_poly, df_point, and df_summary
  prepped <- clean_and_prep_df(df, country_table, year_var, var, region)
    df         <- prepped[["df"]]                # Full data set
    df_poly    <- prepped[["df_poly"]]           # Polygon data only
    df_point   <- prepped[["df_point"]]          # Point data only
    df_summary <- prepped[["df_summary"]]        # One-row-per-NID summary
    rm(prepped)

  # Identify new data & write out results to standard .csv file
  message("Checking to see which data is new...")
  df_summary <- find_new_data(df_summary, indicator, since_date)

  # Grab a country list
  list_output <- get_country_list(df, region)
    reg_title    <- list_output[["reg_title"]]    # Formatted region title
    region_list  <- list_output[["region_list"]]  # Subheadings for regions
    country_list <- list_output[["country_list"]] # ISO3 codes for countries
    rm(list_output)

  # Make a summary data set to use for graphing, subset to region
  df_graph <- make_df_graph(df_summary, country_list)

  # Similarly subset df_point and df_poly to just the region of interest
  df_graph_point <- df_point[country %in% country_list]
  df_graph_poly <- df_poly[country %in% country_list]

  # III. Prepare graphs ###############################################################

  # 1. Setup --------------------------------------------------------------------------

  # This takes a while, so let's do it once only
  if (!("master_shape_all" %in% ls())) {

    if (fast_shapefiles == T) {
      message("\nFast-loading master shapefile... ")
      assign("master_shape_all",
             readRDS('<<<< FILEPATH REDACTED >>>>>master_shape_all.rds'),
             envir = globalenv())
    } else {
    message("\nOpening master shapefile... (good time to go make a cup of coffee)")
    assign("master_shape_all",
            readOGR(dsn = paste0(j_root,"<<<< FILEPATH REDACTED >>>>>/template_shp"),
                    layer = "data_coverage_gaul_template"),
            envir = globalenv())
    }
  }

  # Rename ADM0_CODE field
  names(master_shape_all)[names(master_shape_all) == "ADM0_CODE"] <- "GAUL_CODE"

  # 2X. Make table for the left side of graph-------------------------------------------

  # Deprecated - old make_table
  #summary_tbl_list <- make_table(df, gaul_list = get_gaul_codes(region), order_var = order_var, country_list = country_list, remove_rank = remove_rank)
  #summary_tbl <- summary_tbl_list[[1]]
  #table_data <- summary_tbl_list[[2]]

  table_data <- make_table_new(df_summary, country_list, year_min, year_max)
  polys_total <- sum(table_data$Polygons)
  points_total <- sum(table_data$Points)
  n_total <- sum(table_data$N)

  # 2a. Make plot for left side of graph ------------------------------------------------
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
                                              color_scheme_scatter = color_scheme_scatter)
    g_data     <- data_scatter_list[["g_data"]]                    # scatter for all data
    g_data_new <- data_scatter_list[["g_data_new"]]                # scatter for just new data
    g_data_legend  <- data_scatter_list[["g_data_legend"]]         # legend for g_data
    g_data_new_legend  <- data_scatter_list[["g_data_new_legend"]] # legend for g_data_new
    rm(data_scatter_list)

  # 3. Make plot for right side (map) of figure ------------------------------------------

  # The polygon map works by pulling shapefiles for each polygon in the
  #   data set, and then coloring those by the lastest year. Polygons are
  #   overlaid chronologically, so the most recent of two overlapping
  #   polygons will be visible.

  # Set up point / polygon graphing data sets
  #   will create plots for polygons (latest year of data = color)
  #   and points (color = year of data) separately

  message("Preparing maps...")

  # Bin years to the four standard age bins
  df_graph_poly  <- bin_years(df_graph_poly)
  df_graph_point <- bin_years(df_graph_point)

  # Only do this bit if there are polygons!
  if(nrow(df_graph_poly) > 0) {

    # Check for missing shapefiles & report if present
    shapefile_list <- check_for_missing_shapefiles(df_graph_poly)
    shapefiles <- shapefile_list[["shapefiles"]]
    not_real_shapefiles <- shapefile_list[["not_real_shapefiles"]]

    if (length(not_real_shapefiles) > 0) {
      warning(paste0("DROPPING ", length(not_real_shapefiles), " MISSING SHAPES (in shapefile column but not shapefile dir)!"))
      if (!is.null(log_dir)) {
        warning(paste0("Writing list of missing shapefiles to ", log_dir, "missing_shapes.csv"))
        not_real_shapefiles %>%
          as.data.table %>%
          setnames(., ".", "Missing shapefiles") %>%
          write.csv(., file = paste0(log_dir, "missing_shapes.csv"))
      } else {
        warning(paste0("Missing shapefiles: ", paste(not_real_shapefiles, collapse = ", ")))
      }
    }

    rm(shapefile_list) # clean up

    ## Pull shapes in parallel ----------------------------------------------------------

    # Make a table of shapefiles & associated location codes
    df_shape_loc <- unique(df_graph_poly[, c("shapefile", "location_code")])

    # Pull all polygons in parallel
    poly_list <- pull_polys_in_parallel(shape_loc_list = df_shape_loc,
                                        shapefile_col = "shapefile",
                                        location_code_col = "location_code",
                                        cores = cores,
                                        fast_shapefiles = fast_shapefiles)
    message("Done pulling polys")

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

  } else {

    # If no polygons, set products of the above to NULL
      poly_shapes_all <- NULL
      not_real_shapefiles <- NULL
      broken_shapes <- NULL

  }

    # 3a: Define map basics -----------------------------------------------------------

    # Create a background map for the plot (just the country polygons)
    background_map_list <- make_background_map(region, endemic_gauls, simplify_polys, tolerance, fast_shapefiles, master_shape_all)
    background_outline <- background_map_list[[1]]
    background_map <- background_map_list[[2]]
    background_map_not_endemic <- background_map_list[[3]]
    background_extent <- background_map_list[[4]]
    rm(background_map_list)

    color_list <- get_color_list(color_scheme)
    if (high_is_bad==TRUE) color_list <- rev(color_list)

    # 3b: Make some graphs -------------------------------------------------------------

    for(period in c(2000,2005,2010,2015)) {
     g_map <- suppressMessages(make_a_period_map(period, df, region, poly_shapes_all,
                                                 background_map, background_outline, background_map_not_endemic, background_extent,
                                                 df_graph_poly, df_graph_point, not_real_shapefiles,
                                                 color_list, legend_title, log_dir, base_font_size, cap, cap_type))
     assign(paste0('map_',period), g_map)
    }

  # 4. Set up for save -----------------------------------------------------------------

  if (is.null(out_dir)) out_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator, extra_file_tag, '/')
  if(save_on_share==TRUE) out_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator, extra_file_tag, '/')
  dir.create(out_dir, showWarnings = FALSE, recursive = T)

  # Make the main plots - 4-up

    out_file <- paste0(out_dir, 'data_coverage_', region, '.png')
    message(paste0('Saving ', out_file))
    unlink(out_file) # Delete any old versions

    png(filename=out_file,
        units = "in",
        width = 24.33,
        height = 12,
        pointsize = base_font_size,
        res = 450)

    make_4up_map(g_datamap = g_data,
                 g_data_legend = g_data_legend,
                 summary_tbl = summary_tbl,
                 map_list = list(map_2000, map_2005, map_2010, map_2015),
                 n_countries = length(unique(df_graph$country)),
                 reg_title = reg_title,
                 title = title,
                 base_font_size = base_font_size,
                 n_total = n_total,
                 polys_total = polys_total,
                 points_total = points_total)

    dev.off()

  # Repeat the main plot, but for new data

    out_file <- paste0(out_dir, 'data_coverage_', region, '_new.png')
    message(paste0('Saving ', out_file))
    unlink(out_file) # Delete any old versions

    png(filename=out_file,
        units = "in",
        width = 24.33,
        height = 12,
        pointsize = base_font_size,
        res = 450)

    make_4up_map(g_datamap = g_data_new,
                 g_data_legend = g_data_new_legend,
                 summary_tbl = summary_tbl,
                 map_list = list(map_2000, map_2005, map_2010, map_2015),
                 n_countries = length(unique(df_graph$country)),
                 reg_title = reg_title,
                 title = title,
                 base_font_size = base_font_size,
                 n_total = n_total,
                 polys_total = polys_total,
                 points_total = points_total)

    dev.off()

  # Make the year bin plots - one for each of four years

    for (period in c(2000,2005,2010,2015)) {

      out_file <- paste0(out_dir, 'map_', period, '_', region, '.png')
      message(paste0('Saving ', out_file))
      unlink(out_file) # Delete any old versions

      png(filename=out_file,
          units = "in",
          width = 17,
          height = 10,
          pointsize = 24,
          res = 450)

      print(get(paste0("map_", period)))

      dev.off()

    }

  # 5. Finish up ----------------------------------------------------------------------

  if (prep_shiny == T) {
    # Save some shiny outputs
    prep_data_coverage_shiny(df, df_graph_poly, poly_shapes_all, var, indicator)
  }

  ## Return individual maps if requested
  if(return_maps==TRUE) return(list(map_2000, map_2005, map_2010, map_2015))

}


########################################################################################
##### DATA PREP FUNCTIONS ##############################################################
########################################################################################

# --------------------------------------------------------------------------------------
# Clean and prep your data; merge on country identifiers; get to standard format

clean_and_prep_df <- function(df, country_table, year_var, var, region) {
  # Take df, ensure correct format, and output the relevant bits

  # 1. Preliminary cleaning -----------------------------------------------------------

  df <- as.data.table(df)  # in case passed a data frame

  # Check if latitude and longitude are numeric
  if(!is.numeric(df$latitude) | !is.numeric(df$longitude)) {
    stop("Latitude and longitude not both numeric, please fix.")
  }

  # Truncate survey series names
  df$survey_name <- substr(df$survey_name, 1, 15)

  # Rename the year variable to something standard
  setnames(df, year_var, "year_var")

  # Recode NA attribution for shapefile
  df[shapefile == "", shapefile := NA]

  # Generate pointpoly variable
  df[!is.na(latitude), pointpoly := "Point"]
  df[!is.na(shapefile) & is.na(pointpoly), pointpoly := "Polygon"]
  df <- df[!is.na(pointpoly)]

  # 2. Merge on country identifiers ---------------------------------------------------

  df <- merge(df, country_table, by = "country", all = T)

  if (region == "africa") {
    # Rename "North Africa and Middle East" to just "North Africa"
    df <- df[region_name == "North Africa and Middle East", region_name := "North Africa"]
  } else if (region == "middle_east") {
    # Rename "North Africa and Middle East" to just "Middle East"
    df <- df[region_name == "North Africa and Middle East", region_name := "Middle East"]
  }

  # Drop non-matched countries & notify user
  message(paste0("\nDropping ", nrow(df[is.na(country)]),
                 " rows without matches in GBD country table."))
  message(paste0("Countries affected: ", unique(df[is.na(country)]$location_name)))
  df <- df[!is.na(country)]

  # Tuncate long country names
  df <- df[location_name == "Democratic Republic of the Congo", location_name := "DRC"]
  df <- df[location_name == "Central African Republic", location_name := "CAR"]
  df <- df[location_name == "Sao Tome and Principe", location_name := "STP"]
  df <- df[location_name == "United Arab Emirates", location_name := "UAE"]
  df <- df[location_name == "Equatorial Guinea", location_name := "Eq. Guinea"]
  df <- df[location_name == "Saint Vincent and the Grenadines", location_name := "St. Vin. & Grenadines"]


  # 3. Generate subsets of data for further analysis ----------------------------------

  # Split off a point & polygon data set for use later
  df[, outcome := get(var)]
  df_point <- df[pointpoly == "Point"]
  df_poly  <- df[pointpoly == "Polygon"]

  # Sum over the N of the group - includes rows for countries with no data
  df_summary <- df[, .(n = sum(N), count = .N), by = .(source, country, year_var, pointpoly,
                                           location_name, region_name, svy_id)]

  return(list("df" = df, "df_point" = df_point, "df_poly" = df_poly, "df_summary" = df_summary))

}

# --------------------------------------------------------------------------------------
# Take a df_summary object and check to see if there's new data; mark if so

find_new_data <- function(df_summary, indicator, since_date) {

  # Look for an existing data summary table
  summary_dir <- paste0(j_root, "<<<< FILEPATH REDACTED >>>>>")
  summary_file <- paste0(summary_dir, indicator, "_data_summary.csv")

  df_summary$date <- as.character(Sys.Date())

  if (file.exists(summary_file)) {
    # read in old file
    df_summary_old <- read.csv(summary_file, stringsAsFactors = F) %>% as.data.table

    # grab svy ids & dates from old file
    svy_id_dates_old <- unique(df_summary_old[, c("svy_id", "date")])
    svy_id_dates_old <- svy_id_dates_old[!is.na(svy_id)]
    names(svy_id_dates_old) <- c("svy_id", "date_old")

    # merge
    df_summary$svy_id <- as.numeric(df_summary$svy_id)
    df_summary <- merge(df_summary, svy_id_dates_old, all.x = T)

    # replace if an older date exists
    df_summary[date != date_old, date := date_old]
    df_summary[, date_old := NULL]

  }

  #replace dates for country-rows with no data with "na"
  df_summary[is.na(svy_id), date := NA]
  df_summary <- df_summary[order(location_name)]

  unlink(summary_file)
  write.csv(df_summary, file = summary_file)

  # Mark which data is new with a variable "new_data" (for graphing)
  if (is.null(since_date)) {
    df_summary[date == max(as.Date(df_summary$date), na.rm = T), new_data := 1]
    df_summary[!is.na(date) & is.na(new_data), new_data := 0]
  } else {
    df_summary[as.Date(date) > as.Date(since_date) & !is.na(date), new_data := 1]
    df_summary[!is.na(date) & is.na(new_data), new_data := 0]
  }

  return(df_summary)

}

# --------------------------------------------------------------------------------------
# Pull a country list based on your data & region of choice
# Outputs title for the region, sub-region titles, and vector of countries

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

    remove_countries <- c("AFG", "ARE", "IRN", "IRQ", "JOR", "OMN", "PSE", "SAU", "SYR", "TUR", "YEM","KWT","LBN","QAT","BHR","CPV")

    country_list <- country_list[!(country_list %in% remove_countries)]
  }

  if (region == 'south_asia') {
    reg_title <- "South Asia"
    region_list <- c("South Asia")
    country_list <- unique(df[region_name %in% region_list]$country)

    # add Sri Lanka & the Maldives & Papau New Guinea
    country_list <- unique(c(country_list, "LKA", "MDV"))
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
    country_list <- unique(c(country_list, "PNG"))
  }

  if (region == 'latin_america') {
    reg_title <- "Latin America and Caribbean"
    region_list <- c("Andean Latin America",
                     "Caribbean",
                     "Central Latin America",
                     "Tropical Latin America")
    country_list <- unique(df[region_name %in% region_list]$country)
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
    # Uzbekistan, Turkmenistan, Tajikistan, Kyrgyzstan
    country_list <- unique(c(country_list, "UZB", "TKM", "TJK", "KGZ"))

  }

  return(list("reg_title" = reg_title,
              "region_list" = region_list,
              "country_list" = country_list))
}

# --------------------------------------------------------------------------------------
# Make the df_graph object - subsetted to countries in country_list & formatted

make_df_graph <- function(df_summary, country_list) {
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

}

# --------------------------------------------------------------------------------------
# Take a data frame and bin year_var

bin_years <- function(df_to_bin) {
  df_to_bin <- df_to_bin[, survey := paste0(source,'_',country,'_',year_var)]
  df_to_bin <- subset(df_to_bin, year_var >= 1998)
  df_to_bin <- df_to_bin[year_var > 1997 & year_var <= 2002, plot_year := 2000]
  df_to_bin <- df_to_bin[year_var > 2002 & year_var <= 2007, plot_year := 2005]
  df_to_bin <- df_to_bin[year_var > 2007 & year_var <= 2012, plot_year := 2010]
  df_to_bin <- df_to_bin[year_var > 2012 & year_var <= 2017, plot_year := 2015]

 return(df_to_bin)
}

########################################################################################
##### MISC FUNCTIONS FOR GRAPHING ######################################################
########################################################################################

# --------------------------------------------------------------------------------------
# Extract legend out of a graph to use for final plot
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

# --------------------------------------------------------------------------------------
# Build the data scatterplot (left side of graph)
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
                                 year_max) {

  # Set up table data
  td <- copy(table_data)
  td[, N := NULL] # Don't display N column
  setnames(td, "Country", "country")
  td <- gather(td, type, value, -country) %>% as.data.table

  # Add on regions
  reg_table <- subset(df_graph, select = c("country", "location_name", "region_name")) %>%
                unique
  td <- merge(td, reg_table, by="country")

  base_font_size <- round(base_font_size * 0.8)

  panel_spacing <- 2
  if (region_name == "africa") panel_spacing <- 1


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
                   axis.text.y = element_blank(),,
                   axis.text.x = element_text(size = rel(1)),
                   panel.spacing.y = unit(panel_spacing, "lines")) +
             facet_wrap(~region_name, scales="free_y", ncol = 1)

  # Builds a data scatterplot
  g_datamap <- ggplot(data = df_graph,
                      aes(x = year_var,
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

# --------------------------------------------------------------------------------------
# Adjust the data scatterplot to ensure even distribution of rows
fix_g_data_scatter <- function(g_data_scatter_list, df_graph) {

    # This function takes the data map objects and makes it
    # so that all of the heights are evenly distributed
    # returns a grob object instead of a ggplot object

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

# --------------------------------------------------------------------------------------
# Basically a wrapper for build_g_data_scatter to build both old & new scatterplots
make_data_scatterplots <- function(df_graph, df_summary, title, reg_title, year_min, year_max, table_data, base_font_size, region_name, color_scheme_scatter) {

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
                                        year_max = year_max)

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
                                            year_max = year_max)

    # Pull legends and then remove
    g_data_legend <- gLegend(g_data_list[[1]])
    g_data_list[[1]] <- g_data_list[[1]] + theme(legend.position="none")

    g_data_new_legend <- gLegend(g_data_new_list[[1]])
    g_data_new_list[[1]] <- g_data_new_list[[1]] + theme(legend.position="none")

    g_data <- fix_g_data_scatter(g_data_list, df_graph)
    g_data_new <- fix_g_data_scatter(g_data_new_list, df_graph)

    return(list("g_data_legend" = g_data_legend,
                "g_data_new_legend" = g_data_new_legend,
                "g_data" = g_data,
                "g_data_new" = g_data_new))

}

# --------------------------------------------------------------------------------------
# Make a background map for a given region

make_background_map <- function(region, endemic_gauls, simplify_polys, tolerance, fast_shapefiles, master_shape_all) {

  # Create a background map for the plot (just the country polygons)
  if (region == "africa") {
    gaul_list <- get_gaul_codes('africa')
    gaul_list <- gaul_list[!gaul_list == 269] #remove yemen
    ## Add new background map from Lucas
    if (fast_shapefiles == T) {
      background_map <- readRDS('<<<< FILEPATH REDACTED >>>>>/background_map_africa.rds')
    } else {
      background_map <- readOGR(dsn="<<<< FILEPATH REDACTED >>>>>", layer="africa_ad0")
    }

    names(background_map)[names(background_map) == "ADM0_CODE"] <- "GAUL_CODE"
  }
  if (region == 'south_asia') {
    gaul_list <- get_gaul_codes('south_asia')
    gaul_list <- unique(c(gaul_list, 154, 231)) # add Sri Lanka & Maldives & PNG
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'se_asia') {
    gaul_list <- get_gaul_codes('se_asia')
    gaul_list <- gaul_list[!(gaul_list %in% c(154, 231))] # remove Sri Lanka, Maldives
    add_countries <- c("PNG")
    add_gauls <- get_gaul_codes(add_countries)
    gaul_list <- unique(c(gaul_list, add_gauls))
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'latin_america') {
    gaul_list <- get_gaul_codes('latin_america')
    gaul_list <- c(gaul_list, 86) # Add French Guiana
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
    gaul_list <- get_gaul_codes('middle_east')
    add_countries <- c("UZB", "TKM", "TJK", "KGZ", "YEM")
    add_gauls <- get_gaul_codes(add_countries)
    gaul_list <- unique(c(gaul_list, add_gauls))
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
  }
  if (region == 'eastern_europe') {
    gaul_list <- get_gaul_codes('eastern_europe')
    background_map <- master_shape_all[master_shape_all@data$GAUL_CODE %in% gaul_list, ]
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
  } else {
    background_map_gg_not_endemic <- NULL
  }

  return(list(background_outline_gg, background_map_gg, background_map_gg_not_endemic, extent(background_map)))
}

# --------------------------------------------------------------------------------------
# Make a map of data coverage for a given period

make_a_period_map <- function(period, df, region,
                              poly_shapes_all, background_map,
                              background_outline,
                              background_map_not_endemic,
                              background_extent,
                              df_graph_poly, df_graph_point,
                              not_real_shapefiles,
                              color_list, legend_title,
                              log_dir, base_font_size, cap, cap_type) {

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
    survey_subset <- df_graph_poly[survey == this_survey, ]
    survey_shapes <- try(merge(poly_shapes_all, survey_subset, by=c('location_code','shapefile')))
    if (class(survey_shapes) == "try-error"){
      message(paste(this_survey, "is not merging properly. Please ensure unique location codes in geography codebook."))
      return(NULL)
    }
    survey_shapes <- survey_shapes[!is.na(survey_shapes@data$outcome),]
    return(survey_shapes)
  }

  # error checking - only do this if there is polygon data
  # otherwise, will leave the map as is - empty.

  # set up list of bad polys
  bad_polys <- c()
  poly_outside_list <- list()

  if (nrow(df_period_polys) > 0) {

    all_period_polys <- lapply(unique(df_period_polys[, survey]), function(this_survey) {
      merge_poly_data(this_survey, poly_shapes_all)
    })
    assign(paste0('poly_data_', period), all_period_polys)

    for(i in 1:length(all_period_polys)) {
      survey_spdf <- all_period_polys[[i]]

      if(length(survey_spdf) == 0) {
        bad_polys <- c(bad_polys, unique(df_period_polys[, survey])[i])
      }
      else {
        poly_dt <- fortify(survey_spdf) %>% as.data.table
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
          for(surv in unique(poly_outside$survey)){
            ## Display the unique GAUL_CODEs and surveys to the screen
            warning(paste0("There are polygons outside of the background region range.\n",
                           "  Codebook: ", poly_outside[survey == surv, unique(source)], "\n",
                           "  Survey: ", poly_outside[survey == surv, unique(survey_name)], "\n",
                           "  NID: ", poly_outside[survey == surv, unique(svy_id)], "\n",
                           "  Incorrect GAUL_CODE: ", poly_outside[survey == surv, unique(GAUL_CODE)], "\n"))
          }
          poly_outside_list[[length(poly_outside_list) + 1]] <- poly_outside
        }

        if(nrow(poly_drop)>0) {
          for (grp in unique(poly_drop$group)) {
            warning(paste0("The following polygons are VERY far away from background polygon and will be dropped. \n",
                           "Check the location of these polygons!\n",
                           "  Codebook: ", poly_drop[group == grp, unique(source)], "\n",
                           "  Survey: ", poly_drop[group == grp, unique(survey_name)], "\n",
                           "  NID: ", poly_drop[group == grp, unique(svy_id)], "\n",
                           "  Incorrect GAUL_CODE: ", poly_drop[group == grp, unique(GAUL_CODE)], "\n"))
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
                     size = 0.2)

      }
    }

    # Drop if more than 0.5 degrees away from extent of background map
    if (length(poly_outside_list) > 0) {
      poly_outside_list <- rbindlist(poly_outside_list)
    }



  } else if (nrow(df_period_polys) == 0) {

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
    for(surv in unique(points_outside$survey)){
      ## Display the unique surveys and number of points to the screen
      warning(paste0("There are points outside of the background region range.\n",
                     "  Codebook: ", points_outside[survey == surv, unique(source)], "\n",
                     "  Survey: ", points_outside[survey == surv, unique(survey_name)], "\n",
                     "  NID: ", points_outside[survey == surv, unique(svy_id)], "\n",
                     "  Incorrect Cluster ID: ", points_outside[survey == surv, unique(cluster_id)], "\n"))         }
  }

  if(nrow(point_drop)>0) {
    for (surv in unique(point_drop$surv)) {
      warning(paste0("The following surveys have points that are VERY far away from background polygon and will be dropped. \n",
                     "Check the location of these polygons!\n",
                     "  Codebook: ", points_drop[survey == surv, unique(source)], "\n",
                     "  Survey: ", points_drop[survey == surv, unique(survey_name)], "\n",
                     "  NID: ", points_drop[survey == surv, unique(svy_id)], "\n",
                     "  Incorrect GAUL_CODE: ", points_drop[survey == surv, unique(GAUL_CODE)], "\n"))
    }
  }

  df_period_points <- subset(df_period_points,
                             !(latitude < background_extent@ymin - 0.5 | latitude > background_extent@ymax + 0.5 |
                               longitude < background_extent@xmin - 0.5 | longitude > background_extent@xmax + 0.5))


  ## Set up color scale
  range <- range(df$outcome, na.rm=T)
  if (cap_type == "percentile") range[2] <- quantile(df$outcome, cap/100, na.rm=T)
  if (cap_type == "absolute") range[2] <- cap

  digits <- 10^(floor(log10(diff(range))) - 1)
  limits <- c(digits*floor(range[1]/digits), digits*ceiling(range[2]/digits))
  brks <- round(seq(limits[1], limits[2], length.out = 5), -1*log10(digits))
  brks[c(1, 5)] <- limits
  labs <- format(brks)
  if (cap_type != "none") labs[5] <- paste0(labs[5], "+")

  ## Add points and finishing touches
  g_map <- g_map +
    background_outline +
    geom_point(data = df_period_points,
               aes(x = longitude,
                   y = latitude,
                   fill = outcome),
               alpha = 1,
               size = 0.8,
               color = "black",
               shape = 21,
               stroke = 0.1) +
    coord_equal() +
    theme_empty +
    scale_fill_gradientn(colours=color_list, limits=limits, breaks=brks, labels=labs,
                         na.value = color_list[length(color_list)], name = paste0(legend_title, "\n")) +
    guides(fill = guide_colorbar(barheight = 15, nbin=1000)) +
    labs(fill = "Outcome",
         title = paste0(period-2, "-", period+2)) +
    theme(plot.title = element_text(size = rel(1.5),
                                    face = "bold"))

}


# --------------------------------------------------------------------------------------
# make the table that sits on the left-hand side of the final plot

make_table <- function(df, gaul_list, order_var, country_list, remove_rank) {

  ## Create coverage ranks for table (geo_rank = points/polys per capita, n_rank = sample / population, space_rank = 1 point for points in a year, 0.25 points for polygons)
  source("<<<< FILEPATH REDACTED >>>>>/get_population.R")
  source("<<<< FILEPATH REDACTED >>>>>/get_location_metadata.R")

  loc_meta <- suppressMessages(suppressWarnings(get_location_metadata(location_set_id=2, gbd_round_id=4)))
  loc_meta <- loc_meta[, c('location_id','ihme_loc_id'), with=F]
  pops <- suppressMessages(suppressWarnings(get_population(location_id="-1", sex_id="3", age_group_id="22", year_id="-1", gbd_round_id=4)))
  pops <- merge(pops, loc_meta, by='location_id')
  pops <- pops[year_id %in% c(2000,2005,2010,2015), ]
  setnames(pops, "ihme_loc_id", "country")
  pops <- pops[, list(population=sum(population)), by="country"]

  table_df <- df[, c('country','year_var','pointpoly','N'), with=F]
  table_df <- table_df[country %in% country_list]
  table_df <- table_df[, count := 1]
  table_df <- table_df[year_var >= 1998 & year_var <= 2002, year_id := 2000]
  table_df <- table_df[year_var >= 2003 & year_var <= 2007, year_id := 2005]
  table_df <- table_df[year_var >= 2008 & year_var <= 2012, year_id := 2010]
  table_df <- table_df[year_var >= 2013 & year_var <= 2017, year_id := 2015]
  table_df <- table_df[!is.na(year_id),]
  table_space_rank <- table_df
  table_df <- table_df[, list(total_geo=sum(count), total_n=sum(N)), by = c('country','pointpoly')]
  table_df <- table_df[!is.na(pointpoly),]
  table_df <- merge(pops, table_df, by=c('country'), all.x=TRUE)
  table_df <- table_df[!is.na(pointpoly),]
  table_df <- table_df[, geo_per_capita := total_geo / population]
  table_df <- table_df[, n_per_capita := total_n / population]
  table_df <- table_df[, list(geo_per_capita=sum(geo_per_capita), n_per_capita=sum(n_per_capita)), by='country']

  ## Subset to region
  table_df <- table_df[, gaul := gaul_convert(table_df[, country])]
  table_df <- table_df[gaul == 40764, gaul := 6]
  table_df <- table_df[gaul %in% gaul_list, ]
  table_df <- table_df[, gaul := NULL]

  table_df <- table_df[order(-geo_per_capita)]
  table_df <- table_df[, geo_rank := seq(1:length(table_df[, geo_per_capita]))]
  table_df <- table_df[order(-n_per_capita)]
  table_df <- table_df[, n_rank := seq(1:length(table_df[, n_per_capita]))]
  table_space_rank <- unique(table_space_rank[, c('country','year_id','pointpoly')])
  table_space_rank <- table_space_rank[!is.na(pointpoly),]
  table_space_rank <- table_space_rank[pointpoly=="Polygon", space_rank := 0.01]
  table_space_rank <- table_space_rank[pointpoly=="Point", space_rank := 1]
  table_space_rank <- table_space_rank[, list(space_rank=sum(space_rank)), by='country']
  table_df <- merge(table_df, table_space_rank, by='country')
  table_df <- table_df[, geo_per_capita := NULL]
  table_df <- table_df[, n_per_capita := NULL]
  setnames(table_df, 'country', 'location_name')
  table_df <- table_df[order(space_rank)]

  ### Make table of point/polygon counts
  table_data <- copy(df)
  table_data <- table_data[country %in% country_list]
  table_data <- table_data[!is.na(pointpoly), ]
  table_data[, units := 1]
  table_data <- table_data[, list(units=sum(units)), by=c('country', 'pointpoly')]
  table_data <- data.table::dcast(table_data, country ~ pointpoly, value.var = "units")

  # catch if no point data
  if ("Point" %in% names(table_data)) {
    setnames(table_data, 'Point', 'points')
    } else {
      table_data$points <- rep(0, nrow(table_data))
    }

  # catch if no polygon data
  if ("Polygon" %in% names(table_data)) {
    setnames(table_data, 'Polygon', 'polygons')
    } else {
      table_data$polygons <- rep(0, nrow(table_data))
    }

  # Replace NA with 0
  table_data[is.na(polygons), polygons := 0]
  table_data[is.na(points), points := 0]

  table_data$location_name <- factor(table_data$country,
                                     levels = rev(sort(unique(table_data$country))))
  table_data$country <- NULL

  table_total <- table_data[, list(polygons=sum(polygons, na.rm=T), points=sum(points, na.rm=T)),]
  table_total <- table_total[, Country := "TOTAL:"]

  ## Merge on ranks of data coverage
  table_data <- merge(table_data, table_df, by='location_name')
  table_data <- table_data[order(get(order_var))]

  #table_data <- rbind(table_data, table_total, fill=T)
  setcolorder(table_data, c('location_name','points','polygons','geo_rank','n_rank','space_rank'))
  setnames(table_data, 'location_name', 'Country')
  table_data <- table_data[, c('Country','points','polygons','space_rank')]

  ## Add any totally missing countries
  all_countries <- data.table('Country'=country_list)
  table_data <- merge(all_countries, table_data, by='Country', all.x=TRUE)
  table_data <- table_data[is.na(polygons), polygons := 0]
  table_data <- table_data[is.na(points), points := 0]
  table_data <- table_data[is.na(space_rank), space_rank := 0]
  ## Add back TOTAL
  table_data <- rbind(table_data, table_total, fill=T)

  ## Clean up some names, and make the final TOTAL: row bold by splitting out the table and apprending the grobs.
  table_data <- table_data[, space_rank := as.character(space_rank)]
  table_data <- table_data[Country=="TOTAL:", space_rank := "-"]
  setnames(table_data, "polygons", "Polygons")
  setnames(table_data, "points", "Points")
  setnames(table_data, "space_rank", "Rank")
  if(remove_rank == T) table_data[, Rank := NULL]
  tg1 <- tableGrob(table_data[1:nrow(table_data), 1:ncol(table_data)], rows=NULL,
                  theme = ttheme_minimal(base_size = 7))
  tg2 <- tableGrob(table_data[1, 1:ncol(table_data)], rows=NULL, # can't have empty content
                  cols=as.character(table_data[nrow(table_data), 1:ncol(table_data)]),
                  theme = ttheme_minimal(base_size = 7)) # use 4th row as header
  summary_tbl <- rbind(tg1[-nrow(tg1), ], tg2[1,])

  return(list(summary_tbl, table_data))

}

# --------------------------------------------------------------------------------------
# Make a table to append to the scatter plot

make_table_new <- function(df_summary, country_list, year_min, year_max) {
  td <- subset(df_summary, year_var >= year_min & year_var <= year_max)
  td <- subset(td, select = c("country", "pointpoly", "n", "count"))
  td <- subset(td, country %in% country_list)
  td <- td %>%
    group_by(country, pointpoly) %>%
    summarize(Count = sum(count), N = sum(n)) %>%
    as.data.table

  td[is.na(pointpoly), N := 0]
  td[is.na(pointpoly), Count := 0]

  td_n <- subset(td, select=c("country", "N")) %>%
          group_by(country) %>%
          summarize(N = sum(N)) %>%
          as.data.table

  td_count <- td %>%
                .[, N := NULL] %>%
                spread(pointpoly, Count, fill = 0) %>%
                as.data.table %>%
                merge(., td_n, by = "country") %>%
                setnames("country", "Country")

  # Catch if either no point or no poly data
  if (!("Point" %in% names(td_count))) {
    td_count[, Point := 0]
  } else if (!("Polygon" %in% names(td_count))) {
    td_count[, Polygon := 0]
  }

  setnames(td_count, c("Point", "Polygon"), c("Points", "Polygons"))

  # Add in countries with no data
  no_data_countries <- country_list[!(country_list %in% unique(td_count$Country))]
  no_data_td <- data.table(Country = no_data_countries,
                           Points = 0,
                           Polygons = 0,
                           N = 0)

  td_count <- rbind(td_count, no_data_td)

  if ("<NA>" %in% names(td)) td[, "<NA>" := NULL]

  setcolorder(td_count, c("Country", "Points", "Polygons", "N"))

  return(td_count)
}

# --------------------------------------------------------------------------------------
# Combine all of the map objects into a 4-up map for the right side of the final figure

make_4up_map <- function(g_datamap, g_data_legend, summary_tbl, map_list, n_countries, reg_title, title, base_font_size,
                         n_total, polys_total, points_total) {

  #g_datamap = object for left side of map (lets you add new data map

  # grab your legends using the predefined functions, then state their grid location
  p.legend <- gLegend(map_list[[1]])
  p.legend$vp <- viewport(layout.pos.row = 1:12, layout.pos.col = 9)

  # Note g_data_legend grabbed above while processing g_data
  g_data_legend$vp <- viewport(layout.pos.row = 2:11, layout.pos.col = 4)

  #Legacy form of summary table
  #summary_tbl$vp <- viewport(layout.pos.row = 1:6, layout.pos.col = 1:3)

  # Add a title
  title_grob <- textGrob(paste0(title, ": ", reg_title),
                    gp=gpar(fontsize=base_font_size*1.5))
  title_grob$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1:3)

  # Add notes at bottom
  note_grob <- textGrob(paste0("N: ", formatC(n_total, format="d", big.mark=","), "\n",
                               "Points: ", formatC(points_total, format="d", big.mark=","), "\n",
                               "Polygons: ", formatC(polys_total, format="d", big.mark=",")),
                        gp = gpar(fontsize = base_font_size))
  note_grob$vp <- viewport(layout.pos.row = 11:12, layout.pos.col = 9)

  # Set up based on number of countries
  if(n_countries > 24) g_datamap$vp <-  viewport(layout.pos.row = 2:11, layout.pos.col = 1:3)
  if(n_countries > 12 & n_countries <= 24) g_datamap$vp <-  viewport(layout.pos.row = 3:10, layout.pos.col = 1:3)
  if(n_countries <= 12) g_datamap$vp <-  viewport(layout.pos.row = 4:9, layout.pos.col = 1:3)

  # Initialize plot with master title
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(12, 9)))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y, clip = "off")
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
  #grid.draw(summary_tbl)

}

########################################################################################
##### POLYGON FUNCTIONS ################################################################
########################################################################################

# --------------------------------------------------------------------------------------
# Pull polygons for a single shapefile
# Input `shape_loc` should be a chunk of rows from a larger data.table,
#   with columns `shapefile` (the same for all rows) and
#   a column `location_code` that has individual location codes of interest

pull_polys <- function(shape_loc, fast_shapefiles) {
  # Function to pull in a shapefile and assign each polygon the latest year
  # that any data was collected within it (determines fill color)

  # This is computationally intensive if many shapefiles

  shape <- unique(shape_loc$shapefile)
  loc_codes <- unique(shape_loc$location_code)

  message(paste0("Working on ", shape))

  # Read in the master shapefile
  if (fast_shapefiles == T) {
    master_shape <- fast_load_shapefile(shape)
  } else  {
    master_shape <- readOGR(dsn = paste0(j_root, "<<<< FILEPATH REDACTED >>>>>/Shapefile directory"),
                            layer = shape)
  }

  names(master_shape)[names(master_shape)=="GAUL_Code"] <- "GAUL_CODE"

  # Subsetting will break if NAs in row index (GAUL_CODE)
  master_shape <- master_shape[!is.na(master_shape$GAUL_CODE),]

  # Custom fixes for broken shapefiles
  if (shape == "AZE_DHS_2006" & "GAUL_CODE" %in% names(master_shape) == F) {
    master_shape$GAUL_CODE <- master_shape$REGCODE
  }

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

# --------------------------------------------------------------------------------------
# Pull polygons in parallel. See inside function for description of use

pull_polys_in_parallel <- function(shape_loc_list,
                                    shapefile_col = "shapefile",
                                    location_code_col = "location_code",
                                    cores,
                                    fast_shapefiles = fast_shapefiles) {

  # shape_loc_list:    data table with two columns,
  #                      (one with shapefiles & other with location codes)
  # shapefile_col:     name of column with shapefiles
  # location_code_col: name of column with location codes
  # cores:             number of cores to use in parallelization

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
      .libPaths("<<<< FILEPATH REDACTED >>>>>")
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

# --------------------------------------------------------------------------------------
# Check for missing shapefiles in the df_graph_poly object
# Return a list of those that exist & a list of those that don't

check_for_missing_shapefiles <- function(df_graph_poly) {
  # Pull a list of shapefiles for the polygon data for the region in question
  shapefiles <- unique(df_graph_poly$shapefile) %>% as.character

  ## Check that all shapefile entries are real shapefiles and report bad entries
  real_shapefiles <- gsub('.shp',
                          '',
                          list.files(paste0(j_root, "<<<< FILEPATH REDACTED >>>>>"),
                                     pattern = '.shp'))
  not_real_shapefiles <- shapefiles[!(shapefiles %in% real_shapefiles)]

  shapefiles <- shapefiles[(shapefiles %in% real_shapefiles)]

  return(list("shapefiles" = shapefiles,
              "not_real_shapefiles" = not_real_shapefiles))
}

########################################################################################
##### SHINY FUNCTIONS ##################################################################
########################################################################################


# --------------------------------------------------------------------------------------
# Save the requisite objects for the data coverage shiny

prep_data_coverage_shiny <- function(df, df_graph_poly, poly_shapes_all, var, indicator) {

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
      survey_shapes <- try(merge(poly_shapes_all, survey_subset, by=c('location_code','shapefile')))
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
    all_polygon_data$polygon_survey <- paste0(all_polygon_data$source, '_', all_polygon_data$shapefile, '_', all_polygon_data$year_var)
    just_polygons <- all_polygon_data[c('GAUL_CODE','shapefile')]
    just_polygons$gaul_shape <- paste0(just_polygons$GAUL_CODE, just_polygons$shapefile)
    just_polygons <- just_polygons[which(!duplicated(just_polygons$gaul_shape)), ]
  } else {
    just_polygons <- list(NULL)
  }

  save(list = 'just_polygons', file = paste0('<<<< FILEPATH REDACTED >>>>>', indicator, '_polygons.RData'))
  setnames(df, var, indicator)
  write.csv(df, paste0('<<<< FILEPATH REDACTED >>>>>', indicator, '_points.csv'))

}

simplify_spdf <- function(spdf, tol = tolerance) {
  df_spdf <- data.frame(spdf)
  spdf <- gSimplify(spdf, tol = tol, topologyPreserve = T)
  spdf <- SpatialPolygonsDataFrame(spdf, df_spdf)
  return(spdf)
}

get_color_list <- function(color_scheme) {
  # Set up the color scale
  if (color_scheme == "classic") {
    color_list <- c('#a50026','#d73027','#f46d43', '#fdae61',
                    '#fee090','#ffffbf','#e0f3f8','#abd9e9',
                    '#74add1','#4575b4','#313695')
  }

  if (color_scheme == "darker_middle") {
    color_list <- c("#A50026","#B22E3C","#C05D52","#CD8B69",
                    "#DBBA7F","#E9E996","#C4C595","#9FA195",
                    "#7A7D95", "#555995","#313695")
  }

  if (color_scheme == "red_blue") {
    color_list <- c("#A50026", "#960633", "#880D41", "#79144F",
                    "#6B1B5D", "#5C216B", "#4E2879", "#3F2F87",
                    "#313695")
  }

  if (color_scheme == "carto_red_blue") {
    color_list <- c("#008080","#70a494","#b4c8a8", "#f6edbd",
                    "#edbb8a","#de8a5a", "#ca562c")
  }

  return(color_list)
}

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
