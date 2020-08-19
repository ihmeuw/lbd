input_aggregate_admin <- function(indicator,
                                  indicator_group,
                                  regions = c("cssa", "wssa", "essa", "sssa", "name"),
                                  run_date = NULL,
                                  input_data = NULL,
                                  indicator_family = "binomial",
                                  svy_id = "nid",
                                  sample_column = "sum_of_sample_weights",
                                  subnational_nids = NULL,
                                  shapefile_version = NULL) {

  # Make sure this is run on singularity container in Rstudio or lbd singularity to use sf package
  if (!is_lbd_singularity() & !is_rstudio(check_singularity = TRUE)) {
    stop("must be run in lbd singularity or Rstudio singularity container")
  }

  require(dplyr)
  require(data.table)
  require(sf)
  require(ggplot2)
  require(stringr)

  configurations <- fread(<<<< FILEPATH REDACTED >>>>)
  config_sh_version <- configurations[V1 == "modeling_shapefile_version", V2]

  if (is.null(shapefile_version) == T) {
    shapefile_version <- config_sh_version
    message(paste0("Using modeling_shapefile_version specified in config: ", shapefile_version))
  }

  # If input_data not passed to the function, load in input data from model, or from mbg if the appropriate column name is not specified.
  if (is.null(input_data)) {
    mod_dir <- sprintf(<<<< FILEPATH REDACTED >>>>)
    message(paste0(
      "Input data was not provided, so loading from model directory: ", mod_dir,
      ". \n If this does not exist, consider passing input data as an argument to this function"
    ))
    input_data <- fread(<<<< FILEPATH REDACTED >>>>)

    if (!sample_column %in% names(input_data)) {
      message("Sample weights column not found in model directory, pulling from input data on J drive")
      input_data <- fread(<<<< FILEPATH REDACTED >>>>)
      if (!sample_column %in% names(input_data)) stop("Sample weights column not found on J drive, make sure to add a sum of sample weights column)                                         in collapse code and specify columns name as sample_column argument in function")
    }
  } else {
    if (!sample_column %in% names(input_data)) {
      stop("Sample weights column not found in the provided input data. Make sure that this column is included and specified by the sample_column
            argument in the function call.")
    }
  }

  message(paste0("input_aggregate_admin input_data 1: ", nrow(input_data), " rows"))

  if (!"source" %in% names(input_data)) stop("Need to specify the source of the data under a source column in the input data")
  if (!"point" %in% names(input_data)) stop("You need a column in your input data called 'points' that classifies point data as 1 and polygon as 0")

  # Load in shapefile
  admin_shp <- st_read(get_admin_shapefile(admin_level = 2, version = shapefile_version), quiet = T)

  # Add country name (not just 3 letter abbreviation) to input data
  # Subset input data to given region defined by the regions argument
  gaul_codes <- get_adm0_codes(regions, shapefile_version = shapefile_version)

  message(paste0("gaul_codes: ", gaul_codes))

  gaul_to_loc_id <-
    get_location_code_mapping(shapefile_version = shapefile_version) %>%
    dplyr::select(GAUL_CODE, loc_name, ihme_lc_id) %>%
    dplyr::rename(location_name = loc_name) %>%
    filter(GAUL_CODE %in% gaul_codes)

  if (colnames(input_data)[1] == "V1" & colnames(input_data)[2] == "V1") {
    input_data <- input_data[, 2:ncol(input_data)]
  }
  
  # Join input data to location names to match shapefile. Some of those countries need to be manually changed to fit shapefile names
  input_data <-
    input_data %>%
    left_join(gaul_to_loc_id, by = c("country" = "ihme_lc_id")) %>%
    rowwise() %>%
    ungroup() %>%
    data.table() %>%
    setnames(c(svy_id, sample_column), c("svy_id", "sample_column"))

  message(paste0("input_aggregate_admin input_data 2: ", nrow(input_data), " rows"))

  missing <-
    input_data %>%
    filter(is.na(location_name)) %>%
    nrow()

  message(paste0(round(100 * missing / nrow(input_data), 2), " % of data is from outside of specified regions: ", paste(regions, collapse = " ")))

  input_data <-
    input_data %>% filter(!is.na(location_name))

  message(paste0("input_aggregate_admin input_data 3: ", nrow(input_data), " rows")) ## casch 08/11/2018: inserted for debugging

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

  if ("ADM0_CODE.x" %in% colnames(input_admin)) {
    input_admin <- as.data.table(input_admin)
    input_admin[, ADM0_CODE.x := NULL]
    setnames(input_admin, "ADM0_CODE.y", "ADM0_CODE")
  }

  missing <-
    input_admin %>%
    filter(is.na(ADM0_CODE), is.na(ADM1_CODE), is.na(ADM2_CODE)) %>%
    nrow()

  message(paste0(round(100 * missing / nrow(input_admin), 3), " % of data could not be matched to an admin level"))

  input_admin <-
    input_admin %>%
    filter(!is.na(ADM0_CODE), !is.na(ADM1_CODE), !is.na(ADM2_CODE)) %>%
    data.table() # remove those that are not assigned to an admin0/admin2

  # Make sure they are assigned to the correct country, some are right over the border and so are pushed into a different area.
  # These are excluded (usually not substantial) becuase they then do not nest well with admin1 and admin2 estimates
  wrong_admin0 <- nrow(input_admin[location_name != ADM0_NAME])
  message(paste0(round(wrong_admin0 / nrow(input_admin), 3), " % of input data is matched to a different country.\nThese are usually located on the border and will be dropped for the visualization."))

  input_admin <- input_admin[location_name == ADM0_NAME]

  # If binomial make sure it is prevalence space
  if (indicator_family == "binomial") input_admin[, prev := prev / N]

  # Collapse to admin 0 level
  message("collapsing to admin 0 level")
  input_admin0 <-
    input_admin %>%
    group_by(svy_id, source, point, ADM0_NAME, ADM0_CODE, data_collect_method, diagnostic) %>%
    dplyr::summarise(
      year = floor(median(year, na.rm = T)),
      outcome = weighted.mean(prev, sample_column),
      N = sum(N * weight)
    ) %>%
    ungroup() %>%
    data.table()

  # Collapse to admin 1 level
  message("collapsing to admin 1 level")
  input_admin1 <-
    input_admin %>%
    group_by(svy_id, source, point, ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, data_collect_method, diagnostic) %>%
    dplyr::summarise(
      year = floor(median(year, na.rm = T)),
      outcome = weighted.mean(prev, sample_column),
      N = sum(N * weight)
    ) %>%
    ungroup() %>%
    data.table()

  # Collapse to admin 2 level
  message("collapsing to admin 2 level")
  input_admin2 <-
    input_admin %>%
    group_by(svy_id, source, point, ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE, data_collect_method, diagnostic) %>%
    dplyr::summarise(
      year = floor(median(year, na.rm = T)),
      outcome = weighted.mean(prev, sample_column),
      N = sum(N * weight)
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
      mutate(source = ifelse(nchar(source) > 6, str_trunc(source, 6, ellipsis = ""), source)) %>% # truncate source if it is too long and not specified above
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

  if (!is.null(subnational_nids)) {
    input_admin0 <- subnational_nid_subset(input_admin0)
    input_admin1 <- subnational_nid_subset(input_admin1)
    input_admin2 <- subnational_nid_subset(input_admin2)
  }

  # Final list contains each aggregated admin
  input_admins <- list(ad0 = input_admin0, ad1 = input_admin1, ad2 = input_admin2)
  return(input_admins)
}
