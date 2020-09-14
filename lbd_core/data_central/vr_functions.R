################################# VR Functions ################################
#
# Functions associated with pulling prepped VR data and other related data
#
# Date: 5/2/2019
#
#

#######################################################################################
# Pulling Functions

#' @title vr_get_base_path
#'
#' @description Get base path of the VR folder - defined in one place in case of changes later on
#'
#' @return A file path to the first level of the vr folder with prepped vr data
#'
vr_get_base_path <- function() {
  return("<<<< FILEPATH REDACTED >>>>")
}


#' @title vr_list_sources
#'
#' @description Get a vector of prepped VR sources (typically in the form iso3_admlevel)
#'
#' @return A character vector of available VR sources
#'
vr_list_sources <- function(){
  base_path <- vr_get_base_path()
  sources <- list.files(path = paste0(base_path, "extracts/"))
  sources <- sources[sources != "archive"]
  return(sources)
}


#' @title vr_list_source_files
#'
#' @description Get a character vector of files within a source folder
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param folder_name Boolean, include folder name in returned list? For personal use to help locate files
#'
#' @return A character vector of available VR sources
#'
vr_list_source_files <- function(source, folder_name = F) {
  base_path <- vr_get_base_path()
  source_files <- list.files(path = paste0(base_path, "extracts/", source), recursive = T)
  if(!folder_name){
    source_files <- unlist(lapply(source_files, basename))
  }
  return(source_files)
}


#' @title vr_check_file_exists
#'
#' @description check if a file exists given a source and file. Stops if the source or file does not exist
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param file string, a file name within the `source` directory
#'
#' @return NULL
#'
vr_check_file_exists <- function(source, file) {
  viable_sources <- vr_list_sources()
  if(!(tolower(source) %in% viable_sources)){
    message("When looking for file - ", file, " - from source - ", source)
    stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
  }
  source_files <- vr_list_source_files(source)

  if(!(file %in% source_files)){
    message("When looking for file - ", file, " - from source - ", source)
    stop("file is not valid - use vr_list_source_files() to find available files", call. = F)
  }
}


#' @title vr_pull_loc_meta
#'
#' @description pull a versioned location metadata (table 5) for a source.
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param version string, default `current`, can also take a date in format "YYYYMMDD"
#'
#' @return data.table of location metadata
#'
vr_pull_loc_meta <- function(source, version = "current"){
  base_path <- vr_get_base_path()
  file <- paste0("adm_meta_full_", version, ".csv")
  vr_check_file_exists(source, file)

  filepath <- paste0(base_path, "extracts/", source, "/location_data/", file)

  dt <- fread(filepath)
  return(dt)
}


#' @title vr_pull_shp
#'
#' @description pull a shapefile from the vr_prep directory
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param type string, either `stable` or `annual`
#' @param year string, with different options based on `type`. If `type` is `annual`, `year` can be "first", "last" or a year character, i.e. "2000". first and last will return the first and last available years respectively. If `type` is `stable`, `year` can be "full", which returns the stable shapefile over the longest time period, or a vector of 2 year characters. When years are passed in, shapefiles with these years must be present in the shapefile_single_year (annual) or shapefile_stable (stable) folders.
#' @param return_object string, "sf", "sp", or "filepath". Determines what object to return.
#' @param link boolean, return link objects as well?
#'
#' @return list object with 3 objects: `shapefile`: shapefile object specified by `return_object`, `link_table`: associated link table if `link` = T else NULL, and `id_ras`: associated id raster if `link` = T else NULL
#'
vr_pull_shp <- function(source, type, year, return_object = "filepath", link = T) {
  #checks for correct input
  if(!(source %in% vr_list_sources())){
    stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
  }

  if(!(type %in% c("stable", "annual"))) {
    stop("type must be either stable or annual in vr_pull_shp", call. = F)
  }

  if(type == "annual"){
    if(!(year %in% c("first", "last", 1950:2020))) {
      stop("type specified as annual, year must be 'first', 'last' or a year character, i.e. '2000'", call. = F)
    }
  }
  if(type == "stable"){
    if((!any(year %in% c("full"))) & ((class(year) != "character") | (length(year) != 2) | (!(all(year %in% as.character(c(1950:2020))))))) {
      stop("type specified as stable, year must be 'full', or a vector of 2 year characters i.e. c('2000', '2010')", call. = F)
    }
  }

  if(!(return_object %in% c("sp", "sf", "filepath"))) {
    stop("return object must be 'sp', 'sf', or 'filepath'")
  }

  #finding the right file
  file_list <- vr_list_source_files(source)
  if(type == "annual"){
    shp_names <- grep("shapefile_single_year.*shp", file_list, value = T)
    shp_years <- as.numeric(gsub("[^\\d]+", "", shp_names, perl=TRUE))

    if(year == "first") {
      file_year <- as.character(min(shp_years, na.rm = T))
    } else if(year == "last") {
      file_year <- as.character(max(shp_years, na.rm = T))
    } else {
      file_year <- year
    }

    file_path <- paste0("shapefile_single_year_", file_year, ".shp")
  }

  if(type == "stable") {
    shp_names <- grep("shapefile_stable.*shp", file_list, value = T)
    shp_years <- gsub("[^\\d]+", "", shp_names, perl=TRUE)
    years_list <- lapply((strsplit(shp_years, "(?<=.{4})", perl = TRUE)), as.numeric)

    if(any(year == 'full')){
      years_diff <- unlist(lapply(years_list, function(x){x[2] - x[1]}))
      file_year <- paste(unlist(years_list[which.max(years_diff)]), collapse = "_")
    } else {
      file_year <- paste(year, collapse = "_")
    }

    file_path <- paste0("shapefile_stable_", file_year, ".shp")
  }

  #check the file exists
  vr_check_file_exists(source, file_path)

  #load/convert as specified by return_object param
  base_path <- paste0(vr_get_base_path(), "extracts/", source, if(type == "annual") "/shapefile_single_year/" else "/shapefile_stable/")
  if(return_object == "sp"){
    shp <- readOGR(paste0(base_path, file_path))
  } else if (return_object == "sf"){
    shp <- st_read(paste0(base_path, file_path))
  } else {
    shp <- paste0(base_path, file_path)
  }

  #Read in link tables if specified, otherwise return NULL
  if(link){
    link_filepaths <- vr_convert_link_filepath(paste0(base_path, file_path))
    if (file.exists(link_filepaths$link_path)) {
      link <- readRDS(link_filepaths$link_path)
    } else {
      message("There was no associated link table found for this shapefile, use vr_gen_link_tables() to make one.")
      link <- NULL
    }

    if (file.exists(link_filepaths$id_path)) {
      id <- readRDS(link_filepaths$id_path)
    } else {
      id <- NULL
    }

  } else {
    link <- NULL
    id <- NULL
  }

  return(list("shapefile" = shp, "link_table" = link, "id_raster" = id))
}


#' @title vr_pull_cod
#'
#' @description pull cause of death data for a source
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param stage string, either `formatted` or `redistributed`. Stages of cod prep pipeline.
#' @param cause NOT IMPLEMENTED, causes to subset data by
#'
#' @return dt of cause of death VR data
#'
##TODO update with formal naming schema once data is available, implement cause param
vr_pull_cod <- function(source, stage, cause = "all"){
  if(!(source %in% vr_list_sources())){
    stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
  }
  if(!(stage %in% c("corrections", "disaggregation", "formatted", "misdiagnosiscorrection", "redistribution"))){
    stop('stage must be "corrections", "disaggregation", "formatted", "misdiagnosiscorrection", or "redistribution"', call. = F)
  }
  file_name <- paste0("cod_", stage, ".csv")
  base_path <- paste0(vr_get_base_path(), "extracts/", source, "/cod_output/")

  vr_check_file_exists(source, file_name)

  dt <- fread(paste0(base_path, file_name))

  return(dt)
}


#' @title vr_pull_u5m
#'
#' @description pull prepped under-5 mortality data
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param version string, default `current`, can also take a date in format "YYYYMMDD"
#'
#' @return dt of u5m prepped VR data
#'
vr_pull_u5m <- function(source, version = "current"){
  if(!(source %in% vr_list_sources())){
    stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
  }

  file_name <- if(stage == "current") paste0("u5m_", source, "_current.csv") else paste0("u5m_", source, "_", version, ".csv")
  base_path <- paste0(vr_get_base_path(), "extracts/", source, "/u5m/")

  vr_check_file_exists(source, file_name)

  dt <- fread(paste0(base_path, file_name))

  return(dt)
}


#' @title vr_pull_births
#'
#' @description pull prepped VR birth data
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param version string, default `current`, can also take a date in format "YYYYMMDD"
#'
#' @return dt of prepped VR births data
#'
vr_pull_births <- function(source, version = "current"){
  if(!(source %in% vr_list_sources())){
    stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
  }

  file_name <- if(stage == "current") paste0("births_", source, "_current.csv") else paste0("births_", source, "_", version, ".csv")
  base_path <- paste0(vr_get_base_path(), "extracts/", source, "/births/")

  vr_check_file_exists(source, file_name)

  dt <- fread(paste0(base_path, file_name))

  return(dt)
}


#' @title vr_pull_link_table
#'
#' @description pull a link table from the vr_prep directory. Requires either source/type/year ala vr_pull_shp, or shapefile_path. If shapefile_path is not NULL, it will be prioritized.
#'
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param type string, either `stable` or `annual`
#' @param year string, with different options based on `type`. If `type` is `annual`, `year` can be "first", "last" or a year character, i.e. "2000". first and last will return the first and last available years respectively. If `type` is `stable`, `year` can be "full", which returns the stable shapefile over the longest time period, or a vector of 2 year characters. When years are passed in, shapefiles with these years must be present in the shapefile_single_year (annual) or shapefile_stable (stable) folders.
#' @param shapefile_path default `NULL`, can be used in place of the above 3 options to pull the link table. Takes a string filepath to the .shp file of the shapefile that you want the link table and id raster for.
#' @return list object with 2 objects: `link_table`: associated link table and `id_raster`: associated id raster
#'
#' @example
#' using vr_pull_shp parameters:
#' link <- vr_pull_link_table(source = "chl_adm3",
#'                            type = "stable",
#'                            year = "full")
#'
#'
vr_pull_link_table <- function(source, type, year, shapefile_path = NULL) {
  if(!is.null(shapefile_path)) {
    if(!file.exists(shapefile_path)){
      stop("shapefile_path passed to vr_pull_link_table does not exist")
    }
    link_filepaths <- vr_convert_link_filepath(shapefile_path)

    if(!file.exists(link_filepaths$link_path)){
      stop("shapefile_path passed to vr_pull_link_table does not have corresponding link table")
    }

    link_table <- readRDS(link_filepaths$link_path)
    id_raster <- readRDS(link_filepaths$id_path)
  } else {
    shp_output <- vr_pull_shp(source, type, year, "filepath", link = T)
    link_table <- shp_output$link_table
    id_raster <- shp_output$id_raster
  }

  return(list("link_table" = link_table, "id_raster" = id_raster))
}

#######################################################################################
# Aggregation Functions



#######################################################################################
# Other Functions

#' @title vr_gen_link_tables
#'
#' @description Search through the VR directory to find and make link tables for all shapefiles without them. Saves them to link/ in source directory.
#'
#' @param core_repo string, filepath to lbd_core repo. No "/" at the end.
#' @param cores integer, The number of cores for `build_link_table()`, recommend > 30 for speed
#' @param source string, valid source or list of sources within prepped vr directory (typically in the form iso3_admlevel)
#'
#' @return NULL
#'
vr_gen_link_tables <- function(core_repo, cores, source = NULL) {
  source(paste0(core_repo, "/mbg_central/setup.R"))
  load_mbg_functions(core_repo)

  if(is.null(source)){
    source <- vr_list_sources()
  } else {
    if(any(!(source %in% vr_list_sources()))){
      stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
    }
  }

  #Find shapefiles that need link tables
  need_link <- c()
  for(i in source){
    base_path <- paste0(vr_get_base_path(), "extracts/", i)
    dir.create(file.path(base_path, "link"), showWarnings = FALSE)

    #get all shapefile names
    file_list <- vr_list_source_files(i, folder = T)
    shp_names <- grep("shapefile_stable/shapefile_stable.*shp$|shapefile_single_year/shapefile_single_year.*shp$", file_list, value = T)

    #drop folder and .shp extension from string
    shps <- grep("*.shp", unlist(strsplit(shp_names, split = "/")), value = T)

    #compare to link table folder and add any shps without corresponding files to need_link vector
    #shp_drop are indexes of shpapefiles to drop from shp_names
    shp_drop <- c()
    if(length(shps) > 0) {
      for(j in 1:length(shps)){
        link_files <- grep("link/link_table_shapefile_stable|link/link_table_shapefile_single_year", file_list, value = T)
        if(any(grepl(gsub(".shp", "", shps[j]), link_files))) {
          shp_drop <- c(shp_drop, j)
        }
      }
      shp_to_drop <- shp_names[shp_drop]
      shp_names <- setdiff(shp_names, shp_to_drop)
      if(length(shp_names) > 0){
        shp_names <- paste0(i, "/", shp_names)
      }
    }
    need_link <- c(need_link, shp_names)
  }

  #get full file path
  need_link <- paste0(vr_get_base_path(), "extracts/", need_link)

  #rasterize_check_coverage bug
  modeling_shapefile_version <<- NULL

  #loop through and save link tables for all shapefiles that need them
  for(i in need_link){

    if(i == '<< FILEPATH REDACTED >>>'){
      message("All shapefiles have link tables!")
      return(NULL)
    }

    message("Now building link table for shapefile:")
    message(i)

    output <- build_link_table(shapefile_version = NULL,
                               cores = cores,
                               region = NULL,
                               custom_shapefile_path = i,
                               custom_shapefile_field = "link_id")

    link <- output$link_table
    id_ras <- output$id_raster

    #get filepaths for link table and id raster
    link_filepaths <- vr_convert_link_filepath(i)

    saveRDS(link, link_filepaths$link_path)
    saveRDS(id_ras, link_filepaths$id_path)

    #clean up
    rm(link, id_ras)
    gc()

  }
}

#' @title vr_convert_link_filepath
#'
#' @description Given a shapefile path, return the paths for link table and id raster
#'
#' @param shapefile_path string, full file path to a shapefile in the VR directory
#'
#' @return Named list with the file paths to the link table and id raster (`link_path` and `id_path`)
#'
vr_convert_link_filepath <- function(shapefile_path){
  link_path <- gsub("shapefile_stable/shapefile_stable", "link/link_table_shapefile_stable", shapefile_path)
  link_path <- gsub("shapefile_single_year/shapefile_single_year", "link/link_table_shapefile_single_year", link_path)
  link_path <- gsub(".shp", ".RDS", link_path)
  id_path <- gsub("shapefile_stable/shapefile_stable", "link/id_raster_shapefile_stable", shapefile_path)
  id_path <- gsub("shapefile_single_year/shapefile_single_year", "link/id_raster_shapefile_single_year", id_path)
  id_path <- gsub(".shp", ".RDS", id_path)

  return(list("link_path" = link_path, "id_path" = id_path))
}


#' @title vr_get_population
#' @description: This function uses the frac_agg_covs function to return populations aggregated to
#' shapefile areal units.
#'
#'
#' @param ages Ages of population you wish to aggregate. Currently set up to take ages in five year age groupings,
#'             for example 0, 5, 10, ..., 65+.
#' @param sexes Sex id's for population you wish to aggregate
#' @param years Years of population estimates
#' @param core_repo LBD core repo used to load central functions
#' @param raked Should population be raked to GBD population estimates?
#' @param shapefile_path The filepath to the custom shapefile
#' @param shapfile_field The field in the shapefile that indicates the areal units
#' @param cores How many cores should be used in parallel processing
#'
#' @return Returns either  adata table of population by area, year, age and sex (if raked = F) or
#'         a list containing the raked population by area, year, age and sex and a table of raking factors by year, age and sex.

vr_get_population <- function(ages, sexes, years,
                              core_repo = "<<<< FILEPATH REDACTED >>>>",
                              raked = T,
                              shapefile_path,
                              shapefile_field,
                              cores = 1) {

  # Define functions used
  library(assertthat)
  library(purrr)
  library(tidyr)

  age_sex <- as.list(data.frame(t(expand.grid(ages, sexes))))
  # Define function to convert to worldpop measures
  age_year_to_worldpop_measure <- function(age, sex){
    age_end <- age + 4
    if (nchar(age) == 1) age_start <- paste0(0, age) else age_start <- age
    if (nchar(age_end) == 1) age_end <- paste0(0, age_end)
    measure = paste0("a", age_start, age_end)
    if (age == 65) measure <- "a65pl"
    if (sex == 1)  measure = paste0(measure, "m") else measure = paste0(measure, "f")
    return(measure)
  }
  measures <- map_chr(age_sex, ~ age_year_to_worldpop_measure(.[1], .[2]))
  attributes(measures) <- NULL
  message("Aggregating unraked population estimates")
  pop <- frac_agg_covs(covs = rep("worldpop", length(measures)),
                       years, measures, shapefile_path, shapefile_field,
                       core_repo = core_repo,
                       agg_method = "sum",
                       cores = 5)

  pop_unraked <-
    pop %>%
    dplyr::rename(area = shapefile_field) %>%
    tidyr::gather("age_sex", "pop", -area, -year) %>%
    dplyr::mutate(sex = ifelse(str_sub(age_sex, -1) == "m", 1, 2)) %>%
    dplyr::mutate(age_sex = str_remove(age_sex, "f|m")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(age = ifelse(age_sex != "a65pl", as.numeric(str_sub(age_sex, -2)) - 4, 65)) %>%
    dplyr::ungroup() %>%
    dplyr::select(area, year, age, sex, pop)

  if (raked == F) {
    message("Returning unraked population by age, sex and location")
    return(pop_unraked)
  } else {
    message("Raking to GBD population")
    # Get age link for GBD and age groupings for SAE
    source(paste0(core_repo, "mbg_central/post_estimation_functions.R"))
    source(paste0(CC_ENV_DIR, "/get_age_metadata.R"))
    age_link <-
      get_age_metadata(age_group_set_id = 12, gbd_round_id = 5) %>%
      dplyr::select(age_group_id, age_group_years_start) %>%
      dplyr::rename(age = age_group_years_start)
    gbd_age_groups <- age_link$age_group_id

    # Get country location code
    source(paste0(core_repo, "mbg_central/util/location_metadata_functions.R"))
    source(paste0(core_repo, "mbg_central/shapefile_functions.R"))
    loc_id <-
      get_location_code_mapping("current") %>%
      dplyr::filter(ihme_lc_id == toupper(country)) %>%
      dplyr::pull(loc_id)

    # Get GBD population estimates per age/sex/year for country
    source("<<<< FILEPATH REDACTED >>>>")
    gbd_pop <- get_population(age_group_id = gbd_age_groups,
                              location_id = loc_id,
                              sex = sexes,
                              year_id = years,
                              gbd_round_id = 5)

    # Group into appropriate age groups
    gbd_pop <-
      gbd_pop %>%
      dplyr::left_join(age_link, by = "age_group_id") %>%
      dplyr::rename(year = year_id, sex = sex_id, pop = population) %>%
      dplyr::select(age, year, sex, pop) %>%
      dplyr::mutate(age = case_when(age < 5 ~ 0,
                                    between(age, 5, 65) ~ age,
                                    age > 65 ~ 65)) %>%
      dplyr::group_by(age, year, sex) %>%
      dplyr::summarize(gbd_pop = sum(pop)) %>%
      dplyr::ungroup()

    # generate raking factors
    pop_rf <-
      pop_unraked %>%
      dplyr::group_by(year, age, sex) %>%
      dplyr::summarize(pop = sum(pop)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(gbd_pop, by = c("year", "age", "sex")) %>%
      dplyr::mutate(rf = gbd_pop / pop) %>%
      dplyr::select(year, age, sex, rf)

    # Apply linear raking factor
    pop <-
      pop_unraked %>%
      dplyr::left_join(pop_rf, by = c("year", "age", "sex")) %>%
      dplyr::mutate(pop = pop * rf) %>%
      dplyr::select(area, year, age, sex, pop) %>%
      data.table()

    # Make sure dataframe is square, no missing values
    assert_that(dplyr::select(pop,-pop) %>%
                  dplyr::summarise_all(funs(length(unique(.)))) %>%
                  prod() == nrow(pop))

    message("Returning raked population and raking factors")
    return(list(raked_pop = pop, rf = pop_rf))
  }
}
