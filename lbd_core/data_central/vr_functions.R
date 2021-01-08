################################# VR Functions ################################
#
# Functions associated with pulling prepped VR data and other related data
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
#' @param year string, with different options based on `type`. If `type` is `annual`, `year` can be "first", 
#' "last" or a year character, i.e. "2000". first and last will return the first and last available years respectively.
#' If `type` is `stable`, `year` can be "full", which returns the stable shapefile over the longest time period, 
#' or a vector of 2 year characters. When years are passed in, shapefiles with these years must be 
#' present in the shapefile_single_year (annual) or shapefile_stable (stable) folders.
#' @param return_object string, "sf", "sp", or "filepath". Determines what object to return.
#' @param check_exists boolean, default FALSE. Check that a filepath exists before
#'   returning? This option can only be set to TRUE if the `return_object` is 
#'   set as "filepath".
#' @param link boolean, return link objects as well?
#' @param admin_level int, admin level of shapefile to pull. If NULL will pull admin level specified in source.
#'
#' @return list object with 3 objects: `shapefile`: shapefile object specified by `return_object`, `link_table`: associated link table if `link` = T else NULL, and `id_ras`: associated id raster if `link` = T else NULL
#'
vr_pull_shp <- function(
  source, type, year, return_object="filepath", check_exists=T, link=T, admin_level=NULL
){
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
  if((check_exists==FALSE) & (return_object != 'filepath')){
    stop("check_exists can only be FALSE if the return_object is a 'filepath'.", call.=FALSE)
  }
  
  #finding the right file
  file_list <- vr_list_source_files(source)
  if(type == "annual"){
    shp_names <- grep("shapefile_single_year.*shp", file_list, value = T)
    shp_years <- as.numeric(str_extract(shp_names, "\\d{4}"))
    
    if(year == "first") {
      file_year <- as.character(min(shp_years, na.rm = T))
    } else if(year == "last") {
      file_year <- as.character(max(shp_years, na.rm = T))
    } else {
      file_year <- year
    }
    
    file_path <- paste0("shapefile_single_year_", file_year, ifelse(!is.null(admin_level), paste0("_admin", admin_level), ""), ".shp")
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
  if(check_exists){
    vr_check_file_exists(source, file_path)
  }
  
  #load/convert as specified by return_object param
  base_path <- paste0(vr_get_base_path(), "extracts/", source, if(type == "annual") "/shapefile_single_year/" else "/shapefile_stable/")
  if(return_object == "sp"){
    shp <- readOGR(paste0(base_path, file_path), verbose=FALSE)
  } else if (return_object == "sf"){
    shp <- st_read(paste0(base_path, file_path), quiet=TRUE)
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
vr_pull_u5m <- function(source, version = "current") {
  if(!(source %in% vr_list_sources())){
    stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
  }

  file_name <- if(version == "current") paste0("u5m_", source, "_current.csv") else paste0("u5m_", source, "_", version, ".csv")
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
  
  file_name <- if(version == "current") paste0("births_", source, "_current.csv") else paste0("births_", source, "_", version, ".csv")
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
#' @param year string, with different options based on `type`. If `type` is `annual`, `year` can be "first", 
#' "last" or a year character, i.e. "2000". first and last will return the first and last available years respectively.
#' If `type` is `stable`, `year` can be "full", which returns the stable shapefile over the longest time period, 
#' or a vector of 2 year characters. When years are passed in, shapefiles with these years must be 
#' present in the shapefile_single_year (annual) or shapefile_stable (stable) folders.
#' @param shapefile_path default `NULL`, can be used in place of the above 3 options to pull the link table. Takes a string filepath to the .shp file of the shapefile that you want the link table and id raster for.
#' @return list object with 2 objects: `link_table`: associated link table and `id_raster`: associated id raster
#'
#' @example 
#' using vr_pull_shp parameters:
#' link <- vr_pull_link_table(source = "chl_adm3",
#'                            type = "stable",
#'                            year = "full") 
#'                            
#' using shapefile path:
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


#' @title vr_pull_adjacency_matrix
#'
#' @description pull an adjacency matrix from the vr_prep directory. Requires either source/type/year 
#' ala vr_pull_shp, or shapefile_path. If shapefile_path is not NULL, it will be prioritized.
#' 
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param type string, either `stable` or `annual`
#' @param year string, with different options based on `type`. If `type` is `annual`, `year` can be "first", 
#' "last" or a year character, i.e. "2000". first and last will return the first and last available years respectively.
#' If `type` is `stable`, `year` can be "full", which returns the stable shapefile over the longest time period, 
#' or a vector of 2 year characters. When years are passed in, shapefiles with these years must be 
#' present in the shapefile_single_year (annual) or shapefile_stable (stable) folders.
#' @param shapefile_path default `NULL`, can be used in place of the above 3 options to pull the link table. Takes a string
#' filepath to the .shp file of the shapefile that you want the link table and id raster for.
#' @param suffix string to be added to the end of the file name for the adjacency matrix. Used to distinguish
#' between different types of adjacency matrices.
#' @return adjacency matrix
#'
#' @example 
#' using vr_pull_shp parameters:
#' link <- vr_pull_adjacency_matrix(source = "chl_adm3",
#'                                  type = "stable",
#'                                  year = "full",
#'                                  suffix = "default") 
#'                            
#' using shapefile path:
#'
#' 
vr_pull_adjacency_matrix <- function(source, type, year, shapefile_path = NULL, suffix = "default") {
  #build shapefile path if not passed in
  if(is.null(shapefile_path)) {
    shapefile_path <- vr_pull_shp(source, type, year, "filepath", link=F)$shapefile
  }
  
  #check shapefile exists
  if(!file.exists(shapefile_path)){
    stop("shapefile_path passed to vr_pull_link_table does not exist")
  }
  
  #convert filepath to adj_mat format
  adj_filepath <- gsub(".shp", paste0("_adj_mat_", suffix, ".RDS"), shapefile_path)
  
  #check adjacency matrix file exists
  if(!file.exists(adj_filepath)){
    stop("shapefile_path passed to vr_pull_link_table does not have corresponding link table")
  }

  adj_mat <- readRDS(adj_filepath)
  
  return(adj_mat)
}


#######################################################################################
# Object Generation Functions

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
    system(paste0("chmod -R 777 ", file.path(base_path, "link")))
    
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
    
    if(i == "<<<< FILEPATH REDACTED >>>>"){
      message("All shapefiles have link tables!")
      return(NULL)
    }
    
    if(grepl("admin", i, fixed=TRUE)) {
      custom_field <- "GAUL_CODE"
    } else {
      custom_field <- "link_id"
    }
    
    message("Now building link table for shapefile:")
    message(i)
    
    output <- build_link_table(shapefile_version = NULL,
                               cores = cores,
                               region = NULL,
                               custom_shapefile_path = i,
                               custom_shapefile_field = custom_field)
    
    link <- output$link_table
    id_ras <- output$id_raster
    
    if(is.character(link$link_id)) {
      link$link_id <- as.numeric(link$link_id)
    }
    
    #get filepaths for link table and id raster
    link_filepaths <- vr_convert_link_filepath(i)
    
    saveRDS(link, link_filepaths$link_path)
    saveRDS(id_ras, link_filepaths$id_path)
    
    #clean up
    rm(link, id_ras)
    gc()
  
  }
}


#' @title Generate adjacency matrices
#'
#' @description Search through the VR directory to find any shapefiles without adjacency matrices and build them.
#' The function used to produce the matrices can be passed to the function, with extra args passed in through "...". 
#' Suffix is used to distinguish between different types of adjacency matrices. The function uses the suffix to 
#' identify shapefiles needing matrices, i.e. if all shapefiles have a "default" matrix produced by 
#' \code{\link{vr_build_default_adjacency_matrix}}, but you want to generate adjacency matrices from a different 
#' function, changing the suffix will identify all the shapefiles as needing the new matrix.
#' 
#' @param core_repo string, filepath to lbd_core repo. No "/" at the end.
#' @param source string, valid source or list of sources within prepped vr directory (typically in the form iso3_admlevel). 
#' If NULL, will use all sources
#' @param overwrite boolean, if true, will remake the adjacency matrices for all shapefiles (skips check to see if they already exist)
#' @param adj_function string, the name of the function used to generate the adjacency matrix
#' @param suffix string to be added to the end of the file name for the adjacency matrix. Used to distinguish
#' between different types of adjacency matrices.
#' @param ... additional args to be passed into adj_function. The SpatialPolygonsDataFrame from the source shapefile will
#' always be the first arg passed to the adj function. These arguments must be passed in order of the args to the adj function
#' call because they get unnamed when passed to `do.call`.
#'
#' @return NULL
#' 
#' @examples
#' \dontrun{
#' 
#' #generate default matrices for all shapefiles that don't already have it
#' vr_gen_adjacency_matrices(core_repo)
#' 
#' #generate matrices for "bra_adm2" sources shapefiles using the W style in nb2mat
#' #(passed into vr_build_default_adjacency_matrix)
#' #saves out as shapefile_adj_mat_W.RDS
#' vr_gen_adjacency_matrices(core_repo, source = "bra_adm2", suffix = "W", style = "W", allow_zero_neighbors = TRUE)
#' 
#' }
#'
vr_gen_adjacency_matrices <- function(core_repo,
                                      source = NULL,
                                      overwrite = FALSE,
                                      adj_function = "vr_build_default_adjacency_matrix", 
                                      suffix = "default", 
                                      ...) {
  
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
  need_adj <- c()
  for(i in source){
    base_path <- paste0(vr_get_base_path(), "extracts/", i)
    
    #get all shapefile names
    file_list <- vr_list_source_files(i, folder = T)
    shp_names <- grep("shapefile_stable/shapefile_stable.*shp$|shapefile_single_year/shapefile_single_year.*shp$", file_list, value = T)
    
    #check to see if adj_mat file exists, if not, add to list of files to generate
    if(length(shp_names) > 0) {
      
      #drop folder and .shp extension from string
      shps <- gsub(".shp", "", shp_names)
      shps <- paste0(shps, "_adj_mat_", suffix, ".RDS")
      
      for(j in shps){
        shp_check <- paste0(base_path, "/", j)
        if(!overwrite){
          if(!file.exists(shp_check)) {
            need_adj <- c(need_adj, paste0(i, "/", j))
          }
        } else {
          need_adj <- c(need_adj, paste0(i, "/", j))
        }
      }
    }
  }
  
  #convert need_adj back to shapefile names
  need_adj <- gsub(paste0("_adj_mat_", suffix, ".RDS"), ".shp", need_adj)
  
  error_list <- list()
  message("Building adjacency matrices for -")
  #build and save adjacency matrices
  for(i in need_adj) {
    message("     ", i)
    temp_shp_path <- paste0(vr_get_base_path(), "extracts/", i)
    sp_shp <- readRDS(gsub(".shp", ".RDS", temp_shp_path))
    
    adj_mat_args <- list(sp_shp)
    
    #add in ... args to the adj_function call
    if(length(list(...)) > 0) {
      adj_mat_args <- unname(c(adj_mat_args, list(...)))
    }
    
    result <- 
      tryCatch(adj_mat <- do.call(adj_function, args = adj_mat_args),
               error = function(e) {
                 message("adjacency matrix for ", i, " failed to build, skipping.")
                 return(e)
               })
    if(!("error" %in% class(result))){
      vr_save_adjacency_matrix(adj_mat, shp_path = temp_shp_path, suffix = suffix)
    } else {
      temp_list <- list(result)
      names(temp_list) <- c(eval(i))
      error_list <- c(error_list, temp_list)
    }
  }
  return(error_list)
}


#######################################################################################
# Other Functions

#' @title vr_get_population
#' @description: This function uses the frac_agg_covs function to return populations aggregated to 
#' shapefile areal units for VR data
#' 
#' @param ages Ages of population you wish to aggregate. Set up to take ages in five year age groupings, 
#'             for example 0, 5, 10, ..., 80+, currently only allows ages 0, 5, ..., 80+
#' @param sexes Sex id's for population you wish to aggregate
#' @param years Years of population estimates. If passing multiple years, must be 
#'              a contiguous series of years.
#' @param country VR country you are running vr_get_population for, specified as an ISO3 code.
#' @param admin Admin unit of VR data, the same level as that in the shapefile. Entered
#'              as a string, such as 'adm1' or 'adm2'
#' @param worldpop_release Release of the worldpop being used 
#' @param core_repo LBD core repo used to load central functions
#' @param raked Should population be raked to GBD population estimates? 
#' @param gbd_round_id If raking, the version of GBD to rake against 
#' @param shapefile_path The filepath to the custom shapefile
#' @param shapefile_field The field in the shapefile that indicates the areal units
#' @param cores How many cores should be used in parallel processing
#' 
#' @return Returns either  adata table of population by area, year, age and sex (if raked = F) or
#'         a list containing the raked population by area, year, age and sex and a table of raking factors by year, age and sex. 
vr_get_population <- function(ages = seq(0, 80, 5),
                              sexes = c(1, 2),
                              years,
                              country = NULL,
                              admin = NULL,
                              worldpop_release = "2020_03_20",
                              core_repo = paste0("<<<< FILEPATH REDACTED >>>>/lbd_core/"), 
                              raked = TRUE,
                              gbd_round_id = 5,
                              shapefile_field = "GAUL_CODE",
                              shapefile_path  =  NULL, 
                              link_table = NULL,
                              id_raster = NULL,
                              cores = 1) {
  
  # Load library
  library(purrr)
  library(tidyr)
  library(doParallel)
  
  # Source scripts for central functions
  source(paste0(core_repo, "/mbg_central/post_estimation_functions.R"))
  source(paste0(core_repo, "/mbg_central/prep_functions.R"))
  source(paste0(core_repo, "/mbg_central/shapefile_functions.R"))
  source(paste0(core_repo, "/data_central/sae_functions.R"))
  source("<<<< FILEPATH REDACTED >>>>/get_age_metadata.R")
  source("<<<< FILEPATH REDACTED >>>>/get_population.R")
  
  # Define age measures and sex 
  country <- tolower(country)
  if (is.null(shapefile_path)) shapefile_path <- vr_pull_shp(paste0(country, "_", admin), "stable", "full")$shapefile
  age_sex <- as.list(data.frame(t(expand.grid(ages, sexes))))
  measures <- map_chr(age_sex, ~ vr_age_year_to_worldpop_measure(.[1], .[2]))
  attributes(measures) <- NULL
  
  # Logic checks for raking to GBD; creates map from modeled small areas to gbd location id(s)
  if (raked) {
    if (is.null(country) | is.null(admin)) {
      stop("Country name and admin level needed for pulling location metadata to rake population data, see documentation")
    } else {
      # Grab GBD location id(s)
      gbd_loc_ids <- get_gbd_locs(country, shapefile_version = "current", rake_subnational = T)$location_id
      
      # Grab link between modeling areas and GBD location id(s)
      if (shapefile_field == "uid") {
        area_to_loc_id <-
          vr_pull_loc_meta(paste0(country, "_", admin)) %>%
          filter(level == as.numeric(str_extract(admin, "[0-9]+")),
                 year_end == max(year_end)) %>%
          dplyr::select(area = uid, loc_id = location_parent_id) %>%
          arrange(area, loc_id) %>% 
          unique()
      } else if (shapefile_field == "GAUL_CODE") {
        area_to_loc_id <-
          vr_pull_loc_meta(paste0(country, "_", admin)) %>%
          filter(level == as.numeric(str_extract(admin, "[0-9]+")),
                 year_end == max(year_end)) %>%
          dplyr::select(area = location_id, loc_id = location_parent_id) %>%
          arrange(area, loc_id) %>% 
          unique()
      } else {
        stop("This function was only designed for aggregating uid or GAUL_CODE (admin2 code)")
      }
      
      if (length(gbd_loc_ids) == 1) {
       
         # If only one location id, then GBD is national only
        message("Raking to national GBD Population")
        area_to_loc_id <- unique(mutate(area_to_loc_id, loc_id = gbd_loc_ids))
        rownames(area_to_loc_id) <- NULL # Row numbers can get confusing

      } else {
        message("Raking to subnational GBD Population")
        
        # Make sure area_to_loc_ids contains the correct GBD subnational loc ids
        if (!isTRUE(all.equal(sort(gbd_loc_ids), sort(unique(area_to_loc_id$loc_id))))) {
          stop ("GBD loc ids are not the same as those present in location metadata. Check on VR function vr_pull_loc_meta")
        }
      }
    }
  }

  # Define covariate config for fractionally aggregating covariates
  message("Aggregating unraked population estimates")
  covariate_config <- 
    data.table(covariate = rep("worldpop", length(measures)),
               measure = measures,
               release = rep(worldpop_release, length(measures)),
               agg_method = rep("sum", length(measures)))
  
  if (!is.null(link_table)) {
    message("Building custom link table")
    link_table_list <- build_link_table(shapefile_version = NULL,
                                        cores = cores,
                                        region = NULL,
                                        custom_shapefile_path = shapefile_path,
                                        custom_shapefile_field = shapefile_field)
    link_table <- link_table_list$link_table
    id_raster  <- link_table_list$id_raster
  }
  
  pop <- 
    frac_agg_covs(cov_config = covariate_config,
                  years = years,
                  shapefile_path = shapefile_path,
                  shapefile_field = shapefile_field,
                  core_repo = core_repo,
                  link_table = link_table,
                  id_raster = id_raster,
                  cores = cores)
  
  # Rename dataset to area, year, age, and sex
  pop_unraked <- 
    pop %>% 
    dplyr::rename(area = shapefile_field) %>% 
    tidyr::gather("age_sex", "pop", -area, -year) %>% 
    dplyr::mutate(sex = ifelse(str_sub(age_sex, -1) == "m", 1, ifelse(str_sub(age_sex, -1) == "f", 2, 3))) %>% 
    dplyr::mutate(age_sex = str_remove(age_sex, "f|m|t")) %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(age = ifelse(age_sex != "a80pl", as.numeric(str_sub(age_sex, -2)) - 4, 80)) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(area, year, age, sex, pop) %>% 
    data.table()
  
  if (!raked) {
    message("Returning unraked population by age, sex and location")
    return(pop_unraked)
    
  } else {
    message("Raking to GBD population")
    
    # Link between GBD age group and SAE ages
    age_link <- 
      get_age_metadata(age_group_set_id = 12, gbd_round_id = gbd_round_id) %>% 
      filter(between(age_group_years_start, 5, 75)) %>% 
      dplyr::select(age_group_id, age_group_years_start) %>% 
      dplyr::rename(age = age_group_years_start)
    
    # Grab GBD ids for under 5 and over 80 year olds
    under5 <- data.table(age_group_id = 1, age = 0)
    over80 <- data.table(age_group_id = 21, age = 80)
    age_link <- rbind(under5, age_link, over80)
    
    # Get GBD population estimates per age/sex/year for raking level
    gbd_pop <- 
      get_population(age_group_id = age_link$age_group_id,
                     location_id = gbd_loc_ids,
                     sex = sexes,
                     year_id = years,
                     gbd_round_id = gbd_round_id,
                     #decomp_step = 'step4',
                     status = "best") %>% 
      dplyr::left_join(age_link, by = "age_group_id") %>%
      dplyr::select(age, year = year_id, sex = sex_id, loc_id = location_id, gbd_pop = population)
 
    # generate raking factors 
    rf <- 
      pop_unraked %>%
      left_join(area_to_loc_id, by = "area") %>% 
      dplyr::group_by(year, age, sex, loc_id) %>%
      dplyr::summarize(pop = sum(pop)) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(gbd_pop, by = c("year", "age", "sex", "loc_id")) %>%
      dplyr::mutate(rf = gbd_pop / pop) %>%
      dplyr::select(year, age, sex, loc_id, rf) %>% 
      data.table()
    
    # Apply linear raking factor 
    pop_raked <- 
      pop_unraked %>% 
      dplyr::left_join(area_to_loc_id, by = "area") %>% 
      dplyr::left_join(rf, by = c("year", "age", "sex", "loc_id")) %>% 
      dplyr::mutate(pop = pop * rf) %>% 
      dplyr::select(area, year, age, sex, pop) %>% 
      data.table()
    
    # Make sure dataframe is square
    assert_that(dplyr::select(pop_raked, -pop) %>%
                  dplyr::summarise_all(funs(length(unique(.)))) %>%
                  prod() == nrow(pop_raked))
    
    message("Returning raked population and raking factors")
    return(list(pop_raked = pop_raked, rf = rf))
  }
}


#######################################################################################
# Helper Functions

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


#' @title vr_age_year_to_worldpop_measure
#' @description: Given age in 5 year intervals and sex in numeric, return worldpop measure for population
#' 
#' @param ages Ages of population 5 year intervals (0 represents 0-4, 5 is 5-9, etc)
#' @param sex  Sex of population (1 = men, 2 = women, 3 = both)
#' @return Worldpop measure that can be used to load worldpop rasters of appropriate age and sex
#' 
vr_age_year_to_worldpop_measure <- function(age, sex){
  age_end <- age + 4
  if (nchar(age) == 1) age_start <- paste0(0, age) else age_start <- age
  if (nchar(age_end) == 1) age_end <- paste0(0, age_end)
  measure = paste0("a", age_start, age_end)
  if (age == 80) measure <- "a80pl"
  if (sex == 1) {
    measure = paste0(measure, "m")
  } else if (sex == 2) {
    measure = paste0(measure, "f")
  } else {
    measure = paste0(measure, "t")
  }
  return(measure)
}


#' @title Build default adjacency matrix
#' 
#' @description Builds a matrix describing adjacency between polygons in the 
#'   a spatial object. This function is a wrapper of two functions in the
#'   \code{\link{spdep}} package.
#' 
#' @param poly_sp A SpatialPolygonsDataFrame object
#' @param style One of "B","W","C", or "S". For more information, see the 
#'   documentation in the \code{\link{spdep}} package:
#'   https://www.rdocumentation.org/packages/spdep/versions/1.1-3/topics/nb2mat
#' 
#' @return sparse dsCMatrix representing adjacency between polygons in the
#'   `poly_sp` object
#' 
vr_build_default_adjacency_matrix <- function(poly_sp, style='B', allow_zero_neighbors=TRUE){
  ## Generate adjacency matrix for the polygon
  adjmat <- spdep::nb2mat(
    neighbours = spdep::poly2nb(poly_sp),
    style = style,
    zero.policy = allow_zero_neighbors
  )
}


#' @title Save adjacency matrix
#' 
#' @description Saves out an adjacency matrix for a given source/shapefile into the shapefile folder. 
#' Source and shapefile parameters are the standard format for pulling shapefiles
#' found in \code{\link{vr_pull_shp}}. If a function other than 
#' \code{\link{vr_build_default_adjacency_matrix}} is used to generate the adjacency
#' matrix, a suffix other than "default" can be used.
#' 
#' @param adj_mat An adjacency matrix generate by \code{\link{vr_build_default_adjacency_matrix}}
#' or another function.
#' @param source string, valid source within prepped vr directory (typically in the form iso3_admlevel)
#' @param type string, either `stable` or `annual`
#' @param year string, with different options based on `type`. If `type` is `annual`, `year` can be "first", 
#' "last" or a year character, i.e. "2000". first and last will return the first and last available years respectively.
#' If `type` is `stable`, `year` can be "full", which returns the stable shapefile over the longest time period, 
#' or a vector of 2 year characters. When years are passed in, shapefiles with these years must be 
#' present in the shapefile_single_year (annual) or shapefile_stable (stable) folders.
#' @shp_path string, path to the shapefile getting an adjacency matrix (used for situations where source/type/year are not readily available)
#' @param suffix string to be added to the end of the file name for the adjacency matrix. Used to distinguish
#' between different types of adjacency matrices.
#' 
#' @return NULL
#' 
vr_save_adjacency_matrix <- function(adj_mat,
                                     source, 
                                     type,
                                     year,
                                     shp_path = NULL,
                                     suffix = "default"){
  
  if(is.null(shp_path)) {
    shp_path <- vr_pull_shp(source, type, year, link = F)$shapefile
  }
  adj_path <- gsub(".shp", "", shp_path)
  adj_path <- paste0(adj_path, "_adj_mat_", suffix, ".RDS")
  saveRDS(adj_mat, adj_path)
}


## vr_build_parent_shapefiles ------------------------------------------------->
#' 
#' @title Build VR parent shapefiles
#' @description Given a single-year shapefile corresponding to the most granular
#'   administrative divisions available for VR, build a set of shapefiles
#'   containing the boundaries of the parent administrative divisions up to the
#'   admin0 level. Each shapefile polygon can be uniquely identified by the
#'   GAUL_CODE field, which contains the GBD location IDs for each set of
#'   administrative divisions
#' 
#' @param source [char] valid source within prepped vr directory (typically 
#'   in the form "<ISO3>_<ADMLEVEL>")
#' @param verbose [bool, default TRUE] Should progress messages from this
#'   function be written to screen?
#' 
vr_build_parent_shapefiles <- function(source, verbose=TRUE){
  # Set up logging function
  vbmsg <- function(msg) if(verbose) message(msg)
  # Check that the source exists
  if(any(!(source %in% vr_list_sources()))){
    stop("source is not valid - use vr_list_sources() to find available sources", call. = F)
  }
  # Get full administrative metadata for the source
  loc_meta <- vr_pull_loc_meta(source=source, version='current')
  # Get source max admin level
  max_level <- max(loc_meta$level)

  # Get years for all valid single-year shapefiles in the VR directory
  all_files <- vr_list_source_files(source)
  single_year_shapes <- grep(
    'shapefile_single_year_[0-9]+.shp', all_files, value=TRUE
  )
  shp_years <- gsub("^shapefile_single_year_","",single_year_shapes) %>%
    gsub(".shp$","",.) %>%
    as.numeric
  # Message number of sources found
  if(length(shp_years)==0){
    vbmsg(glue::glue("No annual shapefiles found for source {source}."))
    invisible()
  }
  vbmsg(glue::glue(
    "Found {length(shp_years)} annual shapefiles for source {source} in years ",
    "{paste(shp_years, collapse=', ')}."
  ))
  for(this_year in shp_years){
    vbmsg(glue::glue("  Working on year {this_year}:"))
    this_shp <- vr_pull_shp(
      source=source, type='annual', year=this_year, return_object="sf", 
      check_exists=T, link=F
    )$shapefile
    for(this_level in rev(1:max_level)){
      vbmsg(glue::glue("   - Working on parent_shapefile for level {this_level}"))
      # Merge on parent information based on the location_id
      this_shp <- this_shp[, c('GAUL_CODE','geometry')]
      parent_info <- unique(copy(
        loc_meta[level==this_level, .(location_id, location_parent_id)]
      ))
      merged_shp <- merge(
        x=this_shp, y=parent_info, by.x='GAUL_CODE',by.y='location_id'
      )

      # Group by parent location to create parent shapefile
      parent_shp <- merged_shp %>% group_by(location_parent_id) %>% summarize()
      # Fix any obvious topology issues by buffering with zero distance
      # Trick sourced from https://www.r-spatial.org/r/2017/03/19/invalid.html
      parent_shp <- sf::st_buffer(parent_shp, dist=0.0) %>%
        suppressWarnings %>% suppressMessages
      # Rename the unique ID column of the parent shapefile to "GAUL_CODE"
      names(parent_shp)[names(parent_shp)=='location_parent_id'] <- "GAUL_CODE"

      # Save the parent shapefile using the filepath generated by `vr_pull_shp`
      parent_shp_path <- vr_pull_shp(
        source=source, type='annual', year=this_year, return_object="filepath", 
        check_exists=F, link=F, admin_level=(this_level-1)
      )$shapefile
      sf::st_write(
        obj=parent_shp, dsn=parent_shp_path, delete_dsn=TRUE, quiet=TRUE
      ) %>% suppressWarnings

      # Update so that the new working shapefile is at the parent level
      this_shp <- parent_shp
    }
    vbmsg(glue::glue("    All parent shapefiles created for year {this_year}."))
  }
  # Parent shapefiles have been created for all years
  # Optionally return completion message and finish
  vbmsg(glue::glue("All annual parent shapefiles created for source {source}."))
  invisible()
}
