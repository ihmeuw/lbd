# Simple function to quickly load a shapefile
fast_load_shapefile <- function(shape,
                                fast_shapefile_dir = "<<<<FILEPATH REDACTED>>>>>") {

  wait_for_lock(shape)
  return(readRDS(paste0(fast_shapefile_dir, shape, ".rds")))
}

## get_admin_shape_dir
#' @title Return path to admin shapes directory
#'
#' @description
#' Returns path to the official LBD administrative shape file directory. This
#' actually includes non-shape data that is important for mapping, notably
#' the standard link table and standard id raster.
#'
#' @param version admin shapefile version to pull e.g., '2019_10_10'. Default: 'current'
get_admin_shape_dir <- function(version = "current") {
    paste0("<<<<FILEPATH REDACTED>>>>>", version, "/")
}

## get_admin_shapefile
#' @title Return path to world admin shapefile at specified admin level
#'
#' @description
#' Returns path to world admin shapefile given \code{admin_level}. stop()s if
#' no file exists at that admin_level. Defaults to returning the ".shp" file
#' path, but will substitute any other file \code{suffix} provided.
#'
#' @param admin_level Valid admin level we have a shapefile for. Current 0/1/2.
#'
#' @param suffix '.shp' by default, provide any other suffix to e.g., get the .dbf file
#' associated with the admin shapefile.
#'
#' @param raking boolean, default F. If TRUE pulls subnational raking shapefile.
#'
#' @examples
#' get_admin_shapefile(2)
#' get_admin_shapefile(2, suffix = '.dbf')
#'

get_admin_shapefile <- function(admin_level = 0, suffix = ".shp", type = "admin", version = 'current', raking = F) {
    if (raking) type = "raking"  # backwards compatibility

    base_dir <- get_admin_shape_dir(version)

    if (type == "admin" ) {
        path <- paste0(base_dir, "lbd_standard_admin_", admin_level, suffix)
    } else if (type == "raking") {
        path <- paste0(base_dir, "lbd_standard_raking", suffix)
    } else if (type == "disputed_mask") {
        path <- paste0(base_dir, "lbd_disputed_mask", suffix)
    } else {
        stop(paste("Unknown admin shapefile type '", type, "'"))
    }

    if (!file.exists(path)) {
        warning(sprintf("Could not locate admin shapefile (%s)", path))
    }
    return(path)
}

## is_admin_shapefile_string
#' @title Return logical indicating if string is an admin shapefile version string.
#'
#' @description
#' Checks strings to see if they are equal to "current" or match a YYYY_MM_DD format.
#' If so, returns TRUE. Else FALSE. No checking is done to ensure that the date string
#' has a corresponding subdirectory in the admin_shapefiles directory.
#'
#' @param s String the string to check.
#'
is_admin_shapefile_string <- function(s) {
  if (s == "current") {
      TRUE
  } else if (length(grep("^\\d{4}_\\d{2}_\\d{2}$", s))) {
      TRUE # "2019_02_27" or another admin shapefile release date
  } else {
    FALSE
  }
}



#' @title Detect admin shapefile date and type
#' @description Detects whether the admin shapefile you are using is gaul or gadm
#' and returns the version/date of the shapefile even if the 'current'
#' version was specified
#'
#' @param shpfile_path path to admin shapefile we want to learn about
#'
#' @return two element named list:
#'   1) list element 'shpfile_type' contains string: either 'gadm' or 'gaul'
#'   2) list element 'shpfile_date' contains the actual date of the
#'   shapefile, even if version='current' was used in shpfile_path
#'
#' @examples
#' detect_shapefile_type(get_admin_shapefile(version='current'))

detect_adm_shapefile_date_type <- function(shpfile_path = get_admin_shapefile(version = modeling_shapefile_version)){

  ## resolve the symlink to the path if version='current'
  if(grepl(pattern = 'current', x = shpfile_path)){
    resolve.command <- paste("readlink -f", shpfile_path)
    full.path <- system(resolve.command, intern = TRUE)
  }else{
    full.path <- shpfile_path ## this is faster than resolving if it's an option...
  }

  ## grab the date from the full path
  sf.date <- strsplit(full.path, '/')[[1]][6]

  ## determine if shpfile.date pertains to gaul or gadm

  ## assume versions dated on and after Sept. 1, 2018
  ## (transition.date) are GADM, while dates before transition.date
  ## are GAUL unless otherwise coded in as an exception
  transition.date <- "2018_09_01"
  gaul.exceptions <- c() ## dates after transition.date that are actually GAUL
  gadm.exceptions <- c('2018_08_01') ## dates before transition.date are actually GADMw

  ## determine gaul or gadm
  if(sf.date >= transition.date){
    sf.type <- 'gadm'
  } else{
    sf.type <- 'gaul'
  }

  if(sf.date %in% gaul.exceptions) sf.type <- 'gaul'
  if(sf.date %in% gadm.exceptions) sf.type <- 'gadm'

  return(list(shpfile_type = sf.type,
              shpfile_date = sf.date))
}
