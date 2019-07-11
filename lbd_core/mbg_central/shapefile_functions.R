#####################################################################
## Functions relating to offical world shapefiles
#####################################################################

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

    base_dir <- "<<<< FILEPATH REDACTED >>>>"

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
        stop(paste("Could not locate admin shapefile (", path, ")"))
    }
    return(path)
}

## detect_adm_shapefile_date_type() ################################################

#' Detects whether the admin shapefile you are using is gaul or gadm
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
  
  transition.date <- "<<<< FILEPATH REDACTED >>>>"
  gaul.exceptions <- "<<<< FILEPATH REDACTED >>>>"
  gadm.exceptions <- "<<<< FILEPATH REDACTED >>>>"

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
