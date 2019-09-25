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
get_admin_shapefile <- function(admin_level, suffix = ".shp", raking = F) {
    path <- paste0("<<<< FILEPATH REDACTED >>>>")
    if (!file.exists(path)) {
        stop(paste("Could not locate admin shapefile at admin_level", admin_level, "(", path, ")"))
    }
    if (raking) {
      path <- paste0("<<<< FILEPATH REDACTED >>>>")
    }
    return(path)
}
