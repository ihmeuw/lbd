## #############################################################################
##
## UTILITIES - MAPPING FUNCTIONS
##
## Purpose: Functions for preparing and creating different types of maps.
##
## #############################################################################

library(ggplot2)
library(data.table)
library(dplyr)
library(sf)
library(sp)


## fast_load_shapefile() ------------------------------------------------------>
#'
#' @title Fast Load a Shapefile
#'
#' @description Loads a shapefile using the `sf` library, which is faster than
#'   `readOGR()`. This function also optionally simplifies the polygon as an sf
#'   object.
#'
#' @author Nat Henry, \email{nathaniel.henry@st-hughs.ox.ac.uk}
#'
#' @param shp_path The full filepath to the shapefile, including the '.shp'
#'   extension
#' @param simplify_tol [default NULL] Optional argument giving the tolerance
#'   with which to simplify the input shapefile. In the units of the shapefile's
#'   default projection (for MBG, this is typically decimal degrees).
#'
#' @return A Spatial*DataFrame. For more information, see the `sp` library
#'
fast_load_shapefile <- function(shp_path, simplify_tol=NULL){
  # Load using sf
  shp_sf <- sf::st_read( dsn=shp_path )
  # Optionally simplify the polygons object
  if( !is.null(simplify_tol) ){
    message("Simplifying polygon...")
    shp_sf <- suppressWarnings(
      sf::st_simplify(
        shp_sf,
        preserveTopology=TRUE,
        dTolerance=simplify_tol
      )
    )
    message(" ... Simplification complete")
  }
  # Convert to SpatialPolygonsDataFrame and return
  shp_converted <- as(shp_sf, "Spatial")
  return(shp_converted)
}


## prep_shp_data_for_mapping() ------------------------------------------------>
#'
#' @title Prep shapefile and data for mapping
#'
#' @description Preps a shapefile for mapping in ggplot
#'
#' @author Nat Henry, \email{nathaniel.henry@st-hughs.ox.ac.uk}
#'
#' @param shp SpatialPolygonsDataFrame containing polygon boundaries
#' @param dataset Data.table of relevant data to visualize
#' @param merge_var Variable contained in both datasets to merge on
#' @param merge_var_type [default 'integer'] Type to force both merge variables
#'   to before merging. Should be one of 'integer' or 'character'
#'
#' @return A prepped data.table that can be input to ggplot and visualized using
#'   `geom_polygon()`
#'
prep_shp_data_for_mapping <- function(
  shp,
  dataset,
  merge_var,
  merge_var_type='integer'
  ){
  # Subset to only include observations that will eventually be included in the
  #  merge to speed up processing time
  # Data.tables and data.frames handle subsetting slightly differently
  mvt <- merge_var_type
  if('data.table' %in% class(dataset)){
    shp <- shp[as(shp@data[,merge_var],mvt) %in% as(dataset[,get(merge_var)],mvt),]  
  } else {
    shp <- shp[as(shp@data[,merge_var],mvt) %in% as(dataset[,merge_var],mvt),]  
  }
  # Assign shapefile IDs as a new field
  shp$id <- sapply(shp@polygons, function(x) x@ID) %>% as.numeric
  # Fortify the shapefile, then join all attributes back on using the ID
  shp_fort <- fortify(shp) %>% mutate(id=as.numeric(id)) %>% as.data.table
  setnames(shp_fort, 'order', 'shp_order')
  original_length <- nrow(shp_fort)
  shp_fort <- merge(
    x     = shp_fort,
    y     = shp,
    by    = c('id'),
    all.x = TRUE
  )
  # Merge on new data using the merge_var, ensuring that both variables are the
  #  same type
  dataset <- as.data.table(dataset)
  dataset[,  (merge_var) := as(get(merge_var), mvt) ]
  shp_fort[, (merge_var) := as(get(merge_var), mvt) ]
  shp_final <- merge(
    x     = shp_fort,
    y     = dataset,
    by    = merge_var,
    all.x = TRUE
  )
  # Sort by the 'shp_order' field
  shp_final <- shp_final[order(shp_order)]
  # Ensure that now rows have been added since the original fortify (which would
  #  indicate a bad merge)
  final_length <- nrow(shp_final)
  if(original_length != final_length) stop("Issue with the shapefile merge.")
  # Return the final shapefile for mapping
  return(shp_final)
}


## theme_map() ---------------------------------------------------------------->
#'
#' @title ggplot Map Theme
#'
#' @description A set of themes that can be used for making ggplot maps
#'
#' @author Nat Henry, \email{nathaniel.henry@st-hughs.ox.ac.uk}
#'
#' @details You can pass this function any arguments that would normally be
#'   passed to `ggplot2::theme()` to supplement or override existing settings.
#'
#' @return A ggplot theme object that can be added to ggplot calls to update
#'   the formatting
#'
theme_map <- function(...) {
  theme_minimal() +
    theme(
      text = element_text(color = "#22211d"),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(color = "#ebebe5", size = 0.2),
      panel.grid.minor = element_blank(),
      plot.background = element_rect(fill = "#FFFFFF", color = NA),
      panel.background = element_rect(fill = "#FFFFFF", color = NA),
      legend.background = element_rect(fill = "#FFFFFF", color = NA),
      panel.border = element_blank(),
      ...
    )
}


## coord_map_to_bounds() ------------------------------------------------------>
#'
#' @title Map to bounds
#'
#' @description Create a ggplot2 coordinate map using boundaries that fit the
#'   available dataset, with a buffer
#'
#' @param shp_fort The fortified shapefile that will be included in the ggplot
#'   object: specifically, a data.frame with numeric 'lat' and 'long' fields
#' @param projection [default 'mercator'] Which projection to use? For a full
#'   list, see `mapproj::mapproject()`.
#' @param buffer [default 0.5] Buffer (typically in decimal degrees) defining
#'   the limits of the map
#' @param ... You can add any extra arguments to `coord_map()` to supplement or
#'   override existing behavior of this function
#'
#' @return a `coord_map()` call that can be added to a ggplot object to define
#'   projections.
#'
coord_map_to_bounds <- function(
  shp_fort,
  projection='mercator',
  buffer=0.5,
  ...
  ){
  coord_map(
    projection=projection,
    xlim=c(min(shp_fort$long) - buffer, max(shp_fort$long) + buffer),
    ylim=c(min(shp_fort$lat)  - buffer, max(shp_fort$lat)  + buffer),
    ...
  )
}
