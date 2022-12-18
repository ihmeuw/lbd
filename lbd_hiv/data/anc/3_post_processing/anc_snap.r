## ########################################################################################################################################
##
## ANC SNAP FUNCTIONS
##
## 
## Purpose: Identifies points in codebook that are <10km away from country borders
##          and snaps them into country.
## 
## Sub Functions:
##   format_codebook
##   snap_calculate_points
##
## Principle Function:
##   snap_codebook
##
## ########################################################################################################################################

## SUB FUNCTIONS ---------------------------------------------------------------------------------------------------------------------------
format_codebook <- function(codebook){
  point.locations <- subset(codebook, point == 1)
  point.locations$iso3 <- as.character(point.locations$iso3)
  point.locations$latitude <- as.numeric(as.character(point.locations$latitude))
  point.locations$longitude <- as.numeric(as.character(point.locations$longitude))
  point.locations <- subset(point.locations, !is.na(latitude) & !is.na(longitude))
  return(point.locations)
}

# Identifies out-of-bound points and returns the distance from polygon and nearest point on polygon boundary

snap_calculate_points <- function(codebook, admin0){
  # List countries in codebook
  countries <- sort(unique(codebook$iso3))
  message("Processing ", length(countries), " countries.")
  out_of_bounds <- data.frame()
  # Loop through countries and find points outside of national boundaries
  for (i in 1:length(countries)){
    country <- countries[i]
    
    #subset and prep survey to country
    sub <- unique(codebook[codebook$iso3 == country, ])
    
    #subset shapefile to survey country
    sub.shapefile <- admin0[admin0@data$ISO3 == sub$iso3[1], ] # ihme_lc_id for gbd shpfile
    
    # get points outside boundaries
    outside_pts <- sub[is.na(over(SpatialPoints(sub[, c("longitude", "latitude")], 
                                                CRS(proj4string(sub.shapefile))), 
                                                as(sub.shapefile, "SpatialPolygons"))),]
    sub_CRS <-  CRS(proj4string(sub.shapefile))
    
    # calculate the distance between point and polygon and note the nearest point on the polygon
    if (nrow(outside_pts) != 0) {
      sp <- SpatialPoints(outside_pts[,c("longitude", "latitude")], sub_CRS)
      nearest_point <- dist2Line(sp, as(sub.shapefile, "SpatialPolygons"))
      outside_pts$min_distance_km <- nearest_point[,1]/1000
      outside_pts$target_lat <- nearest_point[,3]
      outside_pts$target_long <- nearest_point[,2]
      
      out_of_bounds <- rbind(out_of_bounds, outside_pts)
    }
  }
  return(out_of_bounds)
}

## PRIMARY FUNCTION -------------------------------------------------------------------------------------------------------------------

# For a given codebook, identifies out-of-bound points <10km from border 
# and changes the lat/long to the nearest point 1km in the country boundary

snap_codebook <- function(codebook) {
  admin0 <- readOGR(FILEPATH)
  codebook$ID <- seq.int(nrow(codebook))
  point.locations <- format_codebook(codebook)
  not_point.locations <- codebook[!(codebook$ID %in% point.locations$ID),]
  bad_points <- subset(point.locations, latitude >90 | latitude < -90 | longitude >180 | longitude < -180)
  bad_iso3 <- subset(point.locations, !(iso3 %in% admin0@data$ISO3))
  point.locations <- subset(point.locations, !(latitude >90 | latitude < -90 | longitude >180 | longitude < -180) 
                                              & (iso3 %in% admin0@data$ISO3))
  to_snap <- point.locations[,c('ID','iso3', 'latitude','longitude')]
  out_of_bounds_over_10 <- data.frame()
  
  #repeat snap until no points are out_of_bounds (up to 10 iterations - 10 iterations were sufficient in testing)
  for (i in 1:10){
    if (dim(to_snap)[1] > 0){
      message("Iteration ", i)
      out_of_bounds_snapped <- snap_calculate_points(to_snap, admin0)
      message(dim(out_of_bounds_snapped)[1], " points out of bounds.")
      if (dim(out_of_bounds_snapped)[1] > 0) {
        if (i == 1){
          out_of_bounds_over_10 <- subset(out_of_bounds_snapped, min_distance_km > 10)
        }
        out_of_bounds_snapped <- subset(out_of_bounds_snapped, min_distance_km <= 10)
        if (dim(out_of_bounds_snapped)[1] > 0) {
          # new lat/long calculated by finding the point 1km inside country boundary along line between original lat/long and target lat/long
          out_of_bounds_snapped$new_lat <- (1/out_of_bounds_snapped$min_distance_km*(out_of_bounds_snapped$target_lat - out_of_bounds_snapped$lat)) + out_of_bounds_snapped$target_lat
          out_of_bounds_snapped$new_long <- (1/out_of_bounds_snapped$min_distance_km*(out_of_bounds_snapped$target_long - out_of_bounds_snapped$long)) + out_of_bounds_snapped$target_long
          
          # the calculated new lat/longs replace the lat/long in point.locations
          point.locations$latitude[match(out_of_bounds_snapped$ID, point.locations$ID)] <- out_of_bounds_snapped$new_lat
          point.locations$longitude[match(out_of_bounds_snapped$ID, point.locations$ID)] <- out_of_bounds_snapped$new_long
          
          out_of_bounds_snapped$latitude <- out_of_bounds_snapped$new_lat
          out_of_bounds_snapped$longitude <- out_of_bounds_snapped$new_long
        }
      } 
      to_snap <- out_of_bounds_snapped 
    }
  }
  
  if (exists('to_snap') && dim(to_snap)[1] >0){
    message(paste0("There are still ", dim(to_snap)[1],"  out of bounds points!"))
    point.locations$snapped[match(to_snap$ID, point.locations$ID)] <- "Flag - still out of bounds"
  }
  if (exists('bad_points') && dim(bad_points)[1]>0){
    bad_points$snapped <- "Flag - bad points"
  }
  if (exists('bad_iso3') && dim(bad_iso3)[1]>0){
    bad_iso3$snapped <- "Flag - bad Iso3"
  }
  if (exists('out_of_bounds_over_10') && dim(out_of_bounds_over_10)[1]>0){
    point.locations$snapped[point.locations$ID %in% out_of_bounds_over_10$ID] <- "Over 10 km away"
  } 
  
  final_codebook <- rbind.fill(list(point.locations,not_point.locations, bad_points, bad_iso3))
  final_codebook <- final_codebook[order(final_codebook$ID),]
  final_codebook$ID <- NULL
  return (final_codebook)
}



