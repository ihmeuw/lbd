## ##############################################################
## Inputs: Extracted, non-geomatched DTAs from Ubcov for BH
## Outputs: Geomatched BH dataset
## Description: This script combines extracted .DTA files 
##    and geomatches them to admin 1 and 2 names using shapefiles
## Dependencies: data.table, haven, plyr, raster, rgdal
## ##############################################################

################################################################
######################### Setup ################################
################################################################
assign("<<<< FILEPATH REDACTED >>>>", "<<<< FILEPATH REDACTED >>>>", envir = .GlobalEnv)
if (Sys.info()["sysname"] != "Linux") stop("This function should be run on the cluster only.")

# Load packages
message("\nLoading packages...")
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                         paste0(j_root,'<<<< FILEPATH REDACTED >>>>'),
                         paste0(j_root,'<<<< FILEPATH REDACTED >>>>'))

package_list <- c('data.table', 'plyr', "haven", "raster", "rgdal")
for(package in package_list) {
  library(package, lib.loc = package_lib, character.only=TRUE)
}

## Choose the working directory
setwd("<<<< FILEPATH REDACTED >>>>")

################################################################
#################### Combine Datasets ##########################
################################################################

## Ensure that all DTA or dta files are selected
file.list <- list.files(, pattern = '*dta$')
file.list2 <- list.files(, pattern = '*DTA$')
file.list <- c(file.list, file.list2)
rm(file.list2)

## Loop through and bind everything together
for(i in 1:length(file.list)) {
  file.name <- file.list[i]
  mydata <- data.table(read_dta(file.name))
  message(paste("processing file ", i, " of ", length(file.list) ))
  if(file.name == file.list[1]){
    combined <- mydata 
  } else{
    combined <- rbind(combined, mydata, fill=TRUE)
  }
}
saveRDS(combined, file = "<<<< FILEPATH REDACTED >>>>")

#####################################################################
#################### Merge Geographies ##############################
#####################################################################

# This following section is to attach geographies from the geography codebooks, and to grab administrative
# polygons for clusters with gps coordinates
# clear workspace

cbh_sbh <- data.frame(readRDS("<<<< FILEPATH REDACTED >>>>"))
rm(list =setdiff(ls(), cbh_sbh))

# read in geography codebook, combined
geographies <- read_dta("<<<< FILEPATH REDACTED >>>>")

# read in admin 1 shapefile
country_a1 <- shapefile('<<<< FILEPATH REDACTED >>>>')

# read in admin 2 shapefile
country_a2 <- shapefile('<<<< FILEPATH REDACTED >>>>')

# match geographies (for point records)
# create a match id, which is identical across the two datasets
cbh_sbh$match_id <- paste(cbh_sbh$nid, cbh_sbh$country, cbh_sbh$cluster_number, sep = "_")
geographies$match_id <- paste(geographies$nid, geographies$iso3, geographies$cluster_number, sep = "_")

cbh_sbh <- merge(cbh_sbh, geographies, by = "match_id")
cbh_sbh$match_id <- NULL

# create subsets based on point/polygon
cbh_sbh_na <- cbh_sbh[is.na(cbh_sbh$point),]
cbh_sbh_nonna <- cbh_sbh[!is.na(cbh_sbh$point),]

cbh_sbh_points <- cbh_sbh_nonna[cbh_sbh_nonna$point == 1, ]
cbh_sbh_non_point <- cbh_sbh_nonna[cbh_sbh_nonna$point == 0, ]

# sample points to grab admin 1 location names
# create matrix of point locations
pts <- cbh_sbh_points[c("longnum", "latnum")]

# using the sp package, where country is the admin 1 shapefile and pts 
# is a two-column (long, lat) matrix of coordinates:
pts_sp <- SpatialPoints(pts,
                        proj4string = CRS(projection(country_a1)))

# get admin1 name codes for all points
admin_1 <- over(pts_sp, country_a1)$NAME

# merge back with background dataset
cbh_sbh_points$admin_1 <- admin_1

# get admin2 names for all points
admin_2 <- over(pts_sp, country_a2)$NAME

# merge back with background dataset
cbh_sbh_points$admin_2 <- admin_2

# merge back macro dhs data
cbh_sbh_matched <- rbind(cbh_sbh_points,
                         cbh_sbh_non_point,
                         cbh_sbh_na)

#####################################################################
########################## Cleaning #################################
#####################################################################

# remove unnecessary cols
cbh_sbh_matched$latnum <- NULL
cbh_sbh_matched$longnum <- NULL
cbh_sbh_matched$point <- NULL
cbh_sbh_matched$location_name <- NULL
cbh_sbh_matched$admin_level <- NULL
cbh_sbh_matched$iso3 <- NULL
cbh_sbh_matched$buffer <- NULL
cbh_sbh_matched$location_code <- NULL
cbh_sbh_matched$shapefile <- NULL
cbh_sbh_matched$uncertain_point<- NULL
cbh_sbh_matched$caseid <- NULL

# correct some positional errors for Egypt surveys ('Administrative unit not available')
cbh_sbh_matched$admin_2[cbh_sbh_matched$admin_2 == "Administrative unit not available"] <- NA
cbh_sbh_matched$admin_2[cbh_sbh_matched$admin_2 == ""] <- NA

#fix some datatyping issues
setnames(cbh_sbh_matched, "nid.x", "nid")
setnames(cbh_sbh_matched, "cluster_number.x", "cluster_number")
cbh_sbh_matched$cluster_number.y <- NULL
cbh_sbh_matched$nid.y <- NULL
cbh_sbh_matched$cluster_number <- as.numeric(cbh_sbh_matched$cluster_number)
cbh_sbh_matched$household_number <- as.numeric(cbh_sbh_matched$household_number)

# write to disk
## Save it as an RDS file
cbh_sbh_matched <- data.table(cbh_sbh_matched)
saveRDS(cbh_sbh_matched, file = "<<<< FILEPATH REDACTED >>>>")


