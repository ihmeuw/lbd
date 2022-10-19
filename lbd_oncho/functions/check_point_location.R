###
### CHECK POINTS IN ADM0
###
### Purpose: Take in a csv of lat/longs and check if the points are inside the designated ADM0 code
###################################################################################################

### Needed Packages ### ----
# library(data.table)
# library(sf)
# library(rgeos)
# library(sp)

### Example Call ### ----
# checkPointLoc(<<<< FILEPATH REDACTED >>>>
#               , <<<< FILEPATH REDACTED >>>>
#               , "Latitude", "Longitude", "Country", "None")

### Function ######################################################################################
### Parameters ### ----
#'@csv - full file path to the csv or data.table/data.frame
#'@o_dir - full directory file path to location to save file
#'@lat_c - name of column in the data that is the latitude column
#'@long_c - name of column in the data that is the longitude column
#'@adm0_Nc - name of column in the data that is the ADM0_Name column, if it exists, not required
#'@adm0_Cc - name of column in the data that is the ADM0_Name column, if it exists, not required
### Notes ### ----
# - You need at least one of these two parameters, not necessarily both: adm0_Nc, adm0_Cc
# - The csv variable can be file location or an already loaded data.table/data.frame    

### Function ### ----
checkPointLoc <- function(csv, o_dir, lat_c, long_c, adm0_Nc = "None", adm0_Cc = "None") {
  
  ### Prep User Provided Data
  if (typeof(csv) == "character") {
    data <- as.data.table(fread(csv))
  } else {
    data <- as.data.table(csv)
  }
  # Rename columns
  setnames(data, old = c(lat_c, long_c), new = c("lat", "long"))
  if (adm0_Cc != "None") { setnames(data, old = c(adm0_Cc), new = c("ADM0_Code")) }
  if (adm0_Nc != "None") { setnames(data, old = c(adm0_Nc), new = c("ADM0_Name")) }
  # Pull shapefile
  shp_f <- <<<< FILEPATH REDACTED >>>>
  shp_all <- st_read(shp_f)
  shp_geom <- st_geometry(shp_all)
  
  # If code doesn't exist, need make it exist
  if (adm0_Cc == "None") {
    codes <- data.table(code = shp_all$ADM0_CODE, name = shp_all$ADM0_NAME)
    d_names <- unique(data$ADM0_Name)
    for (c_name in d_names) {
      # Manually checking for some names that don't exist in the ADM0 file
      if (c_name == "Cote d'Ivoire") {
        new_code <- 45 # Côte d'Ivoire
      } else if (c_name %in% c("Congo (Brazzaville)", "Congo (Kinshasa)")) {
        new_code <- 47 # Democratic Republic of Congo
      } else {
        new_code <- codes[name == c_name, ]$code
      }
      data[ADM0_Name == c_name, ADM0_Code := new_code ]
    }
  }
  
  
  ### Checking each ISO3 for points inside ### ----
  
  # Pull relevant columns, group, and unique them, only want to plot one point at a time 
  dataLoc <- data[, .(lat, long, ADM0_Code)]
  dataLoc <- dataLoc[, .(group = .GRP), by = list(lat, long, ADM0_Code)]
  dataLoc <- as.data.table(unique(dataLoc))
  
  # Prep empty table
  misplacedLatLong <- data.table()
  
  ### looping through all and marking ones that are outside
  for (code in unique(dataLoc$ADM0_Code)) {
    
    # Prep Shape
    c_shp <- shp_geom[which(shp_all$ADM0_CODE == code)]
    c_sp <- as_Spatial(c_shp)
    
    # Prep Points
    c_data <- dataLoc[ADM0_Code == code, ]
    coordinates(c_data) <- ~ long + lat
    
    # Checking
    proj4string(c_data) <- proj4string(c_sp)
    inside_sp <- over(c_sp, c_data, returnList = T)
    inside_sp <- inside_sp[[1]]$group
    
    # Pull ones outside
    c_data2 <- dataLoc[ADM0_Code == code, ]
    c_outside <- c_data2[!(group %in% inside_sp)]
    
    # Add any outside points to one data.table
    if (nrow(c_outside) > 0) {
      if (nrow(misplacedLatLong) == 0) {
        misplacedLatLong <- c_outside
      } else {
        misplacedLatLong <- rbindlist(list(misplacedLatLong, c_outside))
      }
    }
  }
  
  ### Mark in the data provided the points that are not in specified country ### ----
  
  data[, outside := F]
  for (row in 1:nrow(misplacedLatLong)) {
    c_row <- misplacedLatLong[row]
    data[lat == c_row$lat & long == c_row$long, outside := T]
  }
  
  ### Finish up ### ----
  
  message(paste0("There are ", nrow(data[outside == T, ]), "/", nrow(data), " ", "("
                 , round(nrow(data[outside == T, ])/nrow(data) *100, 2)  , "%) points outside of specified ADMIN0"))
  
  # Return names to original
  setnames(data, old = c("lat", "long"), new = c(lat_c, long_c))
  if (adm0_Cc != "None") { setnames(data, old = c("ADM0_Code"), new = c(adm0_Cc)) }
  if (adm0_Nc != "None") { setnames(data, old = c("ADM0_Name"), new = c(adm0_Nc)) }
  
  # Save and return shapefile if user wants it 
  save_f <- paste0(o_dir, "DATA_pointsMarkedOutsideADM0.csv")
  message(paste0("Saving full data with marked problems here: \n...", save_f))
  fwrite(data, save_f)
  return(data)
}








