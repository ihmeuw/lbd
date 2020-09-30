########################################################################
#Switch shapefile of Oxford counterfactual csvs for malaria 
#######################################################################

#(1) Setup ############################################################
rm(list = ls())

#pull arguments
old_shapefile_version <- commandArgs()[4]
new_shapefile_version <- commandArgs()[5]
tridaly_in_date <- commandArgs()[6]
tridaly_out_date <- commandArgs()[7]
pop_release <- commandArgs()[8]
indicator <- commandArgs()[9]
rasterize_input_dat <- as.logical(commandArgs()[10])

tridaly_in_dir <- '<<< FILEPATH REDACTED >>>'
tridaly_out_dir <- '<<< FILEPATH REDACTED >>>'
dir.create(tridaly_out_dir)

#functions set up
# drive locations
commondir      <- '<<< FILEPATH REDACTED >>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/setup.R')
mbg_setup(package_list = package_list, repos = '<<< FILEPATH REDACTED >>>')

library(raster)
library(fasterize)
library(sf)
library(data.table)
library(rgeos)

#(2) Read in shapefile and results ####################################
if (rasterize_input_dat){
  #Shapefile
  old_shapefile <- readRDS('<<< FILEPATH REDACTED >>>')
  
  new_simple_polygon <- load_simple_polygon(
    gaul_list = NULL,
    buffer = 0.4, custom_shapefile = old_shapefile)
  
  new_subset_shape <- new_simple_polygon[["subset_shape"]]
  new_simple_polygon <- new_simple_polygon[["spoly_spdf"]]
  
  message("Loading simple raster\n")
  new_raster_list <- build_simple_raster_pop(new_subset_shape, field = 'ADM2_CODE', modeling_shapefile_version <- old_shapefile_version, raking = FALSE, pop_release = pop_release) #custom function that makes a full pops raster
  new_simple_raster <- new_raster_list[["simple_raster"]]
  pop_raster <- new_raster_list[['pop_raster']]
  
  #results
  mal_input <- fread(paste0(tridaly_in_dir, indicator, '_admin_2_raked_incidence_summary.csv')) %>%
    select(Year, ADM2_CODE, incidence_rate_infants_rmean) %>%
    rename(value = incidence_rate_infants_rmean, year = Year)
  
  #(3) Make a raster of mean results ####################################
  ad2_codes <- unique(mal_input$ADM2_CODE)
  mal_input_raster <- stack(replicate(18, new_simple_raster))
  for (n in 1:18){
    year_raster <- mal_input_raster[[n]]
    yr <- 1999 + n
    for (ad2 in ad2_codes){
      year_raster[year_raster == ad2] <- mal_input[ADM2_CODE == ad2 & year == yr, value]
    }
    
    #NA out any adm2 codes not converted to results value\
    year_raster[year_raster > 1000] <- NA
    mal_input_raster[[n]] <- year_raster
  }
  
  writeRaster(mal_input_raster, paste0(tridaly_out_dir, indicator, '_rasterized_admin2_aggregates.grd'))
  message(paste('DONE WITH', indicator, 'rasterization'))
}

#(4) Run aggregation ###########################################
message(paste('Beginning aggregation for', indicator))

for (admin in c(0:2)){
  #set identifiers by admin_level
  if (admin == 0) {
    admin_info <- c('ADM0_NAME', 'ADM0_CODE')
    admin_code <- 'ADM0_CODE'
  }
  
  if (admin == 1) {
    admin_info <- c('ADM0_NAME', 'ADM0_CODE','ADM1_NAME', 'ADM1_CODE')
    admin_code <- 'ADM1_CODE'
  }
  
  if (admin == 2) {
    admin_info <- c('ADM0_NAME', 'ADM0_CODE','ADM1_NAME', 'ADM1_CODE','ADM2_NAME', 'ADM2_CODE')
    admin_code <- 'ADM2_CODE'
  }
  
  # read in shapefile
  ashp <- shapefile('<<< FILEPATH REDACTED >>>')
  asf <- st_as_sf(ashp)
  
  # read in rasters to be aggregated
  mal_input_raster <- brick(paste0(tridaly_out_dir, indicator, '_rasterized_admin2_aggregates.grd'))
  
  # pop is the weights used to aggregate
  # Standardize rasters to the same extent and same number of pixels
  std_ras <- function(input, standard) {
    a <- crop(input, standard)
    b <- extend(a, standard)
    c <- setExtent(b, standard)
    d <- mask(c, standard)
    return(d)
  }
  
  #load worldpop raster
  pop_dir <- '<<< FILEPATH REDACTED >>>'
  setwd(pop_dir)
  pop_files <- paste0('worldpop_a0004t_1y_', c(2000:2017), '_00_00.tif')
  pop <- lapply(pop_files, raster)
  pop <- stack(pop)
  pop2 <- std_ras(pop, mal_input_raster[[1]])
  
  # make sure ADM2 codes are numeric and make a admin raster
  asf[[admin_code]] <- as.numeric(asf[[admin_code]])
  admin_raster <- fasterize(asf, mal_input_raster[[1]], field = admin_code)
  admin_raster <- mask(admin_raster, mal_input_raster[[1]])
  
  # mask everything by population
  admin_raster <- mask(admin_raster, pop2[[1]])
  mal_input_raster <- mask(mal_input_raster, pop2[[1]])
  
  # miss_admins is the admin2s that are not in the admin raster so we take their 
  # centroids since they are too small to be rasterized
  miss_admins <- setdiff(unique(asf[[admin_code]]), unique(admin_raster))
  centers <- gCentroid(ashp[which(ashp[[admin_code]] %in% miss_admins),], byid = TRUE)
  
  # pixel_spatial is a data.table where pixel_id is the cell numbers where its not-NA 
  # in the admin raster
  a <- as.numeric(t(as.matrix(admin_raster)))
  pixel_id <- which(!is.na(a))
  pixel_spatial<-data.table(pixel_id=pixel_id)
  
  # Looping over layers we are going to aggregate
  results <- list()
  for (j in 1:nlayers(mal_input_raster)) {
    message(paste('layer',j))
    
    # Extract values from the raster of interest by supplying pixel IDs
    values <- data.table(mean = raster::extract(mal_input_raster[[j]], pixel_id))
    
    # add pixel IDs to the data.table
    values[,pixel_id:= pixel_id]
    
    # Extract admin identifier by pixel_id
    values[,admin_x:=raster::extract(admin_raster, pixel_id)]
    
    # Extract population by pixel id
    values[,pop:=raster::extract(pop2[[j]], pixel_id)]
    
    # aggregate values by admin id
    agg <- values[,.(mean=weighted.mean(x = mean, w = pop, na.rm = TRUE),
                     pop = sum(pop, na.rm = TRUE)), by = admin_x]
    
    agg <- plyr::rename(agg, replace = c('admin_x' = admin_code))
    
    # Use admin IDs to get the other spatial identifiers from the shapefile
    ashp_dt <- as.data.table(ashp@data)
    ashp_dt[[admin_code]] <- as.numeric(ashp_dt[[admin_code]])
    ashp_merge <- ashp_dt[,c(admin_info, 'gadm_geoid','iso3'), with = FALSE]
    agg <- merge(agg, ashp_merge, by = admin_code,
                 all.x = TRUE)
    
    # Trim the data set to only variables of interest amd add year column
    agg[,year := j + 1999]
    
    # Extract missing admins
    if (nrow(ashp[which(ashp$ADM2_CODE %in% miss_admins),]@data[,c(admin_info, 'gadm_geoid', 'iso3')]) != 0){
      miss_admins_val <- data.table(mean = raster::extract(mal_input_raster[[j]], centers),
                                    pop = raster::extract(pop2[[j]], centers),
                                    year = 1999 +j,
                                    ashp[which(ashp$ADM2_CODE %in% miss_admins),]@data[,c(admin_info, 'gadm_geoid', 'iso3')])
      
      # Add values for the missing admins to the end of the aggregation data.table
      agg <- rbind(agg, miss_admins_val)
    }
    results[[j]] <- agg
  }
  
  # Rbind list to create final data.table
  results <- do.call(rbind, results)
  write.csv(results, '<<< FILEPATH REDACTED >>>')
  
}