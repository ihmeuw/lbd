## #############################################################################
##
## STAGE 2 PAPER: Supplemental Figures - Admin2 aroc and 2030 projection maps 
##
## Date: January 18, 2019
## Purpose: Make maps of aroc and projected estimates of 2030 
##          for stage 2 admin2s. Rakes to SDG paper projections
## #############################################################################

rm(list=ls())
## IMPORTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)
library(ggplot2)
library(foreign)
library(gridExtra)
library(grid)
library(rgdal)
library(rgeos)
library(raster)
library(matrixStats)
library(parallel)
options(scipen=9999)

source(paste0("<<<< FILEPATH REDACTED >>>>"))
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")

plot <- F

run_date <- "<<<< FILEPATH REDACTED >>>>"
age <- "infant"

modeling_shapefile_version <- "<<<< FILEPATH REDACTED >>>>"
cell_pred_year_list <- c(2000:2017)
pop_year_list <- c(2020)
raking_year_list <- c(2018:2030)

out_dir <- "<<<< FILEPATH REDACTED >>>>"

regions <- c("ansa+trsa-bra", "caca-mex", "cssa", "essa-yem", "mide+yem", "noaf", "ocea+seas-mys", "soas", "sssa", "stan+mng", "wssa")

## READ IN DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#load aroc estimates
e1 <- new.env()
load(sprintf("<<<< FILEPATH REDACTED >>>>"), e1)
mapping_data_full <- get('mapping_data_full', e1)
rm(e1)

if(age == "under5") {
  aroc_raster <- mapping_data_full$rast$aroc$under5_mean_raked_aroc
} else if(age == "infant"){
  aroc_raster <- mapping_data_full$rast$aroc$infant_mean_raked_aroc
} else {
  aroc_raster <- mapping_data_full$rast$aroc$neonatal_mean_raked_aroc
}

#load in 2030 raking factors
rf <- fread("<<<< FILEPATH REDACTED >>>>")
rf <- rf[indicator_short == "Neonatal Mort",]
rf <- rf[year_id %in% raking_year_list, c("location_id", "year_id", "mean_val")]
colnames(rf) <- c("name", "year", "mean")

## REGIONS PREP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

prep_aroc_regions <- function(reg){
  message("loading cell pred for region:", reg," \n")
  cell_pred <- readRDS("<<<< FILEPATH REDACTED >>>>")
  
  #grab the last year of the cell pred
  final_year_cell_pred <- cell_pred[((nrow(cell_pred) - (nrow(cell_pred) / length(cell_pred_year_list))) + 1):nrow(cell_pred),]
  pred_length <- nrow(final_year_cell_pred)
  final_year_cell_pred <- do.call("rbind", replicate(length(raking_year_list), final_year_cell_pred, simplify = FALSE))
  
  #load in GADM simple raster
  adm0_list <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, tolerance = 0.4,
                                        shapefile_version = modeling_shapefile_version)
  subset_shape   <- simple_polygon[['subset_shape']]
  simple_polygon <- simple_polygon[['spoly_spdf']]
  
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  pixel_spatial<-data.table(pixel_id=pixel_id)
  
  pop_raster_annual <- load_and_crop_covariates_annual(covs           = 'worldpop',                
                                                       measures       = 'a0004t',     
                                                       simple_polygon = simple_polygon,
                                                       start_year     = min(pop_year_list),
                                                       end_year       = max(pop_year_list),
                                                       interval_mo    = 12,
                                                       agebin=1)[[1]]
  
  pop_raster_annual <- crop(pop_raster_annual, extent(simple_raster))
  pop_raster_annual <- extend(pop_raster_annual, extent(simple_raster), value = NA)
  pop_raster_annual <- setExtent(pop_raster_annual, simple_raster)
  pop_raster_annual <- mask(pop_raster_annual, simple_raster)
  
  # getting the pixels from the population raster that correspond to the valid pixels in the pop-raster
  pop <- data.table(raster::extract(pop_raster_annual, pixel_id)) 
  pop[,pixel_id:=pixel_id]
  # Melting the dataframe such that it is long () and should match the number of rows in cell_pred
  pop<-melt(pop,id.vars="pixel_id")
  # Converting "worldpop.1" variables to actual years.
  pop[, year := (min(pop_year_list)) + as.numeric(gsub("worldpop.", "", variable)) - 1]
  pop<-pop[,list(pixel_id,year,pop=value)]
  # Setting values where pop is NA or 0 to 0.01 to avoid NAs and NaNs in aggregation
  pop[is.na(pop),pop:=0.01]
  pop[pop == 0 ,pop:=0.01]
  
  
  #Loading and Assigning Admin Units to Pixels
  message("Getting the spatial location (admin unit) of each of the pixel locations.")
  
  message("Rasterizing shapefiles; this may take a while.")
  admin_levels<-list() # Emtpy list of levels that will be filled with admin levels
  for(lvl in c(0,2)){
    fieldname<-paste0("ADM",lvl,"_CODE")
    admin_info<-GetAdmin(admin_level=lvl,simple_raster,adm0_list)
    pixel_spatial[[fieldname]]<-raster::extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
    admin_levels[[as.character(lvl)]]<-admin_info # Add the admin info to the list
    if(sum(is.na(pixel_spatial[[fieldname]]))>0){ # Check to see if any of the pixels don't have a location assigned
      message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level ",lvl))
    }
  }
  
  pop<-merge(pop,pixel_spatial,by="pixel_id",all.x=T) # Merging on the spatial information to the population information.
  pop<-pop[order(year,pixel_id)]
  pop$id <- 1:nrow(pop)
  
  #subset aroc raster to region
  reg_aroc_raster <- crop(aroc_raster, extent(simple_raster))
  reg_aroc_raster <- setExtent(reg_aroc_raster, simple_raster)
  reg_aroc_raster <- mask(reg_aroc_raster, simple_raster)
  
  #get pixels in simple raster that are NA
  simple_raster_extract <- raster::extract(simple_raster, extent(simple_raster))
  simple_raster_subsetter <- !is.na(simple_raster_extract)
  
  #pull pixel values into list and subset based on non-NA pixels in simple raster
  #the number of rows after dropping NAs do not match up so this is necessary
  reg_aroc_extract <- raster::extract(reg_aroc_raster, extent(reg_aroc_raster))
  reg_aroc_extract <- reg_aroc_extract[simple_raster_subsetter]
  
  if(length(reg_aroc_extract) != nrow(pop)) {
    message("aroc pixel estimates do not match up to each pixel")
  }
  
  reg_aroc_extract <- rep(reg_aroc_extract, times = length(raking_year_list))
  pop[, year := NULL]
  year_bind <- rep(raking_year_list, each = pred_length)
  
  #merge on aroc pixel data to cell pred
  aroc_cell_pred <- data.table(cbind(final_year_cell_pred, reg_aroc_extract, pop, year_bind))
  aroc_cell_pred[, n := year_bind - max(cell_pred_year_list)]
  
  aroc_cell_pred[, aroc := exp(reg_aroc_extract * n)]
  
  overs <- paste0("V", 1:ncol(final_year_cell_pred))
  cell_pred_2030 <- aroc_cell_pred[, (overs) := lapply(overs, function(x){get(x) * aroc})]
  
  ## CUSTOM RAKE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #infant has no raking targets, returns unraked, no rasters
  if(age == "infant"){
    message("starting raking process for region:", reg," \n")
    admin2_cell_pred <- cell_pred_2030[, lapply(overs, function(x) weighted.mean(get(x), pop,  na.rm = T)), by = c('ADM2_CODE', 'ADM0_CODE', 'year_bind')]
    admin2_matrix <- as.matrix(admin2_cell_pred[,overs, with = F])
    admin2_cell_pred[, mean := rowMeans(.SD, na.rm =T), .SDcols = overs]
    admin2_cell_pred[, upper := matrixStats::rowQuantiles(admin2_matrix,  probs = 97.5 / 100)]
    admin2_cell_pred[, lower := matrixStats::rowQuantiles(admin2_matrix,  probs = 2.5 / 100)]
    admin2_summary <- admin2_cell_pred[, c("ADM0_CODE", "ADM2_CODE", "year_bind", "mean", "lower", "upper")]
    setnames(admin2_summary, "year_bind", "year")
    
    aroc_summary <- cell_pred_2030[, c("ADM0_CODE", "ADM2_CODE", "reg_aroc_extract", "pop")]
    aroc_summary <- aroc_summary[, aroc := weighted.mean(reg_aroc_extract, pop,  na.rm = T), by = c('ADM2_CODE', 'ADM0_CODE')]
    
    admin2_name <- paste0("admin2_summary_", reg)
    aroc_name <- paste0("aroc_summary_", reg)
    # mean_raster_name <- paste0("mean_raster_", reg)
    # upper_raster_name <- paste0("upper_raster_", reg)
    # lower_raster_name <- paste0("lower_raster_", reg)
    return_list <- list(admin2_summary, aroc_summary) #, mean_raster, upper_raster, lower_raster)
    names(return_list) <- c(admin2_name, aroc_name)#, mean_raster_name, upper_raster_name, lower_raster_name)
    return(return_list)
  }
  else {
    message("starting raking process for region:", reg," \n")
    cell_pred_raker <- cell_pred_2030[, lapply(overs, function(x) weighted.mean(get(x), pop,  na.rm = T)), by = c('ADM0_CODE', 'year_bind')]
    cell_pred_raker[, mean := rowMeans(.SD, na.rm =T), .SDcols = overs]
    cell_pred_raker <- cell_pred_raker[, c('ADM0_CODE', 'mean', 'year_bind')]
    
    stage_list <- fread("<<<< FILEPATH REDACTED >>>>")
    locations <- stage_list[, c("gadm_geoid", "loc_id")]
    
    rf <- merge(rf, locations, by.x = "name", by.y = "loc_id", all.x =T)
    
    raker <- merge(cell_pred_raker, rf, by.x = c("ADM0_CODE", "year_bind"), by.y = c("gadm_geoid", "year"), all.x = T)
    raker[, rf := mean.y / mean.x]
    raker <- raker[, c("ADM0_CODE", "rf", "year_bind")]
    
    
    meta_raked_cell_pred <- merge(cell_pred_2030, raker, by = c("ADM0_CODE", "year_bind"))
    meta_raked_cell_pred <- setorder(meta_raked_cell_pred, year_bind, id)
    
    #make mean/upper/lower raster bricks of raked projected values
    raked_cell_pred <- meta_raked_cell_pred[, lapply(overs, function(x) get(x) * rf)]
    raked_cell_pred_matrix <- as.matrix(raked_cell_pred[, overs, with = F])
    mean_raster <- cell_pred_to_raster_brick(raked_cell_pred_matrix, simple_raster, raking_year_list, "mean")
    upper_raster <- cell_pred_to_raster_brick(raked_cell_pred_matrix, simple_raster, raking_year_list, "upper")
    lower_raster <- cell_pred_to_raster_brick(raked_cell_pred_matrix, simple_raster, raking_year_list, "lower")
    names(mean_raster) <- as.character(raking_year_list)
    names(upper_raster) <- as.character(raking_year_list)
    names(lower_raster) <- as.character(raking_year_list)
    
    #make adm2 summary of raked projected values
    admin2_cell_pred <- meta_raked_cell_pred[, lapply(overs, function(x) get(x) * rf), by = c('ADM2_CODE', 'ADM0_CODE', 'pop', 'year_bind')]
    admin2_cell_pred <- admin2_cell_pred[, lapply(overs, function(x) weighted.mean(get(x), pop,  na.rm = T)), by = c('ADM2_CODE', 'ADM0_CODE', 'year_bind')]
    admin2_matrix <- as.matrix(admin2_cell_pred[,overs, with = F])
    admin2_cell_pred[, mean := rowMeans(.SD, na.rm =T), .SDcols = overs]
    admin2_cell_pred[, upper := matrixStats::rowQuantiles(admin2_matrix,  probs = 97.5 / 100)]
    admin2_cell_pred[, lower := matrixStats::rowQuantiles(admin2_matrix,  probs = 2.5 / 100)]
    admin2_summary <- admin2_cell_pred[, c("ADM0_CODE", "ADM2_CODE", "year_bind", "mean", "lower", "upper")]
    setnames(admin2_summary, "year_bind", "year")
    
    aroc_summary <- cell_pred_2030[, c("ADM0_CODE", "ADM2_CODE", "reg_aroc_extract", "pop")]
    aroc_summary <- aroc_summary[, aroc := weighted.mean(reg_aroc_extract, pop,  na.rm = T), by = c('ADM2_CODE', 'ADM0_CODE')]
    
    admin2_name <- paste0("admin2_summary_", reg)
    aroc_name <- paste0("aroc_summary_", reg)
    mean_raster_name <- paste0("mean_raster_", reg)
    upper_raster_name <- paste0("upper_raster_", reg)
    lower_raster_name <- paste0("lower_raster_", reg)
    return_list <- list(admin2_summary, aroc_summary, mean_raster, upper_raster, lower_raster)
    names(return_list) <- c(admin2_name, aroc_name, mean_raster_name, upper_raster_name, lower_raster_name)
    return(return_list)
  }
}

#prep data
prepped_data <- mclapply(X=regions, FUN=prep_aroc_regions, mc.cores = 13)
prepped_data <- unlist(prepped_data, recursive = F)

#split out 2030 estimates and aroc summary and rbind together
prepped_2030x <- prepped_data[grep("admin2", names(prepped_data))]
prepped_arocx <- prepped_data[grep("aroc_summary", names(prepped_data))]

prepped_2030 <- do.call('rbind', prepped_2030x)
prepped_aroc <- do.call('rbind', prepped_arocx)
prepped_aroc <- unique(prepped_aroc[, c("ADM2_CODE", "ADM0_CODE", "aroc")])

dir.create(sprintf("<<<< FILEPATH REDACTED >>>>"), showWarnings = FALSE, recursive = T)

if(age != "infant") {
  mean_rasterx <- prepped_data[grep("mean_raster", names(prepped_data))]
  upper_rasterx <- prepped_data[grep("upper_raster", names(prepped_data))]
  lower_rasterx <- prepped_data[grep("lower_raster", names(prepped_data))]
  
  mean_raster <- Reduce(function(...)merge(...),mean_rasterx)
  upper_raster <- Reduce(function(...)merge(...),upper_rasterx)
  lower_raster <- Reduce(function(...)merge(...),lower_rasterx)
  
  raster_names <- names(mean_rasterx[[1]])
  names(mean_raster) <- raster_names
  names(upper_raster) <- raster_names
  names(lower_raster) <- raster_names
  
  writeRaster(mean_raster, file= "<<<< FILEPATH REDACTED >>>>", format = "GTiff", overwrite=T)
  writeRaster(upper_raster, file= "<<<< FILEPATH REDACTED >>>>", format = "GTiff", overwrite=T)
  writeRaster(lower_raster, file= "<<<< FILEPATH REDACTED >>>>", format = "GTiff", overwrite=T)
}

saveRDS(prepped_2030, "<<<< FILEPATH REDACTED >>>>")
write.csv(prepped_2030, "<<<< FILEPATH REDACTED >>>>")

## Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if(plot){
  S2_BOUNDS <- list(
    'lat_min'  = -57,
    'lat_max'  = 63,
    'long_min' = -129,
    'long_max' = 165
  )
  
  pull_format_admin_shp <- function(
    adm_level,
    shp_version = 'current',
    simp_prop   = NULL,
    stage2_only = FALSE,
    include_non_gbd = TRUE
  ){
    # Load the shapefile, simplifying Admin 0 shapefiles a bit by default
    if(is.null(simp_prop) & (adm_level==0)) simp_prop <- 0.05
    shp <- fast_load_shapefile(
      shp_path = get_admin_shapefile(admin_level=adm_level, version=shp_version),
      simplify_tol=simp_prop
    )
    # Pull the admin0 lookup table 
    lookup_table <- load_adm0_lookup_table()
    adm_field <- ifelse(
      detect_adm_shapefile_date_type(
        shpfile_path=get_admin_shapefile(admin_level=0, version=shp_version)
      )$shpfile_type=='gaul',
      'GAUL_CODE',
      'gadm_geoid'
    )
    setnames(lookup_table, adm_field, 'ADM0_CODE')
    if(include_non_gbd==FALSE){
      # Drop all non-GBD locations, which will not have a valid location ID
      lookup_table <- lookup_table[ loc_id > 0, ]
    }
    if(stage2_only){
      lookup_table <- lookup_table[ Stage %in% c('1','2a','2b'), ]
    }
    # Fortify for ggplotting, but don't add fields
    shp$admin_code <- shp@data[, paste0('ADM',adm_level,'_CODE')]
    shp@data       <- shp@data[, c('ADM0_CODE','admin_code')]
    shp_fort <- prep_shp_data_for_mapping(
      shp       = shp,
      dataset   = lookup_table[, .(ADM0_CODE, location_name)],
      merge_var ='ADM0_CODE'
    )
    # Subset to the Stage 2 boundaries
    shp_fort[lat  < S2_BOUNDS$lat_min,  lat  := S2_BOUNDS$lat_min  ]
    shp_fort[lat  > S2_BOUNDS$lat_max,  lat  := S2_BOUNDS$lat_max  ]
    shp_fort[long < S2_BOUNDS$long_min, long := S2_BOUNDS$long_min ]
    shp_fort[long > S2_BOUNDS$long_max, long := S2_BOUNDS$long_max ]
    message(sprintf("Done prepping admin%s data.\n",adm_level))
    return(shp_fort)
  }
  
  region_adm0_list <- get_adm0_codes("stage1+stage2", shapefile_version = modeling_shapefile_version)
  simple_polygon <- load_simple_polygon(gaul_list = region_adm0_list, buffer = 1, tolerance = 0.4,
                                        shapefile_version = modeling_shapefile_version)
  subset_shape   <- simple_polygon[['subset_shape']]
  simple_polygon <- simple_polygon[['spoly_spdf']]
  
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  admin_info<-GetAdmin(admin_level=2,simple_raster,region_adm0_list)
  bg_shp <- as(admin_info$spdf, 'SpatialPolygonsDataFrame')
  bg_for_mapping <- fortify(bg_shp)
  
  bg_for_mapping_adm0 <- pull_format_admin_shp(
    adm_level   = 0,
    shp_version = modeling_shapefile_version,
    simp_prop   = NULL,
    stage2_only = FALSE,
    include_non_gbd = TRUE
  )
  
  #color settings for 2030 map
  na_color <- '#E6E6E6'
  col_breaks_abs <- c(0, 5, 25, 50, 100, 200)
  col_labs_abs <- c("", "5", "25", "50", "100", "200")
  grad_vals_abs <- c(
    '#009999','#5aada1','#88c1a9','#b1d6b0','#d8eab8','#ffffbf',"#fff2b6",
    "#fee5ae","#fed8a5","#fecb9c","#fdbe93","#fdb18b","#fda482","#fc9779",
    "#fc8a70","#fc7d68","#fb705f","#fb6356","#fb564d","#fa4945","#fa3c3c",
    "#f53b3a","#f03938","#eb3837","#e63635","#e13533","#db3431","#d6322f",
    "#d1312e","#cc2f2c","#c72e2a","#c22d28","#bd2b26","#b82a25","#b32823",
    "#ae2721","#a8261f","#a3241d","#9e231c","#99211a","#942018"
  )
  
  aroc_col_breaks <- c(-0.08, -0.05, -0.02, 0.0)
  aroc_col_labs <- c('8%','5%', '2%', 'No decline')
  aroc_grad_vals <- c(
    '#002080', '#003866', '#00504C', '#006833', '#008019', '#009900', 
    '#53B639', '#A6D373', '#FAF0AD', '#EDA259', '#E15505'
  )
  
  
  mapping_data <- suppressMessages(prep_shp_data_for_mapping(shp=bg_shp, dataset=prepped_aroc, merge_var='ADM2_CODE'))
  if(age == "under5") {
    legend_title <- 'Annual Percent Rate of\nChange in 5q0\n(per 1000 live births)'
  } else if(age == "infant"){
    legend_title <- 'Annual Percent Rate of\nChange in 1q0\n(per 1000 live births)'
  } else {
    legend_title <- 'Annual Percent Rate of\nChange in neonatal\nmortality\n(per 1000 live births)'
  }
  
  
  aroc_plot <- ggplot() +
    geom_polygon(data = mapping_data, aes(x=long, y=lat, group=group, fill=aroc)) + 
    geom_path(data = bg_for_mapping_adm0,
              aes(x=long, y=lat, group=group),
              color='#000000',
              lwd=0.1) +
    scale_fill_gradientn(limits = c(min(prepped_aroc$aroc, na.rm=T), max(prepped_aroc$aroc, na.rm=T)),
                         colors = aroc_grad_vals,
                         breaks = aroc_col_breaks,
                         labels = aroc_col_labs,
                         na.value = "grey50") +
    labs(#title ='Annual Rate of Change',
      fill = legend_title) +
    theme_map() +
    theme(legend.position = c(0.95, 0.80),
          legend.title = element_text(size=9, color='#222222'),
          legend.text  = element_text(size=6, color='#222222')) +
    coord_map(
      projection = "mollweide",
      xlim = c(-108, 157.5), 
      ylim = c(-35.5,  53)
    )
  
  # Plot
  png(sprintf("<<<< FILEPATH REDACTED >>>>"), 
      height=5.5, width=12, units='in', res=800)
  print(aroc_plot)
  dev.off()
  
  
  if(age == "under5"){
    prepped_2030_data <- prepped_2030[year == 2030,]
    
    mapping_data <- suppressMessages(prep_shp_data_for_mapping(shp=bg_shp, dataset=prepped_2030, merge_var='ADM2_CODE'))
    mapping_data$mean <- mapping_data$mean * 1000
    
    mapping_data <- mapping_data[!(ADM0_NAME %in% c("French Guiana", "Western Sahara")),]
    
    admin2030 <- ggplot() +
      geom_polygon(data = mapping_data,
                   aes(x=long, y=lat, group=group, fill=mean), show.legend = T) + 
      geom_path(data = bg_for_mapping_adm0,
                aes(x=long, y=lat, group=group),
                color='#000000',
                lwd=0.1) +
      scale_fill_gradientn(limits = c(0, 200),
                           colors = grad_vals_abs,
                           breaks = col_breaks_abs,
                           labels = col_labs_abs,
                           na.value = na_color,
                           name = '5q0 per 1000\nLive Births') +
      #labs(title = 'Local Burden of Disease U5M Projection: 2030') +
      theme_map() +
      theme(legend.position = c(0.95, 0.80),
            legend.title = element_text(size=9, color='#222222'),
            legend.text  = element_text(size=6, color='#222222')) +
      coord_map(
        projection = "mollweide",
        xlim = c(-108, 157.5), 
        ylim = c(-35.5,  53)
      )
    
    # Plot
    png(
      sprintf("<<<< FILEPATH REDACTED >>>>"), 
      height=5.5, width=12, units='in', res=800
    )
    print(admin2030)
    dev.off()
  }
}