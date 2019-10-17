## #############################################################################
##
## STAGE 2 PAPER: Supplemental Figures - compare to stage 1 paper estimates 
##
## Date: January 18, 2019
## Purpose: Make scatterplots of stage 1 (previous) estimates and stage 2 
##          (current) estimates for neonatal and under5 age bins, 2005, 2010
##          and 2015. Also makes maps of absolute difference, relative
##          difference, and CI not overlapping. Generates a table of new data
##          sources used in the current estimates, and prints some numbers
##          for numberplugging the SI of the stage 2 paper.
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
options(scipen=9999)

source(paste0("<<<< FILEPATH REDACTED >>>>"))
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")

run_date <- "<<<< FILEPATH REDACTED >>>>"

modeling_shapefile_version <- "<<<< FILEPATH REDACTED >>>>"
#under5 or neonatal
age <- "under5"

out_dir <- "<<<< FILEPATH REDACTED >>>>"

#for the data comparison, the file path to the most up to date inclusion table (must be csv)
inclusion_filepath <- "<<<< FILEPATH REDACTED >>>>"

## READ IN DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#stage 2 africa adm2 aggregates
stage2 <- fread("<<<< FILEPATH REDACTED >>>>")
stage2_upper <- fread("<<<< FILEPATH REDACTED >>>>")
stage2_lower <- fread("<<<< FILEPATH REDACTED >>>>")

stage2 <- merge(stage2, stage2_upper, by = c("ADM2_CODE", "year"))
stage2 <- merge(stage2, stage2_lower, by = c("ADM2_CODE", "year"))

#stage 1 cell pred, saved as "x"
load("<<<< FILEPATH REDACTED >>>>")

#raster to serve in place of simple raster
simple_fake <- raster("<<<< FILEPATH REDACTED >>>>")
simple_list <- raster::extract(simple_fake, extent(simple_fake))
simple_list <- na.omit(simple_list)

#drop NAs to match with simple_fake
cell_pred <- na.omit(x)

#check that cell pred matches simple_fake
if((nrow(cell_pred) /4) != length(simple_list)){
  message("cell pred and raster substitute for simple raster do not align")
}

#load in GADM simple raster
adm0_list <- get_adm0_codes("africa", shapefile_version = modeling_shapefile_version)
simple_polygon <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, tolerance = 0.4,
                                      shapefile_version = modeling_shapefile_version)
subset_shape   <- simple_polygon[['subset_shape']]
simple_polygon <- simple_polygon[['spoly_spdf']]

raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]

simple_raster <- extend(simple_raster, extent(simple_fake), value = NA)

## CROSSWALK ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#convert cell pred to match GADM simple raster, using simple_fake as reference
cross_cell <- crosswalk_cell_pred_add_NA(simple_fake, simple_raster, cell_pred, year_list = 1:4)

#check that crosswalk was successful
new_simple_list <- raster::extract(simple_raster, extent(simple_raster))
new_simple_list <- na.omit(new_simple_list)
if((nrow(cross_cell) / 4) != length(new_simple_list)){
  message("cell pred and new simple raster do not align, problem with crosswalk")
}

## AGGREGATE STAGE 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pixel_id <- seegSDM:::notMissingIdx(simple_raster)
pixel_spatial<-data.table(pixel_id=pixel_id)

year_list <- c(2000,2005,2010,2015)
pop_raster_annual <- load_and_crop_covariates_annual(covs           = 'worldpop',                
                                                     measures       = 'a0004t',     
                                                     simple_polygon = simple_polygon,
                                                     start_year     = 2000,
                                                     end_year       = 2015,
                                                     interval_mo    = 60,
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
pop[, year := (min(year_list)) + as.numeric(gsub("worldpop.", "", variable)) * 5 - 5]
pop<-pop[,list(pixel_id,year,pop=value)]
# Setting values where pop is NA or 0 to 0.01 to avoid NAs and NaNs in aggregation
pop[is.na(pop),pop:=0.01]
pop[pop == 0 ,pop:=0.01]


#Loading and Assigning Admin Units to Pixels
message("Getting the spatial location (admin unit) of each of the pixel locations.")

region_adm0_list<-get_adm0_codes("africa", shapefile_version = modeling_shapefile_version) # Getting the adm0 GAUL codes, we can use this to make sure we don't accidentally include countries from buffers that shouldn't be in this region

message("Rasterizing shapefiles; this may take a while.")
admin_levels<-list() # Emtpy list of levels that will be filled with admin levels
for(lvl in c(0,2)){
  fieldname<-paste0("ADM",lvl,"_CODE")
  admin_info<-GetAdmin(admin_level=lvl,simple_raster,region_adm0_list)
  pixel_spatial[[fieldname]]<-raster::extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
  admin_levels[[as.character(lvl)]]<-admin_info # Add the admin info to the list
  if(sum(is.na(pixel_spatial[[fieldname]]))>0){ # Check to see if any of the pixels don't have a location assigned
    message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level ",lvl))
  }
}

pop<-merge(pop,pixel_spatial,by="pixel_id",all.x=T) # Merging on the spatial information to the population information.
pop<-pop[order(year,pixel_id)]
pop$id <- 1:nrow(pop)

#merge on population and ADM2 codes onto stage1 cell pred
stage1 <- cbind(cross_cell, pop)

#aggregate pixels to ADM2
overs <- paste0("V", c(1:1000))
stage1_agg <- stage1[ADM0_CODE %in% region_adm0_list,lapply(.SD,weighted.mean,w=pop, na.rm=T), by=c("year", "ADM2_CODE","ADM0_CODE"), .SDcols=overs]
#aggregate draws
stage1_agg2 <- stage1_agg[, stage1_mean := rowMeans(.SD, na.rm = T), .SDcols = overs]

#convert draws to matrix to calculate row quantiles faster
stage1_matrix <- as.matrix(stage1_agg2[, 3:1002])
stage1_upper <- rowQuantiles(stage1_matrix, probs = 97.5 / 100)
stage1_lower <- rowQuantiles(stage1_matrix, probs = 2.5 / 100)

stage1_CI <- cbind(stage1_agg2, stage1_upper, stage1_lower)
#drop draws
stage1_fin <- stage1_CI[, c("year", "ADM2_CODE", "ADM0_CODE", "stage1_mean", "stage1_lower", "stage1_upper")]

## AGGREGATE STAGE 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

stage2[year <= 2017, year_bin := 2015]
stage2[year <= 2012, year_bin := 2010]
stage2[year <= 2007, year_bin := 2005]
stage2[year <= 2002, year_bin := 2000]

stage2_fin <- stage2[, .(stage2_mean = mean(mean, na.rm=T), stage2_upper = mean(upper, na.rm = T), stage2_lower = mean(lower, na.rm = T)), by = c("ADM2_CODE", "year_bin")]
setnames(stage2_fin , c("ADM2_CODE", "year", "stage2_mean", "stage2_upper", "stage2_lower"))

## PLOT PREP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

com <- merge(stage1_fin, stage2_fin, by = c("ADM2_CODE", "year"), all.x = T)
com <- com[year != 2000,]
africa_excluded_adm0_list <- get_adm0_codes("dza+tun+lby+cpv+com+stp+mus+esh", shapefile_version = modeling_shapefile_version)
com[ADM0_CODE %in% africa_excluded_adm0_list, c("stage1_mean", "stage1_lower", "stage1_upper", "stage2_mean", "stage2_lower", "stage2_upper")] <- NA

com[ADM0_CODE %in% africa_excluded_adm0_list, ]

com$stage2_mean <- com$stage2_mean * 1000
com$stage1_mean <- com$stage1_mean * 1000
com$stage2_upper <- com$stage2_upper * 1000
com$stage1_upper <- com$stage1_upper * 1000
com$stage2_lower <- com$stage2_lower * 1000
com$stage1_lower <- com$stage1_lower * 1000
com[, abs_diff:= stage2_mean - stage1_mean]
com[, rel_diff:= stage2_mean / stage1_mean]

# Check for non-overlapping CIs
com[, c('stage2_too_high','stage2_too_low','no_overlap') := 0]
com[
  !is.na(stage1_upper) & !is.na(stage2_lower) & (stage1_upper < stage2_lower), 
  stage2_too_high := 1
  ]
com[
  !is.na(stage1_lower) & !is.na(stage2_upper) & (stage1_lower > stage2_upper),
  stage2_too_low := 1
  ]
# Too high --> 1, Too low --> -1
com[, no_overlap := 0 + stage2_too_high - stage2_too_low]
com[ is.na(stage1_mean) | is.na(stage2_mean), no_overlap := NA]

bg_shp <- as(admin_info$spdf, 'SpatialPolygonsDataFrame')
bg_for_mapping <- fortify(bg_shp)

adm0_shp<-GetAdmin(admin_level=0,simple_raster,region_adm0_list)
bg_for_mapping_adm0 <- fortify(adm0_shp$spdf)

## Set color schemes here
col_breaks_abs <- c(0, 5, 15, 50, 100)
col_labs_abs <- c("", "5", "15", "50", "100")
grad_vals_abs <- c(
  "#8436a8","#8436a8","#b14dac","#de64af","#efb2b7","#ffffbf","#fff2b6",
  "#fee5ae","#fed8a5","#fecb9c","#fdbe93","#fdb18b","#fda482","#fc9779",
  "#fc8a70","#fc7d68","#fb705f","#fb6356","#fb564d","#fa4945","#fa3c3c",
  "#f53b3a","#f03938","#eb3837","#e63635","#e13533","#db3431","#d6322f",
  "#d1312e","#cc2f2c","#c72e2a","#c22d28","#bd2b26","#b82a25","#b32823",
  "#ae2721","#a8261f","#a3241d","#9e231c","#99211a","#942018")

col_breaks_adiff <- c(-50, -25, 0, 25, 50, 75)
col_labs_adiff <- c("-50", "-25", "0", "+25", "+50", "+75")
grad_vals_adiff <- c('#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7',
                     '#d8daeb','#b2abd2','#8073ac','#542788')

col_breaks_reldiff <- c(.5, .66, 1, 1.5, 2, 3, 4)
col_labs_reldiff <- c(".5", ".66", "1", "1.5", "2", "3", "4")
grad_vals_reldiff <- c('#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7',
                       '#d8daeb','#b2abd2','#8073ac','#542788')

## SCATTERPLOT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

u5m_upper = ifelse(age == "under5", 300, 100)
fig_scatter <- ggplot(data=com) + 
  geom_point(aes(x=stage2_mean, y=stage1_mean)) +
  facet_wrap('year', ncol=4) +
  geom_abline(slope=1, intercept=0) + 
  guides(fill = guide_colorbar(barheight=15, nbin=100)) +
  labs(#title = 'Comparison between Stage 1 paper and new estimates, ADM2',
       x = sprintf('Current Analysis (%s per 1000 LB)', ifelse(age == "under5", "5q0", "Neonatal Mortality")),
       y = sprintf('Previous Analysis (%s per 1000 LB)',ifelse(age == "under5", "5q0", "Neonatal Mortality")),
       fill = 'Observations') +
  lims(x = c(0,u5m_upper), y = c(0,u5m_upper)) +
  theme_bw()

png(sprintf('%s/%s_stage1_2_scatter.png',out_dir, age),
    height=8, width=8, units='in', res=800)
print(fig_scatter)
dev.off()


## PLOTTING  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(this_year in c(2005,2010,2015)){
  # Check overlap specifically for this year
  this_year_data <- com[year==this_year,]

  # Create mapping dataset
  mapping_data <- suppressMessages(
    prep_shp_data_for_mapping(shp=admin_info$spdf, dataset=this_year_data, merge_var='ADM2_CODE'))
  
  message("printing absolute difference plot for ", this_year, ", ", age)
  #Absolute difference plot
  abs_min <- min(com$abs_diff, na.rm =T)
  abs_max <- max(com$abs_diff, na.rm = T)
  message("absolute mean difference: ", this_year, ", ", age, " - ", mean(this_year_data$abs_diff, na.rm=T))
  fig_adiff_this_year <- ggplot() +
    geom_polygon(data = mapping_data, aes(x=long, y=lat, group=group, fill=abs_diff)) + 
    geom_path(data = bg_for_mapping,
              aes(x=long, y=lat, group=group),
              color='#404040',
              lwd=0.1) +
    geom_path(data = bg_for_mapping_adm0,
              aes(x=long, y=lat, group=group),
              color='#000000',
              lwd=0.175) +
    scale_fill_gradient2(limits = c(abs_min, abs_max),
                         low = '#b35806',
                         high = '#542788',
                         mid = "white",
                         #colors = grad_vals_adiff,
                         breaks = col_breaks_adiff,
                         labels = col_labs_adiff,
                         midpoint = 0,
                         na.value = "grey50") +
    labs(#title = sprintf('Absolute Difference between Current and Previous Analysis: %s',this_year),
         fill = sprintf('Difference in\n%s per 1000\nLive Births\n(Current - Previous)', ifelse(age == "under5", "5q0", "Neonatal Mortality"))) +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data) +
    lims(x = c(-117.5, 154), y = c(-35.5, 99))
  
  # Plot
  png("<<<< FILEPATH REDACTED >>>>"), 
      height=5.5, width=12, units='in', res=800)
  print(fig_adiff_this_year)
  dev.off()
  
  message("printing relative difference plot for ", this_year, ", ", age)
  #Relative difference plot
  rel_min <- min(this_year_data$rel_diff)
  rel_max <- max(this_year_data$rel_diff)
  message("relative mean difference: ", this_year, ", ", age, " - ", mean(this_year_data$rel_diff, na.rm=T))
  fig_reldiff_this_year <- ggplot() +
    geom_polygon(data = mapping_data, aes(x=long, y=lat, group=group, fill=rel_diff)) + 
    geom_path(data = bg_for_mapping,
              aes(x=long, y=lat, group=group),
              color='#404040',
              lwd=0.1) +
    geom_path(data = bg_for_mapping_adm0,
              aes(x=long, y=lat, group=group),
              color='#000000',
              lwd=0.175) +
    scale_fill_gradient2(limits = c(rel_min, rel_max),
                         low = '#b35806',
                         high = '#542788',
                         mid = "white",
                         #colors = grad_vals_adiff,
                         breaks = col_breaks_reldiff,
                         labels = col_labs_reldiff,
                         midpoint = 1,
                         na.value = "grey50") +
    labs(#title = sprintf('Relative Difference between Stage 1 Paper and Stage 2: %s',this_year),
         fill = sprintf('Difference in\n%s per 1000\nLive Births\n(Current / Previous)', ifelse(age == "under5", "5q0", "Neonatal Mortality"))) +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data)
  
  # Plot
  png(sprintf("<<<< FILEPATH REDACTED >>>>"), 
      height=5.5, width=12, units='in', res=800)
  print(fig_reldiff_this_year)
  dev.off()
  
  message("printing CI overlap plot for ", this_year, ", ", age)
  mapping_data[, no_overlap := as.character(no_overlap)]
  total_units <- nrow(this_year_data)
  units_no_overlap <- sum(this_year_data$stage2_too_high, na.rm = T) + sum(this_year_data$stage2_too_low, na.rm = T)
  message("out of ", total_units, " total adm2s, ", units_no_overlap, "(",round(units_no_overlap / total_units,2), "%) did not overlap for: ", this_year, ", ", age)
  fig_uis_this_year <- ggplot() +
    geom_polygon(data = mapping_data,
                 aes(x=long, y=lat, group=group, fill=no_overlap)) + 
    geom_path(data = bg_for_mapping,
              aes(x=long, y=lat, group=group),
              color='#404040',
              lwd=0.1) +
    geom_path(data = bg_for_mapping_adm0,
              aes(x=long, y=lat, group=group),
              color='#000000',
              lwd=0.175) +
    scale_fill_manual(values = c('-1'='#0066cc','0'='#ffffe6','1'='#ff0066', 'NA'="grey50"),
                      labels = c("Current < Previous","Overlap","Current > Previous", "No previous estimates"),
                      na.value = "grey50") +
    labs(#title = sprintf('Areas with non-overlapping UIs: %s',this_year),
         fill = 'Difference\nbetween UIs') +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data) +

  
  # Plot
  png(
    sprintf("<<<< FILEPATH REDACTED >>>>"), 
    height=5.5, width=12, units='in', res=800
  )
  print(fig_uis_this_year)
  dev.off()
} 

## SI numberplug  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for(this_year in c(2005,2010,2015)){
  sink(file="<<<< FILEPATH REDACTED >>>>", append = T)
  this_year_data <- com[year==this_year,]
  print(paste0("absolute mean difference: ", this_year, ", ", age, " - ", round(mean(this_year_data$abs_diff, na.rm=T), 3)))
  print(paste0("relative mean difference: ", this_year, ", ", age, " - ", round(mean(this_year_data$rel_diff, na.rm=T),3)))
  total_units <- nrow(this_year_data)
  units_no_overlap <- sum(this_year_data$stage2_too_high, na.rm = T) + sum(this_year_data$stage2_too_low, na.rm = T)
  print(paste0("out of ", total_units, " total adm2s, ", units_no_overlap, "(",round(units_no_overlap / total_units,2), "%) did not overlap for: ", this_year, ", ", age))
  
  this_year_data <- na.omit(this_year_data)
  print(paste0("correlation for: ", this_year, " ", age, " = ", cor(this_year_data$stage2_mean, this_year_data$stage1_mean)))
  sink()
}
