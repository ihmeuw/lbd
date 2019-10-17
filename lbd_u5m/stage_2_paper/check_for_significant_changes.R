rm(list=ls())


## signficant decrease
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
library(sf)
library(sp)
options(scipen=9999)

source("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
source("<<<< FILEPATH REDACTED >>>>")

run_date <- "<<<< FILEPATH REDACTED >>>>"
modeling_shapefile_version <- "<<<< FILEPATH REDACTED >>>>"
out_dir <- "<<<< FILEPATH REDACTED >>>>"
drop_countries <- c("French Guiana", "Western Sahara", "Malaysia", "Brazil", "Mexico")

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

bg_for_mapping_adm0 <- pull_format_admin_shp(
  adm_level   = 0,
  shp_version = modeling_shapefile_version,
  simp_prop   = NULL,
  stage2_only = FALSE,
  include_non_gbd = TRUE
)

bg_shp <- pull_format_admin_shp(
  adm_level   = 2,
  shp_version = modeling_shapefile_version,
  simp_prop   = NULL,
  stage2_only = TRUE,
  include_non_gbd = TRUE
)

bg_shp <- bg_shp[!(location_name %in% drop_countries),]

for(age in c("under5", "neonatal", "infant")){
  adm2_upper <- fread("<<<< FILEPATH REDACTED >>>>")
  adm2_lower <- fread("<<<< FILEPATH REDACTED >>>>")
  
  adm1_upper <- fread("<<<< FILEPATH REDACTED >>>>")
  adm1_lower <- fread("<<<< FILEPATH REDACTED >>>>")
  
  adm2_2000 <- adm2_upper[year == 2000,]
  adm2_2017 <- adm2_lower[year == 2017,]
  
  adm1_2000 <- adm1_upper[year == 2000,]
  adm1_2017 <- adm1_lower[year == 2017,]
  
  setnames(adm2_2000, "value", "upper_2000")
  setnames(adm2_2017, "value", "lower_2017")
  setnames(adm1_2000, "value", "upper_2000")
  setnames(adm1_2017, "value", "lower_2017")
  
  adm2_comparison <- merge(adm2_2000, adm2_2017, by = "ADM2_CODE")
  adm1_comparison <- merge(adm1_2000, adm1_2017, by = "ADM1_CODE")
  
  adm2 <- adm2_comparison[lower_2017 > upper_2000,]
  adm1 <- adm1_comparison[lower_2017 > upper_2000,]
  
  #print if significant increase in mortality rate from 2000 to 2017
  if(nrow(adm2) > 0){
    print(adm2)
  }
  
  if(nrow(adm1) > 0){
    print(adm1)
  }
  
  adm2_2000 <- adm2_lower[year == 2000,]
  adm2_2017 <- adm2_upper[year == 2017,]
  
  adm1_2000 <- adm1_lower[year == 2000,]
  adm1_2017 <- adm1_upper[year == 2017,]
  
  setnames(adm2_2000, "value", "lower_2000")
  setnames(adm2_2017, "value", "upper_2017")
  setnames(adm1_2000, "value", "lower_2000")
  setnames(adm1_2017, "value", "upper_2017")
  
  adm2_comparison <- merge(adm2_2000, adm2_2017, by = "ADM2_CODE")
  adm1_comparison <- merge(adm1_2000, adm1_2017, by = "ADM1_CODE")
  
  number_adm2s <- nrow(adm2_comparison)
  number_adm1s <- nrow(adm1_comparison)
  
  adm2 <- adm2_comparison[upper_2017 < lower_2000,]
  adm1 <- adm1_comparison[upper_2017 < lower_2000,]
  
  sink("<<<< FILEPATH REDACTED >>>>", append = T)
  #print % and number of significant decreases in mortality rate for adm1 and adm2
  print(sprintf("For the %s age bin, adm%i, %.1f%% (%i of %i) admin units decreased significantly from 2000 to 2017", age, 2, (nrow(adm2) / number_adm2s) * 100, nrow(adm2), number_adm2s))
  print(sprintf("For the %s age bin, adm%i, %.1f%% (%i of %i) admin units decreased significantly from 2000 to 2017", age, 1, (nrow(adm1) / number_adm1s) * 100, nrow(adm1), number_adm1s))
  sink()
  
  adm2_comparison[, decrease := 0]
  adm2_comparison[upper_2017 < lower_2000, decrease := 1]
  adm2_map <- data.table(adm2_comparison[, c("ADM2_CODE", "decrease")])

  
  bg_shp$admin_code <- as.integer(as.character(bg_shp$admin_code))
  #mapping_data <- suppressMessages(prep_shp_data_for_mapping(shp=bg_shp, dataset=adm2_map, merge_var='ADM2_CODE'))
  mapping_data <- merge(adm2_map, bg_shp, by.x = "ADM2_CODE", by.y = "admin_code")
  mapping_data$decrease <- as.factor(mapping_data$decrease)
  
  admin2030 <- ggplot() +
    geom_polygon(data = mapping_data,
                 aes(x=long, y=lat, group=group, fill=as.factor(decrease)), show.legend = T) + 
    geom_path(data = bg_for_mapping_adm0,
              aes(x=long, y=lat, group=group),
              color='#000000',
              lwd=0.1) +
    scale_fill_manual(values = c("grey", "#009999"), name = "", labels = c("Non-significant\ndecrease/increase", "Significant decrease"))+
    theme_map() +
    theme(legend.position = c(0.95, 0.80),
          legend.title = element_text(size=9, color='#222222'),
          legend.text  = element_text(size=6, color='#222222')) +
    # coord_map_to_bounds(shp_fort=mapping_data, projection = "mollweide") +
    # coord_cartesian(xlim = c(-117.5, 157.5), ylim = c(-35.5, 99))
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