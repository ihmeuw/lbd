## ###########################################################################
## 
## PREPARE IRAN MAPS FOR A. KHOSRAVI
## 
## Purpose: Prep Iran SBH data at admin2 level
## 
## ###########################################################################

rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(rgeos)
library(rgdal)
library(raster)
library(RColorBrewer)
# MBG covariate functions
source("<<<< FILEPATH REDACTED >>>>")

## Define inputs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Input and output files
# OUTPUT DIRECTORY AND PDF
main_dir <- '<<<< FILEPATH REDACTED >>>>'
sub_dir <- paste0(main_dir, gsub('-', '_', Sys.Date()) )
dir.create(sub_dir, showWarnings = FALSE)
out_file_base <- paste0(sub_dir, '<<<< FILEPATH REDACTED >>>>')

# Defines larger folder of SBH input data
sbh_run_date <- '<<<< REDACTED >>>>'
# Defines subfolders of SBH input data
year_to_nid <- list(
  '2006' = 39396,
  '2011' = 81291,
  '2016' = 299134
)
# Shapefile locations
shp_dir <- '<<<< FILEPATH REDACTED >>>>'
province_shp_dir <- paste0('<<<< FILEPATH REDACTED >>>>',
                           '<<<< FILEPATH REDACTED >>>>')
province_shp_layer <- '<<<< FILEPATH REDACTED >>>>'
# Location of rasterized polygon
ras_poly_file <- '<<<< FILEPATH REDACTED >>>>'

# MODEL SETTINGS TO LOAD Q AND P
modeled_q_yl <- list(
  '2005'=2005,
  '2010'=2010,
  '2015'=2015
)
rd    <- '<<<< REDACTED >>>>'
group <- 'under5'
model_start_year <- 2000 # MODEL start year
model_end_year   <- 2017 # MODEL end year
model_dir      <- sprintf('<<<< FILEPATH REDACTED >>>>',group,rd)
q_file_unraked <- sprintf('<<<< FILEPATH REDACTED >>>>',
                          model_dir, group, model_start_year, model_end_year)
q_file_raked   <- sprintf('<<<< FILEPATH REDACTED >>>>',
                        model_dir, group, model_start_year, model_end_year)

# Default color scheme limits
col_min <- 0.0
col_max <- 0.08

## Function to pull raw data from flat file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pull_raw_data <- function(keep_nid, dataset){
  # Subset to the given NID
  data_sub <- dataset[nid==keep_nid,]
  # Keep only rows where both CEB and CED is known
  data_sub <- data_sub[!is.na(ceb) & !is.na(ced),]
  # Keep only women under 30 for now (close to present haz value)
  data_sub <- data_sub[age_year <= 30,]
  # Sum CEB and CED by location code
  data_agg <- data_sub[, .(ceb = sum(ceb), ced=sum(ced)), 
                       by=.(nid, country, year, shapefile, location_code)]
  # Keep only rows where CEB > 0
  data_agg <- data_agg[ceb > 0,]
  # haz = CED/CEB
  data_agg[, haz := ced/ceb]
  # Return aggregated data
  return(data_agg)
}

## Function to pull SBH-corrected aggregated data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pull_sbh_corrected_by_admin <- function(nid, run_date){
  # Define model directory
  model_dir <- paste0('<<<< FILEPATH REDACTED >>>>',
                      '<<<< FILEPATH REDACTED >>>>',run_date,'/',nid,'/')
  # Get all files to read in
  in_files <- paste0(
    model_dir, 
    grep(
      '^prepped_yearly_aggregation_summary_',list.files(model_dir), value = TRUE
    )
  )
  # Read in all files and combine row-wise
  all_data <- lapply(as.list(in_files), function(x) readRDS(x)) %>% rbindlist

  # Calculate 5q0 across all age groups, keeping only groups that have all 7
  #  age groups
  all_data <- all_data[, .(numbins = .N, haz=1-prod(1-haz)),
                       by=.(year_entering,shapefile,location_code)]
  all_data <- all_data[numbins == 7,]
  # Return data
  return(all_data)
}


## Function to aggregate Iran MBG estimates to a new admin shapefile ~~~~~~~~~
agg_mbg_results <- function(
  mbg_raster,
  pop_raster,
  rasterized_shp,
  first_modeled_year,
  year
  ){
  # Get index of this year
  year_i <- year - first_modeled_year + 1
  # Get rasters for the year of data
  q_yr <- mbg_raster[[year_i]]
  p_yr <- pop_raster[[year_i]]
  # Mask pop and MBG rasters to the shp admin2 outlines
  q_yr <- raster::mask(x=q_yr, mask=rasterized_shp)
  p_yr <- raster::mask(x=p_yr, mask=rasterized_shp)
  adm_extent <- extent(rasterized_shp)
  q_yr <- setExtent(q_yr, adm_extent)
  p_yr <- setExtent(p_yr, adm_extent)
  # Convert all three to a single data.table
  dt <- data.table(
    q = as.vector(q_yr),
    p = as.vector(p_yr),
    GAUL_CODE = as.integer(as.vector(rasterized_shp))
  )
  dt <- na.omit(dt)
  # Aggregate results by admin2 code, using population-weighted mean
  adm_agg <- dt[, .(haz=weighted.mean(x=q, w=p)), by=.(GAUL_CODE)]
  adm_agg <- adm_agg[order(GAUL_CODE)]
  # Return data.table
  return(adm_agg)
}


## Styling for map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
      plot.background = element_rect(fill = "#f5f5f2", color = NA), 
      panel.background = element_rect(fill = "#f5f5f2", color = NA), 
      legend.background = element_rect(fill = "#f5f5f2", color = NA),
      panel.border = element_blank(),
      ...
    )
}


## Function to make a map of Iran U5M given input data and labels ~~~~~~~~~~~~
iran_u5m_map <- function(
  province_shp,
  county_shp,
  col_min,
  col_max,
  title,
  subtitle,
  legend_title = 'Under-5\nmortality\nprobability\n(5q0)'
  ){
  # Restrict colors to (min, max)
  county_shp[ haz > col_max, haz := col_max]
  county_shp[ haz < col_min, haz := col_min]
  # Make map
  out_map <- ggplot() + 
    theme_map() +
    geom_polygon(
      data=county_shp,
      aes(x=long, y=lat, group=group, fill=haz),
      color='#DDDDDD', lwd=.25
    ) + 
    geom_polygon(
      data=province_shp,
      aes(x=long, y=lat, group=group),
      fill=NA, color='#222222'
    ) + 
    labs(
      title = title,
      subtitle = subtitle,
      fill = legend_title
    ) + 
    scale_fill_gradientn(
      limits   = c(col_min, col_max),
      colors   = rev(brewer.pal(name='Spectral', n=11)),
      guide    = 'colorbar',
      na.value = '#FFFFFF'
    )
  # Return map as ggplot object
  return(out_map)
}




## 1. PREP ALL DATASETS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load in province-level shapefile
province_shp <- readOGR(province_shp_dir,province_shp_layer)
province_fortified <- fortify(province_shp)

## Prep raw data
raw_data_full <- readRDS(
  paste0('<<<< FILEPATH REDACTED >>>>',
         '<<<< FILEPATH REDACTED >>>>')
)
raw_by_year <- lapply(
  year_to_nid,
  function(y) pull_raw_data(keep_nid=y, dataset=raw_data_full)
)

## Prep SBH-corrected data
corrected_data_by_year <- lapply(
  year_to_nid,
  function(y) pull_sbh_corrected_by_admin(nid=y, run_date=sbh_run_date)
)

## Prep modeled data (unraked and raked)
# Prep a single shapefile which will be joined to all years of data
model_shp_layer <- '<<<< FILEPATH REDACTED >>>>'
model_shp <- readRDS(paste0(shp_dir,model_shp_layer,'.rds'))
model_shp@data$GAUL_CODE <- as.integer(model_shp@data$GAUL_CODE)
model_shp$id <- sapply(model_shp@polygons, function(x) x@ID)

# Load rasterized shapefile, population, and 5q0 file
rasterized_shp   <- get(load(ras_poly_file))
q_raster_unraked <- brick(q_file_unraked)
q_raster_raked   <- brick(q_file_raked)

pop_raster <- load_and_crop_covariates_annual(
  covs           = 'worldpop',
  measures       = 'a0004t',
  simple_polygon = q_raster_unraked,
  start_year     = model_start_year,
  end_year       = model_end_year,
  interval_mo    = as.numeric(12),
  agebin=1
)[[1]]

# UNRAKED DATA: Create aggregated q estimates by admin unit
modeled_q_list_unraked <- lapply(
  modeled_q_yl,
  function(y) agg_mbg_results(
    mbg_raster         = q_raster_unraked,
    pop_raster         = pop_raster,
    rasterized_shp     = rasterized_shp,
    first_modeled_year = model_start_year,
    year               = y
  )
)
# RAKED DATA: Create aggregated q estimates by admin unit
modeled_q_list_raked <- lapply(
  modeled_q_yl,
  function(y) agg_mbg_results(
    mbg_raster         = q_raster_raked,
    pop_raster         = pop_raster,
    rasterized_shp     = rasterized_shp,
    first_modeled_year = model_start_year,
    year               = y
  )
)

## 2. PLOT ALL DATA FOR THREE YEARS: 2005, 2010, and 2015 ~~~~~~~~~~~~~~~~~~~~

for(this_year in c(2005, 2010, 2015)){
  # Open a new PDF
  pdf(sprintf(out_file_base, this_year), width=10, height=8)
  
  ## 2A. PLOT RAW DATA FROM THE CENSUS IN YEAR+1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  year_raw <- raw_by_year[[as.character(this_year+1)]]
  year_raw[, GAUL_CODE := as.integer(location_code) ]
  # Read in shapefile
  raw_shp_layer <- gsub("^irn_", "iran_", year_raw[1, shapefile])
  raw_shp <- readRDS(paste0(shp_dir,raw_shp_layer,'.rds'))
  # Merge dataset onto shapefile
  raw_shp@data$GAUL_CODE <- as.integer(raw_shp@data$GAUL_CODE)
  raw_shp$id <- sapply(raw_shp@polygons, function(x) x@ID)
  # "Fortify" shape for mapping
  raw_shp_fortified <- fortify(raw_shp) %>% mutate(id = as.numeric(id)) %>% as.data.table
  setnames(raw_shp_fortified, 'order', 'shp_order')
  raw_shp_fortified[, id := as.character(id)]
  original_length <- dim(raw_shp_fortified)[1]
  raw_shp_fortified <- merge(
    x     = raw_shp_fortified, 
    y     = raw_shp,
    by    = c('id'),
    all.x = TRUE
  )
  raw_shp_fortified <- merge(
    x     = raw_shp_fortified,
    y     = year_raw,
    by    = c('GAUL_CODE'),
    all.x = TRUE
  )
  final_length <- dim(raw_shp_fortified)[1]
  if(final_length > original_length){
    stop("There is an issue with the merges")
  }
  raw_shp_fortified <- raw_shp_fortified[order(shp_order),]
  # CREATE MAP OF RAW DATA
  raw_shp_map <- iran_u5m_map(
    province_shp = province_fortified,
    county_shp   = raw_shp_fortified,
    col_min      = col_min,
    col_max      = col_max,
    title        = paste0("Raw Calculated Child Mortality Data in Iran (",
                          this_year+1," Census)"),
    subtitle     = paste0("Estimated SBH calculated among all women under the age of 30"),
    legend_title = "CED/CEB"
  )
  print(raw_shp_map)
  
  ## 2B. PLOT SBH-CORRECTED DATA FROM ALL FUTURE CENSUSES ~~~~~~~~~~~~~~~~~~~~
  get_census_for_year <- function(dt, which_year){
    dt[,census_year := max(year_entering) + 1]
    dt <- dt[year_entering==which_year,]
    return(dt)
  }
  all_census_data <- lapply(
    corrected_data_by_year,
    function(dt) get_census_for_year(dt, which_year=this_year)
  ) %>% rbindlist
  all_census_data[, GAUL_CODE := as.integer(location_code) ]
  
  for(future_year in unique(all_census_data[,census_year])){
    future_census <- all_census_data[census_year==future_year,]
    # Read in relevant shapefile, fixing name if necessary
    corr_shp_layer <- gsub("^irn_", "iran_", future_census[1, shapefile])
    corr_shp <- readRDS(paste0(shp_dir,corr_shp_layer,'.rds'))
    # Merge dataset onto shapefile
    corr_shp@data$GAUL_CODE <- as.integer(corr_shp@data$GAUL_CODE)
    corr_shp$id <- sapply(corr_shp@polygons, function(x) x@ID)
    # "Fortify" shape for mapping
    corr_shp_fortified <- fortify(corr_shp) %>% mutate(id = as.numeric(id)) %>% as.data.table
    setnames(corr_shp_fortified, 'order', 'shp_order')
    corr_shp_fortified[, id := as.character(id)]
    original_length <- dim(corr_shp_fortified)[1]
    corr_shp_fortified <- merge(
      x     = corr_shp_fortified, 
      y     = corr_shp,
      by    = c('id'),
      all.x = TRUE
    )
    corr_shp_fortified <- merge(
      x     = corr_shp_fortified,
      y     = future_census,
      by    = c('GAUL_CODE'),
      all.x = TRUE
    )
    final_length <- dim(corr_shp_fortified)[1]
    if(final_length > original_length){
      stop("There is an issue with the merges")
    }
    corr_shp_fortified <- corr_shp_fortified[order(shp_order),]
    # CREATE MAP
    corr_shp_map <- iran_u5m_map(
      province_shp = province_fortified,
      county_shp   = corr_shp_fortified,
      col_min      = col_min,
      col_max      = col_max,
      title        = paste0("SBH Adjusted Child Mortality Data in Iran, ",
                            this_year," (adapted from ",future_year,
                            " census)"),
      subtitle     = paste0("Data adjusted using SBH correction method from ",
                            "Burstein et al. (forthcoming)")
    )
    print(corr_shp_map)
  }

  ## 2C. Plot modeled results for this year (unraked and raked) ~~~~~~~~~~~~~~
  q_unraked_year <- modeled_q_list_unraked[[ as.character(this_year) ]]
  q_raked_year   <- modeled_q_list_raked[[ as.character(this_year) ]]
  
  # Standard shapefile fortify for mapping
  model_shp_fortified <- fortify(model_shp) %>% mutate(id = as.numeric(id)) %>% as.data.table
  setnames(model_shp_fortified, 'order', 'shp_order')
  model_shp_fortified[, id := as.character(id)]
  original_length <- dim(model_shp_fortified)[1]
  model_shp_fortified <- merge(
    x     = model_shp_fortified, 
    y     = model_shp,
    by    = c('id'),
    all.x = TRUE
  )
  # Merge with UNRAKED data
  unraked_model_shp_fortified <- merge(
    x     = model_shp_fortified,
    y     = q_unraked_year,
    by    = c('GAUL_CODE'),
    all.x = TRUE
  )
  final_length_unraked <- dim(unraked_model_shp_fortified)[1]
  if(final_length_unraked > original_length){
    stop("There is an issue with the merges")
  }
  unraked_model_shp_fortified <- unraked_model_shp_fortified[order(shp_order),]
  # Merge with RAKED data
  raked_model_shp_fortified <- merge(
    x     = model_shp_fortified,
    y     = q_raked_year,
    by    = c('GAUL_CODE'),
    all.x = TRUE
  )
  final_length_raked <- dim(raked_model_shp_fortified)[1]
  if(final_length_raked > original_length){
    stop("There is an issue with the merges")
  }
  raked_model_shp_fortified <- raked_model_shp_fortified[order(shp_order),]
  
  # Create UNRAKED map
  unraked_shp_map <- iran_u5m_map(
    province_shp = province_fortified,
    county_shp   = unraked_model_shp_fortified,
    col_min      = col_min,
    col_max      = col_max,
    title        = paste0("Modeled Child Mortality Data in Iran, ",this_year),
    subtitle     = paste0("5q0 estimated at the 5x5 km pixel level using a ",
                          "spatio-temporal model adapted from Golding et al. (2017)")
  )
  print(unraked_shp_map)

  # Create RAKED map
  raked_shp_map <- iran_u5m_map(
    province_shp = province_fortified,
    county_shp   = raked_model_shp_fortified,
    col_min      = col_min,
    col_max      = col_max,
    title        = paste0("Modeled and Raked Child Mortality Data in Iran, ",this_year),
    subtitle     = paste0("5q0 estimated at the 5x5 km pixel level and raked ",
                          "to GBD 2017 province-level estimates")
  )
  print(raked_shp_map)

  # Close PDF
  dev.off()
  message(paste0("Mapping finished for year ",this_year,".\n"))
}
