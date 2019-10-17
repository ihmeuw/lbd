## #############################################################################
## CALCULATE LIVES SAVED DUE TO REDUCTION IN MORTALITY BY PIXEL/SUBNATIONAL
##
## Purpose: Calculate the number of deaths that *would have* occurred if the
##   mortality rate remained constant in each pixel since 2000. Using this
##   estimate, determine the number of lives saved in each pixel due to
##   reduction in mortality.
## 
## #############################################################################

rm(list=ls())

## Imports
library(data.table)
library(raster)
library(rgdal)

## Source MBG and U5M functions
repo_path <- '<<<< FILEPATH REDACTED >>>>'
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
source("<<<< FILEPATH REDACTED >>>>")

## DEFINE INPUT VALUES
run_date     <- '<<<< REDACTED >>>>'
group <- 'under5'
model_start_year <- 2000 # MODEL start year
model_end_year   <- 2017 # MODEL end year

save_folder  <- '<<<< FILEPATH REDACTED >>>>'
sub_folder   <- NULL # Default is a timestamp
country_file <- paste0(save_folder,'<<<< FILEPATH REDACTED >>>>')
adm2_file    <- paste0(save_folder,'<<<< FILEPATH REDACTED >>>>')
merged_file  <- paste0(save_folder, '<<<< FILEPATH REDACTED >>>>')
resume       <- TRUE

# Plotting information:
shp_dir <- paste0(
  '<<<< FILEPATH REDACTED >>>>',
  '<<<< FILEPATH REDACTED >>>>'
)
shp_layer <- '<<<< FILEPATH REDACTED >>>>'

plot_min <- -5
plot_max <- 15
plot_min_allyear <- -5
plot_max_allyear <- 50


## Load and prep input dataset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Preliminaries:
# - Prep input folder
if(is.null(sub_folder)) sub_folder <- gsub('-', '_', Sys.Date())
save_folder_full <- paste0(save_folder, sub_folder)
if( !file.exists(save_folder_full) ) dir.create(save_folder_full)
# - Get list of all admin codes in Stage 2
lookup_table <- fread('<<<< FILEPATH REDACTED >>>>')
stage2_gauls <- unique(lookup_table[Stage !='3', GAUL_CODE])


## Function for pulling adm(X) rasters for Stage 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pull_adm_units <- function(save_folder, resume, gauls, compare_raster){
  message("Assembling administrative rasters...")
  out_list <- list('adm0'=NA, 'adm1'=NA, 'adm2'=NA)
  for(lev in 0:2){
    backup_file <- paste0(save_folder,'<<<< FILEPATH REDACTED >>>>',lev,'<<<< FILEPATH REDACTED >>>>')
    if((resume==TRUE) & (file.exists(backup_file))){
      # Get admin file from backup
      this_adm_ras <- get(load(backup_file))
    } else {
      # Create adm(X) rasters
      this_adm_ras <- GetAdmin(0,compare_raster,gauls)[['rast']]
    }
    # Ensure that the comparison raster has the same dimensions and extent
    if( !all.equal(dim(this_adm_ras),dim(compare_raster)) ){
      stop(paste0("Admin",lev," raster does not have the correct dimensions."))
    }
    # Add to list
    out_list[[paste0('adm',lev)]] <- this_adm_ras
  }
  # Return all levels in a list
  return(out_list)
}


if((resume==TRUE) & (file.exists(merged_file))){
  ## The data has already been prepped
  full_data_list <- get(load(merged_file))
  # Extract individual objects
  full_data      <- full_data_list[['full_data']]
  country_raster <- full_data_list[['country_raster']]
  idx            <- full_data_list[['idx']]
} else {
  ## The data needs to be prepped
  ## 1) Load in 5q0 and 5D0 for all years from the model folder
  in_dir <- sprintf('<<<< FILEPATH REDACTED >>>>',group,run_date)
  D <- brick(sprintf('<<<< FILEPATH REDACTED >>>>',
                     in_dir, group, model_start_year, model_end_year))
  Q <- brick(sprintf('<<<< FILEPATH REDACTED >>>>',
                     in_dir, group, model_start_year, model_end_year))
  # Ensure that these bricks have the same extent
  if(extent(D)!=extent(Q)){
    stop("The input raster bricks do not have the same dimensions - please fix.")
  }

  ## 2) Load in the under-5 population covariate
  pop <- load_and_crop_covariates_annual(
    covs           = 'worldpop',
    measures       = 'a0004t',
    simple_polygon = D,
    start_year     = model_start_year,
    end_year       = model_end_year,
    interval_mo    = as.numeric(12),
    agebin=1
  )[[1]]
  
  ## 3) Merge on administrative metadata
  adm_rasters <- pull_adm_units(
    save_folder    = save_folder, 
    resume         = resume, 
    gauls          = stage2_gauls, 
    compare_raster = D
  )
  
  ## 4) Convert 5q0 to 5m0 for the first year of data
  # Get GBD nAx values
  gbd_nax <- get_gbd_crosswalk_values(
    age       = group,
    year_list = model_start_year
  )
  # Using the raster and the GBD nAx values, pull 5m0 for the first year ONLY
  u5mr_firstyear <- crosswalk_qm(
    data      = brick(Q[[1]]), 
    data_type = 'raster', 
    region    = 'stage2',
    year_list = model_start_year,
    age       = group,
    gbd       = gbd_nax,
    simple_raster = country_raster[[1]]
  )
  # Ensure that the 5m0 data has the same dimensions as Q and D
  if( !all.equal(dim(u5mr_firstyear)[1:2], dim(Q)[1:2]) ){
    stop("The U5MR raster does not have the same extent as D and Q.")
  }
  
  ## 5) MERGE EVERYTHING TOGETHER
  full_data <- data.table(
    country = as.vector(adm_rasters[['adm0']]),
    adm1    = as.vector(adm_rasters[['adm1']]),
    adm2    = as.vector(adm_rasters[['adm2']]),
    d_true  = as.vector(D),
    q_true  = as.vector(Q),
    pop     = as.vector(pop),
    m_fy    = rep(as.vector(u5mr_firstyear[[1]]), 
                  length(model_start_year:model_end_year)),
    year    = rep(model_start_year:model_end_year,
                  each=prod(dim(Q)[1:2]))
  )
  # Ensure that the data.table has the expected number of rows
  if( nrow(full_data) != prod(dim(Q))){
    stop("The full merged data.table does not have the correct number of rows.")
  }
  # Drop rows where the simple raster values are missing
  full_data <- full_data[!is.na(country),]
  # Get the cellIdx of the simple raster for easy raster recreation later
  idx <- cellIdx(country_raster)
  
  ## 6) SAVE BACKUP FILE
  full_data_list <- list(
    full_data      = full_data,
    country_raster = country_raster,
    idx            = idx
  )
  save(full_data_list, file=merged_file)
  # Remove intermediate objects from memory
  rm(list=c('D','Q','pop','u5mr_firstyear'))
}

## Prep counterfactual deaths by year ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate the number of deaths that would have occurred if the mortality rate
#  had not decreased since 2000
full_data[, d_nochange := m_fy * pop ]

# Calculate the number of lives saved by a decreasing death rate:
# Lives saved = D(counterfactual) - D(actual)
full_data[, lives_saved := d_nochange - d_true]
# quantile(full_data[,lives_saved], c(.001, .05, .1, .25, .5, .75, .9, .95, .999))

# ISSUE WITH LATER YEARS -- DROP DATA AFTER 2015 FOR NOW
data_sub <- full_data[ (year >= 2001) & (year <= 2015), ]

## #############################################################################
## PLOT LIVES SAVED 
## #############################################################################

# Convert to raster
lives_saved_ras <- insertRaster(
  raster = country_raster[[1]],
  new_vals = matrix(
    data  = as.vector(data_sub[,lives_saved]),
    ncol  = length(2001:2015),
    byrow = FALSE),
  idx = idx
)
lives_saved_ras <- brick(lives_saved_ras)


## Plotting preparation
all_years <- min(data_sub[,year]):max(data_sub[,year])
plotting_brick <- copy(lives_saved_ras)
plotting_brick[plotting_brick < plot_min] <- plot_min
plotting_brick[plotting_brick > plot_max] <- plot_max
# Read shapefile and subset to Stage 2
adm0_shp <- readOGR(shp_dir, shp_layer)
stage2_shp <- adm0_shp[ adm0_shp$ADM0_CODE %in% stage2_gauls , ]

## Plot lives saved *by year* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Define color scale and breaks
breaks <- c(seq(plot_min, 0, length.out=6),
            seq(1, plot_max, length.out=15)) # 21 breaks
worst_col <- "#CC0000"
ramp1 <- colorRampPalette(c("#FF0066","#FFFFCC"))
ramp2 <- colorRampPalette(c("#FFFFCC","#0066CC"))
best_col <- "#003399"
all_cols = c( worst_col, ramp1(4), ramp2(14), best_col ) # 20 colors

# Open PDF and plot
pdf(
  paste0(save_folder_full ,'<<<< FILEPATH REDACTED >>>>'),
  width=15, height=8
)
for(y in 1:length(all_years)){
  plot_year <- all_years[y]
  plot_ras  <- plotting_brick[[y]]

  # Make the plot  
  raster::plot(
    plot_ras,
    axes        = FALSE,
    breaks      = breaks,
    col         = all_cols,
    maxpixels   = length(plot_ras),
    legend      = TRUE,
    main        = paste0("Lives saved due to decreased mortality in ",plot_year),
    axis.args   = list(
      at       = c(plot_min_allyear, seq(0, plot_max, length.out=4)),
      labels   = c(plot_min_allyear, seq(0, plot_max, length.out=4)),
      cex.axis = .6
    ),
    legend.args = list(text='Lives\nsaved\nin\npixel', cex=.8)
  )
  lines(
    stage2_shp,
    lwd=.4,
    col='#666666'
  )
}

dev.off()

## Plot lives saved *across all years* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determine lives saved across all years
all_years_ras <- sum(lives_saved_ras)
all_years_ras[ all_years_ras < plot_min_allyear ] <- plot_min_allyear
all_years_ras[ all_years_ras > plot_max_allyear ] <- plot_max_allyear

# Define color scheme
breaks_allyr   <- c(seq(plot_min_allyear, 0, length.out=6),
                    seq(1, plot_max_allyear, length.out=50)) # 56 breaks
all_cols_allyr <- c( worst_col, ramp1(4), ramp2(49), best_col ) # 55 colors

# Plot the raster
pdf(
  paste0(save_folder_full ,'<<<< FILEPATH REDACTED >>>>'),
  width=15, height=8
)
raster::plot(
  all_years_ras,
  axes        = FALSE,
  breaks      = breaks_allyr,
  col         = all_cols_allyr,
  maxpixels   = length(all_years_ras[[1]]),
  legend      = TRUE,
  main        = paste0("Lives saved due to decreased mortality, 2001-2015"),
  axis.args   = list(
    at       = c(plot_min_allyear, seq(0, plot_max_allyear, length.out=3)),
    labels   = c(plot_min_allyear, seq(0, plot_max_allyear, length.out=3)),
    cex.axis = .6
  ),
  legend.args = list(text='Lives\nsaved\nin\npixel', cex=.8)
)
lines(
  stage2_shp,
  lwd=.4,
  col='#666666'
)
dev.off()