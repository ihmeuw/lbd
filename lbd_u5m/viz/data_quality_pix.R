################################################################################
##                                                                            ## 
## DATA QUALITY PIX                                                           ##                                                  ##
## Purpose: Check the percentage of pixels with underlying data points for    ##
##            checking data quality. Ouputs a line graph by region and year   ##
##            as well as several summary tables.                              ##
##                                                                            ##
##                                                                            ##
################################################################################

################################################################################
## MAIN FUNCTION                                                              ##  
################################################################################

## data_quality_pix() ----------------------------------------------------------
#'
#' @title Data Quality Pix
#' 
#' @description Calculates the percentage of raster pixels with underlying data 
#'   points for each region over year.  
#'   
#' ** ## PARAMETERS ## **
#' @param core_repo (default `'<<<< FILEPATH REDACTED >>>>'`) The directory
#'   used to read in other MBG functions.
#' @param indicator_group The indicator group associated with this data, used to
#'  associate with a filename in the folder.
#' @param indicator The indicator associated with this data, used to associate
#'  with a filename in the folder.
#' @param data_file Prepared data for modelling that includes the latitudes and 
#'  longitudes of all data rows for comparing to raster 
#' @param shapefile_version (default `"current"`) version of the shapefile to
#'  use to load in simple raster
#' @param start_year The earliest year of data to include in plots.
#' @param end_year The latest year of data to include in plots.
#' @param regions A list of regions to check, also used to load in covariate 
#'  raster and training data.
#' @param out_dir Filepath to save summary plots and tables.
#' @param dropResample (default `FALSE`) when false this uses the raw data file 
#'  that is specified. When TRUE this subsets the data file to remove 
#'  pseudoclusters so only points and no resampled polygons are included in the 
#'  calculations.

data_quality_pix <- function(core_repo = "<<<< FILEPATH REDACTED >>>>",
                             indicator_group,
                             indicator,
                             data_file,
                             shapefile_version = "current",
                             start_year,
                             end_year,
                             regions,
                             outdir,
                             dropResample = FALSE){
  
  ##############################################################################
  ## LOAD PACKAGES
  
  message("Loading packages...")
  package_list <- c('data.table', 'raster', 'rgdal', 'rgeos', 'sp',
                    'scales', 'rlist', 'parallel', 'dplyr', 'ggplot2',
                    'RColorBrewer','plyr')
  ## Source package imports and mbg function
  source("<<<< FILEPATH REDACTED >>>>")
  ## Load external packages
  load_R_packages(package_list)
  ## Load mbg_functions
  load_mbg_functions(core_repo)
  
  ## Load data input
  ## NOTE: I do not subset this by region (not sure how to effectively), this should
  ##    be fine since the raster is subset and any data point outside the raster will
  ##    become NA which is then dropped.
  
  data <- readRDS(sprintf('<<<< FILEPATH REDACTED >>>>',
                          indicator_group,indicator,data_file))
  
  ## Remove polygon resampled data
  if (dropResample){
    data <- data[pseudocluster == FALSE]
  }
  
  ## Load in masking raster 
  maskr <- raster("<<<< FILEPATH REDACTED >>>>")
  maskr[is.na(maskr)] <- 0
  maskr[values(maskr)==1] <- NA
  
  ## Calculate time frame
  span <- (end_year - start_year) + 1
  
  region_summary <- NULL
  region_years_summary <- NULL
  all_summary <- NULL
  for (reg in regions){
    message(sprintf("\nNow processing %s...",reg))
    ## SET UP SIMPLE RASTER ####################################################
    ## Loading in the simple raster
    adm0_list      <- get_adm0_codes(reg, shapefile_version = shapefile_version)
    simple_polygon <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, 
                                          tolerance = 0.4,
                                          shapefile_version = shapefile_version)
    subset_shape   <- simple_polygon[['subset_shape']]
    raster_list    <- build_simple_raster_pop(subset_shape)
    simple_raster  <- raster_list[['simple_raster']]
    
    ## mask to remove areas we don't model 
    maskr_crop <- crop(maskr, extent(simple_raster))
    masked_raster <- mask(simple_raster,maskr_crop)  
    
    ## Get pixel count of stage2 simple raster before changing to unique vals
    masked_raster[!is.na(masked_raster)] <- 1
    total_reg_pix <- cellStats(masked_raster,'sum')
    
    ## Multiply to get pix over all years covered 
    total_reg_pix <- total_reg_pix * span 
    
    ## Assign a unique value to each pixel in the simple raster
    ## NOTE: This also assigns values to na cells 
    values(masked_raster) <- 1:ncell(masked_raster)
    
    ## PROCESS INPUT DATA ######################################################
    years_summary <- NULL
    for (i in start_year:end_year) {
      data_sub <- data[year==i,]
      year_pts <- sp::SpatialPoints(data_sub[, .(longitude, latitude)])
      
      ## Get the unique raster values under data points 
      unique_vals <- raster::extract(masked_raster,year_pts, method='simple', 
                                     na.rm=TRUE, df=TRUE)
      ras_vals <- unique_vals$layer 
      ras_vals <- ras_vals[!is.na(ras_vals)]
      
      ## Add year summary to todal list
      add_row <- data.frame(region=reg,
                            year=i,
                            total_reg_ras_pix=length(unique(ras_vals)),
                            total_reg_pix=total_reg_pix/span,
                            percent_cover=(length(unique(ras_vals)) / (total_reg_pix/span))*100)
      years_summary <- rbind(years_summary,add_row)
      region_years_summary <- rbind(region_years_summary,add_row)
      
    }# end loop through years
    
    ## Calculate regions percent cover over all years
    total_reg_ras_pix <- sum(years_summary$total_reg_ras_pix)
    reg_percent_cover <- (total_reg_ras_pix / total_reg_pix) * 100
    
    add_row <- data.frame(region=reg, 
                          total_reg_ras_pix=total_reg_ras_pix,
                          total_reg_pix=total_reg_pix,
                          reg_percent_cover=reg_percent_cover)
    region_summary <- rbind(region_summary,add_row)
    
    ## Calculate the regions percent cover over all years
    pts <- sp::SpatialPoints(data[, .(longitude, latitude)])
    vals <- raster::extract(masked_raster,pts,method='simple',na.rm=T,df=T)
    vals <- vals$layer
    vals <- vals[!is.na(vals)]
    
    all_reg_ras_pix <- length(unique(vals))
    all_reg_percent_cover <- (all_reg_ras_pix / (total_reg_pix/span)) * 100
    
    add_frame <- data.frame(region=reg,
                            all_reg_ras_pix=all_reg_ras_pix,
                            total_reg_pix=total_reg_pix/span,
                            all_percent_cover=all_reg_percent_cover)
    all_summary <- rbind(all_summary,add_frame)
    
  }# end loop through regions
  
  ## Get total for pix split by year 
  total_ras_pix <- sum(region_summary$total_reg_ras_pix)
  total_pix <- sum(region_summary$total_reg_pix)
  percent_cover <- (total_ras_pix / total_pix) * 100
  add_row <- data.frame(region='total',
                        total_reg_ras_pix=total_ras_pix,
                        total_reg_pix=total_pix,
                        reg_percent_cover=percent_cover)
  region_summary <- rbind(region_summary,add_row)
  
  ## Get total for pix across all years
  all_ras_pix <- sum(all_summary$all_reg_ras_pix)
  all_pix <- sum(all_summary$total_reg_pix)
  all_percent_cover <- (all_ras_pix / all_pix) * 100
  add_frame <- data.frame(region='total',
                          all_reg_ras_pix=all_ras_pix,
                          total_reg_pix=all_pix,
                          all_percent_cover=all_percent_cover)
  all_summary <- rbind(all_summary,add_frame)
  
  ## PLOT YEARLY CHANGES BY REGION ###############################################
  
  ## Rename the regions to something understandable in all tables
  
  key <- c('ansa+trsa-bra', 'caca-mex', 'cssa', 'essa-yem', 'mide+yem', 'noaf',
           'ocea+seas-mys', 'soas', 'sssa', 'stan+mng', 'wssa')
  val <- c('South America', 'Central America & Caribbean', 'Central Sub-Saharan Africa',
           'Eastern Sub-Saharan Africa', 'Middle East', 'North Africa',
           'Oceania & Southeast Asia', 'South Asia', 'Southern Sub-Saharan Africa',
           'Central & East Asia', 'Western Sub-Saharan Africa')
  region_years_summary$region <- mapvalues(region_years_summary$region, from=key, to=val)
  region_summary$region <- mapvalues(region_summary$region, from=key, to=val)
  all_summary$region <- mapvalues(all_summary$region, from=key, to=val)
  
  ## Reorder the regions for the map 
  region_years_summary$region <- factor(region_years_summary$region, 
                                        levels=c('Central America & Caribbean',
                                                 'South America', 'North Africa',
                                                 'Western Sub-Saharan Africa',
                                                 'Eastern Sub-Saharan Africa',
                                                 'Central Sub-Saharan Africa',
                                                 'Southern Sub-Saharan Africa',
                                                 'Middle East', 'Central & East Asia',
                                                 'South Asia', 
                                                 'Oceania & Southeast Asia'))
  
  ## Create a color palette 
  numregs <- length(regions)
  getPalette <- colorRampPalette(brewer.pal(n=8, name="Accent"))
  
  fig <- ggplot(region_years_summary, aes(x=year, y=percent_cover, group=region, color=region)) +
    geom_line() + 
    geom_point(size=1.65, shape=18) +
    scale_color_manual(values=getPalette(numregs)) + 
    scale_x_continuous(expand=c(0,0), limits=c(2000,2017.5), breaks=seq(2000,2017,1)) +
    scale_y_continuous(expand=c(0,0), limits=c(0,20), labels = function(x) paste0(x,"%")) +
    labs(title='Data Per Pixel by Region',
         y='Proportion of Pixels with Overlapping Data in Year',
         x='Year Modelled',
         color='Region') +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  ## SAVE ########################################################################
  
  if (dropResample) {
    write.csv(region_summary,sprintf('<<<< FILEPATH REDACTED >>>>',outdir),
              row.names = FALSE)
    write.csv(region_years_summary,sprintf('<<<< FILEPATH REDACTED >>>>',outdir),
              row.names = FALSE)
    write.csv(all_summary,sprintf('<<<< FILEPATH REDACTED >>>>',outdir),
              row.names = FALSE)
    
    ## hi res
    png(sprintf('<<<< FILEPATH REDACTED >>>>',outdir), 
        height=5.5, width=12, units='in', res=1000)
    print(fig)
    dev.off()
  } else {
    write.csv(region_summary,sprintf('<<<< FILEPATH REDACTED >>>>',outdir),
              row.names = FALSE)
    write.csv(region_years_summary,sprintf('<<<< FILEPATH REDACTED >>>>',outdir),
              row.names = FALSE)
    write.csv(all_summary,sprintf('<<<< FILEPATH REDACTED >>>>',outdir),
              row.names = FALSE)
    
    ## hi res
    png(sprintf('<<<< FILEPATH REDACTED >>>>',outdir), 
        height=5.5, width=12, units='in', res=1000)
    print(fig)
    dev.off()
  }#end save if statements
  
}#end main function data_quality_pix()
