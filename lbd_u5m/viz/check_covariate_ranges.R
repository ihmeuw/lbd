################################################################################
##                                                                            ## 
## CHECK COVARIATE RANGES                                                     ##                                                  ##
## Purpose: Summarizes where areas of standard LBD covariates are outide of   ##
##            the range of team specific training data.                       ##
##                                                                            ##
## SECTIONS:                                                                  ##
##   1: MAIN FUNCTION (check_covariate_ranges)                                ##
##   2: DATA PROCESSING FUNCTIONS                                             ##
##                                                                            ##
################################################################################

################################################################################
## 1: MAIN FUNCTION                                                           ##  
################################################################################

## check_covariate_ranges() ----------------------------------------------------
#'
#' @title Check Covariate Ranges
#' 
#' @description Checks covariate raster layer values against training data 
#'   predicter values to determine where and to what extent the raster values 
#'   are out of range of the training data.
#'   
#' @details This function utilizes mclapply to paralilize the work. Cores are 
#'  not specified by use, instead the best number of threads to use is
#'   determined by get_max_forked_threads and the user needs to specifiy the
#'   number of omp threads to use when signing into R (ie. sing -o ##).
#'   
#' ** ## PARAMETERS ## **
#' @param core_repo (default `'<<<< FILEPATH REDACTED >>>>'`) The directory
#'   used to read in other MBG functions.
#' @param indicator_group The indicator group associated with this data, used to
#'  associate with a filename in the folder.
#' @param indicator The indicator associated with this data, used to associate
#'  with a filename in the folder.
#' @param run_date The date and time for the training data, used to associate 
#'  with a filename in the folder. 
#' @param start_year The earliest year of data to include in plots.
#' @param end_year The latest year of data to include in plots.
#' @param interval_mo The number of months in an intervak to load and crop 
#'  covariate raster. Currently load_and_crop_covariates_annual set up to do 
#'  yearly with intervals 12, 24, or 60.
#' @param regions A list of regions to check, also used to load in covariate 
#'  raster and training data.
#' @param remove_cov (default `year_cov`) list of any covariates found in the 
#' model that will not be used in calculations or plotting
#' @param out_dir Filepath to save summary plots and tables.
#' 
#' ** ## PLOTTING ARGUMENTS ## **
#' @param make_region_cov_maps (default `TRUE`) makes a map showing each the 
#'  number of years a pixel was out of range for each covariate for each 
#'  individual region 
#' @param make_regions_cov_hists (default `TRUE`) makes a histogram showing the 
#'  range of each covariate raster used with the range of values in the training 
#'  data used 
#' @param buffer (default `10`)Choose the buffer size for the background map 
#'  to each plot
#' @param point_size (default `0.08`) Choose the point size of geom_point for 
#'  mapping ggplot final all region map
#' @param point_stroke_size (default `0.22`) Choose the stroke suze of geom_point 
#'  for mapping ggplot final all region map
#' @param projection (default `mercator`)Choose the projection for mapping
#' 
#' @return returns `NULL`, all output is saved directly in specified out_dir 
#' 
check_covariate_ranges <- function(core_repo = '<<<< FILEPATH REDACTED >>>>',
                                   indicator_group,
                                   indicator,
                                   run_date,
                                   start_year,
                                   end_year,
                                   interval_mo,
                                   regions,
                                   remove_cov = c('year_cov'),
                                   out_dir,
                                   
                                   make_region_cov_maps = TRUE,
                                   make_region_cov_hists = TRUE,
                                   buffer = 10,
                                   point_size = 0.08,
                                   point_stroke_size = 0.22,
                                   projection = 'mercator') {
    
  ##############################################################################
  ## LOAD PACKAGES AND SETUP                                                  ##
  ##############################################################################
  
  ##############################################################################
  ## LOAD PACKAGES

  message("Loading packages...")
  package_list <- c('data.table', 'raster', 'rgdal', 'rgeos', 'sp', 'ggplot2',
                    'RColorBrewer', 'scales', 'rlist', 'parallel', 'dplyr',
                    'ggthemes', 'tictoc')
  ## Source package imports and mbg function
  source("<<<< FILEPATH REDACTED >>>>")
  ## Load external packages
  load_R_packages(package_list)
  ## Load mbg_functions
  load_mbg_functions(core_repo)
  
  ##############################################################################
  ## MISC SET UP 
  
  ## Create color palette 
  getPalette <- colorRampPalette(brewer.pal(n=9, name="YlGnBu"))
  
  ## Load master shapefile for background map
  message("\nLoading master shapefile... ")
  master_shape_all <- shapefile('<<<< FILEPATH REDACTED >>>>')
  
  ## Get disputed border shapefile for masking
  disputed_shp <- shapefile('<<<< FILEPATH REDACTED >>>>')
  
  ## Create list of region_covariates 
  region_effect_list <- list()
  for (reg in regions) {
    ## load in covariate names
    load(sprintf('<<<< FILEPATH REDACTED >>>>',
                 indicator_group, indicator, run_date, reg))
    effects <- trim(strsplit(fixed_effects, "\\+")[[1]])
    measures <- trim(strsplit(fixed_effects_measures, "\\+")[[1]])
    
    ## Remove covariates if specified
    if (!is.null(remove_cov)) {
      for (cov in remove_cov){
        measures <- measures[-match(cov,effects)]
        effects <- effects[effects != cov]
      }
    }#end segment to remove covs
   
    ## load training data to find effects used in models
    training_data <- fread(sprintf("<<<< FILEPATH REDACTED >>>>",
                                   indicator_group, indicator, run_date, reg))
    for (effect in effects){
      if (effect %in% names(training_data)){
        measure <- measures[match(effect,effects)]
        region_effect_list <- list.append(region_effect_list,c(reg, effect, 
                                                               measure))
      }# if statement to only add effects that are in the training data
    }# end effect loop
  }# end region loop

  ##############################################################################
  ## PROCESS REGION COVARIATES                                                ##
  ##############################################################################
  if (length(region_effect_list)==1) {
    region_output_list <- process_region_effects(region_effect_list[[1]],
                                                 start_year,
                                                 end_year,
                                                 interval_mo,
                                                 indicator_group,
                                                 indicator,
                                                 run_date,
                                                 out_dir,
                                                 master_shape_all,
                                                 getPalette,
                                                 make_region_cov_maps,
                                                 make_region_cov_hists)
    region_output_list <- list(region_output_list)
  } else {
    
    message("\nProcessing region effects and saving to outfile directory...")
    ## Set serial execution, perform mclapply(), and then go back to parallel 
    ## threads
    set_serial_threads()
    ## Check how many threads we should use, based on a function in
    num_threads <- get_max_forked_threads( nobjs = length(region_effect_list) )
    message(sprintf('We will run an mclapply with %s threads',num_threads))
    
    tic(msg="setup mclapply")
    region_output_list <- mclapply(region_effect_list, 
                                   function(x){
                                     process_region_effects(x,
                                                            start_year,
                                                            end_year,
                                                            interval_mo,
                                                            indicator_group,
                                                            indicator,
                                                            run_date,
                                                            out_dir,
                                                            master_shape_all,
                                                            getPalette,
                                                            make_region_cov_maps,
                                                            make_region_cov_hists)
                                     }, mc.cores=num_threads)
    toc()
    
    # Return to multi-threaded processing
    set_original_threads()
  }# end if to allow a group to test a specific cov in single region
  
  ##############################################################################
  ## PULL RASTERS AND REGION NAME FROM REGION OUTPUT LIST
  region_rasters_list <- NULL
  for (reg in regions) {
    region_rasters <- NULL
      for (output in region_output_list) {
        if (reg == as.character(output$region)) {
          ## add cov raster to list of rasters for that region
          region_rasters <- list.append(region_rasters,output[1])
        }# end if statement to find matching region rasters 
      }# end for loop to fo through output
      ## Create list where each index is a list of cpvariate rasters with pixel  
      ## out of range for that region
      assign(paste0(reg,"_rasters"),unlist(region_rasters))
      rasters <- list(c(reg,get(paste0(reg,"_rasters"))))
      region_rasters_list <- c(region_rasters_list,rasters)
  }# end for loop to include all regions
  rasters_meta <- list(start_year, end_year, interval_mo)
  
  ##############################################################################
  ## SET UP ALL REGION PLOTS AND TABLE                                       ##
  ##############################################################################
  
  ##############################################################################
  ## SET UP POPULATION OUT OF RANGE HISTOGRAM
  
  message("\nSetting up population histogram...")
  ## Set serial execution, perform mclapply(), and then go back to parallel 
  ## threads
  if (length(regions)==1){
    pop_output_list <- setup_population_outrange_plot(start_year, 
                                                      end_year, interval_mo, 
                                                      region_rasters_list[[1]])
  } else {
    set_serial_threads()
    
    # Check how many threads we should use, based on a function in
    num_threads <- get_max_forked_threads( nobjs = length(region_effect_list) )
    message(sprintf('We will run an mclapply with %s threads',num_threads))
    pop_output_list <- mclapply(region_rasters_list,
                                function(x){
                                  setup_population_outrange_plot(start_year,
                                                                 end_year,
                                                                 interval_mo,
                                                                 x)
                                }, mc.cores=num_threads)
  
    # Return to multi-threaded processing
    set_original_threads()
  }# end if statement to parallize work if more than one region
  
  ## SEPARATE WORLDPOP AND REGPOP FROM POP OUTPUT LIST #########################
  totalpop <- NULL
  totalpop_outrange <- NULL
  if (length(regions)==1){
    totalpop_outrange <- rbind(totalpop_outrange,pop_output_list[[1]])
    totalpop <- rbind(totalpop,pop_output_list[[2]])
  } else {
    for (i in 1:length(regions)){
      totalpop_outrange <- rbind(totalpop_outrange,pop_output_list[[i]][[1]])
      totalpop <- rbind(totalpop,pop_output_list[[i]][[2]])
    }# end loop to separate pop lists 
  }# end if to handle cases where only 1 region is processed
  
  ## CALC YEAR PERCENTAGE ######################################################
  percentage_pop_out_years <- NULL
  yearly_region_pop_out_years <- NULL
  for (yr in start_year:end_year) {
    ## Find column of year
    index <- match(yr,c(start_year:end_year))
    ## Sum total population out of range for at least 1 covariate in all regions 
    ## and divide by total population in all regions for each year
    percentage_pop_out <- (sum(totalpop_outrange[,index])/
                             sum(totalpop[,index]))*100
    add_frame <- data.frame(year=yr,
                            percentage=percentage_pop_out)
    add_frame_long <- data.frame(year=yr,
                                 totalpop_outrange=totalpop_outrange[,index],
                                 totalpop=totalpop[,index])
    ## Create data frame of year with percentage pop out to plot
    percentage_pop_out_years <- rbind(percentage_pop_out_years,add_frame)
    yearly_region_pop_out_years <- rbind(yearly_region_pop_out_years,add_frame_long)
  }
  
  ##############################################################################
  ## SET UP ALL REGION MAP
  
  message("\nSetting up all region map...")
  if (length(regions)==1){
    map_output_list <- setup_all_region_map(region_rasters_list[[1]])
  } else {
    set_serial_threads()
    
    # Check how many threads we should use, based on a function in
    num_threads <- get_max_forked_threads( nobjs = length(region_effect_list) )
    message(sprintf('We will run an mclapply with %s threads',num_threads))
    map_output_list <- mclapply(region_rasters_list, FUN=setup_all_region_map, 
                                mc.cores=num_threads)
    
    # Return to multi-threaded processing
    set_original_threads()
  }# end if statement to parallize work if more than one region
  
  ## MERGE REGIONS #############################################################
  
  if (length(regions)==1){
    all_regions <- map_output_list
  } else {
    map_output_list$filename <- '<<<< FILEPATH REDACTED >>>>'
    map_output_list$overwrite <- TRUE
    all_regions <- do.call(merge, map_output_list)
  }# end if statement to merge regions if there is more than 1 to map

  ##############################################################################
  ## PULL DATA TABLE FROM REGION OUTPUT LIST
  
  message("\nSetting up summary table...")
  summary_vals <- NULL
  for (output in region_output_list) {
    ## Remove raster and add to list 
    new_summary <- output[-1]
    summary_vals <- rbind(summary_vals,new_summary)
  }# end for each to loop through output and pull data table
  
  ##############################################################################
  ## SAVE PLOTS AND TABLE                                                    ##
  ##############################################################################
  
  ##############################################################################
  ## SAVE POPULATION HISTOGRAM TO OUTDIR AS .PNG
  message("\nSaving histogram...")
  ggplot(percentage_pop_out_years, aes(x=as.factor(year),y=percentage)) +
    geom_col() +
    labs(title="Percentage Population Outside Observed Range for Year",
         y="Percentage Population Outside Observed Range", x="Year") + 
    scale_y_continuous(labels = function(x) paste0(x,"%")) +
    theme_classic()
  
  ggsave(sprintf('<<<< FILEPATH REDACTED >>>>',out_dir))
  
  ## WRITE POP SUMMARY TABLE TO OUTDIR AS .CSV
  message("\nSaving Population Year summary table...")
  write.csv(percentage_pop_out_years, 
            file = sprintf('<<<< FILEPATH REDACTED >>>>',out_dir), 
            row.names = FALSE)
  
  ## WRITE POP SUMMARY TABLE TO OUTDIR AS .CSV
  message("\nSaving Population Year and Region summary table...")
  write.csv(yearly_region_pop_out_years, 
            file = sprintf('<<<< FILEPATH REDACTED >>>>',out_dir), 
            row.names = FALSE)
  
  ##############################################################################
  ## SAVE ALL REGIONS MAP TO OUTDIR AS .PNG
  
  max_val <- max(na.omit(values(all_regions)))
  max_plot <- 10*(max_val%/%10 + as.logical(max_val%%10))
  
  fig <- create_region_map(raster = all_regions,
                           effect = "all",
                           region = "all",
                           getPalette = getPalette,
                           max_plot = max_plot,
                           master_shape_all = master_shape_all,
                           disputed_shp = disputed_shp,
                           buffer = buffer,
                           point_size = point_size,
                           point_stroke_size = point_stroke_size,
                           projection = 'mercator')
  
  ## Hi-res save 
  png(filename=sprintf('<<<< FILEPATH REDACTED >>>>',out_dir), 
      height=5.5, width=12, units='in', res=1000)
  print(fig)
  dev.off()
  
  ##############################################################################
  ## WRITE SUMMARY TABLE TO OUTDIR AS .CSV
  message("\nSaving summary table...")
  write.csv(summary_vals, file = sprintf('<<<< FILEPATH REDACTED >>>>',out_dir), 
            row.names = FALSE)
  
  message("Complete")
}# end main function

################################################################################
## 2: DATA PROCESSING FUNCTIONS                                               ##  
################################################################################

#' Process individual region and effect combinations to create summary information
#' 
#' @description For each region and effect pair it will load in the data to find 
#'   which raster pixels are out of range of the training data and prodcuce summary 
#'   information.
#'   
#' @param input_list Takes a list from check_covariate_ranges that includes the 
#'   region, effect for that region, and measure for the effect to process
#' @param start_year The earliest year of data to include in plots.
#' @param end_year The latest year of data to include in plots.
#' @param interval_mo The number of months in an intervak to load and crop covariate 
#'   raster. Currently load_and_crop_covariates_annual set up to do yearly with 
#'   intervals 12, 24, or 60.
#' @param indicator_group The indicator group associated with this data, used to associate
#'   with a filename in the folder.
#' @param indicator The indicator associated with this data, used to associate
#'   with a filename in the folder.
#' @param run_date The date and time for the training data, used to associate 
#'   with a filename in the folder. 
#' @param out_dir Filepath to save summary plots and tables.
#' @param master_shape_all a combined SPDF with all polygons necessary to make 
#'   maps  
#' @param getPallete univeral color palette for all graphs in 
#'   check_covariate_ranges   
#' @param make_region_cov_maps (default `TRUE`) makes a map showing each the 
#'  number of years a pixel was out of range for each covariate for each 
#'  individual region 
#' @param make_regions_cov_hists (default `TRUE`) makes a histogram showing the 
#'  range of each covariate raster used with the range of values in the training 
#'  data used 
#'      
#' @return returns map of region with number of years out of range for the covariate
#'   directly to out_dir, returns histogram comparing covariate layer and training
#'   data ranges directly to out_dir, returns a list with raster and data.table:
#'     * 'this_pixel-out_range': raster where a pxiel value is 1 if it is out of range.
#'     * 'new_summary': data.table with the summary information for the region
#'         * 'region': region.
#'         * 'covariate': covariate effect.
#'         * 'min_pts': minimum value in the training data.
#'         * 'max_pts': maximum value in the training data.
#'         * 'min_raster': minimum value in the covariate raster.
#'         * 'max_raster': maximum value in the covariate raster. 
#'         * 'pixels_below_range': number of raster values below min_pts.
#'         * 'pixels_above_range': number of raster values above max_pts.
#'         * 'total_pixels_out_range': sum of pixels below and above range.
#'         * 'total_pixels': count of all none na pixels in raster.
#'         * 'percent_out_range': pixels out of range divided by total pixels times 100.
#'         * 'percentage_pop_out_range': sum of people out of range across all years 
#'             divided by sum of people in region across all years times 100.
#'
process_region_effects <- function(input_list,
                                   start_year,
                                   end_year,
                                   interval_mo,
                                   indicator_group,
                                   indicator,
                                   run_date,
                                   out_dir,
                                   master_shape_all,
                                   getPalette,
                                   make_region_cov_maps,
                                   make_region_cov_hists) {
  
  ## parse input list 
  region <- input_list[[1]]
  effect <- input_list[[2]]
  measure <- input_list[[3]]
  
  ## loading in the simple polygon
  adm0_list      <- get_adm0_codes(region, shapefile_version = "current")
  simple_polygon <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, 
                                        tolerance = 0.4,
                                        shapefile_version = "current")
  subset_shape   <- simple_polygon[['subset_shape']]
  simple_polygon <- simple_polygon[['spoly_spdf']]
  
  ## loading in the covariate layers raster 
  cov_layers <- load_and_crop_covariates_annual(covs            = c(effect),
                                                measures        = c(measure),
                                                simple_polygon  = simple_polygon,
                                                start_year      = as.numeric(start_year),
                                                end_year        = as.numeric(end_year),
                                                interval_mo     = as.numeric(interval_mo))
  
  pop_layers <- load_and_crop_covariates_annual(covs            = c("worldpop"),
                                                measures        = c("a0004t"),
                                                simple_polygon  = simple_polygon,
                                                start_year      = as.numeric(start_year),
                                                end_year        = as.numeric(end_year),
                                                interval_mo     = as.numeric(interval_mo))
  
  ## loading in the training data 
  training_data <- fread(sprintf("<<<< FILEPATH REDACTED >>>>",
                                 indicator_group, indicator, run_date, region))
  
  # create background map by subsetting master shape file
  background_map <- master_shape_all[master_shape_all@data$ADM0_CODE %in% 
                                       get_adm0_codes(region,shapefile_version='current'), ]
  
  ## CREATE VARS ###############################################################
  
  ## subset training data to cov
  this_cov_pts <- training_data[[effect]]
  this_cov_raster <- cov_layers[[effect]]
  this_region_pop <- pop_layers[["worldpop"]]
  
  ## remove overlapping borders
  this_cov_raster <- mask(this_cov_raster,background_map)
  this_region_pop <- mask(this_region_pop,background_map)
  
  ## create summary stats across all years
  min_pts <- min(this_cov_pts)
  max_pts <- max(this_cov_pts)
  min_raster <- min(minValue(this_cov_raster))
  max_raster <- max(maxValue(this_cov_raster))
  out_of_range <- sum(na.omit(values(this_cov_raster)) < min_pts |
                        na.omit(values(this_cov_raster)) > max_pts)
  total_pixels <- length(na.omit(values(this_cov_raster)))
  
  ## get the number of years in data
  num_years <- as.numeric(end_year) - as.numeric(start_year) + 1
  
  ## calculate the pixels out of range for the region covariate by setting in 
  ## range pixels to NA then convert to binary so that a pixel is one if out of 
  ## range and 0 if in range
  this_pixel_out_range <- this_cov_raster
  this_pixel_out_range[this_pixel_out_range>=min_pts & 
                         this_pixel_out_range<=max_pts] <- NA
  this_pixel_out_range[!is.na(this_pixel_out_range)] <- 1
  this_pixel_out_range[is.na(this_pixel_out_range)] <- 0
  ## duplicate any single layers that are used all 18 years: will be in range 0 
  ## or 18 times
  if (unlist(dim(this_pixel_out_range))[3] == 1) {
    this_pixel_out_range <- stack(replicate(num_years,this_pixel_out_range))
  }# end if statement to build correct number of yrs cov data
  this_pixel_out_range <- mask(this_pixel_out_range,this_cov_raster)
  
  ## calculate the total pop in the region across all years
  total_pop <- cellStats(this_region_pop, sum)
  total_pop <- Reduce("+", total_pop)
  
  ## calculate the population out of range for the covariate region
  ## multiply pixel out range by worldpop to get total out range for the region 
  ## then divide by total pop for the region across all years
  this_region_pop_out_range <- this_region_pop * this_pixel_out_range
  total_pop_out_range <- cellStats(this_region_pop_out_range, sum)
  total_pop_out_range <- Reduce("+", total_pop_out_range)
  
  ## CREATE OUTPUT TABLE ####################################################### 
  
  ## add summary data to data.table
  new_summary <- data.table(region=region,
                            covariate=effect,
                            min_pts=min_pts,
                            max_pts=max_pts,
                            min_raster=min_raster,
                            max_raster=max_raster,
                            pixels_below_range=sum(na.omit(
                              values(this_cov_raster)) < min_pts),
                            pixels_above_range=sum(na.omit(
                              values(this_cov_raster)) > max_pts),
                            total_pixels_out_range=out_of_range,
                            total_pixels=total_pixels,
                            percentage_out_range=((out_of_range/
                                                     total_pixels)*100),
                            total_pop_out_range=total_pop_out_range,
                            total_pop=total_pop,
                            percentage_pop_out_range=((total_pop_out_range/
                                                         total_pop)*100))
  
  ## PLOT HISTOGRAM OF OBSERVED AND RASTER VALUES ##############################
  
  if (make_region_cov_hists){
    
    ## Create overlapping histogram of point values by density
    training_vals <- data.frame(Value=this_cov_pts, Type="Training Points")
    raster_vals <- melt(na.omit(values(this_cov_raster)))
    raster_vals <- data.frame(Value=raster_vals["value"], Type="Raster Points")
    colnames(raster_vals)[colnames(raster_vals)=='value'] <- 'Value'
    vals_for_plot <- rbind(training_vals,raster_vals)
    
    ggplot(vals_for_plot, aes(x=Value, y = ..density.., fill=Type)) +
      geom_histogram(alpha=0.2, position="identity") +
      scale_fill_manual(values=c("black","red")) + 
      labs(title=sprintf("Point Value Ranges for %s in %s", effect, region)) +
      theme_bw()
    
    ggsave(sprintf('<<<< FILEPATH REDACTED >>>>', out_dir, region, effect))

  }#end make_region_cov_hists
  
  ## CREATE MAP OF REGION WITH NUM TIMES OUT RANGE #############################
  
  if (make_region_cov_maps) {
    ## create map of region with pixels out of range colored
    png(filename=sprintf('<<<< FILEPATH REDACTED >>>>', out_dir, region, effect), 
        units="in", 
        width=20, 
        height=16, 
        pointsize=22, 
        res=72)
    
    ## calc num times the pixel was out of range for this covariate 
    this_pixel_out_range_num <- Reduce("+", unstack(this_pixel_out_range))
    
    ##get max value pixel could be and add 1 for bins 
    max_plot <- num_years + 1
    
    plot(this_pixel_out_range_num+0.01, #add small amount to shift bins
         col=getPalette(max_plot),
         breaks=seq(0,max_plot,n=max_plot), 
         legend=TRUE,
         main = sprintf("Number of Years Pixel was Out of Range for %s in %s", 
                        effect, region)) +
      plot(background_map, add=TRUE)
    
    dev.off()
  }#end make_region_cov_maps
  

  ## RETURN ####################################################################
  output_list <- c(this_pixel_out_range,new_summary)
  return(output_list)
        
}# end process_region_effects


#' Prep output from process_region_effects to plot population histogram
#' 
#' @description Finds the total population out of range and the total population 
#'   living in a single region for each year that is measured. 
#'   
#' @param start_year The earliest year of data to include in plots.
#' @param end_year The latest year of data to include in plots.
#' @param interval_mo The number of months in an intervak to load and crop covariate 
#'   raster. Currently load_and_crop_covariates_annual set up to do yearly with 
#'   intervals 12, 24, or 60.
#' @param reg_rasters Takes a list includes all of the pixel out of range 
#'   rasters, from process_region_effects for a single region and the region 
#'   name as index 1
#'   
#' @return Returns a list of two data tables;
#'  * 'total_region_pop_out_range': total number of people out of range for each year.
#'  * 'total_region_pop': total number of people in the region for each year.
#'  
setup_population_outrange_plot <- function(start_year,
                                           end_year,
                                           interval,
                                           reg_rasters){
   
    region <- reg_rasters[[1]]
    rasters <- reg_rasters[2:length(reg_rasters)]
    
    ## LOAD WORLDPOP LAYER FOR REGION ##########################################
  
    ## Loading in the simple polygon
    adm0_list      <- get_adm0_codes(region, shapefile_version = "current")
    simple_polygon <- load_simple_polygon(gaul_list = adm0_list, buffer = 1, 
                                          tolerance = 0.4,
                                        shapefile_version = "current")
    subset_shape   <- simple_polygon[['subset_shape']]
    simple_polygon <- simple_polygon[['spoly_spdf']]
  
    ## Loading in the covariate layers raster 
    region_pop <- load_and_crop_covariates_annual(covs            = c("worldpop"),
                                                  measures        = c("a0004t"),
                                                  simple_polygon  = simple_polygon,
                                                  start_year      = as.numeric(start_year),
                                                  end_year        = as.numeric(end_year),
                                                  interval_mo     = as.numeric(interval))
  
    ## TABULATE NUM OF PEOPLE OUT OF RANGE THAT YEAR FOR ANY COV ###############
    
    ## create raster stack of layers for each year across all covariates where 0  
    ## means the pixel was not out of range that year for any covariate and 1  
    ## means that the pixel was out of range that year for at least 1 covariate
    ## NOTE: the layer names are defaulted to the first covariate and are no 
    ## longer accurate
    names <- NULL
    num_years <- as.numeric(end_year) - as.numeric(start_year) + 1
    for (i in 1:num_years){
      names <- c(names,paste0("regpop.",i))
    }
    region_pop_out_range <-  Reduce("+", rasters)
    names(region_pop_out_range) <- names
    region_pop_out_range[region_pop_out_range>0] <- 1
    region_pop_out_range <- region_pop_out_range * region_pop[["worldpop"]] 
  
    ## get the total number of people out of range in all pixels for each year 
    ## and the total population in the region for each year
    total_region_pop_out_range <- list(cellStats(region_pop_out_range, sum))
    total_region_pop <- list(cellStats(region_pop[["worldpop"]], sum))
    
    output_list <- c(total_region_pop_out_range,total_region_pop)
    return(output_list)
  
}# end setup_population_outrange_plot function


#' Prep output from process_region_effects for map-making
#' 
#' @description Creates a raster where each pixel value is the percentage of 
#'   covariates that were ever out of range in that region.
#'   
#' @param reg_rasters Takes a list includes all of the pixel out of range 
#'   rasters, from process_region_effects for a single region and the region 
#'   name as index 1.
#'  
#' @return Returns a single raster that contains the percentage of covariates 
#'   out of range in that rgeion.
#'   
setup_all_region_map <- function(reg_rasters){
  
  rasters <- reg_rasters[2:length(reg_rasters)]
  
  ## Create list where each index is a covariate layer for the region where 
  ## 0 means the pixel was never out of range and 1 means the pixel was out 
  ## of range for at least one of the years used 
  ever_out_range <- NULL
  for (layer in rasters){
    layer <- Reduce("+",unstack(layer))
    layer[layer>0] <- 1
    ever_out_range <- list.append(ever_out_range,layer)
  }
  
  ## Create list of rasters for each region where the values of the raster are  
  ## the percentage of covariates that pixel is out of range 
  region_ever_out_range <- Reduce("+", ever_out_range)
  values(region_ever_out_range) <- (values(region_ever_out_range) / 
                                      length(ever_out_range)) * 100
  
  return(region_ever_out_range)
  
}# end setup_all_region_map

#' Create a region map from calculated raster files
#' 
#' @description Takes the raster files computed by other check_covariate_ranges
#'   functions and creates a ggplot map beyond the basic plot functions 
#'   

create_region_map <- function(raster,
                              effect,
                              region,
                              getPalette,
                              max_plot,
                              master_shape_all,
                              disputed_shp,
                              buffer = 10,
                              point_size = 0.08,
                              point_stroke_size = 0.22,
                              projection = 'mercator'){
  
  ## Convert the raster into a data.table for ggplotting
  map_data <- as.data.frame(raster, xy=TRUE, na.rm=TRUE) %>% as.data.table
  names(map_data) <- c('long', 'lat', 'mapvar')
  
  ## set map defaults 
  if (region == "all") {
    col_num <- max_plot -1
    title <- "Percentage of Modelled Covariates Out of Range by Pixel"
  } else {
    col_num = max_plot
    title <- sprintf("Number of Years Pixel was Out of Range for %s in %s", 
                     effect, region)
  }
  
  if (is.null(buffer)) {
    buffer <- 0
  }
  
  ## PLOT GGMAP ################################################################
  
    fig<- ggplot() +
    
    theme_classic(base_size = 12) +
    
    theme(axis.line = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "in"),
          legend.position = "right",
          legend.title = element_text(size=8, color='#222222'),
          legend.text  = element_text(size=6, color='#222222')) +
    
    geom_point(data = map_data, aes(color=mapvar, y = lat, x = long), shape  = 15,
               size   = point_size,
               stroke = point_stroke_size,
               show.legend = T) +
    
    scale_color_gradientn(
      limits = c(0, max_plot),
      colors = c("#E6E6E6",getPalette(col_num)),
      breaks = seq(0, max_plot, 5),
      labels = seq(0, max_plot, 5)) +
    
    geom_path(data = master_shape_all,
      aes(x=long, y=lat, group=group),
      color='#222222',
      lwd=.1) + 
    
    geom_polygon(
      data = disputed_shp,
      aes(x=long, y=lat, group=group),
      fill=NA,
      linetype=3,
      color='#222222',
      lwd=.1) +
    
    labs(title    = title,
      fill     = "Percentage of covariates \nout of range",
      color    = "Percentage of covariates \nout of range",
      x = NULL, y = NULL) + 
    
    xlim(min(map_data$long) - buffer, max(map_data$long) + buffer) +
    
    ylim(min(map_data$lat)  - buffer, max(map_data$lat)  + buffer) +
    
    coord_map(
      projection=projection,
      xlim=c(min(map_data$long) - 0.5, max(map_data$long) + 0.5),
      ylim=c(min(map_data$lat)  - 0.5, max(map_data$lat)  + 0.5)
    ) 
    
  
  return(fig)
  
} # end create region map 

