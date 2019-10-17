
## Make comparison Plots for two run dates
#'
#' @title Make plots comparing two U5M run dates
#' @description Wrapper for `make_comparison_plots`. Loads in admin aggregated data: uses combined csvs if raking_shapefile_version is the same between the two run dates, otherwise uses combined csv for current run date and reaggregates the old_run_date using `custom_aggregate_draws`. 
#'
#' @param current_run_date character string, model run date
#' @param old_run_date character string, model run date to compare with
#' @param current_run_date_name character string, the name of the current dataset to use in the plots
#' @param old_run_date_name character string, the name of the old dataset to use in the plots
#' @param ind mbg indicator 
#' @param ig mbg indicator group
#' @param age string, `under5`, `neonatal`, `infant` acceptable (For U5M use)
#' @param adm_level integer, admin level to compare at, `0`, `1`, or `2`
#' @param out_dir The filepath to save plots to 
#' @param regions A vector of characters, the regions to plot. If `NULL` plots all regions
#' @param plot_adm0 boolean. Plot adm0 boundaries? else plot `adm_level`` boundaries
#' @param cores The number of cores for mclapply to use when aggregating draws by region. Recommended == number of regions.
#' @return NULL
#' 
compare_run_aggregates <- function(current_run_date,
                                   old_run_date,
                                   current_run_date_name,
                                   old_run_date_name,
                                   ind,
                                   ig,
                                   age,
                                   adm_level, 
                                   out_dir,
                                   regions = NULL,
                                   plot_adm0 = T,
                                   cores = get_max_forked_threads(threads = NULL, nobjs = length(regions))) {
  run_date2 <- current_run_date 
  run_date1 <- old_run_date
  
  run_date1_path <- sprintf("<<<< FILEPATH REDACTED >>>>")
  run_date2_path <- sprintf("<<<< FILEPATH REDACTED >>>>")
  
  run_date1_config <- fread(sprintf("<<<< FILEPATH REDACTED >>>>"))
  run_date2_config <- fread(sprintf("<<<< FILEPATH REDACTED >>>>"))
  
  #get year list from configs of run_dates and ensure both are the same
  run_2_year <- run_date2_config[V1 == "year_list",]$V2
  run_1_year <- run_date1_config[V1 == "year_list",]$V2
  
  if(cores == 0){
    #if regions is NULL (pull all regions) and cores are set by default, get number of cores from list.files
    files <- list.files(path = run_date1_path, pattern = "died_under5_cell_draws_eb_bin0_*")
    cores <- get_max_forked_threads(threads = NULL, nobjs = length(files))
  }
  
  #get intersecting years
  yl_intersect <- intersect( eval(parse(text=run_1_year)), eval(parse(text=run_2_year)) )
  if(length(yl_intersect)==0) stop("There are no overlapping years between the two run dates.")
  year_list <- min(yl_intersect):max(yl_intersect)
  
  #get interval mo from configs of run_dates and ensure both are the same
  run_2_int_mo <- run_date2_config[V1 == "interval_mo",]$V2
  run_1_int_mo <- run_date1_config[V1 == "interval_mo",]$V2
  
  if(run_2_int_mo != run_1_int_mo){
    stop("Different interval_mo between run_dates, code is not currently set up to handle properly")
  } else {
    interval_mo <- as.numeric(run_2_int_mo)
  }
  
  #get interval mo from configs of run_dates and ensure both are the same
  run_2_pop_measure <- run_date2_config[V1 == "pop_measure",]$V2
  run_1_pop_measure <- run_date1_config[V1 == "pop_measure",]$V2
  
  if(run_2_pop_measure != run_1_pop_measure){
    stop("Different pop_measure between run_dates, code is not currently set up to handle properly")
  } else {
    pop_measure <- as.character(run_2_pop_measure)
  }
  
  run_2_raking_shapefile <- run_date2_config[V1 == "raking_shapefile_version",]$V2
  run_1_raking_shapefile <- run_date1_config[V1 == "raking_shapefile_version",]$V2
  
  shapefile <- get_admin_shapefile(adm_level, version = run_2_raking_shapefile)
  field <- paste0("ADM", adm_level,"_CODE")
  shapes <- rgdal::readOGR(shapefile)
  shapes[[field]] <- as.integer(as.character(shapes[[field]]))
  
  #TODO: subset to region
  
  if(run_2_raking_shapefile == run_1_raking_shapefile){
    message("raking shapefile versions between run dates are the same, using aggregated draws")
    mean_dt1 <- fread(sprintf("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    lower_dt1 <- fread(sprintf("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    upper_dt1 <- fread(sprintf("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    dt1 <- merge(mean_dt1, lower_dt1, by=c("year", paste0("ADM", adm_level,"_CODE")))
    dt1 <- merge(dt1, upper_dt1, by=c("year", paste0("ADM", adm_level,"_CODE")))
    colnames(dt1) <- c("year", paste0("ADM", adm_level,"_CODE"), "mean", "lower", "upper")
    
    mean_dt2 <- fread("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    lower_dt2 <- fread("<<<< FILEPATH REDACTED >>>>"l), drop = "V1")
    upper_dt2 <- fread("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    dt2 <- merge(mean_dt2, lower_dt2, by=c("year", paste0("ADM", adm_level,"_CODE")))
    dt2 <- merge(dt2, upper_dt2, by=c("year", paste0("ADM", adm_level,"_CODE")))
    colnames(dt2) <- c("year", paste0("ADM", adm_level,"_CODE"), "mean", "lower", "upper")

  } else {
    message("raking shapefile versions between run dates are different, using aggregated draws from current version")
    message("reaggregating old version from draws")
    
    mean_dt2 <- fread("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    lower_dt2 <- fread("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    upper_dt2 <- fread("<<<< FILEPATH REDACTED >>>>"), drop = "V1")
    dt2 <- merge(mean_dt2, lower_dt2, by=c("year", paste0("ADM", adm_level,"_CODE")))
    dt2 <- merge(dt2, upper_dt2, by=c("year", paste0("ADM", adm_level,"_CODE")))
    colnames(dt2) <- c("year", paste0("ADM", adm_level,"_CODE"), "mean", "lower", "upper")
    
    if(is.null(regions)){
      files <- list.files(path = run_date1_path, pattern = "died_under5_cell_draws_eb_bin0_*")
      files <- gsub("died_under5_cell_draws_eb_bin0_", "", files)
      files <- gsub("_0.RDS", "", files)
      custom_regions <- files
    } else {
      custom_regions <- regions
    }
    message("loading and aggregating region draws - this may take awhile depending on the number of cores")
    dt1 <- custom_aggregate_draws(regions = custom_regions,
                                  cores,
                                  age,
                                  run_date = run_date1,
                                  ind,
                                  ig,
                                  shapefile = NULL,
                                  field = NULL,
                                  year_list,
                                  pop_measure,
                                  original_shapefile_version = run_1_raking_shapefile,
                                  standard_shapefile = T,
                                  shapefile_version = run_2_raking_shapefile,
                                  adm_level = adm_level,
                                  aggregation_measure = "mean")
    
    dt1 <- data.table(dt1)
  }
  
  dt1[, mean := mean * 1000]
  dt2[, mean := mean * 1000]
  dt1[, upper := upper * 1000]
  dt2[, upper := upper * 1000]
  dt1[, lower := lower * 1000]
  dt2[, lower := lower * 1000]
  
  #if given regions, subset shapefile and data to desired regions
  if(!is.null(regions)){
    adm0_codes <- get_adm0_codes(regions, shapefile_version = run_2_raking_shapefile)
    shapes <- shapes[shapes$ADM0_CODE %in% adm0_codes,]
    
    adm_codes_to_keep <- data.table(unique(shapes@data[, c("ADM0_CODE", field)]))
    adm_codes_to_keep <- unique(adm_codes_to_keep[ADM0_CODE %in% adm0_codes, field, with =F])
    
    dt1 <- dt1[dt1[[field]] %in% adm_codes_to_keep[[field]]]
    dt2 <- dt2[dt2[[field]] %in% adm_codes_to_keep[[field]]]
  }
  
  if(plot_adm0) {
    adm0_shapefile <- get_admin_shapefile(0, version = run_2_raking_shapefile)
    adm0_shapes <- rgdal::readOGR(adm0_shapefile)
    adm0_shapes$ADM0_CODE <- as.integer(as.character(adm0_shapes$ADM0_CODE))
    if(!is.null(regions)){
      adm0_shapes <- adm0_shapes[adm0_shapes$ADM0_CODE %in% adm0_codes,]
    }
  } else {
    adm0_shapes <- NULL
  }
  
  make_comparison_plots(custom_dt1 = dt1, 
                        custom_dt2 = dt2, 
                        name_1     = old_run_date_name,
                        name_2     = current_run_date_name, 
                        shapefile  = shapes, 
                        field      = field, 
                        out_dir    = out_dir,
                        age        = age,
                        adm0_shapefile = adm0_shapes,
                        cores      = cores)
}


## Make comparison Plots
#'
#' @title Make plots comparing two admin datasets - preparing the data and plotting years in parallel
#' @description Takes two datasets and for each year makes: plot for the values of each dataset, plot of the absolute difference, the relative difference and the CI overlap of the two datasets. The difference plots find custom_dt2 - or / custom_dt1.
#'
#' @param custom_dt1 A data.table with columns `year`, `field` (matching field column below), `mean`, `upper`, `lower`. Can be made with the `custom_aggregate_raster` function
#' @param custom_dt2 A data.table with columns `year`, `field` (matching field column below), `mean`, `upper`, `lower`. Can be made with the `custom_aggregate_raster` function
#' @param name_1 character string, the name of the first dataset to use in the plots
#' @param name_2 character string, the name of the second dataset to use in the plots
#' @param shapefile a filepath to a shapefile or an sp object
#' @param field the field in the shapefile, must match column in `custom_dt1` and `custom_dt2`
#' @param out_dir The filepath to save plots to 
#' @param adm0_shapefile similar to `shapefile`. If null, plots boundaries from shapefile, else plots adm0 boundaries
#' @param cores The number of cores for mclapply to use when plotting (split over years)
#'
#' @return NULL
#'
make_comparison_plots <- function(custom_dt1,
                                  custom_dt2,
                                  name_1,
                                  name_2,
                                  shapefile,
                                  field, 
                                  out_dir,
                                  age,
                                  adm0_shapefile = NULL,
                                  cores = 10){
  
  #rename to avoid renaming out of scope
  c_dt1 <- custom_dt1
  c_dt2 <- custom_dt2
  
  #get shapefile - can take a file path or a spdf object
  if(is.character(shapefile)){
    message("Loading shapefile")
    shapes <- rgdal::readOGR(shapefile)
    shapes[[field]] <- as.integer(as.character(shapes[[field]]))
  } else {
    shapes <- shapefile
  }
  
  if(!is.null(adm0_shapefile)){
    if(is.character(adm0_shapefile)){
      message("Loading adm0 shapefile")
      adm0_shapes <- rgdal::readOGR(shapefile)
      adm0_shapes$ADM0_CODE <- as.integer(as.character(shapes$ADM0_CODE))
    } else {
      adm0_shapes <- adm0_shapefile
    }
    adm0_to_plot <- fortify(adm0_shapes)
  } else {
    adm0_shapes <- NULL
  }
  
  setnames(c_dt1, c("mean", "lower", "upper"), c("mean1", "lower1", "upper1"))
  setnames(c_dt2, c("mean", "lower", "upper"), c("mean2", "lower2", "upper2"))
  
  combined_dt <- merge(c_dt1, c_dt2, by=c(field, "year"))
  
  #plotting variables
  combined_dt[, abs:= mean2 - mean1]
  combined_dt[, rel:= mean2 / mean1]
  # Check for non-overlapping CIs
  combined_dt[, c('two_high','two_low','no_overlap') := 0]
  combined_dt[!is.na(upper1) & !is.na(lower2) & (upper1 < lower2), 
              two_high := 1]
  combined_dt[!is.na(lower1) & !is.na(upper2) & (lower1 > upper2),
              two_low := 1]
  # Too high --> 1, Too low --> -1
  combined_dt[, no_overlap := 0 + two_high - two_low]
  combined_dt[ is.na(mean1) | is.na(mean2), no_overlap := NA]
  
  ## Set color schemes here
  na_color <- '#E6E6E6'
  col_breaks_abs <- c(0, 5, 25, 50, 100, 200)
  col_labs_abs <- c("", "5", "25", "50", "100", ">200")
  grad_vals_abs <- c(
    '#009999','#5aada1','#88c1a9','#b1d6b0','#d8eab8','#ffffbf',"#fff2b6",
    "#fee5ae","#fed8a5","#fecb9c","#fdbe93","#fdb18b","#fda482","#fc9779",
    "#fc8a70","#fc7d68","#fb705f","#fb6356","#fb564d","#fa4945","#fa3c3c",
    "#f53b3a","#f03938","#eb3837","#e63635","#e13533","#db3431","#d6322f",
    "#d1312e","#cc2f2c","#c72e2a","#c22d28","#bd2b26","#b82a25","#b32823",
    "#ae2721","#a8261f","#a3241d","#9e231c","#99211a","#942018"
  )
  #year_list <- c(2012:2017)
  mclapply(year_list, function(i) {
    year = i
    make_comparison_plots_year(year,
                               combined_dt,
                               name_1,
                               name_2,
                               shapes,
                               field,
                               age,
                               na_color,
                               col_breaks_abs,
                               col_labs_abs,
                               grad_vals_abs,
                               out_dir,
                               adm0_shapes,
                               adm0_to_plot)
  }, mc.cores = cores)
}


## Make comparison plots year
#'
#' @title Make plots comparing two admin datasets - plotting 
#' @description Subfunction of `make_comparison_plots`. Called in mclapply looping over years. Not for separate use.
#'
#' @param this_year The year to plot. Comes from mclapply.
#' @param combined_dt A data.table with columns `year`, `field` (matching field column below), `mean`, `upper`, `lower`. Comes from `make_comparison_plots`
#' @param name_1 character string, the name of the first dataset to use in the plots
#' @param name_2 character string, the name of the second dataset to use in the plots
#' @param shapes an sp object to plot
#' @param field the field in `shapes`, must match column in `combined_dt1`
#' @param na_color NA color for plots. Comes from `make_comparison_plots`
#' @param col_breaks_abs legend breaks for run date plot. Comes from `make_comparison_plots`
#' @param col_labs_abs labels for run date plot. Comes from `make_comparison_plots`
#' @param grad_vals_abs color values for run date plot. Comes from `make_comparison_plots`
#' @param out_dir The filepath to save plots to 
#' @param adm0_shapefile similar to `shapefile`. If null, plots boundaries from shapefile, else plots adm0 boundaries
#'
#' @return NULL
#'
make_comparison_plots_year <- function(this_year,
                                       combined_dt,
                                       name_1,
                                       name_2,
                                       shapes,
                                       field,
                                       age,
                                       na_color,
                                       col_breaks_abs,
                                       col_labs_abs,
                                       grad_vals_abs,
                                       out_dir,
                                       adm0_shapes,
                                       adm0_to_plot) {
  message("Making plots for year: ", this_year)
  # Check overlap specifically for this year
  this_year_data <- combined_dt[year==this_year,]
  
  # Create mapping dataset
  mapping_data <- suppressMessages(
    prep_shp_data_for_mapping(shp=shapes, dataset=this_year_data, merge_var=field))
  
  S2_BOUNDS <- list(
    'lat_min'  = -57,
    'lat_max'  = 63,
    'long_min' = -129,
    'long_max' = 165
  )
  
  # Subset to the Stage 2 boundaries
  mapping_data[lat  < S2_BOUNDS$lat_min,  lat  := S2_BOUNDS$lat_min  ]
  mapping_data[lat  > S2_BOUNDS$lat_max,  lat  := S2_BOUNDS$lat_max  ]
  mapping_data[long < S2_BOUNDS$long_min, long := S2_BOUNDS$long_min ]
  mapping_data[long > S2_BOUNDS$long_max, long := S2_BOUNDS$long_max ]
  if(!is.null(adm0_shapes)){
    adm0_to_plot <- as.data.table(adm0_to_plot)
    adm0_to_plot[lat  < S2_BOUNDS$lat_min,  lat  := S2_BOUNDS$lat_min  ]
    adm0_to_plot[lat  > S2_BOUNDS$lat_max,  lat  := S2_BOUNDS$lat_max  ]
    adm0_to_plot[long < S2_BOUNDS$long_min, long := S2_BOUNDS$long_min ]
    adm0_to_plot[long > S2_BOUNDS$long_max, long := S2_BOUNDS$long_max ]
  }
  
  if(age == "under5"){
    age_name <- "5q0"
  } else if(age == "infant"){
    age_name <- "1q0"
  } else {
    age_name <- "NN mortality"
  }
  
  scatterplot <- ggplot(data = this_year_data) +
    geom_point(aes(x=mean1, y=mean2, alpha = .5)) +
    geom_abline(slope=1, intercept=0) +
    labs(title=this_year) +
    xlab(paste0(name_1, " ", age_name, " (per 1000 live births)")) +
    ylab(paste0(name_2, " ", age_name, " (per 1000 live births)")) +
    theme_minimal() +
    theme(legend.position = "none")
  
  png(
    sprintf("<<<< FILEPATH REDACTED >>>>"), 
    height=5.5, width=12, units='in', res=800
  )
  print(scatterplot)
  dev.off()
  
  if(this_year %in% c(2000, 2010, 2017)){
    m <- lm(mean2 ~ mean1, this_year_data)
    
    if(!file.exists("<<<< FILEPATH REDACTED >>>>"))) file.create("<<<< FILEPATH REDACTED >>>>")
    
    sink(sprintf("<<<< FILEPATH REDACTED >>>>"), append = T)
    cat(this_year, ": ", summary(m)$r.squared, "\n")
    sink()
  }
  
  
  
  # Plot values of dataset 1 directly
  name_1_plot <- ggplot() +
    geom_polygon(data = mapping_data,
                 aes(x=long, y=lat, group=group, fill=mean1)) +   
    scale_fill_gradientn(limits = c(0, max(combined_dt$mean1, na.rm =T)),
                         colors = grad_vals_abs,
                         breaks = col_breaks_abs,
                         labels = col_labs_abs,
                         na.value = na_color,
                         name = sprintf('%s per 1000\nLive Births', age_name)) +
    labs(title = sprintf('Local Burden of Disease %s: %s', name_1, this_year)) +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data, projection = "mollweide")
  
  if(!is.null(adm0_shapes)){
    name_1_plot <- name_1_plot +
      geom_path(data = adm0_to_plot,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  } else {
    name_1_plot <- name_1_plot +
      geom_path(data = mapping_data,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  }
  
  png(
    sprintf("<<<< FILEPATH REDACTED >>>>"), 
    height=5.5, width=12, units='in', res=800
  )
  print(name_1_plot)
  dev.off()
  
  
  # Plot values of dataset 2 directly
  name_2_plot <- ggplot() +
    geom_polygon(data = mapping_data,
                 aes(x=long, y=lat, group=group, fill=mean2)) + 
    scale_fill_gradientn(limits = c(0, max(combined_dt$mean2, na.rm =T)),
                         colors = grad_vals_abs,
                         breaks = col_breaks_abs,
                         labels = col_labs_abs,
                         na.value = na_color,
                         name = sprintf('%s per 1000\nLive Births', age_name)) +
    labs(title = sprintf('Local Burden of Disease %s: %s', name_2, this_year)) +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data, projection = "mollweide")
  
  if(!is.null(adm0_shapes)){
    name_2_plot <- name_2_plot +
      geom_path(data = adm0_to_plot,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  } else {
    name_2_plot <- name_2_plot +
      geom_path(data = mapping_data,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  }
  
  png(
    sprintf("<<<< FILEPATH REDACTED >>>>"), 
    height=5.5, width=12, units='in', res=800
  )
  print(name_2_plot)
  dev.off()
  
  
  # Plot absolute difference of dataset 2 - dataset 1
  abs_min <- min(this_year_data$abs, na.rm =T)
  abs_max <- max(this_year_data$abs, na.rm =T)
  abs_val <- max(abs(c(abs_min, abs_max)))
  fig_adiff_this_year <- ggplot() +
    geom_polygon(data = mapping_data, aes(x=long, y=lat, group=group, fill=abs)) + 
    scale_fill_gradient2(limits = c(-abs_val, abs_val),
                         low = '#b35806',
                         high = '#542788',
                         mid = "white",
                         breaks = c(-abs_val, 0, abs_val),
                         labels = as.character(c(round(-abs_val + .05,1), 0, round(abs_val-.05,1))),
                         midpoint = 0,
                         na.value = na_color) +
    labs(title = sprintf('Absolute Difference between %s and %s: %s',name_1, name_2, this_year),
         fill = sprintf('Difference in\n%s per 1000\nLive Births\n(%s - \n%s)',age_name, name_2, name_1)) +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data, projection = "mollweide")
  
  if(!is.null(adm0_shapes)){
    fig_adiff_this_year <- fig_adiff_this_year +
      geom_path(data = adm0_to_plot,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  } else {
    fig_adiff_this_year <- fig_adiff_this_year +
      geom_path(data = mapping_data,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  }
  
  png("<<<< FILEPATH REDACTED >>>>", 
      height=5.5, width=12, units='in', res=800)
  print(fig_adiff_this_year)
  dev.off()
  
  
  # Plot absolute difference of dataset 2 / dataset 1
  rel_min <- min(this_year_data$rel, na.rm = T)
  rel_max <- max(this_year_data$rel, na.rm = T)
  rel_val <- max(abs(c(rel_min, rel_max)))
  fig_reldiff_this_year <- ggplot() +
    geom_polygon(data = mapping_data, aes(x=long, y=lat, group=group, fill=rel)) + 
    scale_fill_gradient2(limits = c(1/(rel_val*1.1), (rel_val*1.1)),
                         low = '#b35806',
                         high = '#542788',
                         mid = "white",
                         breaks = c(1/rel_val, 1, rel_val),
                         labels = as.character(c(round((1/rel_val) + .005,2), 1, round((rel_val) -.005,2))),
                         midpoint = 1,
                         na.value = na_color) +
    labs(title = sprintf('Ratio between %s and %s: %s',name_1, name_2, this_year),
         fill = sprintf('Difference in\n%s per 1000\nLive Births\n(%s / \n%s)', age_name, name_2, name_1)) +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data, projection = "mollweide")
  
  if(!is.null(adm0_shapes)){
    fig_reldiff_this_year <- fig_reldiff_this_year +
      geom_path(data = adm0_to_plot,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  } else {
    fig_reldiff_this_year <- fig_reldiff_this_year +
      geom_path(data = mapping_data,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  }
  
  png("<<<< FILEPATH REDACTED >>>>", 
      height=5.5, width=12, units='in', res=800)
  print(fig_reldiff_this_year)
  dev.off()
  
  
  # Plot confidence interval overlap of dataset 2 and dataset 1
  mapping_data[, no_overlap := as.character(no_overlap)]
  unique_vals <- unique(mapping_data$no_overlap)
  CI_values = c('-1'='#0066cc','0'='#ffffe6','1'='#ff0066', 'NA'="grey50")
  CI_labels = c('-1' = sprintf("%s < %s", name_2, name_1),
                '0' = "Overlap",
                '1' = sprintf("%s > %s", name_2, name_1), 
                'NA' = "No estimates")
  CI_values = CI_values[names(CI_values) %in% unique_vals]
  CI_labels = CI_labels[names(CI_labels) %in% unique_vals]
  
  fig_uis_this_year <- ggplot() +
    geom_polygon(data = mapping_data,
                 aes(x=long, y=lat, group=group, fill=no_overlap)) + 
    scale_fill_manual(values = CI_values,
                      labels = CI_labels,
                      na.value = "grey50") +
    labs(title = sprintf('Areas with non-overlapping UIs between %s and %s: %s',name_1, name_2,this_year),
         fill = 'Difference\nbetween UIs') +
    theme_map() +
    coord_map_to_bounds(shp_fort=mapping_data, projection = "mollweide")
  
  if(!is.null(adm0_shapes)){
    fig_uis_this_year <- fig_uis_this_year +
      geom_path(data = adm0_to_plot,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  } else {
    fig_uis_this_year <- fig_uis_this_year +
      geom_path(data = mapping_data,
                aes(x=long, y=lat, group=group),
                color='#404040',
                lwd=0.1)
  }
  
  png(
    sprintf("<<<< FILEPATH REDACTED >>>>"), 
    height=5.5, width=12, units='in', res=800
  )
  print(fig_uis_this_year)
  dev.off()
}

## Aggregate to custom shapefile
#'
#' @title Aggregate draws to custom shapefile wrapper
#' @description Wrapper for `custom_aggregate_draws_reg` which makes region aggregate draws in parallel
#'
#' @param regions The list of regions to aggregate to
#' @param cores The number of cores for mclapply to use. Recommended == number of regions.
#' @param age string, `under5`, `neonatal`, `infant` acceptable (For U5M use)
#' @param run_date character string, model run date
#' @param ind mbg indicator 
#' @param ig mbg indicator group
#' @param shapefile a filepath to a shapefile or an sp object
#' @param field the field in the shapefile to aggregate to
#' @param year_list numeric list of years present in raster bricks
#' @param pop_measure the population measure to pull in `load_and_crop_covariates_annual()`
#' @param original_shapefile_version standard lbd shapefile versions, the shapefile version used in modelling the cell draws.
#' @param standard_shapefile boolean, if T, builds simple raster, else makes custom admin raster. Must have `shapefile_version` if `T`
#' @param shapefile_version standard lbd shapefile versions - NULL by default, only used if `standard_shapefile` is `T`
#' @param adm_level only used if `standard_shapefile` is `T`. The adm level to aggregate to.
#' @param aggregation_measure type of aggregation. `mean` and `sum` are acceptable
#'
#' @return A data.table with mean, upper, and lower values aggregated to year and field for all `regions`
#'
custom_aggregate_draws <- function(regions,
                                   cores,
                                   age,
                                   run_date,
                                   ind,
                                   ig,
                                   shapefile,
                                   field,
                                   year_list,
                                   pop_measure,
                                   original_shapefile_version,
                                   standard_shapefile = F,
                                   shapefile_version = NULL,
                                   adm_level = 2,
                                   aggregation_measure = "mean") {
  
  agg_list <- mclapply(regions, function(i) {
    reg = i
    custom_aggregate_draws_reg(reg = reg, 
                               age = age,
                               run_date = run_date,
                               ind = ind,
                               ig = ig,
                               shapefile = shapefile,
                               field = field,
                               year_list = year_list,
                               pop_measure = pop_measure,
                               original_shapefile_version = original_shapefile_version,
                               standard_shapefile = standard_shapefile,
                               shapefile_version = shapefile_version,
                               adm_level = adm_level,
                               aggregation_measure = aggregation_measure)
    }, mc.cores = cores)
  
  combined_agg <- do.call(rbind, agg_list)
  
  return(combined_agg)
}

## Aggregate to custom shapefile
#'
#' @title Aggregate draws to custom shapefile units
#' @description Aggregates region draws to shapefile boundaries
#'
#' @details takes in a shapefile and shapefile field to aggregate to. Makes a custom admin raster, then pulls population for that raster. Takes the cell draws, along with population and admin code, and combines into a table. Aggregates to shapefile field. Can take mean for rates or sum for counts.
#'
#' @param reg The region to aggregate to
#' @param age string, `under5`, `neonatal`, `infant` acceptable (For U5M use)
#' @param run_date character string, model run date
#' @param ind mbg indicator 
#' @param ig mbg indicator group
#' @param shapefile a filepath to a shapefile or an sp object
#' @param field the field in the shapefile to aggregate to
#' @param year_list numeric list of years present in raster bricks
#' @param pop_measure the population measure to pull in `load_and_crop_covariates_annual()`
#' @param original_shapefile_version standard lbd shapefile versions, the shapefile version used in modelling the cell draws.
#' @param standard_shapefile boolean, if T, builds simple raster, else makes custom admin raster. Must have `shapefile_version` if `T`
#' @param shapefile_version standard lbd shapefile versions - NULL by default, only used if `standard_shapefile` is `T`. The shapefile version to reaggregate to
#' @param aggregation_measure type of aggregation. `mean` and `sum` are acceptable
#'
#' @return A data.table with mean, upper, and lower values aggregated to year and field
#'
custom_aggregate_draws_reg <- function(reg,
                                       age,
                                       run_date,
                                       ind,
                                       ig,
                                       shapefile,
                                       field,
                                       year_list,
                                       pop_measure,
                                       original_shapefile_version,
                                       standard_shapefile = F,
                                       shapefile_version = NULL,
                                       adm_level = 2,
                                       aggregation_measure = "mean") {
  
  cell_draws <- readRDS("<<<< FILEPATH REDACTED >>>>")
  
  if(standard_shapefile){
    field <- paste0("ADM", adm_level, "_CODE")
  }
  
  #find/calculate interval month
  if(length(year_list) == 1){
    interval_mo <- 12
  } else {
    year_diff <- diff(year_list)
    if(length(unique(year_diff))!=1) {
      stop("Please use annual or 5-year intervals exclusively in year_list")
    } else {
      interval_mo <- year_diff[[1]] * 12
    }
  }
  
  #raster to check against admin raster
  mean_raster <- brick("<<<< FILEPATH REDACTED >>>>")
  mean_raster <- mean_raster[[1]]
  
  message("Loading admin raster\n")
  if(standard_shapefile){
    message("Using lbd standard shapefile")
    if(!is.null(shapefile_version)){
      message('Loading simple polygon')
      simple_polygon <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = shapefile_version),
                                            buffer = 0.4,
                                            shapefile_version = shapefile_version)
      subset_shape   <- simple_polygon[['subset_shape']]
      simple_polygon <- simple_polygon[['spoly_spdf']]
      
      message('Loading simple raster\n')
      raster_list    <- build_simple_raster_pop(subset_shape)
      admin_raster  <- raster_list[['simple_raster']]
    } else {
      stop("Standard_shapefile is T, shapefile must be a standard shapefile date, i.e. "<<<< FILEPATH REDACTED >>>>"")
    }
  } else {
    #get shapefile - can take a file path or a spdf object
    if(is.character(shapefile)){
      message("Loading shapefile")
      shapes <- rgdal::readOGR(shapefile)
    } else {
      shapes <- shapefile
    }
    
    #turns the shapefile into a raster with extent matching mean_raster and values = field
    admin_raster <- load_custom_admin_raster(NULL, field, mean_raster, shapefile=shapes)
  }
  
  #check that the admin raster matches the mean raster (problematic if admin_raster does not match up with
  #cell pred, i.e. different shapefile version than modelling version
  admin_length <- length(seegSDM:::notMissingIdx(admin_raster))
  mean_length <- length(seegSDM:::notMissingIdx(mean_raster[[1]]))
  if((length(mean_length) * length(year_list)) != nrow(cell_draws)){
    message('Recreating simple_raster')
    message("mean raster does not match pixels in cell draws")
    mbg_dir <- "<<<< FILEPATH REDACTED >>>>"
    simple_raster_fp <- sprintf("<<<< FILEPATH REDACTED >>>>")
    if(file.exists(simple_raster_fp)){
      message('Loading simple_raster from file')
      ## If a template file already exists, use that
      attach(simple_raster_fp) # Create a temporary environment from this Rdata object
      mean_raster <- simple_raster # Add an object from the temporary env to your working env
      detach() # Remove the temporary env
    } else {  
      message("loading original simple raster as reference raster for cell draws")
      simple_polygon <- load_simple_polygon(gaul_list = get_adm0_codes(reg, shapefile_version = original_shapefile_version),
                                            buffer = 0.4,
                                            shapefile_version = original_shapefile_version)
      subset_shape   <- simple_polygon[['subset_shape']]
      simple_polygon <- simple_polygon[['spoly_spdf']]
      
      message('Loading simple raster\n')
      raster_list    <- build_simple_raster_pop(subset_shape)
      mean_raster  <- raster_list[['simple_raster']]
    }
    mean_length <- length(seegSDM:::notMissingIdx(mean_raster[[1]]))
  }
  
  if(admin_length != mean_length){
    message("crosswalking cell pred to new admin raster")
    cell_draws <- crosswalk_cell_pred_add_NA(mean_raster, admin_raster, cell_draws, year_list = year_list)
  }

  ## Pull annual population brick using new covariates function
  pop_raster_annual <- load_and_crop_covariates_annual(covs           = 'worldpop_raked',
                                                       measures       = pop_measure,
                                                       simple_polygon = admin_raster,
                                                       start_year     = min(year_list),
                                                       end_year       = max(year_list),
                                                       interval_mo    = as.numeric(interval_mo),
                                                       agebin=1)[[1]]
  
  pop_raster_annual  <- crop(pop_raster_annual, extent(admin_raster))
  pop_raster_annual  <- setExtent(pop_raster_annual, admin_raster)
  pop_raster_annual  <- mask(pop_raster_annual, admin_raster)
  
  message("\nConstructing aggregation table")
  #turn population and pixel_id into a data.table
  pixel_id <- seegSDM:::notMissingIdx(admin_raster)
  pixel_spatial<-data.table(pixel_id=pixel_id)
  pop <- data.table(raster::extract(pop_raster_annual, pixel_id)) 
  pop[,pixel_id:=pixel_id]
  pop<-melt(pop,id.vars="pixel_id")
  pop[, year := (min(year_list) - 1) + as.numeric(gsub("worldpop_raked.", "", variable))]
  pop<-pop[,list(pixel_id,year,pop=value)]
  # Setting values where pop is NA or 0 to 0.01 to avoid NAs and NaNs in aggregation
  pop[is.na(pop),pop:=0.01]
  pop[pop == 0 ,pop:=0.01]
  
  if(standard_shapefile) {
    region_adm0_list<-get_adm0_codes(reg, shapefile_version = shapefile_version)
    
    admin_levels<-list() # Emtpy list of levels that will be filled with admin levels
    for(lvl in adm_level){
      fieldname<-paste0("ADM",lvl,"_CODE")
      modeling_shapefile_version <- shapefile_version
      admin_info<-GetAdmin(admin_level=lvl,admin_raster,region_adm0_list,shapefile_version=modeling_shapefile_version)
      pixel_spatial[[fieldname]]<-raster::extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the pixel_spatial data.table
      admin_levels[[as.character(lvl)]]<-admin_info # Add the admin info to the list
      if(sum(is.na(pixel_spatial[[fieldname]]))>0){ # Check to see if any of the pixels don't have a location assigned
        message(paste0("   Whoah, there are some pixels that are NA, and have not been assigned a location for level ",lvl))
      }
    }
    
    pop<-merge(pop,pixel_spatial,by="pixel_id",all.x=T) # Merging on the spatial information to the population information.
    pop<-pop[order(year,pixel_id)]
  } else {
    #converting admin_raster into vector form
    admin_extract <- raster::extract(admin_raster, extent(value_raster))
    admin_extract  <- na.omit(admin_extract)
    
    #adding on admin codes to pop table
    if((length(admin_extract) * length(year_list)) != nrow(pop)){
      stop("the number of non-na values in the admin raster does not match the number of values in the pop raster")
    } else {
      pop <- cbind(pop, admin_extract)
    }
  }
  
  #adding on draws values to pop table
  if(nrow(cell_draws) != nrow(pop)){
    stop("the number of non-na values in the admin raster does not match the number of values in the cell_draws")
  } else {
    pop <- cbind(cell_draws, pop)
  }
  
  overs <- paste0("V", 1:ncol(cell_draws))
  
  message("Aggregating")
  #aggregate to year, field
  if(aggregation_measure == "mean"){
    admin_cell_pred <- pop[, lapply(overs, function(x) weighted.mean(get(x), pop, na.rm = T)), by = c('year', field)]
  } else if(aggregation_measure == "sum"){
    admin_cell_pred <- pop[, lapply(overs, function(x) sum(get(x) * pop, na.rm = T)), by = c('year', field)]
  } else{
    stop("Aggregation measure must be either mean or sum")
  }
  
  #aggregate draws
  admin_cell_pred <- admin_cell_pred[, mean := rowMeans(.SD, na.rm = T), .SDcols = overs]
  
  #convert draws to matrix to calculate row quantiles faster
  admin_matrix <- as.matrix(admin_cell_pred[, overs, with=F])
  upper <- rowQuantiles(admin_matrix, probs = 97.5 / 100)
  lower <- rowQuantiles(admin_matrix, probs = 2.5 / 100)
  
  agg_pop <- cbind(admin_cell_pred, upper, lower)
  
  setnames(agg_pop, field, "admin_extract")
  agg_pop <- agg_pop[, c("year", "admin_extract", "mean", "upper", "lower")]
  setnames(agg_pop, "admin_extract", field)
  
  return(agg_pop)
}

GetAdmin<-function(admin_level,simple_raster, region_adm0_list, shapefile_version){
  message(paste0("Loading admin level ",admin_level))
  
  # load admin shape file
  admin_shp <- rgdal::readOGR(dsn=get_admin_shapefile(admin_level, version = shapefile_version))
  
  # ensure that the rasterize variable is a numeric
  admin_shp@data[[paste0('ADM', admin_level, '_CODE')]] <- as.numeric(as.character(admin_shp@data[[paste0('ADM', admin_level, '_CODE')]]))
  
  # if it doesn't exist, get areas of polygons. 
  if(is.null(admin_shp$Shape_Area)){
    admin_shp$Shape_Area <- area(admin_shp) / 1e6 ## TODO
  }
  
  message("Rasterizing...")
  # we order by area so small places don't get buried under big places (e.g. Lesotho and S. Africa)
  admin_rast<-rasterize(admin_shp[order(admin_shp$Shape_Area),],simple_raster,paste0("ADM",admin_level,"_CODE"), fun="first")
  
  message("Converted to raster based on simple_raster template. Cropping and masking:")
  admin_rast  <- crop(admin_rast, extent(simple_raster))
  admin_rast  <- setExtent(admin_rast, simple_raster)
  admin_rast  <- mask(admin_rast, simple_raster)
  
  message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
  admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% region_adm0_list,]
  admin_centroids<-SpatialPointsDataFrame(gCentroid(admin_shp, byid=TRUE), admin_shp@data, match.ID=FALSE)
  
  message("Compiling and returning results.")
  admin<-list()
  admin[["spdf"]]<-admin_shp
  admin[["centroids"]]<-admin_centroids
  admin[["rast"]]<-admin_rast
  admin[["attributes"]]<-copy(data.table(admin_shp@data))
  
  return(admin)
}

