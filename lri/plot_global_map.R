####################################################################################################
## Description:   Make maps of estimates in Africa for specified years and year pairs.
##
## Inputs:        mean rasters
##                adm0 and adm1 estimates
##                country outlines, lakes, and population mask
##                optional: custom_data, 
##                          which is a list of aggregated data tables named admin0, admin1, admin2
##                          and requires the following columns: ADM0_CODE, ADM1_CODE, ADM2_CODE, 
##                          and whatever "type" columns your plotting such as "mean", "upper", "lower"
##
## Output:        PDF of maps (/share/geospatial/mbg/[indicator_group]/[indicator]/output/
##                  [run_date]/results_maps/[indicator]_raked_mean.pdf').
####################################################################################################


## make sure required packages are loaded --------------------------
library('viridis')
library('data.table')


## load_map_annotations ----------------------------------------------------------------------------

load_map_annotations <- function(map_region) {
  
  ## Base shapefile (country outlines)
  if (map_region == 'global') {
    message('->loading global country borders')
    stage1 <- shapefile('<<<< FILEPATH REDACTED >>>>')
    stage2 <- shapefile('<<<< FILEPATH REDACTED >>>>')
    adm0 <- bind(stage1, stage2) %>% fortify
  } else {
    message('->loading specific country borders')
    gaul_list           <- get_adm0_codes(map_region, shapefile_version = modeling_shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 0.02, tolerance = 0.01, shapefile_version = modeling_shapefile_version)
    simple_polygon      <- simple_polygon_list[[2]]
    adm0 <- simple_polygon %>% fortify
  }
  
  ## Lakes
  message('-->loading lake mask')
  lakes <- raster('<<<< FILEPATH REDACTED >>>>') %>% 
    rasterToPoints %>% 
    as.data.table %>% 
    setnames(., c("long", 'lat', 'lakes'))
  
  ## Population mask
  message('--->loading population mask')
  mask <- raster('<<<< FILEPATH REDACTED >>>>') %>% 
    rasterToPoints %>% 
    as.data.table %>% 
    setnames(., c("long", 'lat', 'mask'))
  
  ## Stage 3 mask
  message("---->loading stage 3 mask")
  stage3 <- shapefile('<<<< FILEPATH REDACTED >>>>') %>% fortify
  
  return(list(adm0 = adm0, lakes = lakes, mask = mask, stage3 = stage3))
  
}


## load_map_results --------------------------------------------------------------------------------

load_map_results <- function(indicator, indicator_group, run_date, map_region, type, raked, start_year, end_year, 
                             single_year=0, geo_levels = c("raster", "admin1", "admin2"),
                             custom_data, raked_measure = NULL, use.sf = F,
                             cores=1) {
  
  if (missing(custom_data)) message('You have not supplied custom data. It will be loaded for you.')
  
  years <- paste(start_year, end_year, sep="_")
  
  if (single_year == 0) year_list <- start_year:end_year else year_list <- single_year
  message(year_list)
  
  ## Set the input directory
  maps_path <- '<<<< FILEPATH REDACTED >>>>'
  
  ## raster estimates
  if ("raster" %in% geo_levels) {
    message('loading raster data')
    if(missing(custom_data)) {
      raster <- paste0(maps_path, "/", indicator, "_", 
                       ifelse(type == "cirange", "range", type),
                       "_", 
                       ifelse(raked, "raked_", ""), 
                       years, ".tif") %>% brick
    } else raster <- brick(custom_data$raster)
    
    message('-> raster found and bricked')
    
    raster <- mclapply(year_list, 
                       function(y) {
                         
                         message('--> sending to points (year=', y, ')')
                         df <- rasterToPoints(raster[[y - 1999]]) %>% data.table
                         setnames(df, c("long", 'lat', 'outcome'))
                         df[, year := y]
                         
                         return(df)
                         
                       },
                       mc.cores=cores) %>% 
      rbindlist %>% 
      setkey(., long, lat, year)
    
    message('--> converted to dt and keyed')
    
  }
  
  ## admin0 estimates and shape file
  if ("admin0" %in% geo_levels) {
    message('loading admin0 data')
    
    if(missing(custom_data)) {
      pred <- '<<<< FILEPATH REDACTED >>>>' %>% fread
    } else {
      message('using custom data')
      pred <- custom_data$admin0
    }
    message('-> admin0 found and fread')
    
    pred <- pred[year %in% year_list, c('ADM0_CODE', 'year', type), with = F]
    setnames(pred, type, 'outcome')
    
    if (map_region != 'global') {
      cty_code <- get_adm0_codes(map_region, shapefile_version = modeling_shapefile_version)
      pred <- pred[ADM0_CODE == cty_code]
    }
    
    if(use.sf==T) {
      
      admin0 <- get_admin_shapefile(0) %>% st_read
      admin0 <- admin0[admin0$ADM0_CODE %in% pred$ADM0_CODE,]
      admin0 <-merge(admin0, pred, by="ADM0_CODE", allow.cartesian=T)
      
      message('--> admin0 results merged to sf')
      
    } else {
      
      admin0 <- shapefile(get_admin_shapefile(0))
      admin0 <- admin0[admin0@data$ADM0_CODE %in% pred$ADM0_CODE,]
      admin0 <- SpatialPolygonsDataFrame(gSimplify(admin0, tol = 0.1, topologyPreserve = T), data = admin0@data)
      for (i in 1:length(admin0)) admin0@polygons[[i]]@ID <- as.character(admin0@data[i, "ADM0_CODE"])
      admin0 <- data.table(fortify(admin0))
      admin0[, ADM0_CODE := as.numeric(id)]
      admin0 <- merge(admin0, pred, by="ADM0_CODE", allow.cartesian=T)
      setkey(admin0, id, group, order, year)
      
      message('--> admin0 results merged to shapefile')
      
    }
  }
  
  ## admin1 estimates and shape file
  if ("admin1" %in% geo_levels) {
    message('loading admin1 data')
    
    if(missing(custom_data)) {
      pred <- '<<<< FILEPATH REDACTED >>>>' %>% fread
    } else {
      message('using custom data')
      pred <- custom_data$admin1
    }
    message('-> admin1 found and fread')
    
    pred <- pred[year %in% year_list, c('ADM0_CODE', 'ADM1_CODE', 'year', type), with = F]
    setnames(pred, type, 'outcome')
    
    if (map_region != 'global') {
      cty_code <- get_adm0_codes(map_region, shapefile_version = modeling_shapefile_version)
      pred <- pred[ADM0_CODE == cty_code]
    }
    
    if(use.sf==T) {
      
      admin1 <- get_admin_shapefile(1) %>% st_read
      admin1 <- admin1[admin1$ADM1_CODE %in% pred$ADM1_CODE,]
      admin1 <-merge(admin1, pred, by="ADM1_CODE", allow.cartesian=T)
      
      message('--> admin1 results merged to sf')
      
    } else {
      
      admin1 <- shapefile(get_admin_shapefile(1))
      admin1 <- admin1[admin1@data$ADM1_CODE %in% pred$ADM1_CODE,]
      admin1 <- SpatialPolygonsDataFrame(gSimplify(admin1, tol = 0.1, topologyPreserve = T), data = admin1@data)
      for (i in 1:length(admin1)) admin1@polygons[[i]]@ID <- as.character(admin1@data[i, "ADM1_CODE"])
      admin1 <- data.table(fortify(admin1))
      admin1[, ADM1_CODE := as.numeric(id)]
      admin1 <- merge(admin1, pred, by="ADM1_CODE", allow.cartesian=T)
      setkey(admin1, id, group, order, year)
      
      message('--> admin1 results merged to shapefile')
      
    }
  }
  
  ## admin2 estimates and shape file
  if ("admin2" %in% geo_levels) {
    message('loading admin2 data')
    
    if(missing(custom_data)) {
      pred <- '<<<< FILEPATH REDACTED >>>>' %>% fread
    } else {
      message('using custom data')
      pred <- custom_data$admin2
    }
    message('-> admin2 found and fread')
    message(year_list)
    
    pred <- pred[year %in% year_list, c('ADM0_CODE', 'ADM2_CODE', 'year', type), with = F]
    setnames(pred, type, 'outcome')
    
    if (map_region != 'global') {
      cty_code <- get_adm0_codes(map_region, shapefile_version = modeling_shapefile_version)
      pred <- pred[ADM0_CODE == cty_code]
    }
    
    if(use.sf==T) {
      
      admin2 <- get_admin_shapefile(2) %>% st_read
      admin2 <- admin2[admin2$ADM2_CODE %in% pred$ADM2_CODE,]
      admin2 <-merge(admin2, pred, by="ADM2_CODE", allow.cartesian=T)
      
      message('--> admin2 results merged to sf')
      
    } else {
      
      admin2 <- shapefile(get_admin_shapefile(2))
      admin2 <- admin2[admin2@data$ADM2_CODE %in% pred$ADM2_CODE,]
      admin2 <- SpatialPolygonsDataFrame(gSimplify(admin2, tol = 0.1, topologyPreserve = T), data = admin2@data)
      for (i in 1:length(admin2)) admin2@polygons[[i]]@ID <- as.character(admin2@data[i, "ADM2_CODE"])
      admin2 <- data.table(fortify(admin2))
      admin2[, ADM2_CODE := as.numeric(id)]
      admin2 <- merge(admin2, pred, by="ADM2_CODE", allow.cartesian=T)
      setkey(admin2, id, group, order, year)
      
      message('--> admin2 results merged to shapefile')
      
    }
  }
  
  ## combine and return all estimates
  mget(geo_levels) %>% return
  
}


## calc_diff_map -----------------------------------------------------------------------------------

calc_diff_map <- function(pred, diff_years) {
  diff <- lapply(names(pred), function(g) {
    rbindlist(lapply(diff_years, function(y) {
      temp <- pred[[g]][year %in% y, ]
      temp <- temp[, list(outcome = outcome[year == y[2]] - outcome[year == y[1]]), by=setdiff(names(temp), c("outcome", "year"))]
      temp[, years := paste(y, collapse="-")]
      temp
    }))
  })
  names(diff) <- names(pred)
  return(diff)
}


## plot_map ----------------------------------------------------------------------------------------
plot_map <- function(map_data, annotations, title, limits, map_region = map_region,
                     legend_colors, high_bad = T, legend_breaks, legend_labels, legend_title, custom_scale=F,
                     pop.mask=T, lake.mask=T, borders=T, stage3.mask=T,
                     zoom) {
  
  ## Enforce limits
  map_data$plot_var <- pmax(limits[1], pmin(limits[2], map_data$outcome)) #TODO set in to enforce lower limit as well?
  
  if (!custom_scale) {
    
    start_range <- range(map_data$outcome, na.rm = T)
    
    ## Create breaks
    breaks <- pretty(limits, 5)
    if (limits[1] < 0 & limits[2] > 0) breaks <- sort(unique(c(0, breaks)))
    
    ## Create labels
    labels <- format(breaks, nsmall = 2)
    if (min(limits) >= 0) divider <- "-" else divider <- " to "
    if (start_range[1] < limits[1]) {
      labels[1] <- paste0(format(floor(100*start_range[1])/100, nsmall=2), divider, labels[1])
    }
    if (start_range[2] > limits[2]) {
      labels[length(labels)] <- paste0(labels[length(labels)], divider, format(ceiling(100*start_range[2])/100, nsmall=2))
    }
    
  } else {
    map_data$plot_var <- map_data$outcome
    breaks <- legend_breaks
    labels <- legend_labels
  }
  
  ## Plot the base map (this is what shows in places with no estimates and no mask)
  canvas <- ggplot() + 
    geom_polygon(data = annotations$adm0, aes(x = long, y = lat, group = group), color = 'gray90', fill = 'gray90')
  
  ## Zoom
  if (!missing(zoom)) {
    canvas <- canvas + 
      xlim(zoom$x1, zoom$x2)  +
      ylim(zoom$y1, zoom$y2)
  }
  
  ## Plot predictions
  if ("group" %in% names(map_data)) {
    gg <- canvas + geom_polygon(data = map_data, aes(fill = plot_var, y = lat, x = long, group = group)) + 
      coord_equal(ratio = 1)
  } else if (class(map_data)[1] =='sf') {
    gg <- canvas + geom_sf(data = map_data, aes(fill = plot_var), lwd=0) + coord_sf(datum = NA)
  } else {
    gg <- canvas + geom_raster(data = map_data, aes(fill = plot_var, y = lat, x = long)) + 
      coord_equal(ratio = 1)
  }
  
  ## Plot mask, lakes, and adm boarders
  if (pop.mask==T) gg <- gg + annotate(geom = 'raster', x = annotations$mask$long, y = annotations$mask$lat, fill = 'gray70')
  if (lake.mask==T) gg <- gg + annotate(geom = 'raster', x = annotations$lakes$long, y = annotations$lakes$lat, fill = 'lightblue')
  if (borders==T) gg <- gg + geom_path(data = annotations$adm0, aes(x = long, y = lat, group = group), color = 'black', size = 0.2)
  if (stage3.mask==T) gg <- gg + geom_path(data = annotations$stage3, aes(x = long, y = lat, group = group), color = 'gray70', size = 0.2)
  
  ## Scales
  gg <- gg +
    scale_fill_viridis(option = legend_colors, direction = ifelse(high_bad, 1, -1), limits = range(breaks), breaks = breaks, labels = labels, name = legend_title)
  
  ## Labels & aesthetics
  gg <- gg +
    labs(x="", y="", title=title) +
    theme_classic() +
    theme(axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.position = c(0, 0), legend.justification = c(0, 0),
          legend.text=element_text(size=7),
          plot.title = element_text(hjust=0.5), plot.margin=unit(c(0, 0, 0, 0), "in")) +
    guides(fill = guide_colorbar(barwidth = 0.5, barheight = 5), color = guide_legend(override.aes=list(fill=NA)))
  
  return(gg)
  
}


## map_model_results -------------------------------------------------------------------------------

map_model_results <- function(indicator,
                              indicator_group,
                              run_date,
                              custom_data,
                              map_region = 'global',
                              type = 'mean',
                              raked = T,
                              lvl_years = c(2000, 2005, 2010, 2016),
                              lvl_colors = 'magma',
                              high_is_bad = TRUE,
                              lvl_limits = c(0, 0.25),
                              diff_years = list(c(2000, 2005), c(2005, 2010), c(2010, 2016), c(2000, 2016)),
                              diff_colors = 'viridis',
                              diff_limits = c(-0.10, 0.10),
                              include_diff = TRUE,
                              limits_type = "absolute",
                              geo_levels = c("raster", "admin1", "admin2"),
                              plot_by_year = TRUE,
                              plot_combined = TRUE,
                              file_type = "pdf",
                              raked_measure = NULL,
                              file_addin = NULL) {
  
  ## Quick argument checks
  if (!limits_type %in% c("absolute", "quantile")) stop("limits_type must be 'absolute' or 'quantile'")
  if (length(lvl_limits) != 2 | length(diff_limits) != 2) stop("lvl_limits & diff_limits must both be length 2")
  if (sum(!geo_levels %in% c("raster", "admin0", "admin1", "admin2")) > 0) stop("geo_levels can only include 'raster', 'admin0', 'admin1', and 'admin2'")
  if (!file_type %in% c("pdf", "png")) stop("file_type must be 'pdf' or 'png'")
  #if (!type %in% c("mean", "cirange", "cfb", "upper", "lower")) stop("type must be 'mean', 'cirange', or 'cfb', or 'upper', or 'lower'")
  
  ## Create output directory
  out_dir <- '<<<< FILEPATH REDACTED >>>>'
  dir.create(out_dir, showWarnings = F)
  
  ## Load data
  message("Loading data")
  annotations <- load_map_annotations(map_region = map_region)
  pred <- load_map_results(indicator, indicator_group, run_date, map_region, type, raked, 
                           start_year = lvl_years[1], end_year = lvl_years[length(lvl_years)], 
                           single_year = ifelse(length(lvl_years) == 1, lvl_years, 0), geo_levels = geo_levels, 
                           custom_data = custom_data, raked_measure = raked_measure)
  
  ## Make level maps
  message("Make levels maps")
  if (limits_type == "quantile") {
    lvl_limits <- quantile(unlist(lapply(pred, function(x) x$outcome)), probs = lvl_limits, na.rm = T)
    lvl_limits <- c(plyr::round_any(lvl_limits[1], 0.01, floor), plyr::round_any(lvl_limits[2], 0.01, ceiling))
  }
  
  legend_title <- c(mean = "Mean", cirange = "UI range", cfb = "CFB", upper = "Upper", lower = "Lower")[type]
  if (is.na(legend_title) | legend_title == 'NA') legend_title <- gsub('_', ' ', type)
  
  for (g in names(pred)) {
    message(paste0("...", g))
    
    # make maps and plot by year
    plot_list <- lapply(lvl_years, function(y) {
      message(paste0("......", y))
      gg <- plot_map(map_data = pred[[g]][year == y,], annotations = annotations, title = y,
                     legend_title = legend_title, limits = lvl_limits, legend_colors = lvl_colors, high_bad = high_is_bad,
                     pop.mask = ifelse(map_region == 'global', T, F), lake.mask = ifelse(map_region == 'global', T, F), borders=T, stage3.mask = ifelse(map_region == 'global', T, F))
      if (plot_by_year) {
        file_name <- paste0(out_dir, indicator, '_', type, if (raked) '_raked_' else '_unraked_', 
                            ifelse(is.null(raked_measure), '', paste0(raked_measure, '_')), g, '_', y, ifelse(is.null(file_addin), '', paste0('_', file_addin)), '.', file_type)
        if (file_type == "pdf") pdf(file_name, height=4, width=8)
        if (file_type == "png") png(file_name, height=4, width=8, units = "in", res = 1200)
        plot(gg)
        dev.off()
      }
      return(gg)
    })
    
    # plot combined
    if (plot_combined & length(lvl_years) > 1) {
      message("......combined")
      file_name <- paste0(out_dir, indicator, '_', type, if (raked) '_raked_' else '_unraked_', 
                          ifelse(is.null(raked_measure), '', paste0(raked_measure, '_')), g, '_combined', ifelse(is.null(file_addin), '', paste0('_', file_addin)), '.', file_type)
      if (file_type == "pdf") pdf(file_name, height=8, width=16)
      if (file_type == "png") png(file_name, height=8, width=16, units = "in", res = 1200)
      do.call("grid.arrange", plot_list)
      dev.off()
    }
    
    rm(plot_list)
  }
  
  ## Make difference maps
  if (include_diff) {
    message("Make differences maps")
    diff <- calc_diff_map(pred, diff_years)
    if (limits_type == "quantile") {
      diff_limits <- quantile(unlist(lapply(diff, function(x) x$outcome)), probs = diff_limits, na.rm = T)
      if (diff_limits[1] < 0 & diff_limits[2] > 0) diff_limits <- c(-1, 1)*max(abs(diff_limits)) # if crossing zero, ensure symmetry.
      diff_limits <- c(plyr::round_any(diff_limits[1], 0.01, floor), plyr::round_any(diff_limits[2], 0.01, ceiling))
    }
    
    legend_title <- paste0("Change in\n", legend_title)
    
    for (g in names(diff)) {
      message(paste0("...", g))
      
      # make maps and plot by year
      plot_list <- lapply(diff_years, function(y) {
        yrs <- paste(y, collapse = "-")
        message(paste0("......", yrs))
        gg <- plot_map(map_data = diff[[g]][years == yrs,], annotations = annotations, title = yrs,
                       legend_title = legend_title, limits = diff_limits, legend_colors = diff_colors, high_bad = high_is_bad,
                       pop.mask = ifelse(map_region == 'global', T, F), lake.mask = ifelse(map_region == 'global', T, F), borders=T, stage3.mask = ifelse(map_region == 'global', T, F))
        if (plot_by_year) {
          file_name <- paste0(out_dir, indicator, '_diff_', type, if (raked) '_raked_' else '_unraked_', 
                              ifelse(is.null(raked_measure), '', paste0(raked_measure, '_')), g, '_', y[1], '_', y[2], ifelse(is.null(file_addin), '', paste0('_', file_addin)), '.', file_type)
          if (file_type == "pdf") pdf(file_name, height=4, width=8)
          if (file_type == "png") png(file_name, height=4, width=8, units = "in", res = 1200)
          plot(gg)
          dev.off()
        }
        return(gg)
      })
      
      # plot combined
      if (plot_combined & length(diff_years) > 1) {
        message("......combined")
        file_name <- paste0(out_dir, indicator, '_diff_', type, if (raked) '_raked_' else '_unraked_', 
                            ifelse(is.null(raked_measure), '', paste0(raked_measure, '_')), g, '_combined', ifelse(is.null(file_addin), '', paste0('_', file_addin)), '.', file_type)
        if (file_type == "pdf") pdf(file_name, height=8, width=16)
        if (file_type == "png") png(file_name, height=8, width=16, units = "in", res = 1200)
        do.call("grid.arrange", plot_list)
        dev.off()
      }
      
      rm(plot_list)
    }
  }
  
  return("Maps saved!")
}