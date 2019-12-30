shiny_data_and_preds <- function(gaul_list, run_date, indicator, indicator_group, pred_file, layer_name) {
  
  library(ggplot2, lib.loc = package_lib)
  library(rgdal, lib.loc = package_lib)
  
  # Settings
  #color_list <- c("#000000","#00281D","#07425B","#38499A","#8149B9","#C653AF","#EB7190","#EC9F7D","#DCCF91","#DBF0C6")
  color_list <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
  if(indicator=='edu_mean') color_list <- rev(color_list)
  
  # extract admin0
  if(exists("subset_shape")==FALSE) {
    message("Opening master shapefile because not found in global env...")
    master_shape <- shapefile(paste0(root,"DATA/SHAPE_FILES/GBD_geographies/master/GBD_2016/master/shapefiles/GBD2016_analysis_final.shp"))
    subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% gaul_list, ]   
  }
  admin0.dt <- data.table(fortify(subset_shape)) 
  
  ## Logit functions
  logit <- function(x) {
    log(x/(1-x))
  }
  invlogit <- function(x) {
    exp(x)/(1+exp(x))
  }
  
  # Load actual data (df already in memory)
  if('<<< FILEPATH REDACTED >>>>')
  if('<<< FILEPATH REDACTED >>>>')
  plot_dir <- paste0(output_dir,"/plots")
  dir.create(plot_dir, showWarnings = FALSE)
  message(paste0('saving plots to ',plot_dir))
  message(paste0('immediately viewable at http://mbg-viz.duckdns.org:3456/ under your indicator/run_date'))
  
  periods <- data.frame(group = rep(1:length(unique(df$year)),5),years = rep(sort(unique(df$year)),5))
  df$period <- match(df$year, periods$years) # add these to df
  
  # Make quantity of interest
  df <- as.data.table(df)
  if(indicator_family=="binomial") df <- df[, to_map := get(indicator) / N]
  if(indicator_family=="gaussian") df <- df[, to_map := get(indicator)]
  
  # Plot data cluster means for each year
  plot.data <- function(x) {
    df.period <- subset(df, period == x)
    df.year <- df$year[df$period==x][1]
    loop.data.gg <- ggplot() +
      geom_polygon(data=admin0.dt, aes(x=long, y=lat, group=group), fill='grey90', color='grey') +
      geom_point(data=df.period, aes(x=longitude, y=latitude, color=to_map), pch=16, size=1) +
      scale_color_gradientn(colours=rev(color_list), limits=c(min(df[, to_map]), max(df[, to_map])), na.value = "white") +
      guides(fill=guide_colorbar(title=indicator, label=TRUE, ticks=FALSE)) +
      coord_fixed() +
      ggtitle(df.year) +
      guides(size=FALSE) +
      theme_minimal() +
      theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.text.x=element_blank()) + theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines"))
    return(loop.data.gg)
  }
  for(period in sort(unique(df$period))) {
    assign(paste("data.gg", period, sep="."),plot.data(period))
  }
  
  f <- paste0(output_dir,'/',pred_file)
  preds <- brick(f)
  #preds <- setExtent(preds, subset_shape)
  
  # Convert raster to SpatialPointsDataFrame
  preds.sp <- rasterToPoints(preds, spatial=TRUE)
  projection <- proj4string(preds.sp)
  
  # reproject sp object  
  preds.sp <- spTransform(preds.sp, CRS(projection)) 
  preds.sp@data <- data.frame(preds.sp@data, long=coordinates(preds.sp)[,1],lat=coordinates(preds.sp)[,2]) 
  preds.dt <- data.table(preds.sp@data)
  
  ## Plot preds of proportion with 0 years of education
  names(preds.dt)[names(preds.dt) == "lat"] = "latitude"
  names(preds.dt)[names(preds.dt) == "long"] = "longitude" 
  
  # Plot predictions for all periods
  plot.preds <- function(x) {
    period <- period <- strsplit(x, '[.]')[[1]][2]
    df.year <- df$year[df$period==period][1]
    loop.preds.gg <- ggplot(preds.dt,aes(longitude,latitude)) +
      geom_raster(aes(fill=get(x))) +
      coord_fixed() + 
      theme_minimal() +
      geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='white', lwd=.1) +
      scale_fill_gradientn(colours=rev(color_list), limits=c(min(minValue(preds)), max(maxValue(preds))), na.value = "white") +
      guides(fill=guide_colorbar(title=indicator, label=TRUE, ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      ggtitle(df.year) +
      theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) +
      theme(legend.position='bottom', legend.direction='horizontal')
    return(loop.preds.gg)
  }
  for(i.period in sort(unique(df$period))) {
    assign(paste("preds.gg", i.period, sep="."),plot.preds(paste0(layer_name,i.period)))
  }
  
  # Make data and preds pngs for Shiny
  png(paste0(plot_dir,'/data1.png'),width=400)
  print(data.gg.1)
  dev.off()
  png(paste0(plot_dir,'/data2.png'),width=400)
  print(data.gg.2)
  dev.off()
  png(paste0(plot_dir,'/data3.png'),width=400)
  print(data.gg.3)
  dev.off()
  png(paste0(plot_dir,'/data4.png'),width=400)
  print(data.gg.4)
  dev.off()
  png(paste0(plot_dir,'/preds1.png'),width=400)
  print(preds.gg.1)
  dev.off()
  png(paste0(plot_dir,'/preds2.png'),width=400)
  print(preds.gg.2)
  dev.off()
  png(paste0(plot_dir,'/preds3.png'),width=400)
  print(preds.gg.3)
  dev.off()
  png(paste0(plot_dir,'/preds4.png'),width=400)
  print(preds.gg.4)
  dev.off()
  
}

shiny_cov_layers <- function(fixed_effects, gaul_list, run_date, indicator, indicator_group) {
  
  # Plot covariate layers
  library(ggplot2, lib.loc = package_lib)
  
  color_list <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
  
  # Make sure covs are loaded
  load('<<< FILEPATH REDACTED >>>>')
  
  # Load actual data (df already in memory)
  if('<<< FILEPATH REDACTED >>>>')
  if('<<< FILEPATH REDACTED >>>>')
  plot_dir <- paste0(output_dir,"/plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Get templates
  if(exists("subset_shape")==FALSE) {
    message("Opening master shapefile because not found in global env...")
    master_shape <- shapefile(paste0(root,"DATA/SHAPE_FILES/GBD_geographies/master/GBD_2016/master/shapefiles/GBD2016_analysis_final.shp"))
    subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% gaul_list, ]   
  }
  admin0.dt <- data.table(fortify(subset_shape)) 
  
  # Plot varying covariates
  gLegend<-function(a.gplot){
    pdf(NULL) #  
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    graphics.off()
    return(legend)
  }
  selected_covs <- strsplit(fixed_effects," ")
  selected_covs <- selected_covs[[1]][selected_covs[[1]] != "+"]
  library(grid)
  tv_count <- 0
  for(c in selected_covs) {
    if(paste0('tv_',c) %in% grep('tv_*', ls(), value = TRUE)) {
      
      tv_count <- tv_count + 1
      tv_cov <- get(paste0('tv_',c))
      
      # Convert raster to SpatialPointsDataFrame
      preds.sp <- rasterToPoints(tv_cov, spatial=TRUE)
      projection <- proj4string(preds.sp)
      
      # reproject sp object  
      preds.sp <- spTransform(preds.sp, CRS(projection)) 
      preds.sp@data <- data.frame(preds.sp@data, long=coordinates(preds.sp)[,1],lat=coordinates(preds.sp)[,2]) 
      preds.dt <- data.table(preds.sp@data)
      
      ## Plot preds of proportion with 0 years of education
      names(preds.dt)[names(preds.dt) == "lat"] = "latitude"
      names(preds.dt)[names(preds.dt) == "long"] = "longitude" 
      
      # Plot predictions for all periods
      plot.preds <- function(x) {
        period <- gsub(paste0(c,'.'), "", x)
        loop.preds.gg <- ggplot(preds.dt,aes(longitude,latitude)) +
          geom_raster(aes(fill=get(x))) +
          coord_fixed() + 
          theme_minimal() +
          geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='white', lwd=.1) +
          scale_fill_gradientn(colours=rev(color_list), limits=c(min(minValue(tv_cov)), max(maxValue(tv_cov))), na.value = "grey") + 
          guides(fill=guide_colorbar(title=c, label=TRUE, ticks=FALSE)) +
          scale_x_continuous("", breaks=NULL) +
          scale_y_continuous("", breaks=NULL) +
          theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) + 
          theme(legend.position='bottom', legend.direction='horizontal') + 
          ggtitle(paste0("Period ",period))
        return(loop.preds.gg)
      }
      for(i.period in 1:4) {
        assign(paste0("tv_cov_",tv_count,".gg.", i.period),plot.preds(paste0(c,'.',i.period)))
      }
    }    
  }
  
  # Plot all varying covariates
  png(paste0(plot_dir,'/tv_covs.png'),width=1600,height=1600)
  total_rows <- tv_count*5
  # Initialize plot with master title
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(total_rows, 22, heights=c(.25,.25,.25,.25), widths=c(.25,.25,.25,.25))))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  # Plot all data coverage maps
  start_row <- 1
  for(i in 1:tv_count) {
    end_row <- start_row + 4
    print(get(paste0('tv_cov_',i,'.gg.1')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 1:5))
    print(get(paste0('tv_cov_',i,'.gg.2')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 6:10))
    print(get(paste0('tv_cov_',i,'.gg.3')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 11:15))
    print(get(paste0('tv_cov_',i,'.gg.4')) + theme(legend.position="none"), vp = vplayout(start_row:end_row, 16:20))
    p.legend <- gLegend(get(paste0('tv_cov_',i,'.gg.1')))
    p.legend$vp <- viewport(layout.pos.row = start_row:end_row, layout.pos.col = 21:22)
    grid.draw(p.legend)
    start_row <- start_row + 5
  }
  dev.off()
  
  nt_covs <- covs # non-varying
  ntv_plot_function <- function(x) {  
    covs.sp <- rasterToPoints(nt_covs[[x]], spatial=TRUE)
    projection <- proj4string(covs.sp)
    
    # reproject sp object  
    covs.sp <- spTransform(covs.sp, CRS(projection)) 
    covs.sp@data <- data.frame(covs.sp@data, long=coordinates(covs.sp)[,1],lat=coordinates(covs.sp)[,2]) 
    covs.dt <- data.table(covs.sp@data)
    
    ## Plot preds of proportion with 0 years of education
    names(covs.dt)[names(covs.dt) == "lat"] = "latitude"
    names(covs.dt)[names(covs.dt) == "long"] = "longitude" 
    
    covs.gg <- ggplot(covs.dt,aes(longitude,latitude)) +
      geom_raster(aes(fill=get(x))) +
      coord_fixed() + 
      theme_minimal() +
      geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='white', lwd=.1) +
      scale_fill_gradientn(colours=(color_list), limits=c(minValue(nt_covs[[x]]), maxValue(nt_covs[[x]])), na.value = "grey") + 
      guides(fill=guide_colorbar(title=x, label=TRUE, ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) +
      theme(legend.position='bottom', legend.direction='horizontal')
    return(covs.gg)
  }
  
  list_of_ntv_ggs <- lapply(names(nt_covs), ntv_plot_function)
  
  # Plot all non-varying covariates
  png(paste0(plot_dir,'/ntv_covs.png'),width=1200)
  total_columns <- length(names(nt_covs))*6
  # Initialize plot with master title
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(5, total_columns, heights=c(.25,.25,.25,.25), widths=c(.25,.25,.25,.25))))
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  # Plot all data coverage maps
  start_col <- 1
  for(i in 1:length(names(nt_covs))) {
    end_col <- start_col + 4
    print(list_of_ntv_ggs[[i]] + theme(legend.position="none"), vp = vplayout(1:5, start_col:end_col))
    p.legend <- gLegend(list_of_ntv_ggs[[i]])
    p.legend$vp <- viewport(layout.pos.row = 1:5, layout.pos.col = end_col+1)
    grid.draw(p.legend)
    start_col <- end_col+2
  }
  dev.off()
  
}


shiny_raked <- function(gaul_list, run_date, indicator, indicator_group, pred_file, layer_name) {
  
  library(ggplot2, lib.loc = package_lib)
  library(rgdal, lib.loc = package_lib)
  
  # Settings
  #color_list <- c("#000000","#00281D","#07425B","#38499A","#8149B9","#C653AF","#EB7190","#EC9F7D","#DCCF91","#DBF0C6")
  color_list <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
  if(indicator=='edu_mean') color_list <- rev(color_list)
  
  # extract admin0
  if(exists("subset_shape")==FALSE) {
    message("Opening master shapefile because not found in global env...")
    master_shape <- shapefile(paste0(root,"DATA/SHAPE_FILES/GBD_geographies/master/GBD_2016/master/shapefiles/GBD2016_analysis_final.shp"))
    subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% gaul_list, ]   
  }
  admin0.dt <- data.table(fortify(subset_shape)) 
  
  # Load actual data (df already in memory)
  if('<<< FILEPATH REDACTED >>>>')
  if('<<< FILEPATH REDACTED >>>>')
  plot_dir <- paste0(output_dir,"/plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  f <- paste0(output_dir,'/',pred_file)
  preds <- brick(f)

  # Convert raster to SpatialPointsDataFrame
  preds.sp <- rasterToPoints(preds, spatial=TRUE)
  projection <- proj4string(preds.sp)
  
  # reproject sp object  
  preds.sp <- spTransform(preds.sp, CRS(projection)) 
  preds.sp@data <- data.frame(preds.sp@data, long=coordinates(preds.sp)[,1],lat=coordinates(preds.sp)[,2]) 
  preds.dt <- data.table(preds.sp@data)
  
  ## Plot preds of proportion with 0 years of education
  names(preds.dt)[names(preds.dt) == "lat"] = "latitude"
  names(preds.dt)[names(preds.dt) == "long"] = "longitude" 
  
  ## Define maximum value if we rake over 1
  max_raked_value <- max(maxValue(preds))
  if(max_raked_value>1 & indicator_family=='binomial') max_raked_value <- 1
  if(max_raked_value>18 & indicator_family=='gaussian' & indicator=='edu_mean') max_raked_value <- 18
  
  # Plot predictions for all periods
  plot.preds <- function(x) {
    period <- strsplit(x, '[.]')[[1]][2]
    df.year <- df$year[df$period==period][1]
    loop.preds.gg <- ggplot(preds.dt,aes(longitude,latitude)) +
      geom_raster(aes(fill=get(x))) +
      coord_fixed() + 
      theme_minimal() +
      geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='white', lwd=.1) +
      scale_fill_gradientn(colours=rev(color_list), limits=c(min(minValue(preds)), max_raked_value), na.value = "#a6d854") +
      guides(fill=guide_colorbar(title=indicator, label=TRUE, ticks=FALSE)) +
      scale_x_continuous("", breaks=NULL) +
      scale_y_continuous("", breaks=NULL) +
      ggtitle(df.year) +
      theme(panel.margin = unit(0, "lines"), plot.margin = unit(c(0,0,0,0),"lines")) +
      theme(legend.position='bottom', legend.direction='horizontal')
    return(loop.preds.gg)
  }
  for(i.period in sort(unique(df$period))) {
    assign(paste("preds.gg", i.period, sep="."),plot.preds(paste0(layer_name,i.period)))
  }
  
  # Make data and preds pngs for Shiny
  png(paste0(plot_dir,'/raked1.png'),width=400)
  print(preds.gg.1)
  dev.off()
  png(paste0(plot_dir,'/raked2.png'),width=400)
  print(preds.gg.2)
  dev.off()
  png(paste0(plot_dir,'/raked3.png'),width=400)
  print(preds.gg.3)
  dev.off()
  png(paste0(plot_dir,'/raked4.png'),width=400)
  print(preds.gg.4)
  dev.off()
  
}

shiny_data_scatter <- function(df, run_date, indicator, indicator_group, year_var) {
  
  # Make dirs
  load('<<< FILEPATH REDACTED >>>>')
  output_dir <- paste0('<<< FILEPATH REDACTED >>>>')
  plot_dir <- paste0(output_dir,"/plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Make data coverage scatter
  df <- as.data.table(df)
  df[, clusters:=1]
  total_clusters <- df[, list(clusters=sum(clusters)), by=c('original_year','country','source')]
  total_cluster_num <- total_clusters[, list(clusters=sum(clusters))]
  png(paste0(plot_dir,'/data_scatter.png'),width=1200,height=800)
  clusters.gg <- ggplot() +
    geom_point(data=total_clusters, aes(x=original_year, y=country, size=clusters, shape=factor(source), color=factor(source))) +
    guides(color=FALSE) + 
    scale_size(guide = guide_legend(title = "Clusters/polygons"), range=c(1,10)) + 
    ggtitle(paste0("Clusters/polygons by country/year, total points = ", total_cluster_num[1, clusters])) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank()) 
  clusters.gg
  dev.off()
  
}

shiny_raking_map <- function() {
  
  # TRY PLOTTING BY MERGING TO GEO SHAPEFILE
  require(sp)
  library(rgeos)
  library(maptools)
  library(ggplot2)
  
  # extract admin0
  if(exists("subset_shape")==FALSE) {
    message("Opening master shapefile because not found in global env...")
    master_shape <- shapefile(paste0(root,"DATA/SHAPE_FILES/GBD_geographies/master/GBD_2016/master/shapefiles/GBD2016_analysis_final.shp"))
    subset_shape <- master_shape[master_shape@data$GAUL_CODE %in% gaul_list, ]   
  }
  admin0.dt <- data.table(fortify(subset_shape)) 
  
  gLegend<-function(a.gplot){
    pdf(NULL) #  
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    graphics.off()
    return(legend)
  }
  
  if('<<< FILEPATH REDACTED >>>>')
  if('<<< FILEPATH REDACTED >>>>')
  plot_dir <- paste0(output_dir,"/plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  names(rf)[names(rf)=="name"] <- "GAUL_CODE"
  rf$raking_factor[is.infinite(rf$raking_factor)] <- -1
  
  # Loop over periods
  gg.count <- 1
  for(period in unique(rf$year)) {
  
  rf_shape <- merge(subset_shape,rf[rf$year==period,], by="GAUL_CODE")
  rf_shape@data$id = rownames(rf_shape@data)
  rf.pts <- fortify(rf_shape, region="id")
  rf.df = join(rf.pts, rf_shape@data, by="id")
  rf.df$raking_factor[is.na(rf.df$raking_factor)] <- 0
  rf.df$raking_factor <- as.numeric(rf.df$raking_factor)
  
  # Plot
  color_list <- c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
  admin2.gg <- ggplot(rf.df,aes(long,lat,group=group)) + 
    geom_polygon(aes(fill=raking_factor)) + 
    geom_path(data=admin0.dt, aes(x=long, y=lat, group=group), color='black', size=.5) +
    scale_fill_gradientn(colours=color_list, limits=c(0, max(rf$raking_factor[!is.na(rf$raking_factor)])), na.value = "#000000") + 
    guides(fill=guide_colorbar(label=TRUE, ticks=FALSE)) +
    ggtitle(period) +
    scale_x_continuous("", breaks=NULL) +
    scale_y_continuous("", breaks=NULL) +
    coord_equal() 
  # grab your legends using the predefined functions, then state their grid location
  p.legend <- gLegend(admin2.gg)
  assign(paste0("ggplot.",gg.count),admin2.gg)
  gg.count <- gg.count + 1
  
  }
  
  # Make data and preds pngs for Shiny
  png(paste0(plot_dir,'/rf_map1.png'),width=400)
  print(ggplot.1)
  dev.off()
  png(paste0(plot_dir,'/rf_map2.png'),width=400)
  print(ggplot.2)
  dev.off()
  png(paste0(plot_dir,'/rf_map3.png'),width=400)
  print(ggplot.3)
  dev.off()
  png(paste0(plot_dir,'/rf_map4.png'),width=400)
  print(ggplot.4)
  dev.off()

}



