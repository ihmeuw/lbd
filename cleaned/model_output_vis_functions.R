
package_list <- c('rgdal', 'raster', 'ggplot2', 'grid', 'gridExtra', 'RColorBrewer', 'stringr', 'data.table')
suppressPackageStartupMessages(lapply(package_list, library, character.only = TRUE))

.libPaths('FILEPATH')
library(Cairo)
library(scales)

#create wrapper to fit long titles
wrap_text <- function(txt, width) {
  paste(strwrap(txt, width), collapse = "\n")
}


#create capitalization function
capitalize <- function(string){
  first_letter <- substr(string, 1, 1)
  rest <- substr(string, 2, nchar(string))
  fl_upper <- toupper(first_letter)
  capitalized <- paste0(fl_upper, rest)
  return(capitalized)
}



###HISTOGRAM FUNCTION ARGUMENTS###
#dat should be a dat1.preds or dat0.preds object
#out_dir example:'FILEPATH' should *not* end with /
#bg is automatically set to false, if using dat0.preds set to TRUE for appropriate titling
#file_name should be formatted '/what_you_want_your_file_named' it will be automatically saved as jpg this will automatically add to whatever out_dir is set to
#device graphics device for printing, default is 'jpg', but anything used by ggsave will work

#prediction histogram function
#call the datX.preds object and specify bg as TRUE/FALSE
get_pred_hist <- function(dat, 
                          out_dir, 
                          bg=FALSE, 
                          title = NULL, 
                          file_name=NULL,
                          device = 'jpg'){
  #convert data to df to be compatible with ggplot
  preds_dat <- data.frame(dat)
  colnames(preds_dat) <- 'pred_value'
  x_lab <- "Extracted Value"
  if(bg){
    #write y-axis label and title based on T/F of bg if T then set for background if F then occurrence
    y_lab <- "Background Points"
    if(is.null(title)){
      title <- paste0("Extracted Values of Background Data for ", dz)
    } 
    
  }
  else{
    y_lab <- "Occurrence Points"
    if(is.null(title)){
      title <- paste0("Extracted Values of Occurrence Data for ", dz)
    }
    
  }
  if(monthly){
    #alter title if data is monthly
    title <- paste0(title, " - ", substr(month.name[as.numeric(mo)], 1, 3))
  }
  #format ggplot
  pred_hist <- ggplot(data=preds_dat, aes(x=preds_dat$pred_value)) +
    geom_histogram(breaks = seq(0,1,0.05), color="dark grey", fill="light grey") +
    ggtitle(wrap_text(title, width=40))+
    labs(x=x_lab, y=y_lab) +
    theme(plot.title=element_text(size=12, face="bold", vjust=1, hjust=0.5),
          axis.title.x=element_text(size=10, margin = margin(t = 12, r = 0, b = 8, l = 0)),
          axis.title.y=element_text(size=10, margin = margin(t = 0, r = 20, b = 0, l = 0)),
          axis.text.x=element_text(size=8),
          axis.text.y=element_text(size=8),
          panel.background = element_rect(fill="gray95"),
          panel.grid.major = element_line(color="white"),
          panel.grid.minor = element_line(color="white"))
  if(!is.null(out_dir)){
  #decide filepath based on filename and if monthly
  if (monthly){
    file_path <- paste0(out_dir, "/prediction_hist_", mo, ".", device)
  } else {
    file_path <- paste0(out_dir, "/prediction_hist.", device)
  }
  if (!is.null(file_name)){
    file_path <- sub("/prediction_hist", file_name, file_path)
  }
  #ggsave at 300 dpi
  
  
    ggsave(file_path, pred_hist, dpi=350, width = 8.5, height = 6)
    file_path <- sub(paste0('.', device), paste0('_low_res.', device), file_path)
    ggsave(file_path, pred_hist, dpi=100, width = 8.5, height = 6)
    
  } else{
    return(pred_hist)
  }
  
}


###EFFECT HIST PLOT FUNCTION ARGUMENTS###
#point_data should be dat.pts this will be used for the histogram or any other data formatted for a column per covariate
#plot_data should be effect which will be used for the effect plot and there should be a list item named for each covariate
#order should be a vector of the numerical order of covariates e.g. c(4, 3, 1, 6, 5, 2)
#out_dir example:'FILEPATH' should *not* end with /
#file_name defaults to "/FILEPATH" format should be "/file_name" this will automatically add to whatever out_dir is set to
#device graphics device used to print image, default is 'jpg', anything used by ggsave will work

#if there are more than 10 covariates the layout won't pick a numrow or numcol for the grid layout and ask you to specify those numbers manually in the environment


#function for producing effect histogram combo plots
get_effect_hist_plot <- function(point_data, 
                                 plot_data, 
                                 order, 
                                 out_dir, 
                                 main_title = NULL, 
                                 file_name=NULL,
                                 device = 'jpg'){
  plot_list <- list()
  ordered_covs <- cov_names[order]
  for (col_name in ordered_covs){
    print(col_name)
    if (col_name == "bat_frequency_mean"){
      ind <- which(cov_name_list$cov_name=="bat_frequency_mean")
      cov_name_long <- "Mean bat frequency"
    } else if (col_name == 'excess_sum2'){
      ind <- which(cov_name_list$cov_name=="excess_sum2")
      cov_name_long <- "Rainfall (sum prev. 2 months)"
    } else {
      if(!is.null(cov_name_list$title)){
        cov_name_long <- cov_name_list$cov_name_long[grep(col_name, cov_name_list$title)]
      } else {
        cov_name_long <- cov_name_list$cov_name_long[grep(col_name, cov_name_list$cov_name)]
      }
       #match short cov names to long names for titling
    }
    hist_data <- data.frame("cov" = point_data[,col_name]) #select column of specific covariate to be plotted
    effect_data <- data.frame(plot_data[[col_name]][,1:4]) #select columns of covariate effect to be plotted
    names(effect_data) <- c("covariate", "mean", "lower", "upper") #rename for ease of access
    #designate titles and x/y axis labels and transformations based on covariate
    title <- cov_name_long
    inset_text <- paste0("Relative Influence = " ,round(relinf[col_name, 1], digits=2))
    rel_inf_inset <- grobTree(textGrob(inset_text, x=0.38, y=0.94, hjust=0, gp=gpar(fontsize=10, fontface="bold", col="gray35")))
    ylab1 <- "Frequency Density"
    ylab2 <- "Effect"
    if (col_name == "crutstmp"){
      xlab <- "Temperature (°C)"
      hist_data <- (hist_data/10)
      effect_data$covariate <- (effect_data$covariate/10)
    }
    else if (col_name == "crutsard"){
      xlab <- "Precip. / Evap. (mm/mm)"
    }
    else if (col_name == "bat_frequency_mean"){
      xlab <- "Suitability / all species"
    } else{
      if(!is.null(cov_name_list$title)){
        xlab <- capitalize(cov_name_list$units[grep(col_name, cov_name_list$title)])
      } else {
        xlab <- capitalize(cov_name_list$units[grep(col_name, cov_name_list$cov_name)])
      }
      
    }
    
    if(length(xlab > 1)){
      xlab <- xlab[1]
    }
    
    
    
    if(xlab == '100days'){
      xlab <- "Days"
      hist_data <- (hist_data/100)
      effect_data$covariate <- (effect_data$covariate/100)
    }
    
    if(xlab=="Mm" | xlab=="Mm/day" | xlab == 'Km' | xlab == 'M'){
      xlab <- tolower(xlab)
    }
    
   
    
    x_min <- min(effect_data[,"covariate"])
    x_max <- max(effect_data[,"covariate"])
    hist_min <- min(hist_data[,"cov"])
    hist_max <- max(hist_data[,"cov"])
    #plot the histogram
    plt <- ggplot(data=hist_data) +
      geom_histogram(aes(x=cov),
                     breaks = seq(x_min, x_max, (x_max-x_min)/20),
                     color="dark grey",
                     fill="light grey")+
      scale_x_continuous(limits = c(x_min,x_max))


    #set up transformation variables based on max and min values and max bin height of the hist
    max_bin <- max(ggplot_build(plt)$data[[1]]$count)
    print(max_bin)
    if (max_bin>10000){
      max_bin_round <- round(max_bin, digits=-4)
    } else if (max_bin > 1000){
      max_bin_round <- round(max_bin, digits=-3)
    } else if(max_bin > 100){
      max_bin_round <- round(max_bin, digits=-2)
    } else if(max_bin < 100){
      max_bin_round <- round(max_bin)
    }
    max_effect_val <- round(max(effect_data$upper))
    min_effect_val <- round(min(effect_data$lower))
    trans <- (max_bin_round/(max_effect_val-min_effect_val)) * .4
    int <- max_bin_round/2
    resc_int <- int/trans
    #finish by adding the transformed effect plot w/ CI to the histogram
    plt <- plt +
      geom_line(data = effect_data,
                aes(x=covariate, y=(mean*trans)+int),
                color="midnightblue") +
      geom_ribbon(data = effect_data,
                  aes(x=covariate,
                      ymin=(lower*trans)+int,
                      ymax=(upper*trans)+int),
                  alpha=0.6) +
      scale_alpha_continuous(guide=FALSE) +
      scale_y_continuous(sec.axis = sec_axis(trans=~./trans-resc_int, name=ylab2), name=ylab1)+
      ggtitle(title)+
      labs(x=xlab)+
      theme(plot.title=element_text(size=12, vjust=1),
            axis.title.x=element_text(size=12, vjust=-2, margin = margin(t = 12, r = 0, b = 8, l = 0)),
            axis.title.y=element_text(size=12, margin = margin(t = 0, r = 20, b = 0, l = 0)),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10),
            panel.background = element_rect(fill="gray95"),
            panel.grid.major = element_line(color="white"),
            panel.grid.minor = element_line(color="white"),
            plot.margin = unit(c(10,10,10,10),"points"))+
      annotation_custom(rel_inf_inset)
    plt_grobbed <- ggplotGrob(plt)
    plot_list[[col_name]] <- plt_grobbed
  }

  
  #title
 
  if(is.null(main_title)){
    
    main_title <- paste0("Mean Effect Plots and Histograms of Covariate Values for ", dz)
    
  }

  main_title <- textGrob(wrap_text(main_title, 55), gp=gpar(fontface="bold", size=16))
  
  #plot grid
  plot_grid <- marrangeGrob(grobs=plot_list, top=main_title, bottom = quote(textGrob(paste("page", g, "of", npages))), nrow=2, ncol=3)
  
  #save grid page by page if more than one page
  if(!is.null(out_dir)){
    if (length(plot_grid)>1){
      for (i in 1:length(plot_grid)){
        if (monthly){
          file_path <- paste0(out_dir, "/effect_hist_grid_", mo, "_page_", i, ".", device)
        } else {
          file_path <- paste0(out_dir, "/effect_hist_grid_page_", i, ".", device)
        }
        ggsave(file_path, plot_grid[[i]], width=40/3, height = 8, dpi=350)
        file_pathlr <- sub(paste0('.', device), paste0('_low_res.', device), file_path)
        ggsave(file_pathlr, plot_grid[[i]], width=40/3, height = 8, dpi=100)
      }
    } else{
      if (monthly){
        file_path <- paste0(out_dir, "/effect_hist_grid_", mo, ".", device)
      } else {
        file_path <- paste0(out_dir, "/effect_hist_grid.", device)
      }
      if (!is.null(file_name)){
        file_path <- sub("/effect_hist_grid", file_name, file_path)
      }
      ggsave(file_path, plot_grid, width=40/3, height = 8, dpi=350)
      file_path <- sub(paste0('.', device), paste0('_low_res.', device), file_path)
      ggsave(file_path, plot_grid, width=40/3, height = 8, dpi=100)
    }
  } else{
    return(plot_grid)
  }
  
  
}



###Begin prep for raster_map function###

#read in shapefiles for mapping
if(!exists('bg_polygon') | !exists('bg_polygon_df')){
  #read in and convert to df background poly shapefile
  bg_polygon <- shapefile("FILEPATH.shp")
 
  bg_polygon_df <- fortify(bg_polygon)

}

if(!exists('lakes')| !exists('lakes_df')){
  #read in and convert to df lakes shapefile
  lakes <- readOGR(dsn="FILEAPTH", layer="FILEPATH")
  lakes <- lakes[which(lakes@data$POLY_SRC != 'Circle'),]
  lakes_df <- fortify(lakes)
}

if(!exists('disputed') | !exists('disputed_df')){
  disputed <- shapefile("FILEAPTH.shp")
  disputed_df <- fortify(disputed)
}

### MAP FUNCTION ARGUMENTS ###
#raster_layer           should specify a particular layer such as preds_sry[["mean"]]
#out_dir                example:'FILEPATH' should *not* end with /
#color_scheme           should be passed a vector of HEX color codes (three for diverging or two for sequential low - high) e.g. if you want green-yellow-purple -> c("#00681d", "#fbfcc9", "#9b0065") (default)
#bg_pts                 any point data with 'lat' and 'long' columns will work ideally dat0.pts
#occ_pts                any point data with 'lat' and 'long' columns will work ideally dat1.pts
#good_pts               any point data with 'lat' and 'long' columns will work ideally dat1.good
#bad_pts                any point data with 'lat' and 'long' columns will work ideally dat1.bad
##NOTE: bad_pts and good_pts are designed to be used together and the legend for good_pts will not appear if bad_pts not specified ##
#crop                   default is set to TRUE, will crop the raster to the extent of the buffer for the specified experiment
#clip                   default is set to TRUE, will clip the raster to the shape of the buffer for the specified experiment
##NOTE: if clip is used crop should also be used to speed up the time it takes to run the clip function ##
#extents                default is set to the extent of the raster layer after crop/clip occurs, must be entered as an extent(object you want to match extents) object
#legend_position        default is set to c(0.03, 0.03) which is the bottom left corner of the map, format of legend_position should either be c(x_coord, y_coord) or 'top', 'left', 'bottom', 'right'
#file_name              default set to "/prediction_map" format should be "/file_name" this will automatically add to whatever out_dir is set to
#poly_layer             can receive a SpatialPolygonsDataFrame such as those that are produced by readOGR which will add the polygon borders as a final layer to the map
#title                  default set to "Suitability of dz", format "My Map Has a Custom Title"
#lims                   a numeric vector of length 2 setting the limits of values that will be shown, if NULL will be 0-1 ex. c(0,500)
#labs                   a vector of strings to label the colorbar ex. c('Low', 'Medium', 'High')
#breaks                 a vector of the same length as labs which indicates the values the labels correspond to e.g. c(0, 0.5, 1)
#letter                 a vector of length 3 containing the text, x and y coordinates in that order ex. c('A', 45.2, 57.6, 0) and hjust
#na_color               a single color string
#transform              a string containing a function to transform raster values with
#raster2                a raster layer displayed as a semi-transparent population mask
#buffer                 specify a buffer with and spdf
#size                   size specifies number of inches in width to make the plot, height is calculated from that proportionally based on plot dimensions, default is 10
#device                 graphics device string, default is 'jpg', can also be 'pdf', 'png', etc. anything used with Cairo
#resc                   rescale colorbar to be asymmetric default = c(min(raster_df[,3]), max(raster_df[,3], diff(c(min(raster_df[,3]), max(raster_df[,3])/length(colors))))


#function to create map
get_raster_map <- function(raster_layer, 
                           out_dir, 
                           color_scheme = NULL, 
                           bg_pts=NULL, 
                           occ_pts=NULL, 
                           good_pts=NULL, 
                           bad_pts=NULL, 
                           crop = TRUE, 
                           clip=TRUE, 
                           extents = NULL, 
                           expand=TRUE,  
                           legend_position = NULL, 
                           file_name=NULL, 
                           poly_layer=NULL, 
                           title=NULL, 
                           lims = NULL, 
                           labs = NULL, 
                           breaks = NULL, 
                           letter = NULL, 
                           show_leg = TRUE, 
                           plot_title = NULL,
                           na_color = NULL, 
                           transform = NULL, 
                           raster2 = NULL, 
                           buffer = NULL, 
                           size = 10,
                           co_lines = c('gray20', 0.25),
                           device = 'jpg',
                           resc = NULL){

  buffer_smooth <- buffer

  #crop raster if crop=TRUE
  if(crop){
    
    
    if(!is.null(buffer)){
      cext <- extent(buffer_smooth)
    } else if(is.null(buffer) & !is.null(extents)){
      cext <- extents
    } else if (is.null(extents) & is.null(buffer)){
      stop('Nothing to crop raster layer by, please add extents or buffer_smooth sp object')
    }
    
    if(expand){
      cext@xmin <- cext@xmin - 30
      cext@xmax <- cext@xmax + 30
      cext@ymin <- cext@ymin - 30
      cext@ymax <- cext@ymax + 30
    }
    cropped_raster <- crop(raster_layer, cext)
    raster_layer <- cropped_raster
  }

  #clip raster if clip =TRUE
  if(clip){
    if(!is.null(buffer)){
      #clip the raster layer to buffer
      clipped_raster <- raster::mask(cropped_raster, buffer_smooth)
      clip_inverse <- mask(cropped_raster, buffer_smooth, inverse = TRUE)
      raster_layer <- clipped_raster
    } else{
      stop('no buffer provided to clip by')
    }
    
  }


  #convert raster to data.frame and omit NAs
  raster_df <- na.omit(as.data.frame(raster_layer, xy=TRUE))
  names(raster_df)[3] <- 'layer'

  #find map extent based on preds if extents aren't given
  if(is.null(extents)){
    if(is.null(buffer_smooth)){
      extent_ref <- extent(raster_layer)
    } else{
      extent_ref <- extent(buffer_smooth)
    }
    
  } else {
    extent_ref <- extents
  }
  xlims <- c(extent_ref@xmin, extent_ref@xmax)
  ylims <- c(extent_ref@ymin, extent_ref@ymax)
  
  if(!is.null(out_dir)){
    if(grepl('04_rift_valley_fever', out_dir)){
      if(min(ylims) < -38){
        ylims[1] <- -38
      }
      if(min(xlims) < -25){
        xlims[1] <- -25
      }
    }
  }
  

  asp_rat <- diff(ylims)/diff(xlims)
  if(asp_rat > 0.6){
    scl <- asp_rat * 1.6
  } else {
    scl <- 1
  }
  
  if(asp_rat > 1){
    tscl <- asp_rat
  } else{
    tscl <- 1
  }
  
  #position legend in bottom left if specific coords not given
  if(is.null(legend_position)){
    pos <- c(0.03, 0.03)
  } else {
    pos <- legend_position
  }

  #legend title
  if(!is.null(title)){
    legend_title <- title
  } else{
    legend_title <- paste0("Suitability of ", dz)
    if (monthly){
      legend_title <- paste0(legend_title, " - ", substr(month.name[as.numeric(mo)], 1, 3))
    }
  }

  if(is.na(legend_title)){
    legend_title <- NULL
  }
  
  
  if(!is.null(legend_title)){
    legend_title <- paste0(wrap_text(legend_title, 17), '\n')
  }
  
  
  if (is.null(color_scheme)){
    color_scheme <- c("#00681d", "#fbfcc9", "#9b0065")
  }
  
  if(is.null(bg_pts) & is.null(occ_pts) & is.null(good_pts) & is.null(bad_pts)){
    col_bar_w <- 1.1 * scl
    col_bar_l <- 3.5 * scl
  } else {
    col_bar_w <- 1 * scl
    col_bar_l <- 2.7 * scl
  }
  
  if(!is.null(out_dir)){
    if(grepl('04_rift_valley_fever', out_dir)){
      col_bar_w <- 1.8 
      col_bar_l <- 4.5
    }
    
  }
  
  if(!is.null(transform)){
    fun <- function(x){
      eval(parse(text = paste0(transform, '(x)')))
    }
    
    raster_df$layer <- fun(raster_df$layer)
    
  }
  
  if(is.null(lims)){
    lims <- c(min(raster_df$layer), max(raster_df$layer))
  }
  
  if(is.infinite((lims[1]))){
    lims[1] <- 0
  }
  
  if(is.null(labs)){
    labs <- c('Low', 'High')
  }
  
  if(is.null(breaks)){
    breaks <- lims
  }
  
  
  print(toString(co_lines[1]))
  print(as.numeric(co_lines[2]))
  
  

  #plot map
  if(!is.null(na_color)){
    raster_map <- ggplot()+
      geom_polygon(data=bg_polygon_df,aes(x=long, y=lat, group=group), fill="lightgrey", 
                   color=toString(co_lines[1]), size = as.numeric(co_lines[2]))+
      geom_polygon(data = disputed_df, aes(x = long,y = lat, group = group), linetype = 'dashed', fill = 'transparent', size = 0.25, color = 'gray20') +
      geom_tile(data=raster_df, aes(x=x, y=y, fill=layer, linetype = 'No data'), color = 'transparent',
                show.legend = show_leg)
    
  } else {
    raster_map <- ggplot()+
      geom_polygon(data=bg_polygon_df,aes(x=long, y=lat, group=group), fill="lightgrey", 
                   color=toString(co_lines[1]), size = as.numeric(co_lines[2]))+
      geom_tile(data=raster_df, aes(x=x, y=y, fill=layer), color = 'transparent', show.legend = show_leg)
    
  }
  
 
  if(clip){
    cidf <- na.omit(as.data.frame(clip_inverse, xy = TRUE))
    names(cidf)[3] <- 'layer'
    raster_map <- raster_map +
      geom_raster(data = cidf, aes(x = x, y = y, fill = layer), 
                alpha = 0.35, show.legend = show_leg)
  }
  
  if(!is.null(raster2)){
    
    
    raster2_df <- na.omit(as.data.frame(raster2, xy = TRUE))
    names(raster2_df)[3] <- 'layer'
    
    blp <- ifelse(grepl('henipa', out_dir), tolower('barren land,\nlow population'), 'barren land,\nlow population')
    
    raster_map  <- raster_map  +
      geom_raster(data=raster2_df, aes(x=x, y=y, alpha = 'Barren land,\nlow population'), 
                fill = 'gray65', show.legend = show_leg, color = 'transparent') + 
      scale_alpha_manual(name = NULL, values = c('Barren land,\nlow population' = .6), guide = guide_legend(override.aes = list(alpha = .8, color = 'transparent', size = 0), order = 2))
  }
  
  if(clip){
    if(class(buffer_smooth) == 'SpatialPolygonsDataFrame'){
      buff_df <- fortify(buffer_smooth)
      raster_map <- raster_map +
        geom_polygon(data = buff_df, aes(x = long, y = lat, group = group, size = 'Buffer boundary'), color = 'black', fill = 'transparent', show.legend = show_leg) +
        scale_size_manual(name = NULL, values = c('Buffer boundary' = 0.15), guide = guide_legend(override.aes = list(fill = 'transparent', shape = NA)))
    } 
    
  }
  
  if (!is.null(poly_layer)){
    poly_df <- fortify(poly_layer)
    raster_map <- raster_map + 
      geom_polygon(data=poly_df,aes(x=long, y=lat, group=group), fill=NA, color="gray20", size=0.25)
  }
  
  raster_map <- raster_map + 
    geom_polygon(data=bg_polygon_df,aes(x=long, y=lat, group=group), 
                 fill='transparent', color=toString(co_lines[1]), size = as.numeric(co_lines[2]))+
    geom_polygon(data = disputed_df, aes(x = long,y = lat, group = group), linetype = 'dashed', fill = 'transparent', size = 0.25, color = 'gray20') +
    geom_polygon(data=lakes_df,aes(x=long, y=lat, group=group), fill="white")
    
    if(!is.null(letter)){
      
      ldf <- cbind.data.frame(split(letter, rep(1:4, times=length(letter)/4)), stringsAsFactors=F)
      
      raster_map <- raster_map + 
        geom_text(data = ldf, aes(x = as.numeric(ldf[1,2]), 
                                  y = as.numeric(ldf[1,3]), 
                                  label = as.character(ldf[1,1])), size = 11, hjust = ldf[1, 4])
    }
    raster_map <- raster_map + 
    theme(panel.background = element_rect(fill="white"),
          panel.border = element_rect(color="black", fill=NA, size=rel(1)),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 18),
          legend.box.background = element_rect(color='black', fill="white", size=rel(0.75)),
          legend.background=element_blank(),
          legend.position = pos,
          legend.justification = c(0,0),
          # legend.key.height = unit(.75,"cm"),
          legend.key.width = unit(col_bar_w/2, "cm"),
          legend.title=element_text(size=13 * tscl, hjust = 0.5),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size=11 * tscl),
          legend.spacing.y = unit(0, "cm")) 
    if(!is.null(resc) & !is.null(na_color)){
      raster_map <- raster_map +
        scale_fill_gradientn(colours=color_scheme,
                             name=legend_title,
                             limits=lims, breaks=breaks, labels=labs,
                             na.value = na_color,
                             values = rescale(resc),
                             guide=guide_colourbar(barwidth = col_bar_w, 
                                                   barheight = col_bar_l, 
                                                   title.hjust = 0.5, title.vjust = 0.75, order = 1))
      
    } else if(is.null(resc) & !is.null(na_color)){
      raster_map <- raster_map +
        scale_fill_gradientn(colours=color_scheme,
                             name=legend_title,
                             limits=lims, breaks=breaks, labels=labs,
                             na.value = na_color,
                             guide=guide_colourbar(barwidth = col_bar_w, 
                                                   barheight = col_bar_l, 
                                                   title.hjust = 0.5, title.vjust = 0.75, order = 1))
    } else{
      raster_map <- raster_map +
        scale_fill_gradientn(colours=color_scheme,
                             name=legend_title,
                             limits=lims, breaks=breaks, labels=labs,
                             guide=guide_colourbar(barwidth = col_bar_w, 
                                                   barheight = col_bar_l, 
                                                   title.hjust = 0.5, title.vjust = 0.75, order = 1))
    }
       
    if(!is.null(na_color)){
      
      raster_map <- raster_map +
        scale_linetype_manual(values = c('No data' = 'blank'), name = NULL) +              
        guides(linetype=guide_legend(override.aes=list(fill = na_color), order = 2))
    }
    
    
  if(!is.null(plot_title)){
    raster_map <- raster_map + 
      ggtitle(plot_title)
  }
    
  if (expand) {
    raster_map <- raster_map + coord_quickmap(xlim=xlims, ylim=ylims)
  } else{
    raster_map <- raster_map + coord_quickmap(xlim=xlims, ylim=ylims, expand = FALSE)
  }

  


  if (!is.null(bg_pts)){
    raster_map <- raster_map +
      geom_point(data=bg_pts, aes(x=long, y=lat, colour="Background\npoints"), alpha=0.8, size=.01, show.legend = show_leg) +
      scale_color_manual(name=NULL , values = c("Background\npoints" = 'red'))+
      guides(fill = guide_colorbar(order=1), colour=guide_legend(order=0, override.aes = list(size=.01, fill = 'transparent')))
  }
  if (!is.null(occ_pts)){
    raster_map <- raster_map +  
      geom_point(data=occ_pts, aes(x=long, y=lat, colour="Occurrence\npoints"), alpha=0.8, size=.07, show.legend = show_leg) +
      scale_color_manual(name=NULL , values = c( "Occurrence\npoints" = 'blue'))+
      guides(fill = guide_colorbar(order=1), colour=guide_legend(order=0,override.aes = list(size=.07, fill = 'transparent')))
  }

  if(!is.null(good_pts)){
    raster_map <- raster_map +  
      geom_point(data=good_pts, aes(x=long, y=lat, colour="Good fit"), alpha=0.8, size=.01, show.legend = show_leg)
  }
  if(!is.null(bad_pts)){
    raster_map <- raster_map +  
      geom_point(data=bad_pts, aes(x=long, y=lat, colour="Poor fit"), alpha=0.8,  size=.01, show.legend = show_leg)+
      scale_color_manual(name=wrap_text("Occurrences", 15) , values = c("Poor fit" = "red", "Good fit" = "blue"))+
      guides(fill = guide_colorbar(order=1), colour=guide_legend(order=0,  override.aes = list(size=.01, fill = 'transparent'), title.hjust = 0.5))
  }
    
  

  if (monthly){
    file_path <- paste0(out_dir, "/prediction_map_", mo, ".", device)
  } else {
    file_path <- paste0(out_dir, "/prediction_map.", device)
  }
  if (!is.null(file_name)){
    file_path <- sub("/prediction_map", file_name, file_path)
  }
  
  if(is.null(out_dir)){
    return(raster_map)
  } else{
    message('saving...')
    Cairo(file = file_path, dpi=350, width = size, 
           height = size * asp_rat, unit = 'in', type = device)
    plot(raster_map)
    dev.off()
    
    message('saving low res...')
    file_path <- sub(paste0('.', device), paste0('_low_res.', device), file_path)
    Cairo(file = file_path, dpi=100, width = size,
              height = size * asp_rat, unit = 'in', type = device)
    plot(raster_map)
    dev.off()
    
  }
  
  
  
}

### BINARY RASTER MAP FUNCTION ARGUMENTS ###
#raster_layer           should specify a particular layer such as preds_sry[["mean"]]
#levels                 the order in which you want the categorical variable displayed in the legend c('1', '0')
#labels                 the labels you want to be associated with these levels c('High', 'Low')
#out_dir                example:'FILEPATH' should *not* end with /
#color_scheme           should be passed a vector of 2 HEX color codes with designated values i.e. (1 = #234234, 0 = #123190)
#bg_pts                 any point data with 'lat' and 'long' columns will work ideally dat0.pts
#occ_pts                any point data with 'lat' and 'long' columns will work ideally dat1.pts
#good_pts               any point data with 'lat' and 'long' columns will work ideally dat1.good
#bad_pts                any point data with 'lat' and 'long' columns will work ideally dat1.bad
##NOTE: bad_pts and good_pts are designed to be used together and the legend for good_pts will not appear if bad_pts not specified ##
#crop                   default is set to TRUE, will crop the raster to the extent of the buffer for the specified experiment
#clip                   default is set to TRUE, will clip the raster to the shape of the buffer for the specified experiment
##NOTE: if clip is used crop should also be used to speed up the time it takes to run the clip function ##
#extents                default is set to the extent of the raster layer after crop/clip occurs, must be entered as an extent(object you want to match extents) object
#legend_position        default is set to c(0.03, 0.03) which is the bottom left corner of the map, format of legend_position should either be c(x_coord, y_coord) or 'top', 'left', 'bottom', 'right'
#file_name              default set to "/prediction_map" format should be "/file_name" this will automatically add to whatever out_dir is set to
#poly_layer             can receive a SpatialPolygonsDataFrame such as those that are produced by readOGR which will add the polygon borders as a final layer to the map
#title                  default set to "Suitability of dz", format "My Map Has a Custom Title"
#lims                   a numeric vector of length 2 setting the limits of values that will be shown, if NULL will be 0-1 ex. c(0,500)
#device                 graphics devie for printing image, default is 'jpg', can be anything used by Cairo
#letter                 a vector of length 3 containing the text, x and y coordinates in that order ex. c('A', 45.2, 57.6, 0) and hjust


#function to create map
get_binary_raster_map <- function(raster_layer, levels, labels, out_dir, na_color = 'lightgrey',
                                  color_scheme = NULL, bg_pts=NULL, occ_pts=NULL, 
                                  good_pts=NULL, bad_pts=NULL, crop = TRUE, clip=TRUE, 
                                  extents = NULL, expand=TRUE,  legend_position = NULL, 
                                  file_name=NULL, poly_layer=NULL, title=NULL, lims = NULL,
                                  raster2 = NULL, plot_title = NULL, lkh = 0.65, lkw = 0.65,
                                  lw = 17, size = 10, device = 'jpg', co_lines = c('gray20', 0.25),
                                  letter = NULL){
  
  
  #crop raster if crop=TRUE
  if(crop){
    #crop preds raster to extent of buffer
    if(!is.null(extents)){
      cropped_raster <- crop(raster_layer, extents)
    } else if(is.null(extents) & exists('buffer_smooth')){
      cropped_raster <- crop(raster_layer, extent(buffer_smooth))
    } else{
      stop('Nothing to crop raster layer by, please add extents or buffer_smooth sp object')
    }
    
    if(!is.null(raster2)){
      
      if(!is.null(extents)){
        raster2 <- crop(raster2, extents)
      } else if(is.null(extents) & exists('buffer_smooth')){
        raster2 <- crop(raster2, extent(buffer_smooth))
      }
    }
    
    raster_layer <- cropped_raster
  }
  
  
  #clip raster if clip =TRUE
  if(clip){
    #clip the raster layer to buffer
    clipped_raster <- raster::mask(cropped_raster, buffer_smooth)
    raster_layer <- clipped_raster
  }
  
  
  #convert raster to data.frame and omit NAs
  raster_df <- na.omit(as.data.frame(raster_layer, xy=TRUE))
  
  #find map extent based on preds if extents aren't given
  if(is.null(extents)){
    extent_ref <- extent(raster_layer)
  } else {
    extent_ref <- extents
  }
  xlims <- c(extent_ref@xmin, extent_ref@xmax)
  ylims <- c(extent_ref@ymin, extent_ref@ymax)
  
  
  if(!is.null(out_dir)){
    if(grepl('04_rift_valley_fever', out_dir)){
      if(min(ylims) < -38){
        ylims[1] <- -38
      }
      if(min(xlims) < -25){
        xlims[1] <- -25
      }
    }
  }
  
  
  asp_rat <- diff(ylims)/diff(xlims)
  if(asp_rat > 0.6){
    scl <- asp_rat * 1.6
  } else {
    scl <- 1
  }
  
  if(asp_rat > 1){
    tscl <- asp_rat
  } else{
    tscl <- 1
  }
  
  #position legend in bottom left if specific coords not given
  if(is.null(legend_position)){
    pos <- c(0.03, 0.03)
  } else {
    pos <- legend_position
  }
  
  #legend title
  if(!is.null(title)){
    legend_title <- title
  } else{
    legend_title <- paste0("Suitability of ", dz)
    if (monthly){
      legend_title <- paste0(legend_title, " - ", substr(month.name[as.numeric(mo)], 1, 3))
    }
  }
  
  if(is.na(legend_title)){
    legend_title <- NULL
  }
  
  
  
  if(!is.null(legend_title)){
    legend_title <- paste0(wrap_text(legend_title, lw), '\n')
  }
  
  
  if (is.null(color_scheme)){
    color_scheme <- c('0' = "#CBE0CF", '1' = "#9b0065")
  }
  
  if(is.null(bg_pts) & is.null(occ_pts) & is.null(good_pts) & is.null(bad_pts)){
    col_bar_w <- 1.1 * scl
    col_bar_l <- 3.5 * scl
  } else {
    col_bar_w <- 1 * scl
    col_bar_l <- 2.7 * scl
  }
  
  if(!is.null(out_dir)){
    if(grepl('04_rift_valley_fever', out_dir)){
      col_bar_w <- 1.8 
      col_bar_l <- 4.5
    }
  }
  
  
  if(is.null(lims)){
    lims <- c(0,1)
  }
  
  if(is.null(levels)){
    levels <- unique(as.character(raster_df[,3]))
  }
  
  if(is.null(labels)){
    labels <- levels
  }
  
  if(sum(is.na(values(raster_layer))) > 0){
    labels <- append(labels, 'No data')
  }
  
  #plot map
  bin_map <- ggplot()+
    geom_polygon(data=bg_polygon_df,aes(x=long, y=lat, group=group), fill="light grey", color=toString(co_lines[1]), size = as.numeric(co_lines[2]))+
    geom_polygon(data = disputed_df, aes(x = long,y = lat, group = group), linetype = 'dashed', fill = 'transparent', size = 0.25, color = 'gray20') +
    geom_tile(data=raster_df, aes(x=x, y=y, fill=factor(as.character(raster_df[,3]), levels = levels)), color = 'transparent')
    
    if(!is.null(raster2)){
      
      if(clip){
        raster2 <- raster::mask(raster2, buffer_smooth)
      }
      
      raster2_df <- na.omit(as.data.frame(raster2, xy = TRUE))
      names(raster2_df)[3] <- 'layer'
      
      blp <- ifelse(grepl('henipa', out_dir), tolower('Barren low\npopulation'), 'Barren low\npopulation')
      
      bin_map <- bin_map +
        geom_tile(data=raster2_df, aes(x=x, y=y, alpha = blp), color = 'transparent', fill = 'gray65') + 
        scale_alpha_manual(name = NULL, values = c(blp = .8), guide = guide_legend(override.aes = list(alpha = .8), order = 2))
    }
  
  if (!is.null(poly_layer)){
    poly_df <- fortify(poly_layer)
    bin_map <- bin_map + geom_polygon(data=poly_df,aes(x=long, y=lat, group=group), fill=NA, color="gray20", size=0.25)
  }
  
  if(!is.null(letter)){
    
    ldf <- cbind.data.frame(split(letter, rep(1:4, times=length(letter)/4)), stringsAsFactors=F)
    
    bin_map <- bin_map + 
      geom_text(data = ldf, aes(x = as.numeric(ldf[1,2]), 
                                y = as.numeric(ldf[1,3]), 
                                label = as.character(ldf[1,1])), size = 11, hjust = ldf[1, 4])
  }
  
  bin_map <- bin_map +
    geom_polygon(data=bg_polygon_df,aes(x=long, y=lat, group=group), fill=NA, color=toString(co_lines[1]), size=as.numeric(co_lines[2]))+
    geom_polygon(data = disputed_df, aes(x = long,y = lat, group = group), linetype = 'dashed', fill = 'transparent', size = 0.3, color = 'gray20') +
    geom_polygon(data=lakes_df,aes(x=long, y=lat, group=group), fill="white")+
    theme(panel.background = element_rect(fill="white"),
          panel.border = element_rect(color="black", fill=NA, size=rel(1)),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          legend.box.background = element_rect(color='black', fill="white", size=rel(0.75)),
          legend.background=element_blank(),
          legend.position = pos,
          legend.justification = c(0,0),
          plot.title = element_text(size = 18),
          legend.key.height = unit(lkh,"cm"),
          legend.key.width = unit(lkw, "cm"),
          legend.title=element_text(legend_title, size=13 * tscl, hjust = 0.5),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size=11 * tscl),
          legend.spacing.y = unit(0, "cm"))+
    guides(fill = guide_legend(title.hjust = 0.5, title.vjust = 0.75, order = 1)) +
    scale_fill_manual(name = legend_title, values=color_scheme, na.value = na_color,
                      labels = labels, drop = FALSE)
  if (expand) {
    bin_map <- bin_map + coord_quickmap(xlim=xlims, ylim=ylims)
  } else{
    bin_map <- bin_map + coord_quickmap(xlim=xlims, ylim=ylims, expand = FALSE)
  }
  
  
  
  if(!is.null(plot_title)){
    bin_map <- bin_map +
      ggtitle(plot_title)
  }
  
  if (!is.null(bg_pts)){
    bin_map <- bin_map +
      geom_point(data=bg_pts, aes(x=long, y=lat, colour="Background\npoints"), alpha=0.8, size=rel(.2)) +
      scale_color_manual(name=NULL , values = c("Background\npoints" = 'red'))+
      guides(fill = guide_colorbar(order=1), colour=guide_legend(order=0, override.aes = list(size=rel(.7))))
  }
  if (!is.null(occ_pts)){
    bin_map <- bin_map +  geom_point(data=occ_pts, aes(x=long, y=lat, colour="Occurrence\npoints"), alpha=0.8, size=rel(.9)) +
      scale_color_manual(name=NULL , values = c( "Occurrence\npoints" = 'blue'))+
      guides(fill = guide_colorbar(order=1), colour=guide_legend(order=0,override.aes = list(size=rel(.7))))
  }
  
  if(!is.null(good_pts)){
    bin_map <- bin_map +  geom_point(data=good_pts, aes(x=long, y=lat, colour="Good fit"), alpha=0.8, size=rel(.2))
  }
  if(!is.null(bad_pts)){
    bin_map <- bin_map +  geom_point(data=bad_pts, aes(x=long, y=lat, colour="Poor fit"), alpha=0.8,  size=rel(.2))+
      scale_color_manual(name=wrap_text("Occurrence Data", 15) , values = c("Poor fit" = "red", "Good fit" = "blue"))+
      guides(fill = guide_colorbar(order=1), colour=guide_legend(order=0,  override.aes = list(size=rel(.7)), title.hjust = 0.5))
  }
  
  if (monthly){
    file_path <- paste0(out_dir, "/binary_prediction_map_", mo, ".", device)
  } else {
    file_path <- paste0(out_dir, "/binary_prediction_map.", device)
  }
  if (!is.null(file_name)){
    file_path <- sub("/binary_prediction_map", file_name, file_path)
  }
  
  if(is.null(out_dir)){
    return(bin_map)
  } else{
    message('saving...')
    Cairo(file = file_path, dpi=350, width = size, 
          height = size * asp_rat, unit = 'in', type = device)
    plot(bin_map)
    dev.off()
    
    message('saving low res...')
    file_path <- sub(paste0('.', device), paste0('_low_res.', device), file_path)
    Cairo(file = file_path, dpi=100, width = size,
          height = size * asp_rat, unit = 'in', type = device)
    plot(bin_map)
    dev.off()
    
  }
  
  
  
}

### RASTER MAP GRID FUNCTION ARGUMENTS ###
#cov_stack should be a stack of raster layers that correspond to each covariate by name e.g. covs
#order should be a vector of the numerical order of covariates e.g. c(4, 3, 1, 6, 5, 2)
#out_dir example:'FILEPATH' should *not* end with /
#bg_pts any point data with 'lat' and 'long' columns will work ideally dat0.pts
#occ_pts any point data with 'lat' and 'long' columns will work ideally dat1.pts
#legend_position default is set to c(0.03, 0.03) which is the bottom left corner of the map, format of legend_position should either be c(x_coord, y_coord) or 'top', 'left', 'bottom', 'right'
#file_name default set to "/FILEPATH" format should be "/file_name" this will automatically add to whatever out_dir is set to
#device the graphics device used for printing the image - anything used by Cairo is good, default is 'jpg'


#get_raster_grid function
get_raster_grid <- function(cov_stack, 
                            order, 
                            out_dir, 
                            bg_pts=NULL, 
                            occ_pts=NULL, 
                            extents=NULL, 
                            expand=TRUE, 
                            legend_position = NULL, 
                            title = NULL, 
                            file_name=NULL,
                            device = 'jpg'){

  #create map_list
  map_list <- list()
  ordered_covs <- cov_names[order]

  #for loop to manipulate and plot cov data
  for (name in ordered_covs){
    print(name)
    #take data by covariate in order of relinf and get relevant names bat_mean_frequency is a special case
    if (name == "bat_frequency_mean"){
      cov_name_long <- "Mean bat frequency"
    }  else if (name == 'excess_sum2'){
      
        cov_name_long <- "Rainfall (sum prev. 2 months)"
    } else{
      if(!is.null(cov_name_list$title)){
        cov_name_long <- cov_name_list$cov_name_long[grep(name, cov_name_list$title)]
      } else {
        cov_name_long <- cov_name_list$cov_name_long[grep(name, cov_name_list$cov_name)]
      }#match short cov names to long names for titling
    }


    #subset covs by name
    cov_dat <- covs[[name]]

    #crop extent of cov_dat by the extent of the buffer
    if(is.null(extents)){
      extent_ref <- extent(buffer_smooth)
    } else {
      extent_ref <- extents
    }
    
    #find map extent based on extents of cov if extents aren't given
    
    xlims <- c(extent_ref@xmin, extent_ref@xmax)
    ylims <- c(extent_ref@ymin, extent_ref@ymax)

    if(grepl('04_rift_valley_fever', out_dir)){
      if(min(ylims) < -38){
        ylims[1] <- -38
      }
      if(min(xlims) < -25){
        xlims[1] <- -25
      }
    }
    
    asp_rat <- diff(ylims)/diff(xlims)
    if(asp_rat > 0.6){
      scl <- asp_rat * 1.4
    } else {
      scl <- 1
    }
    
    if(asp_rat > 1){
      tscl <- asp_rat
    } else{
      tscl <- 1
    }
    
    if (expand){
      dummy_df <- data.frame(x=numeric(1), y=numeric(1))
      dummy_gg <- ggplot(dummy_df) + geom_point(aes(x=x, y=y)) + coord_quickmap(xlims, ylims)
      
      
      if (packageVersion("ggplot2")=='2.2.1'){
        xrange <- ggplot_build(dummy_gg)$layout$panel_ranges[[1]]$x.range
        yrange <- ggplot_build(dummy_gg)$layout$panel_ranges[[1]]$y.range
      } else if (packageVersion("ggplot2") == '3.0.0'){
        xrange <- ggplot_build(dummy_gg)$layout$panel_params[[1]]$x.range
        yrange <- ggplot_build(dummy_gg)$layout$panel_params[[1]]$y.range
      } else if (packageVersion("ggplot2") == '3.1.0'){
        xrange <- ggplot_build(dummy_gg)$layout$panel_params[[1]]$x.range
        yrange <- ggplot_build(dummy_gg)$layout$panel_params[[1]]$y.range
      } else if (packageVersion("ggplot2") == '3.2.1'){
        xrange <- ggplot_build(dummy_gg)$layout$panel_params[[1]]$x.range
        yrange <- ggplot_build(dummy_gg)$layout$panel_params[[1]]$y.range
      }
      
      
      
      ext_crop <- extent_ref
      
      ext_crop@xmin <- xrange[1]
      ext_crop@xmax <- xrange[2]
      ext_crop@ymin <- yrange[1]
      ext_crop@ymax <- yrange[2]
    } else {
      ext_crop <- extent_ref
    }
    
    cov_dat <- crop(cov_dat, ext_crop)
    
    #convert raster to data.frame and omit NAs
    cov_df <- na.omit(as.data.frame(cov_dat, xy=TRUE))
    names(cov_df) <- c("x", "y", "covar")

    #convert special covariates
    if (name == "crutstmp"){
      unit_lab <- "Temperature (°C)"
      cov_df$covar <- (cov_df$covar/10)
    } else if (name == "crutsard"){
      unit_lab <- "mm/mm"
    } else if (name== "bat_frequency_mean"){
      unit_lab <- "Suitability / all species"
    }  else{
      if(!is.null(cov_name_list$title)){
        unit_lab <- capitalize(cov_name_list$units[grep(name, cov_name_list$title)])
      } else {
        unit_lab <- capitalize(cov_name_list$units[grep(name, cov_name_list$cov_name)])
      }
    }
    
    
    if(length(unit_lab > 1)){
      unit_lab <- unit_lab[1]
    }
    
    if(unit_lab == '10degC'){
      unit_lab <- "Temperature (°C)"
      cov_df$covar <- (cov_df$covar/10)
    }
    
    if(unit_lab == '100days'){
      unit_lab <- "Days"
      cov_df$covar <- (cov_df$covar/100)
    }
    
    if(unit_lab=="Mm" | unit_lab=="Mm/day" | unit_lab == 'Km' | unit_lab == 'M'){
      unit_lab <- tolower(unit_lab)
    }

  
    


    #position legend in bottom left if specific coords not given
    if(is.null(legend_position)){
      pos <- c(0.02, 0.03)
    } else {
      pos <- legend_position
    }
  
   
    col_bar_w <- 1.1
    col_bar_l <- 3.5
    
    r <- range(na.omit(cov_df$covar))

    bks <- seq(r[1], r[2], diff(r)/2)
    
    floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
    
    if(max(bks) < 10){
      bks <- floor_dec(bks, 2)
    } else {
      bks <- floor(bks)
    }
    
    #plot map
    cov_map <- ggplot()+
      geom_tile(data=cov_df, aes(x=x, y=y, fill=covar), color = 'transparent')+
      geom_polygon(data=bg_polygon_df, aes(x=long, y=lat, group=group), fill=NA, color="gray20", size = 0.25)+
      geom_polygon(data = disputed_df, aes(x = long, y = lat, group = group), linetype = 'dashed', fill = 'transparent', size = 0.25, color = 'gray20') +
      scale_fill_distiller(palette = "RdYlBu",
                           name=paste0(wrap_text(unit_lab, 11), '\n'), breaks = bks, limits = c(min(bks), max(bks)),
                           guide=guide_colourbar(barwidth = col_bar_w * scl,
                                                 barheight = col_bar_l * scl,
                                                 title.hjust = 0.5,
                                                 title.vjust=0.75))+
      ggtitle(cov_name_long)+
      theme(plot.title = element_text(size=16, colour="black", hjust = 0),
            panel.background = element_rect(fill="white"),
            panel.border = element_rect(color="black", fill=NA, size=rel(1)),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid = element_blank(),
            legend.box.background = element_rect(color='black', fill="white", size=rel(.6)),
            legend.background=element_blank(),
            legend.position = pos,
            legend.key = element_rect(colour = NA, fill = NA),
            legend.justification = c(0,0),
            legend.title=element_text(size = 11),
            legend.text = element_text(size = 10),
            plot.margin = unit(c(6, 6, 6, 6), "points"))
    
    if(expand){
      cov_map <- cov_map + coord_quickmap(xlim=xlims, ylim=ylims)
    } else {
      cov_map <- cov_map + coord_quickmap(xlim=xlims, ylim=ylims, expand=FALSE)
    }

    #add background points if bg_pts TRUE
    if (!is.null(bg_pts)){
      cov_map <- cov_map + geom_point(data=bg_pts, aes(x=long, y=lat, colour="Background data"), alpha=0.4, size=rel(.01)) +
        scale_color_manual(values = c('Background data' = 'red'))
    }

    #add occurrence points if occ_pts TRUE
    if (!is.null(occ_pts)){
      cov_map <- cov_map + + geom_point(data=occ_pts, aes(x=long, y=lat, colour="Occurrence data"), alpha=0.4, size=rel(.01)) +
        scale_color_manual(values = c('Occurrence data' = 'blue'))
    }

    map_list[[name]] <- cov_map
  }

  ### START creating raster grid ###
  #title
  if(!is.na(title)){
    if(is.null(title)){
      title <- paste0("Maps of Covariate Values for ", dz)
      
    } 

    title <- textGrob(wrap_text(title, 55), gp=gpar(fontface="bold", size=24))
  }
  
  
  #plot grid
  raster_grid <- marrangeGrob(grobs=map_list, top=title, bottom = quote(textGrob(paste("page", g, "of", npages))), nrow=2, ncol=2)

  #save grid page by page if more than one page
  if(!is.null(out_dir)){
    if (length(raster_grid)>1){
    for (i in 1:length(raster_grid)){
      if (monthly){
        file_path <- paste0(out_dir, "/cov_raster_grid_", mo, "_page_", i, ".", device)
      } else {
        file_path <- paste0(out_dir, "/cov_raster_grid_page_", i, ".", device)
      }
      ggsave(raster_grid[[i]], filename = file_path, dpi=350, width = 16, height = 16 * asp_rat)
      file_pathlr <- sub(paste0('.', device), paste0('_low_res.', device), file_path)
      ggsave(raster_grid[[i]], filename = file_pathlr, dpi=100, width = 16, height = 16 * asp_rat)
    }
  } else{
    if (monthly){
      file_path <- paste0(out_dir, "/cov_raster_grid_", mo, ".", device)
    } else {
      file_path <- paste0(out_dir, "/cov_raster_grid.", device)
    }
    if (!is.null(file_name)){
      file_path <- sub("/cov_raster_grid", file_name, file_path)
    }
    ggsave(raster_grid, filename = file_path,  dpi=350, width = 16 , height = (16 * asp_rat))
    file_path <- sub(paste0('.', device), paste0('_low_res.', device), file_path)
    ggsave(raster_grid, filename = file_path,  dpi=100, width = 16 , height = (16 * asp_rat))
    }
  }

  if(is.null(out_dir)){
    return(raster_grid)
  }
  
}


### BUFFER MAP FUNCTION ####

#Args
#color             hex code or color name to color the overlay shapefile
#overlay           shapefile which will be displayed over a background shapefile
#pts               data frame of points to be shown on the map
#extents           map extents to be shown, if none given will snap to overlay
#legend            legend T/F
#legend_position   coordinates within 0-1 plane for bottom left corner of legend
#filename          name of file to be saved
#title             title of plot
#device            graphics device for printing image, default is 'jpg' anything used by Cairo is good


get_buffer_map <- function(color = NULL, 
                           overlay = NULL, 
                           pts=NULL,  
                           extents = NULL, 
                           legend = FALSE, 
                           legend_position = NULL, 
                           file_name=NULL, 
                           title=NULL, 
                           out_dir = NULL,
                           device = 'jpg'){
  
  #find map extent based on preds if extents aren't given
  if(is.null(extents)){
    extent_ref <- extent(overlay)
  } else {
    extent_ref <- extents
  }
  xlims <- c(extent_ref@xmin, extent_ref@xmax)
  ylims <- c(extent_ref@ymin, extent_ref@ymax)
  
  asp_rat <- diff(ylims)/diff(xlims)
  
  #position legend in bottom left if specific coords not given
  if(is.null(legend_position)){
    pos <- c(0.03, 0.03)
  } else {
    pos <- legend_position
  }
  
  if(is.null(color)){
    color <- 'lightgrey'
  }
  
  if(is.null(overlay)){
    stop('no overlay provided, please provide overlay shapefile')
  } else if(class(overlay) == 'SpatialPolygonsDataFrame') {
    poly_df_1 <- fortify(overlay)
    buffer_map <- ggplot()+
      geom_polygon(data=bg_polygon_df,aes(x=long, y=lat, group=group), color="gray20", size = 0.25)+
      geom_polygon(data = disputed_df, aes(x = long,y = lat, group = group), linetype = 'dashed', fill = 'transparent', size = 0.25, color = 'gray20') +
      geom_polygon(data = poly_df_1, aes(x=long, y=lat, group=group),
                   fill = color, color = 'gray20', alpha = 0.6, size = 0.3)
    
  } else if(class(overlay) == 'RasterLayer'){
    poly_df_1 <- na.omit(as.data.frame(overlay, xy = TRUE))
    
    buffer_map <- ggplot()+
      geom_polygon(data=bg_polygon_df,aes(x=long, y=lat, group=group), color="gray20", size=0.3)+
      geom_polygon(data = disputed_df, aes(x = long,y = lat, group = group), linetype = 'dashed', fill = 'transparent', fill = 'transparent', size = 0.3, color = 'gray20') +
      geom_tile(data = poly_df_1, aes(x=x, y=y), color = 'transparent',
                   fill = color, alpha = 0.6)
  }
  
  #plot map
  buffer_map <- buffer_map + 
    geom_polygon(data=lakes_df,aes(x=long, y=lat, group=group), fill="white") +
    theme(panel.background = element_rect(fill="white"),
          panel.border = element_rect(color="black", fill=NA, size=rel(1)),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          legend.box.background = element_rect(color='black', fill="white", size=rel(0.75)),
          legend.background=element_blank(),
          legend.position = pos,
          legend.justification = c(0, 0),
          legend.key.height = unit(.25,"cm"),
          legend.key.width = unit(.4, "cm"),
          legend.title=element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size=rel(.55)),
          legend.spacing.y = unit(0, "cm")) + 
      coord_quickmap(xlim=xlims, ylim=ylims, expand = FALSE)

  
  if (!is.null(pts)){
    buffer_map <- buffer_map +
      geom_point(data=pts, aes(x=long, y=lat, colour="Background\npoints"), alpha=0.8, size=rel(.9)) +
      scale_color_manual(name=NULL , values = c("Background\npoints" = 'black'))
  }
  
  if(!is.null(title)){
    buffer_map <- buffer_map +
      ggtitle(wrap_text(title, 35))
  }
  
  if(is.null(out_dir)){
    message('no out_dir provided, file will not be saved')
  } else{
    if(is.null(file_name)){
      file_name <- paste0('buffer_map.', device)
    } else {
      file_name <- paste0(file_name, '.', device)
    }
    file_name <- paste0(out_dir, '/', file_name)
    
    ggsave(buffer_map, filename = file_name, width = 8, height = 8 * asp_rat, dpi = 350)
    file_name <- sub(paste0('.', device), paste0('_low_res.', device), file_name)
    ggsave(buffer_map, filename = file_name, width = 8, height = 8 * asp_rat, dpi = 100)
  }
  
  if(is.null(out_dir)){
    return(buffer_map)
  }
  
}

