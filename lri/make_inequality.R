###########################################################################################################
# Prepare absolute and relative deviation from country mean and corresponding bar chart 
# for LRI Stg 2 paper
###########################################################################################################
make_inequality <- function(run_date, #lri run date
                            end_year, #end year for modeling (2000-2019 = end year is 2019)
                            no_data_countries, #string of iso3 of countries to drop from bar chart (no data or other exclusion)
                            out_dir, #director in which to save the figure
                            pdf_name){ #name under which to save the figure
  
# (0) Define function to make mini GBD superregion map ################################  
  make_spr_reg_map <- function(spr_reg_colors){
   
     # Load master shapefile
    master_shape <- readOGR(
      '<<<< FILEPATH REDACTED >>>>'
    )
    # Merge on location hierarchy info
    locs <- fread('<<<< FILEPATH REDACTED >>>>')
    locs <- locs[,.(loc_id, ihme_lc_id)]
    superregion_ids <- fread('<<<< FILEPATH REDACTED >>>>')
    super <- superregion_ids[, .(location_id, super_region_name)]
    locs <- merge(
      locs, super, 
      by.x = c("loc_id"), 
      by.y = c("location_id"), 
      all.x=T
    )
    master <- sp::merge(
      x=master_shape, y=locs, 
      by.x = "ADM0_A3", by.y = "ihme_lc_id", 
      all.x=TRUE
    )
    
    # Make Stage 3 countries gray
    stage <- fread('<<<< FILEPATH REDACTED >>>>')
    stage <- stage[,c("loc_id", "Stage")]
    master <- sp::merge(master, stage, by = "loc_id")
    
    #drop antarctica
    master <- master[master$ADM0_A3 != "ATA", ]
    
    # Formatting
    master$ADM0_A3 <- as.character(master$ADM0_A3)
    master$super_region_name <- as.character(master$super_region_name)
    master$Stage <- as.character(master$Stage)
    master[master$ADM0_A3 == "SDS", "ADM0_A3"] <- "SSD"
    master[master$ADM0_A3 == "SSD", "super_region_name"] <- "Sub-Saharan Africa"
    master[master$ADM0_A3 == "SOL", "super_region_name"] <- "Sub-Saharan Africa"
    master[master$ADM0_A3 == "SAH", "super_region_name"] <- "North Africa and Middle East"
    master[master$ADM0_A3 == "PSX", "super_region_name"] <- "North Africa and Middle East"
    master[master$ADM0_A3 == "GUF", "super_region_name"] <- "Latin America and Caribbean"
    master[master$ADM0_A3 == "GRL", "super_region_name"] <- "High-income"
    master[is.na(master$Stage), 'Stage'] <- "4"
    master[master$Stage == "3", 'super_region_name'] <- "High-income"
    
    #merge on superregion info
    x <- fortify(master, region = "super_region_name")
    y <- fortify(master)
    
    #plot
    minimap <- ggplot() + 
      geom_polygon(data = x, aes(x=long, y = lat, group = group, fill = id))+
      theme_void() +
      theme(
        legend.position = "none",
        panel.background = element_rect(fill='#FFFFFF', colour='#000000')
      ) +
      scale_fill_manual(values = spr_reg_colors)
    
    # Return ggplot object
    return(minimap)
  }
  
  #(1) Prep mortality data ##############################################################################
  locs <- read.csv('<<<< FILEPATH REDACTED >>>>')
  
  lri_0 <- read.csv('<<<< FILEPATH REDACTED >>>>')[,-1]
  lri_1 <- read.csv('<<<< FILEPATH REDACTED >>>>')[,-1]
  lri_2 <- read.csv('<<<< FILEPATH REDACTED >>>>')[,-1]
  
  ULocs <- unique(lri_0$ADM0_NAME)
  Which <- rep(NA, length(ULocs))
  for (i in 1:length(ULocs)){
    tmp <- which(as.character(locs$location_name) == as.character(ULocs[i]))
    if (length(tmp)){
      Which[i] <- tmp
    }
  }
  
  #manually fix mismatches due to special characters
  print('Expecting NAs at 12 15 31 73 81')
  print('NAs at')
  print(which(is.na(Which)))
  
  Which[12] <- 896 #Cote d'Ivoire 
  Which[15] <- 796 #Republic of Congo = GBD Congo
  Which[31] <- 897 #Gambia = GBD The Gambia
  Which[73] <- 583 #Palestina = GBD Palestine
  Which[81] <- 906 #Sao Tome and Principe
  
  print('Now expecting no NAs')
  if (length(which(is.na(Which))) > 0) print('NAs remain. Resolve!')
  if (length(which(is.na(Which))) == 0) print('NAs correctly resolved.')

  lri_0$super_region <- as.character(lri_0$region)
  lri_1$super_region <- as.character(lri_1$region)
  lri_2$super_region <- as.character(lri_2$region)
  
  lri_1$abs_lri <- lri_1$rel_lri <- lri_1$a0_lri <- lri_1$mean
  lri_2$abs_lri <- lri_2$rel_lri <- lri_2$a0_lri <- lri_2$mean
  
  #assign GBD superregion to a0,a1,a2 files
  for (year in 2000:end_year){
    for (i in 1:length(ULocs)){
      tmp_a0 <- which(as.character(lri_0$ADM0_NAME) == as.character(ULocs)[i] & lri_0$year == year)
      lri_0$super_region[tmp_a0]<- as.character(locs$super_region_name)[Which[i]]
      
      tmp_a1 <- which(as.character(lri_1$ADM0_NAME) == as.character(ULocs)[i] & lri_1$year == year)
      lri_1$super_region[tmp_a1]  <- as.character(locs$super_region_name[Which[i]])
      
      tmp_a2 <- which(as.character(lri_2$ADM0_NAME) == as.character(ULocs)[i] & lri_2$year == year)
      lri_2$super_region[tmp_a2]  <- as.character(locs$super_region_name)[Which[i]]
      
      #assign absolute and relative deviations
      if (length(tmp_a1)>1){
        lri_1$a0_lri[tmp_a1] <- lri_0$mean[tmp_a0]
        lri_1$abs_lri[tmp_a1] <- (lri_1$mean[tmp_a1] - lri_0$mean[tmp_a0])
        lri_1$rel_lri[tmp_a1] <- (lri_1$mean[tmp_a1] - lri_0$mean[tmp_a0]) / lri_0$mean[tmp_a0]
        
        lri_2$a0_lri[tmp_a2] <- lri_0$mean[tmp_a0]
        lri_2$abs_lri[tmp_a2] <- (lri_2$mean[tmp_a2] - lri_0$mean[tmp_a0])
        lri_2$rel_lri[tmp_a2] <- (lri_2$mean[tmp_a2] - lri_0$mean[tmp_a0]) / lri_0$mean[tmp_a0]
      }
    }
  }
  
  # get iso3s
  isos <- fread('<<<< FILEPATH REDACTED >>>>')
  isos[, grep('V', names(isos)) := NULL]
  isos[, ADM0_CODE := NULL]
  setnames(isos, 'iso3', 'ISO3')
  lri_0 <- data.table(merge(lri_0, isos[, c('ISO3', 'location_name')], by.x = 'ADM0_NAME', by.y = 'location_name', all.x = T))
  lri_1 <- data.table(merge(lri_1, isos[, c('ISO3', 'location_name')], by.x = 'ADM0_NAME', by.y = 'location_name', all.x = T))
  lri_2 <- data.table(merge(lri_2, isos[, c('ISO3', 'location_name')], by.x = 'ADM0_NAME', by.y = 'location_name', all.x = T))
  
  # clean up manually
  lri_0[ADM0_NAME == "Côte d'Ivoire", ISO3 := 'CIV']
  lri_0[ADM0_NAME == "São Tomé and Príncipe", ISO3 := 'STP']
  lri_1[ADM0_NAME == "Côte d'Ivoire", ISO3 := 'CIV']
  lri_1[ADM0_NAME == "São Tomé and Príncipe", ISO3 := 'STP']
  lri_2[ADM0_NAME == "Côte d'Ivoire", ISO3 := 'CIV']
  lri_2[ADM0_NAME == "São Tomé and Príncipe", ISO3 := 'STP']
  lri_0[ADM0_NAME == "Côte d'Ivoire", super_region := 'Sub-Saharan Africa']
  lri_0[ADM0_NAME == "São Tomé and Príncipe", super_region := 'Sub-Saharan Africa']
  lri_1[ADM0_NAME == "Côte d'Ivoire", super_region := 'Sub-Saharan Africa']
  lri_1[ADM0_NAME == "São Tomé and Príncipe", super_region := 'Sub-Saharan Africa']
  lri_2[ADM0_NAME == "Côte d'Ivoire", super_region := 'Sub-Saharan Africa']
  lri_2[ADM0_NAME == "São Tomé and Príncipe", super_region := 'Sub-Saharan Africa']
  lri_0[ISO3 == "COG", super_region := 'Sub-Saharan Africa']
  lri_1[ISO3 == "COG", super_region := 'Sub-Saharan Africa']
  lri_2[ISO3 == "COG", super_region := 'Sub-Saharan Africa']
  lri_0[ISO3 == "GMB", super_region := 'Sub-Saharan Africa']
  lri_1[ISO3 == "GMB", super_region := 'Sub-Saharan Africa']
  lri_2[ISO3 == "GMB", super_region := 'Sub-Saharan Africa']
  lri_0[ISO3 == "PSE", super_region := 'North Africa and Middle East']
  lri_1[ISO3 == "PSE", super_region := 'North Africa and Middle East']
  lri_2[ISO3 == "PSE", super_region := 'North Africa and Middle East']
  
  # remove countries where we have no data
  lri_0 <- lri_0[!grep(paste0(no_data_countries, collapse = '|'), ISO3),]
  lri_1 <- lri_1[!grep(paste0(no_data_countries, collapse = '|'), ISO3),]
  lri_2 <- lri_2[!grep(paste0(no_data_countries, collapse = '|'), ISO3),]
  
  #remove admins with NA estimates
  lri_0 <- lri_0[!is.na(ISO3) & !is.na(mean)]
  lri_1 <- lri_1[!is.na(ISO3) & !is.na(mean)]
  lri_2 <- lri_2[!is.na(ISO3) & !is.na(mean)]
  
  #compare 2000 and end_year values
  old_a2 <- lri_2[which(lri_2$year == 2000),]
  old_a1 <- lri_1[which(lri_1$year == 2000),]
  old_a0 <- lri_0[which(lri_0$year == 2000),]
  new_a2 <- lri_2[which(lri_2$year == end_year),]
  new_a1 <- lri_1[which(lri_1$year == end_year),]
  new_a0 <- lri_0[which(lri_0$year == end_year),]
  
  panels <- vector("list",3)
  
  old_a0 <- old_a0[order(new_a0$super_region, new_a0$mean),]
  new_a0 <- new_a0[order(new_a0$super_region, new_a0$mean),]
  
  #subset to more than 1 death per 10,000
  keep <- which(new_a0$mean > 0.0001)
  
  old_a0 <- old_a0[keep,]
  new_a0 <- new_a0[keep,]
  a0_order <- new_a0$ADM0_NAME[keep]
  a0_iso <- new_a0$ISO3[keep]
  
  stat_fy <- old_a0[,c("ISO3", "mean","lower","upper","year","year", "cirange","super_region")]
  stat_fy$year <- as.numeric(stat_fy$year)
  names(stat_fy)[2:7] <- c("A0_mean","A2_min","A2_max","A2_relmin","A2_relmax","plot_order")
  stat_fy$plot_order <- 1:nrow(stat_fy)
  
  stat_ly <- new_a0[,c("ISO3", "mean","lower","upper","year","year","cirange","super_region")]
  stat_ly$year <- as.numeric(stat_ly$year)
  names(stat_ly)[2:7] <- c("A0_mean","A2_min","A2_max","A2_relmin","A2_relmax","plot_order")
  stat_ly$plot_order <- 1:nrow(stat_ly)
  
  stat_fy$A2_relmax <- as.numeric(stat_fy$A2_relmax)
  stat_ly$A2_relmax <- as.numeric(stat_ly$A2_relmax)
  
  UISO <- unique(stat_ly$ISO3)
  
  for (i in 1:nrow(stat_ly)){
    tmpold <- old_a2[which(old_a2$ISO3 == UISO[i]),]
    stat_fy[i,3] <- range(tmpold$mean,na.rm=TRUE)[[1]]
    stat_fy[i,4] <- range(tmpold$mean,na.rm=TRUE)[[2]]
    stat_fy[i,5] <- range(tmpold$rel_lri,na.rm=TRUE)[[1]]
    stat_fy[i,6] <- range(tmpold$rel_lri,na.rm=TRUE)[[2]]
    tmpnew <- new_a2[which(new_a2$ISO3 == UISO[i]),]
    stat_ly[i,3] <- range(tmpnew$mean,na.rm=TRUE)[[1]]
    stat_ly[i,4] <- range(tmpnew$mean,na.rm=TRUE)[[2]]
    stat_ly[i,5] <- range(tmpnew$rel_lri,na.rm=TRUE)[[1]]
    stat_ly[i,6] <- range(tmpnew$rel_lri,na.rm=TRUE)[[2]]
  }
  
  spr_reg_colors <- c("#6F4070","#CC503E","#E17C05","#73AF48","#0F8554","#1D6996")
  
  #(2) Mortality plots and csvs ##############################################################################
  
  fig1 <- ggplot(stat_ly,
                 aes(x=reorder(ISO3, plot_order), y=A0_mean*100, ymin=A2_min*100, ymax=A2_max*100,
                     color=super_region)
  ) +
    geom_crossbar(
      data=stat_fy, color='#CCCCCC', size=1.8, width=0.05
    ) +
    geom_point(
      data=stat_fy, color='#000000', alpha=.25, shape=18, size=1.8
    ) +
    # Plot last year data colored by GBD super region
    geom_crossbar(size=1.8, width=0.05,
                  position=position_nudge(x=.3)
    ) +
    geom_point(color='#000000', alpha=.5, shape=18, size=1.8,
               position=position_nudge(x=.3)
    )+
    labs(
      title = paste0('LMIC countries with at least 1 LRI death per 10,000 ranked by childhood LRI mortality rate in ', end_year),
      subtitle = 'Highest and lowest second administrative units shown as error bars\n2000 levels shown in grey for comparison',
      x = 'Country',
      y = 'Childhood LRI mortality rate per 1000',
      color = 'GBD Super Region'
    ) +
    # Custom colors for super regions
    scale_color_manual(values=spr_reg_colors) +
    theme_bw() +
    theme(
      legend.position = 'bottom',
      axis.text.x = element_text(angle = 90, hjust = .5),
      text = element_text(size=18)
    )
  
  
  fig2 <- ggplot(stat_ly,
                 aes(x=reorder(ISO3, plot_order), y=0, ymin=A2_relmin, ymax=A2_relmax,
                     color=super_region)
  ) +
    geom_crossbar(
      data=stat_fy, color='#CCCCCC', size=1.8, width=0.05
    ) +
    geom_point(
      data=stat_fy, color='#000000', alpha=.25, shape=18, size=1.8
    ) +
    # Plot last year data colored by GBD super region
    geom_crossbar(size=1.8, width=0.05,
                  position=position_nudge(x=.3)
    ) +
    geom_point(color='#000000', alpha=.5, shape=18, size=1.8,
               position=position_nudge(x=.3)
    )+
    labs(
      title = paste0('Relative deviation from the country mean in ', end_year),
      subtitle = 'Maximum and minimum second administrative unit deviations shown as error bars\n2000 levels shown in grey for comparison',
      x = 'Country',
      y = 'Relative deviation',
      color = 'GBD Super Region'
    ) +
    # Custom colors for super regions
    scale_color_manual(values=spr_reg_colors) +
    theme_bw() +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(angle = 90, hjust = .5),
      text = element_text(size=18)
    )
  
  
  pdf(paste0(out_dir, pdf_name),
      height=14, width=22
  )
  
  #make min superregion map
  spr_reg_colors_high_inc <- c("#6F4070","#D4D2D4","#CC503E","#E17C05","#73AF48","#0F8554","#1D6996")
  spr_reg_map <- make_spr_reg_map(spr_reg_colors = spr_reg_colors_high_inc)
  
  # Draw the high-low plot across the entire image
  grid.arrange(fig1, fig2, nrow = 2, heights=c(8,3.5))
  vp <- grid::viewport(
    x = unit(.03, 'npc'),
    y = unit(.92, 'npc'),
    width = unit(.2, 'npc'),
    height= unit(.16, 'npc'),
    just = c('left','top')
  )
  grid::pushViewport(vp)
  vp_mini <- grid::viewport(width = 1.2, height = 1, x = 0.67, y = 0.45)
  print(spr_reg_map, vp = vp_mini)
  dev.off()
  
  #save out inequality csvs
  save_dir <- '<<<< FILEPATH REDACTED >>>>'
  dir.create(save_dir)
  write.csv(lri_1, '<<<< FILEPATH REDACTED >>>>')
  write.csv(lri_2, '<<<< FILEPATH REDACTED >>>>')
  
  # (3) Prep incidence data ############################################################################
  
  lri_0 <- read.csv('<<<< FILEPATH REDACTED >>>>')[,-1]
  lri_1 <- read.csv('<<<< FILEPATH REDACTED >>>>')[,-1]
  lri_2 <- read.csv('<<<< FILEPATH REDACTED >>>>')[,-1]
  
  ULocs <- unique(lri_0$ADM0_NAME)
  Which <- rep(NA, length(ULocs))
  for (i in 1:length(ULocs)){
    tmp <- which(as.character(locs$location_name) == as.character(ULocs[i]))
    if (length(tmp)){
      Which[i] <- tmp
    }
  }
  
  #manually fix mismatches due to special characters
  print('Expecting NAs at 12 15 31 73 81')
  print('NAs at')
  print(which(is.na(Which)))
  
  Which[12] <- 896
  Which[15] <- 796
  Which[31] <- 897
  Which[73] <- 583
  Which[81] <- 906
  
  lri_0$super_region <- as.character(lri_0$region)
  lri_1$super_region <- as.character(lri_1$region)
  lri_2$super_region <- as.character(lri_2$region)
  
  lri_1$abs_lri <- lri_1$rel_lri <- lri_1$a0_lri <- lri_1$mean
  lri_2$abs_lri <- lri_2$rel_lri <- lri_2$a0_lri <- lri_2$mean
  
  #assign GBD superregion to a0,a1,a2 files
  for (year in 2000:end_year){
    for (i in 1:length(ULocs)){
      tmp_a0 <- which(as.character(lri_0$ADM0_NAME) == as.character(ULocs)[i] & lri_0$year == year)
      lri_0$super_region[tmp_a0]<- as.character(locs$super_region_name)[Which[i]]
      
      tmp_a1 <- which(as.character(lri_1$ADM0_NAME) == as.character(ULocs)[i] & lri_1$year == year)
      lri_1$super_region[tmp_a1]  <- as.character(locs$super_region_name[Which[i]])
      
      tmp_a2 <- which(as.character(lri_2$ADM0_NAME) == as.character(ULocs)[i] & lri_2$year == year)
      lri_2$super_region[tmp_a2]  <- as.character(locs$super_region_name)[Which[i]]
      
      #assign absolute and relative deviations
      if (length(tmp_a1)>1){
        lri_1$a0_lri[tmp_a1] <- lri_0$mean[tmp_a0]
        lri_1$abs_lri[tmp_a1] <- (lri_1$mean[tmp_a1] - lri_0$mean[tmp_a0])
        lri_1$rel_lri[tmp_a1] <- (lri_1$mean[tmp_a1] - lri_0$mean[tmp_a0]) / lri_0$mean[tmp_a0]
        
        lri_2$a0_lri[tmp_a2] <- lri_0$mean[tmp_a0]
        lri_2$abs_lri[tmp_a2] <- (lri_2$mean[tmp_a2] - lri_0$mean[tmp_a0])
        lri_2$rel_lri[tmp_a2] <- (lri_2$mean[tmp_a2] - lri_0$mean[tmp_a0]) / lri_0$mean[tmp_a0]
      }
    }
  }
  
  # get iso3s
  isos <- fread('<<<< FILEPATH REDACTED >>>>')
  isos[, grep('V', names(isos)) := NULL]
  isos[, ADM0_CODE := NULL]
  setnames(isos, 'iso3', 'ISO3')
  lri_0 <- data.table(merge(lri_0, isos[, c('ISO3', 'location_name')], by.x = 'ADM0_NAME', by.y = 'location_name', all.x = T))
  lri_1 <- data.table(merge(lri_1, isos[, c('ISO3', 'location_name')], by.x = 'ADM0_NAME', by.y = 'location_name', all.x = T))
  lri_2 <- data.table(merge(lri_2, isos[, c('ISO3', 'location_name')], by.x = 'ADM0_NAME', by.y = 'location_name', all.x = T))
  
  # clean up manually
  lri_0[ADM0_NAME == "Côte d'Ivoire", ISO3 := 'CIV']
  lri_0[ADM0_NAME == "São Tomé and Príncipe", ISO3 := 'STP']
  lri_1[ADM0_NAME == "Côte d'Ivoire", ISO3 := 'CIV']
  lri_1[ADM0_NAME == "São Tomé and Príncipe", ISO3 := 'STP']
  lri_2[ADM0_NAME == "Côte d'Ivoire", ISO3 := 'CIV']
  lri_2[ADM0_NAME == "São Tomé and Príncipe", ISO3 := 'STP']
  lri_0[ADM0_NAME == "Côte d'Ivoire", super_region := 'Sub-Saharan Africa']
  lri_0[ADM0_NAME == "São Tomé and Príncipe", super_region := 'Sub-Saharan Africa']
  lri_1[ADM0_NAME == "Côte d'Ivoire", super_region := 'Sub-Saharan Africa']
  lri_1[ADM0_NAME == "São Tomé and Príncipe", super_region := 'Sub-Saharan Africa']
  lri_2[ADM0_NAME == "Côte d'Ivoire", super_region := 'Sub-Saharan Africa']
  lri_2[ADM0_NAME == "São Tomé and Príncipe", super_region := 'Sub-Saharan Africa']
  lri_0[ISO3 == "COG", super_region := 'Sub-Saharan Africa']
  lri_1[ISO3 == "COG", super_region := 'Sub-Saharan Africa']
  lri_2[ISO3 == "COG", super_region := 'Sub-Saharan Africa']
  lri_0[ISO3 == "GMB", super_region := 'Sub-Saharan Africa']
  lri_1[ISO3 == "GMB", super_region := 'Sub-Saharan Africa']
  lri_2[ISO3 == "GMB", super_region := 'Sub-Saharan Africa']
  lri_0[ISO3 == "PSE", super_region := 'North Africa and Middle East']
  lri_1[ISO3 == "PSE", super_region := 'North Africa and Middle East']
  lri_2[ISO3 == "PSE", super_region := 'North Africa and Middle East']
  
  # remove countries where we have no data
  lri_0 <- lri_0[!grep(paste0(no_data_countries, collapse = '|'), ISO3),]
  lri_1 <- lri_1[!grep(paste0(no_data_countries, collapse = '|'), ISO3),]
  lri_2 <- lri_2[!grep(paste0(no_data_countries, collapse = '|'), ISO3),]
  
  #remove admins with NA estimates
  lri_0 <- lri_0[!is.na(ISO3) & !is.na(mean)]
  lri_1 <- lri_1[!is.na(ISO3) & !is.na(mean)]
  lri_2 <- lri_2[!is.na(ISO3) & !is.na(mean)]
  
  # (4) Save incidence csvs ############################################################################
  write.csv(lri_1, '<<<< FILEPATH REDACTED >>>>')
  write.csv(lri_2, '<<<< FILEPATH REDACTED >>>>')
}
