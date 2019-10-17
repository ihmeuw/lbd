## ###########################################################################
## 
## PLOT HIGHEST AND LOWEST ADM2 UNITS BY COUNTRY IN 2000 and 2016
## 
## Purpose: Plot highest and lowest admin2 units for under5 mortality by 
##   country in 2000 and 2017
## 
## ###########################################################################

rm(list=ls())

## Imports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(data.table)
library(foreign)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(raster)

core_repo <- '<<<< FILEPATH REDACTED >>>>'
u5m_repo <- '<<<< FILEPATH REDACTED >>>>'
source(paste0(u5m_repo,'<<<< FILEPATH REDACTED >>>>'))
source(paste0(core_repo,'<<<< FILEPATH REDACTED >>>>'))

## Define inputs and filepaths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_date      <- "<<<< REDACTED >>>>"
first_yr      <- 2000
last_yr       <- 2017
high_quantile <- .90
low_quantile  <- .10
raked         <- TRUE # Should these plots use raked or unraked rasters?
resume        <- TRUE # Should the existing admin0 raster be used?

in_dir    <- paste0('<<<< FILEPATH REDACTED >>>>',run_date,'/')
ad0_fp <- '<<<< FILEPATH REDACTED >>>>'
out_dir   <- paste0(
  '<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/'
)
lookup_file <- '<<<< FILEPATH REDACTED >>>>'
dir.create(out_dir, showWarnings = FALSE)

## Prep data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load all input datasets
# Load Q for all pixels
in_q_path <- paste0(
  in_dir,'died_under5_mean_',ifelse(raked,'raked','unraked'),
  '_',first_yr,'_',last_yr,'.tif'
)
q_ras <- brick(in_q_path)

# Load the admin0 dataset
if(resume){
  ad0_ras <- get(load(ad0_fp))
} else {
  ad0_ras <- GetAdmin(
    admin_level = 0,
    simple_raster = q_ras
  )[['rast']]
  save(ad0_ras, file=ad0_fp)
}

## Load lookup table matching countries to stages
loc_meta <- fread(lookup_file)
loc_meta <- loc_meta[ Stage != '3', ]
setnames(
  loc_meta,
  c('GAUL_CODE','location_name'),
  c('country_code','country_name')
)
loc_meta <- loc_meta[,.(country_code, country_name, spr_reg_nm, iso3)]
## Location table fixes
# (Fix location of Uruguay)
loc_meta[iso3=='URY',spr_reg_nm := 'Latin America and Caribbean']
loc_meta <- loc_meta[iso3 != 'GNQ',] # Covariate issue


## Join all data
full_data <- data.table(
  q            = as.vector(q_ras),
  country_code = as.vector(ad0_ras),
  year         = rep(first_yr:last_yr, each=prod(dim(q_ras)[1:2]))
)
# Change to 5q0 PER 1000
full_data[, q := q * 1000]
# Drop any NA values
full_data <- na.omit(full_data)
message(paste("After dropping NA rows,",nrow(full_data),"rows of data remain."))

# Merge onto admin metadata
full_data <- merge(
  x = full_data,
  y = loc_meta,
  by = c('country_code')
)

if(any( is.na(full_data[,country_name]) )) stop("Issue with location meta merge")

## Generate summary statistics
## Get 10th pctile, median, and 90th pctile values by country-year
q_stats <- full_data[, .(q_min = quantile(q, .1, na.rm=TRUE), 
                         q_med = median(  q,     na.rm=TRUE), 
                         q_max = quantile(q, .9, na.rm=TRUE)),
                     by = .(country_code, country_name, iso3, spr_reg_nm, year)]
q_stats[q_min < 0, q_min := 0]
q_stats[q_max > 1000, q_max := 1]

## Create data table to show relative inequality (highest/lowest admX) by year
q_ineq <- copy(q_stats)
q_ineq <- q_ineq[ (year >= first_yr) & (year <= last_yr) ,]
q_ineq[, ratio := q_max / q_min ]
# Get normalized ratio (ratio in year / ratio in first year)
q_norm <- q_ineq[year == first_yr, .(iso3, ratio)]
setnames(q_norm, 'ratio', 'baseline')
q_ineq <- merge(
  x = q_ineq,
  y = q_norm,
  by = c('iso3')
)
q_ineq[, ratio_norm := ratio/baseline ]
# Prep a data.table to use for labels
ineq_labels <- q_ineq[ year==last_yr, .(iso3, spr_reg_nm, year, ratio, ratio_norm)]

## FUNCTION TO MAKE A SUPER-REGION MINIMAP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Returns a ggplot object containing the super region map
make_spr_reg_map <- function(spr_reg_colors){
  # Load master shapefile
  master_shape <- readOGR(
    '<<<< FILEPATH REDACTED >>>>'
  )
  # Merge on location hierarchy info
  locs <- fread("<<<< FILEPATH REDACTED >>>>")
  locs <- locs[,.(loc_id, ihme_lc_id)]
  superregion_ids <- fread(paste0("<<<< FILEPATH REDACTED >>>>",
                                  "<<<< FILEPATH REDACTED >>>>"))
  super <- superregion_ids[, .(location_id, super_region_name)]
  locs <- merge(
    locs, super, 
    by.x = c("loc_id"), 
    by.y = c("location_id"), 
    all.x=T
  )
  master <- sp::merge(master_shape, locs, by.x = "ADM0_A3", by.y = "ihme_lc_id")
  # Make Stage 3 countries gray
  stage <- fread("<<<< FILEPATH REDACTED >>>>")
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
  master[master$ADM0_A3 == "SAH", "super_region_name"] <- "North Africa and Middle East"
  master[master$ADM0_A3 == "PSX", "super_region_name"] <- "North Africa and Middle East"
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


## PLOT FIRST GRAPH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Graph theme: Differences across adminX units within country in 2016

stat_ly <- q_stats[year==last_yr,]
stat_fy <- q_stats[year==first_yr,]
# Set the plotting order (lowest to highest by country)
stat_ly <- stat_ly[order(q_med)]
stat_ly[,plot_order := .I]

# Merge plotting order back onto the first year of data
stat_fy <- merge(
  x = stat_fy,
  y = stat_ly[, .(iso3, plot_order)],
  by = c('iso3')
)
stat_fy <- stat_fy[order(plot_order)]


## Create and print plot

# Custom colors for super regions
spr_reg_colors <- c("#6F4070","#CC503E","#E17C05","#73AF48","#0F8554","#1D6996")
spr_reg_colors_w_high_income <- c("#6F4070","#808080","#CC503E","#E17C05",
                                  "#73AF48","#0F8554","#1D6996")

fig1 <- ggplot(stat_ly,
              aes(x=reorder(iso3, plot_order), y=q_med, ymin=q_min, ymax=q_max,
                  color=spr_reg_nm)
              ) +
  # Plot first year data in grey
  geom_crossbar(
    data=stat_fy, color='#CCCCCC', size=1.8, width=0
  ) +
  geom_point(
    data=stat_fy, color='#000000', alpha=.25, shape=18, size=1.8
  ) +
  # Plot last year data colored by GBD super region
  geom_crossbar(size=1.8, width=0,
                position=position_nudge(x=.3)
  ) + 
  geom_point(color='#000000', alpha=.5, shape=18, size=1.8,
             position=position_nudge(x=.3)
  ) + 
  labs(
    title = paste0('All Stage 2 countries ranked by median 5q0 across all ',
                   'pixels in ',last_yr),
    subtitle = paste0(
      '10th and 90th percentile of 5q0 by country shown as error bars; ',
      first_yr,' 5q0 levels shown in grey for comparison'
    ),
    x = 'Country',
    y = 'Under 5 mortality probability (5q0) per 1000',
    color = 'GBD Super Region'
  ) +
  # Custom colors for super regions
  scale_color_manual(values=spr_reg_colors) + 
  theme_bw() + 
  theme(
    legend.position = 'bottom', 
    axis.text.x = element_text(angle = 90, hjust = .5)
  )

# Create a mini-map of super regions
minimap <- make_spr_reg_map(spr_reg_colors_w_high_income)

pdf(
  paste0(out_dir,'high_low_pixels_ranked.pdf'),
  height=8, width=22
)
# Draw the high-low plot across the entire image
grid.draw( ggplotGrob(fig1) )
# Draw the minimap at the top left of the plot
vp <- grid::viewport(
  x = unit(.03, 'npc'),
  y = unit(.92, 'npc'),
  width = unit(.2, 'npc'),
  height= unit(.22, 'npc'),
  just = c('left','top')
)
grid::pushViewport(vp)
grid.draw( ggplotGrob(minimap) )
# Finish the figure
dev.off()
