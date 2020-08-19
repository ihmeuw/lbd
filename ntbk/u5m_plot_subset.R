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
library(rgdal)
library(dplyr)

run_date  <- "2019_02_28_20_00_15"
location_code_type <- 'gadm' # gaul or gadm, Corresponding to run date
adm_level <- 2
first_yr  <- 2000
last_yr   <- 2017
shp_version <- '2019_02_27'
focus <- FALSE


lookup_file <- "<<<< FILEPATH REDACTED >>>>"
if(location_code_type=='gaul'){
  hierarchy_path <- "<<<< FILEPATH REDACTED >>>>"
} else {
  hierarchy_path <- "<<<< FILEPATH REDACTED >>>>"
}

## Prep data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load all input datasets

# Load hierarchy matching lower admin units to country units
hierarchy <- as.data.table( read.dbf(hierarchy_path) )
# Load lookup table matching countries to stages
lookup_table <- fread(lookup_file)

## Create merge table containing all Stage 2 admin(X) units
# Update column names in hierarchy table
setnames(
  hierarchy,
  c( paste0('ADM',adm_level,c('_CODE','_NAME')), paste0('ADM0_',c('CODE','NAME')) ),
  c('adm_code','adm_name','country_code','country_name')
)
hierarchy <- hierarchy[, .(adm_code, adm_name, country_code, country_name)]
# Update column names in stage lookup table
if(location_code_type == 'gaul'){
  setnames( lookup_table, 'GAUL_CODE', "country_code" )
} else {
  setnames( lookup_table, 'gadm_geoid', "country_code" )
}
# Keep only stage 2
lookup_table <- lookup_table[ Stage %in% c('1','2a','2b'), ]
# Merge tables
loc_meta <- merge(
  x  = lookup_table[, .(country_code, spr_reg_nm, iso3)],
  y  = hierarchy[, .(country_code, country_name, adm_code, adm_name)],
  by = c('country_code')
)
# (Fix location of Uruguay)
loc_meta[iso3=='URY',spr_reg_nm := 'Latin America and Caribbean']
# Remove missing countries from the list
loc_meta <- loc_meta[ !(iso3 %in% c('ESH','GUF', 'CPV', 'CUB',
  'DZA', 'IRN', 'PSE', 'SYR', 'TTO', 'MYS', 'VEN')), ]

# load data function
# load and processdata
load_dat <- function(indi, rd) {
  setwd(paste0("<<<< FILEPATH REDACTED >>>>",
          indi, '/output/', rd, '/'))
  mydat <- do.call(rbind, lapply(list.files(pattern =
              paste0(indi, '_unraked_admin_draws')), function(x) {
    load(x)
    return(admin_2)
  }))

  mydat <- as.data.frame(mydat)
  draws <- mydat[, grep('V', names(mydat))]
  mydat$value <- apply(draws, 1, mean, na.rm = TRUE)
  # mydat$mean <- mydat$mean/3
  mydat <- select(mydat, value, pop, ADM2_CODE, year)
  return(mydat)
}

## Prep Q and merge onto location metadata
load_and_prep_data <- function(x) {
  if (x == 'sani') {
    q <- as.data.frame(load_dat('s_piped', '2019_02_28_20_00_15')) %>%
          rename(piped = value)
    q2 <- as.data.frame(load_dat('s_imp_cr', '2019_02_28_20_00_15'))  %>%
            rename(imp_cr = value)
    q3 <- left_join(q, q2, by = c('ADM2_CODE', 'year', 'pop'))
    q3 <- mutate(q3, imp = piped + ((1-piped)*imp_cr)) %>%
            rename(value = imp) %>%
            select(value, pop, ADM2_CODE, year)
    q <- as.data.table(q3)
    rm(list = c('q2', 'q3'))
  } else {
    q <- as.data.frame(load_dat('w_piped', '2019_02_28_20_00_15')) %>%
          rename(piped = value)
    q2 <- as.data.frame(load_dat('w_imp_cr', '2019_02_28_20_00_15'))  %>%
            rename(imp_cr = value)
    q3 <- left_join(q, q2, by = c('ADM2_CODE', 'year', 'pop'))
    q3 <- mutate(q3, imp = piped + ((1-piped)*imp_cr)) %>%
            rename(value = imp) %>%
            select(value, pop, ADM2_CODE, year)
    q <- as.data.table(q3)
    rm(list = c('q2', 'q3'))
  }

  setnames(q, c(paste0('ADM',adm_level,'_CODE'),'value'), c('adm_code','q') )
  # Change to 5q0 PER 1000
  q[, V1 := NULL]

  q_meta <- merge(
    x = q,
    y = loc_meta,
    by = c('adm_code')
  )


  # Get highest, median, and lowest adm2 values by country-year
  q_stats <- q_meta[, .(q_min = min(q, na.rm=TRUE),
                        q_med = weighted.mean(x = q, w = pop, na.rm=TRUE),
                        q_max= max(q, na.rm=TRUE)),
                    by = .(country_code, country_name, iso3, spr_reg_nm, year)]
  q_stats[q_min < 0, q_min := 0]
  q_stats[q_max > 1000, q_max := 1000]

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
  stat_ly <- q_stats[year==last_yr,]
  stat_fy <- q_stats[year==first_yr,]
  # Set the plotting order (lowest to highest by country)
  stat_ly <- stat_ly[order(iso3)]
  stat_ly[,plot_order := .I]

  # Merge plotting order back onto the first year of data
  stat_fy <- merge(
    x = stat_fy,
    y = stat_ly[, .(iso3, plot_order)],
    by = c('iso3')
  )
  stat_fy <- stat_fy[order(plot_order)]

  return(list(stat_ly, stat_fy))

}

sani_dat <- load_and_prep_data('sani')
sani_ly <- sani_dat[[1]]
sani_ly$indi_group <- 'Improved Sanitation'
sani_fy <- sani_dat[[2]]
sani_fy$indi_group <- 'Improved Sanitation'

water_dat <- load_and_prep_data('water')
water_ly <- water_dat[[1]]
water_ly$indi_group <- 'Improved Water'
water_fy <- water_dat[[2]]
water_fy$indi_group <- 'Improved Water'

plt_order <- select(water_ly, iso3, plot_order)
sani_ly <- select(sani_ly, -plot_order)
sani_ly <- left_join(sani_ly, plt_order)
water_fy <- select(water_fy, -plot_order)
water_fy <- left_join(water_fy, plt_order)
sani_fy <- select(sani_fy, -plot_order)
sani_fy <- left_join(sani_fy, plt_order)

stat_ly <- as.data.table(as.data.frame(rbind(water_ly, sani_ly)))
stat_fy <- as.data.table(as.data.frame(rbind(water_fy, sani_fy)))
# Load Q for all admin2 units

if (focus) {
  stat_fy <- sani_fy
  stat_ly <- sani_ly
  stat_fy <- filter(stat_fy, country_name %in% c('Peru', 'Ethiopia', 'Cambodia', 'Nigeria'))
  stat_ly <- filter(stat_ly, country_name %in% c('Peru', 'Ethiopia', 'Cambodia', 'Nigeria'))
}

make_spr_reg_map <- function(spr_reg_colors){
  # Load master shapefile
  master_shape <- readOGR(
    "<<<< FILEPATH REDACTED >>>>"
  )
  # Merge on location hierarchy info
  locs <- fread("<<<< FILEPATH REDACTED >>>>")
  locs <- locs[,.(loc_id, ihme_lc_id)]
  superregion_ids <- fread("<<<< FILEPATH REDACTED >>>>")
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
  master[master$ADM0_A3 == "SOL", "super_region_name"] <- "Sub-Saharan Africa"
  master[master$ADM0_A3 == "SAH", "super_region_name"] <- "North Africa and Middle East"
  master[master$ADM0_A3 == "PSX", "super_region_name"] <- "North Africa and Middle East"
  master[master$ADM0_A3 == "GUF", "super_region_name"] <- "Latin America and Caribbean"
  master[master$ADM0_A3 == "GRL", "super_region_name"] <- "High-income"
  # Grey out China, Mexico, and Brazil
  master[is.na(master$Stage), 'Stage'] <- "4"
  master[master$Stage == "3", 'super_region_name'] <- "High-income"
  # Fix to keep in Mongolia
  mng <- fortify(master[master$ADM0_A3 == "MNG",], region = "super_region_name")
  #merge on superregion info
  x <- fortify(master, region = "super_region_name")

  #plot
  minimap <- ggplot() +
    geom_polygon(data = x, aes(x=long, y = lat, group = group, fill = id))+
    geom_polygon(data = mng, aes(x=long, y = lat, group = group, fill = id))+
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



## Create and print plot

# Custom colors for super regions
spr_reg_colors <- c("#6F4070","#CC503E","#E17C05","#73AF48","#0F8554","#1D6996")
spr_reg_colors_w_high_income <- c("#6F4070","#808080","#CC503E","#E17C05",
                                  "#73AF48","#0F8554","#1D6996")

# Create a mini-map of super regions
minimap <- make_spr_reg_map(spr_reg_colors_w_high_income)

fig1 <- ggplot() +
          geom_crossbar(data = stat_ly,
              aes(x=reorder(country_name, plot_order), y=q_med, ymin=q_min, ymax=q_max),
                  color="#CC503E",
                  size=2.2, width=0,
                        position=position_nudge(x=.32)) +
          # Plot first year data in grey
          geom_crossbar(aes(x=reorder(country_name, plot_order), y=q_med, ymin=q_min, ymax=q_max),
            data=stat_fy, color='#CCCCCC', size=2.2, width=0
          ) +
          geom_point(
            data=stat_fy, color='#000000', alpha=.25, shape=18, size=3,
            aes(x=reorder(country_name, plot_order), y=q_med)
          ) +
          geom_point(data = stat_ly, color='#000000', alpha=.5, shape=18, size=3,
                     position=position_nudge(x=.32),
            aes(x=reorder(country_name, plot_order), y=q_med)
          ) +
          labs(
            x = 'Country',
            y = 'Access to Improved Sanitation',
            color = 'GBD Super Region'
          ) +
          theme_bw() +
          theme(
            legend.position = 'bottom',
            text = element_text(size=16),
            axis.text.x = element_text(size=16, hjust = .5)
          )
fig1

fig2 <- fig1 + facet_grid(rows = vars(indi_group))

out_dir <- "<<<< FILEPATH REDACTED >>>>"

pdf(
  paste0(out_dir,'high_low_admin',adm_level,'_ranked_focus.pdf'),
  height=8, width=22
)
# Draw the high-low plot across the entire image
grid.draw(ggplotGrob(fig2))
#Draw the minimap at the top left of the plot
vp <- grid::viewport(
  x = unit(.07, 'npc'),
  y = unit(.11, 'npc'),
  width = unit(.07, 'npc'),
  height= unit(.1, 'npc'),
  just = c('left','top')
)
grid::pushViewport(vp)
grid.draw( ggplotGrob(minimap) )
dev.off()



png(
  paste0(out_dir,'high_low_admin',adm_level,'_ranked_focus.png'),
  height=800, width=1600
)
# Draw the high-low plot across the entire image
grid.draw(ggplotGrob(fig2))
#Draw the minimap at the top left of the plot
vp <- grid::viewport(
  x = unit(.07, 'npc'),
  y = unit(.11, 'npc'),
  width = unit(.07, 'npc'),
  height= unit(.1, 'npc'),
  just = c('left','top')
)
grid::pushViewport(vp)
grid.draw( ggplotGrob(minimap) )
dev.off()
