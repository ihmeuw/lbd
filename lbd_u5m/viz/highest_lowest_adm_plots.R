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


## Define inputs and filepaths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_date  <- "<<<< REDACTED >>>>"
location_code_type <- 'gadm' # gaul or gadm, Corresponding to run date
adm_level <- 2
first_yr  <- 2000
last_yr   <- 2017
shp_version <- '<<<< REDACTED >>>>'

in_dir    <- paste0('<<<< FILEPATH REDACTED >>>>',run_date,'/')
out_dir   <- paste0(
  '<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/'
)

lookup_file <- '<<<< FILEPATH REDACTED >>>>'
if(location_code_type=='gaul'){
  hierarchy_path <- paste0("<<<< FILEPATH REDACTED >>>>",
    "<<<< FILEPATH REDACTED >>>>", adm_level, "<<<< FILEPATH REDACTED >>>>", adm_level, "<<<< FILEPATH REDACTED >>>>", 
    adm_level, "<<<< FILEPATH REDACTED >>>>"
  )
} else {
  hierarchy_path <- paste0('<<<< FILEPATH REDACTED >>>>',
    shp_version,'<<<< FILEPATH REDACTED >>>>',adm_level,'.dbf'
  )
}
dir.create(out_dir, showWarnings = FALSE)

## Prep data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Load all input datasets
# Load Q for all admin2 units
q <- fread(paste0(in_dir,'<<<< FILEPATH REDACTED >>>>',adm_level,'.csv'))
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
loc_meta <- loc_meta[ !(iso3 %in% c('ESH','GUF','MEX','CHN','BRA')), ]

## Prep Q and merge onto location metadata
setnames(q, c(paste0('ADM',adm_level,'_CODE'),'value'), c('adm_code','q') )
# Change to 5q0 PER 1000
q[, q := q * 1000]
q[, V1 := NULL]

q_meta <- merge(
  x = q,
  y = loc_meta,
  by = c('adm_code')
)
q_meta[q < 0, q := 0]
q_meta[q > 1000, q := 1000]


# Get highest, median, and lowest adm2 values by country-year
q_stats <- q_meta[, .(q_min = min(q, na.rm=TRUE), 
                      q_med = median(q, na.rm=TRUE), 
                      q_max= max(q, na.rm=TRUE)),
                  by = .(country_code, country_name, iso3, spr_reg_nm, year)]
q_stats[q_min < 0, q_min := 0]
q_stats[q_max > 1000, q_max := 1000]
q_stats[, rel_min := q_min / q_med ]
q_stats[, rel_max := q_max / q_med ]
q_stats[, rel_med := 1 ]

# Merge back onto full q data to get relative q for each district
q_meta_merged <- merge(
  x = q_meta,
  y = q_stats[, .(country_code, year, q_med)],
  by = c('country_code', 'year')
)
q_meta_merged[, q_rel := q / q_med ]
q_meta_merged_fy <- q_meta_merged[year == first_yr,]
q_meta_merged_ly <- q_meta_merged[year == last_yr,]

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
  master[master$ADM0_A3 == "PSX", "super_region_name"] <- "North Africa and Middle East"
  master[master$ADM0_A3 == "GRL", "super_region_name"] <- "High-income"
  # Grey out China, Mexico, and Brazil
  master[master$ADM0_A3 == "SAH", "super_region_name"] <- "High-income"
  master[master$ADM0_A3 == "GUF", "super_region_name"] <- "High-income"
  master[master$ADM0_A3 == "CHN", "super_region_name"] <- "High-income"
  master[master$ADM0_A3 == "MEX", "super_region_name"] <- "High-income"
  master[master$ADM0_A3 == "BRA", "super_region_name"] <- "High-income"
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

stat_ly <- q_stats[year==last_yr,]
stat_fy <- q_stats[year==first_yr,]
# Set the plotting order (lowest to highest by country)
stat_ly <- stat_ly[order(q_med)]
stat_ly[,plot_order := .I]

# Merge plotting order back onto the first year of data and onto other datasets
stat_fy <- merge(
  x = stat_fy,
  y = stat_ly[, .(iso3, plot_order)],
  by = c('iso3')
)
stat_fy <- stat_fy[order(plot_order)]
q_meta_merged_fy <- merge(
  x = q_meta_merged_fy,
  y = stat_ly[, .(iso3, plot_order)],
  by = c('iso3')
)
q_meta_merged_fy <- q_meta_merged_fy[order(plot_order)]
q_meta_merged_ly <- merge(
  x = q_meta_merged_ly,
  y = stat_ly[, .(iso3, plot_order)],
  by = c('iso3')
)
q_meta_merged_ly <- q_meta_merged_ly[order(plot_order)]


## Create and print plot

# Custom colors for super regions
spr_reg_colors <- c("#6F4070","#CC503E","#E17C05","#73AF48","#0F8554","#1D6996")
spr_reg_colors_w_high_income <- c("#6F4070","#808080","#CC503E","#E17C05",
                                  "#73AF48","#0F8554","#1D6996")

# Create a mini-map of super regions
minimap <- make_spr_reg_map(spr_reg_colors_w_high_income)

fig1 <- ggplot(stat_ly,
              aes(x=reorder(iso3, plot_order), y=q_med, ymin=q_min, ymax=q_max,
                  color=spr_reg_nm)
              ) +
  # Plot first year data in grey
  geom_crossbar(
    data=stat_fy, color='#CCCCCC', size=2.2, width=0
  ) +
  geom_point(
    data=stat_fy, color='#000000', alpha=.25, shape=18, size=2.2
  ) +
  # Plot last year data colored by GBD super region
  geom_crossbar(size=2.2, width=0,
                position=position_nudge(x=.32)
  ) + 
  geom_point(color='#000000', alpha=.5, shape=18, size=2.2,
             position=position_nudge(x=.32)
  ) + 
  labs(
    x = 'Country',
    y = 'Under-5 mortality per 1000 live births (U5MR)',
    color = 'GBD Super Region'
  ) +
  # Custom colors for super regions
  scale_color_manual(values=spr_reg_colors) + 
  theme_bw() + 
  theme(
    legend.position = 'bottom', 
    text = element_text(size=20),
    axis.text.x = element_text(size=12, angle = 90, hjust = .5)
  )

pdf(
  paste0(out_dir,'high_low_admin',adm_level,'_ranked.pdf'),
  height=8, width=22
)
# Draw the high-low plot across the entire image
grid.draw( ggplotGrob(fig1) )
# Draw the minimap at the top left of the plot
vp <- grid::viewport(
  x = unit(.04, 'npc'),
  y = unit(.97, 'npc'),
  width = unit(.3, 'npc'),
  height= unit(.33, 'npc'),
  just = c('left','top')
)
grid::pushViewport(vp)
grid.draw( ggplotGrob(minimap) )
# Finish the figure
dev.off()


## Plot second graph: absolute + relative ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fig2a <- ggplot(stat_ly,
              aes(x=reorder(iso3, plot_order), y=q_med, ymin=q_min, ymax=q_max,
                  color=spr_reg_nm)
              ) +
  # Plot first year data in grey
  geom_crossbar(
    data=stat_fy, color='#CCCCCC', size=2.2, width=0
  ) +
  geom_point(
    data=stat_fy, color='#000000', alpha=.25, shape=18, size=2.2
  ) +
  # Plot last year data colored by GBD super region
  geom_crossbar(size=2.2, width=0,
                position=position_nudge(x=.32)
  ) + 
  geom_point(color='#000000', alpha=.5, shape=18, size=2.2,
             position=position_nudge(x=.32)
  ) + 
  labs(
    title = 'Absolute',
    x = 'Country',
    y = '\nUnder-5 mortality per 1000 live births (U5MR)'
  ) +
  # Custom colors for super regions
  scale_color_manual(values=spr_reg_colors) + 
  theme_bw() + 
  theme(
    legend.position = 'none',
    text = element_text(size=20),
    axis.text.x = element_text(size=12, angle = 90, hjust = .5)
  )


fig2b <- ggplot(stat_ly,
              aes(x=reorder(iso3, plot_order), y=rel_med, ymin=rel_min, ymax=rel_max,
                  color=spr_reg_nm)
              ) +
  geom_crossbar(
    data=stat_fy, color='#CCCCCC', size=2.2, width=0
  ) +
  geom_point(
    data=stat_fy, color='#000000', alpha=.25, shape=18, size=2.2
  ) +
  # Plot last year data colored by GBD super region
  geom_crossbar(size=2.2, width=0,
                position=position_nudge(x=.32)
  ) + 
  geom_point(color='#000000', alpha=.5, shape=18, size=2.2,
             position=position_nudge(x=.32)
  ) + 
  labs(
    title='Relative', 
    x='', 
    y='U5MR relative to\nnational mean', 
    color='GBD Super Region'
  ) +
  scale_y_continuous(
    breaks=0:4, 
    limits=c(0,4), 
    labels=c('0','1','2','3','    4')
  ) + 
  scale_color_manual(values=spr_reg_colors) + 
  theme_bw() + 
  theme(
    legend.position = 'bottom',
    text = element_text(size=20),
    axis.text.x = element_text(size=12, angle = 90, hjust = .5)
  )

fig2_layout <- rbind(1,1,2)

pdf(
  paste0(out_dir,'high_low_admin',adm_level,'_with_rel.pdf'),
  height=12, width=22
)
# Draw the high-low plot across the entire image
grid.arrange(fig2a, fig2b, layout_matrix=fig2_layout)
# Draw the minimap at the top left of the plot
vp <- grid::viewport(
  x = unit(.06, 'npc'),
  y = unit(.95, 'npc'),
  width = unit(.3, 'npc'),
  height= unit(.25, 'npc'),
  just = c('left','top')
)
grid::pushViewport(vp)
grid.draw( ggplotGrob(minimap) )
# Finish the figure
dev.off()


## Plot third graph: absolute + relative with dots instead of lines ~~~~~~~~~~~~

fig3a <- ggplot(stat_ly,
              aes(x=reorder(iso3, plot_order), y=q_med, color=spr_reg_nm)
              ) +
  # Plot first year data in grey
  geom_crossbar(
    data=stat_fy, aes(ymin=q_min, ymax=q_max), 
    color='#CCCCCC', size=0, width=.32, alpha=0.5
  ) +
  # Plot last year data colored by GBD super region
  geom_crossbar(
    aes(ymin=q_min, ymax=q_max),
    size=0, width=.32, position=position_nudge(x=.32), alpha=0.5
  ) + 
  # Plot first year data in grey
  geom_point(
    data=q_meta_merged_fy, aes(y=q), color='#CCCCCC', size=1.6, alpha=0.2
  ) +
  geom_point(
    data=stat_fy, color='#000000', alpha=.5, shape=18, size=1.6
  ) +
  # Plot last year data colored by GBD super region
  geom_point(
    data=q_meta_merged_ly, aes(y=q), size=1.6, alpha=0.2, 
    position=position_nudge(x=.32)
  ) + 
  geom_point(color='#000000', alpha=.5, shape=18, size=1.6,
             position=position_nudge(x=.32)
  ) + 
  labs(
    title = 'Absolute',
    x = 'Country',
    y = '\nUnder-5 mortality per 1000 live births (U5MR)'
  ) +
  # Custom colors for super regions
  scale_color_manual(values=spr_reg_colors) + 
  theme_bw() + 
  theme(
    legend.position = 'none',
    text = element_text(size=20),
    axis.text.x = element_text(size=12, angle = 90, hjust = .5),
    panel.grid.major.x = element_blank()
  )


fig3b <- ggplot(stat_ly,
              aes(x=reorder(iso3, plot_order), y=rel_med,
                  color=spr_reg_nm)
              ) +
  # Plot first year data in grey
  geom_crossbar(
    data=stat_fy, aes(ymin=rel_min, ymax=rel_max), 
    color='#CCCCCC', size=0, width=.32, alpha=0.5
  ) +
  # Plot last year data colored by GBD super region
  geom_crossbar(
    aes(ymin=rel_min, ymax=rel_max),
    size=0, width=.32, position=position_nudge(x=.32), alpha=0.5
  ) + 
  geom_point(
    data=q_meta_merged_fy, aes(y=q_rel), color='#CCCCCC', size=1.6, alpha=0.2
  ) +
  geom_point(
    data=stat_fy, color='#000000', alpha=.25, shape=18, size=1.6
  ) +
  # Plot last year data colored by GBD super region
  geom_point(
    data=q_meta_merged_ly, aes(y=q_rel), size=1.6, alpha=0.2, 
    position=position_nudge(x=.32)
  ) + 
  geom_point(color='#000000', alpha=.5, shape=18, size=1.6,
             position=position_nudge(x=.32)
  ) + 
  labs(
    title='Relative', 
    x='', 
    y='U5MR relative to\nnational mean', 
    color='GBD Super Region'
  ) +
  scale_y_continuous(
    breaks=0:4, 
    limits=c(0,4), 
    labels=c('0','1','2','3','    4')
  ) + 
  scale_color_manual(values=spr_reg_colors) + 
  theme_bw() + 
  theme(
    legend.position = 'bottom',
    text = element_text(size=20),
    axis.text.x = element_text(size=12, angle = 90, hjust = .5),
    panel.grid.major.x = element_blank()
  )

fig3_layout <- rbind(1,1,2)

pdf(
  paste0(out_dir,'<<<< FILEPATH REDACTED >>>>',adm_level,'<<<< FILEPATH REDACTED >>>>'),
  height=12, width=22
)
# Draw the high-low plot across the entire image
grid.arrange(fig3a, fig3b, layout_matrix=fig3_layout)
# Draw the minimap at the top left of the plot
vp <- grid::viewport(
  x = unit(.06, 'npc'),
  y = unit(.95, 'npc'),
  width = unit(.3, 'npc'),
  height= unit(.25, 'npc'),
  just = c('left','top')
)
grid::pushViewport(vp)
grid.draw( ggplotGrob(minimap) )
# Finish the figure
dev.off()
