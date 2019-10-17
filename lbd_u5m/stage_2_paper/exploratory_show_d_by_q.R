## #############################################################################
##
## STAGE 2 PAPER: Figure 4 + 5 (d by q by superregion/age)
##
## Date: January 18, 2019
##
## #############################################################################

rm(list=ls())

## IMPORTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)
library(ggplot2)
library(foreign)
library(gridExtra)
library(grid)
library(rgdal)
library(rgeos)
options(scipen=9999)


source(paste0("<<<< FILEPATH REDACTED >>>>"))
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")
load_mbg_functions("<<<< FILEPATH REDACTED >>>>")

## DEFINE INPUT AND OUTPUT INFORMATION
run_date <- "<<<< FILEPATH REDACTED >>>>"
adm_level <- 2

shapefile_version <- "<<<< FILEPATH REDACTED >>>>"
drop_china = T

minimap_save_path <- "<<<< FILEPATH REDACTED >>>>"
figure4_save_path <- "<<<< FILEPATH REDACTED >>>>"
figure5_save_path <- "<<<< FILEPATH REDACTED >>>>"
figure4_proportion_save_path <- "<<<< FILEPATH REDACTED >>>>"

j_head <- "<<<< FILEPATH REDACTED >>>>"
## Superregion minimap ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
master_shape <- readOGR("<<<< FILEPATH REDACTED >>>>")

#get location hierarchy info from somewhere
locs <- fread("<<<< FILEPATH REDACTED >>>>")
locs <- locs[,c("loc_id", "ihme_lc_id")]

superregion_ids <- fread("<<<< FILEPATH REDACTED >>>>")
super <- superregion_ids[, c("location_id", "super_region_name")]

locs <- merge(locs, super, by.x = c("loc_id"), by.y = c("location_id"), all.x=T)

master <- sp::merge(master_shape, locs, by.x = "ADM0_A3", by.y = "ihme_lc_id", all.x=T)

stage <- fread("<<<< FILEPATH REDACTED >>>>"
stage <- stage[,c("loc_id", "Stage")]

master <- sp::merge(master, stage, by = "loc_id")

#drop antarctica, mexico, and brazil
master <- master[master$ADM0_A3 != "ATA", ]

master$ADM0_A3 <- as.character(master$ADM0_A3)
master$super_region_name <- as.character(master$super_region_name)
master$Stage <- as.character(master$Stage)

if(drop_china){
  master[master$ADM0_A3 == "CHN", "super_region_name"] <- "High-income"
}

master[master$ADM0_A3 == "SDS", "ADM0_A3"] <- "SSD"
master[master$ADM0_A3 == "SSD", "super_region_name"] <- "Sub-Saharan Africa"
master[master$ADM0_A3 == "PSX", "super_region_name"] <- "North Africa and Middle East"
master[master$ADM0_A3 == "SOL", "super_region_name"] <- "Sub-Saharan Africa"
master[master$ADM0_A3 == "MNG", "super_region_name"] <- "Central Europe, Eastern Europe, and Central Asia"

#grey out western sahara and french guiana
master[master$ADM0_A3 == "SAH", "super_region_name"] <- "High-income"
master[master$ADM0_A3 == "GUF", "super_region_name"] <- "High-income"

#grey out and malaysia
master[master$ADM0_A3 == "MYS", "super_region_name"] <- "High-income"
master[master$ADM0_A3 == "LKA", "super_region_name"] <- "South Asia"

#grey out mexico and brazil
master[master$ADM0_A3 == "MEX", "super_region_name"] <- "High-income"
master[master$ADM0_A3 == "BRA", "super_region_name"] <- "High-income"

master[is.na(master$Stage), 'Stage'] <- "4"
master[master$Stage == "3", 'super_region_name'] <- "High-income"

mng <- fortify(master[master$ADM0_A3 == "MNG",], region = "super_region_name")
#merge on superregion info
x <- fortify(master, region = "super_region_name")
y <- fortify(master)
#plot

colors_mini <-rev(c('#1D6996', '#0F8554', '#73AF48', '#E17C05', '#CC503E', "#808080", "#6F4070"))

minimap <- ggplot() + 
  geom_polygon(data = x, aes(x=long, y = lat, group = group, fill = id))+
  geom_polygon(data = mng, aes(x=long, y = lat, group = group, fill = id))+
  #geom_polygon(data = y, aes(x=long, y = lat, group = group), fill = NA, color = "black")+
  theme_void() +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors_mini)

pdf(minimap_save_path,
    width=12,
    height=8)
print(minimap)
dev.off()


## PREP DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Function for assembling q and d for an admin level given a MBG run date ~~~~~
assemble_q_d <- function(run_date, adm_level){
  # Define q and d input files
  in_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  fps <- list(
    q = paste0("<<<< FILEPATH REDACTED >>>>"),
    d = paste0("<<<< FILEPATH REDACTED >>>>")
  )
  dts <- list(q=NA, d=NA)
  # Read in and prep data.tables for each entity (q and d)
  for(val in names(dts)){
    dts[[val]] <- fread(fps[[val]])
    setnames(dts[[val]], paste0('ADM',adm_level,"_CODE"), "ADM_CODE")
    dts[[val]] <- dts[[val]][, .(ADM_CODE, year, value)]
    setnames(dts[[val]], 'value', val)
    dts[[val]][, level := adm_level]
  }
  # Merge q and d by admin code and year
  q_d_merged <- merge(
    x   = dts[['q']],
    y   = dts[['d']],
    by  = c("ADM_CODE","year","level"),
    all = TRUE
  )
  # Return prepped data.table
  return(q_d_merged)
}


## Read and prep file
dq <- assemble_q_d(
  run_date = run_date,
  adm_level = adm_level
)

# Clean and format columns
dq[, year_factor := as.factor(year)]
dq[, q := q * 1000]
dq[, d := d / 1E5]
dq <- dq[order(year,q)]
dq[, q_bin := round(q/5,0)*5]
# Aggregate by year and binned q
dq_agg <- dq[, list(d=sum(d, na.rm = T), n_adm=.N), by=c('year','year_factor','q_bin')]

# Subset to relevant years and clean
dq_sub <- dq_agg[year %in% c(2000, 2005, 2010, 2015, 2017)]
dq_sub[, cumul_d := cumsum(d), by=c('year')]
dq_sub[, total_d := sum(d), by=c('year')]
dq_sub[, running_pct_d := cumul_d/total_d ]


# For first and last years only, merge on region
# Load tables connecting admin2 units to countries and countries to regions
dq_firstlast <- dq[year %in% c(2000, 2017)]
suffix      <- '.dbf'
hierarchy_path <- paste0(get_admin_shapefile(admin_level = 2, raking = F, version = shapefile_version))
hierarchy <- readOGR(hierarchy_path)
hierarchy <- as.data.table(hierarchy@data)

setnames(hierarchy, paste0('ADM',adm_level,'_CODE'), "ADM_CODE")
# Subset to only subnational admin units and national admin codes
hierarchy <- hierarchy[, .(ADM_CODE, ADM0_CODE)]
hierarchy$ADM_CODE <- as.integer(as.character(hierarchy$ADM_CODE))
hierarchy$ADM0_CODE <- as.integer(as.character(hierarchy$ADM0_CODE))

# - Admin0 to country name and region
by_stage <- fread("<<<< FILEPATH REDACTED >>>>")[, .(location_name, iso3, Stage, spr_reg_nm, gadm_geoid, mbg_reg_nm)]
setnames(by_stage, 'gadm_geoid','ADM0_CODE')

# - Combine to get merge table from adm2 to all metadata
all_meta <- merge(
  x   = hierarchy,
  y   = by_stage,
  by  = c('ADM0_CODE')
)

all_meta$ADM_CODE <- as.integer(as.character(all_meta$ADM_CODE))
# Merge onto the admin2 table
dq_meta <- merge(
  x     = dq_firstlast,
  y     = all_meta,
  by    = c('ADM_CODE'),
  all.x = TRUE
)

if(drop_china){
  dq_meta <- dq_meta[dq_meta$iso3 != "CHN", ]
}
#hard code french guiana and western sahara
dq_meta[iso3 == "GUF", spr_reg_nm := "Latin America and Caribbean"]
dq_meta[iso3 == "ESH", spr_reg_nm := "North Africa and Middle East"]

total_deaths <- dq_meta[,.(d = sum(d * 1E5, na.rm = T)), by = "year"]
deaths_region <- dq_meta[,.(d = sum(d * 1E5, na.rm = T)), by = c("year", "spr_reg_nm")]
less_than_25 <- dq_meta[q < 25, .(d25 = sum(d * 1E5, na.rm = T)), by = c("year", "spr_reg_nm")]
more_than_100 <- dq_meta[q > 100, .(d100 = sum(d * 1E5, na.rm = T)), by = c("year", "spr_reg_nm")]

deaths_fract <- merge(deaths_region, less_than_25, by = c("year", "spr_reg_nm"))
deaths_fract <- merge(deaths_fract, more_than_100, by = c("year", "spr_reg_nm"))

deaths_fract[, pct_d25 := (d25 / d) * 100]
deaths_fract[, pct_d100 := (d100 / d) * 100]

deaths_fract <- setorder(deaths_fract, spr_reg_nm, year)

##################################################################
##  Superregion distribution plots
##################################################################

## Histograms for first and last year, breaking down by admin2 unit and super-region
# Rename superregion 
dq_meta[spr_reg_nm == "Central Europe, Eastern Europe, and Central Asia", spr_reg_nm := "Eastern Europe and Central Asia"]

dq_meta_first <- dq_meta[year==2000,]
dq_meta_first <- dq_meta_first[Stage != 3,]
dq_meta_first_col <- dq_meta_first[, .(d=sum(d, na.rm = T)), by=c('q_bin', 'spr_reg_nm')]
dq_meta_first_agg <- dq_meta_first[, .(d=sum(d, na.rm = T)), by=c('q_bin')]
dq_meta_last  <- dq_meta[year==2017,]
dq_meta_last <- dq_meta_last[Stage != 3,]
dq_meta_last_col <- dq_meta_last[, .(d=sum(d, na.rm = T)), by=c('q_bin', 'spr_reg_nm')]
dq_meta_last_agg <- dq_meta_last[, .(d=sum(d, na.rm = T)), by=c('q_bin')]

#figure 4 caption number plug
less_than_25_2017 <- dq_meta[q < 25 & year == 2017, .(d25 = sum(d * 1E5, na.rm = T)), by = c("year")]$d25
total_2017 <- dq_meta[year == 2017, .(d25 = sum(d * 1E5, na.rm = T)), by = c("year")]$d25
print(round((less_than_25_2017 / total_2017) * 100, 1))


colors <-rev(c('#1D6996', '#0F8554', '#73AF48', '#E17C05', '#CC503E', "#6F4070"))

fig5_noline_2000 <- ggplot() +
  geom_segment(aes(x = 27.5, xend=27.5, yend=Inf, y = 0)) + 
  geom_bar(data=dq_meta_first_col, aes(x=q_bin, y=d, fill=spr_reg_nm), stat='identity') +
  geom_step(data=dq_meta_last_agg, aes(x=(q_bin - 2.5), y=d), alpha = .3, linetype = 2) +
  labs(title='2000',
       y='Number of Under-5 deaths (100k)',
       fill='GBD Super\nRegion',
       x = "Under-5 Mortality Rate (per 1,000)"
  ) +
  theme_minimal() +
  theme(legend.position='none') +
  scale_fill_manual(values = colors) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_x_continuous(breaks = c(0,27.5,102.5,202.5,302.5), expand=c(0,0), labels = c(0,25,100,200,300))+
  coord_cartesian(xlim = c(0, 315)) +
  annotate("text", x = 275, y = 3.50, label = "--- 2017", colour = "grey35", size = 5) +
  annotate("text", x = 282.7, y = 2.75, label = "SDG 3.2", colour = "grey35", size = 5) +
  annotate("text", x = 282.7, y = 2.25, label = "target", colour = "grey35", size = 5) +
  geom_segment(aes(x = 265.85, xend=270.84, yend=2.50, y = 2.50), color="grey35")
  

fig5_noline_2017 <- ggplot() +
  geom_segment(aes(x = 27.5, xend=27.5, yend=5, y = 0)) + 
  geom_bar(data=dq_meta_last_col, aes(x=q_bin, y=d, fill=spr_reg_nm), stat='identity') +
  geom_step(data=dq_meta_first_agg, aes(x=(q_bin - 2.5), y=d), alpha = .3, linetype = 2) +
  labs(title='2017',
       y='Number of Under-5 deaths (100k)',
       x='Under-5 Mortality Rate (per 1,000)',
       fill='GBD Super\nRegion'
  ) +
  theme_minimal() +
  theme(legend.position='bottom',
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = c(0,27.5,102.5,202.5,302.5), labels = c(0,25,100,200,300), expand=c(0,0)) +
  coord_cartesian(xlim = c(0, 315)) +
  annotate("text", x = 275, y = 3.50, label = "--- 2000", colour = "grey35", size = 5) +
  annotate("text", x = 282.7, y = 2.75, label = "SDG 3.2", colour = "grey35", size = 5) +
  annotate("text", x = 282.7, y = 2.25, label = "target", colour = "grey35", size = 5) +
  geom_segment(aes(x = 265.85, xend=270.84, yend=2.50, y = 2.50), color="grey35")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(fig5_noline_2017)

lay <- rbind(c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,2,2,2,2),
             c(2,2,2,2,2,2,2,2,2),
             c(2,2,2,2,2,2,2,2,2),
             c(4,4,4,4,4,4,4,3,3))

png(filename = figure4_save_path,
    width=12,
    height=8,
    units = "in",
    res = 300)
grid.arrange(fig5_noline_2000, fig5_noline_2017 + theme(legend.position = 'none'), minimap, legend, layout_matrix = lay)
dev.off()

##################################################################
##  Superregion proportional distribution plots
##################################################################

## Histograms for first and last year, breaking down by admin2 unit and super-region
# Rename superregion

dq_meta_first[, d_sum := sum(d, na.rm=T), by = c("q_bin")]
dq_meta_first[, d := (d / d_sum) * 100]
dq_meta_first_col <- dq_meta_first[, .(d_all = sum(d, na.rm=T)), by = c("q_bin", "spr_reg_nm")]
dq_meta_last[, d_sum := sum(d, na.rm=T), by = c("q_bin")]
dq_meta_last[, d := (d / d_sum) * 100]
dq_meta_last_col <- dq_meta_last[, .(d_all = sum(d, na.rm=T)), by = c("q_bin", "spr_reg_nm")]

fig5_noline_2000 <- ggplot() +
  geom_bar(data=dq_meta_first_col, aes(x=q_bin, y=d_all, fill=spr_reg_nm), stat='identity') +
  labs(title='2000',
       y='Percent of Under-5 deaths',
       x='Under-5 Mortality Rate (per 1,000)',
       fill='GBD Super\nRegion'
  ) +
  theme_minimal() +
  theme(legend.position='none') +
  scale_fill_manual(values = colors) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_x_continuous(breaks = c(0,27.5,102.5,202.5,302.5), labels = c(0,25,100,200,300), expand=c(0,0)) +
  coord_cartesian(xlim = c(0, 315)) 


fig5_noline_2017 <- ggplot() +
  geom_bar(data=dq_meta_last_col, aes(x=q_bin, y=d_all, fill=spr_reg_nm), stat='identity') +
  labs(title='2017',
       y='Percent of Under-5 deaths',
       x='Under-5 Mortality Rate (per 1,000)',
       fill='GBD Super\nRegion'
  ) +
  theme_minimal() +
  theme(legend.position='bottom',
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) + 
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = c(0,27.5,102.5,202.5,302.5), labels = c(0,25,100,200,300), expand=c(0,0)) +
  coord_cartesian(xlim = c(0, 315))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

legend <- g_legend(fig5_noline_2017)

lay <- rbind(c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(1,1,1,1,1,1,1,1,1),
             c(2,2,2,2,2,2,2,2,2),
             c(2,2,2,2,2,2,2,2,2),
             c(2,2,2,2,2,2,2,2,2),
             c(4,4,4,4,4,4,4,3,3))

png(filename = figure4_proportion_save_path,
    width=12,
    height=8,
    units = "in",
    res = 300)
grid.arrange(fig5_noline_2000, fig5_noline_2017 + theme(legend.position = 'none'), minimap, legend, layout_matrix = lay)
dev.off()

##################################################################
## Population decomposition plots
##################################################################
pops <- fread(paste0("<<<< FILEPATH REDACTED >>>>"), drop = "V1")

dq_meta[, q_bin := round(q/10,0)*10]


dq_meta_pop <- merge(dq_meta, pops, by.x = c("ADM_CODE", "year"), by.y = c("ADM2_CODE", "year"), all.x =T)
total <- dq_meta_pop[, .(total = sum(pop, na.rm = T)), by = "year"]

adm2_bin <- unique(dq_meta[year == 2000, c("ADM_CODE", "q_bin")])

dq_meta_pop[, q_bin := NULL]
dq_meta_pop <- merge(dq_meta_pop, adm2_bin, by = "ADM_CODE")

dq_meta_pop_agg <- dq_meta_pop[, .(d = sum(d, na.rm=T), pop = sum(pop, na.rm=T)), by = c("q_bin", "year")]

#reshape functions are for losers
dq_meta_pop_2000 <- dq_meta_pop_agg[year == 2000,]
dq_meta_pop_2017 <- dq_meta_pop_agg[year == 2017,]
setnames(dq_meta_pop_2000, c("d", "pop"), c("d_2000", "pop_2000"))
setnames(dq_meta_pop_2017, c("d", "pop"), c("d_2017", "pop_2017"))
dq_meta_pop_2000 <- dq_meta_pop_2000[, !c("year")]
dq_meta_pop_2017 <- dq_meta_pop_2017[, !c("year")]
dq_meta_pop_spread <- merge(dq_meta_pop_2000, dq_meta_pop_2017, by = "q_bin", all = T)
dq_meta_pop_spread <- dq_meta_pop_spread[!is.na(q_bin),]

dq_meta_pop_spread[, d_dis := d_2000 * pop_2017 / pop_2000]

line_colors <- c("2000 deaths" = "#648FFF", "Expected 2017 deaths\nwith population growth"="#DC267F", "2017 deaths"="#FE6100")
fig5_noline_2000 <- ggplot(data = dq_meta_pop_spread) +
  geom_segment(data = dq_meta_pop_spread, aes(x=q_bin-5, xend = q_bin, y=d_2000, yend=d_2000, color="2000 deaths")) +
  geom_segment(data = dq_meta_pop_spread, aes(x=q_bin-5, xend = q_bin+5, y=d_dis, yend=d_dis, color="Expected 2017 deaths\nwith population growth")) +
  geom_segment(data = dq_meta_pop_spread, aes(x=q_bin, xend = q_bin+5, y=d_2017, yend=d_2017, color="2017 deaths")) +
  geom_segment(data = dq_meta_pop_spread, aes(x=q_bin-2.5, xend = q_bin-2.5, y=d_2000, yend=d_dis), arrow=arrow(ends = "last", length = unit(.2, "cm"))) +
  geom_segment(data = dq_meta_pop_spread, aes(x=q_bin+2.5, xend = q_bin+2.5, y=d_dis, yend=d_2017), arrow=arrow(ends = "last", length = unit(.2, "cm"))) +
  labs(y='Number of Under-5 deaths (100k)',
       x = "Under-5 Mortality Rate (per 1,000 LB)") +
  theme_minimal() +
  theme(legend.position='right',
        legend.title = element_blank()) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank())+
  scale_color_manual(values=line_colors, breaks=names(line_colors)) +
  scale_x_continuous(breaks = c(0,27.5,102.5,202.5,302.5, 315), expand=c(0,0), labels = c(0,25,100,200,300, ""))+
  scale_y_continuous(breaks = c(0, 2, 4, 6, 8, 10, 11), expand = c(0,0), labels = c(0, 2, 4, 6, 8, 10, "")) +
  coord_cartesian(xlim = c(0, 315), ylim = c(0, 11))

png(filename = "<<<< FILEPATH REDACTED >>>>",
    width=12,
    height=8,
    units = "in",
    res = 300)
print(fig5_noline_2000)
dev.off()


##################################################################
##  Age distribution plots
##################################################################

pull_age_binned_deaths <- function(run_date, adm_level) {
  infant <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  neonatal <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  under5 <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  
  setnames(neonatal, "value", "neonatal")
  setnames(infant, "value", "infant")
  setnames(under5, "value", "under5")
  
  under1 <- merge(neonatal, infant, by = c("ADM2_CODE", "year"))
  under1[, value := infant - neonatal]
  under1[, age := "Post-Neonatal"]
  child <- merge(infant, under5, by = c("ADM2_CODE", "year"))
  child[, value := under5 - infant]
  child[, age := "Child"]
  
  setnames(neonatal, "neonatal", "value")
  neonatal[,age := "Neonatal"]
  neonatal[,V1 := NULL]
  
  under1 <- under1[, c("ADM2_CODE", "year", "value", "age")]
  child <- child[, c("ADM2_CODE", "year", "value", "age")]
  
  all_ages <- rbind(neonatal, under1, child)
  return(all_ages)
}

d <- pull_age_binned_deaths(run_date, adm_level)
dq <- assemble_q_d(
  run_date = run_date,
  adm_level = adm_level
)
dq[,d:=NULL]
dq <- merge(d, dq, by.x = c("ADM2_CODE", "year"), by.y = c("ADM_CODE", "year"))
setnames(dq, "value", "d")

# Clean and format columns
dq[, year_factor := as.factor(year)]
dq[, q := q * 1000]
dq[, d := d / 1E5]
dq <- dq[order(year,q)]
dq[, q_bin := round(q/5,0)*5]

# For first and last years only, merge on region
# Load tables connecting admin2 units to countries and countries to regions
dq_firstlast <- dq[year %in% c(2000, 2017)]
suffix      <- '.dbf'
hierarchy_path <- paste0(get_admin_shapefile(admin_level = 2, raking = F, version = shapefile_version))
hierarchy <- readOGR(hierarchy_path)
hierarchy <- as.data.table(hierarchy@data)

setnames(hierarchy, paste0('ADM',adm_level,'_CODE'), "ADM_CODE")
# Subset to only subnational admin units and national admin codes
hierarchy <- hierarchy[, .(ADM_CODE, ADM0_CODE)]
hierarchy$ADM_CODE <- as.integer(as.character(hierarchy$ADM_CODE))
hierarchy$ADM0_CODE <- as.integer(as.character(hierarchy$ADM0_CODE))


# - Admin0 to country name and region
by_stage <- fread(
  paste0("<<<< FILEPATH REDACTED >>>>")
)[, .(location_name, iso3, Stage, spr_reg_nm, gadm_geoid, mbg_reg_nm)]
setnames(by_stage, 'gadm_geoid','ADM0_CODE')
# - Combine to get merge table from adm2 to all metadata
all_meta <- merge(
  x   = hierarchy,
  y   = by_stage,
  by  = c('ADM0_CODE')
)

setnames(dq_firstlast, "ADM2_CODE", "ADM_CODE")
all_meta$ADM_CODE <- as.integer(as.character(all_meta$ADM_CODE))
# Merge onto the admin2 table
dq_meta <- merge(
  x     = dq_firstlast,
  y     = all_meta,
  by    = c('ADM_CODE'),
  all.x = TRUE
)

if(drop_china){
  dq_meta <- dq_meta[dq_meta$iso3 != "CHN", ]
}

dq_meta_first <- dq_meta[year==2000,]
dq_meta_first <- dq_meta_first[Stage != 3,]
dq_meta_first_age <- dq_meta_first[, .(d=sum(d, na.rm = T)), by=c('q_bin', "age")]
dq_meta_first_agg <- dq_meta_first[, .(d=sum(d, na.rm = T)), by=c('q_bin')]
dq_meta_last  <- dq_meta[year==2017,]
dq_meta_last <- dq_meta_last[Stage != 3,]
dq_meta_last_age <- dq_meta_last[, .(d=sum(d, na.rm = T)), by=c('q_bin', "age")]
dq_meta_last_agg <- dq_meta_last[, .(d=sum(d, na.rm = T)), by=c('q_bin')]

dq_meta_first_age$age <- factor(dq_meta_first_age$age, levels = c("Child", "Post-Neonatal", "Neonatal"))
dq_meta_last_age$age <- factor(dq_meta_last_age$age, levels = c("Child", "Post-Neonatal", "Neonatal"))

colors <- c("#ff7473", "#ffc952", "#47b8e0")

