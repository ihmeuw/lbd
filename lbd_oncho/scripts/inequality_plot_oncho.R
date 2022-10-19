## ###########################################################################
## 
## PLOT HIGHEST AND LOWEST ADM2 UNITS BY COUNTRY IN 2000 and 2018
## 
## Purpose: Plot highest and lowest admin2 units for onchocerciasis mortality by
##   country in 2000 and 2018
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

location_code_type <- 'gadm' # gaul or gadm, Corresponding to run date
adm_level <- 2
first_yr  <- 1988
last_yr   <- 2018
shp_version <- '2020_05_21'
focus <- FALSE


lookup_file <- <<<< FILEPATH REDACTED >>>>
if(location_code_type=='gaul'){
  hierarchy_path <- <<<< FILEPATH REDACTED >>>>
} else {
  hierarchy_path <- <<<< FILEPATH REDACTED >>>>
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

# Merge tables
loc_meta <- merge(
  x  = lookup_table[, .(country_code, spr_reg_nm, iso3)],
  y  = hierarchy[, .(country_code, country_name, adm_code, adm_name)],
  by = c('country_code')
)

ig <- "oncho"
ind <- "had_oncho_w_resamp"
sharedir <- <<<< FILEPATH REDACTED >>>>
rd <- "2020_08_08_15_01_40"

pathadd <- ""

#### Africa draws
input_dir_africa <- <<<< FILEPATH REDACTED >>>>
output_dir_africa <- <<<< FILEPATH REDACTED >>>>
load(<<<< FILEPATH REDACTED >>>>)
adm2_agg_africa <- fread(<<<< FILEPATH REDACTED >>>>)

admin_2[, c("gam", "gbm", "lasso") := NULL]

all_draws <- copy(admin_2)

mydat <- as.data.frame(all_draws)
draws <- mydat[, grep('V', names(mydat))]
mydat$value <- apply(draws, 1, mean, na.rm = TRUE)
mydat <- select(mydat, value, pop, ADM2_CODE, year)

q <- as.data.table(mydat)

setnames(q, c(paste0('ADM',adm_level,'_CODE'),'value'), c('adm_code','q') )

q_meta <- merge(
  x = q,
  y = loc_meta,
  by = c('adm_code'),
  all.x = TRUE
)

q_meta[country_name == "Djibouti", spr_reg_nm := "North Africa and Middle East"]
q_meta[country_name == "Somalia", spr_reg_nm := "North Africa and Middle East"]
q_meta <- q_meta[!(country_name %in% c("Zambia", "Mauritania"))]

# Get highest, median, and lowest adm2 values by country-year
q_stats <- q_meta[, .(q_min = min(q, na.rm=TRUE), 
                      q_med = weighted.mean(x = q, w = pop, na.rm=TRUE), 
                      q_max= max(q, na.rm=TRUE)),
                  by = .(country_code, country_name, iso3, spr_reg_nm, year)]

first_yr <- 2000

## Create data table to show relative inequality (highest/lowest admX) by year
q_ineq <- copy(q_stats)
q_ineq <- q_ineq[ (year >= first_yr) & (year <= last_yr) ,]
q_ineq[, ratio := q_max / q_min ]
q_ineq <- q_ineq[ratio != "Inf"]

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

stat_fy$spr_reg_nm <- "Sub-Saharan Africa"
stat_ly$spr_reg_nm <- "Sub-Saharan Africa"

# Set the plotting order (lowest to highest by country)
stat_ly <- stat_ly[order(spr_reg_nm, q_med)]
stat_ly[,plot_order := .I]

# Merge plotting order back onto the first year of data
stat_fy <- merge(
  x = stat_fy,
  y = stat_ly[, .(iso3, plot_order)],
  by = c('iso3')
)
stat_fy <- stat_fy[order(plot_order)]


## PLOT FIRST GRAPH ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Graph theme: Differences across adminX units within country in 2016

## Create and print plot

# Custom colors for super regions
spr_reg_colors <- c("#e1aa05","#e1aa05","#e1aa05","#e1aa05","#e1aa05","#e1aa05") # standardize on a single color

stat_fy$q_min <- 100*stat_fy$q_min
stat_fy$q_med <- 100*stat_fy$q_med
stat_fy$q_max <- 100*stat_fy$q_max

stat_ly$q_min <- 100*stat_ly$q_min
stat_ly$q_med <- 100*stat_ly$q_med
stat_ly$q_max <- 100*stat_ly$q_max

## Remove some countries

stat_fy <- stat_fy[!(iso3 %in% c("EGY", "OMN", "SAU", "ERI"))]
stat_ly <- stat_ly[!(iso3 %in% c("EGY", "OMN", "SAU", "ERI"))]

stat_fy[country_name == "Democratic Republic of the Congo", country_name := "Dem. Rep. of Congo"]
stat_ly[country_name == "Democratic Republic of the Congo", country_name := "Dem. Rep. of Congo"]

stat_fy[country_name == "Republic of the Congo", country_name := "Rep. of Congo"]
stat_ly[country_name == "Republic of the Congo", country_name := "Rep. of Congo"]

stat_fy[country_name == "Central African Republic", country_name := "Central African Rep."]
stat_ly[country_name == "Central African Republic", country_name := "Central African Rep."]

stat_fy[country_name == "São Tomé and Príncipe", country_name := "São Tomé & Príncipe"]
stat_ly[country_name == "São Tomé and Príncipe", country_name := "São Tomé & Príncipe"]

fig1 <-
  ggplot() +
  geom_crossbar(data=stat_fy, aes(x=reorder(country_name, plot_order), y=q_med, ymin=q_min, ymax=q_max, color='gray80'), size=2.2, width=0) +
  geom_point(data=stat_fy, aes(x=reorder(country_name, plot_order), y=q_med), color='#000000', alpha=0.25, shape=19, size=2) +
  geom_crossbar(data=stat_ly, aes(x=reorder(country_name, plot_order), y=q_med, ymin=q_min, ymax=q_max, color='orange2'), size=2.2, width=0, position=position_nudge(x=0.32)) +
  geom_point(data=stat_ly, aes(x=reorder(country_name, plot_order), y=q_med), color='#000000', alpha=0.5, shape=19, size=2, position=position_nudge(x=0.32)) +
  labs(x = 'Country', y = 'Onchocerciasis Prevalence (%)', color = 'GBD Super Region') +
  theme_bw() + 
  theme(legend.position = "right", text = element_text(size=12), axis.text.x = element_text(size=12, angle = 90, hjust = 0.95)) +
  scale_color_identity(guide = "legend", name = "Year", breaks = c("gray80", "orange2"), labels = c("2000", "2018"))

fig_lf <- fig1
fig_lf2 <- fig_lf

out_dir <- <<<< FILEPATH REDACTED >>>>

png(
  <<<< FILEPATH REDACTED >>>>,
  height=8, width=18, units = 'in', res = 1200
)
print(fig_lf2)
dev.off()


## Other experimental plots
q_meta <- merge(q_meta, stat_fy[, c("iso3", "q_med", "q_min", "q_max")], by = "iso3", all.x = TRUE)
setnames(q_meta, c("q_med", "q_min", "q_max"), c("q_med_fy", "q_min_fy", "q_max_fy"))
q_meta <- merge(q_meta, stat_ly[, c("iso3", "q_med", "q_min", "q_max", "plot_order")], by = "iso3", all.x = TRUE)
setnames(q_meta, c("q_med", "q_min", "q_max"), c("q_med_ly", "q_min_ly", "q_max_ly"))
q_meta <- q_meta[!(iso3 %in% c("OMN", "SAU", "EGY", "ERI"))]
q_meta$q <- q_meta$q * 100

ggplot() + theme_bw() +
  geom_point(data=q_meta[year == first_yr], aes(x=reorder(country_name, plot_order), y=q), color="#999999", alpha=0.5, shape=19, size=2, position=position_nudge(x=-.16)) +
  geom_point(data=q_meta[year == first_yr], aes(x=reorder(country_name, plot_order), y=q_med_fy), color="#999999", alpha=0.25, shape=19, size=2, position=position_nudge(x=-.16)) +
  geom_point(data=q_meta[year == last_yr], aes(x=reorder(country_name, plot_order), y=q), color="#e1aa05", alpha=0.5, shape=19, size=2, position=position_nudge(x=.16)) +
  geom_point(data=q_meta[year == last_yr], aes(x=reorder(country_name, plot_order), y=q_med_ly), color="#000000", alpha=0.25, shape=19, size=2, position=position_nudge(x=.16)) +
  labs(x = 'Country', y = 'Onchocerciasis Prevalence (%)', color = 'GBD Super Region') +
  scale_color_manual(values="") +
  theme(legend.position = 'none', text = element_text(size=12), axis.text.x = element_text(size=12, angle = 90, hjust = 0.95))

ggplot() + theme_bw() +
  geom_boxplot(data = q_meta[year == first_yr], aes(x=reorder(country_name, plot_order), y=q), color="#999999", alpha=0.5, shape=19, size=0.5, width = 0.2, position=position_nudge(x=-.16)) +
  geom_boxplot(data = q_meta[year == last_yr], aes(x=reorder(country_name, plot_order), y=q), color="#e1aa05", alpha=0.5, shape=19, size=0.5, width = 0.2, position=position_nudge(x=.16)) +
  labs(x = 'Country', y = 'Onchocerciasis Prevalence (%)', color = 'GBD Super Region') +
  scale_color_manual(values="") +
  theme(legend.position = 'none', text = element_text(size=12), axis.text.x = element_text(size=12, angle = 90, hjust = 0.95, vjust=0.5))
