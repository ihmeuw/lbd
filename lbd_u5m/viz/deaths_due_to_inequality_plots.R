## ###########################################################################
## 
## PLOT DEATHS DUE TO INEQUALITY
## 
## Purpose: Plot highest and lowest admin2 units for under5 mortality by 
##   country in 2000 and 2017
## 
## ###########################################################################

rm(list=ls())

## Imports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(colorspace)
library(data.table)
library(dplyr)
library(foreign)
library(ggplot2)
library(ggrepel)
library(grid)
library(gridExtra)
library(raster)
library(rgdal)
library(RColorBrewer)
library(viridisLite)

repo_path <- '<<<< FILEPATH REDACTED >>>>'
# Source function for loading population
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
# Source function for getting admin units
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
source(paste0(repo_path,'<<<< FILEPATH REDACTED >>>>'))
# Source function for getting GBD population
source("<<<< FILEPATH REDACTED >>>>")


## Define inputs and filepaths ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_date  <- "<<<< REDACTED >>>>"
modeling_shapefile_version <- '<<<< REDACTED >>>>' # Corresponding to model run date
adm_level <- 2
first_yr  <- 2000
last_yr   <- 2017

in_dir    <- paste0('<<<< FILEPATH REDACTED >>>>',run_date,'/')
out_dir   <- paste0(
  '<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/'
)
lookup_file <- '<<<< FILEPATH REDACTED >>>>'
dir.create(out_dir, showWarnings = FALSE)

S2_BOUNDS <- list(
  'lat_min'  = -57,
  'lat_max'  = 63,
  'long_min' = -129,
  'long_max' = 165
)

## Load input data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load D for all adminX units
d <- fread(paste0(in_dir,'<<<< FILEPATH REDACTED >>>>',adm_level,'<<<< FILEPATH REDACTED >>>>'))
# Load population data for all years
template_rast <- brick(
  paste0(in_dir,'<<<< FILEPATH REDACTED >>>>',first_yr,'_',last_yr,'.tif')
)
pop_rast <- load_and_crop_covariates_annual(
  covs           = 'worldpop',
  measures       = 'a0004t',
  simple_polygon = template_rast,
  start_year     = first_yr,
  end_year       = last_yr,
  interval_mo    = as.numeric(12),
  agebin=1
)[[1]]
# Load shapefile and raster for admin2 units
adx_spatial <- GetAdmin(
  admin_level   = adm_level,
  simple_raster = template_rast
)
adx_shapefile <- adx_spatial[['spdf']]
adx_raster    <- adx_spatial[['rast']]
# Load lookup table matching countries to stages
lookup_table <- fread(lookup_file)


## Create dataset of all deaths due to inequality in the mortality rate ~~~~~~

## ** Get unique data.table of all adminX units in Stage 2 countries **
# Pre-formatting for lookup table
shapefile_type <- detect_adm_shapefile_date_type()$shpfile_type
if(shapefile_type=='gaul') lt_code <- 'GAUL_CODE' else lt_code <- 'gadm_geoid'
setnames(
  lookup_table,
  c('location_name',lt_code),
  c('country_name','country_code')
)
lookup_table <- lookup_table[, .(country_name, country_code, spr_reg_nm, loc_id, iso3)]
# Make one update in the lookup table for GBD Super-Region
lookup_table[ iso3 == 'URY', spr_reg_nm := "Latin America and Caribbean"]
# Pre-formatting for admx data
adx_data <- as.data.table(adx_shapefile@data)
setnames(
  adx_data, 
  c(paste0('ADM',adm_level,c('_CODE','_NAME')),'ADM0_CODE','ADM0_NAME'),
  c('admin_code','admin_name','country_code','country_name')
)
adx_data <- adx_data[, .(admin_code,admin_name,country_code) ]
# Merge admin data with stage 2 data, keeping only overlapping countries
adx_data <- merge(
  x  = adx_data,
  y  = lookup_table,
  by = c('country_code')
)

## ** Aggregate population by admin2 unit **
# Ensure that the population and admin rasters have the same shape
if( any(dim(pop_rast)[1:2] != dim(adx_raster)[1:2]) ) stop("Issue with raster dimensions")
pop_dt <- data.table(
  pop        = as.vector(pop_rast),
  admin_code = as.vector(adx_raster),
  year       = rep(first_yr:last_yr, each=prod(dim(pop_rast)[1:2]))
)
# Remove NA values
pop_dt <- na.omit(pop_dt)
message(paste("After dropping NA values",nrow(pop_dt),"values remain in pop_dt."))
# Aggregate
pop_dt <- pop_dt[, .(pop = sum(pop)), by=.(admin_code, year)]

## ** Merge deaths and population for all adminX units **
# Pre-formatting for d
setnames(d, c(paste0("ADM",adm_level,"_CODE"),'mean'), c('admin_code','d'))
d <- d[(year >= first_yr) & (year <= last_yr), .(admin_code,year,d) ]
# Merge
full_data <- merge(
  x = adx_data,
  y = d,
  by = c('admin_code')
)
full_data <- merge(
  x = full_data,
  y = pop_dt,
  by = c('admin_code','year')
)
full_data <- na.omit(full_data[pop > 0,])
message(paste("After subsetting and prepping data,",nrow(full_data),"units remain."))

## ** Get lowest admX-level U5MR by country and year **
# Create U5MR column
full_data[, mort_rate := d / pop]
# Get lowest U5MR by country-year
full_data[, min_mort_rate := min(mort_rate), by=.(country_code, year)]
full_data[, min_mort_admin := 0]
full_data[mort_rate == min_mort_rate, min_mort_admin := 1]

# Get mortality difference and lives lost due to in-country inequality by
#  country-year
full_data[, mort_diff := mort_rate - min_mort_rate]
full_data[, d_counterfact := pop * min_mort_rate ] # Counterfactual deaths
full_data[, d_ineq := d - d_counterfact ] # Deaths attributable to inequality
# Correct for rounding errors
full_data[d_ineq < 0, d_ineq := 0]


## Reformat deaths due to inequality and prep for plotting ~~~~~~~~~~~~~~~~~~~

## ** Prep plot 1: Aggregate deaths due to inequality by country-year **

ineq_by_country <- full_data[, .(d = sum(d),
                                 d_counterfact = sum(d_counterfact), 
                                 d_ineq = sum(d_ineq)),
                             by=.(country_name, spr_reg_nm, iso3, year, loc_id)]
ineq_by_country_ly <- ineq_by_country[ year == last_yr, ]
# Get estimated GBD population by country in 2017
gbd_pop <- get_population(
  location_id  = unique( ineq_by_country_ly[,loc_id] ),
  sex_id       = 3,
  age_group_id = 1,
  year_id      = 2017,
  gbd_round_id = 5
) %>% as.data.table
setnames(gbd_pop, c('population','location_id'), c('gbd_pop','loc_id'))
gbd_pop <- gbd_pop[, .(gbd_pop, loc_id)]
# Merge GBD populations onto final year of inequality data  
ineq_by_country_ly <- merge(
  x     = ineq_by_country_ly,
  y     = gbd_pop,
  by    = c('loc_id'),
  all.x = TRUE
)
ineq_by_country_ly <- ineq_by_country_ly[order(-d)]
ineq_by_country_ly[, deaths_rank := .I]
# Keep only the top N countries
n_countries <- 25
top_countries <- ineq_by_country_ly[deaths_rank <= n_countries, ]
other_countries <- ineq_by_country_ly[
  deaths_rank > n_countries,
  .(d=sum(d), d_counterfact=sum(d_counterfact, na.rm=TRUE), d_ineq=sum(d_ineq)),
  by=.(year)
]
other_countries[, country_name := 'All other LMICs' ]
other_countries[, iso3 := 'OTHERS' ]
other_countries[, deaths_rank := n_countries + 1]

ineq_ly_top_n <- rbindlist(
  list(top_countries, other_countries), use.names=TRUE, fill=TRUE
)
ineq_ly_top_n[, d_cf_plot := d_counterfact * -1 ] # For plotting
# Shorten DRC name for conciseness
ineq_ly_top_n[
  country_name == "Democratic Republic of the Congo",
  country_name := "Dem. Repub.\nof the Congo"
]

## ** Prep plot 2: Merge deaths due to inequality with spatial data **

# Pre-prep inequality data
ineq_sub <- full_data[year==last_yr,
                      .(admin_code, min_mort_admin, d_counterfact, d_ineq, iso3)]
ineq_sub <- ineq_sub[, admin_code := as.integer(admin_code)]
# Pre-prep shapefile
adx_shapefile$admin_code <- as.integer(
  adx_shapefile@data[,paste0('ADM',adm_level,'_CODE')]
)
adx_shapefile$id <- sapply(adx_shapefile@polygons, function(x) x@ID)
adx_shapefile$id <- as.integer(adx_shapefile$id)
# Keep only relevant columns to speed up fortify
shp_keep_cols <- c('id','admin_code','NAME_0',paste0('NAME_',adm_level))
adx_shapefile@data <- adx_shapefile@data[, shp_keep_cols]
# 'Fortify' shapefile for mapping
adx_fort <- fortify(adx_shapefile) %>% as.data.table
adx_fort[, id := as.numeric(id)]
setnames(adx_fort,'order','shp_order')
adx_fort <- merge(
  x     = adx_fort,
  y     = adx_shapefile[, shp_keep_cols],
  by    = c('id'),
  all.x = TRUE
)
# Merge on inequality deaths data, keeping only matching adminX codes
adx_fort <- merge(
  x  = adx_fort,
  y  = ineq_sub,
  by = c('admin_code')
)
# Sort by shapefile plotting order
adx_fort <- adx_fort[order(shp_order),]




## PLOT 1: Deaths due to inequality in 2017 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

abs_comma <- function (x, ...) {
  format(abs(x), ..., big.mark = ",", scientific = FALSE, trim = TRUE)
}

# Add a dummy column for plotting
ineq_ly_top_n[, dummy_ineq := 'ineq']
ineq_ly_top_n[, dummy_cf := 'cf']
fig1_cols <- c('ineq'='#FF6666', 'cf'='#6666FF')
# Define ggplot object
fig1 <- ggplot(data = ineq_ly_top_n,
               aes(x=reorder(country_name,deaths_rank))
               ) +
  geom_bar(
    width=.75, stat='identity',
    aes(y=d_ineq, fill=dummy_ineq)
  ) +
  geom_bar(
    width=.75, stat='identity',
    aes(y=d_cf_plot, fill=dummy_cf)
  ) +
  labs(
    title = paste0(''),
    x = 'Country',
    y = 'Total under-5 deaths in 2017',
    fill = 'Deaths Attributable:  '
  ) +
  scale_y_continuous(
    breaks = (-5:9)*1E5,
    labels = abs_comma,
    minor_breaks = FALSE
  ) +
  scale_fill_manual(
    values = fig1_cols,
    breaks = c('cf','ineq'),
    labels = c(
      'Counterfactual (Lowest Administrative Unit)  ',
      'Attributable to Inequality'
    )
  ) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = 'bottom'
  )

# Plot ggplot object
pdf(paste0(out_dir,'<<<< FILEPATH REDACTED >>>>',adm_level,'.pdf'), height=8, width=15)
print(fig1)
dev.off()


## PLOT 2: Deaths due to inequality in 2017 (map) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load an admin0 map for national boundaries
country_shp <- readOGR( dsn=get_admin_shapefile(admin_level=0) )
country_fort <- fortify(country_shp) %>% as.data.table
setnames(country_fort, 'order', 'shp_order')
country_fort[lat  < S2_BOUNDS$lat_min,  lat  := S2_BOUNDS$lat_min  ]
country_fort[lat  > S2_BOUNDS$lat_max,  lat  := S2_BOUNDS$lat_max  ]
country_fort[long < S2_BOUNDS$long_min, long := S2_BOUNDS$long_min ]
country_fort[long > S2_BOUNDS$long_max, long := S2_BOUNDS$long_max ]
country_fort <- country_fort[order(shp_order),]


# Create ggplot object
adx_map <- copy(adx_fort)
# Set color breaks
ad2_color_breaks <- c(50,500,2000,5000)
col_min <- min(ad2_color_breaks)
col_max <- max(ad2_color_breaks)
adx_map[ d_ineq < col_min, d_ineq := col_min ]
adx_map[ d_ineq > col_max, d_ineq := col_max ]
# Extract points for all reference admin units
reference_units <- copy( adx_map[min_mort_admin==1, ] )
reference_units <- reference_units[ ,
  .(lat=mean(lat, na.rm=TRUE), long=mean(long,na.rm=TRUE)),
  by = .(admin_code, NAME_0, NAME_2)
]
# Add a dummy variable for legend plotting
reference_units[, dummy_size := 'ref' ]
fig2_sizes <- c('ref'=2)

# Add labels for select countries
countries_with_labels <- data.table(
  NAME_0 = c(
    'China','India','South Africa','Indonesia','Brazil',
    'Democratic Republic of the Congo','Peru','Tanzania'),
  label_lat = c(
    22.036123, 1.709988, -32.333954, -13.701123, -25.980011, 
    -9.583769, -12.023858, -4.681814
  ),
  label_long = c(
    119.959396, 75.439715, 32.672237, 95.341121, -42.740127,
    -8.052938, -94.063806, 43.311687
  )
)
label_units <- merge(
  x = reference_units,
  y = countries_with_labels,
  by = c("NAME_0")
)
setnames(label_units, paste0("NAME_",adm_level), 'adm_name')
# String formatting for labels
label_units[, adm_name := as.character(adm_name)]
label_units[, ville := as.numeric(endsWith(adm_name, ' (ville)'))]
label_units[ ville==1, adm_name := substr(adm_name, 1, nchar(adm_name)-8)]
label_units[, adm_name := paste0("  ",adm_name)]
# Reshape points long by country
label_units_long <- data.table(
  lat = c(label_units$lat, label_units$label_lat),
  long = c(label_units$long, label_units$label_long),
  adm0 = rep(label_units$NAME_0, 2)
)

xlim <- c(-108, 157.5)
ylim <- c(-35.5,  53)


fig2 <- ggplot() + 
  geom_polygon(
    data = adx_map,
    aes(x=long, y=lat, group=group, fill=d_ineq), color=NA
  ) + 
  geom_polygon(
    data = country_fort,
    aes(x=long, y=lat, group=group),
    fill=NA, color='#222222',
    lwd=.2
  ) + 
  geom_point(
    data = reference_units,
    aes(x=long, y=lat, size=dummy_size),
    shape = 23,
    fill = "#00e673",
    color= "#006633"
  ) +
  labs(
    title = '',
    subtitle = '',
    fill = 'Attributable\nDeaths',
    size = ''
  ) +
  theme_map() + 
  coord_map(
    projection = "mollweide",
    xlim = xlim, 
    ylim = ylim
  ) + 
  scale_fill_gradientn(
    limits = c(col_min, col_max),
    colors = c("#f7f7c0","#ec5f5f","#6d2a80","#400040"),
    breaks = ad2_color_breaks,
    labels = c('50','500','2k','5k+'),
    trans = 'log10',
    guide = guide_colorbar(order=1)
  ) +
  scale_size_manual(
    values = fig2_sizes,
    breaks = c('ref'),
    labels = rep("Reference\nPoint"),
    guide = guide_legend(order = 2)
  ) +
  theme(
    legend.position = c(0.95, 0.7),
    legend.title = element_text(size=9, color='#222222'),
    legend.text  = element_text(size=6, color='#222222'),
    legend.spacing.y = unit(2,'pt')
  )

# # Save map to PDF
# pdf(
#   paste0(out_dir,'deaths_from_ineq_map.pdf'), 
#   height=5.5, width=12
# )
# print(fig2)
# dev.off()

# Save map to PNG
png(
  paste0(out_dir,'<<<< FILEPATH REDACTED >>>>'), 
  height=5.5, width=12, units='in', res=1000
)
print(fig2)
dev.off()
