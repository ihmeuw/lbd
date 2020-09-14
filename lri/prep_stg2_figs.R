###########################################################################################################
# Prepare Stg2 LRI figures
###########################################################################################################
#prep stg2 estimates script need to launch:
  # *aroc
  # *counterfactual

#(0) Setup ################################################################################################
rm(list = ls())

#libraries
library(RColorBrewer)
library(data.table)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(rgdal)

select <- dplyr::select
rename <- dplyr::rename

#custom functions
source('<<<< FILEPATH REDACTED >>>>/lri/move_to_mapping.R')
source('<<<< FILEPATH REDACTED >>>>/lri/make_inequality.R')
source('<<<< FILEPATH REDACTED >>>>/lri/prep_deaths_averted_figures.R')

#user inputs
run_date <- '2020_06_11_11_19_26' #LRI run date
map_date <- run_date #mapping dir date
year_bounds <- '2000_2019' #start and end modeling years for raster names
end_year <- 2019
no_data_countries <- c('BRA','CHN','CPV','ECU','ESH','GUF','LBY','MYS','VEN')
figure_out_dir <- '<<<< FILEPATH REDACTED >>>>'
si_figure_out_dir <- '<<<< FILEPATH REDACTED >>>>'
region_list <- c('dia_afr_horn', 'dia_cssa', 'dia_wssa', 'dia_name-ESH', 'dia_sssa',
                 'dia_mcaca', 'dia_s_america_n', 'dia_s_america_s', 'dia_central_asia',
                 'dia_se_asia', 'dia_malay', 'dia_south_asia-IND', 'dia_mid_east', 'dia_essa', 'IND','MNG')


# Figure (4) Geographic inequality ########################################################################
make_inequality(run_date = run_date,
                end_year = end_year,
                no_data_countries = no_data_countries,
                out_dir = figure_out_dir,
                pdf_name = 'fig_04_ad2_geographic_inequality.pdf')


# Figure (1) Incidence ####################################################################################
# rate/count/aroc: move to mapping dir
move_to_mapping(map_date = map_date,
                model_date = run_date,
                measures = 'incidence',
                include_aroc = TRUE,
                stats = 'mean',
                admin_levels = 2,
                year_bounds = year_bounds,
                rasters = FALSE,
                include_counts = TRUE)

# absolute deviation from country mean: save to mapping dir
map_dir <- '<<<< FILEPATH REDACTED >>>>'

inc_dev <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, year, abs_lri) %>%
  rename(value = abs_lri) %>%
  write.csv('<<<< FILEPATH REDACTED >>>>')

# Figure (2) Mortality ####################################################################################
# rate/count/aroc: move to mapping dir
move_to_mapping(map_date = map_date,
                model_date = run_date,
                measures = 'mortality',
                include_aroc = TRUE,
                stats = 'mean',
                admin_levels = 2,
                year_bounds = year_bounds,
                rasters = FALSE,
                include_counts = TRUE)

# absolute deviation from country mean: save to mapping dir
map_dir <- '<<<< FILEPATH REDACTED >>>>'

mort_dev <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, year, abs_lri) %>%
  rename(value = abs_lri) %>%
  write.csv('<<<< FILEPATH REDACTED >>>>')

# Figure (3) Country-specific mortality relative deviation ##########################################################
# save to mapping dir
map_dir <- '<<<< FILEPATH REDACTED >>>>'

mort_dev <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, year, rel_lri) %>%
  rename(value = rel_lri) %>%
  write.csv('<<<< FILEPATH REDACTED >>>>')

# Figure (6) Bottom 20% deaths in admins ranked by descending mortality rate, by GBD REGION ############################################################################
# create data frame with indicator for bottom 20%--this frame only includes countries NOT excluded by our mask in the rank
locs <- read.csv('<<<< FILEPATH REDACTED >>>>') %>%
  select(ihme_loc_id, super_region_name)

mort_data_global <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  merge(locs, by.x = 'ISO3', by.y = 'ihme_loc_id') %>%
  select(ADM0_NAME, ADM1_NAME, ADM2_NAME, ADM2_CODE, year, mean, super_region_name) %>%
  filter(year %in% c(2000,end_year)) %>%
  as.data.table() %>%
  dcast(ADM0_NAME + ADM1_NAME + ADM2_NAME + ADM2_CODE + super_region_name ~ year, value.var = 'mean') %>%
  rename(mort_2000 = '2000', mort_end_year = as.character(end_year)) %>%
  as.data.table()

#by region, highlight bottom 20%

for (region in unique(mort_data_global$super_region_name)){
  
  #set region tag for file names
  if (region == 'North Africa and Middle East') region_tag <- 'name'
  if (region == 'Sub-Saharan Africa') region_tag <- 'ssa'
  if (region == 'South Asia') region_tag <- 'sa'
  if (region == 'Latin America and Caribbean') region_tag <- 'lac'
  if (region == 'Southeast Asia, East Asia, and Oceania') region_tag <- 'seao'
  if (region == 'Central Europe, Eastern Europe, and Central Asia') region_tag <- 'ceeeca'
  
  mort_data <- filter(mort_data_global, super_region_name == region) %>%
    as.data.table()
  
  #create rank index where 1 is worst mortality
  mort_data <- mort_data[order(-mort_2000)]
  mort_data$mort_rank_2000 <- c(1:nrow(mort_data))
  
  mort_data <- mort_data[order(-mort_end_year)]
  mort_data$mort_rank_end_year <- c(1:nrow(mort_data))
  
  mort_data <- select(mort_data, ADM2_CODE, mort_rank_end_year, mort_rank_2000)
  
  #load in mort counts
  mort_counts <- fread('<<<< FILEPATH REDACTED >>>>') %>%
    select(ADM0_NAME, ADM1_NAME, ADM2_NAME, ADM2_CODE, year, mean) %>%
    filter(year %in% c(2000,end_year)) %>%
    as.data.table() %>%
    dcast(ADM0_NAME + ADM1_NAME + ADM2_NAME + ADM2_CODE ~ year, value.var = 'mean') %>%
    rename(mort_c_2000 = '2000', mort_c_end_year = as.character(end_year)) %>%
    merge(mort_data, by = 'ADM2_CODE') %>%
    as.data.table()
  
  #get n deaths equal to 20% all deaths in 2000, 2017
  all_deaths_2000 <- sum(mort_counts$mort_c_2000)
  bottom_20_2000 <- 0.2 * all_deaths_2000
  all_deaths_end_year <- sum(mort_counts$mort_c_end_year)
  bottom_20_end_year<- 0.2 * all_deaths_end_year
  
  #take cumulative sum with order set by mort_rank
  mort_counts <- mort_counts[order(mort_rank_2000)]
  mort_counts <- mort_counts[, cum_mort_c_2000 := cumsum(mort_c_2000)]
  
  mort_counts <- mort_counts[order(mort_rank_end_year)]
  mort_counts <- mort_counts[, cum_mort_c_end_year := cumsum(mort_c_end_year)]
  
  mort_counts[, value := ifelse(cum_mort_c_2000 <= bottom_20_2000 & cum_mort_c_end_year > bottom_20_end_year, 1, 99)]
  mort_counts[, value := ifelse(cum_mort_c_2000 > bottom_20_2000 & cum_mort_c_end_year <= bottom_20_end_year, 2, value)]
  mort_counts[, value := ifelse(cum_mort_c_2000 <= bottom_20_2000 & cum_mort_c_end_year <= bottom_20_end_year, 3, value)]
  mort_counts[, value := ifelse(cum_mort_c_2000 > bottom_20_2000 & cum_mort_c_end_year > bottom_20_end_year, 4, value)]
  
  map_dir <- '<<<< FILEPATH REDACTED >>>>'
  mort_counts$year <- 2017
  mort_counts <- select(mort_counts, year, ADM2_CODE, value, ADM0_NAME) %>%
    write.csv('<<<< FILEPATH REDACTED >>>>')
}

# Figure (X) inputs for AROC figure ####################################################################################
# output: data frames at admin 0,1,2 with mort rate (2000,2010,2017) and mort aroc (2000-2017, 2010-2017) by draw

# get mort data
load('<<<< FILEPATH REDACTED >>>>')
mort_0 <- admin_0
mort_1 <- admin_1
mort_2 <- admin_2
rm(admin_0, admin_1, admin_2)

#loop over admins
for (admin in c(0:2)){
  if (admin == 0) admin_idx <- c('ADM0_CODE','ADM0_NAME')
  if (admin == 1) admin_idx <- c('ADM0_CODE','ADM0_NAME','ADM1_CODE','ADM1_NAME')
  if (admin == 2) admin_idx <- c('ADM0_CODE','ADM0_NAME','ADM1_CODE','ADM1_NAME','ADM2_CODE','ADM2_NAME')
  
  #get 2000_2017 and 2010_2017 aroc
  aroc_2000_end_year <- readRDS('<<<< FILEPATH REDACTED >>>>') %>%
    as.data.table() %>%
    rename(V1 = V251)
  aroc_2000_end_year$year <- year_bounds
  aroc_2000_end_year$value <- 'mortality_aroc'
  
  aroc_2010_end_year <- readRDS('<<<< FILEPATH REDACTED >>>>') %>%
    as.data.table() %>%
    rename(V1 = V251)
  aroc_2010_end_year$year <- year_bounds
  aroc_2010_end_year$value <- 'mortality_aroc'
  
  mort <- get(paste0('mort_', admin)) %>%
    filter(year %in% c(2000,2010,end_year)) %>%
    select(-region, -pop)
  mort$value <- 'mortality'
  
  #pull pop, admin idx, and iso3 to merge on later
  pop <- select(get(paste0('mort_', admin)), paste0('ADM', admin, '_CODE'), pop, year)
  admin_info <- fread('<<<< FILEPATH REDACTED >>>>') %>%
    select(admin_idx) %>%
    unique.data.frame()
  iso3 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
    select(iso3, gadm_geoid) %>%
    rename(ADM0_CODE = gadm_geoid)
  
  #merge
  merge <- rbind(aroc_2000_end_year, aroc_2010_end_year) %>%
    rbind(mort, fill = TRUE) %>%
    as.data.frame() %>%
    merge(pop, by = c(paste0('ADM', admin, '_CODE'), 'year'), all.x = TRUE) %>%
    merge(admin_info, by = paste0('ADM', admin, '_CODE')) %>%
    merge(iso3, by = 'ADM0_CODE')
  
  assign(paste0('admin_', admin), merge)
}

#save
save(admin_0, admin_1, admin_2, file = '<<<< FILEPATH REDACTED >>>>')

# Figure (5) Alternate: popualtions left behind ####################################################################################################
# absolute deviation from country mortality rate 2000 (Admin 2)
mort_dev_ad2 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(year == 2000) %>%
  select(ADM0_NAME, ADM2_CODE, abs_lri)

# absolute deviation from country mortality AROC 2000
aroc_ad2 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, mean, ADM0_NAME, ADM0_CODE)

aroc_ad0 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(mean, ADM0_NAME, ADM0_CODE) %>%
  rename(country_mean = mean)

aroc_dev_ad2 <- merge(aroc_ad2, aroc_ad0, by = c('ADM0_NAME', 'ADM0_CODE'))
aroc_dev_ad2[, abs_aroc_dev := country_mean - mean]

left_behind_dev <- merge(mort_dev_ad2, aroc_dev_ad2, by = c('ADM0_NAME', 'ADM2_CODE')) %>%
  filter(abs_lri > 0 & abs_aroc_dev < 0)

# countries with deviation > 0 and AROC < 0 --> get count of pop in these admins
#get pops
load('<<<< FILEPATH REDACTED >>>>')
pop <- select(admin_2, ADM2_CODE, pop, year) %>%
  filter(year == 2017) %>%
  filter(ADM2_CODE %in% left_behind_dev$ADM2_CODE)
sum(pop$pop)

############################################################################################
# SI figures
#
##############################################################################################

# Figure S8, S11, S14: Prevalence 2019 (mean, upper, lower) at raster, admin 1, admin 2 ------------------------------------------------------------------
move_to_mapping(map_date = map_date,
                model_date = run_date,
                measures = 'prevalence',
                include_aroc = FALSE,
                stats = c('mean', 'upper', 'lower'),
                admin_levels = c(1,2),
                year_bounds = year_bounds,
                rasters = TRUE,
                include_counts = FALSE)

# Figure S9, S12, S15: Incidence 2019 (mean, upper, lower) at raster, admin 1, admin 2 -------------------------------------------------------------------
move_to_mapping(map_date = map_date,
                model_date = run_date,
                measures = 'incidence',
                include_aroc = FALSE,
                stats = c('mean', 'upper', 'lower'),
                admin_levels = c(1,2),
                year_bounds = year_bounds,
                rasters = TRUE,
                include_counts = FALSE)

# Figure S10, S13, S16: Mortality 2019 (mean, upper, lower) at raster, admin 1, admin 2 ------------------------------------------------------------------
move_to_mapping(map_date = map_date,
                model_date = run_date,
                measures = 'mortality',
                include_aroc = FALSE,
                stats = c('mean', 'upper', 'lower'),
                admin_levels = c(1,2),
                year_bounds = year_bounds,
                rasters = TRUE,
                include_counts = FALSE)

# Figure S6: finite elements mesh --------------------------------------------------------------------------------------------------------------------------
load('<<<< FILEPATH REDACTED >>>>')
pdf(paste0(si_figure_out_dir, '/si_fig_06_mesh.pdf'))
plot(mesh_s)
dev.off()

# Table S7: xgboost parameters by region -------------------------------------------------------------------------------------------------------------------
xgboost_params <- data.frame()
for (region in region_list){
  reg_data <- fread('<<<< FILEPATH REDACTED >>>>')
  reg_data$region <- region
  reg_data$V1 <- NULL
  xgboost_params <- rbind(xgboost_params, reg_data)
}
write.csv(xgboost_params, '<<<< FILEPATH REDACTED >>>>')


