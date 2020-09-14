#############################################################################
# Define countries to include in LRI Stg2 manuscript, and assign to regions
#############################################################################

# (1) Setup -----------------------------------------------------------------
rm(list = ls())
country_specific_models <- c('IND','MNG')

# (2) Read in data ----------------------------------------------------------
# diarrhea regions
base_regions <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(iso3, dia_reg)

#read in stage 2
all_stg2 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(Stage %in% c('1', '2a', '2b')) %>%
  select(location_name, gadm_geoid, iso3) %>%
  rename(ADM0_NAME = location_name, ADM0_CODE = gadm_geoid)

#read in input data
data <- fread('<<<< FILEPATH REDACTED >>>>')

# (3) Subset and assign regions ----------------------------------------------------------
countries_with_data <- filter(all_stg2, (iso3 %in% unique(data$country)))
master <- merge(countries_with_data, base_regions, by = 'iso3')

#custom edits:
master$lri_reg <- master$dia_reg

# dia_chn_mng -> MNG country specific model
master[master$iso3 == 'MNG',]$lri_reg <- 'MNG'

# dia_south_asia -> dia_south_asia-IND + IND
master[master$dia_reg == 'dia_south_asia',]$lri_reg <- 'dia_south_asia-IND'
master[master$iso3 == 'IND',]$lri_reg <- 'IND'

# Palestine (PSE) from dia_NA to dia_mid_east
master[master$iso3 == 'PSE',]$lri_reg <- 'dia_mid_east'

# Trinidad and Tobago (TTO) from dia_NA to dia_mid_east
master[master$iso3 == 'TTO',]$lri_reg <- 'dia_s_america'

# Cuba (CUB) from dia_NA to dia_mcaca
master[master$iso3 == 'CUB',]$lri_reg <- 'dia_mcaca-MEX'

# (4) write out/report results ---------------------------------------------------------------------------
write.csv(master, '<<<< FILEPATH REDACTED >>>>')
no_data_excluded <- filter(all_stg2, !(iso3 %in% countries_with_data$iso3))
write.csv(no_data_excluded, '<<<< FILEPATH REDACTED >>>>')

for (region in unique(master$lri_reg)){
  reg_dat <- filter(master, lri_reg == region)
  print(paste0('iso3 codes for ', region, ':'))
  print(paste0(tolower((unique(reg_dat$iso3))), sep="",collapse="+"))
}