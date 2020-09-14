#################################################################################################################################################################
# Calculates deaths averted for LRI due to total air polution (TAP)
##################################################################################################################################################################

# (1) Setup -----------------------------------------------------------------------------------------------------------------------------------------------------

rm(list = ls())
library(data.table)
library(dplyr)

#user arguments
lri_run_date <- commandArgs()[4]
tap_run_date <- commandArgs()[5]

#sample data set for formatting purposes
save_dir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(save_dir)

#additional arguments for this risk
admin_levels <- c(0,2)

#set directories
in_dir <- '<<<< FILEPATH REDACTED >>>>'

# (2) Read in and bind files by admin -------------------------------------------------------------------------------------------------------------------------------------
for (admin_level in admin_levels){
  
  #set identifiers based on admin level
  if (admin_level == 2){
    admin_code <- 'ADM2_CODE'
    admin_obj <- 'admin_2'
  }
  
  if (admin_level == 1) {
    admin_code <- 'ADM1_CODE'
    admin_obj <- 'admin_1'
  }
  
  if (admin_level == 0) {
    admin_code <- 'ADM0_CODE'
    admin_obj <- 'admin_0'
  }
  
  #read in cause data 
  cause_data_lbd <- fread('<<<< FILEPATH REDACTED >>>>')
  
  admin_files <- list.files(in_dir, pattern = paste0('_ad', admin_level, '_'), full.names = TRUE)
  tap_pafs <- list()
  for (file in admin_files){
    reg_file <- fread(file)
    print(paste0('starting ', file))
    
    if ('hap' %in% names(reg_file)) reg_file <- select(reg_file, -hap)
    if ('lri' %in% names(reg_file)) reg_file <- select(reg_file, -lri)
    if('ADM2_CODE' %in% names(reg_file)) reg_file <- select(reg_file, -ADM0_CODE)
    
    tap_pafs <- rbind(tap_pafs, reg_file)
    print(paste0('success loading ', file))
  }
  
  # (3) Calculate deaths averted from 2000 and 2017 attributable deaths -----------------------------------------------------------------------------------------------------
  #subset to start and end years, then cast by year and clean up for consistency with other risks scripts
  
  if (admin_level == 0){
    tap_pafs <- filter(tap_pafs, year %in% c(2000, 2017)) %>%
      select(-iso3, -location_id) %>%
      melt(id = c(admin_code,'year')) %>%
      dcast(ADM0_CODE ~ variable + year, value.var = 'value') %>%
      as.data.table() %>%
      select(admin_code, tap_paf_2017, tap_paf_2000) %>%
      rename(paf17 = tap_paf_2017, paf00 = tap_paf_2000)
  } else if (admin_level ==2){
    tap_pafs <- filter(tap_pafs, year %in% c(2000, 2017)) %>%
      select(-iso3, -location_id) %>%
      melt(id = c(admin_code,'year')) %>%
      dcast(ADM2_CODE  ~ variable + year, value.var = 'value') %>%
      as.data.table() %>%
      select(admin_code, tap_paf_2017, tap_paf_2000) %>%
      rename(paf17 = tap_paf_2017, paf00 = tap_paf_2000)
    }
  
  #merge on LRI deaths data
  cause_deaths_17 <- cause_data_lbd %>%
    as.data.frame() %>%
    dplyr::filter(year == 2017) %>%
    dplyr::select(mean, admin_code) %>%
    rename(cause_deaths_17 = mean)
  
  cause_deaths_00 <- cause_data_lbd %>%
    as.data.frame() %>%
    dplyr::filter(year == 2000) %>%
    dplyr::select(mean, admin_code) %>%
    rename(cause_deaths_00 = mean)
  
  tap_pafs <- merge(tap_pafs, cause_deaths_17, by = admin_code) %>%
    merge(cause_deaths_00, by = admin_code) %>% 
    mutate(deaths_averted = cause_deaths_17*((1-paf17)/(1-paf00)*paf00-paf17))
  
  write.csv(tap_pafs, paste0(save_dir, 'has_lri_tap_paf_ad', admin_level,'.csv'))
  
}

