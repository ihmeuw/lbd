#Taly tridalys countries
rm(list = ls())
library(dplyr)

#set run dates
tridalys_run_date <- '2019_10_28'

#set countries we will force a mask on
dalys_exclude <- c('South Africa','Algeria','Egypt','Libya','Morocco','Tunisia','Western Sahara') #team decision to exlude SA due to unstable GBD diarrhea estimates, as well as North Africa
dia_exclude <- c('Cape Verde','Libya') #Stg 2 diarrhea exclusions: no data

#link all stage 2 countries
all_stg1 <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(Stage == 1) %>%
  select(location_name, gadm_geoid) %>%
  rename(ADM0_NAME = location_name, ADM0_CODE = gadm_geoid)
  
#countries included in analysis (csvs)
mal_results <- fread('<<< FILEPATH REDACTED >>>') %>%
  select(ADM0_CODE, ADM0_NAME) %>%
  unique.data.frame()

mal_no_results <- filter(all_stg1, !(ADM0_CODE %in% mal_results$ADM0_CODE))

lri_results <- fread('<<< FILEPATH REDACTED >>>') %>%
  select(ADM0_CODE, ADM0_NAME) %>%
  unique.data.frame()

lri_no_results <- filter(all_stg1, !(ADM0_CODE %in% lri_results$ADM0_CODE))

dia_results <- fread('<<< FILEPATH REDACTED >>>') %>%
  select(ADM0_CODE, ADM0_NAME) %>%
  unique.data.frame() %>%
  filter(ADM0_CODE %in% all_stg1$ADM0_CODE)

dia_no_results <- filter(all_stg1, !(ADM0_CODE %in% dia_results$ADM0_CODE))

#countries with NA results
mal_na<- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(is.na(mean_rate)) %>%
  select(ADM0_CODE, ADM0_NAME) %>%
  unique.data.frame()

dia_na<- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(is.na(mean_rate)) %>%
  select(ADM0_CODE, ADM0_NAME) %>%
  unique.data.frame()

lri_na<- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(is.na(mean_rate)) %>%
  select(ADM0_CODE, ADM0_NAME) %>%
  unique.data.frame()

#lri no data
lri_input <- fread('<<< FILEPATH REDACTED >>>')
lri_data <- unique(lri_input$country)

stg1_country_list <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(Stage == 1)

stg1_country <- unique(stg1_country_list$iso3)
  
lri_exclude_table <- filter(stg1_country_list, !(iso3 %in% lri_data)) %>%
  select(location_name, gadm_geoid) %>%
  rename(ADM0_NAME = location_name, ADM0_CODE = gadm_geoid)

lri_exclude <- lri_exclude_table$ADM0_NAME

mask_list <- unique(c(dalys_exclude, dia_exclude, lri_exclude, 
               mal_no_results$ADM0_NAME, lri_no_results$ADM0_NAME, dia_no_results$ADM0_NAME, 
               mal_na$ADM0_NAME, dia_na$ADM0_NAME, lri_na$ADM0_NAME))

mask_list <- data.table(mask_list) %>%
  rename(ADM0_NAME = mask_list)

include_list <- filter(stg1_country_list, !(location_name %in% mask_list$ADM0_NAME)) %>%
  filter(location_name != 'Sao Tome and Principe')

countries_in <- data.table(include_list$location_name)

write.csv(mask_list, '<<< FILEPATH REDACTED >>>')
write.csv(countries_in, '<<< FILEPATH REDACTED >>>')
