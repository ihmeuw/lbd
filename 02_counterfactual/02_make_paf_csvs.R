##############################################################################
# Combine individual PAF/deaths averted csvs for LRI(wasting/stunting) and 
# DIA(wasting/stunting; water/sanitation) into final csv of results for each cause
############################################################################

# (1) Setup ---------------------------------------------------------

rm(list = ls())

run_date <- '2019_10_28' #for tridalys results
shapefile_version <- '2019_09_10'
admin_level <- 0

pop_file <- '<<< FILEPATH REDACTED >>>'

#packages
library(dplyr)
library(data.table)

#set identifiers by admin_level
if (admin_level == 0) {
  admin_info <- c('ADM0_NAME', 'ADM0_CODE')
  admin_code <- 'ADM0_CODE'
  admin_obj <- 'admin_0'
}

if (admin_level == 1) {
  admin_info <- c('ADM0_NAME', 'ADM0_CODE','ADM1_NAME', 'ADM1_CODE')
  admin_code <- 'ADM1_CODE'
  admin_obj <- 'admin_1'
}

if (admin_level == 2) {
  admin_info <- c('ADM0_NAME', 'ADM0_CODE','ADM1_NAME', 'ADM1_CODE','ADM2_NAME', 'ADM2_CODE')
  admin_code <- 'ADM2_CODE'
  admin_obj <- 'admin_2'
}

#link and stage list
link <- readRDS('<<< FILEPATH REDACTED >>>') %>%
  dplyr::select(admin_info) %>%
  unique.data.frame()
stage_list <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(Stage == 1)

# (2) Merge all diarrhea csvs --------------------------------------

#load results
sani <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(sani_paf17 = paf17, 
         sani_paf00 = paf00, 
         sani_deaths_averted = deaths_averted, 
         dia_deaths_17 = cause_deaths_17, 
         dia_deaths_00 = cause_deaths_00) %>%
  dplyr::select(-V1)

water <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(water_paf17 = paf17, 
         water_paf00 = paf00, 
         water_deaths_averted = deaths_averted) %>%
  dplyr::select(-V1, -cause_deaths_00, -cause_deaths_17)

stunting <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(stunting_paf17 = paf17, 
         stunting_paf00 = paf00, 
         stunting_deaths_averted = deaths_averted) %>%
  dplyr::select(-V1, -cause_deaths_00, -cause_deaths_17)

wasting <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(wasting_paf17 = paf17, 
         wasting_paf00 = paf00, 
         wasting_deaths_averted = deaths_averted) %>%
  dplyr::select(-V1, -cause_deaths_00, -cause_deaths_17)

#bind results
dia_results <- merge(sani, water, by = admin_code) %>%
  merge(stunting, by = admin_code) %>%
  merge(wasting, by = admin_code) %>%
  merge(link, by = admin_code)

# (3) Merge all LRI csvs --------------------------------------------
#load results
stunting <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(stunting_paf17 = paf17, 
         stunting_paf00 = paf00, 
         stunting_deaths_averted = deaths_averted, 
         lri_deaths_17 = cause_deaths_17, 
         lri_deaths_00 = cause_deaths_00) %>%
  dplyr::select(-V1)

wasting <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(wasting_paf17 = paf17, 
         wasting_paf00 = paf00, 
         wasting_deaths_averted = deaths_averted) %>%
  dplyr::select(-V1, -cause_deaths_00, -cause_deaths_17)

underweight <- fread('<<< FILEPATH REDACTED >>>') %>%
  rename(underweight_paf17 = paf17, 
         underweight_paf00 = paf00, 
         underweight_deaths_averted = deaths_averted) %>%
  dplyr::select(-V1, -cause_deaths_00, -cause_deaths_17)

#bind results
lri_results <- merge(stunting, wasting, by = admin_code) %>%
  merge(underweight, by = admin_code) %>%
  merge(link, by = admin_code)

# (4) Vet results

# all admins present?
# pull stage master list
stg1_ad0_codes <- stage_list$gadm_geoid
ad_list <- filter(link, ADM0_CODE %in% stg1_ad0_codes)

#pull pops
load(pop_file)
pops <- get(admin_obj) %>%
  filter(year == 2017) %>%
  dplyr::select(pop, admin_code) %>%
  rename(pop_2017 = pop)

dia_missing <- filter(ad_list, !(get(admin_code) %in% pull(dia_results, get(admin_code)))) %>%
  merge(pops, by = admin_code)
lri_missing <- filter(ad_list, !(get(admin_code) %in% pull(lri_results, get(admin_code)))) %>%
  merge(pops, by = admin_code)

dia_missing
lri_missing

lri_to_find <- filter(lri_missing, pop_2017 != 0)
dia_to_find <- filter(dia_missing, pop_2017 != 0)

table(dia_to_find$ADM0_NAME)
table(lri_to_find$ADM0_NAME)

# check that PAFs bounded between 0 and 1
filter(dia_results, 
       sani_paf17 > 1 | 
         water_paf17 > 1 |
         stunting_paf17 > 1 |
         wasting_paf17 > 1 |
         sani_paf00 > 1 |
         water_paf00 > 1 |
         stunting_paf00 > 1 |
         wasting_paf00 > 1 |
         sani_paf17 < 0 | 
         water_paf17 < 0 |
         stunting_paf17 < 0 |
         wasting_paf17 < 0 |
         sani_paf00 < 0 |
         water_paf00 < 0 |
         stunting_paf00 < 0 |
         wasting_paf00 < 0)

filter(lri_results,
         stunting_paf17 > 1 |
         wasting_paf17 > 1 |
         underweight_paf17 > 1 |
         stunting_paf00 > 1 |
         wasting_paf00 > 1 |
         underweight_paf00 > 1 |
         stunting_paf17 < 0 |
         wasting_paf17 < 0 |
         underweight_paf17 < 0 |
         stunting_paf00 < 0 |
         wasting_paf00 < 0 |
         underweight_paf00 < 0)

#are there 0s or NAs?
sapply(dia_results, function(x) sum(is.na(x)))
sapply(lri_results, function(x) sum(is.na(x)))

sapply(dia_results, function(x) sum(x == 0))
sapply(lri_results, function(x) sum(x == 0))

# (5) Save out results ------------------------------------------------

write.csv(dia_results, '<<< FILEPATH REDACTED >>>')
write.csv(lri_results, '<<< FILEPATH REDACTED >>>')
write.csv(dia_missing, '<<< FILEPATH REDACTED >>>')
write.csv(lri_missing, '<<< FILEPATH REDACTED >>>')


