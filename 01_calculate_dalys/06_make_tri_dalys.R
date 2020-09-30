#####################################################################
# Sum dalys for malaria, LRI, and diarrhea by draw
#####################################################################

#(1) Setup ---------------------------------------------------------
rm(list = ls())

#user inputs
tridaly_run_date <- '2019_10_28'
mal_run_date <- '2019_10_28'

library(dplyr)
library(tidyr)
library(data.table)
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/post_estimation_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core_custom/post_estimation_functions.R')

select <- dplyr::select
rename <- dplyr::rename
count <- dplyr::count


#(2) Load in the admin pred objects ----------------------------------------
load('<<< FILEPATH REDACTED >>>')
lri_draws_0 <- admin_0
lri_draws_1 <- admin_1
lri_draws_2 <- admin_2

load('<<< FILEPATH REDACTED >>>')
dia_draws_0 <- admin_0
dia_draws_1 <- admin_1
dia_draws_2 <- admin_2

load('<<< FILEPATH REDACTED >>>')
mal_draws_0 <- admin_0
mal_draws_1 <- admin_1
mal_draws_2 <- admin_2

rm(admin_0, admin_1, admin_2)

#(3) Bind and sum by draw --------------------------------------------------
for (admin in c(0:2)){
  admin_code <- paste0('ADM', admin, '_CODE')
  lri <- select(get(paste0('lri_draws_', admin)), year, starts_with(admin_code), starts_with('V'))
  dia <- select(get(paste0('dia_draws_', admin)), year, starts_with(admin_code), starts_with('V'))
  mal <- select(get(paste0('mal_draws_', admin)), year, starts_with(admin_code), starts_with('V'))
  pop <- select(get(paste0('lri_draws_', admin)), year, starts_with(admin_code), pop)
  
  test <- bind_rows(lri, dia, mal) %>%
    rename('code' = admin_code) %>%
    group_by(year, code)
  
  incomplete_cases <- count(test, code, year) %>%
    filter(n != 3)
  
  #remove any admin that does not appear in malaria, lri, and diarrhea csvs
  incomplete_admins <- unique(incomplete_cases$code)
  test <- filter(test, !(code %in% incomplete_admins))
  
  #merge on pop and remove draws 101-250, since we only have 100 draws for malaria
  extra_draws <- paste0('V', 101:250)
  sums <- summarise_all(test, sum) %>%
    plyr::rename(c('code' = admin_code)) %>%
    merge(pop, by = c(admin_code, 'year')) %>%
    select(-extra_draws)
  
  assign(paste0('admin_', admin), sums)
}

#(4) Save out results --------------------------------------------------

#save in standard format for admin pred object
save(list = c('admin_0', 'admin_1', 'admin_2', 'sp_hierarchy_list'), file = '<<< FILEPATH REDACTED >>>')

#(5) Summarize admins --------------------------------------------------

summarize_admins_tridalys(ind = 'tridalys',
                          ig= 'tridalys',
                          ad_levels = c(0,1,2), 
                          raked = T, 
                          measure = 'daly', 
                          metrics = 'rates',
                          input_dir =  '<<< FILEPATH REDACTED >>>',
                          output_dir = '<<< FILEPATH REDACTED >>>')

#(6) Create color scale at admin 2 --------------------------------------------------
# Key:
# 0 combined dalys below 0.5/child
# 1 LRI > 60%
# 2 DIA > 60%
# 3 MAL > 60%
# 4 LRI 50-60%
# 5 DIA 50-60%
# 6 MAL 50-60%
# 7 all three causes < 50%
# 8 all three causes 20-40%

#summed results
tridaly_sums <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(year == 2017)

admin_2_list <- unique(tridaly_sums$ADM2_CODE)
extra_draws <- paste0('V', 101:250)

#draw results
load('<<< FILEPATH REDACTED >>>')
tridaly_draws_2 <- admin_2 %>%
  filter(year == 2017 & ADM2_CODE %in% admin_2_list) %>%
  select(-pop, -year)

load('<<< FILEPATH REDACTED >>>')
lri_draws_2 <- admin_2 %>%
  filter(year == 2017 & ADM2_CODE %in% admin_2_list) %>%
  select(-pop, -region, -extra_draws, -year)

load('<<< FILEPATH REDACTED >>>')
dia_draws_2 <- admin_2 %>%
  filter(year == 2017 & ADM2_CODE %in% admin_2_list)%>%
  select(-pop, -region, -extra_draws, -year)

load('<<< FILEPATH REDACTED >>>')
mal_draws_2 <- admin_2 %>%
  filter(year == 2017 & ADM2_CODE %in% admin_2_list)%>%
  select(-pop, -region, -year)

#prop by draw -> take mean
lri_draws_2$cause <- 'lri'
dia_draws_2$cause <- 'dia'
mal_draws_2$cause <- 'mal'
tridaly_draws_2$cause <- 'total'

tridaly_prop <- rbind(lri_draws_2, tridaly_draws_2) %>%
  rbind(dia_draws_2) %>%
  rbind(mal_draws_2) %>%
  select('ADM2_CODE','cause','V1','V2','V3') %>%
  melt(id.vars = c('ADM2_CODE','cause')) %>%
  dcast(ADM2_CODE + variable ~ cause, value.var = 'value') %>%
  mutate(prop_lri = lri / total, prop_mal = mal / total, prop_dia = dia / total) %>%
  group_by(ADM2_CODE) %>%
  summarize(lri_prop_mean = mean(prop_lri), mal_prop_mean = mean(prop_mal), dia_prop_mean = mean(prop_dia))

#merge on total tridalys values
tridaly_prop <- merge(tridaly_prop, select(tridaly_sums, ADM2_CODE, ADM2_NAME, ADM1_CODE, ADM1_NAME, ADM0_CODE, ADM0_NAME, year, mean), by = 'ADM2_CODE')

#assign colors
tridaly_prop <- as.data.table(tridaly_prop)

# 1 LRI > 60%
tridaly_prop[, color_code := ifelse(lri_prop_mean > 0.6, 1, 99999)]

# 2 DIA > 60%
tridaly_prop[, color_code := ifelse(dia_prop_mean > 0.6, 2, color_code)]

# 3 MAL > 60%
tridaly_prop[, color_code := ifelse(mal_prop_mean > 0.6, 3, color_code)]

# 4 LRI 50-60%
tridaly_prop[, color_code := ifelse(lri_prop_mean <= 0.6 & lri_prop_mean > 0.5, 4, color_code)]

# 5 DIA 50-60%
tridaly_prop[, color_code := ifelse(dia_prop_mean <= 0.6 & dia_prop_mean > 0.5, 5, color_code)]

# 6 MAL 50-60%
tridaly_prop[, color_code := ifelse(mal_prop_mean <= 0.6 & mal_prop_mean > 0.5, 6, color_code)]

# 7 all three causes < 50%
tridaly_prop[, color_code := ifelse(lri_prop_mean <= 0.5 & mal_prop_mean <= 0.5 & dia_prop_mean <= 0.5, 7, color_code)]

# 8 all three causes 20-40%
tridaly_prop[, color_code := ifelse(lri_prop_mean <= 0.4 & mal_prop_mean <= 0.4 & dia_prop_mean <= 0.4, 8, color_code)]

# 0 combined dalys below 0.5/child
tridaly_prop[, color_code := ifelse(mean < 0.5, 99, color_code)]

#(7) Write results --------------------------------------------------
write.csv(tridaly_prop, paste0('<<< FILEPATH REDACTED >>>'))