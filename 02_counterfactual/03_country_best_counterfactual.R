###########################################################################################################
# Create counterfactual for tridalys paper: how many DALYs averted if every admin 2 below the country mean
# had daly values equal to the country mean for lri, diarrhea, malaria
# by draw
###########################################################################################################

# (1) Setup ##############################################################################################
rm(list = ls())

# user intputs
run_date <- '2019_10_28' #tridalys
mal_run_date <- '2019_10_28' # for malaria admin draw files

# set directories
map_dir <- '<<< FILEPATH REDACTED >>>'
in_dir <- '<<< FILEPATH REDACTED >>>'
mal_in_dir <- '<<< FILEPATH REDACTED >>>'

# (2) Loop over lri, diarrhea, malaria to format draws long and calculate counterfactuals ###############
causes <- c('lri','diarrhea','malaria')

for (cause in causes){
  if (cause == 'lri') indicator <- 'has_lri'
  if (cause == 'diarrhea') indicator <- 'had_diarrhea'
  if (cause == 'malaria') indicator <- 'had_malaria'
  
  if (cause != 'malaria') load(paste0(in_dir, indicator, '_raked_daly_admin_draws_eb_bin0_0.RData'))
  if (cause == 'malaria') load(paste0(mal_in_dir, indicator, '_raked_daly_admin_draws_eb_bin0_0.RData'))
  
  dalys_rate <- admin_2 %>%
    filter(year == 2017) %>%
    merge(sp_hierarchy_list, by = c('ADM2_CODE', 'region')) 
  
  draw_cols <- grep(names(dalys_rate), pattern = 'V', value = TRUE)
  
  min_dalys <- group_by(dalys_rate, ADM0_CODE) %>%
    summarise_at(vars(draw_cols), min) %>%
    melt(id.vars = 'ADM0_CODE') %>%
    rename(draw = variable, counterfactual = value) 
  
  dalys_rate <- melt(dalys_rate, id.vars = c('ADM2_CODE', 'region', 'year', 'ADM0_CODE', 'ADM1_CODE', 'pop', 'ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME')) %>%
    rename(draw = variable, actual = value) %>%
    merge(min_dalys, by = c('ADM0_CODE', 'draw')) %>%
    as.data.table()
  
  # avertable death rate = actual rate - counterfactual rate
  dalys_rate <- dalys_rate[, avertable := actual-counterfactual]
  
  output <- setnames(dalys_rate, old = c('avertable','actual','counterfactual'), new = c(paste0(cause, '_avertable'), paste0(cause, '_actual'), paste0(cause, '_counterfactual')))
  
  assign(paste0(cause, '_dalys'), output)
  
  rm(admin_0, admin_1, admin_2, dalys_rate, sp_hierarchy_list, output, min_dalys)
}


# (4) Merge and assign proportions ################################################################################
malaria_dalys <- select(malaria_dalys, -pop, -region, -ADM1_CODE, -ADM1_NAME, -ADM2_NAME, -ADM0_CODE, -ADM0_NAME)
diarrhea_dalys <- select(diarrhea_dalys, -pop, -region)

dalys <- merge(lri_dalys, diarrhea_dalys, by = c('ADM2_CODE', 'year', 'ADM0_CODE', 'ADM1_CODE', 'ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME','draw')) %>%
  merge(malaria_dalys, by = c('ADM2_CODE', 'year','draw')) %>%
  as.data.table()

dalys <- dalys[, total_avertable_rate := lri_avertable + diarrhea_avertable + malaria_avertable]
dalys <- dalys[, total_avertable_count := total_avertable_rate * pop]

#summarize across draws
dalys <- select(dalys, ADM2_CODE, year, ADM0_NAME, draw, lri_avertable, diarrhea_avertable, malaria_avertable, total_avertable_rate, total_avertable_count) %>%
  group_by(ADM2_CODE, year, ADM0_NAME) %>%
  summarize_at(vars(lri_avertable, diarrhea_avertable, malaria_avertable, total_avertable_rate, total_avertable_count), mean) %>%
  as.data.table()

#color assignments
# lri > 60% 1
dalys <- dalys[, value := ifelse(lri_avertable / total_avertable_rate >= 0.6, 1, 9999)]

# dia > 60% 2
dalys <- dalys[, value := ifelse(diarrhea_avertable / total_avertable_rate >= 0.6, 2, value)]

# mal > 60% 3
dalys <- dalys[, value := ifelse(malaria_avertable / total_avertable_rate >= 0.6, 3, value)]

# lri 50-60% 4
dalys <- dalys[, value := ifelse(lri_avertable / total_avertable_rate < 0.6 & lri_avertable / total_avertable_rate >= 0.5, 4, value)]

# dia 50-60% 5
dalys <- dalys[, value := ifelse(diarrhea_avertable / total_avertable_rate < 0.6 & diarrhea_avertable / total_avertable_rate >= 0.5, 5, value)]

# mal 50-60% 6
dalys <- dalys[, value := ifelse(malaria_avertable / total_avertable_rate < 0.6 & malaria_avertable / total_avertable_rate >= 0.5, 6, value)]

# all under 50% 7
dalys <- dalys[, value := ifelse(lri_avertable / total_avertable_rate < 0.5 & diarrhea_avertable / total_avertable_rate < 0.5 & malaria_avertable / total_avertable_rate < 0.5, 7, value)]

# all 20-40% 8
dalys <- dalys[, value := ifelse(lri_avertable / total_avertable_rate <= 0.4 & diarrhea_avertable / total_avertable_rate <= 0.4 & malaria_avertable / total_avertable_rate <= 0.4, 8, value)]

# less than 10k averted daly: 99
dalys <- dalys[, value := ifelse(total_avertable_count < 10000, 99, value)]


# (5) Save to tridalys dir, format and save to mapping dir ##############################################################################
write.csv(dalys, paste0(in_dir, 'results/tridalys_daly_adm0_counterfactual.csv'))

map_dt <- select(dalys, ADM2_CODE, year, value)
write.csv(map_dt, paste0(map_dir, 'tridaly_daly_mean_averted_proportions_raked_ad2.csv'))
