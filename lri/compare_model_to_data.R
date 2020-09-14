##############################################################
# Compare aggregated model results to aggregated data (admin2)
#############################################################

# (1) Setup -------------------------------------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(data.table)
library(ggplot2)

rename <- dplyr::rename
select <- dplyr::select

run_date <- '2020_01_10_15_18_27'
modeling_shapefile_version <- '2019_09_10'

source('<<<< FILEPATH REDACTED >>>>/aggregate_data_ad2.R')
source('<<<< FILEPATH REDACTED >>>>/aggregate_data.R')
source('<<<< FILEPATH REDACTED >>>>/mbg_central/misc_functions.R')

dir.create('<<<< FILEPATH REDACTED >>>>')

# (2) Pull aggregated model results -------------------------------------------------------------------------

ad1_model_unraked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(model_mean = mean)

ad1_model_raked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(model_mean = mean)

ad2_model_unraked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(model_mean = mean)

ad2_model_raked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(model_mean = mean)


# (3) Aggregate input data -----------------------------------------------------------------------------------
regions <- unique(ad1_model_unraked$region)

ad1_data <- aggregate_input_data(indicator = 'has_lri',
                                 indicator_group = 'has_lri',
                                 run_date,
                                 regions,
                                 modeling_shapefile_version,
                                 file_name = 'has_lri')

ad1_data <- ad1_data$ad1

ad2_data <- aggregate_input_data_ad2(indicator = 'has_lri',
                                     indicator_group = 'has_lri',
                                     run_date,
                                     regions,
                                     modeling_shapefile_version,
                                     file_name = 'has_lri')

ad2_data <- ad2_data$ad2

# remove data not from correct country
orig_input <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(nid, country) %>%
  unique.data.frame() %>%
  rename(survey_country = country)

ad2_data <- merge(ad2_data, orig_input, by = 'nid')
ad1_data <- merge(ad1_data, orig_input, by = 'nid')

# merges
plot_data_ad1_unraked <- merge(ad1_model_unraked, ad1_data, by = c('year', 'ADM1_CODE')) %>%
  filter(survey_country == iso3)

plot_data_ad1_raked <- merge(ad1_model_raked, ad1_data, by = c('year', 'ADM1_CODE')) %>%
  filter(survey_country == iso3)

plot_data_ad2_unraked <- merge(ad2_model_unraked, ad2_data, by = c('year', 'ADM2_CODE')) %>%
  filter(survey_country == iso3)

plot_data_ad2_raked <- merge(ad2_model_raked, ad2_data, by = c('year', 'ADM2_CODE')) %>%
  filter(survey_country == iso3)

# (4) Plot and make correlation table ---------------------------------------------------------------------------------------------------

# all regions
for (type in c('ad1_unraked', 'ad1_raked', 'ad2_unraked', 'ad2_raked')){
  pdf('<<<< FILEPATH REDACTED >>>>', 8, 8)
  
  plot_data <- get(paste0('plot_data_', type))
  
  cor_all <- cor(plot_data$model_mean, plot_data$input_mean, use = 'complete.obs')
  
  plot_all_reg <- ggplot(plot_data, aes(x = input_mean, y = model_mean)) + geom_point() + xlab('aggregated input data') + ylab('aggregated raked prevalence') +
    ggtitle(paste0(type, ' model vs. data for run_date = ', run_date, ' (all regions)')) + labs(subtitle = paste0('R = ', cor_all)) + theme_classic() +
    geom_abline(slope = 1, intercept = 0)
  
  plot(plot_all_reg)
  
  # by region
  results_list <- data_frame('reg' = 'all', 'cor' = cor_all)
  
  for (a_region in regions){
    region_data <- filter(plot_data, region == a_region)
    cor <- cor(region_data$model_mean, region_data$input_mean, use = 'complete.obs')
    results <- data_frame('reg' = a_region,'cor' = cor)
    results_list <- rbind(results_list, results)
    
    plot_reg <- ggplot(region_data, aes(x = input_mean, y = model_mean)) + geom_point() + xlab('aggregated input data') + ylab('aggregated raked prevalence') +
      ggtitle(paste0(type, ' model vs. data for run_date = ', run_date, ' (', a_region, ')')) + labs(subtitle = paste0('R = ', cor)) + theme_classic() + 
      geom_abline(slope = 1, intercept = 0)
    
    plot(plot_reg)
  }
  
  dev.off()

  results_list <- dplyr::rename(results_list, !!paste0(type, '_cor') := cor)
  assign(paste0(type, '_cor_table'), results_list)
}

# merge cor tables
results_table <- merge(ad1_unraked_cor_table, ad1_raked_cor_table, by = 'reg') %>%
  merge(ad2_unraked_cor_table, by = 'reg') %>%
  merge(ad2_raked_cor_table, by = 'reg')

write.csv(results_table, '<<<< FILEPATH REDACTED >>>>')


cor1 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(ad1_unraked_cor_no_xgboost = ad1_unraked_cor) %>%
  rename(ad2_unraked_cor_no_xgboost = ad2_unraked_cor) %>%
  rename(ad2_raked_cor_no_xgboost = ad2_raked_cor) %>%
  rename(ad1_raked_cor_no_xgboost = ad1_raked_cor) %>%
  select(-V1)

cor2 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(ad1_unraked_cor_yes_xgboost = ad1_unraked_cor) %>%
  rename(ad2_unraked_cor_yes_xgboost = ad2_unraked_cor) %>%
  rename(ad2_raked_cor_yes_xgboost = ad2_raked_cor) %>%
  rename(ad1_raked_cor_yes_xgboost = ad1_raked_cor) %>%
  select(-V1)

data <- merge(cor1, cor2, by = 'reg') %>%
  as.data.table()

data[, ad1_raked_results := ifelse(ad1_raked_cor_no_xgboost > ad1_raked_cor_yes_xgboost, 'drop xgboost', 'keep xgboost')]
data[, ad1_unraked_results := ifelse(ad1_unraked_cor_no_xgboost > ad1_unraked_cor_yes_xgboost, 'drop xgboost', 'keep xgboost')]

data[, ad2_raked_results := ifelse(ad2_raked_cor_no_xgboost > ad2_raked_cor_yes_xgboost, 'drop xgboost', 'keep xgboost')]
data[, ad2_unraked_results := ifelse(ad2_unraked_cor_no_xgboost > ad2_unraked_cor_yes_xgboost, 'drop xgboost', 'keep xgboost')]

data_ad1_unraked <- select(data, reg, ad1_unraked_cor_yes_xgboost, ad1_unraked_cor_no_xgboost, ad1_unraked_results)
data_ad1_raked <- select(data, reg, ad1_raked_cor_yes_xgboost, ad1_raked_cor_no_xgboost, ad1_raked_results)
data_ad2_unraked <- select(data, reg, ad2_unraked_cor_yes_xgboost, ad2_unraked_cor_no_xgboost, ad2_unraked_results)
data_ad2_raked <- select(data, reg, ad2_raked_cor_yes_xgboost, ad2_raked_cor_no_xgboost, ad2_raked_results)


write.csv(data_ad1_unraked, '~/ad1_unraked.csv')
write.csv(data_ad1_raked, '~/ad1_raked.csv')
write.csv(data_ad2_unraked, '~/ad2_unraked.csv')
write.csv(data_ad2_raked, '~/ad2_raked.csv')

