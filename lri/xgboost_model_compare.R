cor1 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(ad1_unraked_cor_yes_xgboost = ad1_unraked_cor) %>%
  rename(ad2_unraked_cor_yes_xgboost = ad2_unraked_cor) %>%
  rename(ad2_raked_cor_yes_xgboost = ad2_raked_cor) %>%
  rename(ad1_raked_cor_yes_xgboost = ad1_raked_cor) %>%
  select(-V1)
  
cor2 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  rename(ad1_unraked_cor_no_xgboost = ad1_unraked_cor) %>%
  rename(ad2_unraked_cor_no_xgboost = ad2_unraked_cor) %>%
  rename(ad2_raked_cor_no_xgboost = ad2_raked_cor) %>%
  rename(ad1_raked_cor_no_xgboost = ad1_raked_cor) %>%
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
