setwd("<<<< FILEPATH REDACTED >>>>")

library(dplyr)

reg_list <- readRDS('00_reg_list.rds')
reg_list[['ind_pak_chn']] <- c('IND', 'PAK', 'CHN')
reg_list[['ind_pak']] <- c('IND', 'PAK')
saveRDS(reg_list, '00_reg_list.rds')

water_covs <- readRDS('00_water_selected_covs.rds')
water_covs[['ind_pak']] <- water_covs[['dia_south_asia']]
water_covs[['ind_pak_chn']] <- water_covs[['dia_south_asia']]
saveRDS(water_covs, '00_water_selected_covs.rds')

sani_covs <- readRDS('00_sani_selected_covs.rds')
sani_covs[['ind_pak']] <- water_covs[['dia_south_asia']]
sani_covs[['ind_pak_chn']] <- water_covs[['dia_south_asia']]
saveRDS(sani_covs, '00_sani_selected_covs.rds')

gbm_params <- read.csv('00_gbm_params.csv', stringsAsFactors = FALSE)
new_gbm <- filter(gbm_params, region == 'dia_south_asia') %>%
	mutate(region = 'ind_pak')
new_gbm <- filter(gbm_params, region == 'dia_south_asia') %>%
	mutate(region = 'ind_pak_chn')
gbm_params <- bind_rows(gbm_params, new_gbm)
write.csv(gbm_params, '00_gbm_params.csv', row.names = FALSE)

comp_params <- read.csv('00_comp_params.csv', stringsAsFactors = FALSE)
comp_params <- bind_rows(comp_params, data.frame(region = 'ind_pak',
memory =  200, runtime = 24))
comp_params <- bind_rows(comp_params, data.frame(region = 'ind_pak_chn',
memory =  200, runtime = 24))
write.csv(comp_params, '00_comp_params.csv', row.names = FALSE)

