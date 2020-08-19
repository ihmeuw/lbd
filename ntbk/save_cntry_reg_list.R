# sSet up environment with packages and working directory
library(dplyr)
setwd("<<<< FILEPATH REDACTED >>>>")

# Read region-country dictionary
reg_list <- readRDS('00_reg_list.rds')

# construct data.frame from list
results <- list()
for (i in unique(names(reg_list))) {
  countries <- reg_list[[i]]
  region <- rep(i, length(countries))
  results[[1 + length(results)]] <- data.frame(cbind(countries, region))
}
results <- do.call(rbind, results)
results$countries <- as.character(results$countries)
results$region <- as.character(results$region)

# subset results to regions used
used_reg <- c('dia_afr_horn', 'dia_central_asia', 'dia_chn_mng', 'dia_cssa-cod',
  'cod', 'dia_essa-ken', 'ken', 'dia_malay-idn-png', 'idn', 'png',
  'dia_mcaca-mex-cri', 'mex', 'cri', 'dia_mid_east', 'dia_name',
  'dia_s_america-per', 'per', 'dia_se_asia', 'dia_south_asia', 'dia_sssa',
  'dia_wssa-civ', 'civ')
results <- filter(results, region %in% used_reg)

# output results
setwd("<<<< FILEPATH REDACTED >>>>")
write.csv(results, 'cntry_regs.csv', row.names = FALSE)
