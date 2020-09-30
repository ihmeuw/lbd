#####################################################################
# Produce csvs of raking factors for malaria
#####################################################################

# (1) Setup ---------------------------------------------------------
rm(list = ls())

run_date <- '2019_08_15'
shapefile_version <- '2019_02_27'

library(data.table)
library(dplyr)
library(ggplot2)

source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/util/location_metadata_functions.R')
source('<<< FILEPATH REDACTED >>>/lbd_core/mbg_central/shapefile_functions.R')


# (2) Read in raking factor csvs -------------------------------------

yld_rf <- fread('<<< FILEPATH REDACTED >>>') %>%
  dplyr::rename(mbg_yld = mbg_prev, gbd_yld = gbd_prev, rf_yld = rf) %>%
  dplyr::select(-V1)
yll_rf <- fread('<<< FILEPATH REDACTED >>>') %>%
  dplyr::rename(mbg_yll = mbg_prev, gbd_yll = gbd_prev, rf_yll = rf) %>%
  dplyr::select(-V1)

link <- readRDS('<<< FILEPATH REDACTED >>>')

rf <- merge(yld_rf, yll_rf, by = c('location_id','year'))  

code_map <- get_location_code_mapping(shapefile_version) %>%
  dplyr::select(loc_name, loc_id, ihme_lc_id)

reg_info <- fread('<<< FILEPATH REDACTED >>>') %>%
  dplyr::select(iso3, dia_reg)

rf <- merge(rf, code_map, by.x = 'location_id', by.y = 'loc_id') %>%
  merge(reg_info, by.x = 'ihme_lc_id', by.y = 'iso3')


# (3) Make plots  -------------------------------------
#plot rf by modeling region where color = country, x = year, y = rf
reg_list <- unique(rf$dia_reg)
dir.create('<<< FILEPATH REDACTED >>>')
pdf('<<< FILEPATH REDACTED >>>', width = 10, height = 8)

for (reg in reg_list){
  reg_data <- filter(rf, dia_reg == reg)
  
  inf_data <- filter(reg_data, rf_yld == 'Inf')
  reg_data <- filter(reg_data, rf_yld != 'Inf')
  inf_country <- paste(unique(inf_data$loc_name), sep="", collapse=", ")
  zero_data <- filter(reg_data, rf_yld == 0)
  zero_country <- paste(unique(zero_data$loc_name), sep="", collapse=", ")
  one_data <- filter(reg_data, rf_yld == 1)
  one_country <- paste(unique(one_data$loc_name), sep="", collapse=", ")
  if (inf_country == '') inf_country <- 'no countries'
  if (one_country == '') one_country <- 'no countries'
  if (zero_country == '') zero_country <- 'no countries'
  yld_plot <- ggplot(data = reg_data, aes(x = year, y = rf_yld)) + geom_line(aes(color = loc_name)) + 
    ggtitle(paste0('YLD raking factors for ', reg)) + ylab('YLD raking factor') + labs(caption = paste0('Raking factor = Inf for ', inf_country, ', Raking factor = 0 for ', zero_country, ', Raking factor = 1 for ', one_country))
  plot(yld_plot)
}

dev.off()

# (4) Save out csv -------------------------------------
write.csv(rf, '<<< FILEPATH REDACTED >>>')
