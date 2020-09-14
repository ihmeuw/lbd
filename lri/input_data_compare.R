#####################################################################
# Make scatterplot of two LRI input data sets 
#####################################################################

# (1) Setup ----------------------------------------------------------

#arguments
old_tag <- commandArgs()[4] #name in input data folder. if custom_old_data = T, specify any name for the "old" data to propogate to the file name
new_tag <- commandArgs()[5] #name in input data folder
custom_old_data <- as.logical(commandArgs()[6])
old_data_filepath <- commandArgs()[7]

indicator <- 'has_lri'
indicator_group <- 'lri'
run_date <- '2019_07_23_13_11_22'
regions <- c('cssa','essa','wssa','name','sssa')
modeling_shapefile_version <- '2019_02_27'

#libraries and functions
library(data.table)
library(dplyr)
library(ggplot2)
source('<<<< FILEPATH REDACTED >>>>/lri/post_estimation/functions/aggregate_data.R')
source('<<<< FILEPATH REDACTED >>>>/mbg_central/prep_functions.R')
source('<<<< FILEPATH REDACTED >>>>/mbg_central/shapefile_functions.R')


#link
link <- readRDS('<<<< FILEPATH REDACTED >>>>') %>%
  dplyr::select(c('ADM0_CODE', 'ADM0_NAME')) %>%
  unique()

#reg lookup
reg_table <- list()
for (reg in regions){
  codes <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  table <- data.frame('ADM0_CODE' = codes, 'region' = reg)
  reg_table <- rbind(reg_table, table)
}

# (2) Load in data ----------------------------------------------------
if (custom_old_data){
  old_data <- fread(old_data_filepath)
} else {
  old_data <- fread('<<<< FILEPATH REDACTED >>>>')
}
new_data <- fread('<<<< FILEPATH REDACTED >>>>')

#take out from new data any nids that were added since old collapse date
old_nids <- unique(old_data$nid)
new_data <- filter(new_data, nid %in% old_nids)

#take out from old data any nids that were outliered since old collapse date
new_nids <- unique(new_data$nid)
old_data <- filter(old_data, nid %in% new_nids)

if ('V1' %in% names(old_data)) old_data <- dplyr::select(old_data, -V1)
if ('V1' %in% names(new_data)) new_data <- dplyr::select(new_data, -V1)

write.csv(new_data, '<<<< FILEPATH REDACTED >>>>')
write.csv(old_data, '<<<< FILEPATH REDACTED >>>>')

old_agg <- aggregate_input_data(indicator,
                                indicator_group,
                                run_date,
                                regions,
                                modeling_shapefile_version,
                                file_name = paste0('has_lri_compare_plot_new_', new_tag))

new_agg <- aggregate_input_data(indicator,
                                indicator_group,
                                run_date,
                                regions,
                                modeling_shapefile_version,
                                file_name = paste0('has_lri_compare_plot_old_', old_tag))

old_ad0 <- old_agg$ad0 %>%
  dplyr::rename(old_mean = input_mean, old_ss = input_ss)

old_ad1 <- old_agg$ad1 %>%
  dplyr::rename(old_mean = input_mean, old_ss = input_ss)

new_ad0 <- new_agg$ad0 %>%
  dplyr::rename(new_mean = input_mean, new_ss = input_ss)

new_ad1 <- new_agg$ad1 %>%
  dplyr::rename(new_mean = input_mean, new_ss = input_ss)

#merge ad1, ad0
ad0 <- merge(old_ad0, new_ad0, by = c('nid','year','ADM0_CODE')) %>%
  merge(link, by = 'ADM0_CODE') %>%
  merge(reg_table, by = 'ADM0_CODE', all.x = T)
ad1 <- merge(old_ad1, new_ad1, by = c('nid','year','ADM0_CODE','ADM1_CODE')) %>%
  merge(link, by = 'ADM0_CODE') %>%
  merge(reg_table, by = 'ADM0_CODE', all.x = T)

# (3) Plots ----------------------------------------------------

pdf('<<<< FILEPATH REDACTED >>>>', width = 10, height = 8)
for (reg in regions){
  ad0_reg <- filter(ad0, region == reg)
  ad1_reg <- filter(ad1, region == reg)
  plot0 <- ggplot(ad0_reg, aes(x = old_mean, y = new_mean)) + geom_point(aes(color = ADM0_NAME, size = ((new_ss + old_ss)/2))) + geom_abline(slope = 1, intercept = 0) +
    ggtitle(paste0(reg, ' ADM0 aggregated prevalence')) + xlab(old_tag) + ylab(new_tag)
  plot1 <- ggplot(ad1_reg, aes(x = old_mean, y = new_mean)) + geom_point(aes(color = ADM0_NAME, size = ((new_ss + old_ss)/2))) + geom_abline(slope = 1, intercept = 0) + 
    ggtitle(paste0(reg, ' ADM1 aggregated prevalence')) + xlab(old_tag) + ylab(new_tag)
  plot(plot0)
  plot(plot1)
}
dev.off()                                                       

pdf('<<<< FILEPATH REDACTED >>>>', width = 10, height = 8)
for (reg in regions){
  ad0_reg <- filter(ad0, region == reg)
  ad1_reg <- filter(ad1, region == reg)
  plot0 <- ggplot(ad0_reg, aes(x = old_ss, y = new_ss)) + geom_point(aes(color = ADM0_NAME)) + geom_abline(slope = 1, intercept = 0) +
    ggtitle(paste0(reg, ' ADM0 aggregated sample sizes')) + xlab(old_tag) + ylab(new_tag)
  plot1 <- ggplot(ad1_reg, aes(x = old_ss, y = new_ss)) + geom_point(aes(color = ADM0_NAME)) + geom_abline(slope = 1, intercept = 0) + 
    ggtitle(paste0(reg, ' ADM1 aggregated sample sizes')) + xlab(old_tag) + ylab(new_tag)
  plot(plot0)
  plot(plot1)
}
dev.off() 