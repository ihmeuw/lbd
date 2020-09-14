#################################################################
# Compare two model run dates against each other at admin2 level
#################################################################

# (1) setup ------------------------------------------------------

rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)

#all stg 2 countries
all_stg2 <- fread('<<<< FILEPATH REDACTED >>>>')

#user inputs
run_date_1 <- '2020_01_10_15_18_27'
rd_1_name <- 'with xgboost'

run_date_2 <- '2020_03_09_11_08_00'
rd_2_name <- 'without xgboost'

admin_0_codes <- all_stg2$ADM0_CODE
  
output_pdf <- '<<<< FILEPATH REDACTED >>>>'

# (2) read in necessary data --------------------------------------

rd_1_prev_raked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, ADM0_CODE, year, mean) %>%
  rename(rd_1_mean = mean)

rd_2_prev_raked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, ADM0_CODE, year, mean) %>%
  rename(rd_2_mean = mean)

rd_1_prev_unraked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, ADM0_CODE, year, mean) %>%
  rename(rd_1_mean = mean)

rd_2_prev_unraked <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  select(ADM2_CODE, ADM0_CODE, year, mean) %>%
  rename(rd_2_mean = mean)

data_raked <- merge(rd_1_prev_raked, rd_2_prev_raked, by = c('ADM2_CODE','ADM0_CODE','year'))
data_raked_2000 <- filter(data_raked, year == 2000)
data_raked_2017 <- filter(data_raked, year == 2017)
data_raked_2000_2017 <- filter(data_raked, year %in% c(2000, 2017))

data_unraked <- merge(rd_1_prev_unraked, rd_2_prev_unraked, by = c('ADM2_CODE','ADM0_CODE','year'))
data_unraked_2000 <- filter(data_unraked, year == 2000)
data_unraked_2017 <- filter(data_unraked, year == 2017)
data_unraked_2000_2017 <- filter(data_unraked, year %in% c(2000, 2017))

# (3) plots -------------------------------------------------
pdf(output_pdf, 10, 5)

for (country in admin_0_codes){
  country_data_raked_2000_2017 <- filter(data_raked_2000_2017, ADM0_CODE == country)
  country_data_unraked_2000_2017 <- filter(data_unraked_2000_2017, ADM0_CODE == country)
  
  if (nrow(country_data_raked_2000_2017) == 0){
    print(paste0('No model results in ', all_stg2[ADM0_CODE == country,ADM0_NAME], '; skipping plot'))
  } else {
    plot2 <- ggplot(country_data_raked_2000_2017, aes(x = rd_1_mean, y = rd_2_mean)) + geom_point() + facet_wrap(~year) +
      ggtitle(paste0(all_stg2[ADM0_CODE == country,ADM0_NAME], ' (raked prevalence)')) + xlab(rd_1_name) + ylab(rd_2_name) +
      geom_abline(intercept = 0, slope = 1)
    plot(plot2)
  }
}

dev.off()