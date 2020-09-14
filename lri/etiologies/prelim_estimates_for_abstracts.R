#######################################################################################################
# pull prelim estimates for Hib and P. pneumo abstracts
#######################################################################################################

# (1) Setup ------------------------------------------------------------------------------------------
rm(list = ls())

in_dir <- '<<<< FILEPATH REDACTED >>>>'

countries_to_exclude <- c(33, 38, 44, 52, 66, 69, 94, 128, 158, 239) #iso3 with no data currently

# (2) Pneumococcal pneumonia estimates ----------------------------------------------------------------------------
lri_0 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017)
lri_deaths <- sum(lri_0$mean)

p_pneumo_0 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  as.data.table()

p_pneumo_0_rate <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  as.data.table()

total_deaths <- sum(p_pneumo_0$mean)
total_deaths
lri_deaths

total_deaths/lri_deaths

top_country <- p_pneumo_0[order(-mean)]
head(top_country)

top_country <- p_pneumo_0_rate[order(-mean)]
head(top_country)

#inequality
ad0_lri <- p_pneumo_0_rate %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad0_mean = mean)

ad2_lri <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad2_mean = mean) %>%
  as.data.table()

ad2_lri <- ad2_lri[order(-ad2_mean)]

compare <- merge(ad0_lri, ad2_lri, by = c('ADM0_CODE', 'ADM0_NAME', 'year')) %>%
  filter(ad0_mean > 0.0005) %>%
  as.data.table()
compare[,ad2_over_ad0 := ad2_mean/ad0_mean]
compare <- compare[order(-ad2_over_ad0)]

ad2_lri_count <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad2_mean = mean) %>%
  as.data.table()

ad2_lri_count <- ad2_lri_count[order(-ad2_mean)]

#greatest declines
ad2_lri_2017 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region, -year, -pop) %>%
  dplyr::rename (mean_17 = mean) %>%
  as.data.table()

ad2_lri_2000 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2000) %>%
  dplyr::select(-lower, -upper, -cirange, -region, -year, -pop) %>%
  dplyr::rename (mean_00 = mean) %>%
  as.data.table()

data <- merge(ad2_lri_2000, ad2_lri_2017)
data <- data[, new_over_old := mean_17/mean_00] %>%
  as.data.table()
data <- data[order(new_over_old)]


# (3) Hib pneumonia estimates ----------------------------------------------------------------------------
lri_0 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017)
lri_deaths <- sum(lri_0$mean)

hib_0 <- fread(V) %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  as.data.table()

hib_0_rate <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  as.data.table()

total_deaths <- sum(hib_0$mean)
total_deaths
lri_deaths

total_deaths/lri_deaths

top_country <- hib_0_rate[order(-mean)]
head(top_country)

top_country <- hib_0[order(-mean)]
head(top_country)

#inequality
ad0_lri <- hib_0_rate %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad0_mean = mean)

ad2_lri <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad2_mean = mean)

compare <- merge(ad0_lri, ad2_lri, by = c('ADM0_CODE', 'ADM0_NAME', 'year')) %>%
  as.data.table()
compare[,ad2_over_ad0 := ad2_mean/ad0_mean]
compare <- compare[order(-ad2_over_ad0)]

ad2_lri_count <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad2_mean = mean) %>%
  as.data.table()

ad2_lri_count <- ad2_lri_count[order(-ad2_mean)]

#greatest declines
ad2_lri_2017 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region, -year, -pop) %>%
  dplyr::rename (mean_17 = mean) %>%
  as.data.table()

ad2_lri_2000 <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(!(ADM0_CODE %in% countries_to_exclude)) %>%
  filter(year == 2000) %>%
  dplyr::select(-lower, -upper, -cirange, -region, -year, -pop) %>%
  dplyr::rename (mean_00 = mean) %>%
  as.data.table()

data <- merge(ad2_lri_2000, ad2_lri_2017)
data <- data[, new_over_old := mean_17/mean_00] %>%
  as.data.table()
data <- data[order(new_over_old)]
