################################################################
# pulling various numbers for ASTMH presentation
################################################################

#input data information
inputs <- fread('<<<< FILEPATH REDACTED >>>>')

#surveys
length(unique(inputs$nid))

#countries
length(unique(inputs$country))
unique(inputs$country) #subtract out kenya nid labels
length(unique(inputs$country)) - 3

#DHS
dhs <- filter(inputs, survey_series == 'MACRO_DHS')
length(unique(dhs$nid))

#MICS
mics <- filter(inputs, survey_series == 'UNICEF_MICS')
length(unique(mics$nid))

#points
points <- filter(inputs, point == 1)
nrow(points)

#polygons
poly <- filter(inputs, point == 0)
nrow(poly)

#inequality
ad0_lri <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad0_mean = mean)

ad2_lri <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(year == 2017) %>%
  dplyr::select(-lower, -upper, -cirange, -region) %>%
  dplyr::rename (ad2_mean = mean)

compare <- merge(ad0_lri, ad2_lri, by = c('ADM0_CODE', 'ADM0_NAME', 'year')) %>%
  as.data.table()
compare[,ad2_over_ad0 := ad2_mean/ad0_mean]
compare <- compare[order(-ad2_over_ad0)]

peru <- filter(compare, ADM0_NAME == 'Peru')
myanmar <- filter(compare, ADM0_NAME == 'Myanmar')
