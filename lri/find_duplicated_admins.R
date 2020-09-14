library(dplyr)
library(data.table)
load('<<<< FILEPATH REDACTED >>>>')

length(unique(admin_0$ADM0_CODE))
nrow(admin_0) / 18

length(unique(admin_1$ADM1_CODE))
nrow(admin_1) / 18

length(unique(admin_2$ADM2_CODE))
nrow(admin_2) / 18

admin_0_counts <- admin_0 %>%
  dplyr::count(ADM0_CODE) %>%
  filter(n > 18)

admin_1_counts <- admin_1 %>%
  dplyr::count(ADM1_CODE) %>%
  filter(n > 18)

admin_2_counts <- admin_2 %>%
  dplyr::count(ADM2_CODE) %>%
  filter(n > 18)

ad0_sum <- fread('<<<< FILEPATH REDACTED >>>>')
ad1_sum <- fread('<<<< FILEPATH REDACTED >>>>')
ad2_sum <- fread('<<<< FILEPATH REDACTED >>>>')

filter(ad0_sum, ADM0_CODE %in% admin_0_counts$ADM0_CODE)
filter(ad1_sum, ADM1_CODE %in% admin_1_counts$ADM1_CODE)
filter(ad2_sum, ADM1_CODE %in% admin_1_counts$ADM1_CODE)

load('<<<< FILEPATH REDACTED >>>>')
filter(admin_1, ADM1_CODE == 15039)
