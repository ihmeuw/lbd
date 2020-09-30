missing_ad2 <- lri_to_find$ADM2_CODE

#all are present in lri file
lri_dat <- fread('<<< FILEPATH REDACTED >>>')
lri_present <- filter(lri_dat, ADM2_CODE %in% missing_ad2)
length(missing_ad2)
length(unique(lri_present$ADM2_CODE))

#in cgf file
cgf_dat <- fread('<<< FILEPATH REDACTED >>>')
cgf_present <- filter(cgf_dat, ADM2_CODE %in% missing_ad2)
length(missing_ad2)
length(unique(cgf_present$ADM2_CODE))
cgf_present

cgf_absent <- filter(cgf_dat, !(ADM2_CODE %in% missing_ad2) & ADM2_CODE %in% ad2_list$ADM2_CODE)
filter(lri_dat, ADM2_CODE ==   10045065 )
filter(lri_results, ADM2_CODE ==   10045065 )
filter(cgf_dat, ADM2_CODE ==   10045065 )

under_dat <- fread('<<< FILEPATH REDACTED >>>')
stunt_dat <- fread('<<< FILEPATH REDACTED >>>')
wast_dat <- fread('<<< FILEPATH REDACTED >>>')

under_present <- filter(under_dat, ADM2_CODE %in% missing_ad2) %>%
  as.data.table()
under_present <- under_present[,under_value := ifelse(is.na(mean), 0, 1)] %>%
  select(-region, -upper, -lower, -cirange)

stunt_present <- filter(stunt_dat, ADM2_CODE %in% missing_ad2)  %>%
  as.data.table()
stunt_present <- stunt_present[,stunt_value := ifelse(is.na(mean), 0, 1)] %>%
  select(-region, -upper, -lower, -cirange)

wast_present <- filter(wast_dat, ADM2_CODE %in% missing_ad2)  %>%
  as.data.table()
wast_present <- wast_present[,wast_value := ifelse(is.na(mean), 0, 1)]  %>%
  select(-region, -upper, -lower, -cirange)

nrow(wast_present)
nrow(stunt_present)
nrow(under_present)

compare <- merge(under_present, stunt_present, by = c('ADM0_CODE','ADM0_NAME','ADM1_CODE','ADM1_NAME','ADM2_CODE','ADM2_NAME','year')) %>%
  merge(wast_present, by = c('ADM0_CODE','ADM0_NAME','ADM1_CODE','ADM1_NAME','ADM2_CODE','ADM2_NAME','year')) %>%
  select(-mean.y, -mean.x, -mean)

different <- filter(compare, under_value != stunt_value | under_value != wast_value | wast_value != stunt_value) %>%
  select(-year) %>%
  unique.data.frame()

one_or_more_na <- filter(compare, under_value == 0 | wast_value == 0 | stunt_value == 0 ) %>%
  select(-year) %>%
  unique.data.frame()

#find not present admins
any_cgf_val <- as.data.table(lri_to_find)
any_cgf_val[,in_under := ifelse(ADM2_CODE %in% under_dat$ADM2_CODE, 1, 0)]
any_cgf_val[,in_wast := ifelse(ADM2_CODE %in% wast_dat$ADM2_CODE, 1, 0)]
any_cgf_val[,in_stunt := ifelse(ADM2_CODE %in% stunt_dat$ADM2_CODE, 1, 0)]

not_present <- filter(any_cgf_val, in_under == 0 & in_wast == 0 & in_stunt ==0)


#########################################################
#do we have estimates everywhere pop != 0 for cgf
#########################################################
pop_file <- '<<< FILEPATH REDACTED >>>'

under_dat <- fread('<<< FILEPATH REDACTED >>>')
stunt_dat <- fread('<<< FILEPATH REDACTED >>>')
wast_dat <- fread('<<< FILEPATH REDACTED >>>')

load(pop_file)
pops <- admin_2 %>%
  filter(year == 2017) %>%
  dplyr::select(pop, ADM2_CODE) %>%
  rename(pop_2017 = pop)

s_w_ad2 <- unique(wast_dat$ADM2_CODE)

sw_missing <- filter(lri_dat, !(ADM2_CODE %in% s_w_ad2)) %>%
  select(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME) %>%
  unique.data.frame()

under_dat <- merge(under_dat, pops, by = 'ADM2_CODE') %>%
  filter(pop_2017 > 0)
stunt_dat <- merge(stunt_dat, pops, by = 'ADM2_CODE')  %>%
  filter(pop_2017 > 0)
wast_dat <- merge(wast_dat, pops, by = 'ADM2_CODE') %>%
  filter(pop_2017 > 0)

##############
#check against shapefile
shapefile_version <- '2019_02_27'
link <- readRDS('<<< FILEPATH REDACTED >>>') %>%
  dplyr::select(ADM2_CODE, ADM2_NAME, ADM1_CODE, ADM1_NAME, ADM0_CODE, ADM0_NAME) %>%
  unique.data.frame()
stage_list <- fread('<<< FILEPATH REDACTED >>>') %>%
  filter(Stage == 1)

# pull stage master list
stg1_ad0_codes <- stage_list$gadm_geoid
ad2_list <- filter(link, ADM0_CODE %in% stg1_ad0_codes)

missing_from_ad2_list <- filter(ad2_list, !(ADM2_CODE %in% s_w_ad2)) %>%
  select(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME) %>%
  unique.data.frame()

extra <- filter(under_dat, !(ADM2_CODE %in% ad2_list$ADM2_CODE)) %>%
  select(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME) %>%
  unique.data.frame()

head(extra)
head(missing_from_ad2_list)

table(extra$ADM0_NAME)
