###########################################################################################
# Check to make sure that all nids in vetting tracker made it through collapse, and that no
# nid in collapse was not tracked in vetting tracker
###########################################################################################

#(1) setup --------------------------------------------------------------------------------
rm(list = ls())

library(dplyr)
library(data.table)
library(stringr)

select <- dplyr::select

#list reasons for dropping a survey that cannot be fixed
cant_fix <- c('No DB', 'Intentionally excluded', 'No ages under 5')

#(2) read in data -------------------------------------------------------------

# input and output
input <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  unique()
tab_input <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  unique()
output <- fread('<<<< FILEPATH REDACTED >>>>')

#pull excluded surveys and ISO3 list of stg2 countries
exclude_list <- fread('<<<< FILEPATH REDACTED >>>>')
stg2_list <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(Stage %in% c('1', '2a', '2b'))
stg2_iso3 <- stg2_list$iso3

#pull mat's spreadsheet of post-collapse problems
post_collapse <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  mutate(country = str_trunc(iso3, width = 3, ellipsis = '')) %>%
  filter(country %in% stg2_iso3)

#merge inputs, and exclude all data to stg2 countries
input_total <- rbind(input, tab_input) %>%
  mutate(country = str_trunc(Country, width = 3, ellipsis = '')) %>%
  filter(country %in% stg2_iso3)
output <- filter(output, country %in% stg2_iso3)
output_nids <- select(output, nid, survey_series, year, country) %>%
  unique.data.frame()

# (3) Find extra and dropped surveys -----------------------------------------------------

#find surveys in output that were not in tracking sheet
extra <- filter(output_nids, !(nid %in% input_total$NID))

#dropped
dropped <- filter(input_total, !(NID %in% output_nids$nid)) %>%
  filter(!(NID %in% exclude_list$nid))

#for dropped, merge on known collapse problems list
reason_cols <- select(post_collapse, nid, reason_dropped)

#filter dropped for nids we need to look into
dropped <- merge(dropped, reason_cols, by.x = 'NID', by.y = 'nid', all.x = TRUE)
actionable_dropped <- filter(dropped, !(reason_dropped %in% cant_fix)) %>%
  as.data.table()
actionable_dropped[order(reason_dropped),]

geo_matched <- fread('<<<< FILEPATH REDACTED >>>>')
filter(geo_matched, nid == 2063)

link <- readRDS('<<<< FILEPATH REDACTED >>>>')
cuba <- filter(link, ADM0_NAME == 'Cuba') %>%
  select(ADM1_NAME, ADM2_NAME) %>%
  unique.data.frame() %>%
  as.data.table()
mdg[order(ADM1_NAME),]
