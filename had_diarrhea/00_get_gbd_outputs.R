# Obtain GBD Values of Diarrhea Estimation

# Clear environment and import code to read GBD Databases
rm(list = ls())
source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
library('dplyr')
library('stats')
library('mgcv')

# Read in and subset identifiers to be passed to functions
## REI ID table of GBD Etiologies
eti <- read.csv('<<<< FILEPATH REDACTED >>>>/eti_rr_me_ids.csv')
## GBD Location ID Table
locs <- read.csv('<<<< FILEPATH REDACTED >>>>/ihme_loc_metadata_2017.csv')
## Subset locations to National and Subnational estimates
country <- subset(locs, level >=3)$location_id
## Subset etiology IDs to those relevant for Diarrhea
eti_ids <- subset(eti, cause_id==302)$rei_id

# Get GBD Estimates

## Obtain daly estimates by country and time period
daly <- get_outputs(topic= 'cause', cause_id=302, location_id=country, 
                    year_id=2000:2017,
                    age_group_id=1, gbd_round_id=5, metric_id=3,
                    measure_id = 2,
                    version='best') %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017)))))

daly <- daly[complete.cases(daly),]

daly$val <- ifelse(is.na(daly$val), daly$inter_val, daly$val)

## Obtain yld estimates by country and time period
yld <- get_outputs(topic= 'cause', cause_id=302, location_id=country, 
                   year_id=2000:2017,
                   age_group_id=1, gbd_round_id=5, metric_id=3,
                   measure_id = 3,
                   version='best') %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017)))))

yld <- yld[complete.cases(yld),]

yld$val <- ifelse(is.na(yld$val), yld$inter_val, yld$val)

## Obtain yll estimates by country and time period
yll <- get_outputs(topic= 'cause', cause_id=302, location_id=country, 
                   year_id=2000:2017,
                   age_group_id=1, gbd_round_id=5, metric_id=3,
                   measure_id = 4,
                   version='best') %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017)))))

yll <- yll[complete.cases(yll),]

yll$val <- ifelse(is.na(yll$val), yll$inter_val, yll$val)

## Obtain prevalence estimates by country and time period
prevalence <- get_outputs(topic= 'cause', cause_id=302, location_id=country, 
                          year_id=2000:2017,
                          age_group_id=1, gbd_round_id=5, metric_id=3,
                          measure_id = 5,
                          version='best') %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017)))))

prevalence <- prevalence[complete.cases(prevalence),]

prevalence$val <- ifelse(is.na(prevalence$val), prevalence$inter_val, prevalence$val)

## Obtain deaths estimates by country and time period
deaths <- get_outputs(topic= 'cause', cause_id=302, location_id=country, 
                      year_id=2000:2017,
                      age_group_id=1, gbd_round_id=5, metric_id=3,
                      measure_id = 1,
                      version='best') %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017)))))

deaths <- deaths[complete.cases(deaths),]

deaths$val <- ifelse(is.na(deaths$val), deaths$inter_val, deaths$val)

## Obtain incidence estimates by country and time period
incidence <- get_outputs(topic= 'cause', cause_id=302, location_id=country, 
                         year_id=2000:2017,
                         age_group_id=1, gbd_round_id=5, metric_id=3,
                         measure_id = 6,
                         version='best') %>%
  group_by(location_id) %>%
  mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)),
                                            newdata = data.frame(year_id = 2000:2017)))))

incidence <- incidence[complete.cases(incidence),]

incidence$val <- ifelse(is.na(incidence$val), incidence$inter_val, incidence$val)

## Obtain incidence PAF by etiology, time period, and location
eti_nonfatal <- get_outputs(topic='rei', rei_id=eti_ids, cause_id=302, 
                            metric_id=2,
                            gbd_round_id=5, version='best', age_group_id=1,
                            location_id=country,
                            measure_id=3,
                            year_id = 2000:2017) %>%
  group_by(location_id, rei_id) %>%
  mutate(inter_val = predict(glm(val ~ year_id),
                             newdata = data.frame(year_id = 2000:2017)))

eti_nonfatal <- eti_nonfatal[complete.cases(eti_nonfatal),]

eti_nonfatal$val <- ifelse(is.na(eti_nonfatal$val), eti_nonfatal$inter_val, eti_nonfatal$val)

## Obtain mortality PAF by etiology, time period, and location
eti_fatal <- get_outputs(topic='rei', rei_id=eti_ids, cause_id=302, 
                         metric_id=2,
                         gbd_round_id=5, version='best', age_group_id=1,
                         location_id=country,
                         measure_id=4,
                         year_id = 2000:2017) %>%
  group_by(location_id, rei_id) %>%
  mutate(inter_val = predict(glm(val ~ year_id),
                             newdata = data.frame(year_id = 2000:2017)))

eti_fatal <- eti_fatal[complete.cases(eti_fatal),]

eti_fatal$val <- ifelse(is.na(eti_fatal$val), eti_fatal$inter_val, eti_fatal$val)

## Multiply mortality by etiology PAFs to get etiology mortality
eti_fatal <- merge(data.table(eti_fatal), 
                   data.table(deaths)[, c('location_id', 'year_id', 'val')],
                   by = c('location_id', 'year_id'), allow.cartesian = T)
eti_fatal[, val := val.x*val.y]
eti_fatal[, c('val.x', 'val.y') := NULL]

## Multiply incidence by etiology PAFs to get etiology incidence
eti_nonfatal <- merge(data.table(eti_nonfatal), 
                      data.table(incidence)[, c('location_id', 'year_id', 'val')],
                      by = c('location_id', 'year_id'), allow.cartesian = T)
eti_nonfatal[, val := val.x*val.y]
eti_nonfatal[, c('val.x', 'val.y') := NULL]

# Write CSVs of every GBD Output
metrics <- c('eti_nonfatal', 'eti_fatal', 'deaths', 'daly', 'yld', 'yll',
             'prevalence', 'incidence')

for (m in metrics) {
  write.csv(get(m), paste0('<<<< FILEPATH REDACTED >>>>/gbd_', m, '.csv'), row.names = FALSE)
}