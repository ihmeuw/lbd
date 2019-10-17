## #######################################################
## PURPOSE: generate variables for sbh checks
## INPUTS: sbh data.table
## #######################################################
library(data.table)
library(Biograph)
library(dplyr)

dat_dir <- '<<<< FILEPATH REDACTED >>>>'
d <- setDT(readRDS(sprintf('<<<< FILEPATH REDACTED >>>>',dat_dir)))

#split nid 91740 into separate checks
d$nid <- as.character(d$nid)
d[nid == 91740, nid := paste0(nid,"_",country)]

## #######################################################
## CLEAN DATA AND RUN DATA QUALITY CHECKS SCRIPT

# check for missing geographic identifiers, make geo unit identifier for later on collapse
d[, location_code := as.character(location_code)]
d$location_code[d$location_code=='NA']=''
d[, latnum        := as.numeric(latnum)]
d[, longnum       := as.numeric(longnum)]
d$point[is.na(d$point) & (d$location_code != "" & d$shapefile != "")] = 0
d$point[is.na(d$point) & (!is.na(d$latnum) & !is.na(d$longnum))] = 1
d[point==0, geo_unit := paste0(nid,location_code,shapefile)]
d[point==1, geo_unit := paste0(nid,cluster_number)]

d[, no_geog := (is.na(latnum) | is.na(longnum)) & (shapefile == '' | location_code == '')]

# set maternal age variable names
setnames(d, 'age_group_of_woman', 'mother_age_bin')
setnames(d, 'age_year', 'mothage_atint')

# Point data will not be weighted, so we convert at point row weights to 1
d$weight[d$point == 1] <- 1

# check for missing interview date
d[, yrint := cmc_as_year(interview_date_cmc)]

# start a variable indicating rows to drop, not including weights in this
d[, must_drop := is.na(ceb) | is.na(ced) | ced > ceb | no_geog == TRUE | is.na(mothage_atint)]

# mothage needs to be int
d[, mother_age_bin := round(mother_age_bin, 0)]
d[, mothage_atint := round(mothage_atint, 0)]

# get number of mothers who reported no children ever born by survey
summary_ceb_is_0 <- d %>% 
  group_by(nid, source, country, year) %>% 
  summarise(n = n(),
            ceb_is_0 = sum(ceb==0, na.rm=T))

# get average number of children and number + percent of mothers who reported no children ever born by age bin
summary_ceb_age_bin <- d %>% 
  group_by(nid, mother_age_bin) %>% 
  summarise(total_ceb = n(),
            mean_child = mean(ceb, na.rm=T),
            ceb_0 = sum(ceb == 0, na.rm=T)) %>%
  mutate(ceb_0_percent = (ceb_0 / total_ceb) * 100)
summary_ceb_age_bin <- subset(summary_ceb_age_bin, select=-c(total_ceb))

# drop rows where mother reports no children ever born
d <- subset(d, ceb > 0)

# generate indicators used in warning checks
summary_predrop <- d %>% 
  group_by(nid, source, country, year) %>% 
  summarise(n = n(),
            point_true = sum(point == 1, na.rm=T),
            point_false = sum(point == 0, na.rm=T),
            no_geog = sum(no_geog, na.rm=T),
            point_na = sum(is.na(point)),
            missing_ceb = sum(is.na(ceb)),
            missing_ced = sum(is.na(ced)),
            ced_over_ceb = sum(ced > ced, na.rm=T),
            interview_date_cmc_missing = sum(is.na(interview_date_cmc)),
            mothers_age_missing = sum(is.na(d$mothage_atint)),
            mothers_age_bin_missing = sum(is.na(d$mother_age_bin)),
            weights_missing = sum(is.na(weight)),
            must_drop = sum(must_drop, na.rm=T),
            ceb_is_0 = sum(ceb==0, na.rm=T))

# dropped surveys missing key variables
dr <- subset(d, must_drop == FALSE)

# keep only variables of interest
dr <- dr[, c('nid','source','country','year','strata','weight',
           'point','latnum','longnum','location_code','shapefile',
           'geo_unit','yrint','ceb','ced','mother_age_bin'),
       with = FALSE]

# get final ceb, ced, and raw ced/ceb by age bin for checks table
summary_complete <- dr %>% 
  group_by(nid, mother_age_bin) %>% 
  summarise(n = n(),
            ceb = sum(ceb, na.rm=T),
            ced = sum(ced, na.rm=T)) %>%
  mutate(raw_ced_ceb = (ced / ceb))

## Create SBH checks report - must run on singularity image
rmarkdown::render("<<<< FILEPATH REDACTED >>>>", params = list(
  dt = d,
  summary_predrop = summary_predrop,
  summary_complete = summary_complete,
  summary_ceb_is_0 = summary_ceb_is_0,
  summary_ceb_age_bin = summary_ceb_age_bin,
  outlier_path = "<<<< FILEPATH REDACTED >>>>"))

print("done")
