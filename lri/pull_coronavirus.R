#################################################################
# Pull LRI and S.Pneumo mortality counts in China
# for coronavirus work
##################################################################


# (1) Setup ---------------------------------------------------------------------------------------------------------------------------------

rm(list = ls())
source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
source('<<<< FILEPATH REDACTED >>>>/get_population.R')

# (2) LRI mortality and inc count in China, 2000-2017, by most granular age and sex -------------------------------------------------------------
lri_count <- get_outputs(topic= 'cause',
                         cause_id=322,
                         location_id=6,
                         year_id=c(2000:2017),
                         age_group_id= c(2:22, 28, 30:32, 42, 235),
                         gbd_round_id=5,
                         metric_id= 1,
                         measure_id = c(1,6),
                         sex_id = c(1:2),
                         version = 'latest')

#edit out unnecessary columns and add ID column "lri"
lri_count$cause <- 'lri'
lri_count <- select(lri_count, age_group_name, year_id, location_name, measure_name, metric_name, sex, val, upper, lower, cause, 'age_group_id', 'location_id', 'year_id', 'sex_id')

# (3) S. pneumo mortality count in China, 2000-2017, by most granular age and sex -------------------------------------------------------------

pneumo_count <- get_outputs(topic= 'rei',
                            rei_id=188,
                            cause_id = 322,
                            location_id=6,
                            year_id=c(2000:2017),
                            age_group_id= c(2:22, 28, 30:32, 42, 235),
                            gbd_round_id=5,
                            metric_id= 1,
                            measure_id = 1,
                            sex_id = c(1:2), 
                            version = 'latest')

#edit out unnecessary columns and add ID column "lri"
pneumo_count$cause <- 'lri_pneumo'
pneumo_count <- select(pneumo_count, age_group_name, year_id, location_name, measure_name, metric_name, sex, val, upper, lower, cause, 'age_group_id', 'location_id', 'year_id', 'sex_id')

#bind results
results <- rbind(lri_count, pneumo_count)


# (4) Get population in China, 2000-2017, by most granular age and sex -------------------------------------------------------------
china_pop <- get_population(location_id=6,
                            year_id=c(2000:2017),
                            age_group_id= c(2:22, 28, 30:32, 42, 235),
                            gbd_round_id=5,
                            sex_id = c(1:2))

gbd_lri_chn <- merge(results, china_pop, by = c('age_group_id', 'location_id', 'year_id', 'sex_id'), all.x = TRUE) %>%
  select(-age_group_id, -location_id, -sex_id, -run_id)

# (5) Save out estimates -------------------------------------------------------------------------------------------------------------
saveRDS(gbd_lri_chn, '<<<< FILEPATH REDACTED >>>>')
