########################################################################
#Compare aggregated ad1 results for Stage 1 LRI paper
#######################################################################

#(1) Setup ############################################################
rm(list = ls())

#user arguments
run_date <- '2017_stage_1_updates'
subnats <- c('ZAF','KEN','ETH')
measures <- c('prevalence','incidence','mortality')

source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
source('<<<< FILEPATH REDACTED >>>>/get_location_metadata.R')
library(ggplot2)

#(2) Read in LBD results ##############################################
lbd_prev <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(ADM0_NAME %in% c('South Africa','Kenya','Ethiopia'))
lbd_mort <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(ADM0_NAME %in% c('South Africa','Kenya','Ethiopia'))
lbd_inc <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(ADM0_NAME %in% c('South Africa','Kenya','Ethiopia'))

#(3) Read in GBD results ##############################################
# Subset locations to National and Subnational estimates
locs <- read.csv('<<<< FILEPATH REDACTED >>>>')

subnat_string <- paste(subnats, collapse = '|')
ihme_loc_ids <- locs$ihme_loc_id
subnat_ihme_loc_id <- grep(subnat_string, ihme_loc_ids, value = T)

subnat_loc_ids <- filter(locs, ihme_loc_id %in% subnat_ihme_loc_id)$location_id

for (measure in measures){
  
  print(measure)
  
  if (measure == 'incidence'){
    measure_id <- 6
    lbd_results <- lbd_inc
  } else if (measure == 'prevalence'){
    measure_id <- 5
    lbd_results <- lbd_prev
  } else if (measure == 'mortality'){
    measure_id <- 1
    lbd_results <- lbd_mort
  }
  
  gbd_results <- get_outputs(topic= 'cause',
                          cause_id=322,
                          location_id=subnat_loc_ids,
                          year_id=c(2000:2017),
                          age_group_id=1,
                          gbd_round_id=5,
                          metric_id=3,
                          measure_id = measure_id,
                          sex_id = 3,
                          version='latest')
  
  gbd_results <- dplyr::select(gbd_results, location_id, year_id, val, location_name, upper, lower) %>%
    dplyr::rename('gbd_mean' = 'val', 'year' = 'year_id', 'ADM1_NAME' = 'location_name', 'gbd_upper' = 'upper', 'gbd_lower' = 'lower')
  
  #fix gbd loc names to facilitate merge
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'Addis Ababa'] <- 'Addis Abeba'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'Benishangul-Gumuz'] <- 'Benshangul-Gumaz'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'Gambella'] <- 'Gambela Peoples'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'Harari'] <- 'Harari People'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'HomaBay'] <- 'Homa Bay'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'North-West'] <- 'North West'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'Southern Nations, Nationalities, and Peoples'] <- 'Southern Nations, Nationalities and Peoples'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'TaitaTaveta'] <- 'Taita Taveta'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'TanaRiver'] <- 'Tana River'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'TharakaNithi'] <- 'Tharaka-Nithi'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'TransNzoia'] <- 'Trans Nzoia'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'UasinGishu'] <- 'Uasin Gishu'
  gbd_results$ADM1_NAME[gbd_results$ADM1_NAME == 'WestPokot'] <- 'West Pokot'
  
  #($) Merge to lbd and plot ##############################################
  compare_results <- merge(lbd_results, gbd_results, by = c('ADM1_NAME','year'), all.x = TRUE, all.y = TRUE)
  
  #check success of merge
  unique(filter(compare_results, is.na(gbd_mean))$ADM1_NAME)
  unique(filter(compare_results, is.na(mean))$ADM1_NAME)
  
  pdf('<<<< FILEPATH REDACTED >>>>')
  for (country in c('Ethiopia', 'Kenya', 'South Africa')){
    for (yr in 2000:2017){
      dat <- dplyr::filter(compare_results, ADM0_NAME == country & year == yr)
      plot <- ggplot(dat, aes(x = mean, y = gbd_mean, color = year)) + geom_point() + ggtitle(paste(country, yr))
      plot(plot)
    }
  }
  dev.off()
  
  #make time series with UI of GBD, LBD
  pdf('<<<< FILEPATH REDACTED >>>>')
  for (country in c('Ethiopia', 'Kenya', 'South Africa')){
    dat <- filter(compare_results, ADM0_NAME == country)
    admins <- unique(dat$ADM1_NAME)
    for (admin in admins){
      admin_dat <- filter(dat, ADM1_NAME == admin)
      plot <- ggplot(admin_dat, aes(x = year, y = mean)) + geom_point(color = 'springgreen4') + ggtitle(paste0(admin, ' (', country, ')')) + 
        geom_ribbon(data = admin_dat, aes(x = year, ymin = lower, ymax = upper), na.rm = TRUE, alpha = 0.2, fill = 'springgreen4') +
        geom_point(aes(x = year, y = gbd_mean), color = 'grey36') + 
        geom_ribbon(data = admin_dat, aes(x = year, ymin = gbd_lower, ymax = gbd_upper), na.rm = TRUE, alpha = 0.2, fill = 'grey36')
      plot(plot)
    }
  }
  dev.off()
}


