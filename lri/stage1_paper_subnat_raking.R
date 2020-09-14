########################################################################
#Post-hoc subnational raking for LRI stage 1 paper
#to GBD 2017 results
#######################################################################

#(1) Setup ############################################################
rm(list = ls())

#user arguments
rake_method <- 'logit' #logit or linear
run_date <- '2017_stage_1_updates'
subnats <- c('ZAF','KEN','ETH')
measures <- c('incidence','mortality')
year_list <- c(2000:2017)
modeling_shapefile_version <- '2018_12_04'
raking_shapefile_version <- '2018_12_04'
mort_rasters_dir <- '<<<< FILEPATH REDACTED >>>>'
prev_rasters_dir <- '<<<< FILEPATH REDACTED >>>>'
inc_rasters_dir <- '<<<< FILEPATH REDACTED >>>>'
save_dir <- '<<<< FILEPATH REDACTED >>>>'

source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
source('<<<< FILEPATH REDACTED >>>>/get_location_metadata.R')
source('<<<< FILEPATH REDACTED >>>>/mbg_central/fractional_raking_functions.R')
if (rake_method == 'logit') source('<<<< FILEPATH REDACTED >>>>/lri/logit_rake_subnationals_functions.R')

## drive locations
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source('<<<< FILEPATH REDACTED >>>>/mbg_central/setup.R')
mbg_setup(package_list = package_list, repos = '<<<< FILEPATH REDACTED >>>>')
library(mgcv)

## source any functions with custom edits from lbd_custom folder
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/misc_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/post_estimation_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/stacking_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/fractional_raking_functions.R')
source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')


#(2) Read in LBD results ##############################################
lbd_prev <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(ADM0_NAME %in% c('South Africa','Kenya','Ethiopia'))
lbd_mort <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(ADM0_NAME %in% c('South Africa','Kenya','Ethiopia'))
lbd_inc <- fread('<<<< FILEPATH REDACTED >>>>') %>%
  filter(ADM0_NAME %in% c('South Africa','Kenya','Ethiopia'))

#load the raking raster
shapefile <- readRDS('<<<< FILEPATH REDACTED >>>>')
link <- readRDS('<<<< FILEPATH REDACTED >>>>')

new_simple_polygon <- load_simple_polygon(
  gaul_list = NULL,
  buffer = 0.4, custom_shapefile = shapefile)

new_subset_shape <- new_simple_polygon[["subset_shape"]]
new_simple_polygon <- new_simple_polygon[["spoly_spdf"]]

message("Loading simple raster\n")
new_raster_list <- build_simple_raster_pop(new_subset_shape, field = 'ADM1_CODE') #custom function that makes a full pops raster
new_simple_raster <- new_raster_list[["simple_raster"]]
pop_raster <- new_raster_list[['pop_raster']]

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
    rasters_dir <- inc_rasters_dir
  } else if (measure == 'prevalence'){
    measure_id <- 5
    lbd_results <- lbd_prev
    rasters_dir <- prev_rasters_dir
  } else if (measure == 'mortality'){
    measure_id <- 1
    lbd_results <- lbd_mort
    rasters_dir <- mort_rasters_dir
  }
  
  gbd_results <- get_outputs(topic= 'cause',
                             cause_id=322,
                             location_id=subnat_loc_ids,
                             year_id=year_list,
                             age_group_id=1,
                             gbd_round_id=5,
                             metric_id=3,
                             measure_id = measure_id,
                             sex_id = 3,
                             version='latest') %>%
    
    as.data.frame() %>%
    group_by(location_id) %>%
    mutate(inter_val = as.numeric(try(predict(gam(val ~ s(year_id)), newdata = data.frame(year_id = year_list)))))
  
  gbd_results$val <- ifelse(is.na(gbd_results$val), gbd_results$inter_val, gbd_results$val)
  gbd_results <- gbd_results[complete.cases(gbd_results),]
  
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
  
  #(4) Calculate raking factors ##############################################
  compare_results <- merge(lbd_results, gbd_results, by = c('ADM1_NAME','year'), all.x = TRUE, all.y = TRUE) %>%
    as.data.table()
  
  #check success of merge
  unique(filter(compare_results, is.na(gbd_mean))$ADM1_NAME)
  unique(filter(compare_results, is.na(mean))$ADM1_NAME)
  
  #define which subnats to rake
  ad1_subnat <- unique(compare_results$ADM1_CODE)[!is.na(unique(compare_results$ADM1_CODE))]
  no_ad1_subnat <- unique(filter(link, !(ADM1_CODE %in% ad1_subnat))$ADM1_CODE)
  
  #caculate raking factors: linear
  if (rake_method == 'linear'){
    message('Using linear raking')
    compare_results[,rf := gbd_mean/mean]
    compare_results[, `:=`(mean = mean*rf, upper = upper*rf, lower = lower * rf)]
  }
  
  #calculate raking factors: logit
  if (rake_method == 'logit'){
    message('Using logit raking')
    
    #this is not complete and will break!
    length(values(pop))
    length(values(original))
    length(values(new_simple_raster))
    
    #extract to vectors
    pop_vector <- values(pop)
    lbd_vector <- values(original)
    admin_code_vector <- values(new_simple_raster)
    table <- data.table(pop = pop_vector, lbd = lbd_vector, admin_code = admin_code_vector)
    
    #subset to ad1 in question
    ad1_code <- 1072
    ad1_table <- filter(table, admin_code == ad1_code)
    
    #subset pop and lbd to ad1
    p_i <- ad1_table$pop #populations
    p_N <-compare_results[ADM1_CODE == ad1_code, gbd_mean] #GDB estimates
    N_i <- ad1_table$lbd #raster estimates
    a <- 5 #how close
    
    results <- FindK(p_i, p_N, N_i, a)
  }
  
  #save out new csvs
  lbd_full <- fread('<<<< FILEPATH REDACTED >>>>') %>%
    select(-cirange) %>%
    filter(!(ADM0_NAME %in% c('Ethiopia','Kenya','South Africa')))
  
  results_to_bind <- filter(compare_results, !is.na(ADM1_CODE)) %>%
    select(-cirange, -location_id, -gbd_mean, -gbd_upper, -gbd_lower, -rf)
  
  raked_csv <- rbind(lbd_full, results_to_bind)
  write.csv(raked_csv, paste0(save_dir, 'has_lri_', measure, '_admin_1_raked_summary.csv'))
  
  #(5) Create a raking raster to multiply results rasters by ###################

  #replace admin1s w/ subnat raking with value = rf
  rf_raster <- new_simple_raster
  rf_raster[rf_raster %in% no_ad1_subnat] <- 1
  
  rf_raster <- stack(replicate(18, rf_raster))
  for (n in 1:18){
    year_raster <- rf_raster[[n]]
    yr <- 1999 + n
    for (ad1 in ad1_subnat){
      year_raster[year_raster == ad1] <- compare_results[ADM1_CODE == ad1 & year == yr, rf]
    }
    rf_raster[[n]] <- year_raster
  }
  
  #(6) Create raked results  #################
  
  #load original results rasters
  message(rasters_dir)
  
  for (stat in c('mean','upper','lower')){
    
    #load
    original_file <- paste0(rasters_dir, '/originals/has_lri_', measure, '_', stat, '_raked_2000_2017.tif')
    if (stat == 'mean' & measure == 'incidence') original_file <- paste0(rasters_dir, '/has_lri_', measure, '_', stat, '_raked_2000_2017.tif')
    original <- stack(original_file)
    
    #fix extents
    rf_raster <- extend(rf_raster, original, values = NA)
    rf_raster <- crop(rf_raster, extent(original))
    
    #multiply results rasters by raking factor raster
    raked <- rf_raster * original
    writeRaster(raked, paste0(save_dir, 'has_lri_', measure, '_', stat, '_raked_2000_2017.tif'), overwrite = TRUE)
  }
  
  # (7) Create raked draw files
  load('<<<< FILEPATH REDACTED >>>>')
  ad1_draws <- admin_1
  rm(admin_1, admin_0, admin_2)
  rf_link <- select(compare_results, year, ADM1_CODE, rf)
  ad1_draws <- merge(ad1_draws, rf_link, by = c('ADM1_CODE','year'), all.x = TRUE)
  ad1_draws$rf[is.na(ad1_draws$rf)] <- 1
  
  col_names <- names(select(ad1_draws, starts_with("V")))
  raked_results <- mutate_at(ad1_draws, .vars = col_names, funs(.*rf))
  admin_1 <- raked_results
  
  save(admin_1, file = paste0(save_dir, 'has_lri_', measure, '_raked_admin_draws_eb_bin0_0.RData'))
}