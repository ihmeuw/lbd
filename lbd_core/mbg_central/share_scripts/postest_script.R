# Parallel script for post-estimation
#
# Note that prep_postest is used to pass objects to the parallel script
# (saved in a temporarly RData file).  If you want to use logit raking,
# pass "rake_transform" as a string that equals "logit".  You can set this
# up easily in the config file with a row for rake_transform with value "logit"
#
# Sample use:
#
# # Prepare for parallel post-estimation - save file with objs to re-load
# prep_postest(indicator = indicator,
#              indicator_group = indicator_group,
#              run_date = run_date,
#              save_objs = c("core_repo", "gbd", "year_list", "summstats", "rake_transform"))
#
# ## Parallelized post-estimation over region
# postest_script <- "postest_script"
#
# for (s in strata) {
#  qsub <- make_qsub_postest(code = postest_script,
#                            stratum = s,
#                            log_location = 'sharedir',
#                            memory        = mem,
#                            cores         = cores)
#  system(qsub)
# }
#
# ## check to make sure post-est done before continuing
# waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

## SETUP ################################################################################

stratum=as.character(commandArgs()[4])
run_date=as.character(commandArgs()[5])
indicator=as.character(commandArgs()[6])
indicator_group=as.character(commandArgs()[7])
geos_node = as.logical(commandArgs()[8])
reg = as.character(commandArgs()[9])
interval_mo <- 12

# Define directories
main_dir <- paste0('<<< FILEPATH REDACTED >>>>')
temp_dir <- paste0(main_dir, "temp_post_est/")

# Load objects from convenience temp file
load(paste0(temp_dir, "post_est_temp_objs.RData"))

# Print some settings to console
if(grepl("rescaled", indicator)) {indicator <- substr(indicator, 1, nchar(indicator)-9)}
message(main_dir)
message(indicator)
message(indicator_group)
message(run_date)
message(stratum)
message(pop_measure)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

user <- Sys.info()['user']
message(user)
message(sessionInfo())
core_repo <- paste0('<<< FILEPATH REDACTED >>>>')
indic_repo <- paste0('<<< FILEPATH REDACTED >>>>')
indicator_group <- 'education'
rerun_run_date <- NA

root           <- ifelse(Sys.info()[1]=='Windows', ''<<< FILEPATH REDACTED >>>>'/', ''<<< FILEPATH REDACTED >>>>'')
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                         sprintf('<<< FILEPATH REDACTED >>>>'),
                         sprintf('<<< FILEPATH REDACTED >>>>')
#'<<< FILEPATH REDACTED >>>>')
sharedir       <- sprintf('<<< FILEPATH REDACTED >>>>')
commondir      <- sprintf('<<< FILEPATH REDACTED >>>>')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
load_R_packages(c('data.table', 'assertthat'))
mbg_setup(package_list = package_list, repos = core_repo, load_from_sing = TRUE)

# #delete oceania from regions since it has no data and didn't run
# Regions <- c("latin_america", "asia_png", "south_asia", "name", "essa", "wssa", "cssa", "sssa", "middle_east")
# strata <- Regions


# #CUSTOM RAKING
if(indicator=='edu_mean_15_49_male' | indicator == 'edu_mean_20_24_male') {
  raking_sex <- 1
  raking_pop <- 'total'
}
if(indicator=='edu_mean_15_49_female' | indicator=='edu_mean_20_24_female') {
  raking_sex <- 2
  raking_pop <- 'wocba'
}
# # 
# # 
if(indicator=='edu_mean_15_49_female' | indicator == 'edu_mean_15_49_male') {
  year_ids  <- c(min(year_list):max(year_list))
  gaul_list <- get_gaul_codes(stratum)
  # #   
  library('dplyr')
  metadata <- get_covariate_metadata()
  source('<<< FILEPATH REDACTED >>>>')
  #gbd_estimates <- get_covariate_estimates(covariate_name_short = 'education_yrs_pc')
  gbd_estimates <- get_covariate_estimates(covariate_id = metadata[covariate_name_short == 'education_yrs_pc', covariate_id], gbd_round_id = 5, year_id = paste(as.character(year_ids), collapse = " "))
  
  #   # Merge pops to collapse age-specific to 10-54 women
  source('<<< FILEPATH REDACTED >>>>')
  gbd_pops <- get_population(age_group_id = paste(as.character(unique(gbd_estimates[, age_group_id])), collapse = " "),
                             sex_id = "1 2",
                             location_id = paste(as.character(unique(gbd_estimates[, location_id])), collapse = " "),
                             year_id = paste(as.character(unique(gbd_estimates[, year_id])), collapse = " "))
  gbd_pop_estimates <- merge(gbd_estimates, gbd_pops, by=c('age_group_id','sex_id','location_id','year_id'))
  # Subset to 10-54 women
  gbd_pop_estimates <- gbd_pop_estimates[age_group_id %in% c(7:15) & sex_id == raking_sex, ]
  # Pop-weighted collapse to mean for 10-54
  gbd_means <- gbd_pop_estimates[,.(gbd_mean=weighted.mean(x=mean_value,w=population)), by=.(location_id,year_id)]
  # Convert to GAUL_CODE
  gbd_estimates <- gbd_means
  #gaul_to_loc_id <- fread('<<< FILEPATH REDACTED >>>>')
  #names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"
  gbd_estimates <- gbd_estimates[year_id %in% year_ids,]
  # gbd_estimates <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
  # gbd_estimates <- gbd_estimates[GAUL_CODE %in% gaul_list,]
  names(gbd_estimates)[names(gbd_estimates)=="location_id"] <- "name"   
  names(gbd_estimates)[names(gbd_estimates)=="year_id"] <- "year"
  names(gbd_estimates)[names(gbd_estimates)=="gbd_mean"] <- "mean"
  gbd_estimates <- gbd_estimates[, c('name', 'year', 'mean'), with = FALSE]
  gbd <- gbd_estimates
}
if(indicator=='edu_mean_20_24_female' | indicator=='edu_mean_20_24_male') {
  year_ids  <- c(min(year_list):max(year_list))
  gaul_list <- get_gaul_codes(stratum)
  # #   
  library('dplyr')
  metadata <- get_covariate_metadata()
  source('<<< FILEPATH REDACTED >>>>')
  #gbd_estimates <- get_covariate_estimates(covariate_name_short = 'education_yrs_pc')
  gbd_estimates <- get_covariate_estimates(covariate_id = metadata[covariate_name_short == 'education_yrs_pc', covariate_id], gbd_round_id = 5, year_id = paste(as.character(year_ids), collapse = " "))
  
  #   # Merge pops to collapse age-specific to 10-54 women
  source('<<< FILEPATH REDACTED >>>>')
  gbd_pops <- get_population(age_group_id = paste(as.character(unique(gbd_estimates[, age_group_id])), collapse = " "),
                             sex_id = "1 2",
                             location_id = paste(as.character(unique(gbd_estimates[, location_id])), collapse = " "),
                             year_id = paste(as.character(unique(gbd_estimates[, year_id])), collapse = " "))
  gbd_pop_estimates <- merge(gbd_estimates, gbd_pops, by=c('age_group_id','sex_id','location_id','year_id'))
  # Subset to 10-54 women
  gbd_pop_estimates <- gbd_pop_estimates[age_group_id == 9  & sex_id == raking_sex, ]
  # Pop-weighted collapse to mean for 10-54
  gbd_means <- gbd_pop_estimates[,.(gbd_mean=weighted.mean(x=mean_value,w=population)), by=.(location_id,year_id)]
  # Convert to GAUL_CODE
  gbd_estimates <- gbd_means
  #gaul_to_loc_id <- fread('<<< FILEPATH REDACTED >>>>')
  #names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"
  gbd_estimates <- gbd_estimates[year_id %in% year_ids,]
  # gbd_estimates <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
  # gbd_estimates <- gbd_estimates[GAUL_CODE %in% gaul_list,]
  names(gbd_estimates)[names(gbd_estimates)=="location_id"] <- "name"   
  names(gbd_estimates)[names(gbd_estimates)=="year_id"] <- "year"
  names(gbd_estimates)[names(gbd_estimates)=="gbd_mean"] <- "mean"
  gbd_estimates <- gbd_estimates[, c('name', 'year', 'mean'), with = FALSE]
  gbd <- gbd_estimates
}
#   gaul_list <- get_gaul_codes(model_domain)
#   
#   library('<<< FILEPATH REDACTED >>>>')
#   metadata <- get_covariate_metadata()
#   source('<<< FILEPATH REDACTED >>>>')
#   #gbd_estimates <- get_covariate_estimates(covariate_name_short = 'education_yrs_pc')
#   gbd_estimates <- get_covariate_estimates(covariate_id = metadata[covariate_name_short == 'education_yrs_pc', covariate_id], gbd_round_id = 4)
#   
#   # Merge pops to collapse age-specific to 10-54 women
#   source('<<< FILEPATH REDACTED >>>>')
#   gbd_pops <- get_population(age_group_id = paste(as.character(unique(gbd_estimates[, age_group_id])), collapse = " "),
#                              sex_id = "1 2",
#                              location_id = paste(as.character(unique(gbd_estimates[, location_id])), collapse = " "),
#                              year_id = paste(as.character(unique(gbd_estimates[, year_id])), collapse = " "))
#   gbd_pop_estimates <- merge(gbd_estimates, gbd_pops, by=c('age_group_id','sex_id','location_id','year_id'))
#   # Subset to 10-54 women
#   gbd_pop_estimates <- gbd_pop_estimates[age_group_id == 9 & sex_id == raking_sex, ]
#   # Pop-weighted collapse to mean for 10-54
#   gbd_means <- gbd_pop_estimates[,.(gbd_mean=weighted.mean(x=mean_value,w=population)), by=.(location_id,year_id)]
#   # Convert to GAUL_CODE
#   gbd_estimates <- gbd_means
#   gaul_to_loc_id <- fread('<<< FILEPATH REDACTED >>>>')
#   names(gaul_to_loc_id)[names(gaul_to_loc_id)=="loc_id"] <- "location_id"
#   gbd_estimates <- gbd_estimates[year_id %in% year_ids,]
# gbd_estimates <- merge(gaul_to_loc_id, gbd_estimates, by="location_id")
#    gbd_estimates <- gbd_estimates[GAUL_CODE %in% gaul_list,]
#    names(gbd_estimates)[names(gbd_estimates)=="GAUL_CODE"] <- "name"
#   names(gbd_estimates)[names(gbd_estimates)=="year_id"] <- "year"#
#   names(gbd_estimates)[names(gbd_estimates)=="gbd_mean"] <- "mean"
#    gbd_estimates <- gbd_estimates[, c('name', 'year', 'mean'), with = FALSE]
#    gbd <- gbd_estimates
#  }
# 
# #postest script
# if(grepl('prop',indicator)) gbd <- NULL
# summstats <- 'mean'
# prep_postest(indicator = indicator,
#              indicator_group = indicator_group,
#              run_date = run_date,
#              save_objs = c("core_repo", "gbd", "year_list", "summstats", 
#                            "rake_transform", "pop_measure"))
# 
# # For now just assign stratum to reg (will need to modify the below for strata beyond reg)
# reg <- stratum
# 
# ## Load libraries and miscellaneous MBG project functions.
# setwd(core_repo)
# root <- ifelse(Sys.info()[1]=="Windows", "'<<< FILEPATH REDACTED >>>>'/", "'<<< FILEPATH REDACTED >>>>'")
# 
# ## drive locations
# root           <- ifelse(Sys.info()[1]=='Windows', ''<<< FILEPATH REDACTED >>>>'/', ''<<< FILEPATH REDACTED >>>>'')
# package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
#                          paste0('<<< FILEPATH REDACTED >>>>'),
#                          paste0('<<< FILEPATH REDACTED >>>>')
# sharedir       <- sprintf('<<< FILEPATH REDACTED >>>>')
# commondir      <- sprintf('<<< FILEPATH REDACTED >>>>')
# 
# 
# ## Load libraries and  MBG project functions.
# .libPaths(package_lib)
# package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
# 
# if(Sys.info()[1] == 'Windows'){
#   stop('STOP! you will overwrite these packages if you run from windows\n
#        STOP! also, lots of this functions wont work so get on the cluster!')
# } else {
#   for(package in package_list)
#     require(package, lib.loc = package_lib, character.only=TRUE)
#   for(funk in list.files(recursive=TRUE,pattern='functions')){
#     if(length(grep('central',funk))!=0){
#       message(paste0('loading ',funk))
#       source(funk)
#     }
#   }
# }
## PREPARE RASTERS, ETC. ################################################################

# Load cell draws
message('Loading Data...')
# if(grepl('prop', indicator) & !grepl('zero', indicator)){
#   cell_pred <- load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
#   cell_pred<- get(cell_pred)
# }
# if(grepl('mean', indicator) | grepl('zero', indicator)){
load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
# }
if(grepl('mean', indicator)) cell_pred[cell_pred <0] <- 0
message('Loading Map Stuff...')
# Get aggregated estimates for all admin0. Aggregate to level you rake to
r <- reg
if(reg == 'trsa') r <- 'trsa-pry-guf'
if(reg == 'eaas') r <- 'eaas-prk'
simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(r), buffer = 0.4, subset_only = FALSE)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]
# 
# ## Pull 2000-2015 annual population brick using new covariates function
# pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',
#                                                      measures = pop_measure,
#                                                      simple_polygon = simple_polygon,
#                                                      start_year  = min(year_list),
#                                                      end_year    = max(year_list),
#                                                      interval_mo = as.numeric(interval_mo),
#                                                      agebin=1)
# pop_raster_annual  <- pop_raster_annual[[1]]
# pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
# pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
# pop_raster_annual  <- mask(pop_raster_annual, simple_raster)

## POPULATION WEIGHTS ###################################################################

## Create population weights using the annual brick and feed custom year argument to aggregation function
# pop_wts_adm0 <- make_population_weights(admin_level   = 0,
#                                         simple_raster = simple_raster,
#                                         pop_raster    = pop_raster_annual,
#                                         gaul_list     = get_gaul_codes(reg))
# 
# cond_sim_raw_adm0 <- make_condSim(pop_wts_object = pop_wts_adm0,
#                                   gaul_list      = get_gaul_codes(reg),
#                                   admin_level    = 0,
#                                   cell_pred      = cell_pred,
#                                   summarize      = TRUE,
#                                   years          = year_list)

if (!is.null(gbd)) {
  
  message("Raking...")
  
  ## RAKING ###############################################################################
  
  rake_to <- gbd
  ifelse(reg %in% c('ansa+pry+sur+guy+guf', 'noaf'), cross <- T, cross <- F) 
  #run raking
  if(reg %in% c('trsa', 'eaas')){
    raked_outputs <- rake_cell_pred(cell_pred,
                                    rake_to,                 #dataframe or table with columns name, year, mean where name is the admin0/1 gbd location id
                                    reg,
                                    year_list,
                                    pop_measure,             
                                    rake_method = rake_transform,   #linear or logit
                                    crosswalk = cross,           #keep true if the model was run before 8/3/18 when the shapefiles changed
                                    rake_subnational = F,    #T or F
                                    MaxJump = 10,            #optional param for logit raking
                                    MaxIter = 80,            #optional param for logit raking
                                    FunTol = 1e-5,           #optional param for logit raking
                                    iterate = F,             #optional param for logit raking, try larger jumps and iter if not converging
                                    if_no_gbd = "return_unraked",
                                    shapefile_path = get_admin_shapefile(0, raking = T))
  } else{
    raked_outputs <- rake_cell_pred(cell_pred,
                                    rake_to,                 #dataframe or table with columns name, year, mean where name is the admin0/1 gbd location id
                                    reg,
                                    year_list,
                                    pop_measure,             
                                    rake_method = rake_transform,   #linear or logit
                                    crosswalk = cross,           #keep true if the model was run before 8/3/18 when the shapefiles changed
                                    rake_subnational = F,    #T or F
                                    MaxJump = 10,            #optional param for logit raking
                                    MaxIter = 80,            #optional param for logit raking
                                    FunTol = 1e-5,           #optional param for logit raking
                                    iterate = F,             #optional param for logit raking, try larger jumps and iter if not converging
                                    if_no_gbd = "return_unraked",
                                    shapefile_path = get_admin_shapefile(0, raking = T),
                                    simple_polygon = simple_polygon,
                                    simple_raster = simple_raster,
                                    pop_raster = pop_raster) #"return_unraked" or "return_na" for admin units without GBD raking targets
  }
  raked_cell_pred <- raked_outputs[[1]]
  simple_raster <- raked_outputs[[2]]
  rf <- raked_outputs[[3]]
  ## Get raking factors
  # if (rake_transform == "logit") {
  #   
  #   rf <- calc_logit_raking_factors(pop_wts_object = pop_wts_adm0,
  #                                   gaul_list      = get_gaul_codes(reg),
  #                                   admin_level    = 0,
  #                                   cell_pred      = cell_pred,
  #                                   rake_to        = gbd,
  #                                   years           = year_list)
  # 
  #   raked_cell_pred <- rake_predictions(raking_factors = rf,
  #                                       pop_wts_object = pop_wts_adm0,
  #                                       cell_pred      = cell_pred,
  #                                       logit_rake     = TRUE)
  #   
  #   ## Get rf into standard format if logit raking - get aggregated values to append to rf
  #   aggs <- get_aggs(agg_geo_est = cond_sim_raw_adm0,
  #                    gaul_list   = get_gaul_codes(reg),
  #                    rake_to     = gbd)
  #   
  #   aggs$name <- as.character(aggs$name)
  #   aggs$year <- as.character(aggs$year)
  #   rf <- merge(rf, aggs, by=c("name", "year"))
  #   
  # } else {
  #   
  #   rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
  #                               gaul_list   = get_gaul_codes(reg),
  #                               rake_to     = gbd)
  #   
  #   raked_cell_pred <- rake_predictions(raking_factors = rf,
  #                                       pop_wts_object = pop_wts_adm0,
  #                                       cell_pred      = cell_pred,
  #                                       logit_rake     = FALSE)
  # }
  # 
  
  
} else if (is.null(gbd)) {
  
  # For unraked models
  rf                <- NULL
  raked_cell_pred   <- NULL
  
}

## SAVE THE RESULTS #####################################################################
message('Saving results...')

## save RF
save_post_est(rf, "csv", paste0(reg, "_rf"))

## save raked cell preds
save(raked_cell_pred, file = paste0(sharedir, "/output/", run_date, "/",
                                    indicator, "_raked_cell_draws_eb_bin0_", reg, "_0.RData" ))

# make and save summaries
if(grepl('mean', indicator)){
  save_cell_pred_summary <- function(summstat, raked, ...) {
    message(paste0('Making unraked summmary raster for: ',summstat, " (", raked, ")"))
    if (raked=="unraked") cpred <- "cell_pred"
    if (raked=="raked")   cpred <- "raked_cell_pred"
    ras <- make_cell_pred_summary(
      draw_level_cell_pred = get(cpred),
      mask                 = simple_raster,
      return_as_raster     = TRUE,
      summary_stat         = summstat,
      ...)
    save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ""),'_',summstat,'_raster'))
  }
  
  # Do this as lapply to not fill up memory in global env with big obs
  if (is.null(gbd))  rake_list <- c("unraked")
  if (!is.null(gbd)) rake_list <- c("unraked", "raked")
  
  summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)
  
  lapply(1:nrow(summ_list), function(i) {
    summstat <- as.character(summ_list[i, 1])
    raked <- as.character(summ_list[i, 2])
    save_cell_pred_summary(summstat, raked)
  })
  
  ## Can't pass additional params in the above framework, so will code by hand here
  for (r in rake_list) {
    if ("p_below" %in% summstats) {
      save_cell_pred_summary(summstat = "p_below",
                             raked = r,
                             value = 0.8,
                             equal_to = F)
    }
  }
}

# Write a file to mark done
output_dir <- paste0('<<< FILEPATH REDACTED >>>>')
pathaddin <- paste0('_bin0_',reg,'_0') # To allow us to use waitformodelstofinish()
write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

# All done
message(paste0("Done with post-estimation for ", stratum))

