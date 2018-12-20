# Parallel script for post-estimation

## SETUP ################################################################################

stratum=as.character(commandArgs()[3])
run_date=as.character(commandArgs()[4])
indicator=as.character(commandArgs()[5])
indicator_group=as.character(commandArgs()[6])
geos_node = as.logical(commandArgs()[7])
interval_mo <- 12

# Define directories
main_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date, '/')
temp_dir <- paste0(main_dir, "temp_post_est/")

# Load objects from convenience temp file
load(paste0(temp_dir, "post_est_temp_objs.RData"))

# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(stratum)
message(pop_measure)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
reg <- stratum

## Load libraries and miscellaneous MBG project functions.
setwd(repo)
root <- ifelse(Sys.info()[1]=="Windows", "'<<<< FILEPATH REDACTED >>>>>", "'<<<< FILEPATH REDACTED >>>>>")

## drive locations
root           <- ifelse(Sys.info()[1]=='Windows', '<<<< FILEPATH REDACTED >>>>>', '<<<< FILEPATH REDACTED >>>>>')
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                          paste0(root,'temp/geospatial/geos_packages'),
                          paste0(root,'temp/geospatial/packages'))
sharedir       <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s',indicator_group,indicator)
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>>/mbg/common_inputs')


## Load libraries and  MBG project functions.
.libPaths(package_lib)
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

if(Sys.info()[1] == 'Windows'){
  stop('STOP! you will overwrite these packages if you run from windows\n
        STOP! also, lots of this functions wont work so get on the cluster!')
} else {
  for(package in package_list)
    require(package, lib.loc = package_lib, character.only=TRUE)
  for(funk in list.files(recursive=TRUE,pattern='functions')){
    if(length(grep('central',funk))!=0){
        message(paste0('loading ',funk))
        source(funk)
    }
  }
}

## PREPARE RASTERS, ETC. ################################################################

# Load cell draws
message('Loading Data...')
load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))

message('Loading Map Stuff...')
# Get aggregated estimates for all admin0. Aggregate to level you rake to
simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = FALSE)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]

## Pull 2000-2015 annual population brick using new covariates function
pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',                
                                                     measures = pop_measure,          
                                                     simple_polygon = simple_polygon,
                                                     start_year  = min(year_list),
                                                     end_year    = max(year_list),
                                                     interval_mo = as.numeric(interval_mo),
                                                     agebin=1)
pop_raster_annual  <- pop_raster_annual[[1]]
pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
pop_raster_annual  <- mask(pop_raster_annual, simple_raster)

## POPULATION WEIGHTS ###################################################################

## Create population weights using the annual brick and feed custom year argument to aggregation function
pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                        simple_raster = simple_raster,
                                        pop_raster    = pop_raster_annual,
                                        gaul_list     = get_gaul_codes(reg))

cond_sim_raw_adm0 <- make_condSim(pop_wts_object = pop_wts_adm0,
                                  gaul_list      = get_gaul_codes(reg),
                                  admin_level    = 0,
                                  cell_pred      = cell_pred,
                                  summarize      = TRUE,
                                  years          = year_list)

if (!is.null(gbd)) {

  message("Raking...")

  ## RAKING ###############################################################################

  ## Get raking factors
  if (rake_transform == "logit") {
  
    rf <- calc_logit_raking_factors(pop_wts_object = pop_wts_adm0,
                                    gaul_list      = get_gaul_codes(reg),
                                    admin_level    = 0,
                                    cell_pred      = cell_pred,
                                    rake_to        = gbd,
                                    year           = 2000:2015)
      
    raked_cell_pred <- rake_predictions(raking_factors = rf,
                                    pop_wts_object = pop_wts_adm0,
                                    cell_pred      = cell_pred,
                                    logit_rake     = TRUE)

    ## Get rf into standard format if logit raking - get aggregated values to append to rf
    aggs <- get_aggs(agg_geo_est = cond_sim_raw_adm0,
                     gaul_list   = get_gaul_codes(reg),
                     rake_to     = gbd)

    aggs$name <- as.character(aggs$name)
    aggs$year <- as.character(aggs$year)
    rf <- merge(rf, aggs, by=c("name", "year"))

  } else {
  
    rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
                                gaul_list   = get_gaul_codes(reg),
                                rake_to     = gbd)  

    raked_cell_pred <- rake_predictions(raking_factors = rf,
                                    pop_wts_object = pop_wts_adm0,
                                    cell_pred      = cell_pred,
                                    logit_rake     = FALSE)
  }



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

# Write a file to mark done
output_dir <- paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator, '/output/', run_date)
pathaddin <- paste0('_bin0_',reg,'_0') # To allow us to use waitformodelstofinish()
write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

# All done
message(paste0("Done with post-estimation for ", stratum))
