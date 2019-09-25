
## SETUP ################################################################################
stratum=as.character(commandArgs()[4])
run_date=as.character(commandArgs()[5])
indicator=as.character(commandArgs()[6])
indicator_group=as.character(commandArgs()[7])
geos_node = as.logical(commandArgs()[8])
prev_only = as.logical(commandArgs()[9])
rerun = T

interval_mo <- 12

# Define directories
main_dir <- '<<<< FILEPATH REDACTED >>>>'
temp_dir <- '<<<< FILEPATH REDACTED >>>>'

# Load objects from convenience temp file
load('<<<< FILEPATH REDACTED >>>>')

# Print some settings to console
message(names(rake_targets))
message(indicator)
message(indicator_group)
message(run_date)
message(stratum)
message(pop_measure)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

rake_targets = rake_targets[c('prevalence' ,'incidence', 'mortality')]


# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
reg <- stratum

## Load libraries and miscellaneous MBG project functions.
repo = '/ihme/code/geospatial/lbd_core/'

## drive locations
root           <- '<<<< FILEPATH REDACTED >>>>'
package_lib    <- '<<<< FILEPATH REDACTED >>>>'
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'


#get slots
slots = as.numeric(Sys.getenv('NSLOTS'))
cores_to_use = ifelse(grepl('Intel', system("cat /proc/cpuinfo | grep \'name\'| uniq", inter = T)), floor(slots * .86), floor(slots*.64))


## Load libraries and  MBG project functions.
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>')))

if(Sys.info()[1] == 'Windows'){
  stop('STOP! you will overwrite these packages if you run from windows\n
       STOP! also, lots of this functions wont work so get on the cluster!')
} else {
  for(package in package_list)
    require(package, character.only=TRUE)
  for(funk in list.files(repo, recursive=TRUE,pattern='functions',full.names = T)){
    if(length(grep('central',funk))!=0){
      #message(paste0('loading ',funk))
      source(funk)
    }
  }
}

# Custom load indicator-specific functions

#Load raking functions
source('<<<< FILEPATH REDACTED >>>>')

## PREPARE RASTERS, ETC. ################################################################

# Load cell draws
message('Loading Data...')
if(rerun) load(paste0('<<<< FILEPATH REDACTED >>>>'))

message('Loading Map Stuff...')
# Get aggregated estimates for all admin0. Aggregate to level you rake to
simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = FALSE)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]

if(class(year_list)=='character'){
  year_list = eval(parse(text = year_list))
}

## Pull 2000-2015 annual population brick using new covariates function
pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',
                                                     measures = pop_measure,
                                                     simple_polygon = simple_polygon,
                                                     start_year  = min(year_list),
                                                     end_year    = max(year_list),
                                                     interval_mo = as.numeric(interval_mo),
                                                     agebin=1)$worldpop
#make sure the dimensions are the same
pop_raster_annual = raster::crop(pop_raster_annual, simple_raster)

cell_pred_to_raster = function(cell_pred, summstat, raked, measure, ...){
  ras <- make_cell_pred_summary(draw_level_cell_pred = cell_pred, mask = simple_raster, #simple raster comes from the parent environment
                                return_as_raster = TRUE, summary_stat = summstat, ...)
  save_post_est(ras,'raster',paste0(reg,'_', measure, '_', ifelse(raked == "raked", "_raked", ""),'_',summstat,'_raster'))
}



rake_me = function(rrr, cell_pred, weight_brick, logit_rake = F, return_raked = F, suffix = "", calc_summaries = T, cores = 1, fractional_pixel = F){

  message(rrr)

  rake_method = ifelse(logit_rake, 'logit','linear')
  #year list and rake targets are borrowed from the upper environment
  rf = calculate_raking_factors(cell_pred = cell_pred, rake_targets = rake_targets[[rrr]],
                                lyv = c('name','year','mean'), year_list = year_list,
                                simple_raster = simple_raster,
                                weight_brick = weight_brick,
                                rake_method = rake_method,
                                fractional_pixel = fractional_pixel)

  index = ifelse(fractional_pixel, 2, 1 )

  #rake the predictions
  raked_cell_pred <- apply_raking_factors(cell_pred = cell_pred,
                                          simple_raster = simple_raster,
                                          rake_dt = rf[[index]], logit_rake = logit_rake,
                                          fractional_pixel = fractional_pixel,
                                          force_simp_ras_dt_match=F)
  #save results
  save_post_est(rf[[1]], "csv", paste0(reg,'_',rrr,suffix, "_rf"))

  saveRDS(raked_cell_pred, file = '<<<< FILEPATH REDACTED >>>>')

  #save some raked summary estimates
  if(calc_summaries){
    summ_list <- summstats[summstats != "p_below"]
    mclapply(summ_list, function(x) cell_pred_to_raster(raked_cell_pred, x, 'raked', rrr), mc.cores =cores)
  }
  gc()

  if(return_raked){
    return(raked_cell_pred)
  } else {
    invisible()
  }


}

#Deal with prevalence
#save results for unraked

#calc linear raked prev
summ_list <- summstats[summstats != "p_below"]
#calculate logit raked prevalence
if(rerun){
  lapply(summ_list, function(x) cell_pred_to_raster(cell_pred, x, 'unraked', 'prevalence'))

  raked_prev = rake_me(rrr = 'prevalence', cell_pred = cell_pred, weight_brick = pop_raster_annual, logit_rake = T,
                       return_raked = T, suffix ='', calc_summaries = T, cores = 1)
}else{
  raked_prev = readRDS('<<<< FILEPATH REDACTED >>>>')
  lapply(summ_list, function(x) cell_pred_to_raster(raked_prev, x, 'raked', 'prevalence'))
}

if(!prev_only){
  #remove regular cell pred to save memory
  rm(cell_pred)

  #calculate for other metrics

  rakes = names(rake_targets)
  rakes = rakes[!rakes %in% 'prevalence']
  for(rake_measure in rakes){

    if(!rerun){
      raked_prev = readRDS('<<<< FILEPATH REDACTED >>>>')
      lapply(summ_list, function(x) cell_pred_to_raster(raked_prev, x, 'raked', rake_measure))

    }else{
      message('no re re')
      rake_me(rrr = rake_measure, cell_pred = raked_prev, weight_brick = pop_raster_annual, logit_rake = F, return_raked = F,
              suffix ='', calc_summaries = T, cores = 1)
    }

  }
}
# All done
message(paste0("Done with post-estimation for ", stratum))
