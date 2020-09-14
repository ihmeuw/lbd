######################################################################
## Script to split raster-level MBG outputs by age and/or sex groups
######################################################################

# ----CONFIG------------------------------------------------------------------------------------------------------------
# clear memory
rm(list=ls())

# runtime configuration
if (Sys.info()["sysname"] == "Linux") {
  package_lib    <- '<<< FILEPATH REDACTED >>>'
  ## Load libraries and  MBG project functions.
  .libPaths(package_lib)
  
  # necessary to set this option in order to read in a non-english character shapefile on a linux system (cluster)
  Sys.setlocale(category = "LC_ALL", locale = "C")
  
} 

#load external packages
pacman::p_load(magrittr, mgsub)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

if (interactive) {
  
  ## Set repo location, indicator group, and some arguments
  repo <-'<<< FILEPATH REDACTED >>>'
  indicator_group <- 'ort'
  indicator <- 'had_diarrhea'
  measure <- 'deaths'
  run_date <- '2019_09_17_14_12_54'
  region <- 'dia_afr_horn'
  modeling_shapefile_version <- '2019_09_10'
  split_by <- 'age_sex'
  
} else {

  # Get MBG arguments
  repo            <- commandArgs()[4]
  indicator_group <- commandArgs()[5]
  indicator       <- commandArgs()[6]
  measure         <- commandArgs()[7]
  run_date        <- commandArgs()[8]
  region          <- commandArgs()[9]
  modeling_shapefile_version <- commandArgs()[10]
  split_by        <- commandArgs()[11]
}

# Set GBD arguments
cause_id        <- ifelse(indicator=='has_lri', 322, 302) #if not LRI, diarrhea
year_id         <- 2000:2019
metric_id       <- c(1,3)
gbd_round_id    <- 6

if (split_by == 'age'){
  age_group_id    <- c(2,3,4,5) 
  sex_id          <- 3
} else if(split_by == 'sex'){
  age_group_id    <- 1 
  sex_id          <- c(1,2)
} else if(split_by == 'age_sex'){
  age_group_id    <- 1:5 
  sex_id          <- 1:3
}

# Load MBG packages and functions
package_list <- c(t(read.csv('<<< FILEPATH REDACTED >>>', header=FALSE)))
source(paste0(repo, 'lbd_core/mbg_central/setup.R'))
mbg_setup(package_list = package_list, paste0(repo, 'lbd_core'))

# Load custom functions
source(paste0(repo, 'lri/y4_age_sex_splitting/functions/format_draw_objects.R'))
source(paste0(repo, 'lri/y4_age_sex_splitting/functions/format_gbd_results.R'))

#use your own diacritics fx, due to inscrutable error
#note: requires mgsub pkg
fix_diacritics <<- function(x) {
  
  require(mgsub)
  
  #first define replacement patterns as a named list
  defs <-
    list('??'='S', '??'='s', '??'='Z', '??'='z', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', '??'='A', 
         '??'='C', '??'='E', '??'='E','??'='E', '??'='E', '??'='I', '??'='I', '??'='I', '??'='I', '??'='N', '??'='O', 
         '??'='O', '??'='O', '??'='O', '??'='O', '??'='O', '??'='U','??'='U', '??'='U', '??'='U', '??'='Y', '??'='B', 
         '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='a', '??'='c','??'='e', '??'='e', '??'='e', 
         '??'='e', '??'='i', '??'='i', '??'='i', '??'='i', '??'='o', '??'='n', '??'='o', '??'='o', '??'='o', '??'='o',
         '??'='o', '??'='o', '??'='u', '??'='u', '??'='u', '??'='y', '??'='y', '??'='b', '??'='y', '??'='Ss')
  
  #then force conversion to UTF-8 and replace with non-diacritic character
  enc2utf8(x) %>% 
    mgsub(., pattern=enc2utf8(names(defs)), replacement = defs) %>% 
    return
  
}

# Define share directory
share_dir <- '<<< FILEPATH REDACTED >>>'

# Create output directory
outdir <- '<<< FILEPATH REDACTED >>>'
dir.create(outdir, showWarnings = FALSE)

# Get GBD measure_id based on measure
measure_id <- get_gbd_measure_id(measure)

# Report arguments
message('\nStarting age and/or sex splitting at the raster level for:')
message('indicator_group = ', indicator_group)
message('indicator = ', indicator)
message('measure = ', measure)
message('run_date = ', run_date)
message('modeling_shapefile_version = ', modeling_shapefile_version)
message('region = ', paste0(region, '\n'))


## Get GBD scaling factors ------------------------------------------------------------------

# Get scaling factors
message('Getting GBD scaling factors')
results <- get_gbd_group_sf(region = region,
                            modeling_shapefile_version = modeling_shapefile_version,
                            year_id = year_id,
                            cause_id = cause_id,
                            age_group_id = age_group_id,
                            sex_id = sex_id,
                            measure_id = measure_id,
                            gbd_round_id = gbd_round_id)


## Pull raster-level outputs ------------------------------------------------------------------

# read in region mean raster brick
message('Reading in raster estimates')
mean_raster <- brick(paste0(share_dir, indicator, '_prediction_', measure,  
                            '_eb_bin0_', region, '_0.grd'))

# read in region upper raster brick
upper_raster <- brick(paste0(share_dir, indicator, '_upper_', measure, 
                             '_eb_bin0_', region, '_0.grd'))

# read in region lower raster brick
lower_raster <- brick(paste0(share_dir, indicator, '_lower_', measure, 
                             '_eb_bin0_', region, '_0.grd'))

# read in GADM codes raster
gadm_raster <- load_admin_raster(admin_level = 0, 
                                 simple_raster = mean_raster,
                                 shapefile_version = modeling_shapefile_version)

# create data table for faster calculations
dt <- data.table(mean = NA, upper = NA, lower = NA, ADM0_CODE = NA, year = NA)
for (i in 1:length(year_id)) {
  dt <- rbindlist(list(dt,
                       data.table(
                         mean = as.vector(mean_raster[[i]]),
                         upper = as.vector(upper_raster[[i]]),
                         lower = as.vector(lower_raster[[i]]),
                         ADM0_CODE = as.vector(gadm_raster[[i]]),
                         year = rep(year_id[[i]], length(as.vector(mean_raster[[i]])))
                       )
  ),
  use.names = T
  )
}

# remove raster values outside of region
dt <- dt[!is.na(ADM0_CODE)]


## Split raster-level outputs ------------------------------------------------------------------

# loop over groups
for (g in unique(results$group)) {
  message('Performing splits for group: ', g)
  
  # subset GBD results to group
  results_by_group <- results[group == g]
  
  # merge GBD results with raster estimates
  dt_group <- copy(dt)
  dt_group <- merge(dt_group, results_by_group, by = c('ADM0_CODE', 'year'), 
                    all.x = T, sort = F, allow.cartesian = T)
  
  # multiply rasters to get group rates
  cols <- c('mean', 'upper', 'lower')
  dt_group[, (cols) := lapply(.SD, function(x, y) {x * y}, y = scaling_factor), 
           .SDcols = cols, by = .I]
  
  # insert split rasters by year
  ras_mean <- ras_upper <- ras_lower <- mean_raster
  for (i in 1:length(year_id)) {
    ras_mean[[i]] <- insertRaster(mean_raster[[i]], 
                                  matrix(dt_group[year == year_id[[i]], mean], ncol = 1))
    ras_upper[[i]] <- insertRaster(mean_raster[[i]],
                                   matrix(dt_group[year == year_id[[i]], upper], ncol = 1))
    ras_lower[[i]] <- insertRaster(mean_raster[[i]], 
                                   matrix(dt_group[year == year_id[[i]], lower], ncol = 1))
  }
  
  # save rasters
  writeRaster(ras_mean, 
              file = paste0(outdir, indicator, '_', measure,
                            '_mean_', g, '_', region), 
              overwrite = TRUE)
  writeRaster(ras_upper, 
              file = paste0(outdir, indicator, '_', measure,
                            '_upper_', g, '_', region), 
              overwrite = TRUE)
  writeRaster(ras_lower, 
              file = paste0(outdir, indicator, '_',  measure,
                            '_lower_', g, '_', region),
              overwrite = TRUE)
  
} # End loop over groups

# done!
message('Finished age and/or sex group splits for raster-level estimates. Saved at\n', outdir)