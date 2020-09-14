################################################################################
## Script to split aggregated admin-level MBG outputs by age and/or sex groups
################################################################################

## Setup -------------------------------------------------------------------------

# Clear environment
rm(list = ls())

pacman::p_load(magrittr)

#detect if running interactively
interactive <- F  %>% #manual override
  ifelse(., T, !length(commandArgs())>2) %>%  #check length of arguments being passed in
  ifelse(., T, !(is.na(Sys.getenv("RSTUDIO", unset = NA)))) #check if IDE

if (interactive) {

  repo <- '<<< FILEPATH REDACTED >>>'
  indicator       <- 'has_lri'
  indicator_group <- 'lri'
  measure <- 'mortality'
  run_date        <- '2020_06_11_11_19_26'
  split_by        <- 'age_sex' #run both for Y5
  raked           <- T
  modeling_shapefile_version <- '2019_09_10'

} else {
  
  # Get MBG arguments
  repo            <- commandArgs()[4]
  indicator_group <- commandArgs()[5]
  indicator       <- commandArgs()[6]
  raked           <- commandArgs()[7]
  measure         <- commandArgs()[8]
  run_date        <- commandArgs()[9]
  modeling_shapefile_version <- commandArgs()[10]
  split_by        <- commandArgs()[11]

}
# Set GBD arguments
cause_id        <- 322
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
mbg_setup(package_list = package_list, repos = paste0(repo, 'lbd_core'))

# Load custom functions
source(paste0(repo, '/lri/y4_age_sex_splitting/functions/format_draw_objects.R'))
source(paste0(repo, '/lri/y4_age_sex_splitting/functions/format_gbd_results.R'))
source(paste0(repo, '/ort/mbg_central/misc_functions.R'))

# Define share directory
share_dir <- '<<< FILEPATH REDACTED >>>'

# Create output directory
outdir <- '<<< FILEPATH REDACTED >>>'
dir.create(outdir, showWarnings = FALSE)

# Get modeling regions
region <- get_output_regions(share_dir)
if (NA %in% region) region <- region[!is.na(region)]

# Get GBD measure_id based on measure
measure_id <- get_gbd_measure_id(measure)

# Report arguments
message('\nStarting age and/or sex splitting at the admin level for:')
message('indicator_group = ', indicator_group)
message('indicator = ', indicator)
message('raked = ', raked)
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

## Pull draw-level aggregates ------------------------------------------------------------

#first combine
#combine and summarize etiology admin draws
use_inla_country_fes <- FALSE
combine_aggregation(rd       = run_date,
                    indic    = indicator,
                    ig       = indicator_group,
                    ages     = 0,
                    regions  = region,
                    holdouts = 0,
                    raked    = T,
                    measure  = measure,
                    delete_region_files = F,
                    metrics = c('rates','counts'))


# pull estimate aggregates
dt <- format_admin_results(ind_gp = indicator_group,
                           ind = indicator,
                           rd = run_date,
                           measure = paste0('_', measure),
                           suffix = '_eb_bin0_0',
                           rk = raked)

# clean up
dt[, name := NULL]
setnames(dt, indicator, 'b')
dt[, b := as.numeric(b)]

# add admin 1 and 2 codes to GBD results
results <- merge(results, sp_hierarchy_list[, c('ADM0_CODE', 'ADM1_CODE', 'ADM2_CODE')], 
                 by = 'ADM0_CODE', allow.cartesian = T)


## Split draw-level aggregates ------------------------------------------------------------

# loop over admin levels
for (a in 0:2) {
  message('Performing splits for admin level: ', a)
  
  # create a copy and subset
  dt_admin <- copy(dt)
  dt_admin <- dt_admin[agg_level == paste0('ADM', a)]
  dt_admin[, agg_level := NULL]
  
  # create admin column name objects
  ad_code <- paste0('ADM', a, '_CODE')
  all_ad_codes <- paste0('ADM', 0:a, '_CODE')
  all_ad_names <- paste0('ADM', 0:a, '_NAME')
  
  # merge MBG results with with GBD estimates
  setnames(dt_admin, 'code', c(ad_code))
  dt_admin <- merge(dt_admin, 
                    unique(results[, c(ad_code, 'year', 'group', 'scaling_factor'), with = FALSE]), 
                    by = c(ad_code, 'year'), allow.cartesian = T)

  # loop over groups
  for (g in unique(results$group)) {
    message('Saving splits for group: ', g)
    
    # create a copy and subset
    dt_group <- copy(dt_admin)
    dt_group <- dt_group[group == g]
    dt_group[, group := NULL]
    
    # multiply to get split estimate and clean up
    message('~>multiplying by scaling factor')
    dt_group[, split := lapply(.SD, function(x, y) {x * y}, 
                               y = scaling_factor),
             .SDcols = 'b', by = c('year', 'draw', ad_code)]
    dt_group[, c('b', 'scaling_factor') := NULL]
    
    # get mean, upper, and lower
    message('~~>calculating and saving mean, upper, and lower estimates')
    dt_group[, mean := lapply(.SD, mean, na.rm = T), 
             .SDcols = 'split', by = c(ad_code, 'year')]
    dt_group[, upper := lapply(.SD, upper), 
             .SDcols = 'split', by = c(ad_code, 'year')]
    dt_group[, lower := lapply(.SD, lower), 
             .SDcols = 'split', by = c(ad_code, 'year')]
    
    # save summary and free up space
    ad_summary <- merge(unique(sp_hierarchy_list[, c(all_ad_codes, all_ad_names), with = FALSE]), 
                        unique(dt_group[, c(ad_code, 'year', 'mean', 'upper', 'lower'), with = FALSE]), 
                        by = ad_code, allow.cartesian = T)
    write.csv(ad_summary, paste0(outdir, indicator, '_', 
                                 ifelse(raked, measure, 'unraked'),
                                 '_adm', tolower(a), '_summary_', g, '.csv'))
    rm(ad_summary)
    
    # reshape draw object wide
    message('~~~>reshaping wide')
    dt_group <- dcast(dt_group, ... ~ draw, value.var = 'split')
    
    # check for NAs
    message('TESTING: Percent of NA rows per column is: ', mean(is.na(dt_group[, V1])), '%')
    
    # save draw objects
    outpath <- paste0(outdir, indicator, '_', 
                      ifelse(raked, measure, 'unraked'),
                      '_adm', tolower(a), '_draws_', g, '.rds')
    message('-- finished making splits across admin draws. now saving at \n', outpath)
    saveRDS(dt_group, outpath)
    rm(dt_group)
    
  } # End loop over groups
  
  rm(dt_admin)
  
} # End loop over admins

# done!
message('Finished all age and/or sex group splits for admin 0, 1, and 2 draws')