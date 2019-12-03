## Wrapper script to submit mbg_collapse.R by indicator.

library(data.table)

## Set general arguments for this model set.
indicator_group <- 'education'
repo <- "<<<< FILEPATH REDACTED >>>>"
collapse_script <- 'mbg_collapse.R'
launch_slots <- 25

## Log and environment options.
geo_nodes <- FALSE
log_tag <- gsub("-", "_", Sys.Date())
log_dir <- paste0("<<<< FILEPATH REDACTED >>>>", log_tag)
dir.create(log_dir, showWarnings = TRUE)

## make sure to process the data before you run - this reads in all the single year and binned data and saves it as
## "<<<< FILEPATH REDACTED >>>>",all_{by|sy}_{YYYY_MM_DD}.rds 
## to be read in by each collapse job
source(paste0(repo, "/education/process_data.r"))

# if you've added a bunch of new data, recalculate what surveys should be dropped for missingness using drop_for_missingness.R

## Toggle the indicators for which you want to launch models. Right now we have 16 total indicators = (1 mean + 3 props) * (4 age/sex groups).
mean_inds <- c('edu_mean_15_49_female', 'edu_mean_15_49_male', 'edu_mean_20_24_female', 'edu_mean_20_24_male')
prop_inds <- c('edu_zero_prop_15_49_female', 'edu_zero_prop_15_49_male', 'edu_zero_prop_20_24_female', 'edu_zero_prop_20_24_male',
               'edu_no_primary_prop_15_49_female', 'edu_no_primary_prop_15_49_male', 'edu_no_primary_prop_20_24_female', 'edu_no_primary_prop_20_24_male',
               'edu_primary_prop_15_49_female', 'edu_primary_prop_15_49_male', 'edu_primary_prop_20_24_female', 'edu_primary_prop_20_24_male')

all_inds <- c(mean_inds, prop_inds)

#all_inds <- c('edu_zero_prop_15_49_female','edu_mean_15_49_female')

## Loop over indicators to submit collapse script.
node.flag <- ''
shell     <- "shell_prod.sh"  
proj      <- "proj_geospatial" 
if(geo_nodes==TRUE) {
  node.flag <- '-l geos_node=TRUE -v sing_image=default'
  shell     <- "shell_sing.sh"  
  proj      <- "proj_geo_nodes" 
}
for(indicator in all_inds) {
  
  ind_log_dir <- paste0(log_dir, '/', indicator)
  dir.create(ind_log_dir, showWarnings = TRUE)
  qsub <- paste0('qsub -e ', ind_log_dir, ' -o ', ind_log_dir,' -l m_mem_free=200G -q geospatial.q -l h_rt=03:00:00:00 -l fthread=4 -P ', proj, ' -N ', indicator, ' ', '-v sing_image=default /ihme/code/geospatial/nathenry/lbd_core/mbg_central/share_scripts/shell_sing.sh ', repo, '/', indicator_group,
                 '/', collapse_script, ' ', indicator) 
  system(qsub)
  
}

