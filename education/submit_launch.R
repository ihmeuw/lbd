## Wrapper script to submit launch_mean.R by indicator.
## Also includes the option to edit the configs systematically for each model submission,
## this is helpful if you want to submit multiple models per indicator with different configs.

library(data.table)

## Set general arguments for this model set.
indicator_group <- 'education'
repo <- "<<<< FILEPATH REDACTED >>>>"
core_repo <- "<<<< FILEPATH REDACTED >>>>"
launch_script <- 'launch_mean.R'
launch_slots <- 5


## Log and environment options.
geo_nodes <- FALSE
log_tag <- 'stage2_inla_1yearknots'
log_dir <- paste0("<<<< FILEPATH REDACTED >>>>", log_tag)
dir.create(log_dir, showWarnings = TRUE)

## Toggle the indicators for which you want to launch models. Right now we have 16 total indicators = (1 mean + 3 props) * (4 age/sex groups).
mean_inds <- c('edu_mean_15_49_female', 'edu_mean_15_49_male', 'edu_mean_20_24_female', 'edu_mean_20_24_male')
prop_inds <- c('edu_zero_prop_15_49_female', 'edu_zero_prop_15_49_male', 'edu_zero_prop_20_24_female', 'edu_zero_prop_20_24_male',
               'edu_no_primary_prop_15_49_female', 'edu_no_primary_prop_15_49_male', 'edu_no_primary_prop_20_24_female', 'edu_no_primary_prop_20_24_male',
               'edu_primary_prop_15_49_female', 'edu_primary_prop_15_49_male', 'edu_primary_prop_20_24_female', 'edu_primary_prop_20_24_male')

#Groups for launching proportional models through single launch
prop_group_15_49_female <- c('edu_zero_prop_15_49_female', 'edu_no_primary_prop_15_49_female', 'edu_primary_prop_15_49_female')
prop_group_15_49_male <- c('edu_zero_prop_15_49_male', 'edu_no_primary_prop_15_49_male', 'edu_primary_prop_15_49_male')
prop_group_20_24_female <- c('edu_zero_prop_20_24_female', 'edu_no_primary_prop_20_24_female', 'edu_primary_prop_20_24_female')
prop_group_20_424_male <- c('edu_zero_prop_20_24_male', 'edu_no_primary_prop_20_24_male', 'edu_primary_prop_20_24_male')

all_inds <- c(mean_inds, prop_inds)
## Use job names to toggle different config options if you want to submit multiple models per indicator.
job_names <- ''

#########################################################################
########## Submit jobs by indicator #####################################
#########################################################################
## Loop over indicators to submit launch script.
node.flag <- ''
shell     <- "shell_sing.sh"  
proj      <- "proj_geospatial_edu" 
if(geo_nodes==TRUE) {
  node.flag <- '-l geos_node=TRUE'
  shell     <- "shell_sing.sh"  
  proj      <- "proj_geo_nodes_edu" 
}
for(indicator in all_inds) {
  for(job_name in job_names) {
    
    target_sex <- ifelse(grepl('female', indicator), "female", "male")
    age_min <- ifelse(grepl('15_49', indicator), 15, 20)
    age_max <- ifelse(grepl('15_49', indicator), 49, 24)
    if(length(indicator) > 1) indicator <- paste0('prop_group_', age_min, "_", age_max, "_", target_sex)[1] 
    
    # ## Default to basic config file, edit below prior to submission.
    if(grepl('prop', indicator)){
      config_file <- paste0(repo, '/', indicator_group, '/config_props_base.csv')
      config <- fread(config_file, header= FALSE)
    }else {
      config_file <- paste0(repo, '/', indicator_group, '/config_means_base.csv')
      config <- fread(config_file, header= FALSE)
    }
    #check_config()
    if("V3" %in% colnames(config)){
      config$V1 <- NULL
      colnames(config) <- c("V1", "V2")
    } 
    
    ### Make all variants of config, saving as "indicator_[i].csv".

    write.table(config, paste0(repo, '/', indicator_group, '/config_', indicator, '_', job_name, '.csv'), row.names=FALSE, col.names=FALSE, sep=",")
    
    ind_log_dir <- paste0(log_dir, '/', indicator)
    dir.create(ind_log_dir, showWarnings = TRUE)
    dir.create(paste0(ind_log_dir, '/output'), showWarnings = TRUE)
    dir.create(paste0(ind_log_dir, '/errors'), showWarnings = TRUE)
    qsub <- paste0('qsub -e ', ind_log_dir, '/errors -o ', ind_log_dir, '/output -cwd -q all.q -l fthread=2 -l archive=TRUE -l m_mem_free=10G -l h_rt=03:00:00:00', ' ',
                   node.flag, ' -P ', proj, " -v sing_image=default -N ", indicator, ' ', core_repo,  '/mbg_central/share_scripts/', shell, ' ', repo, '/', indicator_group,
                   '/', launch_script, ' config_',indicator,'_',job_name,' ', indicator) 
    system(qsub)
    
  }
}


