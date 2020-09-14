##################################################
###           LRI COLLAPSE SCRIPT             ####
#####      Written By: Mathew Baumann        #####
##################################################
rm(list=ls())

##############################################
# General setup
##############################################
cores <- 5
repo <- '<<<< FILEPATH REDACTED >>>>/02_collapse' 
setwd(repo)

module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)
folder_in <- '<<<< FILEPATH REDACTED >>>>'
folder_out <- '<<<< FILEPATH REDACTED >>>>'

source('collapse_functions/general_functions.R')
source("<<<< FILEPATH REDACTED >>>>/mbg_central/setup.R")
package_list <- c('mgcv', 'nlme', 'lme4', 'survey', 'foreign', 'rgeos', 'data.table','raster',
                  'rgdal','seegSDM','seegMBG','plyr','dplyr', 'foreach', 'doParallel', 'feather','boot')
mbg_setup(package_list = package_list, repos='<<<< FILEPATH REDACTED >>>>')

##############################################
# Load and Clean Data
##############################################
latest_postextraction <- get_latest_file(folder_in, '*.feather')

#load most recent geomatched data
lri_data <- data.table(read_feather(paste0(folder_in, latest_postextraction)))

#Read in supplementary information
locs <- fread('<<<< FILEPATH REDACTED >>>>') #location metadata
locs <- locs[,c('location_name','location_id','super_region_name','super_region_id',
                'region_name','region_id','parent_id', 'ihme_loc_id')]
setnames(locs, old = c('ihme_loc_id'), new = c('country'))
scalar <- fread('<<<< FILEPATH REDACTED >>>>') #GBD seasonality scalars
mm_scalar <- fread('<<<< FILEPATH REDACTED >>>>') #nids missing monthly information for seasonality
dm.coeffs <- fread('<<<< FILEPATH REDACTED >>>>') #DisMod crosswalk coefficients
duration <- fread('<<<< FILEPATH REDACTED >>>>') #LRI duration draws

#Initial cleaning: Fixing/renaming columns, subsetting ages, prep case definitions
lri_data <- initial_cleaning(lri_data)
#Assign case definitions for surveys that include and do not include chest symptoms and fever
lri_data <- mclapply(unique(lri_data$nid), assign_case_defs, mc.cores=cores)
lri_data <- rbindlist(lri_data)


##############################################
# Seasonality and Definition Crosswalks
##############################################
lri_data <- lri_data[,-'location_name']
source('collapse_functions/seasonality_missing_month_prep.R')
lri_data <- apply_crosswalks(lri_data)

##############################################
# Collapse
##############################################
lri_data <- collapse(lri_data)

report_data <- fread('<<<< FILEPATH REDACTED >>>>') #tabulated data not run through normal pipeline until this point
report_data <- report_data[,-'optional']
report_data$cluster_id <- c((nrow(lri_data) + 1): (nrow(lri_data) + nrow(report_data)))
lri_data <- rbind(lri_data, report_data)

##############################################
# Age Crosswalks
##############################################
coverage_data <- copy(lri_data)
names(coverage_data)[names(coverage_data) == 'nid'] <- 'svy_id'
#save unadjusted prevalence
coverage_data[,unadj := has_lri]
#age crosswalk
source('collapse_functions/crosswalk/lri_age_crosswalk_insert_code.R')
coverage_data[,age_adj := has_lri]
names(coverage_data)[names(coverage_data) == 'svy_id'] <- 'nid'

write.csv(coverage_data, '<<<< FILEPATH REDACTED >>>>',row.names = F)

##############################################
# Resampling
##############################################
#Resample in parallel by shapefile, outputs csvs of each shapefile
#Code will stall until qsubs finish
source('collapse_functions/resample/parent.R')

##############################################
# Final Cleanup
##############################################
#Append resampled data, adjust for recall period
all_collapsed <- final_cleanup()
write.csv(all_collapsed, '<<<< FILEPATH REDACTED >>>>',row.names = F)
nids_diag[nid %in% all_collapsed$nid, reason_dropped := 'NOT DROPPED']
write.csv(nids_diag, '<<<< FILEPATH REDACTED >>>>', row.names = F)
