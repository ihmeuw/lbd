###############################################################################
###############################################################################
## Correlations before and after calibration ("raking")
##
## Purpose: Create scatter plots of pre- and post-calibration estimates
###############################################################################
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific code begins here ##########################################

indicator_group <- 'vaccine'
indicator <- "dpt3_cov"
run_date <- NULL # Change this to reflect the model run date
input_date <- NULL # Change this to reflect the date stamp for the input data  

indicators <- c("dpt3_cov")
regions <- c("wssa", "cssa", "name", "sssa", "essa")
raked_vals <- c("_raked", "")

corr_table <- data.table(expand.grid(indicators, c(T,F), stringsAsFactors = F))
names(corr_table) <- c("indicator", "raked")
corr_table[, pearson := as.numeric()]
corr_table[, spearman := as.numeric()]

for (i in 1:nrow(corr_table)) {
  ind <- corr_table[i, indicator]
  rr <- corr_table[i, raked]
  rake <- ifelse(rr == T, "_raked", "")
  
  message(ind, " | ", rake)
  
  in_dir <- paste0("<<<< FILEPATH REDACTED >>>>/",
                   indicator_group, "/", 
                   ind, "/output/",
                   run_date, "/dhs_compare_data/")
  
  # Grab all data for raked/unraked
  all_df <- lapply(regions, function(reg) {
    reg_df <- fread(paste0(in_dir, reg, "_", rake, ".csv"))
    reg_df[, region := reg]
    return(reg_df)                      
  })
  all_df <- rbindlist(all_df)
  
  corr_table[i, ]$pearson <- cor(all_df$outcome, all_df$geo_mean, method = "pearson")
  corr_table[i, ]$spearman <- cor(all_df$outcome, all_df$geo_mean, method = "spearman")
  
}