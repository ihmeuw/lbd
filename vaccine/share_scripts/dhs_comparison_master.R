###############################################################################
###############################################################################
## DHS launch (master) script
##
## Purpose: Launch DHS comparison code (dhs_comparison_parallel.R) by region
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
indicator <- "dpt3_cov"
run_date <- NULL    # Change this to reflect the model run date
input_date <- NULL  # Change this to reflect the input data date stamp

indicators <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout")
regions <- c("wssa", "cssa", "name", "sssa", "essa")
raked_vals <- c("_raked", "")

dhs_lv <- expand.grid(indicators,
                      regions,
                      raked_vals,
                      stringsAsFactors = F)
dhs_lv <- as.data.table(dhs_lv)
names(dhs_lv) <- c("indicator", "region", "use_raked")
 
dhs_qsub_output <- parallelize(script = "dhs_comparison_parallel",
                               script_dir = paste0(indic_repo, "share_scripts"),
                               log_location = "sgeoutput",
                               lv_table = dhs_lv,
                               save_objs = c('indicator_group', 'run_date', 'input_date'),
                               prefix = "dhs",
                               slots = 15,
                               memory = 20,
                               geo_nodes = TRUE, 
                               singularity = 'default')

monitor_jobs(dhs_qsub_output)

for (ind in indicators) {
  
  # Set up the directory with DHS compare data
  in_dir <- paste0("<<<< FILEPATH REDACTED >>>>/",
                   indicator_group, "/", 
                   ind, "/output/",
                   run_date, "/dhs_compare_data/")
  
  # Do this once each for raked & unraked 
  for (rake in raked_vals) {

    message(ind, " | ", rake)
  
    # Grab all data for raked/unraked
    all_df <- lapply(regions, function(reg) {
      reg_df <- fread(paste0(in_dir, reg, "_", rake, ".csv"))
      reg_df[, region := reg]
      return(reg_df)                      
    })
    all_df <- rbindlist(all_df)

  if (ind == "dpt3_cov") title_ind <- "DPT3 Coverage"
  if (ind == "dpt1_cov") title_ind <- "DPT1 Coverage"
  if (ind == "dpt1_3_abs_dropout") title_ind <- "DPT 1-3 Absolute Dropout"

  plot_title <- paste0(title_ind, " at the Admin 1 Level")
  
  all_df[region == "wssa", region := "Western Sub-Saharan Africa"]
  all_df[region == "cssa", region := "Central Sub-Saharan Africa"]
  all_df[region == "sssa", region := "Southern Sub-Saharan Africa"]
  all_df[region == "essa", region := "Eastern Sub-Saharan Africa"]
  all_df[region == "name", region := "Northern Africa"]

  all_df[,V1 := NULL]

  gg <- ggplot(all_df, aes(x=outcome,y=geo_mean))+
    geom_abline(intercept=0,slope=1,colour='red')+
    geom_point(aes(size = N), colour='black', alpha = 0.2, pch = 16)+
    theme_bw()+
    xlab('Data Estimate (DHS)') +
    ylab('Mean Prediction (MBG)')+
    theme(strip.background = element_rect(fill="white"))+
    geom_errorbar(aes(ymin=geo_lower, ymax=geo_upper), colour="black", width=0, size=.5, alpha = 0.2) +
    geom_abline(intercept=0,slope=1,colour='red') +
    scale_size_area() +
    coord_equal() +
    scale_x_continuous(expand = c(0, 0), limits = c(0,1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,1)) +
    labs(title = plot_title,
         subtitle = paste0("DHS vs ", 
                           ifelse(rake == "", 
                                  "uncalibrated", 
                                  "calibrated"), 
                           " MBG estimates"),
         size = "N",
         color = "Modeling region")

  png(file = paste0(in_dir, "compare_dhs_mbg_", ind, rake, ".png"),
      width = 8, 
      height = 7,
      units = "in",
      res = 300)
  print(gg)
  dev.off()

  # By region
  gg <- gg + facet_wrap(~region)
  png(file = paste0(in_dir, "compare_dhs_mbg_", ind, rake, "_by_region.png"),
      width = 8, 
      height = 7,
      units = "in",
      res = 300)
  print(gg)
  dev.off()
  }
}
