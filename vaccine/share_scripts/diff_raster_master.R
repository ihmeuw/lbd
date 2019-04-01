###############################################################################
###############################################################################
## Generate difference rasters - Master Script
##
## Purpose: Launch parallel jobs to calculate draw-level difference between
##			first and last modeling year and then generate rasters of 
##			mean as well as upper and lower CIs. Launches 
#			`diff_raster_parallel.R` by region.
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

### OPTIONS DEFINED HERE ################################################################

## Options
indicator_group <- "vaccine"
run_date <- NULL # Model run date
year_list <- c(2000:2016)
summstats <- c("mean", "upper", "lower")

indicators <- c("dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout", "dpt1_3_rel_dropout")
regions <- c("wssa", "cssa", "name", "sssa", "essa")
raked_vals <- c(TRUE)


diff_lv <- expand.grid(indicators,
                       regions,
                       raked_vals,
                       stringsAsFactors = F)
diff_lv <- as.data.table(diff_lv)
names(diff_lv) <- c("indicator", "region", "raked")
 
diff_qsub_output <- parallelize(script = "diff_raster_parallel",
                                script_dir = paste0(indic_repo, "share_scripts"),
                                log_location = "sgeoutput",
                                lv_table = diff_lv,
                                save_objs = c('indicator_group', 'run_date', 'year_list', 'summstats'),
                                prefix = "diff",
                                slots = 15,
                                memory = 20,
                                geo_nodes = FALSE, 
                                c2_nodes = TRUE,
                                singularity = 'default')

monitor_jobs(diff_qsub_output)

# Combine everything together
message("\nCombining rasters...")
ig <- indicator_group

for (indic in indicators) {
	sharedir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", indic, "/output/", run_date, "/")

	for (rake in raked_vals) {
		rake_addin <- ifelse(rake, "_raked", "")
		message(paste0("indic: ", indic))
		
  		for(ss in summstats){
	      message(paste0('  ',ss))
	      rlist <- lapply(regions, function(reg) {
	      	message(paste0('    ',reg))
	        ras <- brick(paste0(sharedir, indic, "_", reg, rake_addin, "_diff_", 
	                            min(year_list), "-", max(year_list), "_", ss, "_raster.tif"))
	        return(ras)
	      })	

	      if(ss=='cirange') ssname = 'range' else ssname = ss # naming convention
	      rbrick <- do.call(raster::merge, unname(rlist))
	      out_file <- paste0(sharedir, indic, rake_addin, "_diff_", 
	                         min(year_list), "-", max(year_list), "_", ss, "_raster.tif")
		    writeRaster(rbrick,
			            file = out_file,
			            format='GTiff',
			            overwrite = TRUE)			
		}
	}
}

message("Done with master script")