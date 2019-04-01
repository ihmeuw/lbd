###############################################################################
###############################################################################
## MBG raster plotting script
##
## Purpose: Create quick maps of rasters for model vetting
###############################################################################
###############################################################################

###############################################################################
## SETUP
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

# Load from qsub
# indicator, indicator_group, vaccine, run_date, summstat, raked(T/F), year_list
load_from_parallelize()

## Read config file (from sharedir) and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,           
                      run_date        = run_date)

## Ensure you have defined all necessary settings in your config
check_config()

## Create a few objects from options above
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
gaul_list <- get_gaul_codes(Regions)
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
test <- as.logical(test)

#######################################################################################
## Make some raster maps
#######################################################################################
raked <- ifelse(raked, "raked", "unraked")

make_maps_of_raster(summstats = summstat,
                    rake      = raked,
                    min_value = 0,
                    mid_value = 0.5,
                    max_value = 1,
                    highisbad = F,
                    cores     = as.numeric(slots))

if (summstat == "mean" & indicator == paste0(vaccine, "3_cov")) {
  # Send a pushover map!

  img_file <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", indicator, "/output/", 
                     run_date, "/plots/", indicator, "_", summstat, "_", 
                     raked, "_", max(year_list), ".png")

  pushover_notify(title = paste0("Mapping done for ", indicator, " | ", summstat, " | ", raked),
                  message = paste0("Map of ", indicator, " | ", summstat, " | ", raked, " | ", max(year_list)),
                  image = img_file)
}