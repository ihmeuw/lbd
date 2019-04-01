###############################################################################
###############################################################################
## MBG Merge Regions Script
##
## Purpose: Merge together regional outputs to create unified tifs and csvs
##          for all of Africa together, for easier subsequent analysis
###############################################################################
###############################################################################
####
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

## Convenience
sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>/%s/%s',indicator_group, indicator)
strata <- Regions

#######################################################################################
## Set up options
#######################################################################################

if (raked == T) {
  rr <- c("raked")
  rf_tab <- T
} else if (raked == F) {
  rr <- c("unraked")
  rf_tab <- F
}

if (modeled == T) {
  run_summary <- T
} else if (modeled == F) {
  run_summary <- F
}

if (("p_below" %in% summstats) & !(indicator %in% c(paste0(vaccine, "1_cov"), paste0(vaccine, "3_cov")))) {
  # Only run this for 1st and 3rd doses
  summstats <- summstats[summstats != "p_below"]
}

#######################################################################################
## Merge regions and save
#######################################################################################

post_load_combine_save(regions    = Regions,
                       summstats  = summstats,
                       raked      = rr,
                       rf_table   = rf_tab,
                       run_summ   = run_summary,
                       indic      = indicator,
                       ig         = indicator_group,
                       sdir       = sharedir)

# Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)