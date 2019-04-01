###############################################################################
###############################################################################
## MBG Launch Script
##
## Purpose: Launch MBG models for each modeled indicator
###############################################################################
###############################################################################

###########################################################################
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
load_from_parallelize() #indicator, vaccine, run_date, use_gn

# Definitions
indicator_group <- 'vaccine'
sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>/%s/%s',indicator_group, indicator)

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

## Set up for individual country runs if needed
if (individual_countries == TRUE) {
  # Convert all Regions to individual countries
  Regions <- get_individual_countries(gaul_list)

  # Turn off all FEs
  use_child_country_fes <- F
  use_inla_country_fes <- F
  use_inla_country_res <- F
} 

###############################################################################
## Make Holdouts
###############################################################################
if(as.logical(makeholdouts) & !(file.exists(paste0(sharedir,'/output/',run_date,"/stratum.rds")))){

  message("Recreating holdouts ...")
  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        simple      = NULL,
                        removeyemen = TRUE,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share),
                        yl          = year_list)

  # add in location information
  df <- merge_with_ihme_loc(df)

  # make a list of dfs for each region, with 5 qt folds identified in each

  stratum_ho <- recreate_holdouts(data = df,
                                  row_id_col = "row_id",
                                  load_from_indic = paste0(vaccine, "3_cov"),
                                  rd = run_date,
                                  ig = indicator_group)

  # Save stratum_ho object for reference
  save_stratum_ho(indic       = indicator,
                  ig          = indicator_group,
                  rd          = run_date,
                  stratum_obj = stratum_ho)

}

###############################################################################
## Launch Parallel Script
###############################################################################

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(as.logical(makeholdouts)) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)

## Set up memory requirements
mem <- 100
cores <- 16

## loop over each row (one row per job), save images and submit qsubs
for(i in 1:nrow(loopvars)){

  message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))

  # make a qsub string
  qsub <- make_qsub_share(age           = loopvars[i,2],
                          reg           = as.character(loopvars[i,1]),
                          holdout       = loopvars[i,3],
                          test          = test,
                          indic         = indicator,
                          saveimage     = TRUE,
                          coderepo      = core_repo,
                          memory        = mem,
                          cores         = cores,
                          geo_nodes     = use_gn)

  system(qsub)
}
## check to make sure models are done before continuing
waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

##############################################################################
## Summarize model results
##############################################################################
if (use_gp == T & use_stacking_covs == T) {
  message("Summarizing model results")
  clean_model_results_table()
}

##############################################################################
## Indicate that this model is done
##############################################################################

message(paste0("Done with launch script for ", indicator))

## END OF FILE
###############################################################################