###############################################################################
###############################################################################
## MBG Launch Script
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- 'USERNAME'
core_repo          <- sprintf("FILEPATH")
indic_repo         <- sprintf("FILEPATH")

## sort some directory stuff
setwd(core_repo)
commondir      <- sprintf("FILEPATH")
package_list <- c(t(read.csv(sprintf("FILEPATH"),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0("FILEPATH"))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0("FILEPATH"))

## Script-specific code begins here ##########################################

# Load from qsub
load_from_parallelize()

# Definitions
indicator_group <- 'vaccine'
sharedir <- sprintf("FILEPATH")

## Read config file (from FILEPATH) and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,      
                      run_tests       = FALSE,     
                      run_date        = run_date)

## Ensure you have defined all necessary settings in your config
check_config()

## Create a few objects from options above
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
gaul_list <- get_adm0_codes(Regions, shapefile_version=modeling_shapefile_version, subnational_raking=subnational_raking)
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
## Launch Parallel Script
###############################################################################

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(as.logical(makeholdouts)) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)

if (use_gn == F) {
  use_c2_nodes <- TRUE
} else {
  use_c2_nodes <- FALSE
}

# If planning to resume broken models, filter loopvars to the correct ones
if (resume_broken == T) {
  message("Checking for completed models and running only broken ones...")
  rerun_all <- data.table(loopvars)
  rerun_all[, rdata_exists := file.exists(paste0(sharedir, "/output/", run_date, "/",
                                                  indicator, "_cell_draws_eb_bin0_", 
                                                  Var1, "_", Var3, ".RData"))]

  rerun_all[, chunk_exists := file.exists(paste0(sharedir, "/output/", run_date, "/", 
                                                 "inla_draws_" ,Var1, "_holdout_", Var3, 
                                                 "_agebin_0_chunk_1.rds"))]
  loopvars <- subset(rerun_all, 
                     rdata_exists == F & chunk_exists == F, 
                     c("Var1", "Var2", "Var3"))
  loopvars <- data.frame(loopvars)
}

## loop over them, save images and submit qsubs
## first check to make sure that there are models to run (for resum_broken)
if (nrow(loopvars) > 0) { 

  for(i in 1:nrow(loopvars)){

    message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))

    # Try to load profiling information from csv

    mem_cores <- load_profiling(vax = vaccine,
                                draws = samples,
                                reg = as.character(loopvars[i, 1]))
    cores <- 4
    # make a qsub string
    qsub <- make_qsub_share(age               = loopvars[i,2],
                            reg               = as.character(loopvars[i,1]),
                            holdout           = loopvars[i,3],
                            test              = test,
                            indic             = indicator,
                            saveimage         = TRUE,
                            corerepo          = core_repo,
                            memory            = as.numeric(mem_cores["mem"]),
                            cores             = as.numeric(mem_cores["slots"]),
                            singularity       = "FILEPATH",
                            singularity_opts  = list(SET_OMP_THREADS = cores, SET_MKL_THREADS = cores),
                            geo_nodes         = TRUE,
                            use_c2_nodes      = FALSE)

    system(qsub)
  
  }
  ## check to make sure models are done before continuing
  waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)
}

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