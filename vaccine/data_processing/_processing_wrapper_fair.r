#----HEADER------------------------------------------------------------------------------------------------------------
# Purpose: Process extracted unit record data from UbCov
#          Step 1: Launch processing code; saves dataset by NID
#               Creates indicators of interest
#               Apply restrictions (date formatting, age restrictions, vaccine definitions)
#               Explore card versus recall
#               Address missingness
#          Step 2: Launch tabulation code; saves dataset by NID
#               LBD - Collapse by lat/long or admin unit
#               GBD - Collapse by GBD location ids
# Inputs:  nid  ->  NID on which to launch processing
#                   if nid=="all", parent script will launch all NIDs that have not already been processed
#***********************************************************************************************************************


#----SETUP-------------------------------------------------------------------------------------------------------------
username        <- 'USERNAME'
core_repo       <- paste0("FILEPATH")
vaccines_repo   <- paste0("FILEPATH")
extraction_root <- "FILEPATH"

# Load packages available to the LBG Singularity image
mbg_setup(c("proto", "findpython", "getopt", "argparse", "data.table", "magrittr", "survey", "parallel", "plyr", "dplyr", "haven",
            "rgeos", "raster", "rgdal", "dismo", "gbm", "foreign", "doParallel", "grid", "gridExtra", "gtools", "ggplot2",
            "sf", "assertthat", "INLA", "seegSDM", "seegMBG", "pacman", "glmnet", "RMySQL", "tictoc"), repos=core_repo)

# Load packages directly from file path for packages not yet available to LBD Singularity image
library("fasterize", lib.loc="FILEPATH", character.only=TRUE)
library("binom", lib.loc="FILEPATH", character.only=TRUE)


# Introduced to record data dropped at each stage of data filtering - used for publications table
record_data_drop <- FALSE
if(record_data_drop) {
  filter_table_path <- "FILEPATH"
  vaccine_to_record <- "mcv"
}

#----SETUP--------------------------------------------------------------------------------------------------------------

#**********************************************************************************************************************


#----LOAD ARGUMENTS----------------------------------------------------------------------------------------------------
### setup arguments
parser <- ArgumentParser()
parser$add_argument("--action", help="'launch' to launch a job to process each NID listed; 'process' to actually run processing", default="launch", type="character")
parser$add_argument("--nid", help="NID(s) to process", type="character")
parser$add_argument("--project", help="prod cluster project", default="proj_covariates", type="character")
parser$add_argument("--slots", help="how many slots?", default=5, type="integer")
parser$add_argument("--date", help="leave missing", default=as.character(Sys.Date()), type="character")
parser$add_argument("--lbd_only", help="logical; only run LBD tabulation (i.e. not GBD tabulation or processing)?", default="FALSE", type="character")
parser$add_argument("--gbd_only", help="logical; run everything but LBD tabulation (i.e. processing and GBD tabulation)?", default="FALSE", type="character")


### save job args as objects
args <- parser$parse_args()
message(args)
list2env(args, environment()); rm(args)
if (nid=="all_missing") {
  extracted <- gsub(".csv", "", gsub(".*_\\s*|_.*", "", list.files(file.path(extraction_root, "raw")))) %>% unique %>% as.numeric
  processed <- gsub(".csv", "", list.files(file.path(extraction_root, "processed"))) %>% unique %>% as.numeric
  nids      <- extracted[!extracted %in% processed] %>% sort
} else if (nid=="all") {
  nids      <- gsub(".csv", "", gsub(".*_\\s*|_.*", "", list.files(file.path(extraction_root, "raw")))) %>% unique %>% as.numeric %>% sort
} else { 
  nids      <- as.numeric(strsplit(gsub(",", " ", nid), split=" +")[[1]])
}
lbd_only <- as.logical(lbd_only)
gbd_only <- as.logical(gbd_only)

### source function
source(paste0(vaccines_repo, "/process.r"))

#**********************************************************************************************************************


#----LAUNCH------------------------------------------------------------------------------------------------------------
if (action=="launch") {
  print("step 1")
  ### launch one job per NID
  for (NID in nids) {
    
    # set date
    date <- format(lubridate::with_tz(Sys.time(), tzone="America/Los_Angeles"), "%Y-%m-%d")
    
    # launch
    job_name <- paste0("CHILD_", NID, "_vax_processing")
    shell <- "/FILEPATH"
    job <- paste0("qsub -N ", job_name,
                  " -l m_mem_free=8G ", "-P ", project,
                  " -q geospatial.q",
                  " -l fthread=8",
                  " -l archive=TRUE",
                  " -o FILEPATH",
                  " ", shell, " -e s ",
                  vaccines_repo, "FILEPATH/_processing_wrapper_fair.r",
                  " --action process --nid ", NID, " --project ", project, " --date ", date, " --gbd_only ", gbd_only, " --lbd_only ", lbd_only)
    system(job); print(job)
  } 
  ### check for success/failure
  source("FILEPATH/ubcov_tools.r")
  # wait

  message("Beginning job_hold()")
  job_hold(paste0("CHILD_", nids, "_vax_processing"))
  # check each job
  message("Beginning lapply()")
  invisible(lapply(nids, check_nids))
  
} else if (action=="process") {
  
  ### Step 1: process that NID
  
  # launch
  if (!lbd_only) success_process <- process(nid) else success_process <- TRUE
  
  ### Step 2: launch tabulation
  message("Launch tabulation")
  print("Step 2")
  if (success_process) {
    
    # GBD tabulation
    if (!lbd_only) source(paste0(vaccines_repo, "FILEPATH/tabulate_gbd.r"))
    if (!lbd_only) tabulate_gbd(nid)
    
    # LDB tabulation
    if (!gbd_only) {
      
      # tabulate
      source(paste0(vaccines_repo, "FILEPATH/tabulate_lbd.r"))
      tabulate_lbd(nid)
      
      # resample
      source(paste0(vaccines_repo, "FILEPATH/resample.r"))
      resample_geo_polys(nid, extraction_root=extraction_root, prefixes=c(unique(gsub("[0-9]", "", vaccines))[!unique(gsub("[0-9]", "", vaccines)) %in% c("pent", "tetra", "mmr")], "hib3_dpt3_ratio",
      "hepb3_dpt3_ratio", "dpt3_under_1", "mcv2_mcv1_ratio", "pcv3_dpt3_ratio", "rotac_dpt3_ratio"))
      
    }
  }
}

#**********************************************************************************************************************
