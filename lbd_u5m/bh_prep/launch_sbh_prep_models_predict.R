## #######################################################
## PURPOSE: launches run_sbh_prep_models_predict.R per NID
## #######################################################

####################################################################
# USER INPUT #######################################################

existing_run_date <- NULL

####################################################################

# Assumes you have your repo cloned 
user      <- Sys.info()['user']
gitbranch <- 'master'

# git pull, if you want
system(sprintf('<<<< FILEPATH REDACTED >>>>', user, gitbranch))

# root dir for outputs
root <- '<<<< FILEPATH REDACTED >>>>'

# load custom functions
setwd(sprintf('<<<< FILEPATH REDACTED >>>>', user))
source('<<<< FILEPATH REDACTED >>>>')

core_repo <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c("data.table", "splitstackshape")
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '<<<< FILEPATH REDACTED >>>>'))
mbg_setup(package_list = package_list, repos = core_repo)

# get a run date
if(is.null(existing_run_date)){
  run_date <- make_time_stamp()
  make_output_dirs(root = root, run_date = run_date)
  message(sprintf('New run date: %s',run_date))
} else {
  run_date <- existing_run_date
  message(sprintf('Using existing run date: %s',run_date))
}

# make directory for the rundate
rd_dir <- sprintf('%s/%s/',root,run_date)

# grab all nids of SBH surveys in our database
nids <- fread('<<<< FILEPATH REDACTED >>>>')
nids <- subset(nids, type == 'SBH') # keep only those that are SBH
nids <- subset(nids, year >= 1998 ) # keep only those after 1998

#######################################################
# QSUB all NIDS

# Set the below to TRUE  if you want estimates at the survey level (use this for SBH-CBH paper)
# Set the below to FALSE if you want estimates at the smallest available geography level
aggregate_to_nid <- FALSE

for(n in nn){ 
  
  # make qsub
  qsub <- make_qsub(rtdir       = root,
                    rd          = run_date,
                    cores       = 10, 
                    memory      = 150,
                    proj        = 'proj_geospatial',
                    coderepo    = '<<<< FILEPATH REDACTED >>>>',
                    codescript  = sprintf('<<<< FILEPATH REDACTED >>>>', user),
                    geo_nodes   = TRUE,
                    singularity = TRUE, # set to 10 cores right now
                    cntry       = n,
                    adl_args    = aggregate_to_nid)

  system(qsub)
}

system('qstat')

### NOTE: ONCE THESE MODELS HAVE FINISHED, YOU NEXT RUN SBH_COLLAPSE.R TO COMBINE THEM
