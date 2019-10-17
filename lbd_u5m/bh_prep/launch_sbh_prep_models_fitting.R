## #######################################################
## PURPOSE: makes dated config and launches ./u5m/run_sbh_prep_models.R per region
## #######################################################

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

# get a run date
run_date <- make_time_stamp()

# make directory for the rundate
make_output_dirs(root = root, run_date = run_date)
rd_dir <- sprintf('%s/%s',root,run_date)

#######################################################
## ADD TO NOTES FILE FOR THIS RUN_DATE
write_rd_notes(note = 'Running predictions with censored fix in place plus svy weights in the model fit. ')

#######################################################
## SET OPTIONS AND WRITE CONFIGS
# CONFIG Formula
zzz_formula <- cbind('formula',
                     'died ~ -1 + s(yrborn, sdi, by = factor(ab), k = 9) + s(cdceb100, birthorder, mothage_atbirth, k = 9) + factor(ab) + s(nid, bs = \'re\') + s(country_ab, bs = \'re\')')

# CONFIG Cores
zzz_cores <- cbind('cores', '10')

# CONFIG base year, typically 1970
zzz_by <- cbind('base_year','1970')

# place options into config and write
config <- matrix(NA,ncol=2,nrow=0)
for(o in ls()[grepl('zzz',ls())]) config <- rbind(config, get(o))
write.csv(config,file=sprintf('<<<< FILEPATH REDACTED >>>>',rd_dir))

## Age bins, seperate config
ab_times <- data.frame(
  tstart = c(0, 1,  6, 12, 24, 36, 48),
  tstop  = c(1, 6, 12, 24, 36, 48, 60),
  ab     = c('NN','PNN1','PNN2','1yr','2yr','3yr','4yr') )

ab_times$ab <- paste0(1:length(ab_times$ab),"_",ab_times$ab)
write.csv(ab_times,file=sprintf('<<<< FILEPATH REDACTED >>>>',rd_dir))

#######################################################
## Submit the script for each regional run

# Loop through regions
regions <- c('NAME','SSAWC','ASIA','SSASE','LAC')

for(region in regions){ 
  # make qsub
  qsub <- make_qsub(rtdir       = root,
                    rd          = run_date,
                    cores       = as.numeric(zzz_cores[1,2]),
                    memory      = 150,
                    proj        = 'proj_geospatial',
                    singularity = TRUE,
                    coderepo    = '<<<< FILEPATH REDACTED >>>>',
                    codescript  = sprintf('<<<< FILEPATH REDACTED >>>>', user),
                    geo_nodes   = FALSE,
                    cntry       = region  )
  
  system(qsub)
}

system('qstat')



