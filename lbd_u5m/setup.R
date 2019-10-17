
# save some other important drive locations
root           <- '<<<< FILEPATH REDACTED >>>>'
datadrive      <- '<<<< FILEPATH REDACTED >>>>'
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'

#unlink(sprintf('%s/model_image_history/',sharedir))

# check if in singularity image, if not then set up package locations
if(!'SINGULARITY_NAME' %in% names(Sys.getenv())){
  package_lib    <- ifelse(grepl("geos", Sys.info()[4]),
                           paste0(root,'<<<< FILEPATH REDACTED >>>>'),
                           paste0(root,'<<<< FILEPATH REDACTED >>>>'))
  .libPaths(package_lib)
} 

## Load libraries and MBG project functions.
package_list<-c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
# 
# if(Sys.info()[1] == "Windows"){
#   message("STOP! you will overwrite these packages if you run from windows")
#   message("STOP! also, lots of this functions won't work so get on the cluster!")
# } else {
#   for(package in package_list)
#     require(package, character.only=TRUE)
#   for(funk in list.files(path = c(ig_repo,core_repo), recursive=TRUE, pattern='functions', full.names=TRUE)){
#    message(funk)
#     source(funk) 
#   }
# }
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)
load_mbg_functions(ig_repo)


## Read config file and save all parameters in memory
config <- set_up_config(
  repo=ig_repo,
  core_repo=core_repo,
  indicator_group='u5m',
  indicator='died',
  config_file=paste0(ig_repo,config_name,'.csv'),
  covs_file=paste0(ig_repo,covs_config_name,'.csv'),
  post_est_only=FALSE,
  run_date=run_date,
  push_to_global_env=TRUE,
  run_tests=FALSE
)

slots      <- as.numeric(slots)
inla_cores <- as.numeric(inla_cores)


## Create directory structure for this model run, including an output for each agebin
create_dirs(indicator_group = indicator_group,
            indicator       = indicator)

dir.create(sprintf('%s/output/%s', sharedir, run_date))

for(a in 1:numagebins){
   dir.create(sprintf('/%s_age%i/', sharedir, a))
   dir.create(sprintf('/%s_age%i/output', sharedir, a))
   dir.create(sprintf('/%s_age%i/model_image_history', sharedir, a))
   dir.create(sprintf('/%s_age%i/output/%s', sharedir, a, run_date))
}

for(outputage in c('neonatal','infant','under5'))
  dir.create(sprintf('/%s_%s/output/%s', sharedir, outputage,run_date))


