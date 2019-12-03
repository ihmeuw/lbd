###########################################################################################
###########################################################################################
## Launch MBG model for education
###########################################################################################
###########################################################################################

rm(list=ls())

## Set repo locations and indicator group
user <- Sys.info()['user']
message(user)
message(sessionInfo())
core_repo <- paste0("<<<< FILEPATH REDACTED >>>>")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")
indicator_group <- 'education'
rerun_run_date <- NA

config_arg <- as.character(commandArgs()[4]); message(config_arg)
indicator <- as.character(commandArgs()[5]); message(indicator)

if(grepl("prop_group_", indicator)){
  target_sex <- ifelse(grepl('female', indicator), "female", "male")
  age_min <- ifelse(grepl('15_49', indicator), 15, 20)
  age_max <- ifelse(grepl('15_49', indicator), 49, 24)
  i <- c("edu_zero_prop", "edu_primary_prop", "edu_no_primary_prop")
  indicator <- paste0(i, "_", age_min, "_", age_max, "_", target_sex)
}

root           <- ifelse(Sys.info()[1]=='Windows', 'J:/', '/home/j/')
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                         sprintf("<<<< FILEPATH REDACTED >>>>"),
                         sprintf("<<<< FILEPATH REDACTED >>>>"))

commondir      <- sprintf("<<<< FILEPATH REDACTED >>>>")

## Load libraries and  MBG project functions.
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)#, load_from_sing = FALSE)
library(assertthat)
## Load libraries and miscellaneous MBG project functions.

## Read config file and save all parameters in memory
config <- load_config(repo = indic_repo,
                      core_repo = core_repo,
                      indicator_group = indicator_group,
                      indicator = indicator,
                      config_name = config_arg,
                      covs_name = 'covs_edu_mean',
                      run_tests = FALSE)


Regions <- strsplit(region_list," ")
Regions <- Regions[[1]][Regions[[1]] != "+"]
message(fixed_effects)
year_list <- eval(parse(text=year_list))

## Create run date in correct format
run_date <- make_time_stamp(time_stamp)
if(!is.na(rerun_run_date)) {
  ## Replace run_date with existing model date and wipe any "fin" files so waitformodelstofinish works on the reruns
  run_date <- rerun_run_date
  fin_files <- list.files("<<<< FILEPATH REDACTED >>>>", pattern = '^fin_', full.names = TRUE)
  for(f in fin_files) unlink(f)
}
run_date <- paste0("<<<< DATE REDACTED >>>>", config_arg)
## Create directory structure for this model run
for(i in indicator){
  create_dirs(indicator_group = indicator_group,
              indicator = i)
}

## Load full gaul list to model over
gaul_list <- get_adm0_codes(unlist(Regions), shapefile_version = modeling_shapefile_version)

new_folds <- F
## Run function to create holdouts (returns list of data.frames with an additional "folds" column)
if(crossval==TRUE & new_folds == T) {
  
  ## Load simple polygon template to model over
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 0.4, shapefile_version = modeling_shapefile_version)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  ## Load input data and make sure it is cropped to modeling area
  df <- load_input_data(indicator   = indicator,
                        simple      = simple_polygon,
                        agebin      = 0,
                        removeyemen = FALSE,
                        years       = yearload,
                        use_share   = TRUE,
                        yl          = year_list,
                        update_run_date = FALSE)
  
  ## Add GAUL_CODE and region to df given the lat/longs. We need this to stratify in holdout function.
  df <- add_gauls_regions(df = df,
                          simple_raster = simple_raster)
    table(df$region, df$year)
  long_col = 'longitude'
  lat_col = 'latitude'
  n_folds = as.numeric(n_ho_folds)
  stratum_ho <- make_folds(data = df, n_folds = n_folds, spte_strat = 'nids', strat_cols = "region",
                           ts = 20, mb = 10)
} 
if(crossval == T & new_folds == F){
  ## Load simple polygon template to model over
  load("<<<< FILEPATH REDACTED >>>>")
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  load("<<<< FILEPATH REDACTED >>>>")
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  ## Load input data and make sure it is cropped to modeling area
  df <- load_input_data(indicator   = indicator,
                        simple      = simple_polygon,
                        agebin      = 0,
                        removeyemen = FALSE,
                        years       = yearload,
                        use_share   = TRUE,
                        yl          = year_list,
                        update_run_date = FALSE)
  
  ## Add GAUL_CODE and region to df given the lat/longs. We need this to stratify in holdout function.
  df <- add_gauls_regions(df = df,
                          simple_raster = simple_raster)
  
  from_date <- "<<<< DATE REDACTED >>>>"
  stratum_ho <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
  saveRDS(stratum_ho, file = paste0("<<<< FILEPATH REDACTED >>>>"))
}

## Set strata as character vector of each strata (in my case, just stratifying by region whereas U5M stratifies by region/age)
strata <- unlist(Regions)

## ~~~~~~~~~~~~~~~~~~~~~~~~  Parallel MBG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~ Submit job by strata/holdout  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for(r in strata){
  
  fin <- 'FALSE'
  if(!is.na(rerun_run_date)) {
    # Only qsub if that strata has any failed holdout jobs ("fin" file doesn't exist)
    fins <- 1
    for(holdout in c(0:n_folds)) {
      holdout_fin <- file.exists(paste0("<<<< FILEPATH REDACTED >>>>"))
      fins <- c(fins, holdout_fin)
    }
    if(sum(fins)==length(fins)) {
      # If this strata was completed finished, mark to not rerun
      fin <- 'TRUE'
    }
    # If we are rerunning this strata because one or more holdouts broke, delete fin files for all holdouts.
    # This is because if one broke, we have to relaunch all because we have to reorganize holdout ids.
    if(fin==FALSE) {
      delete_fin_files <- list.files("<<<< FILEPATH REDACTED >>>>",
                                     pattern = paste0('fin__bin0_', r, '_'),
                                     full.names = TRUE)
      lapply(delete_fin_files, unlink)
    }
  }
  
}


ram_per_draw_list <- c("noaf" = 150/500,
                       "cssa" = 90/500, 
                       "essa" = 140/500,
                       "wssa" = 150/500,
                       "sssa" = 65/500,
                       "stan" = 90/500,
                       "mide" = 50/500,
                       "seas" = 50/500,
                       "soas" = 70/500,
                       "eaas" = 230/500,
                       "ocea" = 70/500,
                       "caca" = 70/500,
                       "ansa+pry+sur+guy+guf" = 120/500,
                       "bra+ury" = 160/500,
                       "bra" = 160/500,
                       'uga' = 70/500)

samples <- eval(parse(text=samples))
ram_list <- round(ram_per_draw_list * samples)

if(!is_new_cluster()){
  if(use_geos_nodes){
    proj <- 'proj_geo_nodes_edu'
    slot_list <- ram_list/18
  } else{
    proj <- 'proj_geospatial_edu'
    slot_list <- ram_list/9
    use_c2_nodes <- TRUE
  }
}

if(crossval) loopvars <- expand.grid(strata, 0, 0:n_ho_folds) else loopvars <- expand.grid(strata, 0, 0)

for(indic in indicator){
  ## loop over them, save images and submit qsubs
  if(fin == FALSE) {
    
    for(i in 1:nrow(loopvars)){
      region <- as.character(loopvars[i,1])
      sharedir <- "<<<< FILEPATH REDACTED >>>>"
      message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))
      
      # make a qsub string
      qsub <- make_qsub_share(age           = loopvars[i,2],
                              reg           = as.character(loopvars[i,1]),
                              holdout       = loopvars[i,3],
                              test          = FALSE,
                              memory        = ram_list[[region]],
                              queue         = "long.q",
                              run_time      = "01:00:00:00",
                              indic         = indic,
                              saveimage     = TRUE,
                              cores         = ifelse(is_new_cluster(), 4, slot_list[[region]]),
                              corerepo = core_repo,
                              singularity =  "<<<< FILEPATH REDACTED >>>>",
                              proj = 'proj_geospatial',
                              singularity_opts = list(SET_OMP_THREADS=4, SET_MKL_THREADS = 1))
      system(qsub)
      
      
    }
  }
}
ifelse(grepl('prop', indicator), tag <- 'prop', tag <- 'mean')

waitformodelstofinish(sleeptime=60, lv=loopvars)

####################################################
####                  RAKING                    ####
####################################################

## Parallelized post-estimation over region
print(indicator)
summstats <- c('mean', 'upper', 'lower')
postest_script <- "postest_script"
setwd(core_repo)
raked <- T
year_list <- c(2000:2017)
rake_transform <- 'logit'

for(i in indicator){
  run_date <- paste0("<<<< DATE REDACTED >>>>")
  if(raked == T & !grepl('mean', i)){
    gbd <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  } else{
    gbd <- NULL
  }
  ifelse(grepl('mean', i), rake_transform <- 'linear', rake_transform <- 'logit')
  
  prep_postest(indicator = i,
               indicator_group = indicator_group,
               run_date = run_date,
               save_objs = c("core_repo", "year_list", "summstats", "gbd",
                             "rake_transform", "pop_measure", "pop_release"))
  for (s in strata) {
    
    qsub <- make_qsub_postest(code = postest_script,
                              stratum = s,
                              log_location = 'sharedir',
                              memory        = ifelse(s %in% c('bra', 'eaas', 'noaf', 'wssa'),300,150),
                              run_time = "03:00:00:00",
                              queue = 'geospatial.q',
                              cores         = 2,
                              geo_nodes = F,
                              indic = i,
                              singularity = 'default',
                              proj = 'proj_geo_nodes',
                              raking_shapefile_version = raking_shapefile_version,
                              modeling_shapefile_version = modeling_shapefile_version,
                              subnat_raking = subnational_raking)
    qsub <- paste0(qsub, ' ', s)
    system(qsub)
  }
}

waitforpostesttofinish(strata = strata, sleeptime=100)

#####################################################
####                  RESCALE                    ####
#####################################################
if(!grepl('mean', indicator)){
  for(age_sex in c('15_49_male','20_24_male',  '15_49_female','20_24_female')){
    for(r in strata) {
      run_date <- paste0('st_test')
      message(paste0('Submitting rescale ', r))
      use_raked <- 1
      sharedir <- paste0("<<<< FILEPATH REDACTED >>>>")
      rescale_qsub <- paste0('qsub -e ', sharedir,  '/errors -o ', sharedir, '/output -cwd -l m_mem_free=500G -P proj_geo_nodes_edu -q geospatial.q -l h_rt=00:24:00:00 -l fthread=5 -v sing_image=default -N rescale_',
                             r, '_', age_sex, ' ', core_repo, '/mbg_central/share_scripts/shell_sing.sh ', indic_repo, '/education/parallel_post_est_props.R ', core_repo, " ", indic_repo, " ", age_sex, " ", run_date, " ", use_raked, " ", r) 
      system(rescale_qsub)
    }
  }
  
  waitforrescaletofinish(strata = strata)
  
  
  i <- c( "edu_primary_prop", "edu_secondary_prop", "edu_zero_no_primary_prop")
  indicator <- paste0(i, "_", age_sex)
  for(i in indicator){
    if(!grepl('zero_prop', i)) sharedir <- paste0("<<<< FILEPATH REDACTED >>>>")
    if(grepl('zero', i)) sharedir <- paste0("<<<< FILEPATH REDACTED >>>>")
    if(grepl('secondary|zero_no_primary', i)){
      post_load_combine_save(indic = i, summstats = c("mean", "lower", "upper"), raked = "raked", rf_table = FALSE, run_summ = F, regions = strata, sdir = sharedir)
    }else {
      post_load_combine_save(indic = i, summstats = c("mean", "lower", "upper"), raked = "raked", rf_table = TRUE, run_summ = F, regions = strata, sdir = sharedir)
    }
  }
} else{
  for(i in indicator){
    run_date <- paste0("<<<< DATE REDACTED >>>>")
    sharedir <- paste0("<<<< FILEPATH REDACTED >>>>")
    post_load_combine_save(regions = strata, indic = i, summstats = c('mean', 'lower', 'upper'), raked = c('raked'), sdir = sharedir, rf_table = T, run_summ = F)
  }
}


######################################################
####                  AGGREGATE                   ####
######################################################

if(!grepl('mean', indicator)){
  dir.create("<<<< FILEPATH REDACTED >>>>")
  file.copy("<<<< FILEPATH REDACTED >>>>", 
            "<<<< FILEPATH REDACTED >>>>")
}

if(!grepl('mean', indicator)) rak <- T
if(grepl('mean', indicator)) rak <- c(T,F)
indicator <- c(indicator, paste0('edu_secondary_prop_', c('15_49_female', '15_49_male', '20_24_female', '20_24_male')),
               paste0('edu_zero_no_primary_prop_', c('15_49_female',  '20_24_female', '20_24_male','15_49_male')))
## Launch jobs to aggregate to admin units
for(i in indicator){
  for(s in strata){
    run_date <- paste0()
    log_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
    dir.create(log_dir)
    submit_aggregation_script(indicator=i,
                              indicator_group=indicator_group,
                              run_date=run_date,
                              raked=rak,
                              pop_measure=pop_measure,
                              run_time = "03:00:00:00",
                              queue = "geospatial.q",
                              slots = 2,
                              memory = ifelse(s %in% c('bra', 'eaas'),500,200),
                              overwrite=TRUE,
                              ages=0,
                              holdouts=0,
                              regions=s,
                              corerepo=core_repo,
                              log_dir=log_dir,
                              singularity = 'default',
                              proj = 'proj_geo_nodes',
                              modeling_shapefile_version = modeling_shapefile_version,
                              raking_shapefile_version = raking_shapefile_version)
  }
}


waitforaggregation()

#launch combine aggregation
for(i in indicator){
  for(r in c('raked', 'unraked')){
    run_date <- paste0("<<<< DATE REDACTED >>>>")
    dir.create("<<<< FILEPATH REDACTED >>>>")
    dir.create("<<<< FILEPATH REDACTED >>>>")
    combine_agg_qsub <- paste0('qsub -e ', ' "<<<< FILEPATH REDACTED >>>>"',i ,'/output/', run_date, '/errors -o ', ' "<<<< FILEPATH REDACTED >>>>"',i , '/output -cwd -l m_mem_free=50G -l h_rt=00:03:00:00 -l fthread=1 -q geospatial.q -v sing_image=default -P proj_geo_nodes_edu -N combine_',
                               i, ' ', core_repo, '/mbg_central/share_scripts/shell_sing.sh ', indic_repo, 'education/combine_admin_draws.R ', run_date, " ", i, " ", r) 
    system(combine_agg_qsub)
  }
}

