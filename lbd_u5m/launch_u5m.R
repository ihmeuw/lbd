## #####################################################################################
## AUTHOR:  Roy Burstein
## DATE:    2018
## PURPOSE: U5M Launch Script
## #####################################################################################


## #####################################################################################
## SETUP

# some bash commands to clean up system and pull any new code
rm(list=ls())

doannual <- TRUE # if false will do a 5-year run (used for testing)
testing  <- FALSE # run on a subset of the data
testgauss <- FALSE # Test SBH Data As Gaussian

## Set repo location and indicator group
user            <- Sys.info()['user']                                  # user running this code
core_repo       <- '<<<< FILEPATH REDACTED >>>>'
ig_repo         <- '<<<< FILEPATH REDACTED >>>>'                       # u5m specific repo
core_remote     <- 'origin'                                            # gits name of lbd_core upstream repo (should stay origin)
core_branch     <- 'dev-u5m'                                           # working branch to pull lbd_core from stash
ig_remote       <- 'origin'                                            # gits name of u5m upstream repo (should stay origin)
ig_branch       <- 'develop'                                           # woring branch to pull u5m from stash  
indicator_group <- 'u5m'                                               # indicator group is u5m, this shouldn't change
indicator       <- 'died'                                              # indicator is died, this shouldn't change 
Regions <- regions <- c(
  'soas',
  'ocea+seas-mys',
  'stan+mng',
  'caca-mex', # Removing Mexico from Stage 2 paper
  'ansa+trsa-bra', # Removing Brazil from Stage 2 paper
  'noaf',
  'mide+yem',
  'cssa',
  'essa-yem',
  'sssa',
  'wssa'
)

makeholdouts    <- TRUE   # set to true if this is an out of sample run
pullgit         <- FALSE  # set to true if you will pull into the repo on share before running coe
numagebins      <- 7      # this depends on data prep, right now its 5, but it will likely change to 7
custom          <- FALSE  # specific for post estimation in India as MC set up

## IF YOU WANT TO PRESET TO AN OLD RUN DATE DO IT HERE, ELSE SET IT NULL AND IT WILL SPIN UP A NEW RUN DATE(see run_date_notes.csv)
preset_run_date <- '2019_06_07_rtr'

## If you want to use an OOS fold table from a different run, list it here
## WARNING: This will fail in the share script if new sources have been added 
##  since the reference model run date.
preset_holdout_date <- '2019_06_07_rtr'

## sort some directory stuff and pull newest code into share
if(pullgit) {
  system(sprintf('cd %s\ngit pull %s %s', core_repo, core_remote, core_branch))
  system(sprintf('cd %s\ngit pull %s %s', ig_repo,   ig_remote,   ig_branch))
}



## Create run date in correct format
if(is.null(preset_run_date)){
  run_date <- make_time_stamp(TRUE)
  message(sprintf('New run_date: %s', run_date))
} else {
  run_date <- preset_run_date
  message(sprintf('Preset run_date: %s', run_date))
}

## Create config and covariate config names
config_name <- paste0('config_died_', run_date) 
covs_config_name <- paste0('covs_died_', run_date)

# run some setup scripts
source(sprintf('%s/setup.R',ig_repo)) 

## #####################################################################################




## #####################################################################################
## LOAD UP DATA AND MAKE VALIDATION FIELDS
## TODO UPDATE THIS LATER, CURRENTLY HOLD OUT IS NOT FUNCATIONAL

if(makeholdouts){
  if( !is.null(preset_holdout_date) ){
    # If a previous OOS run date has been selected, use that
    message(sprintf("Using OOS table from run date %s.",preset_holdout_date))
    oos_fold_table <- fread('<<<< FILEPATH REDACTED >>>>')
  } else {
    # Make holdouts, by NID, save as a table of NIDs to be held out in the parallel script
    d <- readRDS('<<<< FILEPATH REDACTED >>>>')
    if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
    d <- subset(d, year %in% year_list)
    # Create a new holdout table, using sampling of holdouts WITHOUT replacement
    #  within each modeling region. For more information, see documentation for
    #  `make_holdouts_by_region()` in "u5m_specific_functions.R"
    oos_fold_table <- make_holdouts_by_region(
      in_data = d, 
      regions = Regions,
      num_holdouts = as.numeric(n_ho_folds)
    )
    rm('d')
  }
  write.csv(
    oos_fold_table,
    '<<<< FILEPATH REDACTED >>>>',
    row.names=FALSE
  )

  # Looping variables for qsub
  loopvars <- expand.grid(Regions,0, 0:as.numeric(n_ho_folds))
} else {
  loopvars <- expand.grid(Regions,0,0)
}

# FOR COMBINED AGE BIN MODEL AGE NEEDS TO BE ZERO


## #####################################################################################
## Copy config file and covariates config file to the run directory
## `config_name`, `sharedir`, `run_date`, and `covs_config_name` are all defined
##  in ./setup.R
message("Copying config and covariate config files to the run directory.")
file.copy(
  from = sprintf('%s/%s.csv', ig_repo, config_name), 
  to = sprintf('%s/output/%s/config.csv', sharedir, run_date), 
  overwrite = TRUE
)
file.copy(
  from = sprintf('%s/%s.csv', ig_repo, covs_config_name), 
  to = sprintf('%s/output/%s/covariates_config.csv', sharedir, run_date), 
  overwrite = TRUE
)


## #####################################################################################
##  LAUNCH MBG  
## FOR NOW RUN AGES SEPERATELY WHILE WORKING ON UPDATING A SIMULATENOUS MODEL IN TMB
## #####################################################################################

queue <- 'geospatial.q'
samples <- as.numeric(samples)
# multiplier for smaller sample size run
msampleMult <- .25^(samples <= 100)


loopvars$pred_slots <- 49
# We want region specific pred mem requests
loopvars$pred_mem <- case_when(
    loopvars$Var1 == 'soas' ~ 520 * msampleMult,
    loopvars$Var1 == 'ocea+seas-mys' ~ 666 * msampleMult,
    loopvars$Var1 == 'stan+mng' ~ 800 * msampleMult,
    loopvars$Var1 == 'caca-mex'~ 150 * msampleMult,
    loopvars$Var1 == 'ansa+trsa-bra' ~ 750 * msampleMult, 
    loopvars$Var1 == 'noaf' ~ 1000 * msampleMult,
    loopvars$Var1 == 'mide+yem' ~ 500 * msampleMult,
    loopvars$Var1 == 'cssa' ~ 700 * msampleMult,
    loopvars$Var1 == 'essa-yem' ~ 1000 * msampleMult,
    loopvars$Var1 == 'sssa' ~ 420 * msampleMult,
    loopvars$Var1 == 'wssa' ~ 1000 * msampleMult,
    TRUE ~ 1000 * msampleMult
)

# loopvars[ loopvars$Var1 =='caca-mex', 'pred_slots'] <- 20
# loopvars[ loopvars$Var1 =='caca-mex', 'pred_mem'] <- 300


for(i in 1:nrow(loopvars)){
  
  # setup stratification variables
  reg         <- as.character(loopvars[i,1])
  age         <- loopvars[i,2]
  holdout     <- loopvars[i,3]
  pred_slots  <- unname(loopvars[i, 'pred_slots'])
  pred_mem    <- unname(loopvars[i, 'pred_mem'])
  message(paste('\n',age,reg,holdout))
  
  make_qsub_split_jobs(
    reg              = reg,
    age              = age,
    holdout          = holdout,
    corerepo         = core_repo,
    indic_repo       = ig_repo,
    cores_fit        = 10,
    cores_predict    = pred_slots,
    memory_fit       = 200,
    memory_predict   = pred_mem,
    proj             = 'proj_geo_nodes_u5m',
    indic            = ifelse(age==0,indicator,paste0(indicator,'_age',age)), # 0 = age combined model,
    rd               = run_date,
    log_location     = 'sharedir',
    saveimage        = TRUE,
    test             = FALSE,
    geo_nodes        = TRUE,
    use_c2_nodes     = FALSE,
    singularity_opts = list(SET_OMP_THREADS=7, SET_MKL_THREADS=7),
    predict_only     = TRUE
  )

}







## #######################################################
## POST ESTIMATION
## Launch Post estimation once the model jobs have finished
## Using misc_functions.R::parallelize() to launch post estimation by region

# what regions have completed running?
compregs <- completed_regions(rd=run_date)
for_postest <- regions_needing_postest(rd=run_date)

# some options to pass to post estimation script
ageasindic                  <- TRUE   # should always be true, false is now depreciated
abs                         <- 1:7    # age bins in use for this model run
ab_fracs                    <- list(neonatal = 1,     # reporting age bins composed of which analytical age groups?
                                    infant   = 1:3,
                                    under5   = 1:7)
agegroups                   <- names(ab_fracs)  # names of reporting age bins
admin_levels_to_aggregate   <- 0:2 # can be any combination of 0,1,2
get_summary_parameter_table <- TRUE
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (class(z_list) == "character") z_list <- eval(parse(text=z_list))

parallelize(user             = user,
            slots            = 10,
            memory           = 150,
            script           = 'postestimation_u5m.R', #_deathsonly
            ig               = indicator_group,
            indic            = indicator,
            rd               = run_date,
            expand_vars      = list(region = for_postest),
            save_objs        = c('run_date','ageasindic','year_list','doannual','testing',
                                 'abs','ab_fracs','agegroups','admin_levels_to_aggregate','interval_mo','get_summary_parameter_table',
                                 'config_name', 'covs_config_name'),
            script_dir       = ig_repo,
            prefix           = 'job',
            log_location     = '<<<< FILEPATH REDACTED >>>>',
            corerepo         = core_repo,
            proj             = 'proj_geo_nodes_u5m',
            geo_nodes        = TRUE,
            singularity      = 'default',
            singularity_opts = list(SET_OMP_THREADS=5, SET_MKL_THREADS=1)) # MKL must be 1 for now for mclapply bug



## #######################################################
## COMBINE AND MERGE OUTPUTS TO GLOBAL
## Once Region wise post estimation is done:
## Run a script that summarizes and merges together outputs across regions for global cohesive outputs
source(sprintf('%s/combine_summarize_outputs.R',ig_repo))



## #######################################################
## CHECK FOR MISSING DATA IN PARALLEL MODEL
## check to see if any data was dropped between the start of parallel_model and model run

#any countries that shouldn't be compared, not really necessary normally because data is subset to the countries in the region
countries_to_drop <- NULL

#input data for the model run
input_data_path <- '<<<< FILEPATH REDACTED >>>>'

#place for missing data to be saved
outpath <- '<<<< FILEPATH REDACTED >>>>'

#used for parallelized reading in data
cores <- length(regions)

parallelize(user             = user,
            slots            = cores,
            memory           = 50,
            script           = 'check_for_dropped_surveys_parallel_launch.R',
            ig               = indicator_group,
            indic            = indicator,
            rd               = run_date,
            expand_vars      = list(level = c(1)),
            save_objs        = c('core_repo', 'ig_repo', 'run_date', 'regions', 'cores', 'countries_to_drop', 'input_data_path', 'outpath'),
            script_dir       = sprintf("%s/launch_scripts", ig_repo),
            prefix           = 'missing_data',
            log_location     = '<<<< FILEPATH REDACTED >>>>',
            corerepo         = core_repo,
            proj             = 'proj_geo_nodes',
            geo_nodes        = TRUE,
            singularity      = 'default',
            run_time         = "3:00",
            threads          = cores,
            singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1))


## #######################################################
## PREDICTIVE VALIDITY METRICS

# IS, OOS, or BOTH
oos_fold_table<-fread('<<<< FILEPATH REDACTED >>>>')
if(0 %in% unique(loopvars[,3]) & length(loopvars[,3])>1){
  data_sample <- c('IS','OOS')
} else if(unique(loopvars[,3])==0) {
  data_sample <- 'IS'
  oos_fold_table <- NULL
} else if(unique(loopvars[,3])!=0) {
  data_sample <- 'OOS'
} else {
  stop('Not sure what to do with OOS given this loopvars input')
}

## First get in sample out of sample draws at each data point
## For now doing this by region and reporting age bin
## Should take 20 minutes or so if all goes well
# could just do under5 then break up stuff in pv table?
parallelize(user             = user,
            slots            = 10,
            memory           = 100,
            script           = 'get_is_oos_data_draws.R', #_deathsonly
            ig               = indicator_group,
            indic            = indicator,
            rd               = run_date,
            expand_vars      = list(region = Regions, agebin = 'under5'), # always do under5 so it has all age bins within it
            save_objs        = c('run_date','year_list','ab_fracs','makeholdouts','doannual','testing','numagebins','oos_fold_table','data_sample'),
            script_dir       = ig_repo,
            prefix           = 'j_isoosdraws_',
            log_location     = '<<<< FILEPATH REDACTED >>>>',
            corerepo         = core_repo,
            proj             = 'proj_geo_nodes_u5m',
            geo_nodes        = TRUE,
            singularity      = 'default',
            singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1)) # MKL must be 1 for now for mclapply bug


## Once those have run, get the pv table outputs for this run

# run for each analytical age group. save to under5 9since this abfrac holds all age bins
get_pv_table_u5m(agebin = 'under5', rd = run_date, collapse_to_agebin = FALSE, aggregate_on = c('ad0','ad1','ad2'))

# run one for each reporting age bin ( this could be made more efficient, its repeating a lot of data loading.)
for(a in c('under5','infant','neonatal'))
  get_pv_table_u5m(agebin = a, rd = run_date, collapse_to_agebin = TRUE, aggregate_on = c('ad0','ad1','ad2'))




## #######################################################
## MAKE ADMIN LEVEL VIZ PLOTS
## Once Region wise post estimation is done:
## Run a script that produces admin0, 1,and 2 plots with data
source(sprintf("%s/misc/subnational_viz_prep.R", ig_repo))
source(sprintf("%s/mbg_central/setup.R", core_repo))

#T to make raked viz plots, F to make unraked, c(T,F) to make both
rake <- c(FALSE, TRUE)

#get admin levels to make plots of
#tied to admin_levels_to_aggregate because this data from post-est is needed to make the plots
plot_levels <- c("ad0")
if(1 %in% admin_levels_to_aggregate) {
  plot_levels <- c(plot_levels, "ad1")
}
if(2 %in% admin_levels_to_aggregate) {
  plot_levels <- c(plot_levels, "ad2")
}

#make the plots; note this may take several hours
for(type in rake) {
  run_u5m_subnat_vis(reg           = compregs,
                     agegroups     = c('under5','neonatal','infant'),
                     run_date      = run_date,
                     plot_levels   = plot_levels,
                     outdir        = NULL,
                     legend_max    = NULL,
                     raked         = type,
                     with_data     = TRUE,
                     tol           = 0.0003)
}


## #######################################################
## MAKE MAPPING_DATA_FULL
## Once Region wise post estimation and combination is done:
## Run a script that produces mapping data full - contains
## pred derivative rasters, including aroc rasters used 
## for projecting estimates below

# Should countries not included in GBD (ESH Western Sahara and GUF French Guiana) be included in the map?
include_non_gbd <- TRUE

prep_mex_bra <- TRUE

#takes 17G memory with 100 draws

parallelize(user             = user,
            slots            = 1,
            memory           = 100,
            script           = 'make_mapping_data_full.R',
            ig               = indicator_group,
            indic            = indicator,
            rd               = run_date,
            expand_vars      = list(mapping=c(1)),
            save_objs        = c('core_repo', 'ig_repo', 'run_date', 'modeling_shapefile_version', 'include_non_gbd', 'prep_mex_bra'),
            script_dir       = sprintf("%s/launch_scripts", ig_repo),
            prefix           = 'mapping_data',
            log_location     = '<<<< FILEPATH REDACTED >>>>',
            corerepo         = core_repo,
            proj             = 'proj_geo_nodes',
            geo_nodes        = TRUE,
            singularity      = 'default',
            threads          = 1,
            run_time         = "25:00",
            singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1))


## #######################################################
## MAKE PROJECTED ESTIMATES
## Once Region wise post estimation and combination is done:
## Run a script that produces raked projections up to 2030

#years in the cell pred
cell_pred_year_list <- c(2000:2017)
#years to project
proj_year_list <- c(2018:2030)
#for mclapply over regions
cores <- 13
#population to use to aggregate - base worldpop goes up to 2030, worldpop_raked only goes up to 2017
pop_cov <- "worldpop"

regions <- Regions

#for 100 draws, memory = 280G. Using less is ok, mclapply does not go over memory limit.

parallelize(user             = user,
            slots            = cores,
            memory           = 1000,
            script           = 'make_projections.R',
            ig               = indicator_group,
            indic            = indicator,
            rd               = run_date,
            expand_vars      = list(age = c("under5", "infant", "neonatal")),
            save_objs        = c('core_repo', 'ig_repo', 'run_date', 'modeling_shapefile_version', 'cell_pred_year_list', 'proj_year_list', 'regions', 'cores', 'pop_cov'),
            script_dir       = sprintf("%s/launch_scripts", ig_repo),
            prefix           = 'projection',
            log_location     = '<<<< FILEPATH REDACTED >>>>',
            corerepo         = core_repo,
            proj             = 'proj_geo_nodes',
            geo_nodes        = TRUE,
            singularity      = 'default',
            threads          = cores,
            run_time         = "40:00",
            singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1))


## #######################################################
## MAKE POPULATION AGGREGATES
## Once Region wise post estimation and combination is done:
## Run a script that produces worldpop population aggregates

# worldpop measure to aggregate (probably won't change but "total" or other valid measures will work)
measures <- "a0004t"

# one of "raked", "unraked", "combined" - "raked" uses worldpop_raked (raked to GBD), which as of 5/10/2019 does not have any years after 2017. "unraked" uses plain worldpop, which has projected estimates out to 2030. "Combined" uses worldpop_raked up until 2017, then worldpop afterwards. Once worldpop_raked has been updated to have projections, the default can be updated to "raked"
rake_type <- "combined"

# level to aggregate to, 0, 1, 2, or any combination. Splits out into separate jobs if multiple.
agg_level <- 2

# years to get population for
year_list <- c(2000:2020)

#no parallelizing at the moment, but is an argument in necessary function
cores <- 1

parallelize(user             = user,
            slots            = cores,
            memory           = 35,
            script           = 'make_worldpop_aggregates.R',
            ig               = indicator_group,
            indic            = indicator,
            rd               = run_date,
            expand_vars      = list(level = agg_level),
            save_objs        = c('core_repo', 'ig_repo', 'run_date', 'modeling_shapefile_version', 'measures', 'year_list', 'cores', 'rake_type'),
            script_dir       = sprintf("%s/launch_scripts", ig_repo),
            prefix           = 'pop_agg',
            log_location     = '<<<< FILEPATH REDACTED >>>>',
            corerepo         = core_repo,
            proj             = 'proj_geo_nodes',
            geo_nodes        = TRUE,
            singularity      = 'default',
            run_time         = "3:00",
            threads          = cores,
            singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1))


