## Set repo location and indicator group
user            <- Sys.info()['user']
makeholdouts    <- F
run_models      <- F
run_postest     <- T
postest_rd      <- '<<<< RUN DATE REDACTED >>>>'
rerun           <- F
rerun_rd        <- '<<<< RUN DATE REDACTED >>>>'
rerun_sp        <- 3
possible_permutations    <- c("", 'covs', 'stacking', 'gpcovs','gponly', '1000_draws') #""
desired_permutations <- c('covs', 'stacking', 'gpcovs','gponly')
mod_slots = 30; pe_slots = 35;
country_models  <- F
geos_nodes <- F
prev_only = F
user_repo = '<<<< FILEPATH REDACTED >>>>'
plib = '<<<< FILEPATH REDACTED >>>>'

#rerun_broken
only_broken = #boolean
if(rerun){
  source('postproc_functions.R')
  source(paste0('model_status.R'))
}

#launch helper files
source(paste0(user_repo, 'shared/helpers_launch.R'))

config_num = 17
brt_num = 2
setup_models('lri',paste0(config_num, '_', brt_num),'8')

#stacking instructions
.libPaths(plib)
library('mbgstacking')
gam_model = init_gam(model_name = 'gam', arguments = list(spline_args = list(bs = 'ts', k = 6)))
enet_model = init_penalized('enet', arguments = list(alpha = .5))
brt_model1 = init_xgboost(model_name = 'gbm1', params_arg = list(eta = .007, max_depth = 8, subsample = .7, colsample_bytree = .5),
                          nrounds = 1000, binomial_evaluation = 'prev', weight_column = 'alt_weight')
brt_model2 = init_xgboost(model_name = 'gbm2', params_arg = list(eta = .06, max_depth = 16, min_child_weight = 1, subsample = .7, colsample_bytree = 1),
                          nrounds = 10000, binomial_evaluation = 'prev', weight_column = 'alt_weight')

#summary statistics
summstats <- c('mean', 'cirange', 'upper', 'lower')


## drive locations
root           <- '<<<< FILEPATH REDACTED >>>>'
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'

#load basic functions and what not
core_repo = '<<<< FILEPATH REDACTED >>>>'
library(data.table)
if(Sys.info()[1] == 'Windows'){
  stop('STOP! you will overwrite these packages if you run from windows\n
       STOP! also, lots of this functions wont work so get on the cluster!')
} else {
  library('data.table')
  library('magrittr')
  source(paste0(core_repo,'mbg_central/mbg_functions.R'))
  source(paste0(core_repo,'mbg_central/prep_functions.R'))
  source(paste0(core_repo,'mbg_central/covariate_functions.R'))
  source(paste0(core_repo,'mbg_central/misc_functions.R'))
  source(paste0(core_repo,'mbg_central/post_estimation_functions.R'))
  source(paste0(core_repo,'mbg_central/gbd_functions.R'))
  source(paste0(core_repo,'mbg_central/holdout_functions.R'))
  source(paste0(core_repo,'mbg_central/validation_functions.R'))
  source(paste0(core_repo,'mbg_central/validation_report_functions.R'))
  source(paste0(core_repo,'mbg_central/seegMBG_transform_functions.R'))
  source('get_outputs.R') #shared fucntion to get GBD estimates
}

if(run_models & !rerun){
  ## Load libraries and  MBG project functions.
  #.libPaths(package_lib)
  package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))


  ## Read config file and save all parameters in memory
  config <- load_config(repo            = user_repo,
                        indicator_group = indicator_group,
                        indicator       = indicator,
                        config_name     = config_name,
                        covs_name       = c_list)


  if(as.logical(brts2)){
    stacking_models = list(gam_model, enet_model, brt_model1, brt_model2)
  }else{
    stacking_models = list(gam_model, enet_model, brt_model1)
  }
  ## Create run date in correct format
  run_date <- make_time_stamp(TRUE)

  #make new folder
  dir.create(sprintf('%s/output/%s',sharedir,run_date))

  # Load region list from config file
  if(!country_models){
    region_list = unlist(strsplit(region_list, "+", fixed = T))
    region_list = sapply(region_list, trimws)

  } else{
    df <- load_input_data(indicator   = indicator,
                          simple      = NULL,
                          removeyemen = TRUE,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share))
    df[, .N, by = country]
    region_list = as.character(df[, .N, by = country][N>10, country])
    print(region_list)
  }
  Regions <- region_list
  strata = Regions

  ## Ensure you have defined all necessary settings in your config
  #check_config()

  #set things properly
  #constrain_coefs = as.logical(constrain_coefs)

  ###############################################################################
  ## Make Holdouts
  ###############################################################################
  if(makeholdouts){
    # load the full input data
    for(package in package_list){
      library(package, character.only=TRUE)
    }

    df <- load_input_data(indicator   = indicator,
                          simple      = NULL,
                          removeyemen = TRUE,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share))

    # add in location information
    df <- merge_with_ihme_loc(df)

    #make period the year col vis a vis hold outs
    if(yr_col == 'ppp'){
      df[year<=2003, ppp:= 1]
      df[is.na(ppp) & year<= 2007, ppp:= 2]
      df[is.na(ppp) & year<= 2011, ppp:= 3]
      df[is.na(ppp) , ppp:= 4]
    }

    shape_dir <- '' #path to shapefiles

    if(!exists('spatial_ho')) spatial_ho = 'qt'
    if(spatial_ho == 'nid'){
      stratum_ho = lapply(unique(df$region[!is.na(df$region)]), function(x){

        dat = copy(df[region == x,])
        nids = unique(dat[,nid])

        #randomly rearrange
        nids = sample(nids, size = length(nids))

        #chop into 5 folds
        foldz = split(nids, ceiling(seq_along(nids)/(length(nids)/5)))

        iter = 1
        for(fdz in foldz){
          dat[nid %in% fdz, fold:=iter]
          iter = iter + 1
        }

        dat[, ho_id := fold]

        return(dat)

      })

      names(stratum_ho) = paste0('region__',unique(df$region[!is.na(df$region)]))
    }else{
      stratum_ho <- make_folds(data       = df,
                               n_folds    = as.numeric(n_ho_folds),
                               spat_strat = spatial_ho,
                               temp_strat = 'prop',
                               strat_cols = 'region',
                               ts         = as.numeric(ho_ts),
                               mb         = as.numeric(ho_mb),
                               lat_col    = "latitude",
                               long_col   = "longitude",
                               ss_col     = "N",
                               admin_shps   = paste0(shape_dir, 'ad1_raster.grd'),
                               admin_raster = paste0(shape_dir, 'africa_ad1.shp'),
                               mask_shape   = paste0(shape_dir, 'africa_simple.shp'),
                               mask_raster  = paste0(shape_dir, 'ad0_raster'))
    }
  }

  ###############################################################################
  ## Format shared files/variables/etc.
  ###############################################################################

  ## Save to use in producing aggregated fit statistics
  dir.create(paste0(sharedir, '/output/', run_date, '/fit_stats'),recursive = T)
  save(strata, file = '<<<< FILEPATH REDACTED >>>>')

  g2l <- fread("/snfs1/WORK/11_geospatial/10_mbg/gaul_to_loc_id.csv")
  locs = unique(g2l[GAUL_CODE %in% get_gaul_codes('africa'),loc_id])

  ## Load GBD Estimates for this indicator which will be used in raking
  severe <- .151

  fetch_raking_targets = function(measure_id, g2l){
    #takes year list and gbd id from the global environment because I'm lazy
    dat = get_outputs(topic = 'cause', measure_id = measure_id, location_id = locs, year_id = eval(parse(text = year_list)), sex_id = 3, age_group_id = 1, cause_id = gbd_id,
                      gbd_round_id = 4, metric = 3)

    #convert to severe incidence; stolen from the global environment
    if(measure_id == 6){
      dat[,val:= val * severe]
    }

    dat = merge(dat, g2l, by.x= 'location_id', by.y = 'loc_id', all.x = T)
    setnames(dat, c('GAUL_CODE', 'year_id', 'val'), c('name','year', 'mean'))

    dat = dat[ ,.(name, year, mean, measure_id, metric_id)]

    return(dat)

  }

  rake_targets = lapply(c(1:6), function(x) fetch_raking_targets(x,g2l))
  names(rake_targets) = c('mortality','daly','yld','yll', 'prevalence', 'incidence')




  ## Make loopvars aka strata grid (format = regions, ages, holdouts)
  if(makeholdouts){
    loopvars <- expand.grid(reg = Regions, age = 0, holdout = 0:n_ho_folds, test = F, start_point = 1, stringsAsFactors = F, permutations = possible_permutations)
  } else {
    loopvars <- expand.grid(reg = Regions, age = 0, holdout = 0, test = F, start_point = 1, stringsAsFactors = F, permutations = possible_permutations)
  }

  loopvars$id = 1:nrow(loopvars)

  #save some stuff for post estimation
  dir.create('<<<< FILEPATH REDACTED >>>>',recursive = T)
  prep_postest(indicator = indicator,
               indicator_group = indicator_group,
               run_date = run_date,
               save_objs = c("rake_targets", "year_list", "summstats", "rake_transform", 'pop_measure'))

  #stacking instructions moved to the top of the script

  #cross val rounds
  nfc = 1
  nf = 5

  r_path = '<<<< FILEPATH REDACTED >>>>' #path to r shell

  #things to save in the temp image
  keepme = unique(c('config', 'df', 'g2l', 'loopvars','rake_targets',ifelse(exists('stratum_ho'), 'stratum_ho', ""),config$V1, 'run_date', summstats,
                    'stacking_models', 'nfc','nf','r_path'))
  keepme = keepme[!keepme %in% ""]
  #save temprun
  save(list = keepme,
       file = '<<<< FILEPATH REDACTED >>>>')
} else{
  load('') #load prerun temp image

  if(!any(grepl('permutations', names(loopvars)))){
    loopvars$permutations = ""
  }

  if(rerun){
    if(all(desired_permutations %in% "")){
      loopvars[, c('permutations','original_rundate', 'id')] = NULL
      loopvars = unique(loopvars)
      expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))
      loopvars = expand.grid.df(loopvars, data.frame(permutations = possible_permutations, stringsAsFactors = F),data.frame(original_rundate = rerun_rd, stringsAsFactors = F))
    }

    loopvars$id = 1:nrow(loopvars)

    loopvars$start_point <- rerun_sp
    #things to save in the temp image
    keepme = unique(c('config', 'df', 'g2l', 'loopvars','rake_targets',ifelse(exists('stratum_ho'), 'stratum_ho', ""),config$V1, 'run_date', summstats,
                      'stacking_models', 'nfc','nf','r_path'))
    keepme = keepme[!keepme %in% ""]
    save(list = keepme, file = '<<<< FILEPATH REDACTED >>>>')
  }


}
#prep things

#make qsub
share_scripts = '<<<< FILEPATH REDACTED >>>>'
shelldir = '<<<< FILEPATH REDACTED >>>>'
#determine start points
dir.create(paste0(sharedir, '/output/', run_date, '/checkpoints/'),recursive = T)
loopvars = loopvars[loopvars$permutations %in% desired_permutations,]
loopvars$full_rd = paste0(loopvars$original_rundate,ifelse(loopvars$permutations == "", "", paste0('_', loopvars$permutations)))
if(only_broken){
  mod_stat = rbindlist(lapply(unique(loopvars$full_rd), function(x) model_status(x, 0:5, rake_targets = 'prevalence', indicator, indicator_group, T)))

  #keep where rawcellpred didn't get created
  mod_stat = mod_stat[rawcellpred == F,]

  loopvars = merge(loopvars, mod_stat, by.x=c('full_rd','reg','holdout'), by.y = c('run_date', 'region','holdout'))
}


#for each possible model
#to do, change to be keyed off loopvars id
source('misc_functions.R')
for(i in unique(loopvars$id)){
  therd = ifelse(run_models | rerun, run_date, postest_rd)
  if(loopvars$permutations[loopvars$id == i] != ""){
    therd = paste0(therd, '_', loopvars$permutations[loopvars$id == i])
  }
  if(run_models | rerun){
    model_qsub = build_qsub(job_name = paste(loopvars$reg[loopvars$id == i],paste0(indicator, '_mbg'), loopvars$holdout[loopvars$id == i], sep = "_"),
                            output_folder = '<<<< FILEPATH REDACTED >>>>',
                            error_folder ='<<<< FILEPATH REDACTED >>>>',
                            make_folders = T,
                            shell_path = paste0(shelldir, 'shell_sing.sh'),
                            script_path = paste0(share_scripts, 'parallel_model_newstacking.R'),
                            additional_options = ' -v sing_image=default',
                            project = '<<<< PROJECT NAME REDACTED >>>>', #project name
                            geos_node = as.logical(geos_nodes),
                            slots = mod_slots,
                            num_tasks = 0,
                            arguments_string = paste(indicator_group, indicator, run_date, i, 'blarg'))
    #print(model_qsub)
    job = as.numeric(system(model_qsub, intern = T))
    print(job)
  }

  if(loopvars$holdout[loopvars$id == i]==0 & run_postest){ #& loopvars$permutations[loopvars$id == i] %in% c("", '1000_draws')
    #launch post estimation for full models
    postest_qsub = build_qsub(job_name = paste0('postest_',indicator, paste(loopvars$reg[loopvars$id == i])),
                              output_folder = '<<<< FILEPATH REDACTED >>>>',
                              error_folder = '<<<< FILEPATH REDACTED >>>>',
                              make_folders = T,
                              shell_path = paste0(shelldir, 'shell_sing.sh'),
                              script_path = paste0(share_scripts, 'postest_script.R'),
                              additional_options = paste(switch(exists('job'), paste0('-hold_jid ', job), NULL),
                                                         paste0('-v SET_MKL_THREADS=',floor(4/2), ' -v sing_image=default')),
                              project = '<<<< PROJECT NAME REDACTED >>>>', #project name
                              geos_node = geos_nodes,
                              slots = pe_slots,
                              num_tasks = 0,
                              arguments_string = paste(loopvars$reg[loopvars$id == i], therd, indicator, indicator_group, F, prev_only, 'blarg'))
    #print(postest_qsub)
    system(postest_qsub)

  }

}
print(run_date)
