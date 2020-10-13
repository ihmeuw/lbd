rm(list = ls())

source('FILEPATH.R')
source('FILEPATH.R')
source('FILEPATH.R')

print('setting config global variables... ')
configs <- read.csv('FILEPATH.csv', as.is=T)

# read in relevant year bins
raw_dat <- read.csv(paste0(configs$dat_dir[1], '/', configs$data_file[1]))
raw_dat.yrs <- na.omit(raw_dat$year)
yr_min <- min(raw_dat.yrs)
yr_max <- max(raw_dat.yrs)
if (yr_max > 2016) yr_max <- 2016
all_yrs <- c(yr_min:yr_max)

## MAKE PREDICTION BRICK RUN: skip to "LAUNCH MODELS..." if already made

# master list of covariates to use
cov_name_list <- read.csv('FILEPATH.csv', as.is=T)

# list of dirs
in_dir <- toString(configs$in_dir[1])
cov_dir <- "FILEPATH/"
stack_out_dir <- 'FILEPATH/'
project <- toString(configs$project[1])
queue <- toString(configs$queue[1])

# launch job to make stacks / covariate in master list
for (j in 1:nrow(cov_name_list)) {
  static <- as.logical(cov_name_list$static[j])
  if (!static){
    cov_name <- cov_name_list$cov_name[j]
    measure <- cov_name_list$measure_yrly[j]
    l_lim <- cov_name_list$l_lim[j]
    u_lim <- cov_name_list$u_lim[j]
    cov_stack.run <- 'FILEPATH.R'
    qsub(paste0("FILEPATH_", cov_name, "_", substr(make_time_stamp(T),6,10)),
         cov_stack.run,
         pass=c(toString(l_lim), toString(u_lim), cov_name, measure, 'FALSE', cov_dir, stack_out_dir, '0', cov_name, '0'),
         log=T,
         fthreads=1,
         m_mem_free='50G',
         archive=T,
         submit=T,
         q=queue,
         proj=project,
         user='USERNAME')
  }
}

# launch job to make prediction raster stacks for latest available year
cov_pred.run <- paste0(in_dir, '/FILEPATH.R')
qsub(paste0("FILEPATH", substr(make_time_stamp(T),6,10)),
     cov_pred.run,
     pass=c(toString(yr_max), in_dir, 'FALSE'),
     log=T,
     fthreads=1,
     m_mem_free='4G',
     archive=T,
     submit=T,
     q=queue,
     proj=project,
     user='USERNAME')

## LAUNCH BRT MODEL RUNS
for (i in 1:nrow(configs)) {
  # file location variables
  in_dir <- toString(configs$in_dir[i])
  data_file <- toString(configs$data_file[i])
  runGBM_file <- toString(configs$runGBM_file[i])
  smry_file <- toString(configs$smry_file[i])
  dat_dir <- toString(configs$dat_dir[i])
  # modeling options
  #  index of model permutations -- determined by user
  index <- toString(configs$index[i])
  #  type of model-based hyperparameter selection to perform: rf, brt, or gp
  opt_type <- toString(configs$opt_type[i])
  #  number of BRT jobs to run
  njobs <- as.numeric(configs$njobs[i])
  # qsub options
  brt_threads <- as.numeric(configs$brt_threads[i])
  brt_mem <- toString(configs$brt_mem[i])
  smry_threads <- as.numeric(configs$smry_threads[i])
  smry_mem <- toString(configs$smry_mem[i])
  project <- toString(configs$project[i])
  queue <- toString(configs$queue[i])
  
  idx_dir <- paste0(in_dir, '/', max(list.files(in_dir)[grep(index, list.files(in_dir))]))

  print('checking for or making directory tree... ')
  ## CREATE DIRECTORIES
  out_dir <- paste0(idx_dir, '/', opt_type)
  out_dirs <- c(idx_dir,
                out_dir,
                paste0(out_dir, '/FILEPATH'),
                paste0(out_dir, '/FILEPATH'),
                paste0(out_dir, '/FILEPATH'),
                paste0(out_dir, '/FILEPATH'),
                paste0(out_dir, '/FILEPATH'),
                paste0(out_dir, '/FILEPATH'),
                paste0(out_dir, '/FILEPATH'))
  for (j in 1:length(out_dirs)){
    if(!dir.exists(out_dirs[j])){
      dir.create(out_dirs[j])
    }
  }

  # Check for results - makes sure models are running and allows time for them to run
  print('checking output directories for model results... ')
  Sys.sleep(3600)
  check_loc_results(c(1:njobs), paste0(out_dir, '/FILEPATH'), prefix="FILEPATH",postfix=".tif")

  ## SUMMARIZE BRT RUNS
  print('summarizing model runs... ')
  sry.run <- paste0(in_dir, '/', smry_file)
  qsub("ADDRESS",
       sry.run,
       pass=c(in_dir, out_dir, yr_max),
       log=T,
       fthreads=smry_threads,
       m_mem_free=smry_mem,
       archive=F,
       submit=T,
       q=queue,
       proj=project,
       user="USERNAME")
}

