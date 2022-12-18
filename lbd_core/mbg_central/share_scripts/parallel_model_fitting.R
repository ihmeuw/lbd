#################################################################################
## Generic parallel model fitting script for running MBG models             ##
#################################################################################

sessionInfo()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
doannual <- TRUE # if false will do a 5-year run (used for testing)
testing  <- FALSE # run on a subset of the data
testgauss <- FALSE # Test SBH Data As Gaussian
clamp_covs <- FALSE
## grab arguments
## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

print(reg)
print(age)
print(run_date)
print(test)
print(holdout)
print(indicator)
print(indicator_group)
## make a pathaddin that get used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

## load an image of the main environment
load(paste0(<<<< FILEPATH REDACTED >>>>'/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))

## In case anything got overwritten in the load, reload args
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
outputdir <- file.path(<<<< FILEPATH REDACTED >>>>'output',run_date,'/')
dir.create(outputdir, showWarnings = FALSE)

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

# print out session info so we have it on record
sessionInfo()

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)



source(paste0(<<<< FILEPATH REDACTED >>>>'/lbd_hiv/mbg/hiv_prev_disagg/3_functions/read_inla_prior.r'))


# We need to be in the singularity image, and specifically the LBD one if using TMB
if(!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if(as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Print the core_repo hash and check it
message("Printing git hash for 'core_repo' and checking against LBD Core Code master repo")
#record_git_status(core_repo = core_repo, check_core_repo = TRUE)

## Make sure this inla patch is implemented if running on geos
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

## cores to use
cores_to_use <- Sys.getenv("SGE_HGR_fthread")
#if (cores_to_use=="") 
cores_to_use=1
message(paste("Model set to use", cores_to_use, "cores"))

inla.setOption("pardiso.license", <<<< FILEPATH REDACTED >>>>)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PID <- Sys.getpid()
tic("Entire script") # Start master timer
## some set up
if (class(year_list)     == "character") year_list    <- eval(parse(text=year_list))
if (class(z_list)        == "character") z_list       <- eval(parse(text=z_list))
if (class(sex_list)      == "character") sex_list     <- eval(parse(text=sex_list))
if (class(interacting_gp_1_effects) == "character" & length(interacting_gp_1_effects) == 1) interacting_gp_1_effects <- eval(parse(text=interacting_gp_1_effects))
if (class(interacting_gp_2_effects) == "character" & length(interacting_gp_2_effects) == 1) interacting_gp_2_effects <- eval(parse(text=interacting_gp_2_effects))

message(paste0('interacting gp 1 effects: ', interacting_gp_1_effects))
message(paste0('interacting gp 2 effects: ', interacting_gp_2_effects))



## reload data an prepare for MBG
load(paste0(<<<< FILEPATH REDACTED >>>> '/model_image_history/', run_date, pathaddin, '.RData'))

#Load input data
input_data<-readRDS( ## save this here in case predict dies
  file = sprintf('<<<< FILEPATH REDACTED >>>>/output/%s/%s_data_input_list_%s_holdout_%i_agebin_%i.RDS',
               run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Try doing some environment cleanup before launching into fitting
if(exists('stop_after_fit')) {
  try(remove(cov_list))
  try(remove(df))
  try(remove(mesh_int))
  try(remove(mesh_s))
  try(remove(simple_raster))
  try(remove(MbgStandardCovariateLoader))
  try(remove(CovariatePathHelper))
  gc()
}
if(!exists('skipped_inla_rundate')) skipped_inla_rundate = run_date
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fit MBG model

tic("MBG - fit model") ## Start MBG - model fit timer

## Set the number of cores to be equal to input;
## If missing, then revert to cores_to_use value
if(Sys.getenv("OMP_NUM_THREADS") != "") {
  setompthreads(Sys.getenv("OMP_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input OpenMP threads.")
  setompthreads(cores_to_use)
}

if(Sys.getenv("MKL_NUM_THREADS") != "") {
  setmklthreads(Sys.getenv("MKL_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input MKL threads.")
  setmklthreads(cores_to_use)
}



if(!as.logical(skipinla)) {
  if(fit_with_tmb == FALSE) {
    message('Fitting model with R-INLA')
    
    model_fit <- fit_mbg(indicator_family = indicator_family,
                         stack.obs        = stacked_input,
                         spde             = spde_int,
                         cov              = outcome,
                         N                = N,
                         int_prior_mn     = intercept_prior,
                         f_mbg            = mbg_formula,
                         run_date         = run_date,
                         keep_inla_files  = keep_inla_files,
                         cores            = cores_to_use,
                         wgts             = weights,
                         intstrat         = intstrat,
                         verbose_output = TRUE,
                         fe_sd_prior      = 1 / 9, 
                         sparse_ordering  = as.logical(sparse_ordering))    
    
  } else {
    message('Fitting model with TMB')
    message(sprintf('%i Data points and %i mesh nodes',nrow(df),length(input_data$Parameters$Epsilon_stzx)))
    
    use_ar2=F
    cpp_filetag <- '_cre_z_x_2'
    
    message(paste0('using ',cpp_filetag))
    
    # run the model
    system.time(
      model_fit <- fit_mbg_tmb( 
        lbdcorerepo     = core_repo,
        cpp_template    = paste0('mbg_tmb_model', cpp_filetag),
        tmb_input_stack = input_data,
        control_list    = list(trace=1, eval.max=8000, iter.max=8000, abs.tol=1e-20),
        optimizer       = 'nlminb',
        ADmap_list      = NULL,
        sparse_ordering = as.logical(sparse_ordering),
        int_gp_1_effs   = interacting_gp_1_effects,
        int_gp_2_effs   = interacting_gp_2_effects)
    )
  }
  
  saveRDS(object = model_fit, ## save this here in case predict dies
          file = sprintf('<<<< FILEPATH REDACTED >>>>/%s/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                        run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
  
} else {
  model_fit <- readRDS( file = sprintf('<<<< FILEPATH REDACTED >>>>/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                                       skipped_inla_rundate, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
  saveRDS(object = model_fit, ## save this here in case predict dies
          file = sprintf('<<<< FILEPATH REDACTED >>>>/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                    run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
}
  
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up fitting~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# save training data
write.csv(
  df,
  file = (paste0(outputdir, '/', indicator,'_trainingdata',pathaddin,'.csv')),
  row.names = FALSE
)

message('Done with model fitting.')

toc(log = T) # End master timer

## Format timer
ticlog   <- tic.log(format = F)
df_timer <- generate_time_log(ticlog)
df_timer[, region := reg]
df_timer[, holdout := holdout]
setcolorder(df_timer, c("region", "holdout", "step", "time"))

## Pull run time for this run
run_time_all <- df_timer[step == "Entire script", time]

## Write to a run summary csv file in the output directory
output_file <- paste0(outputdir, "run_summary_fitting_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}


## Write CSV
write.csv(df_timer,file = output_file, row.names = FALSE)

toc(log = T) ## End MBG - model fit timer
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Launch Predict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Launch prediction launch job
  qsub <- make_qsub_share(code_path      = sprintf("%s/mbg_central/share_scripts/parallel_model_launch_prediction.R", core_repo),
                          addl_job_name  ='pred_launch',
                          age            = age,
                          reg            = reg,
                          holdout        = holdout,
                          test           = testing,
                          indic          = indicator,
                          saveimage      = FALSE,
                          memory         = 5,
                          cores          = 5,
                          proj           = 'proj_geo_nodes',
                          geo_nodes      = FALSE,
                          use_c2_nodes   = FALSE,
                          singularity    = '<<<< FILEPATH REDACTED >>>>/lbd_full_20200128.simg',
                          singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1),
                          queue          = 'geospatial.q',
                          run_time       = '05:00:00:00')
  
  system(qsub)
  message('finished with model fitting, prediction has been launched. Done with this job!')
  