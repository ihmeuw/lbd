#################################################################################
## Generic parallel model prediction script for running MBG models             ##
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
load(paste0(<<<< FILEPATH REDACTED >>>> '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))

## In case anything got overwritten in the load, reload args
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
outputdir <- file.path(<<<< FILEPATH REDACTED >>>> 'output',run_date,'/')
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



source(paste0(<<<< FILEPATH REDACTED >>>> '/lbd_hiv/mbg/hiv_prev_disagg/3_functions/read_inla_prior.r'))
source(paste0(indic_repo, 'mbg/', indicator, '/3_functions/waitformodelstofinish_prediction.r'))

if(use_gbd_covariate==T & !exists('create_gbd_covariate_raster')) create_gbd_covariate_raster <- T



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

if(!exists('agesex_in_stacking')) agesex_in_stacking = F
if(!exists('remake_adult_prev_stackers')) remake_adult_prev_stackers = F
if(!exists('skipped_inla_rundate')) skipped_inla_rundate = run_date
if(!exists('use_adult_prev_stackers'))  use_adult_prev_stackers = FALSE

message(paste0('interacting gp 1 effects: ', interacting_gp_1_effects))
message(paste0('interacting gp 2 effects: ', interacting_gp_2_effects))



## reload data an prepare for MBG
load(paste0(<<<< FILEPATH REDACTED >>>> '/model_image_history/', run_date, pathaddin, '.RData'))


toc(log = T) ## End MBG - model fit timer

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Predict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tic("MBG - predict model") ## Start MBG - model predict timer

## Run predict_mbg on chunks of 50 samples (to avoid memory issues)
message('Making predictions in 50 draw chunks.')

max_chunk <- 50
samples   <- as.numeric(samples)

## Create vector of chunk sizes
chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)


memory <- c(cssa = 400, 'essa_sdn-COM' = 600, sssa = 280, 'wssa-CPV-STP-MRT' = 600)[reg]
if(is.na(memory)|z_list==0) memory = 300


#Launch all the chunks
for(chunk in 1:length(chunks)){

    qsub <- make_qsub_share(code_path      = sprintf("%s/mbg_central/share_scripts/parallel_model_chunk_prediction.R", core_repo),
                            addl_job_name  = paste0('pred_chunk_', chunk),
                            age            = age,
                            reg            = reg,
                            holdout        = holdout,
                            test           = testing,
                            indic          = indicator,
                            saveimage      = FALSE,
                            memory         = memory,
                            cores          = 5,
                            pred_chunk     = chunk,
                            proj           = 'proj_geo_nodes',
                            geo_nodes      = FALSE,
                            use_c2_nodes   = FALSE,
                            singularity    = '<<<< FILEPATH REDACTED >>>>/lbd_full_20200128.simg',
                            singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1),
                            queue          = 'geospatial.q',
                            run_time       = '05:00:00:00')
    
    system(qsub)
  }


#wait for the prediction chunks to finish
waitformodelstofinish_prediction(lv = 1:length(chunks), sleeptime = 60)
############################################

memory <- c(cssa = 160, 'essa_sdn-COM' = 250, sssa = 100, 'wssa-CPV-STP-MRT' = 200)[reg]
if(is.na(memory)|z_list==0) memory = 300

#Launch the combination job
qsub <- make_qsub_share(code_path      = sprintf("%s/mbg_central/share_scripts/parallel_model_combine_prediction.R", core_repo),
                        addl_job_name  = 'pred_combine',
                        age            = age,
                        reg            = reg,
                        holdout        = holdout,
                        test           = testing,
                        indic          = indicator,
                        saveimage      = FALSE,
                        memory         = memory,
                        cores          = 5,
                        pred_chunk     = chunk,
                        proj           = 'proj_geo_nodes',
                        geo_nodes      = FALSE,
                        use_c2_nodes   = FALSE,
                        singularity    = '<<<< FILEPATH REDACTED >>>>/lbd_full_20200128.simg',
                        singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1),
                        queue          = 'geospatial.q',
                        run_time       = '05:00:00:00')

system(qsub)

message('Prediction combining stage launched!')
#########################################