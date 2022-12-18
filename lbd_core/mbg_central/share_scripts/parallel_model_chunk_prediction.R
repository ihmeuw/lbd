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
chunk                    <- as.character(commandArgs()[11])

print(reg)
print(age)
print(run_date)
print(test)
print(holdout)
print(indicator)
print(indicator_group)
print(chunk)
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
chunk                    <- as.character(commandArgs()[11])
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
pathaddin_chunk <- paste0('_bin',age,'_',reg,'_',holdout, '_', chunk)
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



source(paste0(<<<< FILEPATH REDACTED >>>> '/lbd_hiv/mbg/hiv_prev_disagg/3_functions/read_inla_prior.r'))

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
if (class(chunk)      == "character") chunk <- as.numeric(chunk)

if(!exists('agesex_in_stacking')) agesex_in_stacking = F
if(!exists('remake_adult_prev_stackers')) remake_adult_prev_stackers = F
if(!exists('skipped_inla_rundate')) skipped_inla_rundate = run_date
if(!exists('use_adult_prev_stackers'))  use_adult_prev_stackers = FALSE

message(paste0('interacting gp 1 effects: ', interacting_gp_1_effects))
message(paste0('interacting gp 2 effects: ', interacting_gp_2_effects))



## reload data, fitting, and covariates to prepare for prediction
load(paste0(<<<< FILEPATH REDACTED >>>> '/model_image_history/', run_date, pathaddin, '.RData'))

if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)){
  message('Converting stackers to logit space')
  
  if(as.logical(agesex_in_stacking)) stacked_rasters <- cov_list[1:(length(z_list)*length(child_model_names)*length(sex_list))] else stacked_rasters <- cov_list
  
  if(use_gbd_covariate==T){
    stacked_rasters[['gbd_est']] <- gbd_raster
    remove(gbd_raster)
  }
  
  ## transform the rasters
  for (ii in 1:length(stacked_rasters)) { 
    
    ## Preserve variable names in the raster first
    tmp_rastvar <- names(stacked_rasters[[ii]])
    
    message(paste0('working on transforming covar ', tmp_rastvar))
    
    ## Logit
    stacked_rasters[[ii]] <- logit(stacked_rasters[[ii]])
    
    ## Reassign names
    names(stacked_rasters[[ii]]) <- tmp_rastvar
    rm(tmp_rastvar)
  }
  
  cov_list <- stacked_rasters
  remove(stacked_rasters)
}

#Load input data
input_data<-readRDS( ## save this here in case predict dies
  file = sprintf('<<<< FILEPATH REDACTED >>>>/output/%s/%s_data_input_list_%s_holdout_%i_agebin_%i.RDS',
                run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))

#Load mbg model
model_fit <- readRDS( file = sprintf('<<<< FILEPATH REDACTED >>>>/output/%s/%s_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                                  run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Predict ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tic("MBG - predict model") ## Start MBG - model predict timer

## Run predict_mbg on chunks of 50 samples (to avoid memory issues)
message('Making predictions in 50 draw chunks.')

max_chunk <- 50
samples   <- 1000

if(!as.logical(use_subnat_res)) simple_raster2 <- NULL

## Create vector of chunk sizes
chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)

outputdir_chunks <- file.path(outputdir, 'prediction_chunks')
dir.create(outputdir_chunks, showWarnings = FALSE)

#Set a seed that is unique to this chunk
seed = holdout + chunk*10 + as.numeric(gsub('_', '', run_date))/100000 + c(cssa = 500, 'essa_sdn-COM' = 700, sssa = 400, 'wssa-CPV-STP-MRT' = 700)[reg]


    pm<-  predict_mbg_tmb(samples              = 50,
                          seed                 = seed,
                          tmb_input_stack      = input_data,
                          model_fit_object     = model_fit,
                          fes                  = all_fixed_effects,
                          sr                   = simple_raster,
                          yl                   = year_list,
                          zl                   = z_list,
                          xl                   = sex_list,
                          covs_list            = cov_list, 
                          clamp_covs           = clamp_covs,
                          cov_constraints = covariate_constraint_vectorize(config),
                          use_space_only_gp = as.logical(use_space_only_gp),
                          use_time_only_gmrf = as.logical(use_time_only_gmrf),
                          use_age_only_gmrf = as.logical(use_age_only_gmrf),
                          use_sex_only_gmrf = as.logical(use_sex_only_gmrf),
                          int_gp_1_effect =            as.logical(use_gp),
                          use_sz_gp = as.logical(use_sz_gp),
                          use_sx_gp = as.logical(use_sx_gp),
                          use_tx_gp = as.logical(use_tx_gp),
                          use_zx_gp = as.logical(use_zx_gp),
                          use_cre_z_gp = as.logical(use_cre_z_gp),
                          use_cre_x_gp = as.logical(use_cre_x_gp),
                          use_gbd_covariate = as.logical(use_gbd_covariate),
                          int_gp_2_effect =   as.logical(any(interacting_gp_2_effects != '')),
                          #use_covs         = FALSE,
                          int_gp_1_effs   = interacting_gp_1_effects,
                          int_gp_2_effs   = interacting_gp_2_effects,
                          mesh_int = mesh_int,
                          mesh_s = mesh_s)  
    
    for(x in 1:length(sex_list)) {  
      for(z in 1:length(z_list)) {
        if (length(z_list)>1){
          if (length(sex_list)==1){
            pm_zx <- pm[[z]] 
          } else if (length(sex_list)>1)
            pm_zx <- pm[[x]][[z]]
        }  else {
          pm_zx<- pm
        }
        
        pathaddin <- paste0('_bin', z_list[z], '_sex', sex_list[x], '_chunk',chunk, '_', reg,'_',holdout)
        message(paste0('Saving draws for chunk ', chunk, ' for agebin ', z_list[z], ' sex ', sex_list[x]))
        save(
          pm_zx,
          file     = paste0(outputdir_chunks, '/', indicator,'_cell_draws_eb',pathaddin,'.RData'),
          compress = TRUE
        )
        
        rm(pm_zx)
        gc()
      }
    }
    rm(pm)
    gc()


# Write a an empty file to indicate done with this parallel script
write(NULL, file = paste0(outputdir, "/prediction_chunks/fin_", pathaddin_chunk))

message(paste("Done with chunk prediction for", reg))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

