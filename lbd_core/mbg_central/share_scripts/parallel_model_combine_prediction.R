#################################################################################
## Generic parallel model prediction chunk combine script                      ##
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

if(!as.logical(use_subnat_res)) simple_raster2 <- NULL

## Create vector of chunk sizes
chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)

outputdir_chunks <- file.path(outputdir, 'prediction_chunks')
dir.create(outputdir_chunks, showWarnings = FALSE)


  for (x in sex_list) {
    for (z in z_list) {
      cell_pred<-NULL
      message('Reload and combine temporary chunk files agebin ', z, ' sex ', x)
      for (chunk in 1:length(chunks)) {
        pathaddin <- paste0('_bin',z, '_sex', x, '_chunk',chunk, '_', reg,'_',holdout)
        load(paste0(outputdir_chunks, '/', indicator,'_cell_draws_eb',pathaddin,'.RData'))
        
        cell_pred <- cbind(cell_pred, pm_zx)
        
        rm(pm_zx)
        gc()
      }
      
      pathaddin <- paste0('_bin',z, '_sex', x, '_',reg,'_',holdout) # new pathaddin
      
      # make a mean raster
      library(matrixStats)
      mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))
      sd_ras    <- insertRaster(simple_raster,matrix(  rowSds(cell_pred),ncol = max(period_map$period)))
      
      # save z specific objects
      writeRaster(
        mean_ras,
        file      = paste0(outputdir, '/', indicator,'_prediction_eb',pathaddin),
        overwrite = TRUE
      )
      
      save(
        cell_pred,
        file = (paste0(outputdir, "/", indicator, "_cell_draws_eb", pathaddin, ".RData")),
        compress = TRUE
      )
      
      #Create mean & sd pdfs
      pdf(paste0(outputdir,'mean_raster', pathaddin, '.pdf'))
      plot(mean_ras,main='mean',maxpixel=1e6)
      plot(sd_ras,main='sd',maxpixel=1e6)
      dev.off()
      
      #Clean up
      rm(cell_pred)
      rm(mean_ras)
      rm(sd_ras)
      gc()
      
    }
  }
  
  #Remove temporary chunk files
  message('Remove temporary chunk files')
  for(x in sex_list) {
    for (z in z_list) {
      for (chunk in 1:length(chunks)) {
        pathaddin <- paste0('_bin',z,'_sex', x, '_chunk', chunk, '_', reg,'_',holdout)
        file.remove(paste0(outputdir_chunks, '/', indicator,'_cell_draws_eb',pathaddin,'.RData'))
      }
    }
  }

message("wrapping up")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Save log of config file
write.csv(config, paste0(outputdir, "/config.csv"), row.names = FALSE)

## timer stuff
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
output_file <- paste0(outputdir, "run_summary_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}


# Write a an empty file to indicate done with this parallel script
write(NULL, file = paste0(outputdir, "/fin_", pathaddin))

## Write CSV
write.csv(df_timer,file = output_file, row.names = FALSE)


message(paste("Done with model prediction for", reg))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

