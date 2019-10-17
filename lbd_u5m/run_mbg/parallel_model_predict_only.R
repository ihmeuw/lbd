#####################################################################
## Generic parallel script for predict only u5m
## RB

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## grab arguments
## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

## make a pathaddin that get used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

## load an image of the main environment
load('<<<< FILEPATH REDACTED >>>>')

## In case anything got overwritten in the load, reload args
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
outputdir <- file.path('<<<< FILEPATH REDACTED >>>>')
dir.create(outputdir, showWarnings = FALSE)
# Check that a small file has written to the model fit directory
#  This indicates that the model has properly fit
check_model_fit_fp <- paste0(outputdir,'/model_fit_progress/fin_',pathaddin)
if(!file.exists(check_model_fit_fp)){
  stop("ISSUE: Model has not fit correctly. Please re-run TMB fitting to proceed.")
}


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
mbg_setup(package_list = c(package_list, 'matrixcalc'), repos = core_repo)

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
check_config()

# We need to be in the singularity image, and specifically the LBD one if using TMB
if(!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if(as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Print the core_repo hash and check it
message("Printing git hash for 'core_repo' and checking against LBD Core Code master repo")
record_git_status(core_repo = core_repo, check_core_repo = TRUE)

## Make sure this inla patch is implemented if running on geos
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

## cores to use
cores_to_use <- round(as.numeric(slots)*.5)

## some set up
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (class(z_list)    == "character") z_list    <- eval(parse(text=z_list))

## Load simple polygon template to model over
gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
                                           shapefile_version = modeling_shapefile_version)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']] # TODO DELETE LINE


## reload data an prepare for MBG
load('<<<< FILEPATH REDACTED >>>>')

input_data <- readRDS('<<<< FILEPATH REDACTED >>>>')
model_fit <- readRDS('<<<< FILEPATH REDACTED >>>>')

clamp_covs <- FALSE

## Create a temporary directory to save MBG outputs, just in case
preds_temp_dir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(preds_temp_dir, showWarnings=FALSE)


## ** COMMENTING OUT FOR NOW - BREAKS ON LARGE MATRICES **
## Force the joint precision matrix to be positive definite. This will only
##  be run in cases where the user has decided to override the default model
##  behavior (which is to fail if the joint precision matrix is not positive
##  definite) by writing a `fin__(job info)` file in the `model_fit_progress/`
##  folder.
# if( !matrixcalc::is.positive.definite(as.matrix(model_fit$sdrep$jointPrecision)) ){
#   message("WARNING: Joint Precision matrix is not positive definite.")
#   message("Finding closest positive definite matrix...\n")
#   model_fit$sdrep$jointPrecision <- Matrix(
#     Matrix::nearPD(model_fit$sdrep$jointPrecision)$mat,
#     sparse=TRUE
#   )
# }

# Define filepath where the list of cell pred chunks will be saved
cellpred_list_fp <- sprintf("%s/cp_list_%s.RDS",preds_temp_dir,pathaddin)
# If the cell pred chunks file exists, ALWAYS skip to age combination for now
skip_drawgen_if_done <- TRUE


if((skip_drawgen_if_done == TRUE) & (file.exists(cellpred_list_fp))){
  message("Draws have already been generated but not split out by age.")
  message("Skipping to the age group splitting step now.")
  max_chunk <- 50
  # U5M MANUAL OVERRIDE FOR NOW
  samples   <- 500

  ## Create vector of chunk sizes
  chunks <- rep(max_chunk, samples %/% max_chunk)
  if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
  chunks_idx <- as.list(1:length(chunks))

} else {  
  tic("MBG - predict model") ## Start MBG - model predict timer

  ## Run predict_mbg on chunks of 50 samples (to avoid memory issues)
  message('Making predictions in 50 draw chunks.')

  max_chunk <- 50
  samples   <- as.numeric(samples)

  ## Create vector of chunk sizes
  chunks <- rep(max_chunk, samples %/% max_chunk)
  if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
  chunks_idx <- as.list(1:length(chunks))

  # NOTE
  # We currently run this code in a for loop in order to save on memory
  # that is required for doing this in a lapply
  pm <- lapply(chunks_idx, function(x) NULL)

  for(i in 1:length(chunks_idx)){
    samp_idx <- chunks_idx[[i]]
    samp <- chunks[[samp_idx]]
    if(fit_with_tmb == FALSE){
      pm[[i]] <- predict_mbg(
        res_fit       = model_fit,
        cs_df         = cs_df,
        mesh_s        = mesh_s,
        mesh_t        = mesh_t,
        cov_list      = cov_list,
        samples       = samp,
        simple_raster = simple_raster,
        transform     = transform,
        coefs.sum1    = coefs_sum1,
        pred_gp       = as.logical(use_gp),
        shapefile_version = modeling_shapefile_version
      )[[3]]
    } else {
      pm[[i]] <- predict_mbg_tmb(
        samples              = samp,
        seed                 = NULL,
        tmb_input_stack      = input_data,
        model_fit_object     = model_fit,
        fes                  = all_fixed_effects, # TODO use input_data or model_fit object for this (in case its changed due to checks)
        sr                   = simple_raster,
        yl                   = year_list,
        zl                   = z_list,
        covs_list            = cov_list,
        clamp_covs           = clamp_covs # TODO ADD TO CONFIG
        # cov_constraints      = covariate_constraint_vectorize(config)
      )
    }
    # gc(full=TRUE)
    
  }

  # The command above can fail for large objects
  message("Memory has successfully been allocated for the full preds list")
  # Clear up some memory
  rm(list=c("input_data","model_fit"))
  gc(full=TRUE)

  # message(sprintf("Saving cell pred list to %s for now...",cellpred_list_fp))
  # saveRDS(pm, file=cellpred_list_fp)
  # message("  ...Finished saving cell pred list.")
  
  toc(log = T) # Stop MBG - model predict timer
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read cell pred back into memory if it does not already exist
if(not('pm' %in% ls())){
  message("Reading cell preds object from file...")
  pm <- readRDS(cellpred_list_fp)
  message("  ...Read from file successfully.")
}

# if z dimension has more than one level, then save each z as a different indicator
if(length(z_list) > 1){
  
  # reorder pm list, right now its z within each chunk. rbind all z's together
  message("Reordering z list in cell pred:")
  for(z in z_list){ # z_list must be integers starting with 1
    message(sprintf("  Age group %s",z))
    if(length(chunks) > 1){
      for(ch in 2:length(chunks)){
        message(sprintf("    cbinding chunk %s of %s to chunk 1...",ch,length(chunks)))
        pm[[1]][[z]] <- cbind(pm[[1]][[z]], pm[[ch]][[z]])
        pm[[ch]][[z]] <- integer(0)
      }
    }
  } 
  pm <- pm[[1]] # pm is now a list of cell_preds by z
      
  # loop over z and save as an indicator each one
  orig_indic  <- indicator
  orig_paddin <- pathaddin
  orig_outdir <- outputdir
  
  message('Wrapping up')
  
  for(z in z_list) {
    cptmp <- pm[[z]]
    
    indicator <- sprintf('%s_%s%i',orig_indic,zcol,z)  # new indicator name
    pathaddin <- paste0('_bin',z,'_',reg,'_',holdout) # new pathaddin
    outputdir <- '<<<< FILEPATH REDACTED >>>>' # new outputdir
    dir.create(outputdir)
    message(sprintf('New indicator: %s',indicator))
    
    # make a mean raster
    library(matrixStats)
    mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cptmp),ncol = max(period_map$period)))
    sd_ras    <- insertRaster(simple_raster,matrix(  rowSds(cptmp),ncol = max(period_map$period)))
    
    # save z specific objects
    writeRaster(
      mean_ras,
      file      = paste0(outputdir, '/', indicator,'_prediction_eb',pathaddin),
      overwrite = TRUE
    )
    
    save(
      cptmp,
      file     = paste0(outputdir, '/', indicator,'_cell_draws_eb',pathaddin,'.RData'),
      compress = TRUE
    )
    
    pdf(paste0(outputdir,'mean_raster', pathaddin, '.pdf'))
    plot(mean_ras,main='mean',maxpixel=1e6)
    plot(sd_ras,main='sd',maxpixel=1e6)
    dev.off()
    
    rm(cptmp)
    gc(full=TRUE)
  }
  
  indicator <- orig_indic
  pathaddin <- orig_paddin
  outputdir <- orig_outdir
  
  # save training data
  write.csv(
    df,
    file = (paste0(outputdir, '/', indicator,'_trainingdata',pathaddin,'.csv')),
    row.names = FALSE
  )
  
  message('done saving indicator-specific outputs by z')
      
}  else { # if no z colums (most peoples cases)
  
  ## THIS IS NOT THE CASE FOR U5M -- SKIP FOR NOW SO WE CAN DELETE THE 
  ##  model_fit OBJECT AT AN EARLIER STAGE
  stop("This section is not yet implemented for U5M.")
  # ## Make cell preds and a mean raster
  # cell_pred <- do.call(cbind, pm)
  # mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))
  # toc(log = T) # Stop MBG - model predict timer
    
  # message('Wrapping up')
  # save_mbg_preds(config     = config,
  #                time_stamp = time_stamp,
  #                run_date   = run_date,
  #                mean_ras   = mean_ras,
  #                sd_ras     = NULL,
  #                res_fit    = model_fit,
  #                cell_pred  = cell_pred,
  #                df         = df,
  #                pathaddin  = pathaddin)
  
  
  # # plot the mean raster
  # pdf(paste0(outputdir,'mean_rasterXX', pathaddin, '.pdf'))
  # plot(mean_ras,maxpixel=1e6)
  # dev.off()
}

# Write a an empty file to indicate done with this parallel script
write(NULL, file = paste0(outputdir, "/fin_", pathaddin))


## Clean up the temporary backup RDS files for the draw chunks
cleanup_files <- grep(
  pathaddin, 
  list.files(preds_temp_dir),
  value=TRUE
)
if(length(cleanup_files) > 0) file.remove(paste0(preds_temp_dir, cleanup_files))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
