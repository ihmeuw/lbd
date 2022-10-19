#####################################################################
## Run predict_mbg() in parallel
#####################################################################

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source(<<<< FILEPATH REDACTED >>>>)
load_from_parallelize()
message(paste0("Using ", core_repo))

user <- Sys.info()[["user"]]

## Set repo locations
indic_repo <- <<<< FILEPATH REDACTED >>>>

## Load central libraries, packages, and miscellaneous MBG project functions.
commondir <- <<<< FILEPATH REDACTED >>>>
package_list <- c(t(read.csv(<<<< FILEPATH REDACTED >>>>, header = FALSE)))

message("Loading in required R packages and MBG functions")
source(<<<< FILEPATH REDACTED >>>>)

mbg_setup(package_list = package_list, repos = core_repo)

library(fasterize)
library(matrixStats)
library(sf)

## Focal 3 specific workflow: Pull in custom scripts.
setwd(indic_repo)
for (funk in list.files(<<<< FILEPATH REDACTED >>>>, recursive = TRUE)) {
  message(funk)
  source(<<<< FILEPATH REDACTED >>>>)
}

rownumber <- as.numeric(commandArgs()[5])
message(paste0("\nrownumber: ", rownumber, "\n"))

config <- set_up_config_focal_3(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator, config_file = paste0(indic_repo, "/config/", config_file, ".csv"), run_tests = FALSE)

year_list <- eval(parse(text = year_list))
predict_years <- eval(parse(text = predict_years))

if (!exists("restart_predict_years")) {
  restart_predict_years <- NULL
} else {
  restart_predict_years <- eval(parse(text = restart_predict_years))
}

if (!is.null(restart_predict_years)) {
  predict_years <- eval(restart_predict_years)[rownumber]
} else {
  predict_years <- predict_years[rownumber]
}

source(<<<< FILEPATH REDACTED >>>>)

reg <- region

message(paste0("\npredict_years: ", predict_years, "\n"))
message(paste0("\nperiod_map_predict: ", period_map_predict, "\n"))

outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)
dir.create(outputdir, showWarnings = FALSE)

## make a pathaddin that gets used widely
pathaddin <- paste0('_bin', 0, '_', region, '_', holdout)

# We need to be in the singularity image, and specifically the LBD one if using TMB
if (!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if (as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Make sure this inla patch is implemented if running on geos
INLA:::inla.dynload.workaround()

## cores to use
cores_to_use <- Sys.getenv("SGE_HGR_fthread")
message(paste("Model set to use", cores_to_use, "cores"))

# print out session info so we have it on record
sessionInfo()

### Load saved draws object
run_dir <- <<<< FILEPATH REDACTED >>>>
message(paste0("Loading saved draws object"))
load(<<<< FILEPATH REDACTED >>>>)

samples <- as.numeric(samples)

# Run predict_mbg on chunks of X samples (to avoid memory issues)
if (exists("max_chunk")) {
  max_chunk <- as.integer(max_chunk)
} else {
  max_chunk <- 50
}

for (l in 1:length(cov_list)) {
  message(sprintf("On cov %i out of %i", l, length(cov_list)))
  cov_list[[l]] <- raster::crop(cov_list[[l]], raster::extent(simple_raster))
  cov_list[[l]] <- raster::setExtent(cov_list[[l]], simple_raster)
  cov_list[[l]] <- raster::mask(cov_list[[l]], simple_raster)
}

## Create vector of chunk sizes
chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
c <- 0
pm <- lapply(chunks, function(samp) {
  c <<- c + 1
  message(paste0("Starting prediction batch ", c, " of ", length(chunks), "..."))
  if (fit_with_tmb == FALSE) {
    draws_start <- sum(chunks[1:c]) - chunks[c] + 1
    draws_end <- sum(chunks[1:c])
    draws_current <- draws[draws_start:draws_end]
    predict_mbg(res_fit       = model_fit,
                cs_df         = cs_df,
                mesh_s        = mesh_s,
                mesh_t        = mesh_t,
                cov_list      = cov_list,
                samples       = samp,
                region        = region,
                simple_raster = simple_raster,
                transform     = transform,
                coefs.sum1    = coefs_sum1,
                pred_gp       = as.logical(use_gp),
                shapefile_version = modeling_shapefile_version,
                predict_years = predict_years,
                return_mean_sd = FALSE,
                draws = draws_current,
                use_space_only_gp = use_space_only_gp,
                use_time_only_gmrf = use_time_only_gmrf,
                use_timebyctry_res = as.logical(use_timebyctry_res),
                predict_diagnostic = predict_diagnostic,
                predict_age_start = predict_age_start,
                predict_age_end = predict_age_end,
                simple_raster_subnats = simple_raster2,
                rw1_raw_covar_list = rw1_raw_covar_list,
                subnat_country_to_get = subnat_country_to_get)
  } else {
    predict_mbg_tmb(samples              = samp,
                    seed                 = NULL,
                    tmb_input_stack      = input_data,
                    model_fit_object     = model_fit,
                    fes                  = all_fixed_effects, # TODO use input_data or model_fit object for this (in case its changed due to checks)
                    sr                   = simple_raster,
                    yl                   = year_list,
                    zl                   = z_list,
                    covs_list            = cov_list,
                    clamp_covs           = clamp_covs,
                    cov_constraints = covariate_constraint_vectorize(config),
                    use_full_interacting_effect = as.logical(use_gp),
                    use_space_only_gp = as.logical(use_space_only_gp),
                    use_time_only_gmrf = as.logical(use_time_only_gmrf),
                    use_age_only_gmrf = as.logical(use_age_only_gmrf),
                    coefs.sum1           = coefs_sum)
  }
})

# Save predictions
run_dir <- <<<< FILEPATH REDACTED >>>>
message(paste0("Now saving predictions for years ", min(unlist(predict_years)), "-", max(unlist(predict_years))))
dir.create(run_dir, showWarnings = FALSE)

if(fit_with_tmb == FALSE){
  preds <- vector("list", length(pm))
  
  outputs <- vector("list", length(pm[[1]]))
  
  for (i in 1:length(pm)) {
    preds[[i]] <- pm[[i]]$cell_pred
    
    for (j in 1:length(outputs)) {
      if (!is.null(pm[[i]][[j]])) {
        if (length(pm[[i]][[j]]) == 1) {
        } else {
          outputs[[j]][[i]] <- pm[[i]][[j]]
          names(outputs)[[j]] <- names(pm[[i]])[[j]]
        }
      }
    }
  }
}

pm <- preds
cell_pred <- do.call(cbind, pm)

message(paste0("Saving cell_preds for ", predict_years))
save(cell_pred, <<<< FILEPATH REDACTED >>>>, compress = TRUE)

# cell_pred <- do.call(cbind, pm)
mean_ras <- insertRaster(simple_raster, matrix(rowMeans(cell_pred), ncol = length(predict_years)))
median_ras <- insertRaster(simple_raster, matrix(rowQuantiles(cell_pred, probs = 0.5), ncol = length(predict_years)))

message("Wrapping up")

## Plot the mean raster
cols <- c("#810f7c", "#9ebcda", "#ffffcc", "#ffeda0", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#800026")
cut <- c(0, 0.01, 0.02, 0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 1)
cuts <- c(seq(0, .3, .005))
cuts2 <- c(seq(0, 1, .05))

names(mean_ras) <- paste0("mean_pred_year.", predict_years)
names(median_ras) <- paste0("median_pred_year.", predict_years)

outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)
pdf(<<<< FILEPATH REDACTED >>>>)
plot(mean_ras, maxpixel = 1e6, breaks = cut, col = cols)
dev.off()

outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)
pdf(<<<< FILEPATH REDACTED >>>>)
plot(median_ras, maxpixel = 1e6, breaks = cut, col = cols)
dev.off()

### Save remaining cell preds
for (i in 1:length(outputs)) {
  if (!is.null(outputs[[i]])) {
    
    cell_pred <- do.call(cbind, cell_pred)
    mean_ras <- insertRaster(simple_raster, matrix(rowMeans(cell_pred), ncol = length(predict_years)))
    median_ras <- insertRaster(simple_raster, matrix(rowQuantiles(cell_pred, probs = 0.5), ncol = length(predict_years)))
    
    message("Wrapping up")

    names(mean_ras) <- paste0("mean_pred_year.", predict_years)
    names(median_ras) <- paste0("median_pred_year.", predict_years)
    
    outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)
    pdf(<<<< FILEPATH REDACTED >>>>)
    plot(mean_ras, maxpixel = 1e6)
    dev.off()
    
    outputdir <- file.path(<<<< FILEPATH REDACTED >>>>)
    pdf(<<<< FILEPATH REDACTED >>>>)
    plot(median_ras, maxpixel = 1e6)
    dev.off()
    
    writeRaster(mean_ras, file = <<<< FILEPATH REDACTED >>>>, overwrite = TRUE, format = "GTiff")
    writeRaster(median_ras, file = <<<< FILEPATH REDACTED >>>>, overwrite = TRUE, format = "GTiff")
  }
}

### Now save summary prediction rasters
pm <- preds
cell_pred <- do.call(cbind, pm)

s <- paste0(0, "_", holdout)

message("Starting mean raster summary")
mean_raster <- insertRaster(simple_raster, matrix(rowMeans(cell_pred), ncol = length(predict_years)))

message("Saving mean raster summary")
writeRaster(mean_raster, file = <<<< FILEPATH REDACTED >>>>, format='GTiff', overwrite = TRUE)

message("Starting upper 95 raster summary")
upper_95_raster <- insertRaster(simple_raster, matrix(rowQuantiles(cell_pred, probs = 0.975), ncol = length(predict_years)))

message("Saving upper 95 raster summary")
writeRaster(upper_95_raster, file = <<<< FILEPATH REDACTED >>>>, format='GTiff', overwrite = TRUE)

message("Starting lower 95 raster summary")
lower_95_raster <- insertRaster(simple_raster, matrix(rowQuantiles(cell_pred, probs = 0.025), ncol = length(predict_years)))

message("Saving lower 95 raster summary")
writeRaster(lower_95_raster, file = <<<< FILEPATH REDACTED >>>>, format='GTiff', overwrite = TRUE)

message("Starting range raster summary")
range_raster <- upper_95_raster - lower_95_raster

message("Saving range raster summary")
writeRaster(range_raster, file = <<<< FILEPATH REDACTED >>>>, format='GTiff', overwrite = TRUE)
