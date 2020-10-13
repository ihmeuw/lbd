in_dir <- commandArgs()[4]
out_dir <- commandArgs()[5]
yr_max <- commandArgs()[6]

source('FILEPATH.R')
package_list <- c('seegSDM','dismo', 'raster', 'maptools', 'sp', 'rgdal', 'ggplot2', 'rgeos')
load_R_packages(package_list)
source('FILEPATH.R')
source('FILEPATH.R')

# read in covariates and data
cov_name_list <- read.csv(paste0(in_dir, '/FILEPATH.csv'), as.is=T)
cov_names <- cov_name_list$cov_name

covs <- stack(paste0(in_dir, '/FILEPATH', yr_max, '.grd'))
names(covs) <- cov_names

## SUMMARY CALCULATIONS

# make lists of all files in output directories
data_files <- list.files(path = paste0(out_dir, '/FILEPATH'), full.names = TRUE)
dat_all <- c()
stats_files <- list.files(path = paste0(out_dir, '/FILEPATH'), full.names = TRUE)
stats_list <- c()
model_files <- list.files(path = paste0(out_dir, '/FILEPATH'), full.names = TRUE)
model_list <- c()
for (i in 1:length(data_files)){
  model_list[[i]] <- get(load(model_files[[i]]))
  stats_list[[i]] <- get(load(stats_files[[i]]))
  dat_all <- rbind(dat_all, read.csv(data_files[[i]]))
}

# write data files - all bootstrap data
dat_all.bg <- dat_all[dat_all$PA==0,]
write.csv(dat_all.bg, paste0(out_dir, '/FILEPATH.csv'), row.names=F)
dat_all.occ <- dat_all[dat_all$PA==1,]
write.csv(dat_all.occ, paste0(out_dir, '/FILEPATH.csv'), row.names=F)

# get prediction raster and save it as GeoTiff
preds_files <- list.files(path = paste0(out_dir, '/FILEPATH'), full.names = TRUE)
preds_list <- lapply(preds_files, raster)
preds_list <- stack(preds_list)

preds_sry <- combinePreds(preds_list)
names(preds_sry) <- c('mean', 'median', 'lowerCI', 'upperCI')

# get synoptic fit statistics for summarized model...
# ...mean
syn_stats <- get_fit_stats(dat_all, preds_sry$mean, get_ROC=T)
fit_stats <- syn_stats$fit_stats
fpr <- syn_stats$fpr
sens <- syn_stats$sensitivity
opt_thresh <- syn_stats$opt_thresh
write.csv(data.frame(opt_thresh), paste0(out_dir, '/FILEPATH.csv'))
write.csv(fit_stats, paste0(out_dir, '/FILEPATH.csv'))
write.csv(cbind(fpr, sens), paste0(out_dir, '/FILEPATH.csv'))
# ...upper CI
syn_stats_uppCI <- get_fit_stats(dat_all, preds_sry$upperCI)
write.csv(syn_stats_uppCI, paste0(out_dir, '/FILEPATH.csv'))
#... lower CI
syn_stats_lowCI <- get_fit_stats(dat_all, preds_sry$lowerCI)
write.csv(syn_stats_lowCI, paste0(out_dir, '/FILEPATH.csv'))

# calculate uncertainty in the predictions 
preds_sry$uncertainty <- preds_sry[[4]] - preds_sry[[3]]

# classify predictions based on optimal threshold of mean predictions
preds_sry$binary <- preds_sry$mean
values(preds_sry$binary)[values(preds_sry$binary) >= opt_thresh] <- 1
values(preds_sry$binary)[values(preds_sry$binary) < opt_thresh] <- 0

# write mean prediction, binary prediciton, and uncertainty rasters as GeoTiffs
writeRaster(preds_sry$mean, file = paste0(out_dir, '/FILEPATH.tif'), format = 'GTiff', overwrite=T)
writeRaster(preds_sry$uncertainty, file = paste0(out_dir, '/FILEPATH.tif'), format = 'GTiff', overwrite=T)
writeRaster(preds_sry$binary, file = paste0(out_dir, '/FILEPATH.tif'), format = 'GTiff', overwrite=T)

writeRaster(preds_sry,
            file = paste0(out_dir, '/FILEPATH.tif'),
            format = 'GTiff',
            overwrite=T)

write.csv(names(preds_sry), paste0(out_dir, '/FILEPATH.csv'))

# save matrix of summary statistics across all bootstraps
stats <- do.call('rbind', stats_list)
write.csv(stats, file = paste0(out_dir, '/FILEPATH.csv'))

# get and save relative influence scores of covariates; make a box plot
relinf <- getRelInf(model_list, plot=F) 
write.csv(relinf, file = paste0(out_dir, '/FILEPATH.csv'))

# get MESS raster for preds_sry
MESS <- mess(preds_sry, dat_all.occ[,c('longitude', 'latitude')])
writeRaster(x=MESS, filename=paste0(out_dir, '/FILEPATH'), format="GTiff", overwrite=T)

## COMPUTE OoS BIN:
# code for computing the probability of a single pixel w/in a given IU has an...
# ...environmental suitability value GT a threshold found by optimizing on an ROC curve

probGT <- function(vals, thresh){
  dens <- density(vals)
  if (thresh*max(dens$x) < min(dens$x)){
    l.lim <- min(dens$x)
  } else {
    l.lim <- thresh*max(dens$x)
  }
  return(integrate.xy(dens$x, dens$y, l.lim, max(dens$x)))
}

ius <- shapefile('FILEPATH.shp')

for(i in 1:length(index_list)){
  bin_dir <- paste0(out_dir, '/FILEPATH')
  if (!dir.exists(bin_dir)) dir.create(bin_dir)

  mean_preds <- preds_sry$mean
  lowCI <- preds_sry$lowerCI
  uppCI <- preds_sry$upperCI
  # load model threshold
  thresh <- read.csv(paste0(in_dir, '//FILEPATH.csv'))
  
  ## COMPUTE MEAN, UPPER CI AND LOWER CI PROBABILITIES AT THE OPTIMIZED THRESHOLD VALUE
  
  # compute probs for mean preds
  probs <- c()
  for (i in 1:nrow(ius@data)){
    vals <- extract(x = mean_preds, y = ius[i,])[[1]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0){
      probs <- append(probs, NA)
    } else if (length(vals) < 2){
      if (vals[[1]] >= opt_thresh){
        probs <- append(probs, 1)
      } else {
        probs <- append(probs, 0)
      }
    } else {
      probs <- append(probs, probGT(vals, opt_thresh))
    }
  }
  write.csv(data.frame(probs), paste0(bin_dir, '/FILEPATH', toString(opt_thresh*100), '.csv'), row.names = F)
  
  # compute probs for upper CI preds
  uppCI_probs <- c()
  for (i in 1:nrow(ius@data)){
    vals <- extract(x = uppCI, y = ius[i,])[[1]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0){
      uppCI_probs <- append(uppCI_probs, NA)
    } else if (length(vals) < 2){
      if (vals[[1]] >= opt_thresh){
        uppCI_probs <- append(uppCI_probs, 1)
      } else {
        uppCI_probs <- append(uppCI_probs, 0)
      }
    } else {
      uppCI_probs <- append(uppCI_probs, probGT(vals, opt_thresh))
    }
  }
  write.csv(data.frame(uppCI_probs), paste0(bin_dir, '/FILEPATH', toString(opt_thresh*100), '.csv'), row.names = F)
  
  # compute probs for lower CI preds
  lowCI_probs <- c()
  for (i in 1:nrow(ius@data)){
    vals <- extract(x = lowCI, y = ius[i,])[[1]]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0){
      lowCI_probs <- append(lowCI_probs, NA)
    } else if (length(vals) < 2){
      if (vals[[1]] >= opt_thresh){
        lowCI_probs <- append(lowCI_probs, 1)
      } else {
        lowCI_probs <- append(lowCI_probs, 0)
      }
    } else {
      lowCI_probs <- append(lowCI_probs, probGT(vals, opt_thresh))
    }
  }
  write.csv(data.frame(lowCI_probs), paste0(bin_dir, '/FILEPATH', toString(opt_thresh*100), '.csv'), row.names = F)
  
  endm <- shapefile('FILEPATH.shp')
  ius <- cbind(endm$IU_ID, probs, lowCI_probs, uppCI_probs)
  write.csv(ius, paste0(bin_dir, '/FILEPATH', toString(opt_thresh*100), '.csv'), row.names = F)
  print(paste0('...is done'))
}


## PLOTTING

source('\FILEPATH.R')

# PLOT OCCURRENCE AND BACKGROUND HISTOGRAMS

get_pred_hist(dat0.preds, out_dir=out_dir, TRUE, file_name="/FILEPATH")
get_pred_hist(dat1.preds, out_dir=out_dir, file_name="/FILEPATH")

## PLOT EFFECT CURVES

# get the covariate mean effects; make a plot of each mean effect curve
effect <- getEffectPlots(model_list, plot=T)

# get the order of plots (all except relinf)! - customize list to covs of interest
order <- match(rownames(relinf), cov_names)

#covariate histogram data
dat.pts <- dat_all[dat_all$PA==1, cov_names[order]]
get_effect_hist_plot(dat.pts, effect, order, out_dir=out_dir, file_name="/FILEPATH")

## RASTER DATA PLOTTING

get_raster_grid(covs, order, extent=ext, out_dir=out_dir, file_name="/FILEPATH")

# PLOT SUITABILITY MAPS

# choose disease extent to display manually
expmt.pts <- rbind(dat0.pts, dat1.pts)
zoom <- c(abs(min(expmt.pts$long))/180 - 0.20,
          abs(max(expmt.pts$long))/180 + 0.10,
          abs(min(expmt.pts$lat))/180 - 0.1,
          abs(max(expmt.pts$lat))/180 + 0.27)

ext <- zoom*extent(preds_sry)
legend_position <- c(0.8, 0.03)  

# choose legend_position manually
# legend_position <- c(0.83, 0.03)

blue_red <- c("#313695", "#ffffbf", "#d73027")
green_purple_bin <- c("#CBE0CF", "#9b0065")

# mean pred map
get_raster_map(preds_sry$mean,
               out_dir=out_dir,
               crop=FALSE,
               clip=FALSE,
               extent=ext,
               legend_position=legend_position,
               title=paste0("Env. Suitability"),
               file_name="/FILEPATH")

# mean pred map - occ points
get_raster_map(preds_sry$mean,
               out_dir=out_dir,
               occ_pts = dat1.pts,
               crop=FALSE,
               clip=FALSE,
               extent=ext,
               legend_position=legend_position,
               title=paste0("Env. Suitability"),
               file_name="/FILEPATH")

# mean pred map - bg points
get_raster_map(preds_sry$mean,
               out_dir=out_dir,
               bg_pts = dat0.pts,
               crop=FALSE,
               clip=FALSE,
               extent=ext,
               legend_position=legend_position,
               title=paste0("Env. Suitability"),
               file_name="/FILEPATH")

# uncertainty map
get_raster_map(preds_sry$uncertainty, 
               color_scheme=blue_red, 
               out_dir=out_dir, 
               crop=FALSE,
               clip=FALSE,
               extent=ext,
               legend_position=legend_position,
               title=paste0("Uncertainty"), 
               file_name="/FILEPATH")

# threshold-classified mean pred map
get_raster_map(preds_sry$binary, 
               color_scheme=green_purple_bin, 
               out_dir=out_dir, 
               crop=FALSE,
               clip=FALSE,
               extent=ext,
               legend_position=legend_position,
               title=paste0("Threshold Env. Suitability"),
               file_name="/FILEPATH")



