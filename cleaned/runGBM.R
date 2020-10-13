jobnum <- commandArgs()[4]
opt_type <- commandArgs()[5]
in_dir <- commandArgs()[6]
out_dir <- commandArgs()[7]
data_file <- commandArgs()[8]
stack_out_dir <- commandArgs()[9]
dat_dir <- commandArgs()[10]
index <- commandArgs()[11]
cov_dir <- commandArgs()[12]

source('FILEPATH.R')
source('FILEPATH.R')
source('FILEPATH.R')
source('FILEPATH.R')

source('FILEPATH.R')
package_list <- c('seegSDM', 'gbm', 'dismo', 'maptools', 'rgdal')
load_R_packages(package_list)

## CONSTRUCT DATA FOR BRT STEP
cov_name_list <- read.csv(paste0(in_dir, '/FILEPATH.csv'), as.is=T)
cov_names <- cov_name_list$cov_name

if (!file.exists(paste0(out_dir, '/FILEPATH', jobnum, '.csv'))){
  # read in data
  dat.orig <- read.csv(paste0(dat_dir, '/', data_file), as.is=T)
  
  covs <- brick('FILEPATH.grd')
  # make sure bias grid exists
  bias_fname <- paste0(in_dir, '/FILEPATH.tif')
  if (!file.exists(bias_fname)) {
    ref <- raster('FILEPATH.tif')
    IU_shapefile <- shapefile('FILEPATH.shp')
    base <- shapefile("FILEPATH.shp")
    field <- 'Endem_MDA'
    exclude <- c("Endemic (under MDA)", "Endemic (MDA not started)", "Endemic, Post-MDA Surveillance")
    bias_grid_final <- make_bias_grid(IU_shapefile, base, field, exclude, ref)
    writeRaster(bias_grid_final, filename = bias_fname)
  } else {
    bias_grid_final <- raster(bias_fname)
  }
  
  # fill in missing years in data for annual covariate sampling
  year <- na.omit(dat.orig$year)
  yr_freqs <- data.frame(table(year))
  yr_freqs$prob <- yr_freqs$Freq / sum(yr_freqs$Freq)
  for (r in 1:nrow(dat.orig)){
    if (is.na(dat.orig$year[r])) {
      s <- sample(yr_freqs$year, size=1, prob=yr_freqs$prob)
      dat.orig$year[r] <- as.numeric(levels(s))[s]
    }
  }
  
  ## get lat-longs
  dat.pts <- dat.orig[!(dat.orig$MDA_status == 0 & dat.orig$presence == 0),]
  dat.pts <- dat.pts[,c("long", "lat", "year")]
  dat.pts <- na.omit(dat.pts)
  names(dat.pts) <- c("long", "lat", "year")
  dat.pts <- cbind(PA = rep(1, nrow(dat.pts)), dat.pts)
  
  pts_buffer <- shapefile(paste0(in_dir, '/FILEPATH', index, '.shp'))
  
  bg <- data.frame(bgSample(n = nrow(dat.pts),
                            raster = mask(bias_grid_final, pts_buffer),
                            prob=T))
  bg_year <- c()
  for (r in 1:nrow(bg)){
    s <- sample(yr_freqs$year, size=1, prob=yr_freqs$prob)
    bg_year[r] <- as.numeric(levels(s))[s]
  }
  bg$year <- bg_year
  
  names(bg) <- c('long', 'lat', 'year')
  
  dat <- rbind(dat.pts, cbind(PA = rep(0, nrow(bg)), bg))
  
  ## get lat-longs and cov values for polygon data
  
  polys <- dat.orig[dat.orig$shapefile != "",]
  polys <- polys[polys$presence == 1,]
  
  path_converter <- function(winpath, disc){
    if (disc == "ADDRESS"){
      path <- gsub("FILEPATH", "FILEPATH", gsub("\\\\", "/", winpath))
    }
    return(path)
  }
  
  smp.pts <- data.frame()
  smp.bg <- data.frame()
  for (i in 1:nrow(polys)){
    print(i)
    shapeDF <- readOGR(path_converter(polys$shapefile[i], "J"))
    currGaul <- polys$poly_id[i]
    if ("GAUL_CODE" %in% colnames(shapeDF@data)) {
      if (currGaul == 0) {
        subs <- subset(shapeDF, (is.na(GAUL_CODE) | GAUL_CODE == 0))
      } else {
        subs <- subset(shapeDF, GAUL_CODE == currGaul)
      }
    } else {
      subs <- shapeDF
    }
    if(is.na(polys$villages_inside_polygon[i])) {
      polys$villages_inside_polygon[i] <- 1
    } else if(polys$villages_inside_polygon[i] == ""){
      polys$villages_inside_polygon[i] <- 1
    }
    pts <- data.frame(spsample(subs, polys$villages_inside_polygon[i], "random", iter=10))
    pts$year <- rep(polys$year[i], polys$villages_inside_polygon[i])
    crs <- CRS(proj4string(subs))
    bg_buffer <- spTransform(buffer(spTransform(subs, CRS=CRS("+proj=merc +ellps=GRS80")), 1e+05), CRS=crs)
    bg <- data.frame(bgSample(n = polys$villages_inside_polygon[i],
                              raster = mask(bias_grid_final, bg_buffer, updatevalue=0),
                              replace=F, prob=T))
    names(pts) <- c("long", "lat", "year")
    names(bg) <- c("long", "lat")
    smp.pts <- rbind(smp.pts, pts)
    smp.bg <- rbind(smp.bg, bg)
  }
  
  smp.pts <- cbind(PA = rep(1, nrow(smp.pts)), smp.pts)
  smp.bg <- cbind(PA = rep(0, nrow(smp.bg)), smp.bg)
  bg_year <- c()
  for (r in 1:nrow(smp.bg)){
    s <- sample(yr_freqs$year, size=1, prob=yr_freqs$prob)
    bg_year[r] <- as.numeric(levels(s))[s]
  }
  smp.bg$year <- bg_year
  
  dat_all <- rbind(dat, rbind(smp.pts, smp.bg))
  
  dat_all <- dat.all(dat_all, buf.dat=data.frame(), pt.dat=data.frame(), bg.dat=data.frame(), cov_name_list, cov_dir, stack_out_dir)
  
  write.csv(dat_all, paste0(out_dir, '/FILEPATH', jobnum, '.csv'), row.names=F)
} else {
  dat_all <- read.csv(paste0(out_dir, '/FILEPATH', jobnum, '.csv'), as.is=T)
}



## MAKE .CSV OF OPTIMIZED HPARS

# split to make OoS hold-out data
indx <- sample(nrow(dat_all), 0.80*nrow(dat_all))
data_train <- dat_all[indx,]
data_test <- dat_all[-indx,]
write.csv(data_train, file = paste0(out_dir, '/FILEPATH', jobnum, '.csv'), row.names=F)

run_optimizerPy(python <- 'FILEPATH',
                funcs.file_path <- 'FILEPATH',
                funcs.file <- '/FILEPATH.py',
                bounds.file_path <- in_dir,
                bounds.file <- '/FILEPATH.csv',
                data.file_path <- out_dir,
                data.loc <- paste0('/FILEPATH', jobnum, '.csv'),
                optimizer <- paste0(opt_type),
                learner <- 'brtR',
                cv_folds <- '10',
                n_calls <- '150',
                jobnum <- jobnum,
                col_start <- toString(which(names(data_train) == cov_names[1])))

# delete extraneous file
file.remove(paste0(out_dir, '/FILEPATH', jobnum, '.csv'))

# read in optimized hyperparameters
par <- read.csv(paste0(out_dir, '/FILEPATH', jobnum, '.csv'))

# get performance for selected hyperparameter values on training set
model_train <- run_brt_model(data_train, par, covs, cov_names, final=F)
train_stats <- model_train$stats
save(train_stats, file = paste0(out_dir,'/FILEPATH', jobnum,'.Rdata'))

# get performance for selected hyperparameter values on hold out set
model_test <- run_brt_model(data_test, par, covs, cov_names, final=F)
test_stats <- model_test$stats
save(test_stats, file = paste0(out_dir,'/FILEPATH', jobnum,'.Rdata'))

# get model outputs and performance for selected hyperparameter values on whole set
model <- run_brt_model(dat_all, par, covs, cov_names, final=T)
model_list <- model$model.list
stats <- model$stats
pred_rast <- model$pred.raster

save(model_list, file = paste0(out_dir,'/FILEPATH', jobnum,'.Rdata'))
save(stats, file = paste0(out_dir,'/FILEPATH', jobnum,'.Rdata'))
writeRaster(pred_rast,
            paste0(out_dir, '/FILEPATH', jobnum, '.tif'),
            format = 'GTiff', overwrite=T)

# get ROC-optimimal env. suit. threshold for this bootstrap
opt_thresh <- get_fit_stats(dat_all, pred_rast)$opt_thresh
if(!dir.exists(paste0(out_dir, '/FILEPATH'))) dir.create(paste0(out_dir, '/FILEPATH'))
write.csv(opt_thresh, paste0(out_dir, '/FILEPATH', jobnum, '.csv'))
  

