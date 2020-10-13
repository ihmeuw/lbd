source('FILEPATH.R')
package_list <- c('seegSDM')
load_R_packages(package_list)

yr <- as.numeric(commandArgs()[4])
in_dir <- commandArgs()[5]
monthly <- as.logical(commandArgs()[6])

cov_dir <- "ADDRESS"
cov_name_list <- read.csv(paste0(in_dir, '/FILEPATH.csv'), as.is=T)
out_dir <- paste0(in_dir, '/FILEPATH')
if(!dir.exists(out_dir)) dir.create(out_dir)

# make/write covariate stack(s)
if (monthly) {
  for (mo in 1:12) {
    cov_stack <- stack()
    if (mo < 10) mo <- paste0(0, mo)
    if (!file.exists(paste0(out_dir, '/pred_stack_', yr, '_', mo, '.grd'))) {
      for (i in 1:nrow(cov_name_list)) {
        cov <- cov_name_list$cov_name[i]
        measure <- cov_name_list$measure[i]
        static <- as.logical(cov_name_list$static[i])
        # .tif file path
        if (static) {
          f_dir <- paste0(cov, '/', measure, '/synoptic/', cov, '_', measure, '_synoptic.tif')
        } else {
          f_dir <- paste0(cov, '/', measure, '/1m/', cov, '_', measure, '_1m_', yr, '_', mo, '_00.tif')
        }
        cov_stack <- stack(cov_stack, raster(paste0(cov_dir, f_dir)))
      }
      writeRaster(cov_stack, paste0(out_dir, '/pred_stack_', yr, '_', mo, '.grd'))
    }
  }
} else {
  cov_stack <- stack()
  if (!file.exists(paste0(out_dir, 'pred_stack_', yr, '.grd'))){
    for (i in 1:nrow(cov_name_list)) {
      cov <- cov_name_list$cov_name[i]
      measure <- cov_name_list$measure_yrly[i]
      static <- as.logical(cov_name_list$static[i])
      # .tif file path
      if (static) {
        f_dir <- paste0(cov, '/', measure, '/synoptic/', cov, '_', measure, '_synoptic.tif')
      } else {
        f_dir <- paste0(cov, '/', measure, '/1y/', cov, '_', measure, '_1y_', yr, '_00_00.tif')
      }
      cov_stack <- stack(cov_stack, raster(paste0(cov_dir, f_dir)))
    }
    writeRaster(cov_stack, paste0(out_dir, '/pred_stack_', yr, '.grd'))
  }
}

