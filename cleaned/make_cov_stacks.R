invisible(source('FILEPATH.R'))
source('FILEPATH.R')
package_list <- c('seegSDM')
load_R_packages(package_list)

yr_min <- as.numeric(commandArgs()[4])
yr_max <- as.numeric(commandArgs()[5])
cov_name <- toString(commandArgs()[6])
measure <- toString(commandArgs()[7])
monthly <- as.logical(commandArgs()[8])
cov_dir <- toString(commandArgs()[9])
stack_out_dir <- toString(commandArgs()[10])
months_lag <- as.numeric(commandArgs()[11])
title <- toString(commandArgs()[12])
mo_min <- as.numeric(commandArgs()[13])

all_yrs <- c(yr_min:yr_max)

cov_stack <- stack()

# make/write covariate stacks for each year

if (monthly){
  for (yr_orig in all_yrs){
    stack_name <- paste0(stack_out_dir, title, '_monthly_', yr_orig, '.grd')
    if (!file.exists(stack_name)){
      # in bootstrap, fitting to distribution of months => to impute monthly data, need all monthly covariates
      if (yr_orig!=min(all_yrs)) mo_min <- "01"
      for (mo_orig in mo_min:12){
        yr <- yr_orig + floor((mo_orig - months_lag - 1)/12)
        yr <- as.character(yr)
        mo <- (mo_orig - months_lag) - 12*floor((mo_orig - months_lag - 1)/12)
        mo <- as.character(mo)
        if (as.numeric(mo) <= 9) mo <- paste0(0, mo)
        if (as.numeric(mo_orig) <= 9) mo_orig <- paste0(0, mo_orig)
        f_dir <- paste0(cov_name, '/', measure, '/1m/', cov_name, '_', measure, '_1m_', yr, '_', mo, '_00.tif')
        if (file.exists(paste0(cov_dir, f_dir))){
          cov_stack_new <- raster(paste0(cov_dir, f_dir))
          names(cov_stack_new) <- paste0(title, '_', measure, '_1m_', yr_orig, '_', mo_orig, '_00')
          cov_stack <- addLayer(cov_stack, cov_stack_new)
        }
      }
      # save monthly covariate data stack
      if (nlayers(cov_stack)!=0) writeRaster(cov_stack, paste0(stack_name))
    }
  }
} else {
  for (yr_orig in all_yrs){
    print(yr_orig)
    stack_name <- paste0(stack_out_dir, cov_name, '_yearly.grd')
    # .tif file path
    f_dir <- paste0(cov_name, '/', measure, '/1y/', cov_name, '_', measure, '_1y_', yr_orig, '_00_00.tif')
    # stack and save yearly covariates
    if (file.exists(paste0(cov_dir, f_dir))) cov_stack <- addLayer(cov_stack, raster(paste0(cov_dir, f_dir)))
  }
}

writeRaster(cov_stack, paste0(stack_name), overwrite=T)