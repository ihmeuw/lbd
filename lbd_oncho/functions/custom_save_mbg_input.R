save_mbg_input <- function(indicator = indicator, indicator_group = indicator_group, df = df , simple_raster = simple_raster, simple_raster2 = NULL, simple_raster_raf1 = simple_raster, mesh_s = mesh_s, mesh_t = mesh_t, cov_list = all_cov_layers, run_date = NULL, pathaddin="", child_model_names = NULL, all_fixed_effects = NULL, period_map = NULL, centre_scale = T, rw1_raw_covar_list = NULL) {
  
  period_map <- copy(period_map) ## to avoid global scoping issues

  if (!is.null(rw1_raw_covar_list)) {
    rw1_raw_covs <- extract_covariates(df, cov_list[which(names(cov_list) %in% rw1_raw_covar_list)], return_only_results = T, centre_scale = FALSE, period_var = 'year', period_map = period_map)
    just_covs <- extract_covariates(df, cov_list, return_only_results = T, centre_scale = centre_scale, period_var = 'year', period_map = period_map)
    if(centre_scale==TRUE){
      just_covs <- just_covs$covs
    }
    keep_cols <- which(!(colnames(just_covs) %in% rw1_raw_covar_list))
    just_covs <- just_covs[, ..keep_cols]
    just_covs <- just_covs[, year := NULL]
    just_covs <- just_covs[, period_id := NULL]
    
    if (length(just_covs) > 0) {
      just_covs <- cbind(just_covs, rw1_raw_covs)
    } else {
      just_covs <- rw1_raw_covs
    }

    keep_cols <- which(!(colnames(df) %in% rw1_raw_covar_list))
    df <- df[, ..keep_cols]
    keep_cols <- which(!(colnames(just_covs) %in% colnames(df)))
    
    df <- cbind(df, just_covs[, ..keep_cols])
  } else {
    just_covs <- extract_covariates(df, cov_list, return_only_results = T, centre_scale = centre_scale, period_var = 'year', period_map = period_map)
    if(centre_scale==TRUE){
      just_covs <- just_covs$covs
    }
    just_covs <- just_covs[, year := NULL]
    just_covs <- just_covs[, period_id := NULL]

    keep_cols <- which(!(colnames(just_covs) %in% colnames(df)))

    df <- cbind(df, just_covs[, ..keep_cols])
  }
  
  #create a period column
  if(is.null(period_map)) {
    period_map <- make_period_map(c(2000,2005,2010,2015))
  }
  df[, period := NULL]
  setnames(period_map, 'data_period', 'year')
  setnames(period_map, 'period_id', 'period')
  df <- merge(df, period_map, by = 'year', sort =F)
  
  ## Now that we've extracted covariate values to our data in the buffer zone, clip cov list to simple_raster instead of simple_polygon
  ##    (simple_raster is area we actually want to model over)
  for(l in 1:length(cov_list)) {
    cov_list[[l]]  <- crop(cov_list[[l]], extent(simple_raster_raf1))
    cov_list[[l]]  <- setExtent(cov_list[[l]], simple_raster_raf1)
    cov_list[[l]]  <- raster::mask(cov_list[[l]], simple_raster_raf1)
  }
  
  ## Save all inputs
  to_save <- c("df", "simple_raster", "mesh_s", "mesh_t", 'cov_list', 'child_model_names', 'all_fixed_effects','period_map')
  if(!is.null(simple_raster2)) to_save <- c(to_save, "simple_raster2", "simple_raster_raf1")
  save(list = to_save, file = <<<< FILEPATH REDACTED >>>>)
  message(<<<< FILEPATH REDACTED >>>>)
  
}
