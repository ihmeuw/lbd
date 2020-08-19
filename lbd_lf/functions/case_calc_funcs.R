### Functions for creating viz inputs
##########################################################################################################################################################

make_case_rasters_parallel <- function(diag = "anti", ind = indicator, ig = indicator_group, rd = run_date, pop_measure = NULL, output_dir = <<<< FILEPATH REDACTED >>>>, indic_repo) {
  # define directories
  message(rd)
  sharedir <- <<<< FILEPATH REDACTED >>>>
  input_dir <- <<<< FILEPATH REDACTED >>>>
  
  # get config args
  config <- fread(<<<< FILEPATH REDACTED >>>>)
  if (is.null(pop_measure) == T) pop_measure <- config[V1 == "pop_measure", V2]
  raster_agg_factor <- config[V1 == "raster_agg_factor", V2] %>% as.numeric()
  region <- reg <- config[V1 == "region_list", V2]
  shapefile_version <- config[V1 == "modeling_shapefile_version", V2]
  predict_years <- config[V1 == "predict_years", V2]
  if (class(predict_years) == "character") predict_years <- eval(parse(text = predict_years))
  fixed_effects_config <- fread(<<<< FILEPATH REDACTED >>>>)
  
  # load spatial templates
  load(<<<< FILEPATH REDACTED >>>>)
  
  # Pulling Population
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Pull annual population brick using new covariates function
  message(paste0("Pulling ", pop_measure, " population raster."))
  
  worldpop_current <- fixed_effects_config[covariate == "worldpop",]
  worldpop_current$year_lag <- 0
  
  pop <- load_lagged_covariates(covariate_config = worldpop_current,
                                template = simple_raster,
                                start_year = min(predict_years),
                                end_year = max(predict_years),
                                raster_agg = as.integer(raster_agg_factor))$worldpop
  
  pop <- raster::crop(pop, extent(simple_raster))
  pop <- setExtent(pop, simple_raster)
  pop <- raster::mask(pop, simple_raster)
  
  lower_prev <- brick(<<<< FILEPATH REDACTED >>>>)
  mean_prev <- brick(<<<< FILEPATH REDACTED >>>>)
  upper_prev <- brick(<<<< FILEPATH REDACTED >>>>)
  
  for (i in 2:length(predict_years)) {
    lower_prev[[i]] <- brick(<<<< FILEPATH REDACTED >>>>)
    mean_prev[[i]] <- brick(<<<< FILEPATH REDACTED >>>>)
    upper_prev[[i]] <- brick(<<<< FILEPATH REDACTED >>>>)
  }
  
  lower_cases <- lapply(1:length(predict_years), FUN = function(y) lower_prev[[y]] * pop[[y]]) %>% unlist()
  mean_cases <- lapply(1:length(predict_years), FUN = function(y) mean_prev[[y]] * pop[[y]]) %>% unlist()
  upper_cases <- lapply(1:length(predict_years), FUN = function(y) upper_prev[[y]] * pop[[y]]) %>% unlist()
  
  ## make a mask and NA lakes
  if (reg == "lf_endem_afr") {
    msk <- raster(<<<< FILEPATH REDACTED >>>>)
    lks <- raster(<<<< FILEPATH REDACTED >>>>)
    
    ## set prepped simple_raster object to be 1 where it's not NA
    ## crop mean raster
    a <- crop(simple_raster, mean_prev)
    a <- extend(a, mean_prev)
    a2 <- copy(a)
    a2[!is.na(a2)] <- 1
    if (exists("ctry_to_NA")) {
      a[!a %in% g$code[!g$iso3 %in% ctry_to_NA]] <- NA
    }
    
    msk[is.na(msk)] <- 0
    msk[msk == 1] <- NA
    msk[msk == 0] <- 1
    
    ## crop and extend the mask
    msk <- crop(msk, mean_prev)
    msk <- extend(msk, mean_prev)
    
    msk <- msk * a2
    a[!is.na(a)] <- 1
    
    ## crop and extend the mask
    msk <- crop(msk, mean_prev)
    msk <- extend(msk, mean_prev)
    
    ## save the mask to the output dir
    writeRaster(
      msk,
      file = paste0(<<<< FILEPATH REDACTED >>>>),
      format = "GTiff",
      datatype = "INT1U",
      NAflag = 0,
      overwrite = TRUE
    )
  } else {
    lks <- raster(<<<< FILEPATH REDACTED >>>>)
  }
  
  lks <- extend(lks, mean_prev)
  lks <- crop(lks, mean_prev)
  extent(lks) <- extent(mean_prev)
  lks[lks > 1] <- NA
  lks[is.na(lks)] <- 0
  lks[lks == 1] <- NA ## set lake locations to be NA
  lks[lks == 0] <- 1 ## otherwise set to 1
  
  output_ras_list <- c("lower_prev", "mean_prev", "upper_prev", "lower_cases", "mean_cases", "upper_cases")
  output_ras <- list(lower_prev, mean_prev, upper_prev, lower_cases, mean_cases, upper_cases)
  names(output_ras) <- output_ras_list
  
  for (ras_name in output_ras_list) { # for each prev/cases raster, save one rasterlayer per year
    for (i in 1:length(predict_years)) {
      out <- output_ras[[ras_name]][[i]]
      out <- out * lks # apply lakes mask
      if (length(grep("prev", ras_name)) > 0) out <- out * 100 # multiply 100 if is prevalence
      if (length(grep("cases", ras_name)) > 0) out <- round(out, 0) # round to integer if cases
      out[is.na(out)] <- -999999
      writeRaster(
        out,
        file = paste0(<<<< FILEPATH REDACTED >>>>),
        format = "GTiff",
        datatype = "FLT4S",
        NAflag = -999999,
        overwrite = TRUE
      )
    }
  }
}

make_case_rasters <- function(ind = indicator, ig = indicator_group, rd = run_date, pop_measure = NULL, output_dir = <<<< FILEPATH REDACTED >>>>, indic_repo) {
  # define directories
  message(rd)
  sharedir <- sprintf(<<<< FILEPATH REDACTED >>>>)
  input_dir <- paste0(<<<< FILEPATH REDACTED >>>>)
  
  # get config args
  config <- fread(<<<< FILEPATH REDACTED >>>>)
  if (is.null(pop_measure) == T) pop_measure <- config[V1 == "pop_measure", V2]
  raster_agg_factor <- config[V1 == "raster_agg_factor", V2] %>% as.numeric()
  region <- reg <- config[V1 == "region_list", V2]
  shapefile_version <- config[V1 == "modeling_shapefile_version", V2]
  predict_years <- config[V1 == "predict_years", V2]
  if (class(predict_years) == "character") predict_years <- eval(parse(text = predict_years))
  fixed_effects_config <- fread(paste0(indic_repo, config[V1 == "covs_name", V2], ".csv"))
  
  # load spatial templates
  load(<<<< FILEPATH REDACTED >>>>)

  # Pulling Population
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## Pull annual population brick using new covariates function
  message(paste0("Pulling ", pop_measure, " population raster."))

  worldpop_current <- fixed_effects_config[covariate == "worldpop_raked",]
  worldpop_current$year_lag <- 0
  
  pop <- load_lagged_covariates(covariate_config = worldpop_current,
                                template = simple_raster,
                                start_year = min(predict_years),
                                end_year = max(predict_years),
                                raster_agg = as.integer(raster_agg_factor))$worldpop
  
  pop <- raster::crop(pop, extent(simple_raster))
  pop <- setExtent(pop, simple_raster)
  pop <- raster::mask(pop, simple_raster)

  lower_prev <- brick(<<<< FILEPATH REDACTED >>>>)
  mean_prev <- brick(<<<< FILEPATH REDACTED >>>>)
  upper_prev <- brick(<<<< FILEPATH REDACTED >>>>)
  
  lower_cases <- lapply(1:length(predict_years), FUN = function(y) lower_prev[[y]] * pop[[y]]) %>% unlist()
  mean_cases <- lapply(1:length(predict_years), FUN = function(y) mean_prev[[y]] * pop[[y]]) %>% unlist()
  upper_cases <- lapply(1:length(predict_years), FUN = function(y) upper_prev[[y]] * pop[[y]]) %>% unlist()

  ## make a mask and NA lakes
  if (reg == "lf_endem_afr") {
    msk <- raster(<<<< FILEPATH REDACTED >>>>)
    lks <- raster(<<<< FILEPATH REDACTED >>>>)

    ## set prepped simple_raster object to be 1 where it's not NA
    ## crop mean raster
    a <- crop(simple_raster, mean_prev)
    a <- extend(a, mean_prev)
    a2 <- copy(a)
    a2[!is.na(a2)] <- 1
    if (exists("ctry_to_NA")) {
      a[!a %in% g$code[!g$iso3 %in% ctry_to_NA]] <- NA
    }

    msk[is.na(msk)] <- 0
    msk[msk == 1] <- NA
    msk[msk == 0] <- 1

    ## crop and extend the mask
    msk <- crop(msk, mean_prev)
    msk <- extend(msk, mean_prev)

    msk <- msk * a2
    a[!is.na(a)] <- 1

    ## crop and extend the mask
    msk <- crop(msk, mean_prev)
    msk <- extend(msk, mean_prev)

    ## save the mask to the output dir
    writeRaster(
      msk,
      file = paste0(<<<< FILEPATH REDACTED >>>>),
      format = "GTiff",
      datatype = "INT1U",
      NAflag = 0,
      overwrite = TRUE
    )
  } else {
    lks <- raster(<<<< FILEPATH REDACTED >>>>)
  }

  lks <- extend(lks, mean_prev)
  lks <- crop(lks, mean_prev)
  extent(lks) <- extent(mean_prev)
  lks[lks > 1] <- NA
  lks[is.na(lks)] <- 0
  lks[lks == 1] <- NA ## set lake locations to be NA
  lks[lks == 0] <- 1 ## otherwise set to 1

  output_ras_list <- c("lower_prev", "mean_prev", "upper_prev", "lower_cases", "mean_cases", "upper_cases")
  output_ras <- list(lower_prev, mean_prev, upper_prev, lower_cases, mean_cases, upper_cases)
  names(output_ras) <- output_ras_list

  for (ras_name in output_ras_list) { # for each prev/cases raster, save one rasterlayer per year
    for (i in 1:length(predict_years)) {
      out <- output_ras[[ras_name]][[i]]
      out <- out * lks # apply lakes mask
      if (length(grep("prev", ras_name)) > 0) out <- out * 100 # multiply 100 if is prevalence
      if (length(grep("cases", ras_name)) > 0) out <- round(out, 0) # round to integer if cases
      out[is.na(out)] <- -999999
      writeRaster(
        out,
        file = paste0(<<<< FILEPATH REDACTED >>>>),
        format = "GTiff",
        datatype = "FLT4S",
        NAflag = -999999,
        overwrite = TRUE
      )
    }
  }
}

make_case_adm_tables <- function(ind, ig, rd, pop_measure = NULL) {

  sharedir <- sprintf(<<<< FILEPATH REDACTED >>>>)
  input_dir <- paste0(<<<< FILEPATH REDACTED >>>>)

  config <- fread(<<<< FILEPATH REDACTED >>>>)
  if (is.null(pop_measure) == T) pop_measure <- config[V1 == "pop_measure", V2]
  raster_agg_factor <- config[V1 == "raster_agg_factor", V2] %>% as.numeric()
  region <- reg <- config[V1 == "region_list", V2]
  shapefile_version <- config[V1 == "modeling_shapefile_version", V2]
  predict_years <- config[V1 == "predict_years", V2]
  if (class(predict_years) == "character") predict_years <- eval(parse(text = predict_years))

  load(<<<< FILEPATH REDACTED >>>>)

  adm_tables <- list(admin_0, admin_1, admin_2)
  adm_tables <- lapply(adm_tables, FUN = function(tb) {
    mean_prev <- tb[, grep("V", names(tb)), with = F] %>% rowMeans()
    lower_prev <- as.matrix(tb[, grep("V", names(tb)), with = F]) %>% rowQuantiles(probs = .025)
    upper_prev <- as.matrix(tb[, grep("V", names(tb)), with = F]) %>% rowQuantiles(probs = .975)
    keep_names <- grep("ADM", names(tb), value = T)
    keep_names <- c(keep_names, "year", "pop")
    out_tb <- cbind(tb[, keep_names, with = F], lower_prev, mean_prev, upper_prev)
    out_tb[, c("lower_cases", "mean_cases", "upper_cases") := list(lower_prev * pop, mean_prev * pop, upper_prev * pop)]
    setnames(out_tb, grep("ADM", names(tb), value = T), "gaul_code")
    return(out_tb)
  })

  summary_list <- c("lower_prev", "mean_prev", "upper_prev", "lower_cases", "mean_cases", "upper_cases")
  summary_out <- lapply(summary_list, FUN = function(summ) {
    adm_return <- lapply(adm_tables, FUN = function(tb) {
      out_tb <- tb[, c("gaul_code", "year", summ), with = F]
      setnames(out_tb, summ, "value")
      out_tb <- out_tb[!is.na(value) & !is.nan(value)]
      if (length(grep("cases", summ)) > 0) out_tb$value <- round(out_tb$value, 0) # round to integer if cases
      if (length(grep("prev", summ)) > 0) out_tb$value <- out_tb$value * 100 # multiply 100 if is prevalence

      return(out_tb)
    }) %>% rbindlist()
  })
  names(summary_out) <- summary_list
  return(summary_out)
}
