### Wrapper for set_up_config(), but with covs file specified in config file, and with loading of covariate lags
set_up_config_focal_3 <- function(repo, 
                                  core_repo = repo,
                                  indicator_group, 
                                  indicator, 
                                  config_file = NULL, 
                                  run_tests = TRUE,
                                  return_list = FALSE) {

  ## Load config
  config <- data.table::fread(config_file, header = FALSE)
  
  ## For parsimony, let's make sure that the config column names are V1 and V2
  config <- data.table(config)
  if(colnames(config)[1] != "V1" & colnames(config)[2] != "V2") {
    warning("Renaming config column names to V1 and V2. Please verify that 'config' is properly built")
    colnames(config) <- c("V1", "V2")
  }
  
  ## Check for covs file specified in config
  if (!is.null(config[V1 == "covs_name", V2]) & !(config[V1 == "covs_name", V2] == "")) {
    covs_file <- paste0(config[V1 == "covs_name", V2], ".csv")
  }
  
  ## Now call set_up_config()
  config <- set_up_config(repo = indic_repo, core_repo = core_repo, indicator_group = indicator_group, indicator = indicator,
                          config_file = config_file, covs_file = covs_file, run_tests = run_tests, return_list = return_list)
  
  fixed_effects_lags <- paste0(fixed_effects_config$year_lag, collapse = " + ")
  assign("fixed_effects_lag", fixed_effects_lags, envir = globalenv())
  
  assign("fixed_effects", fixed_effects[[2]], envir = globalenv())
  
  return(config)
}
