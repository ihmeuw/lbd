  ##############################################################################
  ## Variance inflation factor (VIF) launcher script
  #REDACTED
  
  ## Setup -------------------------------------------------------------------------
  
  # clear environment
  rm(list = ls())
  
  # set general arguments
  #REDACTED
  
  # Load MBG packages and functions
  message('Loading in required R packages and MBG functions')
  package_list <- c(t(read.csv('#REDACTED/package_list.csv',header=FALSE)))
  source(paste0(repo, '/mbg_central/setup.R'))
  mbg_setup(package_list = package_list, repos = repo)
  
  # set node preference
  use_geos_nodes <- TRUE
  if (use_geos_nodes) {
    proj <- 'proj_geo_nodes'
    r_shell <- 'shell_geos.sh'
  } else {
    proj <- 'proj_geospatial'
    r_shell <- 'shell_prod.sh'
  }
  
  # set config and covariate files
  config_par   <- 'hap_sp_fine'
  config_file <- file.path(indicator_group, 'model/configs/')
  cov_par      <- 'hap_standard'
  cov_file <- config_file
  
  # indicate whether to use old run date
  use_old_run_date <- FALSE
  old_run_date_input <- ''
  
  # set run date
  if (use_old_run_date == FALSE) {
    run_date <- make_time_stamp(TRUE)
  } else {
    run_date <- old_run_date_input
  }
  
  # set whether running individual countries
  individual_countries <- FALSE
  
  # indicate threshold parameters
  threshold_min <- 2
  threshold_max <- 5
  threshold_step <- 1
  
  # indicate whether or not to crop covariates (must give .RData file with cropped covariate file if crop_covs = FALSE)
  crop_covs <- TRUE
  cropped_covs_file <- ''
  
  # list indicators
  indics <- c('cooking_fuel_solid')
  
  # list regions
  if (indicator_group=='cooking') {
    
    regions <- c('essa-ERI-DJI-YEM', "ERI+DJI+YEM", 'essa',
                 'sssa-ZAF', 'ZAF', 'sssa',
                 'cssa-AGO-GNQ', 'AGO', 'cssa',
                 'wssa-CPV-NGA', 'NGA', 'wssa',
                 'noaf-ESH',
                 'caca-CUB',
                 'ansa-VEN', 'trsa-GUF',
                 'stan-TKM', 'stan',
                 'CHN', 'MNG', 
                 'ocea-MYS', 
                 'seas-VNM-THA', 'VNM', 'THA', 'seas',
                 'mide+TKM', 'mide', 'soas')
    
  } else {
    regions <- c('dia_afr_horn', 'dia_name', 'dia_sssa', 
                 'dia_mcaca', 'dia_s_america', 'dia_central_asia', 'dia_chn_mng', 
                 'dia_se_asia', 'dia_malay', 'dia_mid_east',
                 'ZWE', 'KEN', 'NGA', 'COD', 'IND', 'PAK',
                 'dia_essa-zwe-ken', 'dia_wssa-nga', 'dia_cssa-cod', 'dia_south_asia-ind-pak')
  }
  
  
  
  
  ## Run vif selection script -------------------------------------------------------------------------
  
  for (i in indics) {
    
    for (reg in regions) {
      
      # set specific arguments
      indicator       <- i
      jname           <- paste0(indicator, '_', reg, '_vif_selection')
      mymem <- '100G'
      sys.sub <- paste0('qsub -e #REDACTED',
                        '-l m_mem_free=', mymem, ' -P ', proj, ifelse(use_geos_nodes, ' -q geospatial.q ', ' '),
                        '-l fthread=1 -l h_rt=00:12:00:00 -v sing_image=default -N ', jname, ' ')
      r_shell <- file.path(repo, '#REDACTED')
      script <- paste0(repo, 'modules/covariate_selection/02_vif_selection_script.R')
      args <- paste(user, repo, indicator_group, indicator, config_par, config_file, cov_par, cov_file, 
                    run_date, reg, threshold_min, threshold_max, threshold_step, individual_countries,
                    crop_covs, cropped_covs_file)
      
      # run launch script
      paste(sys.sub, r_shell, script, args) %>% system
    
    }
    
  }
  #***********************************************************************************************************************