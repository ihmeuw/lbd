  ##############################################################################
  ## MBG launch, aggregate results, and diagnostics launcher script 
  ## Created 2018/02/23
  ##############################################################################
  ## Setup -------------------------------------------------------------------------
  
  # clear environment
  rm(list = ls())
  
  #REDACTED
  
  # Load MBG packages and functions
  message('Loading in required R packages and MBG functions')
  package_list <- c(t(read.csv('/#REDACTED/package_list.csv',header=FALSE)))
  source(paste0(core_repo, '/mbg_central/setup.R'))
  mbg_setup(package_list = package_list, repos = core_repo)
  
  # set cluster arguments
  use_geos_nodes  <- F
  proj_arg        <- ifelse(use_geos_nodes, 'proj_geo_nodes', 'proj_geospatial')
  proj            <- ifelse(use_geos_nodes, paste0(' -P ', proj_arg, ' -l gn=TRUE '), paste0(' -P ', proj_arg, ' '))
  
  # set script arguments
  skip_entry <- F #use to skip to descent stage
  
  # set covariate arguments
  plot_covariates <- TRUE
  covariate_plotting_only <- FALSE
  
  # indicate whether to use old run date
  use_old_run_date <- T
  old_run_date_input <- '2020_09_01_11_42_52'
  
  # set run date
  if (use_old_run_date == FALSE) {
    run_date <- make_time_stamp(TRUE)
  } else {
    run_date <- old_run_date_input
  }

  # custom region list
  regions <- c('essa-ERI-DJI-YEM', "ERI+DJI+YEM",
               'sssa-ZAF', 'ZAF',
               'cssa-AGO-GNQ', 'AGO',
               'wssa-CPV-NGA', 'NGA',
               'noaf-ESH',
               'caca-CUB',
               'ansa-VEN', 'trsa-GUF',
               'stan-TKM',
               'CHN', 'MNG',
               'ocea-MYS',
               'seas',
               'mide+TKM', 'soas')

  #REDACTED

if(!skip_entry) {  
  for (reg in regions) {
  
    # set memory based on region
    if (reg %in% c('trsa-GUF')) { mymem <- '750G'
    } else if (reg %in% c('wssa-CPV-NGA', 'CHN', 'soas', 'ansa-VEN', 'ocea-MYS')) { mymem <- '500G'
    } else if (reg %in% c('seas-VNM-THA', 'VNM', 'THA', "ERI+DJI+YEM", 'sssa-ZAF')) { mymem <- '200G'
    } else mymem <- '350G'
  
        #name job
        jname           <- paste('eDL', reg, indicator, sep = '_')
  
        #setup covars file
        cov_par <- paste(indicator_group, reg, sep='_')
  
        # set up qsub
        sys.sub <- paste0('qsub -e #REDACTED/', user,'/errors -o #REDACTED/', user, '/output ',
                          '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                          '-l fthread=1 -l h_rt=', ifelse(use_geos_nodes, '16:00:00:00', '3:00:00:00'),
                          ' -v sing_image=default -N ', jname, ' -l archive=TRUE ')
        r_shell <- file.path(core_repo, '#REDACTED/shell_sing.sh')
        script <- file.path(my_repo, indicator_group, 'post/2_entry.R')
        args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, run_date, measure, holdout, my_repo)
  
        # run launch script
        paste(sys.sub, r_shell, script, args) %>%
          system
  
  }

} else {
  #use to launch descent instead if necessary
  
  for (reg in regions) {
  
    # set memory based on region
    if (reg %in% c('CHN', 'trsa-GUF', 'ansa-VEN')) { mymem <- '750G'
    } else if (reg %in% c('wssa-CPV-NGA', 'trsa-GUF', 'CHN', 'soas', 'ocea-MYS')) { mymem <- '500G'
    } else if (reg %in% c('seas-VNM-THA', 'VNM', 'THA', "ERI+DJI+YEM", 'sssa-ZAF')) { mymem <- '100G'
    } else mymem <- '300G'
  
      #name job
      jname           <- paste('EdL', reg, indicator, sep = '_')
  
      #setup covars file
      cov_par <- paste(indicator_group, reg, sep='_')
      
      # set up qsub
      sys.sub <- paste0('qsub -e /#REDACTED', user,'/errors -o /#REDACTED', user, '/output ',
                        '-l m_mem_free=', mymem, ' -P ', proj_arg, ifelse(use_geos_nodes, ' -q geospatial.q ', ' -q all.q '),
                        '-l fthread=1 -l h_rt=', ifelse(use_geos_nodes, '16:00:00:00', '3:00:00:00'),
                        ' -v sing_image=default -N ', jname, ' -l archive=TRUE ')
      r_shell <- file.path(core_repo, '#REDACTED/shell_sing.sh')
      script <- file.path(my_repo, indicator_group, 'post/3_descent.R')
      args <- paste(user, core_repo, indicator_group, indicator, config_par, cov_par, reg, run_date, measure, holdout, my_repo)
  
      # run launch script
      paste(sys.sub, r_shell, script, args) %>%
        system
  
  }
}