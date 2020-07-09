# --------------------------------------------------------------------------------------------------
# Script designed to facilitate covariate selection using variance inflation factor (VIF) analysis
# The purpose of the VIF analysis is to test for multicollinearity in covariates
#
# (1) load covariate data cropped to MBG input data by region
# (2) perform VIF analysis and drop covariates above a threshold by region
# (3) take results and summarize in a table that is human-readable
# (4) take results and summarize in MBG input format
#
# References: 
# Faraway's Linear Models in R, Chapter 4: Diagnostics
# https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/
# --------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
# Setup

# load arguments from qsub
user            <- commandArgs()[4]
repo            <- commandArgs()[5]
indicator_group <- commandArgs()[6]
indicator       <- commandArgs()[7]
config_par      <- commandArgs()[8]
config_file     <- commandArgs()[9]
cov_par         <- commandArgs()[10]
cov_file        <- commandArgs()[11]
run_date        <- commandArgs()[12]
regions         <- commandArgs()[13]
threshold_min   <- commandArgs()[14]
threshold_max   <- commandArgs()[15]
threshold_step  <- commandArgs()[16]
individual_countries <- commandArgs()[17]
crop_covs <- commandArgs()[18]
cropped_covs_file <- commandArgs()[19]
covs_to_keep    <- commandArgs()[20]
message(indicator)

# set additional arguments
core_repo <- repo
cov_repo <- paste0(repo, '/3_modeling/vif_selected_covs/')
age <- 0
holdout <- 0
threshold <-  c(seq(as.numeric(threshold_min), as.numeric(threshold_max), by = as.numeric(threshold_step)))

# list regions
if (regions == 'all') {
  regions <- c('dia_afr_horn-eth-yem', 'ETH', 'YEM', 'dia_name', 'dia_sssa', 
               'dia_mcaca', 'dia_central_asia', 'MNG', 
               'dia_se_asia', 'dia_malay', 'dia_mid_east',
               'ZWE', 'KEN', 'NGA', 'COD', 'IND', 'PAK',
               'dia_essa-zwe-ken', 'dia_cssa-cod', 'dia_south_asia-ind-pak',
               'dia_s_america_n', 'dia_s_america_s')
}

# create directories to save results in
save_dir <- paste0('<<<< FILEPATH REDACTED >>>>/covariate_selection_', run_date, '/')
dir.create(save_dir)
mbg_input_dir <- paste0(save_dir, 'mbg_inputs/')
dir.create(mbg_input_dir)

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## Read config file and save all parameters in memory
config <- set_up_config(repo            = core_repo,
                        indicator_group = '',
                        indicator       = '',
                        config_name     = paste0(config_file, config_par),
                        covs_name       = paste0(cov_file, cov_par))

## Create proper year list object
if (class(year_list) == 'character') year_list <- eval(parse(text=year_list))

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
check_config()

# If running individual countries make sure all country FEs and REs off
if (individual_countries) {
  use_child_country_fes <- FALSE
  use_inla_country_fes  <- FALSE
  use_inla_country_res  <- FALSE
}
# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Crop covariates to data

if (crop_covs == TRUE) {

  # set up blank list to store cropped covariate data
  all_cov_data <- list()
  
  # start loop over regions
  for (region in regions) {
    
    
    # --------------------------------------------------------------------------------------------------
    # Prep MBG inputs and load data
    
    # report region and indicator
    message(paste0('\n', indicator))
    reg <- region
    message(reg)
    
    # make pathaddin
    pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
    
    ## some set up
    if (class(year_list) == 'character') year_list <- eval(parse(text=year_list))
    if (class(z_list)    == 'character') z_list    <- eval(parse(text=z_list))
    
    ## Load simple polygon template to model over
    gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version) 
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4)
    subset_shape        <- simple_polygon_list[[1]]
    simple_polygon      <- simple_polygon_list[[2]]
    
    ## Load list of raster inputs (pop and simple)
    raster_list        <- build_simple_raster_pop(subset_shape)
    simple_raster      <- raster_list[['simple_raster']]
    pop_raster         <- raster_list[['pop_raster']]
    
    ## Load input data
    df <- load_input_data(indicator   = indicator,
                          simple      = simple_polygon,
                          agebin      = age,
                          pathaddin   = pathaddin,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list)
    
    ## Remove any data outside of region for ORS
    reg_list <- fread('<<<< FILEPATH REDACTED >>>>/dia_region_iso.csv')
    reg_list[, gadm_adm0_code := get_adm0_codes(iso3, shapefile_version = modeling_shapefile_version), by = iso3]
    iso3_list <- filter(reg_list, gadm_adm0_code %in% gaul_list)
    iso3_list <- unique(iso3_list$iso3)
    df <- filter(df, country %in% iso3_list)
    df <- as.data.table(df)
    
    # remove unecessary folder and data file
    unlink(paste0('<<<< FILEPATH REDACTED >>>>'), recursive = TRUE)
    
    ## Some built in data checks that cause known problems later on
    if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
      stop('You have binomial data where k > N. Check your data before proceeding')
    if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
      stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')
    # --------------------------------------------------------------------------------------------------
    
    
    # --------------------------------------------------------------------------------------------------------------
    # Pull covariates
    
    ## Define modeling space. In years only for now.
    if(yearload=='annual') period_map <-
      make_period_map(modeling_periods = c(min(year_list):max(year_list)))
    if(yearload=='five-year') period_map <-
      make_period_map(modeling_periods = seq(min(year_list),max(year_list),by=5))
    
    ## Make placeholders for covariates
    cov_layers <- gbd_cov_layers <- NULL
    
    ## Pull all covariate bricks/layers
    if (nchar(fixed_effects) > 0) {
      message('Grabbing raster covariate layers')
      
      effects <- trim(strsplit(fixed_effects, '\\+')[[1]])
      measures <- trim(strsplit(fixed_effects_measures, '\\+')[[1]])
      cov_layers <- load_and_crop_covariates_annual(covs            = effects,
                                                    measures        = measures,
                                                    simple_polygon  = simple_polygon,
                                                    start_year      = min(year_list),
                                                    end_year        = max(year_list),
                                                    interval_mo     = as.numeric(interval_mo))
    }
    
    ## Pull country level gbd covariates
    if (nchar(gbd_fixed_effects) > 0) {
      message('Grabbing GBD covariates')
      
      # can't pull gbd covs past 2017 currently
      if (max(year_list) > 2017) {
        gbd_yr_list <- year_list[-c(which(year_list > 2017))]
      } else {
        gbd_yr_list <- year_list
      }
      
      effects <- trim(strsplit(gbd_fixed_effects, '\\+')[[1]])
      measures <- trim(strsplit(gbd_fixed_effects_measures, '\\+')[[1]])
      gbd_cov_layers <- load_gbd_covariates(covs     = effects,
                                            measures = measures,
                                            year_ids = gbd_yr_list,
                                            age_ids  = gbd_fixed_effects_age,
                                            template = cov_layers[[1]][[1]],
                                            simple_polygon = simple_polygon,
                                            interval_mo = interval_mo)
      
      # copy 2017 covs to later years
      if (max(year_list) > 2017) {
        for (c in names(gbd_cov_layers)) {
          for (i in which(year_list > 2017)) {
            gbd_cov_layers[[c]][[i]] <- gbd_cov_layers[[c]][[17]]
          }
        }
      }
      
    }
  
    ## Combine all covariates
    all_cov_layers <- c(cov_layers, gbd_cov_layers)
    
    ## regenerate all fixed effects equation from the cov layers
    all_fixed_effects <- paste(names(all_cov_layers), collapse = ' + ')
    
    ## Make stacker-specific formulas where applicable
    all_fixed_effects_brt <- all_fixed_effects
    
    ## Set Up Country Fixed Effects
    if(use_child_country_fes == TRUE | use_inla_country_fes == TRUE) {
      message('Setting up country fixed effects')
      fe_gaul_list <- unique(c(gaul_convert(unique(df[, country])), gaul_list))
      fe_template  <- cov_layers[[1]][[1]]
      simple_polygon_list <- load_simple_polygon(gaul_list   = fe_gaul_list,
                                                 buffer      = 0.4,
                                                 subset_only = TRUE)
      fe_subset_shape     <- simple_polygon_list[[1]]
      gaul_code <- rasterize(fe_subset_shape, fe_template, field = 'GAUL_CODE')
      gaul_code <- setNames(gaul_code,'gaul_code')
      gaul_code <- create_categorical_raster(gaul_code)
      
      ## update covlayers and add country fixed effects to the
      all_cov_layers = update_cov_layers(all_cov_layers, gaul_code)
      all_fixed_effects_cfes = paste(all_fixed_effects,
                                     paste(names(gaul_code)[1:length(names(gaul_code))],
                                           collapse = ' + '), sep=' + ')
      
      ## update specific stacker formulas (for now we just want country effects in BRT)
      all_fixed_effects_brt <- all_fixed_effects_cfes
    }
    
    ## Add these to the fixed effects if we want them in stacking
    if(use_child_country_fes == TRUE) {
      gaul_fes <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = ' + ')
      all_fixed_effects = paste(all_fixed_effects, gaul_fes, sep = ' + ')
    }
    
    # get cov names
    the_covs <- format_covariates(all_fixed_effects)
    
    ## copy the dataset to avoid unintended namespace conflicts
    the_data <- copy(df)
    
    ## add a row id column
    the_data[, a_rowid := seq(1:nrow(the_data))]
    
    ## extract covariates to the points and subset data where its missing covariate values
    cs_covs <- extract_covariates(the_data,
                                  all_cov_layers,
                                  id_col              = 'a_rowid',
                                  return_only_results = TRUE,
                                  centre_scale        = TRUE,
                                  period_var          = 'year',
                                  period_map          = period_map)
    
    
    ## Check for data where covariate extraction failed
    rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
    if (rows_missing_covs > 0) {
      pct_missing_covs <- round((rows_missing_covs/nrow(the_data))*100, 2)
      warning(paste0(rows_missing_covs, ' out of ', nrow(the_data), ' rows of data ',
                     '(', pct_missing_covs, '%) do not have corresponding ',
                     'covariate values and will be dropped from child models...'))
      if (rows_missing_covs/nrow(the_data) > 0.1) {
        stop(paste0('Something has gone quite wrong: more than 10% of your data does not have ',
                    'corresponding covariates.  You should investigate this before proceeding.'))
      }
    }
    
    the_data <- merge(the_data, cs_covs[[1]], by = 'a_rowid', all.x = F, all.y = F)
    
    ## store the centre scaling mapping
    covs_cs_df  <-  cs_covs[[2]]
    
    ## this will drop rows with NA covariate values
    the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))
    
    ## stop if this na omit demolished the whole dataset
    if(nrow(the_data) == 0) stop('You have an empty df, make sure one of your covariates was not NA everywhere.')
    
    ## make sure we're working with a data table
    the_data <- data.table(the_data)
    
    ## create a data table with just covariates cropped by locations in region with data
    cov_cols <- names(all_cov_layers)
    cov_data <- the_data[, cov_cols, with = FALSE]
    
    ## save data to list
    all_cov_data[[region]] <- cov_data
    # --------------------------------------------------------------------------------------------------------------
    
  } # end loop over regions
  
  
  # save data (for easier loading in the future)
  cropped_covs_file <- paste0(save_dir, 'covs_cropped_to_data_', Sys.Date(), '.RData')
  save(all_cov_data, file = cropped_covs_file)
  # ------------------------------------------------------------------------------------

  
} else {
  
  # load previously cropped covariate data
  load(cropped_covs_file)
}


# ---------------------------------------------------------------
# Stepwise VIF selection function

# Start function
stepwise_vif_selection <- function(thresh, 
                                   covariate_data, 
                                   reg,
                                   not_to_drop = NULL,
                                   trace = TRUE
) {
  # ---------------------------------------------------------------
  
  
  # ------------------------------------------------------------------------------------------------------
  # Data checks
  
  # make sure input data is a data table
  if(any(!'data.table' %in% class(covariate_data))) covariate_data<-data.table(covariate_data)
  
  # remove any covariate values that do not vary accross the region
  var_names <- names(covariate_data)
  for (val in var_names) {
    if (!is.na(val) & length(unique(covariate_data[[val]])) < 2) {
      message(paste0('\nRemoving covariate "', val, '" because it is constant across the region.'))
      covariate_data[, (val) := NULL]
    }
  }
  # ------------------------------------------------------------------------------------------------------
  
  
  # ------------------------------------------------------------------------------------------------------------------------
  # VIF selection
  
  # load vif analysis package
  library(fmsb, lib.loc = '<<<< FILEPATH REDACTED >>>>')
  
  # get initial vif value for all comparisons of variables
  vif_init<-NULL
  
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = covariate_data))))
  }
  
  vif_max<-max(as.numeric(vif_init[,2]), na.rm = TRUE)
  
  if(vif_max < thresh){
    
    if(trace==T){ # print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      
      message(paste('\nAll variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    
    cat(var_names)
    
  } else {
    
    # backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(covariate_data)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = covariate_data))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      
      # if there are any covs that we don't want dropped, set their vif to 0
      if (!is.null(not_to_drop)) {
        for (k in not_to_drop) {
          if (trace==T) message(paste0('Making sure we keep ', k))
          vif_vals[which(vif_vals[,1] == k),2] <- 0
        }
      }
      
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2]), na.rm = TRUE))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ # print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        message(paste0('removed: ',vif_vals[max_row,1], ' ', vif_max,'\n\n'))
        flush.console()
      }
      
      keep_cols <- var_names[-which(var_names == vif_vals[max_row,1])]
      covariate_data<-covariate_data[, keep_cols, with = FALSE]
      
    }
    
    message(paste0('\nCovariates selected for ', indicator, ' in ', reg, ' with a threshold of ', thresh, ':'))
    cat(names(covariate_data))
    
  }
  
  # vector of selected covariates
  selected_covs <- c(names(covariate_data))
  # ------------------------------------------------------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------------------
  # Add selection results to table
  
  # add 0/1 for whether or not covariate was selected to a data table
  d1 <- data.table(cov_names)
  d2 <- data.table(selected_covs)
  d2[, included := TRUE]
  dt <- merge(d1, d2, by.x = 'cov_names', by.y = 'selected_covs', all = T)
  dt[is.na(included), included := FALSE]
  dt <- dcast(melt(dt, id.vars = 'cov_names'), variable ~ cov_names)
  dt[, variable := NULL]
  
  # add threshold, region, number of covariates selected, and run date that corresponds to data
  dt[, region := reg]
  dt[, vif_threshold := thresh]
  dt[, num_covs_selected := length(selected_covs)]
  dt[, data_run_date := run_date]
  
  # bind to original data table
  selection_results <- rbindlist(list(selection_results, dt), use.names = TRUE)
  # --------------------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------------------------------------
  # Save mbg input files
  
  # create mbg input files
  selected_covs <- data.table(covariate = selected_covs)
  cov_config <- merge(covs, selected_covs, by = 'covariate')
  
  # save mbg input files
  dir.create(paste0(mbg_input_dir, 'vif_', i, '/'), showWarnings = FALSE)
  write.csv(cov_config, paste0(mbg_input_dir, 'vif_', i, '/covs_', indicator_group, '_', region, '.csv'))
  # ----------------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------
  # End function
  
  # return a vector of covariates selected
  return(selection_results)
}
# ---------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
# Stepwise covariate removal using VIF

# loop over thresholds
for (i in threshold) {
  
  # report threshold
  message(paste0('\nThreshold: ', i))
  
  # load covariate data
  load(cropped_covs_file)
  
  # load input covariate names
  covs <- fread(paste0(core_repo, cov_file, cov_par, '.csv'))
  covs <- covs[include == TRUE]
  cov_names <- covs[, covariate]
  
  # set up table
  col_names <- c(cov_names, 'region', 'vif_threshold', 'num_covs_selected', 'data_run_date')
  selection_results <- setNames(data.table(matrix(nrow = 0, 
                                                  ncol = length(col_names))), 
                                col_names)
  
  # run vif selection analysis
  for (region in regions) {
    message(paste0('\n\n', region))
    
    # run stepwise vif selection function
    selection_results <- stepwise_vif_selection(thresh = i, 
                                                covariate_data = all_cov_data[[region]],
                                                reg = region,
                                                trace = T,
                                                not_to_drop = covs_to_keep)
  }
  
  # save selection results
  write.csv(selection_results, paste0(save_dir, 'selection_results_', indicator_group, '_vif_', i, '.csv'))
  
  # reset for next threshold
  message('\nSelection complete. Results saved.')
  
}
# ---------------------------------------------------------------------------------------------------------------

