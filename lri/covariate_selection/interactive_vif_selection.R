# --------------------------------------------------------------------------------------------------
# Script designed to facilitate covariate selection using variance inflation factor (VIF) analysis
# The purpose of the VIF analysis is to test for multicollinearity in covariates
#
# (1) load covariate data cropped to MBG input data by region
# (2) perform VIF analysis and drop covariates above a threshold by region
# (3) take results and summarize in a table that is human-readable
# --------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------------
# Setup

# clear environment
rm(list = ls())

# set arguments
user <- Sys.info()['user']
repo <- core_repo <- '<<<< FILEPATH REDACTED >>>>'
cov_repo <- paste0(repo, '/3_modeling/vif_selected_covs/')
config_file <- '3_modeling/config_ort_enet'
cov_file <- '3_modeling/covs_ort_standard'
indicator <- 'ors'
indicator_group <- 'ort'
run_date <- '2018_09_12_16_17_25'
age <- 0
holdout <- 0
threshold <- c(seq(2, 20, by = 0.5))
regions <- c('s2_afr_horn', 's2_cssa', 's2_wssa', 's2_name', 's2_sssa', 
             's2_mcacaf', 's2_s_america', 's2_central_asia', 's2_chn_mng', 
             's2_se_asia', 's2_malay', 's2_south_asia', 's2_mid_east', 's2_essa')

# indicate whether to use short list of covariates
short_list <- TRUE
select_cov_file <- '3_modeling/covs_ort_select'

# create directory to save results in
save_dir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(save_dir)

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>',header=FALSE)))
source(paste0(repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = repo)

## Read config file and save all parameters in memory
config <- load_config(repo            = core_repo,
                      indicator_group = '',
                      indicator       = '',
                      config_name     = config_file,
                      covs_name       = cov_file)
# ------------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Crop covariates to data

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
  gaul_list           <- get_gaul_codes(reg)
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
                        removeyemen = TRUE,
                        pathaddin   = pathaddin,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share),
                        yl          = year_list)
  
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
    
    effects <- trim(strsplit(gbd_fixed_effects, '\\+')[[1]])
    measures <- trim(strsplit(gbd_fixed_effects_measures, '\\+')[[1]])
    gbd_cov_layers <- load_gbd_covariates(covs     = effects,
                                          measures = measures,
                                          year_ids = year_list,
                                          age_ids  = gbd_fixed_effects_age,
                                          template = cov_layers[[1]][[1]],
                                          simple_polygon = simple_polygon,
                                          interval_mo = interval_mo)
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
save(all_cov_data, file = paste0(save_dir, 'covs_cropped_to_data_', Sys.Date(), '.RData'))
# ------------------------------------------------------------------------------------


# ---------------------------------------------------------------
# Stepwise VIF selection

# Start function
stepwise_vif_selection <- function(thresh, 
                                   covariate_data, 
                                   reg,
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
  
  # add threshold, region, and run date that corresponds to data
  dt[, region := reg]
  dt[, vif_threshold := thresh]
  dt[, data_run_date := run_date]
  
  # bind to original data table
  selection_results <- rbindlist(list(selection_results, dt), use.names = TRUE)
  # --------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------
  # End function
  
  # return a vector of covariates selected
  return(selection_results)
}
# ---------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Stepwise covariate removal using VIF

# loop over thresholds
for (i in threshold) {
  
  # report threshold
  message(paste0('\nThreshold: ', i))
  
  # load covariate data
  load(paste0(save_dir, 'covs_cropped_to_data_2018-09-28.RData'))
  
  # load input covariate names
  covs <- fread(paste0(core_repo, cov_file, '.csv'))
  original_covs <- covs[include == TRUE, covariate]
  
  # crop to just select covariates
  if (short_list) {
    # get select covariate names
    short_list_covs <- fread(paste0(core_repo, select_cov_file, '.csv'))
    short_list_names <- short_list_covs[include == TRUE, covariate]
    original_covs <- short_list_names
    # remove from covariate data
    regs <- names(all_cov_data)
    for (r in regs) {
      all_cov_data[[r]] <- all_cov_data[[r]][, short_list_names, with = FALSE]
    }
  }
  
  # set up for looping
  mincov <- 0
  cov_names <- original_covs

  # perform selection until all regions have VIFs less than the threshold for each covariate
  while (mincov < length(regions)) {
    
    # set up table
    col_names <- c(cov_names, 'region', 'vif_threshold', 'data_run_date')
    selection_results <- setNames(data.table(matrix(nrow = 0, 
                                                    ncol = length(col_names))), 
                                  col_names)
    
    # run vif selection analysis
    for (region in regions) {
      selection_results <- stepwise_vif_selection(thresh = i, 
                                                  covariate_data = all_cov_data[[region]],
                                                  reg = region,
                                                  trace = F)
    }
    
    # get first minimum value
    temp <- copy(selection_results)
    temp <- melt(temp, id.vars = 'region', measure.vars = cov_names)
    temp[, total := sum(value), by = 'variable']
    temp <- unique(temp[, c('region', 'value') := NULL])
    mincov <- temp[which.min(temp$total), total]
    mincov_name <- as.character(temp[which.min(temp$total), variable])
    
    # stop if all covariates selected
    if (mincov == length(regions)) break
      
    # remove the first minimum value covariate
    message(paste0('\n\nDropping ', mincov_name, ' which was selected ', mincov, ' times'))
    for (region in regions) {
      all_cov_data[[region]][, (mincov_name) := NULL]
    }
    
    # update covariate names
    cov_names <- names(all_cov_data[[1]])
    rm(selection_results)
    rm(temp)
    
  }
  
  # create and save covariate input file
  standard_covs <- fread('<<<< FILEPATH REDACTED >>>>')
  selected_covs <- data.table(covariate = cov_names)
  cov_config <- merge(standard_covs, selected_covs, by = 'covariate')
  write.csv(cov_config, paste0(save_dir, 'covs_', indicator_group, '_vif_', i, '.csv'))
  
  # reset for next threshold
  message('Selection complete. Results saved.')
  rm(selection_results)
  rm(temp)
  
}
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Summarize configs and remove duplicate selections

# load saved files
setwd(save_dir)
cov_files <- list.files(pattern = paste0('covs_', indicator_group, '_vif_'))
cov_list <- lapply(cov_files, fread)
vifs <- gsub('covs_ort_', '', cov_files)
vifs <- gsub('.csv', '', vifs)
names(cov_list) <- vifs

# create data table of covariate strings
cov_dt <- data.table(threshold = vifs, covariates = '')
for (i in vifs) {
  cov_dt[threshold == i, covariates := paste(cov_list[[i]][['covariate']], collapse = ', ')]
}
cov_dt[, threshold := as.numeric(gsub('vif_', '', threshold))]
setorderv(cov_dt, cols = 'threshold')

# find unique covariate sets
unique_covs <- copy(cov_dt)
unique_covs$threshold <- NULL
unique_covs <- unique(unique_covs)
unique_covs[, cov_group := .I]

# merge back into data table
cov_dt <- merge(cov_dt, unique_covs, by = 'covariates')
setorderv(cov_dt, cols = 'threshold')

# get threshold group names
cov_dt[, thresholds := paste(threshold, collapse = ', '), by = cov_group]
cov_dt <- unique(cov_dt[, threshold := NULL])
# ------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Save for upload into MBG model

# re-load standard cov list
standard_covs <- fread('<<<< FILEPATH REDACTED >>>>')

# loop over thresholds and save
i <- 1
for (i in 1:nrow(cov_dt)) {
  selected_covs <- data.table(covariate = strsplit(cov_dt$covariates[[i]], ', ')[[1]])
  cov_config <- merge(standard_covs, selected_covs, by = 'covariate')
  write.csv(cov_config, paste0(cov_repo, 'covs_', indicator_group, '_group_', cov_dt$cov_group[i], '.csv'))
  i <- i + 1
}

# save record of what thresholds correspond to what covariate groups
write.csv(cov_dt, paste0(cov_repo, 'cov_group_thresholds.csv'))
# ------------------------------------------------------------------------------------

