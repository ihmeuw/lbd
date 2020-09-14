#####################################################################
# MBG launch script: 1) launch MBG models
# 2) rake to incidence, prevalence, mortality
#####################################################################

#--------------------------------------------------------------------
# (1) Launch MBG models
#--------------------------------------------------------------------

# Setup -------------------------------------------------------------

## set arguments
user                <- commandArgs()[4]
use_run_date        <- commandArgs()[5]
indicator_group     <- commandArgs()[6]
indicator           <- commandArgs()[7]
core_repo           <- commandArgs()[8]
indicator_repo      <- commandArgs()[9]
config_par          <- commandArgs()[10]
cov_par             <- commandArgs()[11]
gbm_par             <- commandArgs()[12]
region_specific_cov <- commandArgs()[13]
Regions             <- unlist(strsplit(commandArgs()[14], '~'))
go_to_inla          <- commandArgs()[15]
rake_measures       <- unlist(strsplit(commandArgs()[16], '~'))
skip_to_rake        <- commandArgs()[17]
etiology            <- commandArgs()[18]
priority            <- commandArgs()[19]
go_to_inla_from_rundate <- commandArgs()[20]
make_holdouts       <- commandArgs()[21]
holdout_type    <- commandArgs()[22]

## drive locations
sharedir       <- '<<<< FILEPATH REDACTED >>>>'
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

## Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)
library(mgcv)

## source any functions with custom edits from lbd_custom folder
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/misc_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/post_estimation_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/stacking_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/fractional_raking_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/prep_functions.R')
source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')

## Set up config
config <- set_up_config(repo = '<<<< FILEPATH REDACTED >>>>',
                        core_repo = core_repo,
                        indicator_group,
                        indicator,
                        config_name = paste0('config_', config_par),
                        covs_name = paste0('covs_', cov_par),
                        run_tests = FALSE)

## Set project
proj <- ifelse(as.logical(use_geos_nodes), 'proj_geo_nodes', 'proj_geospatial')

## Create run date in correct format
if (is.na(use_run_date) | use_run_date == 'NA'){
  run_date <- make_time_stamp(TRUE)
} else {
  run_date <- use_run_date
}

## Create output folder with the run_date
outputdir <- '<<<< FILEPATH REDACTED >>>>'

dir.create(outputdir)

## Make sure year object is in the correct format
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))


## Set optimized BRT parameters if using GBM
if (stacked_fixed_effects %like% 'gbm') {
  gbm_params <- fread('<<<< FILEPATH REDACTED >>>>', stringsAsFactors = F)
} else {
  gbm_cv <- NA
  gbm_tc <- NA
  gbm_bf <- NA
}


## Record model parameters in a google sheet model tracker
if (is.na(use_run_date) | use_run_date == 'NA'){
  library('googlesheets')
  region_as_char <- paste(Regions, sep=" ", collapse=" + ")
  if (stacked_fixed_effects %like% 'gbm') {
    gbm_bf <- paste(gbm_params$gbm_bf, sep=" ", collapse=" + ")
    gbm_cv <- paste(gbm_params$gbm_cv, sep=" ", collapse=" + ")
    gbm_tc <- paste(gbm_params$gbm_tc, sep=" ", collapse=" + ")
  } else {
    gbm_bf <- 'NA'
    gbm_cv <- 'NA'
    gbm_tc <- 'NA'
    gbm_par <- 'NA'
  }
  model_params <- c(indicator_group, indicator, run_date, outputdir, region_as_char, 
                    config_par, cov_par, gbm_par, region_specific_cov, fixed_effects, gbd_fixed_effects,
                    stacked_fixed_effects, makeholdouts, use_inla_country_fes, gbm_cv, gbm_bf, gbm_tc,
                    mesh_t_knots, rho_prior, nugget_prior, ctry_re_prior, modeling_shapefile_version)
  
  ttt <- readRDS('<<<< FILEPATH REDACTED >>>>')
  gs_auth(token = ttt)
  lri_tracker <- gs_title('LRI model tracker')
  gs_add_row(lri_tracker, ws = 'model runs', input = model_params)
  rm(lri_tracker)
}else if(is.null(gbm_par)){
  gbm_par <- 'NA'
}

## Make holdouts -------------------------------------------------------------------------

if(as.logical(make_holdouts) & !as.logical(skiptoinla)){
  message('Making holdouts')
  
  set.seed(98112)
  
  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        simple      = NULL,
                        removeyemen = FALSE,
                        years       = yearload,
                        yl          = year_list,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share))
  
  # add in location information
  df <- merge_with_ihme_loc(df, shapefile_version = modeling_shapefile_version)
  
  # remove data for countries outside of region
  df <- df[!is.na(region)]
  
  # make holdouts based on nids, if selected
  if (holdout_type == 'nids') {
    
    # make sure we have enough nids to do holdouts
    if (!is.logical(skipinla)) {
      if (length(unique(df$nid)) < 5) {
        n_ho_folds <- length(unique(df$nid))
      }
    }
    
    # make a list of dfs for each region, with 5 nid folds identified in each
    stratum_ho <-   make_folds(data       = df,
                               n_folds    = as.numeric(n_ho_folds),
                               spte_strat = 'nids',
                               strat_cols = 'region',
                               ss_col     = ss_col,
                               yr_col     = yr_col,
                               seed       = 98112,
                               save.file = '<<<< FILEPATH REDACTED >>>>')
  } 
  
  # make holdouts based on admin 1, if selected
  if (holdout_type == 'admin') {
    
    # Load simple polygon template by region
    gaul_list           <- get_adm0_codes(Regions, shapefile_version = modeling_shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
    subset_shape        <- simple_polygon_list[[1]]
    
    # Load simple raster by region
    raster_list        <- build_simple_raster_pop(subset_shape)
    simple_raster      <- raster_list[['simple_raster']]
    
    # load admin 2 shapefile and crop by region
    shapefile_admin <- shapefile(get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))
    shapefile_admin <- gBuffer(shapefile_admin, byid = TRUE, width = 0)
    shapefile_admin <- crop(shapefile_admin, extent(subset_shape))
    
    # load mask raster and crop by region
    raster_mask <- raster('<<<< FILEPATH REDACTED >>>>')
    raster_mask <- crop(raster_mask, extent(simple_raster))
    
    # make admin 2 holdouts
    stratum_ho <- make_folds(data = df,
                             n_folds = as.numeric(n_ho_folds),
                             spat_strat = 'poly',
                             temp_strat = 'prop',
                             strat_cols = 'region',
                             ts = as.numeric(ho_ts),
                             mb = as.numeric(ho_mb),
                             admin_shps = shapefile_admin,
                             admin_raster = simple_raster,
                             mask_shape = subset_shape,
                             mask_raster = raster_mask,
                             shape_ident = 'ADM2_CODE',
                             lat_col = lat_col,
                             long_col = long_col,
                             ss_col = ss_col,
                             yr_col = yr_col,
                             seed = 98112,
                             save.file = '<<<< FILEPATH REDACTED >>>>')
  } 
}

## Launch parallel script -------------------------------------------------------------------------

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(make_holdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)

if (skip_to_rake == FALSE){
  ## loop over strata, save images and submit qsubs
  
  ## set region-specific covariates if needed
  for(i in 1:nrow(loopvars)){
    
    if (region_specific_cov){
      reg <- loopvars[i,1]
      cov_par_reg <- paste0('covs_', cov_par, '_', reg)
      
      config <- set_up_config(repo = '<<<< FILEPATH REDACTED >>>>',
                              core_repo = core_repo,
                              indicator_group,
                              indicator,
                              config_name = paste0('config_', config_par),
                              covs_name = paste0('covs_', cov_par, '_', reg),
                              run_tests = FALSE)
      
      ##double check pop
      if (pop_measure != 'a0004t'){
        print('check population measure!')
      }
      
      message(paste('Used region-specific covariate file', cov_par_reg, 'for region', reg))
    }
    
    if (gbm_par != 'NA'){
      reg <- loopvars[i,1]
      gbm_tc <- gbm_params[region == reg, gbm_tc]
      gbm_lr <- gbm_params[region == reg, gbm_lr]
      gbm_bf <- gbm_params[region == reg, gbm_bf]
      gbm_nminobs <- gbm_params[region == reg, gbm_nminobs]
      gbm_ntrees <- gbm_params[region == reg, gbm_ntrees]
      gbm_cv <- gbm_params[region == reg, gbm_cv]
      
      message(paste0('Used region-specific gbm parameters ', 'gbm_params_', gbm_par, '.csv for region ', reg))
    }
    
    # set cores by region
    region_cores <- 8
    reg <- loopvars[i,1]
    if(reg == 'dia_central_asia') region_cores <- 4
    if(reg == 'sssa' | reg == 'dia_se_asia' | reg == 'dia_sssa') region_cores <- 6
    if(reg == 'dia_malay' | reg == 'dia_name' | reg == 'name') region_cores <- 10
    if(reg == 'dia_wssa' | reg =='dia_south_asia' | reg == 'wssa' | reg == 'ansa' | reg == 'caca' | reg == 'dia_afr_horn') region_cores <- 12
    if(reg == 'dia_chn_mng') region_cores <- 40
    if(reg == 'dia_s_america-GUF') region_cores <- 12
    
    # convert cores approximately to memory and run time
    region_rt <- '06:00:00:00'
    if (region_cores < 9) region_rt <- '03:00:00:00'
    if (region_cores < 6) region_rt <- '01:12:00:00'
    if (region_cores > 9) region_rt <- '16:00:00:00'
    region_mem <- region_cores*15
    
    if (go_to_inla == TRUE){
      skiptoinla <- TRUE
      skiptoinla_from_rundate <- go_to_inla_from_rundate
    }
    
    #message age, region, holdout
    message(paste(loopvars[i,2],as.character(loopvars[i,1]),loopvars[i,3]))
    
    # make a qsub string
    qsub <- make_qsub_share(age           = loopvars[i,2],
                            reg           = as.character(loopvars[i,1]),
                            holdout       = loopvars[i,3],
                            test          = F,
                            indic         = indicator,
                            user          = indicator_group,
                            rd            = run_date,
                            saveimage     = TRUE,
                            use_c2_nodes  = FALSE,
                            memory        = region_mem,
                            cores         = region_cores,
                            run_time      = region_rt,
                            proj          = proj,
                            geo_nodes     = as.logical(use_geos_nodes),
                            corerepo      = core_repo,
                            code_path          = '<<<< FILEPATH REDACTED >>>>', # custom parallel model
                            addl_job_name = paste0(indicator, '_parallel'),
                            singularity   = '<<<< FILEPATH REDACTED >>>>')
        # submit job
    system(paste0(qsub, ' -p ', priority))
  }
  
  waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)
  
} else if (skip_to_rake == TRUE){
  message(paste0('Skipping to raking without launching new models for run_date ', run_date))
}

#--------------------------------------------------------------------
# (2) Rake to GBD using fractional raking
#--------------------------------------------------------------------
#loop over regions, holdouts, measures
for(i in 1:nrow(loopvars)){
  for (reg in loopvars[i,1]) {
    for (holdout in loopvars[i,3]) {
      for (measure in rake_measures) {
        
        # Make qsub string and submit -----------------------------------------
        mem <- '200G'
        if (reg == 'dia_chn_mng' | reg == 'dia_s_america' | reg == 'dia_s_america-GUF') mem <- '300G'
        rt <- '04:00:00'
        queue <- 'geospatial.q'
        proj <- 'proj_geo_nodes'
        name <- paste0('rake_',reg, '_', measure, '_', holdout)
        shell <- '<<<< FILEPATH REDACTED >>>>'
        code <- '<<<< FILEPATH REDACTED >>>>'
        
        args <- paste(user,
                      core_repo,
                      indicator_group,
                      indicator,
                      config_par,
                      cov_par,
                      reg,
                      run_date,
                      measure,
                      holdout,
                      etiology)
        
        qsub <- paste0('qsub -l m_mem_free=', mem,
                       ' -l fthread=1 -l h_rt=', rt, 
                       ' -v sing_image=default -q ', queue, 
                       ' -P ', proj,
                       ' -p ', priority,
                       ' -N ', name, 
                       ' ', shell, 
                       ' ', code, 
                       ' ', args)
        
        system(qsub)
        
      }
    }
  }
}

