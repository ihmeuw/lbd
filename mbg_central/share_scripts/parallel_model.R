#####################################################################
## Generic parallel script for running MBG models                  ##
#####################################################################


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

## grab arguments
## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

## make a pathaddin that get used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

## load an image of the main environment
load(paste0('<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator,
            '/model_image_history/pre_run_tempimage_', run_date,
            pathaddin,'.RData'))

## In case anything got overwritten in the load, reload args
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s',indicator_group,indicator)
outputdir <- file.path('<<<< FILEPATH REDACTED >>>>>',indicator_group,
                       indicator,'output',run_date,'/')
dir.create(outputdir, showWarnings = FALSE)

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

## load packages and custom functions
## there is a difference depending on which nodes you are on
root        <- ifelse(Sys.info()[1]=="Windows", "<<<< FILEPATH REDACTED >>>>>", "<<<< FILEPATH REDACTED >>>>>")
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                         paste0(root,'temp/geospatial/geos_packages'),
                         paste0(root,'temp/geospatial/packages'))

commondir   <- sprintf('<<<< FILEPATH REDACTED >>>>>/mbg/common_inputs')
.libPaths(package_lib)
message(paste0('Loading packages from ',package_lib))
if(Sys.info()[1] == "Windows"){
  stop("STOP! you will overwrite these packages if you run from windows")
} else {
  for(package in c(package_list)) {
    message(paste0('loading ',package))
    library(package, lib.loc = package_lib, character.only=TRUE)
  }
}

## workaround for old INLA on geos nodes
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

# Ensure raster library is loaded and configure it.
(function() {
    require(raster)
    # set temporary file dir. INFR does not regularly delete files from here
    rasterOptions(tmpdir="<<<< FILEPATH REDACTED >>>>>/geospatial-tempfiles")
    # delete files older than 1 week. INFR may take over this duty at some point
    removeTmpFiles(h=24*7)
})()

# print out session info so we have it on record
sessionInfo()

## looks for any R script in main with name function
setwd(repo)
message(paste0('Loading local functions from ./mbg_central and ./',indicator_group,' ... '))
for(funk in list.files(recursive=TRUE,pattern='functions')){
  if(length(grep('central',funk))!=0){
    message(funk)
    source(funk)
  }
}

## Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
check_config()

## cores to use
cores_to_use <- round(as.numeric(slots)*.5)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PID <- Sys.getpid()
tic("Entire script") # Start master timer
psample_file <- paste0("psample_", indicator, "_", run_date, ".csv")
system(sprintf('dos2unix -n %s/mbg_central/share_scripts/psample %s/mbg_central/share_scripts/psample',repo,repo))
system(paste0(repo, '/mbg_central/share_scripts/psample --showthreads --pid=',PID,' ',
              sprintf('/share/geospatial/mbg/%s/%s/output/%s/%s ',
                      indicator_group, indicator, run_date, psample_file)),
       wait=FALSE) ## start profiling script

## skip a large chunk if requested in config
if(as.logical(skiptoinla) == FALSE){
message('You have chosen to not skip directly to inla.')

## some set up
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

## Load simple polygon template to model over
gaul_list           <- get_gaul_codes(reg)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, use_premade = T)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']]

## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
# For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
if(holdout!=0) {
  message(paste0('Holdout != 0 so loading holdout data only from holdout ',holdout))
  message('Please be sure you have a list object called stratum_ho in your environment.')
  # if strateifies by age then make sure loads correctly
  if(age!=0) df <- as.data.table(stratum_ho[[paste('region',reg,'_age',age,sep='__')]])
  if(age==0) df <- as.data.table(stratum_ho[[paste('region',reg,sep='__')]])
  df <- df[fold != holdout, ]
}
if(holdout==0) {
  message('Holdout == 0 so loading in full dataset using load_input_data()')
  df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                        simple      = simple_polygon,
                        agebin      = age,
                        removeyemen = TRUE,
                        pathaddin   = pathaddin,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share))
}

## if testing, we only keep 1000 or so observations
if(test == 1){
  test_pct <- as.numeric(test_pct)
  
  message(paste0('Test option was set on and the test_pct argument was found at ',test_pct,'% \n
                 ... keeping only ', round(nrow(df)*(test_pct/100),0),' random rows of data.'))
  df <- df[sample(nrow(df),round(nrow(df)*(test_pct/100),0)),]
  
  message('Also, making it so we only take 100 draws')
  samples <- 100
}


 ## for u5m in particular,  make sure indicator is properly named here, wont affect others
  df[[indicator]] <- df[[gsub(paste0('_age',age),'',indicator)]]

  ## if there is another weight column, multiply it with weight now
  if(exists('other_weight')) if(other_weight!='') {
                               message(paste0('Multiplying weight and ',other_weight))
                               df[['weight']] <- df[['weight']]*df[[other_weight]]
                             }

  ## Some built in data checks that cause known problems later on
  if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
    stop('You have binomial data where k > N. Check your data before proceeding')
  if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
    stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')

  ## make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

  ## Save distribution of data for this region
  png(paste0(outputdir, reg, '.png'))
  if(indicator_family=='binomial') hist(df[, get(indicator)]/df$N) else hist(df[, get(indicator)])
  dev.off()

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ## Define modeling space. In years only for now.
  if(yearload=='annual') period_map <-
                         make_period_map(modeling_periods = c(min(year_list):max(year_list)))
  if(yearload=='five-year') period_map <-
                            make_period_map(modeling_periods = seq(min(year_list),max(year_list), by=5))

  ## make covariates conditional
  cov_layers <- gbd_cov_layers <- mbg_cov_layers <- NULL

  ## Pull all covariate bricks/layers
  if(nchar(fixed_effects)> 0){
    message('Grabbing raster covariate layers')
    selected_fixed_effects <- strsplit(fixed_effects," ")[[1]][strsplit(fixed_effects," ")[[1]] != "+"]
    selected_measures     <- strsplit(fixed_effects_measures," ")[[1]][strsplit(fixed_effects_measures," ")[[1]] != "+"]
    cov_layers <- load_and_crop_covariates_annual(
      covs            = selected_fixed_effects,
      measures        = selected_measures,
      simple_polygon  = simple_polygon,
      start_year      = min(year_list),
      end_year        = max(year_list),
      interval_mo     = as.numeric(interval_mo) )
  }

  ## pull country level gbd covariates
  if(nchar(gbd_fixed_effects)>0){
    message('Grabbing GBD covariates')
    selected_gbd_fixed_effects <-
    strsplit(gbd_fixed_effects," ")[[1]][strsplit(gbd_fixed_effects," ")[[1]] != "+"]

    ## get full gaul list, including buffer zone
    gbd_gaul_list <- unique(c(gaul_convert(unique(df[, country])), gaul_list))
                                        # rasterize gbd info
    gbd_cov_layers <- load_gbd_covariates(gbd_fixed_effects = selected_gbd_fixed_effects,
                                          year_ids          = year_list,
                                          gaul_list         = gbd_gaul_list,
                                          template          = cov_layers[[1]][[1]])
                                        # combine raster and gbd covariates
    all_cov_layers <- c(cov_layers, gbd_cov_layers)
  } else {
    all_cov_layers <- cov_layers
  }


  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
  ## regenerate all fixed effects equation from the cov layers
  all_fixed_effects <- all_raw_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")

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
                                         collapse = " + "), sep=" + ")

    ## update specific stacker formulas (for now we just want country effects in BRT)
    all_fixed_effects_brt <- all_fixed_effects_cfes
  }

  ## Add these to the fixed effects if we want them in stacking
  if(use_child_country_fes == TRUE) {
    gaul_fes <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
    all_fixed_effects = paste(all_fixed_effects, gaul_fes, sep = " + ")
  }


  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Stacking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  tic("Stacking - all") ## Start stacking master timer

  ## Figure out which models we're going to use
  child_model_names <- stacked_fixed_effects        %>%
  gsub(" ", "", .)          %>%
  strsplit(., "+", fixed=T) %>%
  unlist
  message(paste0('Child stackers included are: ',paste(child_model_names,collapse=' // ')))

  the_covs <- format_covariates(all_fixed_effects)
  
  ## copy the dataset to avoid unintended namespace conflicts
  the_data <- copy(df)
  
  ## shuffle the data into six folds
  the_data <- the_data[sample(nrow(the_data)),]
  the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=as.numeric(n_stack_folds),labels=FALSE)]

  ## add a row id column
  the_data[, a_rowid := seq(1:nrow(the_data))]
  
  ## extract covariates to the points and subset data where its missing covariate values
  cs_covs <- extract_covariates(the_data,
                                all_cov_layers,
                                id_col              = "a_rowid",
                                return_only_results = TRUE,
                                centre_scale        = TRUE,
                                period_var           = 'year',
                                period_map          = period_map)
  
  ## Check for data where covariate extraction failed
  rows_missing_covs <- nrow(the_data) - nrow(cs_covs[[1]])
  if (rows_missing_covs > 0) {
    pct_missing_covs <- round((rows_missing_covs/nrow(the_data))*100, 2)
    warning(paste0(rows_missing_covs, " out of ", nrow(the_data), " rows of data ",
                   "(", pct_missing_covs, "%) do not have corresponding ", 
                   "covariate values and will be dropped from child models..."))
    if (rows_missing_covs/nrow(the_data) > 0.1) {
      stop(paste0("Something has gone quite wrong: more than 10% of your data does not have ", 
                  "corresponding covariates.  You should investigate this before proceeding."))
    }
  }

  ## merge back on the extractes covariates
  the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)               

  ## store the centre scaling mapping
  covs_cs_df  <-  cs_covs[[2]]
  
  ## this will drop rows with NA covariate values
  the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))
  
  if(use_stacking_covs){
   

    ## seperated out into a different script
    message('Fitting Stackers')
    source('./mbg_central/share_scripts/run_child_stackers.R')

    ## combine the children models
    the_data  <- cbind(the_data, do.call(cbind, lapply(lapply(child_model_names, 'get'), function(x) x[[1]])))
    child_model_objs <- setNames(lapply(lapply(child_model_names, 'get'), function(x) x[[2]]), child_model_names)

    ## fit GAM stacker (just for fun? Daniel?)
    stacked_results <- gam_stacker(the_data,
                                   model_names      = child_model_names,
                                   indicator        = indicator,
                                   indicator_family = indicator_family)

    ## return the stacked rasters
    stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                          period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                          child_models     = child_model_objs,
                                          stacker_model    = stacked_results[[2]],
                                          indicator_family = indicator_family,
                                          return_children  = TRUE,
                                          centre_scale_df  = covs_cs_df)

    ## plot stackers
    pdf(paste0(outputdir, 'stacker_rasters', pathaddin, '.pdf'))
    for(i in 1:length(stacked_rasters))
      plot(stacked_rasters[[i]],main=names(stacked_rasters[[i]]),maxpixel=ncell(stacked_rasters[[i]]))
    dev.off()

    message('Stacking is complete')
  }

  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## copy things back over to df
  df <- copy(the_data)

  ## remove the covariate columns so that there are no name conflicts when they get added back in
  df <- df[,paste0(the_covs) := rep(NULL, length(the_covs))]

  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
  df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]

  
  ## create a full raster list to carry though to the shiny/next steps
  if(use_stacking_covs){
    cov_list         <- c(unlist(stacked_rasters), unlist(all_cov_layers))
    child_mod_ras    <- cov_list[child_model_names]
  }else{
    cov_list <- unlist(all_cov_layers)
  }
  

  toc(log = T) ## End stacking master timer

  ## make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d           = df,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset)

  ## Build temporal mesh (standard for now)
  mesh_t <- build_time_mesh(periods=eval(parse(text=mesh_t_knots)))


}else{ ## END !SKIPTOINLA
  message(paste0('You have chosen to skip directly to INLA. Picking up objects from run_date ',
                 skiptoinla_from_rundate))
  message('Now loading saved MBG inputs from that chosen run_date.')

  load(paste0('<<<<< FILEPATH REDACTED >>>>>/', indicator_group, '/',
              indicator, '/model_image_history/',
              skiptoinla_from_rundate, pathaddin, '.RData'))  
}


## setup our fixed effects 
if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
  all_fixed_effects <- ''
}
if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
  all_fixed_effects <- stacked_fixed_effects
}
if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
  all_fixed_effects <- fixed_effects ## from config
}
if(!as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
  all_fixed_effects <- paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + ")
}
if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & !as.logical(use_inla_country_fes)){
  all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, sep = " + ")
}
if(as.logical(use_stacking_covs) & !as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
  all_fixed_effects <- paste(stacked_fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
}
if(!as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
  all_fixed_effects <- paste(fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
}
if(as.logical(use_stacking_covs) & as.logical(use_raw_covs) & as.logical(use_inla_country_fes)){
  all_fixed_effects <- paste(stacked_fixed_effects, fixed_effects, paste(names(gaul_code)[2:length(names(gaul_code))], collapse = " + "), sep = " + ")
}

## Save all inputs for MBG model into correct location on /share
if(as.logical(use_stacking_covs) & !as.logical(skiptoinla)){
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 mesh_s            = mesh_s,
                 mesh_t            = mesh_t,
                 cov_list          = cov_list,
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = TRUE)  
}else if(!as.logical(skiptoinla)){
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 mesh_s            = mesh_s,
                 mesh_t            = mesh_t,
                 cov_list          = cov_list,
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = '',
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = FALSE)
}


## reload data an prepare for MBG
if(!as.logical(skiptoinla)){

  load(paste0('<<<<< FILEPATH REDACTED >>>>>', indicator_group, '/', indicator,
              '/model_image_history/', run_date, pathaddin, '.RData'))
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


tic("MBG - all") # Start MBG master timer

## for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
if(as.logical(use_stacking_covs)){
  df <- df[,paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]
}

if(as.logical(use_stacking_covs)){ all_fixed_effects <- 'gam + gbm + lasso'}

## Generate MBG formula for INLA call
mbg_formula <- build_mbg_formula_with_priors(fixed_effects   = all_fixed_effects,
                                             add_nugget      = use_inla_nugget,
                                             add_ctry_res    = use_inla_country_res,
                                             no_gp = !as.logical(use_gp),
                                             coefs.sum1 = as.logical(coefs.sum1))

if(constrain_children_0_inf)
  mbg_formula <- build_mbg_formula_with_priors(fixed_effects = all_fixed_effects,
                                               add_nugget    = use_inla_nugget,
                                               positive_constrained_variables = child_model_names,
                                               no_gp = !as.logical(use_gp))

## we may want to use covariates in logit space un-centre-scaled
if(as.logical(coefs.sum1) | as.logical(logit.stack.cov)){
  message("using stackers in logit space and non-center-scaled")
  
  f.e.v <- base::strsplit(all_fixed_effects, ' + ', fixed = TRUE)[[1]]

  df[, (f.e.v) := lapply(.SD, 'logit'), .SDcols = f.e.v]    
  
  ## Create SPDE INLA stack
  input_data <- build_mbg_data_stack(df            = df,
                                     fixed_effects = all_fixed_effects,
                                     mesh_s        = mesh_s,
                                     mesh_t        = mesh_t,
                                     use_ctry_res  = use_inla_country_res,
                                     use_nugget    = use_inla_nugget,
                                     coefs.sum1    = as.logical(coefs.sum1), 
                                     exclude_cs    = f.e.v) ## dont center-scale
}else{
  ## Create SPDE INLA stack
  input_data <- build_mbg_data_stack(df            = df,
                                     fixed_effects = all_fixed_effects,
                                     mesh_s        = mesh_s,
                                     mesh_t        = mesh_t,
                                     use_ctry_res  = use_inla_country_res,
                                     use_nugget    = use_inla_nugget,
                                     coefs.sum1    = as.logical(coefs.sum1))
}

## combine all the inputs
stacked_input <- input_data[[1]]
spde          <- input_data[[2]]
cs_df         <- input_data[[3]]

## save spde object (useful later for transforming INLA hyperparameters)
saveRDS(object = spde,
        file   = sprintf('%s/output/%s/inla_spde_object_%s_holdout_%i.RDS',
                         sharedir, run_date, reg, holdout))

## Generate other inputs necessary
outcome <- df[[indicator]] ## N+_i - event obs in cluster
N       <- df$N            ## N_i - total obs in cluster
weights <- df$weight

## catch in case there is no weight column
if(is.null(weights)){
  weights = rep(1,nrow(df))
}

tic("MBG - fit model") ## Start MBG - model fit timer

## Fit MBG model
source('mbg_central/mbg_functions.R')

if(!as.logical(skipinla)){
  model_fit <- fit_mbg(indicator_family = indicator_family,
                       stack.obs        = stacked_input,
                       spde             = spde,
                       cov              = outcome,
                       N                = N,
                       int_prior_mn     = intercept_prior,
                       f_mbg            = mbg_formula,
                       run_date         = run_date,
                       keep_inla_files  = keep_inla_files,
                       cores            = cores_to_use,
                       wgts             = weights,
                       verbose          = TRUE, 
                       intstrat         = intstrat,
                       fe_sd_prior      = 1 / 9) ## this actually sets precision!. prec=1/9 -> sd=3


  saveRDS(object = model_fit, ## save this here in case predict dies
          file = sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/inla_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                         indicator_group, indicator, run_date, reg, holdout, age))
}else{
  ## skipped fitting INLA so just load model and move to predict
  model_fit <- readRDS( file = sprintf('<<<< FILEPATH REDACTED >>>>>/%s/%s/output/%s/inla_model_fit_pre_preds_%s_holdout_%i_agebin_%i.RDS',
                         indicator_group, indicator, run_date, reg, holdout, age))
}

toc(log = T) # End MBG - model fit timer

tic("MBG - predict model") # Start MBG - model predict timer

# Run predict_mbg on chunks of 50 samples (to avoid memory issues)

message('Making predictions in 50 draw chunks.')

max_chunk <- 50
samples   <- as.numeric(samples)

## Create vector of chunk sizes
chunks <- rep(max_chunk, samples %/% max_chunk)
if (samples %% max_chunk > 0) chunks <- c(chunks, samples %% max_chunk)
pm <- lapply(chunks, function(samp) {
  predict_mbg(res_fit       = model_fit,
              cs_df         = cs_df,
              mesh_s        = mesh_s,
              mesh_t        = mesh_t,
              cov_list      = cov_list,
              samples       = samp,
              simple_raster = simple_raster,
              transform     = transform,
              coefs.sum1    = as.logical(coefs.sum1),
              logit.covs    = as.logical(logit.stack.cov), 
              pred_gp       = as.logical(use_gp)
              )[[3]]
})

## Make cell preds and a mean raster
cell_pred <- do.call(cbind, pm)
mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),
                                               ncol = max(period_map$period)))

toc(log = T) ## Stop MBG - model predict timer


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message('Wrapping up')
save_mbg_preds(config     = config,
               time_stamp = time_stamp,
               run_date   = run_date,
               mean_ras   = mean_ras,
               sd_ras     = NULL,
               res_fit    = model_fit,
               cell_pred  = cell_pred,
               df         = df,
               pathaddin  = pathaddin)

## plot the mean raster
pdf(paste0(outputdir,'mean_raster', pathaddin, '.pdf'))
plot(mean_ras,maxpixel=ncell(mean_ras))
dev.off()

# shut off sampling script
fname <- paste0("$(hostname).", PID, ".pid")
cmd <- paste0("kill $(cat $HOME/", fname, ")")
system(cmd)

## timer stuff
toc(log = T) # End master timer

## Format timer
ticlog   <- tic.log(format = F)
df_timer <- generate_time_log(ticlog)
df_timer[, region := reg]
df_timer[, holdout := holdout]
setcolorder(df_timer, c("region", "holdout", "step", "time"))

## Pull run time for this run
run_time_all <- df_timer[step == "Entire script", time]

## Write to a run summary csv file in the output directory
output_file <- paste0(outputdir, "run_summary_", indicator, "_", run_date, ".csv")

## Pull in file contents from other reg/holdouts (if exists)
if (file.exists(output_file)) {
  file_contents <- read.csv(output_file, stringsAsFactors = F) %>% as.data.table
  df_timer      <- rbind(file_contents, df_timer)
}

## Write CSV
write.csv(df_timer, file = output_file, row.names = FALSE)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
