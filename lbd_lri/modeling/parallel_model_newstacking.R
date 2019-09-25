## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# grab arguments
indicator_group          <- as.character(commandArgs()[4])
indicator                <- as.character(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
task_id                  <- as.numeric(commandArgs()[7])

load('<<<< FILEPATH REDACTED >>>>') #prerun temp image

#reload run_date
loopvars = loopvars[loopvars$id == task_id,]
run_date <- as.character(commandArgs()[6])
reg = loopvars$reg
age = loopvars$age
holdout = loopvars$holdout
test = loopvars$test
start_point = loopvars$start_point

#pull out the permutation, if using.
if(loopvars$permutations!=""){
  permutation = loopvars$permutations[loopvars$id == task_id]
  permute = T
  run_date = paste0(run_date,'_', loopvars$permutations)
}else{
  permute = F
  permutation = ""
}

path_to_save_file <- '<<<< FILEPATH REDACTED >>>>'
# echo $! writes the PID of the process out. intern=TRUE causes the PID to be returned as a string
meminfo_pid <- system('<<<< FILEPATH REDACTED >>>>', intern=TRUE)

useelogit = T
message(useelogit)

#sort out constraints and things
coefssum1 = as.logical(coefssum1)
coefsmore0 = as.logical(coefsmore0)

message(paste(coefssum1, coefsmore0))

#create pathaddin
pathaddin <- '<<<< FILEPATH REDACTED >>>>'
outputdir <- '<<<< FILEPATH REDACTED >>>>'

# print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir','permutation', 'coefssum1','coefsmore0', 'start_point'))
  message(paste0(arg, ' = ',get(arg)))

## load packages and custom functions
# there is a difference depending on which nodes you are on
root        <- '<<<< FILEPATH REDACTED >>>>'
package_lib    <- '<<<< FILEPATH REDACTED >>>>'

commondir   <- '<<<< FILEPATH REDACTED >>>>'

package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
if(Sys.info()[1] == "Windows"){
  stop("STOP! you will overwrite these packages if you run from windows")
} else {
  for(package in c(package_list)) {
    suppressMessages(suppressWarnings(library(package, character.only=TRUE)))
  }
}
.libPaths(package_lib)
library('mbgstacking')
# Ensure raster library is loaded and configure it.
rasterOptions(tmpdir='<<<< FILEPATH REDACTED >>>>')
# delete files older than 1 week.
removeTmpFiles(h=24*7)


# looks for any R script in main with name function
user            <- Sys.info()['user']
repo            <- sprintf('<<<< FILEPATH REDACTED >>>>')
for(funk in list.files(repo, recursive=TRUE,pattern='functions', full.names = T)){
  if(length(grep('central',funk))!=0){
    #message(funk)
    source(funk)
  }
}

# cores to use
cores_to_use <- as.numeric(Sys.getenv('NSLOTS')) -1
if(is.na(cores_to_use)) cores_to_use = 5
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(start_point <= 1){
  # skip a large chunk if requested in config

  ## some set up
  year_list    <- eval(parse(text=year_list))

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

  #drop outliered data
  if(file.exists('<<<< FILEPATH REDACTED >>>>')){
    outliers = read.csv('<<<< FILEPATH REDACTED >>>>')
    df = df[!nid %in% outliers$nid,]
  }


  #drop years outside the period of interest
  df = df[year<=max(year_list) & year>=min(year_list), ]

  # for u5m in particular,  make sure indicator is properly named here, wont affect others
  df[[indicator]] <- df[[gsub(paste0('_age',age),'',indicator)]]

  # if there is another weight column, multiply it with weight now
  if(exists('other_weight')) if(other_weight!='') {
    message(paste0('Multiplying weight and ',other_weight))
    df[['weight']] <- df[['weight']]*df[[other_weight]]
  }

  df[, alt_weight := weight * N]

  ## Some built in data checks that cause known problems later on
  if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
    stop('You have binomial data where k > N. Check your data before proceeding')
  if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
    stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')

  # make sure this inla patch is implemented if running on geos
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
    make_period_map(modeling_periods = seq(min(year_list),max(year_list),by=5))

  # make covariates conditional
  cov_layers <- gbd_cov_layers <- mbg_cov_layers <- NULL

  # Pull all covariate bricks/layers
  if(nchar(fixed_effects)> 0){
    message('Grabbing raster covariate layers')
    selected_fixed_effects <-
      strsplit(fixed_effects," ")[[1]][strsplit(fixed_effects," ")[[1]] != "+"]
    selected_measures      <-
      strsplit(fixed_effects_measures," ")[[1]][strsplit(fixed_effects_measures," ")[[1]] != "+"]
    cov_layers <- load_and_crop_covariates_annual(
      covs            = selected_fixed_effects,
      measures        = selected_measures,
      simple_polygon  = simple_polygon,
      start_year      = min(year_list),
      end_year        = max(year_list),
      interval_mo     = as.numeric(interval_mo))
  }

  # pull country level gbd covariates (assuming not running a single country model)
  if(nchar(gbd_fixed_effects)>0 & length(gaul_list)>1){
    message('Grabbing GBD covariates')
    selected_gbd_fixed_effects <-
      strsplit(gbd_fixed_effects," ")[[1]][strsplit(gbd_fixed_effects," ")[[1]] != "+"]
    # get full gaul list, including buffer zone
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

  # regenerate all fixed effects equation from the cov layers
  all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")


  #save checkpoint
  save.image('<<<< FILEPATH REDACTED >>>>')
} else{
  sp = start_point
  hold_rd = run_date
  hold_permute = permute
  hold_permutation = permutation
  load('<<<< FILEPATH REDACTED >>>>')
  start_point =  sp
  run_date = hold_rd
  permute = hold_permute
  permutation = hold_permutation
}



#country value
loc = 47 #cape verde

if(loc %in% unique(simple_raster) & 'stunting_mod_b_raked' %in% names(all_cov_layers)){

  #add country levels
  fixme = all_cov_layers[['stunting_mod_b_raked']]
  simp_ply_ras = fixme[[1]]
  simp_ply_ras[] = 1:length(simp_ply_ras)
  simp_ply_ras_crop = crop(simp_ply_ras,simple_raster)

  val_list = read.csv('<<<< FILEPATH REDACTED >>>>')
  setDT(val_list)
  setorder(val_list, +year_id)
  val_list = val_list[location_name == 'Cape Verde', mean]

  all_cov_layers[['stunting_mod_b_raked']] = brick(lapply(seq(year_list), function(x){

    #replace all cells in the covariate layer as corrosponding to simple raster
    ras = all_cov_layers[['stunting_mod_b_raked']][[x]]

    ras[simp_ply_ras_crop[simple_raster[] == loc]] = val_list[x]
    return(ras)

  }))
}


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Stacking ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(start_point <= 2){
  tic("Stacking - all") # Start stacking master timer
  wf = '<<<< FILEPATH REDACTED >>>>'
  dir.create(wf, recursive = T)
  logloc = '<<<< FILEPATH REDACTED >>>>'
  oot = paste0(logloc, 'output/')
  eet = paste0(logloc, 'errors/')

  #check for omp
  if(!grepl('SINGULARITYENV_OMP_NUM_THREADS', r_path)){
    r_path = paste('SINGULARITYENV_OMP_NUM_THREADS=1', r_path)
  }

  sgeset =init_sge(working_folder = wf, r_path = r_path, repeat_iterations = 2, output_files= oot, error_files = eet,
                   project_name = '<<<< PROJECT NAME REDACTED >>>>', other_options = NULL, slots_per_job = c(3,6),
                   package_location = '<<<< FILEPATH REDACTED >>>>')

  steak = init_stacker(stacking_models,
                       inlist = T,
                       data = df,
                       indicator = indicator,
                       indicator_family = indicator_family,
                       fe_equation = paste(names(all_cov_layers), sep = ' + '),
                       covariate_layers = all_cov_layers,
                       centre_scale = T,
                       time_var = 'year',
                       time_scale = year_list,
                       weight_col = 'weight',
                       num_fold_cols = nfc,
                       num_folds = nf,
                       cores = cores_to_use,
                       sge_parameters = sgeset)

  #run child models
  model_results = run_stacking_child_models(steak)

  #make rasters
  child_ras = make_all_children_rasters(st = steak, model_objects = model_results[[2]], time_points = year_list)


  #save checkpoint
  save(list = c('steak', 'model_results', 'child_ras'),
       file = '<<<< FILEPATH REDACTED >>>>')
} else{
  load('<<<< FILEPATH REDACTED >>>>')
}

if(start_point <= 3){
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  child_model_names <- names(steak$models)

  ## Start prepping for some INLA stuff
  #if running a normal model, or doing stacking, no gp
  if(permute == F | permutation %in% c('stacking', '1000_draws')){
    all_fixed_effects <- paste(child_model_names, collapse = ' + ')
  }

  ## copy things back over to df
  #note, stacking will remove data with NA values, hence why nrow(df) may not equal nrow(steak$data)
  df <- copy(merge(steak$data,model_results$preds, by = 'rid'))

  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
  df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]

  #if we want to constrain the coefs and keep the logited probs as covariates
  if(coefssum1 & (permute == F | permutation %in% c('stacking', '1000_draws'))){

    the_names = as.vector(outer(child_model_names,c('_full_pred','_cv_pred'), paste, sep=""))
    #fix the data frame
    df = df[, (the_names) := lapply(the_names, function(x) logit(get(x)))]


    #do the rasters now
    stacked_rasters = lapply(child_ras, logit)
  }else{
    stacked_rasters = child_ras
  }


  ## create a full raster list to carry though to the shiny/next steps
  cov_list         <- c(unlist(stacked_rasters),unlist(all_cov_layers))
  child_mod_ras    <- cov_list[child_model_names]

  toc(log = T) # End stacking master timer

  # make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

  ## Build spatial mesh over modeling area
  mesh_s <- build_space_mesh(d           = df,
                             simple      = simple_polygon,
                             max_edge    = mesh_s_max_edge,
                             mesh_offset = mesh_s_offset)

  ## Build temporal mesh (standard for now)
  mesh_t <- build_time_mesh(periods=eval(parse(text=mesh_t_knots)))

  #reconfigure save mbg input step
  #map the time step to period
  df = merge(df, period_map, by.x = 'year', by.y = 'data_period')
  #adjust cov_list to match the extent of simple raster and not polygon
  cov_list = lapply(cov_list, function(x) mask(setExtent(crop(x,simple_raster),simple_raster),simple_raster))

  ## Save all inputs for MBG model into correct location on /share
  to_save <- c("df", "simple_raster", "mesh_s", "mesh_t", 'cov_list', 'child_model_names', 'all_fixed_effects','period_map')
  save(list = to_save, file = '<<<< FILEPATH REDACTED >>>>')

}
## reload data an prepare for MBG
load('<<<< FILEPATH REDACTED >>>>')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if(start_point <= 4){
  tic("MBG - all") # Start MBG master timer

  ## for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
  df <- df[,paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]

  #if we ask for constrained coefficients, and its an appropriate time to do so
  if(coefssum1 & any(permutation %in% c('','stacking', '1000_draws'))){
    coefs.sum1 = T
  }else{
    coefs.sum1 = F
  }

  no_gp = F
  if(any(permutation %in% c('covs','stacking'))){
    no_gp = T
  }

  if(permutation == 'gponly'){
    usenull = T
  } else{
    usenull = F
  }


  #overrule random effect option if a single country model
  use_inla_country_res = as.logical(use_inla_country_res)
  if(length(get_gaul_codes(reg)) == 1){
    use_inla_country_res = F
  }

  if(coefs.sum1 | coefsmore0){
    exclude_these = strsplit(all_fixed_effects, " + ", fixed = T)[[1]]
  } else{
    exclude_these = ""
  }

  if(coefsmore0){
    pcv = all_fixed_effects
  }else{
    pcv = NULL
  }

  ## Generate MBG formula for INLA call
  mbg_formula <- build_mbg_formula_with_priors(fixed_effects   = all_fixed_effects,
                                               add_nugget      = as.logical(use_inla_nugget),
                                               add_ctry_res    = as.logical(use_inla_country_res),
                                               coefs.sum1      = coefs.sum1,
                                               no_gp           = no_gp,
                                               nullmodel       = usenull,
                                               positive_constrained_variables = pcv)

  ## Create SPDE INLA stack
  input_data <- build_mbg_data_stack(df            = df,
                                     fixed_effects = all_fixed_effects,
                                     mesh_s        = mesh_s,
                                     mesh_t        = mesh_t,
                                     use_ctry_res  = as.logical(use_inla_country_res),
                                     use_nugget    = as.logical(use_inla_nugget),
                                     coefs.sum1    = coefs.sum1,
                                     exclude_cs    = exclude_these)

  ## combine all the inputs
  stacked_input  <- input_data[[1]]
  spde           <- input_data[[2]]
  cs_df          <- input_data[[3]]

  ## Generate other inputs necessary
  outcome <- df[[indicator]] # N+_i - event obs in cluster
  N       <- df$N                  # N_i - total obs in cluster
  weights <- df$weight

  #catch in case there is no weight column
  if(is.null(weights)){
    weights = rep(1,nrow(df))
  }

  tic("MBG - fit model") # Start MBG - model fit timer

  ## Fit MBG model
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
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
                       intstrat         = intstrat)

  save(list = c('model_fit','cs_df', 'coefs.sum1', 'no_gp', 'mbg_formula'),
       file = '<<<< FILEPATH REDACTED >>>>')

} else{
  load('<<<< FILEPATH REDACTED >>>>')
}


toc(log = T) # End MBG - model fit timer

tic("MBG - predict model") # Start MBG - model predict timer

# Run predict_mbg on chunks of 50 samples (to avoid memory issues)
message('Making predictions in 50 draw chunks.')

max_chunk <- 50
samples   <- as.numeric(samples)


#reload functions until the old versions clear checkpoints
for(funk in list.files(recursive=TRUE,pattern='functions')){
  if(length(grep('central',funk))!=0){
    #message(funk)
    source(funk)
  }
}

# Create vector of chunk sizes
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
if(permutation == '1000_draws') samples = 1000
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
              coefs.sum1    = coefs.sum1,
              pred_gp       = !no_gp)[[3]] #we've already logited the appropriate covs
})

## Make cell preds and a mean raster
cell_pred <- do.call(cbind, pm)
mean_ras  <- insertRaster(simple_raster,matrix(rowMeans(cell_pred),ncol = max(period_map$period)))

toc(log = T) # Stop MBG - model predict timer


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Finish up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

message('Wrapping up')
message(run_date)
summary(model_fit)
save_mbg_preds(config     = config,
               time_stamp = time_stamp,
               run_date   = run_date,
               mean_ras   = mean_ras,
               sd_ras     = NULL,
               res_fit    = model_fit,
               cell_pred  = cell_pred,
               df         = df,
               pathaddin  = pathaddin)

system(paste("kill", meminfo_pid))

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
