## join_modeling_regions() ---------------------------------------------------->
#' 
#' @title Join MBG Modeling Region
#' 
#' @description Given the input MBG dataset and a list of MBG modeling regions,
#'   join each NID to the modeling region that it is expected to fully overlap.
#'   This function relies on get_adm0_codes to identify overlapping country
#'   codes between regions and the individual ISO3 codes associated with each
#'   NID.
#'   
#'   ** ASSUMPTIONS (function will break if these are not met) ** 
#'   1. Each NID corresponds to a single modeling region
#'   2. The various modeling regions are non-overlapping by country (ie. there
#'      are no countries fully within multiple modeling regions).
#'   
#'   If a country does not fall within any known modeling region, either because
#'   it is not in the location lookup table or because it is not in the list of
#'   modeled countries, the assigned modeling region will be 'None'.
#'   
#' @param in_data Full dataset read from MBG input_data folder
#' @param regions Vector of all modeling regions
#' 
#' @return A data.table with two fields, 'nid' and 'model_region'
#' 
join_modeling_regions <- function(in_data, regions){
  # Get admin0 lookup table that will be used for all calls to get_adm0_codes
  lookup_table <- load_adm0_lookup_table()

  # Get admin0 codes associated with each modeling region
  reg_dt <- suppressMessages(lapply(
    as.list(regions), 
    function(rr) {
      data.table(
        adm0_code = get_adm0_codes(rr, lookup_table=lookup_table),
        model_region = rr
      )
    }
  )) %>% rbindlist

  # Check that modeling regions are non-overlapping
  check_reg_dt <- copy(reg_dt)
  check_reg_dt[, ct := 1 ]
  reg_dt_errors <- check_reg_dt[, 
    .(ct=sum(ct), model_regions=sort(unique(model_region))),
    by=.(adm0_code)
  ][ct > 1,][order(adm0_code)]
  if(nrow(reg_dt_errors) > 0){
    message("There are countries listed within multiple modeling regions.")
    print(reg_dt_errors)
    stop("Please resolve duplicate country issues before continuing.")
  }

  # Get all unique NIDs and ISO codes
  nid_iso <- unique(in_data[, .(nid, country)])
  nid_iso[, country := tolower(country) ]
  setnames(nid_iso, 'country', 'iso3')

  # Merge on adm0 codes, then use adm0 codes to merge onto modeling regions
  adm0_type <- detect_adm_shapefile_date_type(
    get_admin_shapefile(version = 'current', raking = subnational_raking)
  )$shpfile_type
  adm0_field <- ifelse(adm0_type=='gadm', 'gadm_geoid', 'GAUL_CODE')
  lookup_for_merge <- copy(lookup_table)
  setnames(lookup_for_merge, adm0_field, 'adm0_code')
  nid_to_admin <- merge(
    x = nid_iso,
    y = lookup_for_merge[, .(iso3, adm0_code)],
    by = c('iso3'),
    all.x = TRUE
  )
  # Merge onto modeling regions
  nid_to_reg <- merge(
    x = nid_to_admin,
    y = reg_dt,
    by = c('adm0_code'),
    all.x = TRUE
  )

  # Check the number of unique modeling regions by NID
  nid_to_reg <- unique(nid_to_reg[, .(model_region, nid)])
  nid_to_reg[, ct := 0 ]
  nid_to_reg[!is.na(model_region), ct := 1 ]
  reg_by_nid <- nid_to_reg[, 
    .(ct=sum(ct), model_region=paste(sort(unique(model_region)),collapse=';;')),
    by=.(nid)
  ]
  reg_by_nid[ ct==0, model_region := 'None' ]
  # Throw an error if a single NID matches to more than one modeling region
  reg_errors <- reg_by_nid[ ct > 1, ]
  if(nrow(reg_errors) > 0){
    message("Some NIDs in your dataset match with multiple modeling regions.")
    print(reg_errors)
    stop("Please resolve region uniqueness issues before continuing.")
  }

  # Return the dataset
  return(reg_by_nid[, .(nid, model_region)])
}


## make_holdouts_by_region() -------------------------------------------------->
#' 
#' @title Make holdouts by region
#' 
#' @description Given a table of all input data for MBG, create holdouts so
#'   that held-out NIDs are each sampled evenly across regions.
#' 
#' @param in_data Full dataset read from MBG input_data folder
#' @param regions Vector of all modeling regions
#' @param num_holdouts Number of holdout runs
#' 
#' @return data.table with three following fields: 'nid', 'model_region', and
#'   'fold'
#' 
make_holdouts_by_region <- function(in_data, regions, num_holdouts=5){
  # Associate model regions with each unique NID
  nids_by_region <- join_modeling_regions(
    in_data = in_data, 
    regions = regions
  )
  # Get list of all unique regions, including the "None" region signifying
  #  a source that does not fall squarely within any modeling region
  data_regions <- sort(unique(nids_by_region[, model_region]))
  # Helper function to evenly sample across NIDs in a given region
  sample_region <- function(sub_data){
    # Create a vector that has the same length as the number of NIDs and is the 
    #  vector 1:n_holdouts repeating. For an example with 12 NIDs and 5 holdouts,
    #  this would produce the vector c(1,2,3,4,5,1,2,3,4,5,1,2)
    num_nids <- nrow(sub_data)
    holdouts <- rep_len(1:num_holdouts, length.out=num_nids)
    # Shuffle the vector using `sample()` and return as the holdout number
    sub_data$fold <- sample(x=holdouts, size=length(holdouts), replace=FALSE)
    # Return the dataset with a new field, `fold`
    return(sub_data)
  }
  # Sample the data across each region
  sampled_by_region <- lapply(
    as.list(data_regions),
    function (rr) sample_region(nids_by_region[model_region==rr,])
  ) %>% rbindlist
  # Return the holdout table, now with the additional `fold` field
  return(sampled_by_region)
}



## make_qsub_split_jobs ------------------------------------------------------->
#' 
#' @title Make Qsub with split fit and predict jobs
#' 
#' @description This function is essentially a clone of `make_qsub_share()` that
#'   is designed to submit separate TMB fitting ('fitting') and draw-generation 
#'   ('predict') programs. This is a short-term solution that will eventually
#'   be fixed in an updated MBG pipeline.
#' 
#' For information about parameters, see `make_qsub_share()` in the lbd_core
#'   repository in "<<<< FILEPATH REDACTED >>>>"
#' 
make_qsub_split_jobs <- function(
  user             = Sys.info()['user'],
  code             = NULL,
  cores_fit        = slots,
  cores_predict    = slots,
  memory_fit       = 100,
  memory_predict   = 100,
  proj             = NULL,
  ig               = indicator_group,
  indic            = indicator,
  reg              = "test",
  age              = 0,
  rd               = run_date,
  log_location     = 'sharedir',
  addl_job_name    = '',
  saveimage        = FALSE,
  test             = FALSE,
  holdout          = 0,
  corerepo         = core_repo,
  indic_repo       = NULL, # Indicator-specific repo
  geo_nodes        = FALSE,
  use_c2_nodes     = FALSE,
  singularity      = 'default',
  singularity_opts = NULL,
  predict_only     = FALSE # When you ONLY want to submit predict jobs
){
  # The indicator-specific repository MUST be specified
  if(is.null(indic_repo)) stop("The indicator-specific repository must be specified!")

  # save an image
  if(saveimage==TRUE) save.image(pre_run_image_path(ig, indic, rd, age, reg, holdout))

  # Define project first (necessary to validate node options)
  proj <- get_project(proj, use_geo_nodes=geo_nodes)

  # Validate arguments
  validate_singularity_options(singularity, singularity_opts)
  validate_node_option(geo_nodes, use_c2_nodes, proj)

  # Create sharedir
  sharedir = get_model_output_dir(ig, indic, rd)
  dir.create(sharedir, showWarnings = FALSE)

  # Determine where stdout and stderr files will go
  output_err = setup_log_location(log_location, user, indic, ig, rd)
  output_log_dir = output_err[[1]]
  error_log_dir = output_err[[2]]

  # Define remaining job attributes
  job_name <- paste0("job_", addl_job_name, "_", reg, "_", age, "_", holdout)
  queue <- get_queue(use_geo_nodes=geo_nodes, use_c2_nodes=use_c2_nodes)
  shell <- paste0("<<<< FILEPATH REDACTED >>>>")
  sing_image <- get_singularity(image = singularity)

  # resources are all the -l qsub arguments
  resources_fit <- get_resources(
    use_geo_nodes=geo_nodes, 
    cores=cores_fit, 
    ram_gb=memory_fit,
    runtime='3:00:00:00'
  )
  resources_predict <- get_resources(
    use_geo_nodes=geo_nodes, 
    cores=cores_predict, 
    ram_gb=memory_predict,
    runtime='3:00:00:00'
  )

  code_fit <- sprintf("<<<< FILEPATH REDACTED >>>>")
  code_predict <- sprintf("<<<< FILEPATH REDACTED >>>>")

  if(predict_only==FALSE){
    # RUN BOTH FIT AND PREDICT
    qsub_fit <- generate_qsub_command(
      # qsub-specific arguments
      stderr_log=error_log_dir,
      stdout_log=output_log_dir,
      project=proj,
      resources=resources_fit,
      job_name=paste0(job_name,'_fit'),
      singularity_str=qsub_sing_envs("", singularity_opts, sing_image),
      cores=cores_fit,
      queue=queue,
      # Command to qsub
      shell, code_fit, reg, age, rd, as.numeric(test), holdout, indic, ig, "fin")

    # Submit job
    message(qsub_fit)
    qsub_output <- system(qsub_fit, intern=TRUE)
    message(qsub_output)

    # Parse output from the first job submission to get the hold_jid for the
    #  predict job
    qsub_output <- qsub_output[ grepl('has been submitted', qsub_output) ]
    if (length(qsub_output) == 0) stop("Issue with submission of the fitting job.")

    fit_job_id <- gsub('Your job ', '', qsub_output[1]) %>% gsub(" (.)+$", "", .)
    if(!(grepl('^([[:digit:]])+$', fit_job_id))){
      stop(sprintf("ID of previously submitted job, %s, is invalid.", fit_job_id))
    }

    qsub_predict <- generate_qsub_command(
      # qsub-specific arguments
      stderr_log=error_log_dir,
      stdout_log=output_log_dir,
      project=proj,
      resources=resources_predict,
      job_name=paste0(job_name,'_predict'),
      singularity_str=qsub_sing_envs("", singularity_opts, sing_image),
      cores=cores_predict,
      queue=queue,
      # Hold Job ID (for fitting job)
      '-hold_jid', fit_job_id,
      # Command to qsub
      shell, code_predict, reg, age, rd, as.numeric(test), holdout, indic, ig, "fin")
    message(qsub_predict)
    system(qsub_predict)    
  } else {
    # RUN PREDICT ONLY
    qsub_predict <- generate_qsub_command(
      # qsub-specific arguments
      stderr_log=error_log_dir,
      stdout_log=output_log_dir,
      project=proj,
      resources=resources_predict,
      job_name=paste0(job_name,'_predict'),
      singularity_str=qsub_sing_envs("", singularity_opts, sing_image),
      cores=cores_predict,
      queue=queue,
      # Command to qsub
      shell, code_predict, reg, age, rd, as.numeric(test), holdout, indic, ig, "fin")
    message(qsub_predict)
    system(qsub_predict)    
  }

  return(NULL)
}



#################################################################################
### Combine cell pred files for age bins


combine_agebins <- function(age_bins    = 1:5,
                            reg,
                            ageasindic  = TRUE,    # is each age bin its own indicator? 
                           monthlyprob = FALSE,   # modeled monthly or bin probability?
                            fractions   = list(neonatal = 1,     # which age bins make up the bigger ones
                                               infant   = 1:3,
                                               under5   = 1:5)){
  
  if(monthlyprob == TRUE)
    stop('Model Conditional Probabilities, modelling monthly probabilities is depreciated.')
  
  message('\nLoading Survival Data for:')
  for (bin in age_bins) {
    message(sprintf('age bin %i', bin))
    ## object name
    name <- sprintf('surv_all_bin%i',bin)
    ## load the object (named cell_preds)
    assign(
      name,
      1 - load_cell_preds(
        indicator_group = indicator_group,
        indicator       = 'died',
        rd              = run_date,
        region          = reg,
        agebin          = bin,
        u5m             = TRUE,
        other           = '',
        ageasindic      = ageasindic
      )
    )
  }
  gc(full=TRUE)
  
  message('\nConstructing mortality probability samples for larger bins of interest:')
  # Create empty list to populate
  res <- lapply(as.list(names(fractions)), function(x) NULL)
  names(res) <- paste0(names(fractions),'_all')
  for (f in 1:length(fractions)) {
    message(names(fractions)[f])
    ## get bins required
    bins <- fractions[[f]]
    n_bins <- length(bins)    
    ## names of all survival probability objects in this fraction
    names_from <- sprintf('surv_all_bin%i',bins)
    ## output mortality fraction object name
    name_to <- sprintf('%s_all', names(fractions)[f])
    ## multiply mortalities for all objects in names_from    
    ## Construct survival probability statement to evaluate
    res[[name_to]] <- get(names_from[1])
    ## if more than 1, multiply by all subsequent ones
    if (n_bins > 1) {
      for (j in 2:n_bins) {
        res[[name_to]] <- res[[name_to]] * get(names_from[j])
        gc(full=TRUE)
      }
    }
    # Mortality = 1 - prod(survival by age group)
    res[[name_to]] <- 1 - res[[name_to]]
  }
    
  return(res)
}





#################################################################################
### FULL POST EST FOR U5M

u5mpostest <- function(reg,
                       subnational_condsim = TRUE,
                       ageasindic          = TRUE,
                       yl                  = year_list,
                       agegroups           = c('neonatal','infant','under5')){
  
  
  message(reg)
  
  # initiate output results list
  outputlist <- list()
  
  # Combine cell_pred age groups in each region
  cell_pred <- combine_agebins(reg=reg)
  gc(full=TRUE)
  
  ## get aggregated estimates for all admin0. Aggregate to level you rake to
  message('Loading spatial data')
  simple_polygon <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4)
  subset_shape   <- simple_polygon[['subset_shape']]
  simple_polygon <- simple_polygon[['spoly_spdf']]
  
  raster_list    <- build_simple_raster_pop(subset_shape) #,u5m=TRUE)
  simple_raster  <- raster_list[['simple_raster']]
  pop_raster     <- raster_list[['pop_raster']]
  
  ## Pull 2000-2015 annual population brick using new covariates function
  if (class(yl) == "character") yl <- eval(parse(text=yl))
  pop_raster_annual <- load_and_crop_covariates_annual(covs           = 'worldpop',                
                                                       measures       = pop_measure,# from config      
                                                       simple_polygon = simple_polygon,
                                                       start_year     = min(yl),
                                                       end_year       = max(yl),
                                                       interval_mo    = as.numeric(interval_mo),
                                                       agebin=1)[[1]]
  
  pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
  pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
  pop_raster_annual  <- mask(pop_raster_annual, simple_raster)
  
  ## Create population weights using the annual brick and feed custom year argument to aggregation function
  pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                          simple_raster = simple_raster,
                                          pop_raster    = pop_raster_annual,
                                          gaul_list     = get_gaul_codes(reg))
  
  
  # for both child and neonatal
  for(group in agegroups){
    message(paste0(group,group,group,group,group))
    
    outdir <- sprintf("<<<< FILEPATH REDACTED >>>>")
    
    # draw level cond sim at national level
    message('  National Cond Sim')
    cond_sim_draw_adm0 <- make_condSim(admin_level    = 0,
                                       pop_wts_object = pop_wts_adm0,
                                       cell_pred      = cell_pred[[sprintf('%s_all',group)]],
                                       gaul_list      = get_gaul_codes(reg),
                                       summarize      = FALSE,
                                       years          = yl)
    
    # save some raw country estimates to compare with GBD in a plot later on
    cond_sim_raw_adm0 <- apply(cond_sim_draw_adm0   , 1, mean)
    adm0_geo          <- cbind(mean=cond_sim_raw_adm0,
                               lower=apply(cond_sim_draw_adm0, 1, quantile, probs=.025),
                               upper=apply(cond_sim_draw_adm0, 1, quantile, probs=.975))
    outputlist[[sprintf('%s_%s_adm0_geo',reg,group)]] <- data.table(split_geo_names(adm0_geo),adm0_geo)
    
    # get gbd estimates for raking
    gbd <- load_u5m_gbd(gaul_list = get_gaul_codes(reg), 
                        years     = yl,
                        age_group = group, 
                        exact     = TRUE)
    
    ## Get raking factors
    message('  Establishing raking factors')
    rf   <- calc_raking_factors(agg_geo_est = cond_sim_raw_adm0,
                                rake_to     = gbd)
    
    
    # rake the cell preds
    raked_cell_pred <- rake_predictions(raking_factors = rf,
                                        pop_wts_object = pop_wts_adm0,
                                        cell_pred      = cell_pred[[sprintf('%s_all',group)]])
    
    ## summarize raked predictions for each cell
    message('  Summarizing unraked cell preds:\nmean')
    outputlist[[sprintf('%s_mean_unraked_%s_raster',reg,group)]] <-
      make_cell_pred_summary( draw_level_cell_pred = cell_pred[[sprintf('%s_all',group)]],
                              mask                 = simple_raster,
                              return_as_raster     = TRUE,
                              summary_stat         = 'mean')
    
    message('  Summarizing raked cell preds:')
    for(summeasure in c('mean','cirange','lower','upper')){
      message(sprintf('    %s',summeasure))
      outputlist[[sprintf('%s_%s_raked_%s_raster',reg,summeasure,group)]] <-
        make_cell_pred_summary( draw_level_cell_pred = raked_cell_pred,
                                mask                 = simple_raster,
                                return_as_raster     = TRUE,
                                summary_stat         = summeasure)
    }
    
    
    # do subnational 
    if(subnational_condsim){
      message('  Subnational aggregation:')
      
      for(ad in 1:2){
        message(sprintf('    Admin Level %i',ad))
        
        pop_wts <- make_population_weights(admin_level   = ad,
                                           simple_raster = simple_raster,
                                           pop_raster    = pop_raster_annual,
                                           gaul_list     = get_gaul_codes(reg))
        message('     raked')
        condsim <-    make_condSim(admin_level   = ad,
                                   pop_wts_object = pop_wts,
                                   cell_pred      = raked_cell_pred,
                                   gaul_list      = get_gaul_codes(reg),
                                   summarize      = FALSE,
                                   years          = yl)
        
        message('     unraked')
        condsim2 <-    make_condSim(admin_level   = ad,
                                    pop_wts_object = pop_wts,
                                    cell_pred      = cell_pred[[sprintf('%s_all',group)]],
                                    gaul_list      = get_gaul_codes(reg),
                                    summarize      = FALSE,
                                    years          = yl)
        
        outputlist[[paste0('cond_sim_raked_adm',ad,'_',group,'_',reg)]] <- 
          data.table(cbind(split_geo_names(as.matrix(condsim)),mean=unname(condsim) ))
        outputlist[[paste0('cond_sim_unraked_adm',ad,'_',group,'_',reg)]] <- 
          data.table(cbind(split_geo_names(as.matrix(condsim2)),mean=unname(condsim2) ))
        
        rm(condsim);rm(condsim2)
        
      }
    }
    
    # save cell preds and rf file
    message('  saving cell preds')
    save(raked_cell_pred, file=sprintf("<<<< FILEPATH REDACTED >>>>")
    cell_draws <- cell_pred[[sprintf('%s_all',group)]]
    save(cell_draws, file=sprintf("<<<< FILEPATH REDACTED >>>>")
    rm(cell_draws)
    
    write.csv(rf, file=sprintf("<<<< FILEPATH REDACTED >>>>")
    save_post_est(rf[,.(ADM0_CODE = name, year = year, mean = rake_to_mean)], 
                  'csv', 'mean_raked_adm0',   indic = sprintf('died_%s',group))
    
  }
  
  # return output list
  return(outputlist)
  
}


load_u5m_gbd <-function(gaul_list,
                        age_group='under5',
                        years = c(2000,2005,2010,2015),
                        exact = TRUE, # if folse, takes a 5 year average
                        getci = FALSE){
  
  require(plyr)
  
  
  if(! age_group %in% c('neonatal','infant','under5'))
    stop('Group must be either "neonatal", "infant", or "under5". ')
  
  # load in gbd estimates
  gbd <- fread("<<<< FILEPATH REDACTED >>>>")
  message('WARNING: PULLING GBD FILE FROM J DRIVE, CONFIRM WITH MORTALITY TEAM THIS IS UP DATE. ')
  
  
  # rename some variables. 
  setnames(gbd, c('year_id','sex_id','age_group_id','ihme_loc_id'), c('year', 'sex', 'age', 'iso3'))
  
  ## 1 is under-5, 2 is early neonatal, 3 is late neonatal, 4 is post-neonatal, 5 is ages 1-4, 28 is under-1, and 42 is neonatal
  ## both sex only
  ## study years only
  if(exact==TRUE){
    # subset to exact years and u5 ages needed
    gbd <- subset(gbd,sex == 3 &
                    age  %in% c(1,28,42) &
                    year %in% years)
  } else {
    # if not exact, get a mean over years
    yrmap <- data.table(year=1998:2017,yrmap=rep(years,each=5))
    gbd   <- merge(gbd,yrmap,by='year')
    gbd   <- subset(gbd,sex == 3 &
                      age  %in% c(1,28,42))
    
    # CIs are just an approx this way, but not using anyway for analysis just plots so ok
    gbd <- gbd[,.(qx_mean=mean(qx_mean),
                  qx_lower=mean(qx_lower),
                  qx_upper=mean(qx_upper)),
               by=.(sex,age,iso3,location_id,yrmap)]
    gbd <- rename(gbd, c('yrmap'='year'))
    
  }
  
  # both sexes for now, so sex is not an important variable here
  gbd$sex=NULL
  
  # load in country id info an merge it
  gaul_to_loc_id <- fread("<<<< FILEPATH REDACTED >>>>")
  gaul_to_loc_id <- gaul_to_loc_id[GAUL_CODE%in%gaul_list,]
  gbd <- merge(gbd,gaul_to_loc_id,by.x='location_id',by.y='loc_id',all.y=T)
  
  # keep only the year for this age group of interest
  if(age_group=='neonatal')  gbd <- gbd[age==42,]
  if(age_group=='infant')    gbd <- gbd[age==28,]
  if(age_group=='under5')    gbd <- gbd[age==1,]
  
  # rename some variables
  setnames(gbd, c('qx_mean','qx_lower','qx_upper','GAUL_CODE'), 
           c('mean',   'lower',   'upper',   'name'))
  
  # keep upper lower, or dont.
  if(getci)  gbd = gbd[,c('name','year','mean','lower','upper'),with=FALSE]
  if(!getci) gbd = gbd[,c('name','year','mean'),with=FALSE]
  
  # return the gbd dataset
  return(gbd)
}


make_alpha_raster <-  function(ig,
                               ind,
                               rd,
                               alpha = 0.1,
                               reg   = NULL,
                               year_list,
                               scope,  # country or regional
                               simple_raster = NULL, # can be supplied to skip the step of pulling it if use has it in env already
                               raked      = TRUE
){
  
  # load the data - must be in the proper lucas format: "<<<< FILEPATH REDACTED >>>>"
  rasterbrick <- sprintf("<<<< FILEPATH REDACTED >>>>")
  r <- brick(rasterbrick)
  
  # keep only first and last year
  r <- brick(r[[1]],r[[dim(r)[3]]])
  
  # check for a weird args combination
  if(length(get_gaul_codes(reg))==1 & scope == 'regional'){
    scope <- 'country'
    message('Gaul length is one, so defaulting back to country from regional scope.')
  }
  
  # 
  if(is.null(reg)){
    message('Region was null so defaulting to the results extent from this run date')
    reg <- 'all'
  }
  
  # get admin0 raster
  if(is.null(simple_raster)){
    simple_polygon <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4)
    simple_raster <- build_simple_raster_pop(simple_polygon[['subset_shape']])[['simple_raster']]
  }
  
  # get population rasters
  p <- load_and_crop_covariates_annual(covs           = 'worldpop',                
                                       measures       = pop_measure, # from config      
                                       simple_polygon = r[[1]],
                                       start_year     = min(year_list),
                                       end_year       = max(year_list),
                                       interval_mo    = as.numeric(interval_mo),
                                       agebin=1)[[1]]
  p  <- brick(p[[1]],p[[dim(p)[3]]])
  
  
  a  <- crop(simple_raster,r) # gauls
  
  # check rasters are all good
  if(!all(dim(r)[1:2] == dim(a)[1:2] & dim(a)[1:2] == dim(p)[1:2]))
    stop('Admin, Results, Population Rasters do not seem to conform')
  
  # helper function to get the locations
  getalphalocs<-function(gaul,       # country gaul code(s)
                         alpha,      # alpha level
                         mr          # master raster (used if iteratively adding countries to the regional map)
  ){
    # tmp rasters
    ta <- a
    tp <- p 
    tr <- r
    
    # subset to desired gauls
    ta[!ta %in% gaul] <- NA
    tp <- raster::mask(tp,ta)
    tr <- raster::mask(tr,ta)
    tr1 <- tr[[1]]
    tr2 <- tr[[2]]
    
    # identify worst of alpha% of the pop
    require(Hmisc)
    #require(reldist)
    cutoff1 <- Hmisc::wtd.quantile(as.vector(tr[[1]]),weights=as.vector(tp[[1]]),probs=1-alpha)
    cutoff2 <- Hmisc::wtd.quantile(as.vector(tr[[2]]),weights=as.vector(tp[[2]]),probs=1-alpha)
    
    # categorize them and add to the master raster
    tr[[1]][tr[[1]] <  cutoff1]  <- 0
    tr[[2]][tr[[2]] <  cutoff2]  <- 0
    tr[[1]][tr[[1]] >= cutoff1]  <- 1
    tr[[2]][tr[[2]] >= cutoff2]  <- 2
    tr[[1]][is.na(tr[[1]])]      <- 0
    tr[[2]][is.na(tr[[2]])]      <- 0
    mr <- mr + tr[[1]] + tr[[2]]
    
    return(mr)
  }
  
  all_gauls <- unique(a)
  master_raster <- a
  values(master_raster) <- 0
  
  # if scope is country, then loop over every country in this region and add them to the master shapefile
  if(scope == 'country'){
    for(gaul in all_gauls){
      message(gaul)
      master_raster <- getalphalocs(gaul = gaul, mr = master_raster, alpha = alpha)
    }
  }
  
  # if scope is region, then do one run of the alphalocs function with all gaul codes
  if(scope == 'region'){
    master_raster <- getalphalocs(gaul = all_gauls, mr = master_raster, alpha = alpha)
  }
  
  return(master_raster)
  
}


aroc_need  <-  function(ig,
                        ind,
                        rd,
                        goal,
                        goalyr,
                        fromyr,
                        yl,
                        raked      = TRUE
){
  
  # load the data - must be in the proper lucas format: "<<<< FILEPATH REDACTED >>>>"
  rasterbrick <- sprintf("<<<< FILEPATH REDACTED >>>>")
  
  # keep only from year
  r <- brick(rasterbrick)[[which(year_list==fromyr)]]
  
  # figure the goal aroc needed
  rd  <- na.omit(values(r))
  gl  <- rep(goal,length(rd))
  res <- log(gl/rd)/(goalyr-fromyr)
    
  # inserraster
  res <- seegMBG::insertRaster(r,cbind(res))
  
  return(res)
  
}


lucas_plot <- function(ig,
                       ind,
                       rd,
                       regs,
                       year_list,
                       years_to_plot = c(2000,2005,2010,2015),
                       adm1_borders = FALSE,
                       ind_territories=FALSE,
                       raked = TRUE){
  
  if(length(years_to_plot) > 4)
    stop('years_to_plot must be 4 or less in length')

  # load the data - must be in the proper lucas format: "<<<< FILEPATH REDACTED >>>>"
  rasterbrick <- sprintf("<<<< FILEPATH REDACTED >>>>")
  r <- brick(rasterbrick)
  idx <- which(year_list %in% years_to_plot)
  r <- r[[idx]]
  
  # gauls
  gauls <- c()
  for(re in regs)
    gauls <- c(gauls, get_gaul_codes(re))
  
  # load in adm0 and adm1
  gaulpath <- "<<<< FILEPATH REDACTED >>>>"
  sad0 <- shapefile("<<<< FILEPATH REDACTED >>>>")
  if(ind_territories==TRUE){
    sad0[[paste0('ADM', 0,'_CODE')]][sad0[[paste0('ADM', 0,'_CODE')]]==52   ]=115
    sad0[[paste0('ADM', 0,'_CODE')]][sad0[[paste0('ADM', 0,'_CODE')]]==40781]=115
    sad0[[paste0('ADM', 0,'_CODE')]][sad0[[paste0('ADM', 0,'_CODE')]]==15   ]=115  
  }
  
  sad0 <- sad0[sad0@data$ADM0_CODE %in% gauls,]
  sad0 <- gSimplify(sad0, tol = 0.02, topologyPreserve = TRUE)

  if(adm1_borders==TRUE){
    sad1 <- shapefile("<<<< FILEPATH REDACTED >>>>")
    sad1 <- sad1[sad1@data$ADM0_CODE %in% gauls,]
    sad1 <- gSimplify(sad1, tol = 0.02, topologyPreserve = TRUE)
  }

  # scale values to 1000
  rr         <- r
  values(rr) <- as.vector(rr)*1000
  rr[rr>200] <- 200
  
  # define color breaks
  col.f1 <- colorRampPalette(c("#e58bba", "#f2e8b5"))
  col.f2 <- colorRampPalette(c("#f2e8b5", "#ed152e"))
  if(ind == 'died_under5'){
    breaks <- c(0, 25,
                26:50,
                51:200,
                max(c(200,max(as.vector(rr)))))
    col   <- c("#74039E",
               col.f1(25),
               col.f2(150),
               "#ED152E")
    arg <- list(at=c(25,50,200), labels=c('<=25','50','200+'))
    
  } else {
    breaks <- c(0, 12,
                13:50,
                51:200,
                max(c(200,max(as.vector(rr)))))
    col   <- c("#74039E",
               col.f1(38),
               col.f2(150),
               "#ED152E")
    arg <- list(at=c(12,50,200), labels=c('<=12','50','200+'))
    
  }
  

  # save the plot
  pdf(sprintf("<<<< FILEPATH REDACTED >>>>"), 
      height = ifelse(length(years_to_plot) >= 3, 8, 5), 
      width  = ifelse(length(years_to_plot) >= 2, 8, 5))
  if(length(years_to_plot) >= 3) 
    par(mfrow=c(2,2))
  if(length(years_to_plot) == 2 )
    par(mfrow=c(1,2))
  
  for(y in 1:length(years_to_plot)){
    raster::plot(rr[[y]], 
                 axes      = FALSE,
                 breaks    = breaks,
                 col       = col,
                 maxpixels = length(rr[[y]]),
                 legend    = FALSE,
                 main      = years_to_plot[y])
    box(col='white')
    lines(sad0, lwd = 0.8)
    if(adm1_borders==TRUE) 
      lines(sad1, lwd = 0.1)
  }
  plot(rr[[1]], legend.only=TRUE, col=col, breaks= breaks, legend.width=1, legend.shrink=0.75,
       smallplot = c(0.03,0.07, 0.3, 0.5),  axis.args = arg)
  par(mar = par("mar"))
  dev.off()

}


#' Title
#'
#' @param file_path Full file path to csv containing GBD estimates for U5M for India Subnationals 
#' @param age_group "Infant", "neonatal", or "under5".
#' @param years List of years.
#' @param getci Include upper and lower values in dataframe
#'
#' @return Dataframe with columns "name", "year", and "mean", corresponding to the gbd_loc_id, year, and mean value for that location-year
#' @export
#'
#' @examples rake_to <- load_u5m_gbd_india_state(age_group = "neonatal")
#' @note for u5m GBD india subnational raking for lalit only
load_u5m_gbd_india_state <- function(file_path = "<<<< FILEPATH REDACTED >>>>",
                                     age_group ='under5',
                                     years = c(2000:2015),
                                     getci = FALSE){
  
  if(! age_group %in% c('neonatal','infant','under5'))
    stop('Group must be either "neonatal", "infant", or "under5". ')
  
  #read in gbd variables
  gbd <- fread(file_path)
  
  # rename some variables. 
  setnames(gbd, c('year_id','sex_id','age_group_id', 'location_id'), c('year', 'sex', 'age', "name"))
  
  
  # subset to exact years and u5 ages needed
  gbd <- subset(gbd,sex == 3 &
                  age  %in% c(1,28,42) &
                  year %in% years)
  
  # keep only the year for this age group of interest
  if(age_group=='neonatal')  gbd <- gbd[age==42,]
  if(age_group=='infant')    gbd <- gbd[age==28,]
  if(age_group=='under5')    gbd <- gbd[age==1,]
  
  # keep upper lower, or dont.
  if(getci)  gbd = gbd[,c('name','year','mean','lower','upper'),with=FALSE]
  if(!getci) gbd = gbd[,c('name','year','mean'),with=FALSE]
  
  # return the gbd dataset
  return(gbd)
}

# Defining a function that will get the raster versions of each Admin level:
GetAdmin<-function(admin_level,simple_raster, region_adm0_list, shapefile_version){
  message(paste0("Loading admin level ",admin_level))
  
  # load admin shape file
  admin_shp <- rgdal::readOGR(dsn=get_admin_shapefile(admin_level, version = shapefile_version))
  
  # ensure that the rasterize variable is a numeric
  admin_shp@data[[paste0('ADM', admin_level, '_CODE')]] <- as.numeric(as.character(admin_shp@data[[paste0('ADM', admin_level, '_CODE')]]))
  
  # if it doesn't exist, get areas of polygons. 
  if(is.null(admin_shp$Shape_Area)){
    admin_shp$Shape_Area <- area(admin_shp) / 1e6
  }
  
  message("Rasterizing with the custom function...")
  # we order by area so small places don't get buried under big places (e.g. Lesotho and S. Africa)
  admin_rast<-rasterize_check_coverage(admin_shp[order(admin_shp$Shape_Area),],simple_raster,paste0("ADM",admin_level,"_CODE"), fun="first")
  
  message("Converted to raster based on simple_raster template. Cropping and masking:")
  admin_rast  <- crop(admin_rast, extent(simple_raster))
  admin_rast  <- setExtent(admin_rast, simple_raster)
  admin_rast  <- mask(admin_rast, simple_raster)
  
  message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
  admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% region_adm0_list,]
  admin_centroids<-SpatialPointsDataFrame(gCentroid(admin_shp, byid=TRUE), admin_shp@data, match.ID=FALSE)
  
  message("Compiling and returning results.")
  admin<-list()
  admin[["spdf"]]<-admin_shp
  admin[["centroids"]]<-admin_centroids
  admin[["rast"]]<-admin_rast
  admin[["attributes"]]<-copy(data.table(admin_shp@data))
  
  return(admin)
}


#' get the raster versions of a given admin level
#'
#' Returns a list of objects used by fing_missing_adms()
#' 
#' @param admin_level 0, 1, 2
#' @param simple_raster simple raster
#' @param region_adm0_list region name
#'
GetAdmin_fast<-function(admin_level,simple_raster, region_adm0_list=NULL){
  # This function should be run in the LBD singularity image
  library(sp)
  library(rgeos)
  library(sf)
  library(fasterize)
  message(paste0("Loading admin level ",admin_level))
  
  # load admin shape file
  admin_shp <- sf::st_read(get_admin_shapefile(admin_level = admin_level, version = modeling_shapefile_version))
  
  # ensure that the rasterize variable is a numeric
  adm_field <- paste0('ADM', admin_level, '_CODE')
  admin_shp[[adm_field]] <- as.numeric(as.character(admin_shp[[adm_field]]))
  
  message("Rasterizing...")
  admin_rast <- fasterize::fasterize(
    sf     = admin_shp,
    raster = raster(simple_raster), # Ensures that output is a rasterLayer
    field  = adm_field,
    fun    = 'first'
  )
  
  message("Converted to raster based on simple_raster template. Cropping and masking:")
  admin_rast  <- crop(admin_rast, extent(simple_raster))
  admin_rast  <- setExtent(admin_rast, simple_raster)
  admin_rast  <- mask(admin_rast, mask=simple_raster)

  # Convert admin_rast to rasterBrick (expected output format)
  if("RasterLayer" %in% class(admin_rast)) admin_rast <- brick(admin_rast)
  
  # Convert df object back to SpatialPolygonsDataFrame
  admin_shp <- readOGR(get_admin_shapefile(admin_level = admin_level, version = modeling_shapefile_version))
  message("Subsetting polygon and point objects to only contain the relevant ADM0 codes; calculating centroids.")
  if(!is.null(region_adm0_list)) admin_shp<-admin_shp[admin_shp@data$ADM0_CODE %in% region_adm0_list,]
  admin_centroids<-SpatialPointsDataFrame(gCentroid(admin_shp, byid=TRUE), admin_shp@data, match.ID=FALSE)
  
  message("Compiling and returning results.")
  admin<-list()
  admin[["spdf"]]<-admin_shp
  admin[["centroids"]]<-admin_centroids
  admin[["rast"]]<-admin_rast
  admin[["attributes"]]<-copy(data.table(admin_shp@data))
  
  return(admin)
}


#' Looks for missing admins that wont be picked up in the simple raster method of doing things
#' 
#' Returns a list of pixels to use for the missing admins, if there are any. 
#' Also retuns a rasterized admin for use later on.
#'
#' @param lvl admin level, 0, 1, 2
#' @param simple_raster simple raster for this region
#' @param reg region name
#' @param adm0_list list of adm0 codes for this region
#' @param cp cell pred object that missing adm values will be taken from
#' @param year_list numeric list of years in cp
#'

find_missing_adms <- function(lvl,simple_raster,reg,adm0_list, cp, year_list) {
  
  # get needed objects
  pixel_id <- seegSDM:::notMissingIdx(simple_raster)
  pixel_spatial<-data.table(pixel_id=pixel_id)
  message("Generating a raster of the pixel_id values")
  pixel_id_raster <- insertRaster(simple_raster, cbind(pixel_id))
  
  # get admin shapefile
  fieldname  <- paste0("ADM",lvl,"_CODE")
  admin_info <-GetAdmin(admin_level=lvl,simple_raster,region_adm0_list = adm0_list, modeling_shapefile_version)
  pixel_spatial[[fieldname]]<-raster::extract(admin_info[["rast"]],pixel_spatial$pixel_id) # Generate a field based on the admin boundary unit that has the ID code in the 
  
  #1st row of cell pred, used to check for NA admin units
  simple_extract <- raster::extract(simple_raster, pixel_id)
  est <- cp[,1]
  est <- data.table(cbind(rep(pixel_spatial$pixel_id,18), unname(unlist(pixel_spatial[,2])), rep(simple_extract,18), rep(year_list, each = (length(est) / length(year_list))), est))
  setnames(est, c("V1", "V2", "V3", "V4"), c("pixel_id" , paste0("ADM",lvl,"_CODE"), "ADM0_CODE", "year_list"))
  if(lvl == 0) {
    #duplicate ADM0_CODE column when level is 0
    est[,3] <- NULL
    agg_est <- est[,mean(est, na.rm=T), by=c("ADM0_CODE", "year_list")]
    missing_agg_est <- unique(agg_est[is.na(V1) & !is.na(ADM0_CODE), 1])
  } else if(lvl == 1) {
    agg_est <- est[,mean(est, na.rm=T), by=c("ADM0_CODE", "ADM1_CODE", "year_list")]
    missing_agg_est <- unique(agg_est[
      is.na(V1) & !is.na(ADM0_CODE) & !is.na(ADM1_CODE), c(1,2)
    ])
  } else if(lvl == 2){
    placeholder <-GetAdmin(admin_level=1,simple_raster,region_adm0_list = adm0_list, modeling_shapefile_version)
    placeholder_adm1<-raster::extract(placeholder[["rast"]],pixel_spatial$pixel_id)# Generate a field based on
    est <- cbind(est, rep(placeholder_adm1, 18))
    setnames(est, "V2", "ADM1_CODE")
    agg_est <- est[,mean(est, na.rm=T), by=list(ADM2_CODE, ADM1_CODE, ADM0_CODE, year_list)]
    missing_agg_est <- unique(agg_est[
      is.na(V1) & !is.na(ADM0_CODE) & !is.na(ADM1_CODE) & !is.na(ADM2_CODE) , c(1:3)
    ])
  }
  missing_admins<-list() # This will contain the cell_preds for otherwise missing admin units, by level.
  
  # First, we discover what GAUL codes are missing from the pixel_spatial table compared to the shapefile:
  shpfile_attributes <- admin_info[["attributes"]] # Get the attribute table from the admin level's shapefile
  admin_cols         <- names(shpfile_attributes)[names(shpfile_attributes)%in%c("ADM0_CODE","ADM1_CODE","ADM2_CODE")] # Any columns in ADM0, ADM1, and ADM2_CODE that should be included as merges.
  shpfile_attributes <- subset(shpfile_attributes,ADM0_CODE != 40762)
  attribute_codes    <- shpfile_attributes[[fieldname]] # Get list of codes for that level from the admin level's attribute table
  pixel_codes        <- pixel_spatial[[fieldname]] # Get list of codes based on what's in the pixel_spatial object

  missing_table<-shpfile_attributes[!(attribute_codes %in% pixel_codes),admin_cols,with=F]
  missing_table <- unique(rbind(missing_table, missing_agg_est))
  
  #drop all rows where ADM0 is not in the region
  sanity_check <- get_adm0_codes(reg)
  missing_table <- missing_table[ADM0_CODE %in% sanity_check,]
  
  #pull all unique missing codes at lvl
  missing_codes <- unname(unlist(na.omit(unique(
    missing_table[, c(paste0("ADM", lvl, "_CODE")), with=F],
    na.rm=T
  ))))  
  if(length(missing_codes)==0){
    message(paste0("No missing codes at level ",lvl))
    return(list(admin_info = admin_info))
  }else{
    message(paste0("Missing codes found for level ",lvl,":"))
    # Strategy 1: Assign a pixel location based on centroid location
    # Develop a raster of the pixel-IDs:
    message("  Discovering centroid location")
    points <- admin_info[["centroids"]] # SpatialPointsDataFrame of that level
    missing_points<-points[(points@data[[fieldname]] %in% missing_codes),] # getting only the missing points
    missing_centroid_locs<-data.table(raster::extract(pixel_id_raster,missing_points), missing_points[[fieldname]]) # Extracting the missing location values as centroid points...
    names(missing_centroid_locs)<-c("point_method", fieldname)
    if(nrow(missing_centroid_locs) == 0) {
      missing_centroid_locs <- data.table(fn = missing_codes)
      setnames(missing_centroid_locs, 'fn', fieldname)
      missing_centroid_locs$point_method <- NA 
    }
    
    missing_admins_centroids<-missing_centroid_locs
    
    point_check <- est[pixel_id %in% missing_admins_centroids$point_method,]
    point_check <- point_check[is.na(est),]
    missing_admins_centroids[point_method %in% point_check$pixel_id, point_method := NA]
    
    # Strategy 2: Assign a pixel location based on polygon extraction
    # Develop a raster of the pixel-IDs:
    message("  Discovering first raster pixel touching polygon")
    polys<-admin_info[["spdf"]] # SpatialPointsDataFrame of that level
    missing_admins_polys <- tryCatch(
      {
        missing_polys<-polys[(polys@data[[fieldname]] %in% missing_codes),] # getting only the missing points
        # Extracting the missing location values as polygons, pulling the first raster cell that it touches...
        missing_poly_locs<-data.table(
          raster::extract(x=pixel_id_raster,y=missing_polys,small=T,fun=function(x,...)first(x)), 
          missing_polys[[fieldname]]
        )
        names(missing_poly_locs)<-c("poly_method", fieldname)
        missing_poly_locs # Return the polygon locations
      },
      error = function(e){
        missing_poly_ids <- intersect(polys@data[[fieldname]], missing_codes)
        if(length(missing_poly_ids)==0) missing_poly_ids <- "NA"
        missing_poly_locs <- data.table(
          poly_method=NA,
          fn=missing_poly_ids
        )
        setnames(missing_poly_locs, 'fn', fieldname)
        return(missing_poly_locs)
      }
    )
    
    # Merging strategies together: centroids and polygons; adding to list.
    missing_locs<-merge(missing_admins_polys,missing_admins_centroids,by=fieldname)
    missing_locs<-merge(missing_table,missing_locs,by=fieldname, all.x=TRUE) # Add in admin 0, 1 levels if relevant
    missing_locs[,pixel_id:=point_method] # If centroid produced something, go with centroid
    missing_locs[is.na(point_method),pixel_id:=poly_method] # Otherwise, go with the polygon method.
    
    point_check <- est[pixel_id %in% missing_locs$pixel_id,]
    point_check <- point_check[is.na(est),]
    missing_locs[pixel_id %in% point_check$pixel_id, pixel_id := NA]
    
    # For those still missing, assign them to nearest non-NA pixel using centroid
    if(NA %in% missing_locs$pixel_id){
      message('  After centroids and polygon methods, there are still some missing admins.')
      message('  Now sampling to find nearest non-NA pixel and using that')
      
      for(rr in which(is.na(missing_locs$pixel_id))){
        message(sprintf('  -- finding nearest pixel to gaul_code: %i',
                        as.numeric(missing_locs[rr, sprintf('ADM%i_CODE', lvl), with = F])))
        
        ## get centroid of chape
        mp <- polys[polys[[sprintf('ADM%i_CODE', lvl)]] == missing_locs[[sprintf('ADM%i_CODE', lvl)]][rr],]
        cent <- getSpPPolygonsLabptSlots(mp)
        
        ## loop through withh an increasing radius and see if nearby points are non-NA
        found <- 0
        radius <- .005
        while(found != 1 & radius < 10 & prod(dim(cent)) > 0){ ## stop for max radius or match found
          ## sample 1000 nearby locs
          near <- matrix(runif(2000, -radius, radius), ncol = 2)
          near[, 1] <- near[, 1] + cent[1, 1]
          near[, 2] <- near[, 2] + cent[1, 2]
          
          ## extract raster pixels
          near <- data.table(raster::extract(pixel_id_raster, near), near)
          colnames(near) <- c('pixel_id', 'x', 'y')
          
          if(mean(is.na(near[,pixel_id])) < 1){ ## then we've found a non-NA neighbor
            found <- 1 ## end while loop

            
            ## find closest neighbor
            near <- na.omit(near) ## those with non-NA pixels
            dist <- sqrt((near[, x] - cent[1, 1]) ^ 2 +
                           (near[, y] - cent[1, 2]) ^ 2)
            min.ind <- which.min(dist) ## in case of tie, this returns 1st
            
            pix <- near[min.ind, pixel_id]
            
            if(is.na(est[pixel_id == pix, est][[1]])){
              #pixel is na in cell pred, try again
              found <- 0
            } else {
              ## take the pixel id of the nearest sampled neighbor
              message(sprintf('  ---- found neighbor using radius: %f', radius))
              missing_locs[rr, pixel_id := pix]
            }
          }
          
          ## increase radius in case we didn't catch anything
          radius <- radius + .005
          
        } ## end while loop
      } ## for each admin with NA pixel_id
    } ## if any NA pixel ids after centroid and poly methods
    
    point_check <- est[pixel_id %in% missing_locs$pixel_id,]
    point_check <- point_check[is.na(est),]
    missing_locs[pixel_id %in% point_check$pixel_id, pixel_id := NA]
    
    # Check to see if NAs still exist in missing locations, make a warning about missing areas.
    if(NA %in% missing_locs$pixel_id){
      message( "The following admin units appear to be getting lost forever:")
      print(missing_locs[is.na(pixel_id),c("pixel_id",admin_cols),with=F])
    }
    
    # Merging on locations with pixel IDs
    missing_locs<-missing_locs[!is.na(pixel_id),c("pixel_id",admin_cols),with=F]
    
    return(list(missing_locs  = missing_locs,
                pixel_spatial = pixel_spatial,
                admin_info    = admin_info))
    
  }
}


#' Grab cell pred for missing admin areas
#'
#' @param cp cell pred object to grab from
#' @param mo object output from find_missing_adms()
#' @param yl year list
#' @param ad admin level
#'
get_missing_admin_cell_pred <- function(cp,mo,yl,ad){  
  if(!'missing_locs' %in% names(mo)){
    message('No missing admins to deal with here')
    return(NULL)
  } else {
    mo$missing_locs
    
    if(nrow(cp) != (nrow(mo$pixel_spatial) * length(yl))){
      stop('CELL PREP NOT SEARCHABLE TO FILL IN MISSING ADS: nrow(cp) != nrow(mo$pixel_spatal) * length(yl)')
    }
    
    added_admin_cp <- list()
    for(missad in mo$missing_locs[[paste0('ADM',ad,'_CODE')]]){
      # get the missing add pixel from the shapefile (by row from the pixel spatial object)
      px                <- mo$missing_locs[['pixel_id']][mo$missing_locs[[paste0('ADM',ad,'_CODE')]]==missad]
      cprow_base        <- which(mo$pixel_spatial$pixel_id == px)
      cprows            <- c(cprow_base + (nrow(cp)/length(yl)) * c(0:(length(yl) -1)))
      miss_cp           <- cp[cprows,]
      rownames(miss_cp) <- c(paste0(missad,'_',yl))
      added_admin_cp[[length(added_admin_cp)+1]] <- miss_cp
    }
    
    added_admin_cp <- do.call('rbind',added_admin_cp)
    
    return(added_admin_cp)
  }
}


#' Get names of completed regions based on output director
#'
#' @param rd run date
#'
completed_regions <- function(rd=run_date){
  dir  <- sprintf("<<<< FILEPATH REDACTED >>>>")
  fns  <- grep('*_0$', list.files(dir,pattern='fin__bin0_'), value=TRUE)
  regs <- gsub(gsub(fns,pattern='^fin__bin0_',replacement=''),pattern='_0$',replacement='')
  return(regs)
}


## regions_needing_postest() -------------------------------------------------->
#' 
#' Returns a vector of regions that have finished modeling but have NOT finished
#'   U5M postestimation
#' 
#' @param rd run date
#' 
regions_needing_postest <- function(rd=run_date){
  # Check regions have finished modeling
  compregs <- completed_regions(rd=rd)
  # Check regions that have finished postestimation based on the last file that
  #  gets saved in postest
  need_postest <- character(0)
  dir_template <- paste0("<<<< FILEPATH REDACTED >>>>")
  for(ag in c('neonatal','infant','under5')){
    dir <- sprintf(dir_template, ag)
    if(dir.exists(dir)){
      fns <- grep("*.csv$", list.files(dir,pattern=sprintf('fr_%s_raking_factors',ag)),value=TRUE)
      post_fin_regs <- gsub(
        gsub(fns, pattern=sprintf('^fr_%s_raking_factors_',ag), replacement=''),
        pattern='.csv$', replacement=''
      )
      # Return regions that finished prediction but not postestimation
      need_postest <- unique(c(need_postest, compregs[!(compregs %in% post_fin_regs)]))
    } else {
      need_postest <- compregs
    }
  }
  return(need_postest)
}


#' Pull gbd estimates for q for raking
#'
#' @param age either `under5`, `infant`, or `neonatal`
#' @param year_list list of years to pull from gbd
#'
get_gbd_q <- function(age,
                      year_list) {
  
  library(ini, lib="<<<< FILEPATH REDACTED >>>>")
  library(slackr, lib="<<<< FILEPATH REDACTED >>>>")
  
  library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")
  
  if(age == "under5"){
    age_group_id <- c("1")
  } else if(age == "infant") {
    age_group_id <- c("28")
  } else if(age == "neonatal") {
    age_group_id <- c("2", "3")
  }
  
  est_5q0_final <- get_mort_outputs("with shock life table", "estimate", run_id = "best", life_table_parameter_id = 3, age_group_id = age_group_id, gbd_year = 2017)
  
  #clean up
  est_5q0_final <- est_5q0_final[sex_id == 3,]
  lt <- est_5q0_final[year_id %in% year_list,]
  
  if(age == "neonatal"){
    lt <- lt[, .(mean = 1 - prod(1 - mean, na.rm=T)), by = c( "year_id", "location_id")]
    
    #pull neonatal from flat file for now, leaving old code for future when NN is added to db
    lt <- fread("<<<< FILEPATH REDACTED >>>>")
    lt <- lt[year_id %in% year_list,]
    lt <- lt[age_group_id == 42,]
    lt <- lt[sex_id == 3,]
  }
  
  #clean up
  setnames(lt, c("location_id", "year_id"), c("name", "year"))
  lt <- lt[,c("name", "year", "mean")]
  
  return(lt)
}


#' Fractionally aggregate populations to shapefile
#'
#' @description: Builds a link table out of the given shapefile and field, and aggregates the given covariates to that shapefile. 
#' @depends sp, sf, raster, rgeos, rgdal, data.table, and dplyr packages
#'
#' @param shapefile_path Path to shapefile to aggregate to
#' @param field Field in shapefile to aggregate to (like GAUL_CODE)
#' @param covs a vector of covariates to pass to `load_and_crop_covariates_annual`. Must be the same length as `measures`
#' @param measures a vector of measures to pass to `load_and_crop_covariates_annual`. Must be the same length as `covs`
#' @param year_list years to pull covariates for
#' @param interval_mo interval month to pass to `load_and_crop_covariates_annual`. Multiple values not supported. 
#' @param cores number of cores to run mclapply over - not implemented yet
#' @param link_table a link table that matches the shapefile_path and field. Used to bypass `build_custom_link_table`.`
#' @param id_raster an id_raster that matches the shapefile_path and field. Used to bypass `build_custom_link_table`.
#' 
#' @return returns the aggregated covariate data (agg_table), link table (link_table), and the id raster(id_raster) in a named list 
#'
fr_aggregate_pop_custom_shapefile <- function(shapefile_path,
                                              field,
                                              covs,
                                              measures,
                                              year_list,
                                              interval_mo,
                                              cores,
                                              link_table = NULL,
                                              id_raster = NULL) {
  
  #build link table if not passed in to function call
  if(is.null(link_table) & is.null(id_raster)) {
    message("Building link table for provided shapefile. Region or Global shapefiles may take several hours to run.")
    link_output <- build_link_table(shapefile_path=NULL, field, cores)
    link_table <- link_output$link_table
    id_raster <- link_output$id_raster
  }
  
  #loop through covariates, because load_and_crop can't do multiple measures of the same covariate
  agg_list = list()
  for(i in 1:length(covs)) {
    pop_raster_annual <- load_and_crop_covariates_annual(covs           = covs[i],                
                                                         measures       = measures[i],      
                                                         simple_polygon = id_raster,
                                                         start_year     = min(year_list),
                                                         end_year       = max(year_list),
                                                         interval_mo    = as.numeric(interval_mo),
                                                         agebin=1)[[1]]
    
    #get pixel_id, year, and covariate values into vectors of the same length
    id_raster_extract <- raster::extract(id_raster, extent(id_raster))
    year_extract <- rep(year_list, each = length(id_raster_extract))
    id_raster_extract <- rep(id_raster_extract, times = length(year_list))
    pop_raster_extract <- as.vector(raster::extract(pop_raster_annual, extent(pop_raster_annual)))
    
    if(length(id_raster_extract) != length(pop_raster_extract)) {
      message("cov raster and id raster have different length, this should never happen")
    }
    
    #make a table of pixel_id, year, and covariate values
    pop_table <- data.table(cbind(id_raster_extract, year_extract, pop_raster_extract))
    
    #merge onto link table by pixel_id (this duplicates rows for each year)
    link_pop <- merge(link_table, pop_table, by.x = "pixel_id", by.y = "id_raster_extract", all.x = T)
    
    #check the duplication didn't go wrong
    if(nrow(link_table) != (nrow(link_pop) / length(year_list))) {
      message("problem merging covariate and link table, extra rows added")
    }
    
    #calculate fractional covariate value
    link_pop[, frac_pop := pop_raster_extract * area_fraction]
    #aggregate to field and year
    agg_table <- link_pop[, .(pop = sum(frac_pop, na.rm = T)), by = c(field, "year_extract")]
    setnames(agg_table, "pop", paste0(covs[i], "_", measures[i]))
    agg_list[[paste0(covs[i], "_", measures[i])]] <- agg_table
  }
  
  #combine covariates into a single table
  if(length(agg_list) > 1){
    agg_table <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = c(field, "year_extract")), agg_list)
  } else {
    agg_table <- agg_list[[1]]
  }
  
  setnames(agg_table, c("year_extract"), "year")
  
  return(list("agg_table" = agg_table, "link_table" = link_table, "id_raster" = id_raster))
}


get_is_oos_draws_u5m_old <- function(ind_gp,
                                 ind,
                                 rd,
                                 agebin, # 'neonatal','infant', 'under5'
                                 regions,
                                 abfracs = 'ab_fracs',
                                 ind_fm = 'binomial',
                                 yrs = 2000:2017,
                                 write.to.file = FALSE,
                                 get.oos       = FALSE,
                                 year_col      = 'year'){
  
  
  if (class(yrs) == "character")  yrs  <- eval(parse(text=yrs))
  nperiod <- length(yrs)
  
  if(length(regions)>1)
    stop('NOT MORE THAN ONE REGION PLEASE')
  
  ###############
  ## Load data ##
  ###############
  
  message("Load input data used in model")
  
  # load regions that were used in modeling
  mod.dir <- sprintf("<<<< FILEPATH REDACTED >>>>")
  
  # load raw data
  df <- data.table()
  if(!get.oos) {
    for(r in regions){
      tmp <- fread("<<<< FILEPATH REDACTED >>>>")
      tmp$region <- r
      tmp <- tmp[age %in% abfracs[[agebin]],]
      df <- rbind(df, tmp)
    }
    df <- merge_with_ihme_loc(df)
    df <- df[GAUL_CODE %in% get_gaul_codes(r)] # removing countries not in region
  } else {
    stop('OOS NOT YET SET UP')
    #df <- rbindlist(lapply(readRDS(paste0(mod.dir, 'stratum.rds')), data.table))
  }
  
  message(sprintf('DF has %i rows',nrow(df)))
  # rename year column for convenience
  setnames(df, year_col, "the_year_col")
  
  ###################
  ## Assign admins ##
  ###################
  
  message('Identify ad1 and ad2 membership')
  
  # load admin2 shapefile (also contains admin1)
  message('WARNING USING GAUL STILL....')
  admin2_shapefile <- rgdal::readOGR(get_admin_shapefile(admin_level=2))
  for (v in grep("CODE", names(admin2_shapefile@data))) admin2_shapefile@data[[v]] <- as.numeric(as.character(admin2_shapefile@data[[v]]))
  
  # identify the admin2 (and by extension, the admin1) each point belongs to
  locs <- SpatialPoints(cbind(df$longitude, df$latitude), proj4string = CRS(proj4string(admin2_shapefile)))
  adm.df <- sp::over(locs, admin2_shapefile)
  
  # for those that don't fall inside a polygon, assign them to the nearest polygon (this tends to happen on coastlines)
  # do this by country to speed things up and to make sure the point ends up at least in the correct admin0
  for (ctry in unique(df[is.na(adm.df[,1]), GAUL_CODE])) {
    ii <- which(is.na(adm.df[,1]) & df$GAUL_CODE == ctry)
    temp_shape <- admin2_shapefile[admin2_shapefile@data$ADM0_CODE == ctry,]
    distmat <- gDistance(locs[ii], temp_shape, byid = T)
    jj <- apply(distmat, 2, which.min)
    adm.df[ii,] <- temp_shape@data[jj,]
    rm(ii, jj, temp_shape, distmat)
  }
  
  # copy admin information into df
  df$ad1 <- adm.df$ADM1_CODE
  df$ad2 <- adm.df$ADM2_CODE
  df$ad0 <- df$GAUL_CODE # for consistency...
  
  
  
  ###############
  ## Get draws ##
  ###############
  
  message('Get draws')
  
  # loop over regions
  #  df_all <- rbindlist(lapply(all.regions, function(rr) {
  rr <- regions
  
  message(paste('...Region:', rr))
  
  # load the simple raster template
  message('......load simple raster template')
  gaul_list <- get_gaul_codes(rr)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 0.4, subset_only = TRUE)
  subset_shape <- simple_polygon_list[[1]]
  raster_list <- build_simple_raster_pop(subset_shape)
  template <- raster_list[['simple_raster']]
  
  # subset data to this region
  df.r <- df[region == rr,]
  loc.r <- as.matrix(df.r[, list(longitude, latitude)])
  
  # make a raster of cell_pred row ids
  id_raster <- insertRaster(template, matrix(1:(length(cellIdx(template)) * nperiod), ncol = nperiod))
  
  # get cell_pred row ids for each data point
  for (yy in 1:nperiod) {
    this_year <- which(df.r$the_year_col == yrs[yy])
    if (length(this_year) == 0) next
    df.r[this_year, cell_pred_id := raster::extract(id_raster[[yy]], loc.r[this_year,])]
  }
  
  df.r <- df.r[!is.na(cell_pred_id),]
  
  
  # if out of sample metrics are requested, duplicate the data to create separate rows for in and out of sample
  if (get.oos) {
    df.r <- rbind(df.r, cbind(df.r[, -"fold", with=F], fold = 0), use.names = T)
  } else {
    df.r[, fold := 0]
  }
  
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  # loop over holdouts
  for(a in abfracs[[agebin]]){
    for (this_fold in sort(unique(df.r$fold))) {
      
      # load cell pred objects
      message(sprintf('......load cell preds for age %i and holdout %i', a, this_fold))
      tmpmoddir <- sprintf("<<<< FILEPATH REDACTED >>>>")
      cell_pred <- loadRData("<<<< FILEPATH REDACTED >>>>")
      
      # extract draws
      df.r[fold == this_fold & age == a, paste0("draw", 1:ncol(cell_pred)) := as.data.table(cell_pred[cell_pred_id,])]
      
    }
  }
  
  # rename year column back to original
  setnames(df.r, "the_year_col", year_col)
  
  # save combined data and draws to file and return
  outputfilename <- sprintf("<<<< FILEPATH REDACTED >>>>")
  sprintf('Saving file as %s',outputfilename)
  if (write.to.file) saveRDS(df.r, outputfilename)
  return(df.r)
  
}



get_is_oos_draws_u5m <- function(ind_gp,
                                 ind,
                                 rd,
                                 region,
                                 data_sample, # can be 'IS', 'OOS', or a vector of both
                                 abfracs, # list describing which analytical ages go into each agebin (produced in launch_u5m.R)                                 
                                 agebin        = 'under5', # 'neonatal','infant', 'under5' # JUST DO UNDER5 to get all age bins. Need all later for PV table anyway

                                 ind_fm        = 'binomial',
                                 yrs           = 2000:2017,
                                 write.to.file = TRUE,
                                 fold_table    = NULL,
                                 year_col      = 'year'){
  
  
  if (class(yrs) == "character")  yrs  <- eval(parse(text=yrs))
  nperiod <- length(yrs)
  
  if(length(region)>1)
    stop('NOT MORE THAN ONE REGION PLEASE')
  
  # detect fold table and if doing IS, OOS, or both
  if(is.null(fold_table) & 'IS' %in% data_sample & length(data_sample) == 1) {
    message('IN SAMPLE ONLY')
  } else if (all(data_sample == c('IS','OOS')) & !is.null(fold_table)) {
    message(sprintf('IN SAMPLE AND OUT OF SAMPLE WITH FOLDS: %s', paste0(sort(unique(fold_table$fold)),collapse=', ')))
  } else if (!is.null(fold_table) & 'OOS' %in% data_sample & length(data_sample) == 1) {
    message(sprintf('OUT OF SAMPLE ONLY WITH FOLDS: %s', paste0(sort(unique(fold_table$fold)),collapse=', ')))
  } else {
    stop('data_sample AND fold_table ARGUMENTS DO NOT SEEM TO AGREE ON IF THIS IS "IS","OOS",or c("IS,"OOS"). CHECK INPUTS')
  }
  
  if(is.null(fold_table) & 9999 %in% fold_table$fold) fold_table[fold==9999, fold := NA]
  
  ###############
  ## Load data ##
  ###############
  
  message("Load input data used in model")
  
  # load regions that were used in modeling
  mod.dir <- sprintf("<<<< FILEPATH REDACTED >>>>")

  # get admin codes mapping
  adm_to_loc_id <- get_location_code_mapping(shapefile_version = modeling_shapefile_version)
  
  ### load raw data
  df <- data.table()
  
  if('IS' %in% data_sample) {
    message('LOADING IN-SAMPLE RAW DATA')

    tmp <- fread(sprintf("<<<< FILEPATH REDACTED >>>>") # use training data for in-sample
    tmp$region <- region
    tmp <- tmp[age %in% abfracs[[agebin]],]
    
    # removing countries not in region (BORDER DATA)
    tmp <- tmp[, ihme_lc_id := as.character(country)]
    tmp <- merge(tmp, adm_to_loc_id, by='ihme_lc_id',all.x=T)
    tmp <- tmp[ADM_CODE %in% get_adm0_codes(region, shapefile_version = modeling_shapefile_version)] 
    
    # in sample gets fold zero
    tmp$fold <- 0
    
    # rbind it on to the full dataset
    tmp <- tmp[,c('nid','ihme_lc_id',year_col,'age','data_type',ind,'N','weight','fold','latitude','longitude','ADM_CODE'),with=FALSE]
    df  <- rbind(df, tmp)
  } 
  
  if('OOS' %in% data_sample) {
    message('LOADING OUT-OF-SAMPLE RAW DATA')
    
    tmp <- fread("<<<< FILEPATH REDACTED >>>>") # use full dataset for oos (all the same at this point, pre drops from folds)
    # Drop fold 9999 (workaround for single-holdout runs)
    tmp[ fold == 9999, fold := NA ]
    tmp$region <- region
    tmp <- tmp[age %in% abfracs[[agebin]],]

    # removing countries not in region (BORDER DATA)
    tmp <- tmp[, ihme_lc_id := as.character(country)]
    tmp <- merge(tmp, adm_to_loc_id, by='ihme_lc_id',all.x=T)
    tmp <- tmp[ADM_CODE %in% get_adm0_codes(region, shapefile_version = modeling_shapefile_version)] 
    
    # assing folds
    tmp <- merge(tmp, fold_table, by = 'nid', all.x = TRUE)
    
    # rbind it on to the full dataset
    tmp <- tmp[,c('nid','ihme_lc_id',year_col,'age','data_type',ind,'N','weight','fold','latitude','longitude','ADM_CODE'),with=FALSE]
    df  <- rbind(df, tmp)
  } 
  
  
  message(sprintf('DF has %i rows',nrow(df)))
  # rename year column for convenience
  setnames(df, year_col, "the_year_col")
  
  ###################
  ## Assign admins ##
  ###################
  
  message('Identify ad1 and ad2 membership')
  
  # load admin2 shapefile (also contains admin1)
  admin2_shapefile <- rgdal::readOGR(get_admin_shapefile(admin_level=2,version=modeling_shapefile_version))
  for (v in grep("CODE", names(admin2_shapefile@data))) admin2_shapefile@data[[v]] <- as.numeric(as.character(admin2_shapefile@data[[v]]))
  

  # identify the admin2 (and by extension, the admin1) each unique point belongs to
  ucoords  <- unique(data.table(longitude=df$longitude, latitude=df$latitude))[,ucid:=1:.N] # do unique for sppedups on the sp::over call
  locs     <- SpatialPointsDataFrame(cbind(ucoords$longitude, ucoords$latitude),
                            data        = ucoords,
                            proj4string = CRS(proj4string(admin2_shapefile)))
  adm.df   <- sp::over(locs, admin2_shapefile) # MAKE SURE THIS RETURNS IN ORDER OF locs OBJECT
  
  
  # for those that don't fall inside a polygon, assign them to the nearest polygon (this tends to happen on coastlines)
  # do this by country to speed things up and to make sure the point ends up at least in the correct admin0
 # for (ctry in unique(df[is.na(adm.df[,1]), ADM0_CODE])) {
  ii          <- which(is.na(adm.df[,1]))
  temp_shape  <- admin2_shapefile[admin2_shapefile@data$ADM0_CODE %in% unique(df$ADM_CODE),] # keep only the region
  distmat     <- gDistance(locs[ii,], temp_shape, byid = TRUE)
  jj          <- apply(distmat, 2, which.min)
  adm.df[ii,] <- temp_shape@data[jj,]
  rm('distmat')

  # merge admins back on to df based on coordinates
  adm.df     <- cbind(adm.df, ucoords)
  
  # CONFIRM ALL ORDERING IS CORRECT
  #png("<<<< FILEPATH REDACTED >>>>")
  #ggplot(df[sample(1:nrow(df),50000)],aes(x=longitude,y=latitude,color=ihme_lc_id))+theme_minimal()+geom_point()
  #dev.off()
  
  adm.df.tmp <- data.table(adm.df)[,c('ADM0_CODE','ADM1_CODE','ADM2_CODE','ADM1_NAME',
                                      'ADM2_NAME','longitude','latitude'),with=FALSE]
  df <- merge(df,adm.df.tmp,by=c('latitude','longitude'),all.x=TRUE)

  # copy admin information into df
  setnames(df,paste0('ADM',0:2,'_CODE'),paste0('ad',0:2))


  ###############
  ## Get draws ##
  ###############
  
  message('DONE GETTING RAW DATA, NOW GETTING DRAWS')
  
  rr <- region # to avoid namespace pollution
  message(paste('...Region:', region))
  
  # load the raster template to base pixel id off of
  # any output raster from this region will do
  # search for file to open through holdouts
  aaa <- 0
  ffexist <- FALSE
  while(ffexist == FALSE & aaa <= 5){
    ff <- sprintf("<<<< FILEPATH REDACTED >>>>")
    ffexist <- file.exists(ff)
    aaa <- aaa + 1
  }
  template <- brick(ff)[[1]]
  

  # coords
  loc.r <- as.matrix(df[, list(longitude, latitude)])
  
  # make a raster of cell_pred row ids
  id_raster <- insertRaster(template, matrix(1:(length(cellIdx(template)) * nperiod), ncol = nperiod))
  
  # get cell_pred row ids for each raw data point in the main df
  for (yy in 1:nperiod) {
    this_year <- which(df$the_year_col == yrs[yy])
    if (length(this_year) == 0) next
    df[this_year, cell_pred_id := raster::extract(id_raster[[yy]], loc.r[this_year,])]
  }
  
  # drop any data rows with no estimate
  df <- df[!is.na(cell_pred_id),]
  
  
  loadRData <- function(fileName){
    #loads an RData file, and returns it
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  
  # loop over holdouts
  for(a in abfracs[[agebin]]){
    for (this_fold in na.omit(sort(unique(df$fold)))) {
      
      # load cell pred objects
      message(sprintf('..load cell preds for age %i and holdout %i', a, this_fold))
      tmpmoddir <- sprintf("<<<< FILEPATH REDACTED >>>>")
      cell_pred <- loadRData("<<<< FILEPATH REDACTED >>>>")
      
      cell_pred <- na.omit(cell_pred) # remove NAs since those are from the simple_raster, and here we are using model output as the template, so NA rows must be removed. 
      
      message(sprintf('Number of rows in cell_pred: %i, number of pixels in pixel_id raster: %i', 
                      nrow(cell_pred), max(as.vector(id_raster),na.rm=TRUE)))
      
      # extract draws
      df[fold == this_fold & age == a, paste0("draw", 1:ncol(cell_pred)) := as.data.table(cell_pred[cell_pred_id,])]
      
      # save the mem, drop cell_pred
      rm(cell_pred)
    }
  }
  
  # rename year column back to original
  setnames(df, "the_year_col", year_col)
  
  # save combined data and draws to file and return
  outputfilename <- sprintf("<<<< FILEPATH REDACTED >>>>")
  sprintf('Saving file as %s',outputfilename)
  if (write.to.file) saveRDS(df, outputfilename)
  return(df)
  
}







## get_pv_table_u5m ################################################

#' Generate plots and tables of metrics of predictive validity
#'
#' @param age_bin to look for the folder containing all the regional runs 
#' @param indicator indicator
#' @param indicator_group indicator_group
#' @param rd run date
#' @param aggregate_on column of spatial aggregation ("country", "ad1", or "ad2" usually) - can pass vector
#' @param collapse_to_agebin collapse to age bin based on abfracs
#' @param coverage_probs probability at which you want to assess coverage of predictive intervals
#'                       can be a single number or a vector; if a vector can produce calibration plots
#'                       assessing x% coverage at various cutoffs
#' @param result_agg_over what should the final table aggregate over (vector)? For instance
#'                        c("year", "oos") are the defaults.  If you include "region" here
#'                        the function will produce a separate set of validation plots,
#'                        one for each region
#' @param weighted do you want to weight PV metrics on sample size? (Boolean)
#' @param family distribution to use ("binomial" or "gaussian")
#' @param plot produce plots? (Boolean)
#' @param plot_by produce separate plots for an element of `result_agg_over`?
#'                for instance, `plot_by = region` with `result_agg_over = c("year", "oos", "region"`
#'                will produce a separate plot for each region/oos combo
#' @param plot_by_title title for captions for "plot_by" object, i.e. "Region" (defaults to plot_by)
#' @param plot_ci plot CIs (Boolean)
#' @param plot_ci_level what ci level to plot (numeric, i.e. "95" is the 95% UI)
#' @param ci_color what color do you want the CI lines to be?
#' @param point_alpha how transparent do you want the points themselves to be (numeric, [0-1])
#' @param point_color what color do you want the points to be in your validation plots?
#' @param plot_title main title for your plots (defaults to "indicator")
#' @param plot_ncol how many columns do you want your plots to have?
#' @param save_csv save a csv table of your predictive validity metrics?
#' @param abfracs list of ages that goes into the agebins
#' @param out.dir where do you want these results to be saved?
#'
#' @return data.table of PV metrics
#' @export plots and tables of PV metrics depending on the combinations of options above
#'
#' @examples
#'
get_pv_table_u5m <- function(agebin, # if want to collapse (neonatal, infant, or under5)
                             rd,
                             collapse_to_agebin,
                             indicator = 'died',
                             indicator_group = 'u5m',
                             aggregate_on = c('ad0','ad1','ad2'),
                             coverage_probs = seq(5,95,by=5),
                             result_agg_over =  c('oos'), #'year'
                             weighted = TRUE,
                             family = 'binomial',
                             plot = TRUE,
                             plot_by = NULL,
                             plot_by_title = NULL,
                             plot_ci = FALSE,
                             plot_ci_level = 95, 
                             ci_color = "grey",
                             point_alpha = 1,
                             point_color = "black",
                             plot_title = indicator,
                             plot_ncol = 4,
                             save_csv = T,
                             abfracs = ab_fracs) {
  
  
  out.dir <- sprintf("<<<< FILEPATH REDACTED >>>>")
  
  require(boot)
  str_match <- stringr::str_match
  
  # load in the data
  message('Combining regional files')
  datdir <- sprintf("<<<< FILEPATH REDACTED >>>>") # always pull under5 because it has all the age bins in it
  files  <- list.files(path=datdir,pattern ='output_draws_data_')
  d <- rbindlist(
    mclapply(files, function(f) readRDS("<<<< FILEPATH REDACTED >>>>")[,1:168], mc.cores=4), 
    use.names=T, fill=T
  )
    
   
  #number of draws, detect
  draws <- sum(grepl('draw',names(d)))
  
  ## Collapse data and draws spatially
  message('Collapse data and draws spatially')
  d[ is.na(fold), fold := 0 ]
  d[, oos := (fold != 0)]
  d[, p := get(indicator) / N]
  #d[data_type=='sbh',N:=N*.1] old heuristic
  d[, exposure := N]
  #d[, exposure := N * weight] # all weights should be 1 now anyway 

  # Create list of pv metrics for all levels of aggregation in `aggregate_on`
  message('Outputting predictive validity metrics for your models at each level of aggregation')
  pv_table_list <- lapply(aggregate_on, function(agg_on) {
    message(paste0("  ", agg_on))
    by_vars <- unique(c('oos', 'year', result_agg_over, agg_on,'age','nid')) # 21 jan 2019 RB: adding in NID here
    bv2 <- by_vars[by_vars!='age']
    collapse_vars <- c('p', paste0('draw', 1:draws)) # RB 8/29/2018 note for now removing  paste0('clusters_covered_', coverage_probs)
    
    # for coverage, we need to simulate predictive draws of deaths here. (about a minute)
    # note if using 100, should use more draws here in future
    Y_sim <- data.table(sapply(1:draws, function(i) rbinom(nrow(d), round(d$exposure), d[[paste0('draw', i)]])))
    # copy d for aggregation but with Ys
    colnames(Y_sim) <- paste0('draw',1:draws)
    Y_sim <- cbind(Y_sim,d[,c('exposure','died','p',by_vars),with=FALSE])
    
    
    # collapse clusters to the agg on level (P)
    resP <- d[!is.na(draw1),
             c(list(total_clusters = .N, exposure = sum(exposure), died = sum(died)),
               lapply(.SD, function(x) weighted.mean(x, exposure))),
             keyby = by_vars, .SDcols = collapse_vars]
    
    resY <- Y_sim[!is.na(d$draw1),
              c(list(total_clusters = .N, exposure = sum(exposure), died = sum(died)),
                lapply(.SD, function(x) sum(x))),
              keyby = by_vars, .SDcols = collapse_vars]
    
    # HERE COLLAPSE TO THE AGE BIN (ie collapse 1:3 to be infant) # using under5 prepped data can do all age bins. can also just run parallel over age bins
    if(collapse_to_agebin == TRUE){
      resP <- resP[age %in% abfracs[[agebin]],]
      resP <- resP[,c(list(A = .N, total_clusters = sum(total_clusters), exposure = max(exposure)), 
                 lapply(.SD, function(x) { 1 - prod(1-x) })), # note max exposures give number of births
               by = bv2, .SDcols = collapse_vars]
      resP <- resP[A==length(abfracs[[agebin]])] # only keep ones that actually evaluated the correct number of ages
      resP$age <- agebin
      
      resY <- resY[age %in% abfracs[[agebin]],]
      resY <- resY[,c(list(A = .N, total_clusters = sum(total_clusters), exposure = max(exposure), died=sum(died)), 
                      lapply(.SD, function(x) { sum(x) })), # note max exposures give number of births
                   by = bv2, .SDcols = collapse_vars]
      resY <- resY[A==length(abfracs[[agebin]])] # only keep ones that actually evaluated the correct number of ages
      resY$age <- agebin
      resY$p <- NULL # meaningless as a sum
      
    }
    
    # Get coverage for each aggregated adm point
    ## Coverage
    message('Calculate coverage')
    x <- as.matrix(resY[,paste0('draw',1:draws),with=FALSE])
    
    one_side <- (1 - coverage_probs/100)/2
    ui <- t(apply(x, 1, quantile, c(one_side, 1 - one_side), na.rm = T)) # most are 0, 0. this is a problem...
    for (c in coverage_probs) {
      c_id <- which(c == coverage_probs)
      resY[, paste0('clusters_covered_', c) := died %between% list(ui[, c_id], ui[, c_id + length(coverage_probs)])]
    }
    
    # merge on coverage to resP for finalizing table
    resP <- cbind(resP, resY[,paste0('clusters_covered_',coverage_probs),with=FALSE])
    
    # metrics
    resP[, mean_draw := rowMeans(.SD, na.rm=TRUE), .SDcols = paste0('draw', 1:draws)]
    resP[, error := p - mean_draw]
    resP[, abs_error := abs(error)]
    
   
    ## Collapse to calculate predictive validity metrics
    weighted.rmse <- function(error, w) sqrt( sum(  (w/sum(w)) * ((error)^2) ) )
    if (weighted) resP$weight <- resP$exposure else resP$weight <- 1
    res2 <- resP[, c(lapply(.SD, function(x) weighted.mean(x, weight)),
                    rmse = weighted.rmse(error, weight),
                    median_SS = median(exposure),
                    cor = corr(cbind(p, mean_draw), weight)),
                by = c(result_agg_over, 'age'), #'age'
                .SDcols = c("error", "abs_error", "mean_draw", "p")]
    res2 <- cbind(res2,
              resY[, c(lapply(.SD, mean(x))),
                               by = c(result_agg_over, 'age'), #'age'
                               .SDcols = c(paste0("clusters_covered_", coverage_probs))])
    
    ## Clean up names
    setnames(res2, c("error", "abs_error", "p", paste0("clusters_covered_", coverage_probs)),
             c("me", "mae", "mean_p",paste0("coverage_", coverage_probs)))
    
    return(list(resP = resP, res2 = res2))
  })
  
  
  names(pv_table_list) <- aggregate_on
 
  if(collapse_to_agebin==FALSE) plot_by = 'age'
  
  # Make plots
  if(plot==TRUE){
    message('Making plots of aggregated data and estimates')
    
    # Get unique levels of `plot_by` and set up a plot_by title if needed
    if (!is.null(plot_by)) {
      plot_by_levels <- unique(pv_table_list[[1]][["resP"]][, get(plot_by)])
    }  else {
      plot_by_levels <- NULL
    }
    
    if (is.null(plot_by_title)) plot_by_title <- plot_by
    
    # Create a table of things to plot
    plot_table <- CJ(aggregate_on = aggregate_on,
                     oos = unique(pv_table_list[[1]][["resP"]]$oos),
                     plot_by_value = if (is.null(plot_by_levels)) NA else plot_by_levels)
    
    
    message("...saving plots here: ", out.dir)
    # Loop over plots
    for (i in 1:nrow(plot_table)) {
      
      # Grab items from plot table
      agg_on <- plot_table[i, aggregate_on]
      oosindic <- plot_table[i, oos]
      pb_val <- plot_table[i, plot_by_value]
      
      # Set up titles
      if (agg_on == "ad0") agg_title <- "Admin 0 (National)"
      if (agg_on == "ad1") agg_title <- "Admin 1"
      if (agg_on == "ad2") agg_title <- "Admin 2"
      
      res <- pv_table_list[[agg_on]][["resP"]]
      res2 <- pv_table_list[[agg_on]][["res2"]]
      
      # Make a validation plot -----------------------------------------------------
      
      # Set up filename and file
      plot_filename <- paste0(indicator, '_validation_plot_',
                              paste(c(as.character(agg_on),
                                      setdiff(result_agg_over, "oos")),
                                    collapse="_"), "_",
                              ifelse(oosindic, "OOS", "IS"),
                              ifelse(is.na(pb_val), "", paste0("_", pb_val)),
                              '.png')
      message(paste('    ', plot_filename))
      png("<<<< FILEPATH REDACTED >>>>", width = 12, height = 12, units = 'in', res = 350)
      
      # Subset data
      fdata <- res[oos == oosindic,]
      if (!is.na(pb_val)) {
        setnames(fdata, plot_by, "plot_by_column") # convenience
        fdata <- fdata[plot_by_column == pb_val,]
      }
      
      
      # Set up CI bar limits; range as defaults
      if (plot_ci) {
        fdata[, upper := apply(.SD, 1, quantile, p = 0.01*(plot_ci_level + (100 - plot_ci_level)/2), rm.na=TRUE), .SDcols = paste0('draw', 1:draws)]
        fdata[, lower := apply(.SD, 1, quantile, p = 0.01*((100 - plot_ci_level)/2), rm.na=TRUE), .SDcols = paste0('draw', 1:draws)]
        limits <- c(0,quantile(c(fdata$p, fdata$mean_draw, fdata$lower, fdata$upper),c(0.98))) #fdata[, range(c(p, mean_draw, lower, upper))]
      } else {
        limits <- c(0,quantile(c(fdata$p, fdata$mean_draw),c(0.98)))
      }
      
      # The plot code itself
      fdata$N <- fdata$weight
      gg <- ggplot(fdata, aes(x = p, y = mean_draw, size=sqrt(N))) + #, size = log(N))) +
        geom_abline(intercept=0, slope=1, color = 'red') +
        geom_point(colour = point_color, aes(alpha = sqrt(N))) +
        scale_size_area() +
        xlim(limits) +
        ylim(limits) +
        coord_equal() +
        theme_bw() +
        theme(strip.background = element_rect(fill="white")) +
        labs(x = 'Data Estimate',
             y = 'Mean Prediction',
             size = "sqrt(N=Births)",
             alpha = 'sqrt(N=Births)',
             title = paste0("Validation Plot for ", plot_title, " by ", agg_title),
             subtitle = paste0("OOS: ", oosindic, ifelse(is.na(pb_val), "", paste0(" | ", plot_by_title, ": ", pb_val))))
      
      if(plot_ci) {
        gg <- gg +  geom_errorbar(aes(ymin = lower, ymax = upper), colour = point_color, width = 0, size = .3, alpha = min(point_alpha, 0.2))
      }
      if (length(setdiff(result_agg_over, 'oos')) > 0) {
        gg <- gg + facet_wrap(as.formula(paste("~", paste(setdiff(result_agg_over, c('oos', plot_by)), collapse = "+"))),
                              ncol = plot_ncol)
      }
      
      plot(gg)
      dev.off()
    }
    
    
    # Make a COVERAGE calibration plot -----------------------
    if (plot == T & length(coverage_probs) > 1) {
      
      # Set up a subset of the data for plotting
      fdata <- res2[, unique(c("oos", result_agg_over, plot_by,paste0("coverage_", coverage_probs))), with=F]
      #fdata <- fdata[oos == oosindic]
      #if (!is.na(pb_val)) {
      if(collapse_to_agebin==FALSE) {
        setnames(fdata, plot_by, "plot_by_column") # convenience
       # fdata <- fdata[plot_by_column == pb_val,]
      } 
      fdata <- melt(fdata, id.vars = names(fdata)[!grepl("coverage", names(fdata))],
                    value.name = "observed_coverage",
                    variable.name = "coverage")
      fdata[, expected_coverage := as.numeric(gsub("coverage_", "", coverage))]
      fdata[, observed_coverage := observed_coverage * 100]
      #fdata$group <- apply(fdata[, result_agg_over[!(result_agg_over %in% c('oos', 'plot_by_column', plot_by))], with=F], 1, paste, collapse=" ")
      #if (sum(!is.na(fdata$group)) == 0) fdata[, group := "All"]
      if(collapse_to_agebin==FALSE) {
        fdata[, agebin := factor(plot_by_column)]
      } else {
        fdata$agebin <- agebin
      }
      
      # Create filename
      
      cplot_filename <- paste0(indicator, '_calibration_plot_',
                               paste(c(as.character(agg_on),
                                       setdiff(result_agg_over, "oos")),
                                     collapse="_"), "_",
                               ifelse(oosindic, "OOS", "IS"),'_',
                               ifelse(collapse_to_agebin==FALSE, "allages", agebin),
                               '.png')
      message(paste('    ', cplot_filename))
      png("<<<< FILEPATH REDACTED >>>>", width = 12, height = 12, units = 'in', res = 350)
      
      # Set limits and plot
      limits <- fdata[, range(c(observed_coverage, expected_coverage))]
      gg <- ggplot(fdata, aes(x = expected_coverage, y = observed_coverage, group = agebin, color = agebin)) +
        geom_abline(intercept = 0, slope = 1, color = "red") +
        geom_point() +
        geom_line(alpha = 0.2) +
        scale_color_discrete(name = "") +
        coord_equal() +
        xlim(limits) +
        ylim(limits) +
        theme_bw()+ facet_wrap(~oos) +
        labs(x = "Expected coverage", y = "Observed coverage")
      
      plot(gg)
      dev.off()
      
    } # End calibration plot if statement
  } # End if (plot==T) loop
  
  # Format, save (if desired), and return `res2` objects
  output_list <- lapply(aggregate_on, function(agg_on) {
    res2 <- copy(pv_table_list[[agg_on]][["res2"]])
    setorderv(res2, result_agg_over)
    #setnames(res2, result_agg_over, ifelse(result_agg_over == 'oos', 'OOS', gsub("(^.)", "\\U\\1", result_agg_over, perl=T)))
    setnames(res2, c("me", "mae", "mean_draw", "mean_p", "rmse", "median_SS", "cor", paste0("coverage_", coverage_probs)),
             c("Mean Err.", "Mean Abs. Err.", "Mean Pred.", "Mean Obs.", "RMSE", "Median SS", "Corr.", paste0(coverage_probs, "% Cov.")))
    
    # Save final tables if desired
    if (save_csv) {
      a_pathaddin <- ifelse(length(setdiff(result_agg_over, 'oos')) > 0,
                            paste0("_by_", paste0(setdiff(result_agg_over, 'oos'), collapse = "_")),
                            "")
      b_pathaddin <- ifelse(collapse_to_agebin,agebin,'allages')
      filename <- paste0("<<<< FILEPATH REDACTED >>>>")
      
      message(paste0("Saving csv to ", filename, "..."))
      write.csv(res2, file = filename)
    }
    return(res2)
  })
  
  names(output_list) <- aggregate_on
  
  return(output_list)
}


#' Generate plots and table of mean survey-level random effect bias adjustment
#'
#' @description uses the training data and the values of the survey-level random effects to calculate the mean effect of the random effects on the mean q of each survey. Calculates mean q for each survey, converts to logit space, adds on the random effect value, converts back to linear space, then takes (new/old - 1) * 100 to get percent increase or decrease. 
#' 
#' @param regions character vector, list of regions in a model run
#' @param run_date character, model run date
#' @param plot boolean, make plots showing distribution of mean effect?
#' @param plot_path character, place to save plots
#' @param cores int, number of cores, recommend length(regions)
#' 
#' @return data.table of mean source bias adjustment
#' @export plots given `plot` is T
#'
pull_u5m_source_bias_adjustment <- function(regions, 
                                            run_date,
                                            plot = F,
                                            plot_path = sprintf("<<<< FILEPATH REDACTED >>>>"),
                                            cores = 11) {
  
  data_list <- fread("<<<< FILEPATH REDACTED >>>>")
  stage_list <- fread("<<<< FILEPATH REDACTED >>>>")
  stage_list <- stage_list[, .(iso3, gadm_geoid)]
  setnames(stage_list, c("iso3", "gadm_geoid"), c("country", "code"))
  data_list <- merge(data_list, stage_list, by = "country")
  
  pull_u5m_source_bias_adjustment_reg <- function(reg,
                                                  run_date,
                                                  plot,
                                                  plot_path) {
    message("Now working on reg: ", reg)
    #training data going into model
    train <- fread(sprintf("<<<< FILEPATH REDACTED >>>>"))
    #tmb object with random effect values
    tmb <- readRDS(sprintf("<<<< FILEPATH REDACTED >>>>"))
    
    # make an nid_re mapping table NID NID - Pulled from TMB code
    nidEff <- unique(select(train, country, nid)) %>% # get unique nid country combos
      group_by(country) %>% # group by country to see...
      mutate(nidCount=n()) %>% # the number of nids per country
      filter(nidCount!=1) %>% # ignore where we only have 1 nid in a country
      ungroup %>% select(nid) # simplify
    
    nidEff <- nidEff %>%
      mutate(re_id=1:n()) %>% # generate a random id number
      right_join(select(train, nid), by="nid") %>% # merge on the og data to index re
      mutate(re_id=ifelse(is.na(re_id), 0, re_id)) # replace nans with zero effect
    
    nidEff <- unique(as.data.frame(nidEff))
    
    #tmb$sdrep$par.random
    random_effects <- tmb$sdrep$par.random
    random_effects <- random_effects[names(random_effects) == "nid_re"]
    rand <- data.table("re_id" = 1:length(random_effects), "random_effect" = random_effects)
    
    #merge oto nidEff
    calc_table <- merge(nidEff, rand, all = T)
    
    #train, by nid, get weighted average haz, weight is weight
    agg_data <- train[, .(input_mean = sum(died)/sum(N)), by = c("nid", "ab")]
    #merge onto one table
    calc_table <- data.table(merge(calc_table, agg_data, by = "nid", all = T))
    
    #qlogis on average haz
    calc_table[, log_input_mean := qlogis(input_mean)]
    #that + nid level random effect
    calc_table[, log_adj_mean := log_input_mean - random_effect]
    #plogis out of logit space
    calc_table[, output_mean := plogis(log_adj_mean)]
    
    random_eff <- calc_table[, c("nid", "random_effect")]
    
    calc_table <- calc_table[, .(haz = 1-prod(1 - input_mean), haz_adj = 1-(prod(1 - output_mean))), by = "nid"]
    
    
    
    calc_table <- merge(calc_table, random_eff, by = "nid", all.x = T)
    #(transformed / original) - 1
    calc_table[, source_bias_adj := ((haz_adj / haz) - 1) * 100]
    
    #make formatted column
    calc_table[, formatted_source_bias := as.character(round(source_bias_adj, 1))]
    calc_table[!grep("-", formatted_source_bias), formatted_source_bias := paste0("+", formatted_source_bias)]
    calc_table[is.na(source_bias_adj), formatted_source_bias := "No source adjustment"]
    
    calc_table[, region := reg]
    
    #figure out valid nids within region to prevent cross-region nid duplication
    reg_codes <- get_adm0_codes(reg)
    valid_nids <- data_list[code %in% reg_codes]$nid
    calc_table <- calc_table[nid %in% valid_nids]
    
    if(plot){
      #make plot
      plot_table <- calc_table[!is.na(source_bias_adj)]
      # plot e ^ tmb$sdrep$par.fixed$log_nidre_sigma, mean 0, everything else as vert lines
      sd = exp(tmb$sdrep$par.fixed["log_nidre_sigma"]) * 100
      
      fake_data<-data.frame(value=rnorm(1000000, mean=0, sd=sd))

      source_bias <- ggplot() +
        #stat_function(fun = dnorm, geom = line, args = list(mean = 0, sd = sd)) +
        theme_classic() +
        scale_x_continuous(limits = c(-100, 100)) +
        scale_y_continuous(limits = c(0, .05), breaks = NULL) + 
        geom_vline(xintercept = 0) +
        xlab("Source bias adjustment from survey-level random effects, %") + 
        ylab(NULL) +
        ggtitle(reg)
      
      for(i in 1:nrow(plot_table)) {
        intercept <- plot_table$source_bias_adj[i]
        source_bias <- source_bias + geom_vline(xintercept = intercept)
      }
      
      source_bias <- source_bias + geom_density(data = fake_data, mapping=aes(value, colour = "darkred", fill = "darkred", alpha = .5)) +
        theme(legend.position = "none")
      
      png("<<<< FILEPATH REDACTED >>>>")
      print(source_bias)
      dev.off()
      
      saveRDS(source_bias, file = sprintf("<<<< FILEPATH REDACTED >>>>")
    }
    
    return(calc_table)
  }
  region_list <- mclapply(X=regions, FUN=pull_u5m_source_bias_adjustment_reg, run_date = run_date, plot = plot, plot_path = plot_path, mc.cores = cores)
  
  combined_table <- do.call('rbind', region_list)
  
  return(combined_table)
}


#' check for any data dropped in parallel_model
#' 
#' @description Looks for any surveys that were dropped between the start of parallel model (a.k.a) and the final modeling dataset.
#'
#' @param regions character vector of regions with model estimates/input data.
#' @param run_date model run date
#' @param cores number of threads to parallelize reading in data
#' @param countries_to_drop character vector of all caps iso3 codes. Does not compare these locations - only necessary if country is within `regions` and want to ignore for whatever reason
#' @param input_data_path path to data that feeds into parallel model, e.g. died_global_7ab_annual_weightedSS. Will be standardized as 6/2019
#'
#' @return data.table of surveys that have been dropped in parallel model (subset from the input dataset)
#'
check_for_dropped_surveys_parallel <- function(regions,
                                               run_date,
                                               cores,
                                               countries_to_drop = NULL,
                                               input_data_path = sprintf("<<<< FILEPATH REDACTED >>>>"){
  
  #data that goes into parallel_model, e.g. died_global_7ab_annual_weightedSS
  input_nids <- readRDS(input_data_path)
  
  #data that feeds directly into the model
  train <- mclapply(X=regions, FUN = function(x) {fread("<<<< FILEPATH REDACTED >>>>")}, mc.cores = cores)
  train2 <- do.call("rbind", c(train, list(fill = T)))
  
  #drop countries if specified
  if(length(countries_to_drop) > 0) {
    input_nids <- input_nids[!(country %in% countries_to_drop)]
    train2 <- train2[!(country %in% countries_to_drop)]
  }
  
  #drop all data from before 2000
  input_nids <- input_nids[year >= 2000,]
  train2 <- train2[year >= 2000,]
  
  stage_list <- fread("<<<< FILEPATH REDACTED >>>>")
  stage_list <- stage_list[, .(iso3, gadm_geoid)]
  names(stage_list) <- c("country", "code")
  
  input_nids <- merge(input_nids, stage_list, by = "country", all.x = T)
  train2 <- merge(train2, stage_list, by = "country", all.x = T)
  
  country_codes <- get_adm0_codes(regions)
  
  #drop all data from countries not in regions
  input_nids <- input_nids[code %in% country_codes,]
  train2 <- train2[code %in% country_codes,]
  
  prep_nids <- unique(train2$nid)
  missing_data <- input_nids[!(nid %in% prep_nids),]
  
  if(nrow(missing_data) > 0){
    message("The following surveys were dropped in the parallel_model script: ")
    message(paste0(unique(missing_data$nid), collapse = ", "))
  } else {
    message("No surveys were dropped in parallel_model!")
  }
  
  return(missing_data)
}

