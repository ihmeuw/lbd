#####################################################################
## Generic parallel script for running MBG models                  ##
#####################################################################

sessionInfo()
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
doannual <- TRUE # if false will do a 5-year run (used for testing)
testing  <- FALSE # run on a subset of the data
testgauss <- FALSE # Test SBH Data As Gaussian
clamp_covs <- FALSE
## grab arguments
## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])

print(reg)
print(age)
print(run_date)
print(test)
print(holdout)
print(indicator)
print(indicator_group)
## make a pathaddin that get used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

## load an image of the main environment
load(paste0(<<<< FLIEPATH REDACTED >>>> '/model_image_history/pre_run_tempimage_', run_date, pathaddin,'.RData'))

## In case anything got overwritten in the load, reload args
reg                      <- as.character(commandArgs()[4])
age                      <- as.numeric(commandArgs()[5])
run_date                 <- as.character(commandArgs()[6])
test                     <- as.numeric(commandArgs()[7])
holdout                  <- as.numeric(commandArgs()[8])
indicator                <- as.character(commandArgs()[9])
indicator_group          <- as.character(commandArgs()[10])
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)
outputdir <- file.path(<<<< FLIEPATH REDACTED >>>>'output',run_date,'/')
dir.create(outputdir, showWarnings = FALSE)

## print run options
message("options for this run:\n")
for(arg in c('reg','age','run_date','test','holdout',
             'indicator','indicator_group','pathaddin','outputdir'))
  message(paste0(arg,':\t',get(arg),'\t // type: ',class(get(arg))))

# print out session info so we have it on record
sessionInfo()

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)



source(paste0(<<<< FLIEPATH REDACTED >>>> '/lbd_hiv/mbg/hiv_prev_disagg/3_functions/read_inla_prior.r'))
source(paste0(<<<< FLIEPATH REDACTED >>>> '/lbd_hiv/mbg/hiv_prev_disagg/3_functions/load_input_data_for_ap_stackers.r'))


# We need to be in the singularity image, and specifically the LBD one if using TMB
if(!is_singularity()) {
  stop('YOU MUST USE THE SINGULARITY IMAGE TO FIT YOUR MODELS.')
}

if(as.logical(fit_with_tmb) & !is_lbd_singularity()) {
  stop('YOU MUST USE THE LBD SINGULARITY IMAGE IF YOU WANT TO FIT YOUR MODEL USING TMB.')
}

## Print the core_repo hash and check it
message("Printing git hash for 'core_repo' and checking against LBD Core Code master repo")
#record_git_status(core_repo = core_repo, check_core_repo = TRUE)

## Make sure this inla patch is implemented if running on geos
if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()

## cores to use
cores_to_use <- Sys.getenv("SGE_HGR_fthread")
#if (cores_to_use=="") 
cores_to_use=1
message(paste("Model set to use", cores_to_use, "cores"))

inla.setOption("pardiso.license", <<<< FLIEPATH REDACTED >>>>)

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PID <- Sys.getpid()
tic("Entire script") # Start master timer
## some set up
if (class(year_list)     == "character") year_list    <- eval(parse(text=year_list))
if (class(z_list)        == "character") z_list       <- eval(parse(text=z_list))
if (class(sex_list)      == "character") sex_list     <- eval(parse(text=sex_list))
if (class(interacting_gp_1_effects) == "character") interacting_gp_1_effects <- eval(parse(text=interacting_gp_1_effects))
if (class(interacting_gp_2_effects) == "character") interacting_gp_2_effects <- eval(parse(text=interacting_gp_2_effects))

if(!exists('stop_before_fit'))stop_before_fit = FALSE
if(!exists('stop_after_fit'))stop_after_fit = FALSE
if(!exists('drop_low_weight')) drop_low_weight = F
if(!exists('low_weight_multiplier')) low_weight_multiplier=2
if(!exists('agesex_in_stacking')) agesex_in_stacking = F
if(!exists('remake_adult_prev_stackers')) remake_adult_prev_stackers = F
if(!exists('skipped_inla_rundate')) skipped_inla_rundate = run_date
if(!exists('use_adult_prev_stackers')){
  use_adult_prev_stackers = FALSE
} else if(use_adult_prev_stackers==TRUE){
  load(paste0('<<<< FLIEPATH REDACTED >>>>/model_image_history/', adult_prev_run_date, pathaddin, '.RData'))
  adult_prev_period_map <- period_map
  names(adult_prev_period_map) = c('data_period', 'period_id')
  stacked_rasters <- cov_list[child_model_names]
  remove(df)
  remove(mesh_s)
  remove(mesh_t)
  remove(period_map)
  remove(simple_raster)
  
}
if(exists('holdout_run_date')) holdout_outputdir = file.path(<<<< FLIEPATH REDACTED >>>>'output',holdout_run_date,'/') else holdout_outputdir = ''

message(paste0('interacting gp 1 effects: ', interacting_gp_1_effects))
message(paste0('interacting gp 2 effects: ', interacting_gp_2_effects))


###############################################################
## skip a large chunk if requested in config
if(as.logical(skiptoinla) == FALSE){
  
  message('You have chosen to not skip directly to inla.')
  
  
  ## Load simple polygon template to model over
  gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
  simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4,
                                             shapefile_version = modeling_shapefile_version)
  subset_shape        <- simple_polygon_list[[1]]
  simple_polygon      <- simple_polygon_list[[2]]
  
  ## Load list of raster inputs (pop and simple)
  raster_list        <- build_simple_raster_pop(subset_shape, pop_release=pop_release)
  simple_raster      <- raster_list[['simple_raster']]
  pop_raster         <- raster_list[['pop_raster']]
  
  ## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
  ## For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
  if(makeholdouts==T & holdout !=0){
    message("Loading holdout data")
    if(create_holdout_files==T) holdout_file = paste0(outputdir, '/holdout_file_', n_ho_folds, 'f_', reg, '.rda')
    if(create_holdout_files==F) holdout_file = paste0(holdout_outputdir, '/holdout_file_', n_ho_folds, 'f_', reg, '.rda')
    load(holdout_file)
  }
  
  if((method != "point only" & makeholdouts == T & holdout !=0) | makeholdouts == F | holdout==0) {
    message("Loading in data")
    # Load data and subset to non point
    df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                          agebin      = age,
                          removeyemen = FALSE,
                          pathaddin   = pathaddin,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list,
                          region      = reg)
    
    #Divide up data that could have been included in the holdout process and those that couldn't
    if(makeholdouts == T & holdout != 0){
      df.other.tmp <- df[df$point==0 | is.na(agebin) | type=='ANC', ]
      df.other.tmp$fold <- 0
      
      # add in location information
      df.other.tmp <- merge_with_ihme_loc(df.other.tmp, shapefile_version = modeling_shapefile_version)
      
      
      df.other.tmp <- data.frame(df.other.tmp)
      df.other.tmp <- df.other.tmp[,colnames(df.other.tmp) %in% colnames(df.points.tmp)]
      df <- rbind(df.points.tmp[,!colnames(df.points.tmp) %in% c("t_fold","ho_id")],
                  df.other.tmp)
    }
    
  } else {
    message("Loading in point data which uses holdouts")
    df <- df.points.tmp
  }
  
  df <- setDT(df)
  
  if(makeholdouts == T & holdout !=0){
    df <- df[fold!=holdout, ]
  }
  
  #---------------------------------
  #Remove non-point data if called for
  if(!exists('point_only')) point_only <- FALSE
  if(point_only==TRUE){
    df <- df[point==1 & !is.na(agebin)]
  }

  #--------------------------------------------
  ##Set up ANC variables for data weighting and ANC fixed-effect
  if(use_anc_corrections==TRUE){
    df[type=='ANC', ANC := TRUE]
    df[type!='ANC', ANC := FALSE]
    
    df[ANC == TRUE, ANC_01 := 1]
    df[ANC == FALSE, ANC_01 := 0]
    
    df[ANC == TRUE, agebin_ag:='1:8']
    
    
  } else{
    df[,ANC:=FALSE]
    df[,ANC_01:=0]
  }
  
  #---------------------------------------------
  
  ## Create df with pixel ids (and polygon centroids and age group weights if aggregation)
  if(!exists("loadinputdata_from_rundate")) {
    df <- process_input_data(df, 
                             pop_release                = pop_release, 
                             interval_mo                = interval_mo, 
                             modeling_shapefile_version = modeling_shapefile_version,
                             poly_ag                    = as.logical(poly_ag), 
                             zcol                       = zcol, 
                             zcol_ag                    = zcol_ag, 
                             zcol_ag_id                 = zcol_ag_id, 
                             z_map                      = z_map, 
                             z_ag_mat                   = z_ag_mat,
                             auto_disaggregate          = F,
                             drop_low_weight=as.logical(drop_low_weight),
                             low_weight_multiplier=as.numeric(low_weight_multiplier))
    save(df, file=paste0(outputdir,"processed_input_data", pathaddin,".rda"))
  } else {
    load(paste0(<<<< FILEPATH REDACTED >>>> '/output/',loadinputdata_from_rundate,"/","processed_input_data", pathaddin,".rda"))
    save(df, file=paste0(outputdir,"processed_input_data", pathaddin,".rda")) # save to current rundate dir
  }
  
  
  ##Calculate age/country/year-specific Fertility-Rate Ratios (FRRs)-------------------------------------------------
  if(!exists('load_data_with_frrs_from_date')){
    df$hiv_frr <- 1 #default value
    df[is.na(agebin), agebin:= 0]
    df[is.na(sex_id), sex_id:= 0]
    if(use_anc_corrections == T){
      #Split data into anc/non-anc
      if(length(z_list)>1){
        subs <- df[type == 'ANC']
        df   <- df[type != 'ANC']
        
        if(nrow(subs) > 0){
          ##Get subnat adm1 info for ZAF, ETH, KEN, NGA & IND
          
          #get adm0 codes for countries to get subnat REs
          countries_to_get_subnat_frr <- get_adm0_codes(c('zaf', 'eth', 'nga'), shapefile_version = modeling_shapefile_version)
          
          
          #load and subset standard admin1 shape to countries with subnational REs 
          subnat_full_shp    <- readRDS(get_admin_shapefile( admin_level = 1, raking = F, suffix = '.rds', version = modeling_shapefile_version ))
          subnat_shapefile <- raster::subset(subnat_full_shp, 
                                             ADM0_CODE %in% countries_to_get_subnat_frr)
          
          simple_polygon_list2 <- load_simple_polygon(gaul_list = NULL, 
                                                      buffer = 1, tolerance = 0.4, 
                                                      custom_shapefile = subnat_shapefile)
          subset_shape2        <- simple_polygon_list2[[1]]
          
          ## Load list of raster inputs (pop and simple)
          raster_list2        <- build_simple_raster_pop(subset_shape2, field = "ADM1_CODE")
          simple_raster2      <- raster_list2[['simple_raster']]
          
          #simple_raster2 has to be the same size as the simple raster for predict to work correctly
          simple_raster2 <- raster::extend(simple_raster2, extent(simple_raster))
          simple_raster2 <- raster::crop(simple_raster2, extent(simple_raster))
          
          ## Merge ADM0/1 codes to df
          adm1_subset_lox <- over(SpatialPoints(subs[,.(long = longitude, lat = latitude)], 
                                                CRS(proj4string(subnat_shapefile))), subnat_shapefile)
          subs[, subnat_frr_ADM1_CODE := as.numeric(as.character(adm1_subset_lox$ADM1_CODE))]
          subs[!(country %in%  c('ZAF', 'ETH', 'NGA')), subnat_frr_ADM1_CODE := NA]
          
          ##get locs to match lbd & gbd codes  
          locs <- get_location_code_mapping(shapefile_version=shapefile_version)
          
          for(iso in unique(subs[ANC==TRUE]$country)){
            print(iso)
            
            for(i in unique(subs[country==iso]$subnat_frr_ADM1_CODE)){
              print(i)
              if(is.na(i) & iso %in% c('ZAF', 'ETH', 'NGA')){
                next()
              } else if(is.na(i)){
                adm1 <- ''
              } else{
                adm1 <- locs[ADM_CODE==i]$loc_id
                adm1 <- paste0('_', adm1)
              }
              print(paste0('calling ', iso, adm1))
              #Read in and organize epp data for fertility etc
              epp_dat <- fread(paste0(<<<< FILEPATH REDACTED >>>>, iso, adm1, ".csv"))
              epp_dat[, age := 5*floor(age/5)] 
              epp_dat <- epp_dat[sex == "female" & age >= 15 & age < 50,] 
              epp_dat <- epp_dat[year >= min(year_list) & year <= max(year_list),]
              
              #Calculate frr's
              epp_dat <- epp_dat[, list(pop = sum(pop), 
                                        pop_neg = sum(pop_neg),  
                                        total_births = sum(total_births), 
                                        hiv_births = sum(hiv_births)),  
                                 by = 'age,year,run_num'] 
              
              epp_dat[, hiv_fert := hiv_births / (pop - pop_neg)] 
              epp_dat[, nohiv_fert := (total_births - hiv_births) / pop_neg] 
              epp_dat[, hiv_frr := hiv_fert / nohiv_fert] 
              
              epp_dat <- epp_dat[, list(hiv_frr = mean(hiv_frr)),  
                                 by = 'age,year'] 
              
              #Aggregated again by age if running an age-aggregated model
              if(zcol=="z_column_default_blank" | is.null(zcol_ag)) {
                epp_dat <- epp_dat[, list(hiv_frr = mean(hiv_frr)),  
                                   by = 'year'] 
              }
              
              #Convert epp ages to agebins & set up 50-54 FRR to equal 45-49 FRR
              epp_dat[,agebin:=as.factor(age)]
              levels(epp_dat$agebin)<-c(1:7)
              epp_dat$agebin<-as.integer(epp_dat$agebin)
              if(zcol=="z_column_default_blank" | is.null(zcol_ag)) epp_dat$agebin <- 0
              
              #Apply frr to data in particular years
              for(y in year_list) {
                for(zed in 1:7) {
                  frr <- epp_dat[year == y & agebin == zed]$hiv_frr
                  if(!is.na(i)){
                    subs[year == y & agebin == zed & country == iso & subnat_frr_ADM1_CODE== i]$hiv_frr <-frr
                  } else{
                    subs[year == y & agebin == zed & country == iso]$hiv_frr <-frr 
                  }
                if(zed==7)
                  if(!is.na(i)){
                    subs[year == y & agebin == 8 & country == iso & subnat_frr_ADM1_CODE== i]$hiv_frr <-frr
                  } else{
                    subs[year == y & agebin == 8 & country == iso]$hiv_frr <-frr 
                  }
                  }
              }
            }
          }
        }
        #Recombine data when finished
        df[,subnat_frr_ADM1_CODE:=NA]
        df <- rbind(df, subs)
        df <- df[!is.na(agg_weight)]
      } 
      
    }
    save(df, file=paste0(outputdir,"processed_input_data", pathaddin,"_with_frrs.rda"))
  } else {
    load(paste0(<<<< FLIEPATH REDACTED >>>> '/output/',load_data_with_frrs_from_date,"/processed_input_data", pathaddin,"_with_frrs.rda"))
  }
  
  
  #----------------------------------------------------------------------------------------------------
  
  if(as.logical(use_subnat_res)) {
    
    #get adm0 codes for countries to get subnat REs
    if("all" %in% subnat_country_to_get){
      countries_to_get_subnat_res <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
    } else {
      countries_to_get_subnat_res <- get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version)[get_adm0_codes(subnat_country_to_get, shapefile_version = modeling_shapefile_version) %in% get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)]
    }
    
    #load and subset standard admin1 shape to countries with subnational REs 
    subnat_full_shp    <- readRDS(get_admin_shapefile( admin_level = 1, raking = F, suffix = '.rds', version = modeling_shapefile_version ))
    subnat_shapefile <- raster::subset(subnat_full_shp, 
                                       ADM0_CODE %in% countries_to_get_subnat_res)
    
    simple_polygon_list2 <- load_simple_polygon(gaul_list = NULL, 
                                                buffer = 1, tolerance = 0.4, 
                                                custom_shapefile = subnat_shapefile)
    subset_shape2        <- simple_polygon_list2[[1]]
    
    ## Load list of raster inputs (pop and simple)
    raster_list2        <- build_simple_raster_pop(subset_shape2, field = "ADM1_CODE")
    simple_raster2      <- raster_list2[['simple_raster']]
    
    #simple_raster2 has to be the same size as the simple raster for predict to work correctly
    simple_raster2 <- raster::extend(simple_raster2, extent(simple_raster))
    simple_raster2 <- raster::crop(simple_raster2, extent(simple_raster))
    
    ## Merge ADM0/1 codes to df
    adm1_subset_lox <- over(SpatialPoints(df[,.(long = longitude, lat = latitude)], 
                                          CRS(proj4string(subnat_shapefile))), subnat_shapefile)
    df[, subnat_re_ADM1_CODE := as.numeric(as.character(adm1_subset_lox$ADM1_CODE))]
    df[, subnat_re_ADM0_CODE := as.numeric(as.character(adm1_subset_lox$ADM0_CODE))]
    
    #create new ADM1 columns for each country in subnat_country_to_get so data can be fit separately
    for(i in 1:length(unique(na.omit(df$subnat_re_ADM0_CODE)))) {
      df[subnat_re_ADM0_CODE == unique(na.omit(df$subnat_re_ADM0_CODE))[i], (paste0("SUBNAT", i)) := subnat_re_ADM1_CODE]
    }
  } else {
    simple_raster2 <- NULL
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
  
  ## if there is another weight column, multiply it with weight now
  if(exists('other_weight')) if(other_weight!='') {
    message(paste0('Multiplying weight and ',other_weight))
    df[['weight']] <- df[['weight']]*df[[other_weight]]
  }
  
  ## Some built in data checks that cause known problems later on
  if (indicator_family == 'binomial' & any(df[, get(indicator)] / df$N > 1)) {
    stop('You have binomial data where k > N. Check your data before proceeding')
  }
  
  if(any(df[['weight']] %in% c(Inf, -Inf) | any(is.na(df[['weight']])))) {
    stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding') 
  }
  
  ## for u5m in particular,  make sure indicator is properly named here, wont affect others
  df[[indicator]] <- df[[gsub(paste0('_age',age),'',indicator)]]
  
  ## Save distribution of data for this region
  png(paste0(outputdir, reg, '.png'))
  if(indicator_family=='binomial') hist(df[df$first_entry==1, get(indicator)]/df$N[df$first_entry==1]) else hist(df[df$first_entry==1, get(indicator)])
  dev.off()
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## Define modeling space. In years only for now.
  if(yearload=='annual') period_map <- make_period_map(modeling_periods = c(min(year_list):max(year_list)))
  if(yearload=='five-year') period_map <- make_period_map(modeling_periods = seq(min(year_list),max(year_list),by=5))
  
  ## Make placeholders for covariates
  cov_layers <- gbd_cov_layers <- age_rasters <- sex_rasters <- NULL
  
  ## Pull all covariate bricks/layers
  if (nrow(fixed_effects_config) > 0) {
    message('Grabbing raster covariate layers')
    loader <- MbgStandardCovariateLoader$new(start_year = min(year_list),
                                             end_year = max(year_list),
                                             interval = as.numeric(interval_mo),
                                             covariate_config = fixed_effects_config)
    cov_layers <- loader$get_covariates(simple_polygon)
  }
  
  ## Pull country level gbd covariates
  if (nchar(gbd_fixed_effects) > 0) {
    message('Grabbing GBD covariates')
    
    effects <- trim(strsplit(gbd_fixed_effects, "\\+")[[1]])
    measures <- trim(strsplit(gbd_fixed_effects_measures, "\\+")[[1]])
    gbd_cov_layers <- load_gbd_covariates(covs     = effects,
                                          measures = measures,
                                          year_ids = year_list,
                                          age_ids  = gbd_fixed_effects_age,
                                          template = cov_layers[[1]][[1]],
                                          simple_polygon = simple_polygon,
                                          interval_mo = interval_mo)
  }
  
  ###  
  
  if(agesex_in_stacking==T){
    age_raster <- cov_layers[[1]][[1]]
    age_rasters <- age_raster
    names(age_rasters) = 'age.1'
    for(i in 2:length(z_list)){
      age_rasters <-stack(age_rasters, age_raster)
      names(age_rasters)[i] = paste0('age.', i)
      age_rasters[[i]] = i
    }
    age_rasters[[1]] = 1
    
    sex_rasters <-stack(age_raster, age_raster)
    sex_rasters[[1]] = 1
    sex_rasters[[2]] = 2
    names(sex_rasters) = c('sex.1', 'sex.2')
    
    agesex_layers <- c(age_rasters, sex_rasters)
    names(agesex_layers) <- c('age_cov', 'sex_cov')
  }
  
  
  
  ## Combine all covariates
  all_cov_layers <- c(cov_layers, gbd_cov_layers) #Add in agesex layers later
  
  ## regenerate all fixed effects equation from the cov layers
  all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")
  
  ## Make stacker-specific formulas where applicable
  all_fixed_effects_brt <- all_fixed_effects
  
  ## Set Up Country Fixed Effects
  if(use_child_country_fes == TRUE | use_inla_country_fes == TRUE) {
    message('Setting up country fixed effects')
    fe_gaul_list <- unique(c(gaul_convert(unique(df[, country]),
                                          shapefile_version =
                                            modeling_shapefile_version),
                             gaul_list))
    fe_template  <- cov_layers[[1]][[1]]
    simple_polygon_list <- load_simple_polygon(gaul_list   = fe_gaul_list,
                                               buffer      = 0.4,
                                               subset_only = TRUE,
                                               shapefile_version = modeling_shapefile_version)
    fe_subset_shape     <- simple_polygon_list[[1]]
    gaul_code <- rasterize_check_coverage(fe_subset_shape, fe_template, field = 'ADM0_CODE')
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
  
  df[, period := year - 1999]
  df[, period_id := year - 1999]
  
  ## copy the dataset to avoid unintended namespace conflicts
  new_data <- copy(df)
  
  if(!as.logical(use_adult_prev_stackers)){ #If you don't plan to load stackers from another source
    
    ## only use data where we know what age group or point
    if(remake_adult_prev_stackers==F | z_list == 0){
      ag_data <- new_data[new_data$agg_weight!=1 | ANC_01==1, ]
      the_data <- new_data[new_data$agg_weight==1 & ANC_01==0, ]
      remove(new_data)
    }  else{
      
      the_data <- load_input_data_for_ap_stackers(indicator   = gsub(paste0('_age',age),'',indicator),
                                                  agebin      = age,
                                                  removeyemen = FALSE,
                                                  pathaddin   = pathaddin,
                                                  years       = yearload,
                                                  withtag     = as.logical(withtag),
                                                  datatag     = '_age_aggregated_pointonly_survey_only',
                                                  use_share   = as.logical(use_share),
                                                  yl          = year_list,
                                                  region      = reg)
      ag_data <- copy(new_data)
      remove(new_data)
    }
    
    #remove data that's being held out in the true data
    the_data <- match_df(the_data, df[point==1], on=c('latitude','longitude','year'))
    
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
                                  period_var          = 'year',
                                  period_map          = period_map)
    
    # A check to see if any of the variables do not vary across the data. This could break model later so we check and update some objects
    covchecklist <- check_for_cov_issues(check_pixelcount = check_cov_pixelcount,
                                         check_pixelcount_thresh = ifelse(exists("pixelcount_thresh"), as.numeric(pixelcount_thresh), 0.95))
    for(n in names(covchecklist)){
      assign(n, covchecklist[[n]])
    }
    
    # plot covariates as a simple diagnostic here
    pdf(sprintf('%s/raw_covariates_%s.pdf',outputdir,pathaddin), height=12, width=12)
    for(covname in names(all_cov_layers)){
      plot(all_cov_layers[[covname]],main=covname,maxpixel=1e6)
    }
    dev.off()
    
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
    
    the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)
    
    
    ######################################################
    #centre_scale the age-sex variables
    if(agesex_in_stacking==T){  
      the_data[, age_cov := agebin/max(agebin)]
      the_data[, sex_cov := sex_id]
      
      design_matrix = data.frame(the_data[,names(agesex_layers), with =F])
      cs_the_data <- getCentreScale(design_matrix)
      design_matrix <- centreScale(design_matrix, df = cs_the_data)
      
      #replace the the_data columns with the design matrix
      the_data[, names(agesex_layers) := NULL]
      the_data = cbind(the_data, design_matrix)
    }
    #######################################################
    
    
    ## store the centre scaling mapping
    covs_cs_df  <-  cs_covs[[2]]
    
    #Add age & sex to various effects mapping objects
    if(agesex_in_stacking==T){  
      mz = mean(the_data$age_cov)
      sz = sd(the_data$age_cov)
      
      z <- as.data.table('age_cov')
      z <- cbind(z, mz, sz)
      
      ms = mean(the_data$sex_cov)
      ss = sd(the_data$sex_cov)
      
      s <- as.data.table('sex_cov')
      s <- cbind(s, ms, ss)
      
      q <- rbind(z, s, use.names=F)
      
      covs_cs_df  <- rbind(covs_cs_df, q, use.names = F)
      
      the_covs <- c(the_covs, 'age_cov', 'sex_cov')
      all_fixed_effects = paste0(all_fixed_effects, ' + age_cov + sex_cov')
      all_fixed_effects_brt = all_fixed_effects
    }
    
    if(as.logical(use_raw_covs) | as.logical(use_stacking_covs)) {
      ## this will drop rows with NA covariate values
      the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))
    }
    
    ## stop if this na omit demolished the whole dataset
    if(nrow(the_data) == 0) stop('You have an empty df, make sure one of your covariates was not NA everywhere.')
    
    ## seperated out into a different script
    if(as.logical(use_stacking_covs)){
      message('Fitting Stackers')
      
      # Run the child stacker models 
      child_model_run <- run_child_stackers(models = child_model_names, input_data = the_data)
      
      # Bind the list of predictions into a data frame
      child_mods_df <- do.call(cbind, lapply(child_model_run, function(x) x[[1]]))
      
      ## combine the children models with the_data
      the_data  <- cbind(the_data, child_mods_df)
      
      ## Rename the child model objects into a named list
      child_model_objs <- setNames(lapply(child_model_run, function(x) x[[2]]), child_model_names)
      
      
      
      ## return the stacked rasters
      ###ADDING OPTION TO PRODUCE AGE/SEX-SPECIFIC STACKERS###
      if(agesex_in_stacking == TRUE){
        stacked_rasters <- list()
        stacked_rasters[[1]] <- list()
        stacked_rasters[[2]] <- list()
        for(x in 1:length(sex_list)){
          for(z in 1:length(z_list)){
            stacked_rasters[[x]][[z]] <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                                            period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                                            child_models     = child_model_objs,
                                                            indicator_family = indicator_family,
                                                            centre_scale_df  = covs_cs_df,
                                                            age=z,
                                                            sex=x)
            for(i in 1:length(stacked_rasters[[x]][[z]])){
              names(stacked_rasters[[x]][[z]][[i]]) = paste(names(stacked_rasters[[x]][[z]][[i]]), z, x, sep = '_')
            }
          }
        }
      } else{
        stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                              period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                              child_models     = child_model_objs,
                                              indicator_family = indicator_family,
                                              centre_scale_df  = covs_cs_df)
      }
      
      
      ## plot stackers
      message('plot stackers')
      pdf(paste0(outputdir, 'stacker_rasters', pathaddin, '.pdf'))
      if(agesex_in_stacking == TRUE){
        for(x in 1:length(sex_list)){
          for(z in 1:length(z_list)){
            for(i in 1:length(stacked_rasters[[x]][[z]])){
              plot(stacked_rasters[[x]][[z]][[i]],main=names(stacked_rasters[[x]][[z]][[i]]),maxpixel=ncell(stacked_rasters[[x]][[z]][[i]]))
            }
          }
        }
      } else {
        for(i in 1:length(stacked_rasters))
          plot(stacked_rasters[[i]],main=names(stacked_rasters[[i]]),maxpixel=ncell(stacked_rasters[[i]]))
      }
      dev.off()
      
      message('Stacking is complete')
    } ## if(use_stacking_covs)
    
    ## add aggregate data back in, with stacking predictions from the full model
    if (nrow(ag_data) > 0) {
      print(paste0('nrow ag data is ', nrow(ag_data)))
      ag_data[, a_rowid := 1:.N + max(the_data$a_rowid)]
      if(as.logical(use_stacking_covs)) {
        if(agesex_in_stacking == T){
          all_ag_data = data.table()
          for(x in sort(unique(ag_data$sex_id))){
            print(paste0('working on sex ', x))
            for(z in sort(unique(ag_data$agebin))){
              print(paste0('working on age ', z))
              as_ag_data  =  ag_data[sex_id==x & agebin ==z]
              as_ag_data  =  as_ag_data[!is.na(latitude) & !is.na(longitude)]
              if(nrow(as_ag_data) ==0) next(paste0('No ag data for sex = ', x,' age = ', z, ', moving to next'))
              ag_stackers <- extract_covariates(as_ag_data,
                                                stacked_rasters[[x]][[z]],
                                                id_col              = "a_rowid",
                                                return_only_results = TRUE,
                                                centre_scale        = FALSE,
                                                period_var          = "year",
                                                period_map          = period_map)
              
              ag_stackers <- ag_stackers[, c("a_rowid", child_model_names, child_model_names), with = F]
              
              stacker_names <- c(paste0(child_model_names, "_full_pred"), paste0(child_model_names, "_cv_pred"))
              setnames(ag_stackers, c("a_rowid", stacker_names))
              
              as_ag_data <- merge(as_ag_data, ag_stackers, by = "a_rowid", all.x=T)
              
              all_ag_data <- rbind(all_ag_data, as_ag_data)
              
              all_ag_data[is.na(period_id), period_id := year - 1999]
            }
          }
          
          ag_data <- all_ag_data
          remove(all_ag_data)
          
          print(paste0('nrow ag data is ', nrow(ag_data)))
          
          if(any(is.na(ag_data[, ..stacker_names]))) {
            warnings("There are NAs in predictions from stackers for aggregated data. 
               Please contact the core code team if you encounter this problem.")
            ag_data <-ag_data[!is.na(gam_full_pred)]
          }
          
        }else{
          ag_stackers <- extract_covariates(ag_data,
                                            stacked_rasters,
                                            id_col              = "a_rowid",
                                            return_only_results = TRUE,
                                            centre_scale        = FALSE,
                                            period_var          = "year",
                                            period_map          = period_map)
          
          ag_stackers <- ag_stackers[, c("a_rowid", child_model_names, child_model_names), with = F]
          
          stacker_names <- c(paste0(child_model_names, "_full_pred"), paste0(child_model_names, "_cv_pred"))
          setnames(ag_stackers, c("a_rowid", stacker_names))
          
          ag_data <- merge(ag_data, ag_stackers, by = "a_rowid")
        }
        
        
        
        if(any(is.na(ag_data[, ..stacker_names]))) {
          warnings("There are NAs in predictions from stackers for aggregated data. 
               Please contact the core code team if you encounter this problem.")
        }
      } else {
        ag_covs <- extract_covariates(ag_data,
                                      all_cov_layers,
                                      id_col              = "a_rowid",
                                      return_only_results = TRUE,
                                      centre_scale        = FALSE,
                                      period_var          = 'year',
                                      period_map          = period_map)
        cov_names <- names(all_cov_layers)
        
        cs_ag_covs <- centreScale(ag_covs[, ..cov_names], df = covs_cs_df)
        ag_covs <- cbind(ag_covs[,- ..cov_names], cs_ag_covs)
        
        ag_data <- merge(ag_data, ag_covs, by = "a_rowid")
        
        
        if(as.logical(use_raw_covs)) {
          if(any(is.na(ag_data[, ..cov_names]))) {
            stop("There are NAs in covariates for aggregated data. 
               Please contact the core code team if you encounter this problem.")
          }
        }
      }
      
      if(remake_adult_prev_stackers==T & z_list != 0){
        the_data <- data.table()
      }
      
      the_data <- rbind(the_data, ag_data, fill = T)
      
      #Check if any data points got missed in stacking extraction
      if(use_stacking_covs==T){
        if(!any(the_data[is.na(gam_full_pred)]$point==0)){
          message(paste0(nrow(the_data[is.na(gam_full_pred)]), ' points being removed as they fall outside simple raster'))
          the_data <- the_data[!is.na(gam_full_pred)]
        } else if(any(the_data[is.na(gam_full_pred)]$point==0)){
          warning('there is polygon data where stackers were not extracted, investigate')
          the_data <- the_data[!is.na(gam_full_pred)]
        }
      }
    }
  } else { #Load in adult prev stackers and extract those  
    ## add a row id column
    the_data[, a_rowid := seq(1:nrow(the_data))]
    the_data[, period := year - 1999]
    the_data[, period_id := year - 1999]
    the_stackers <- extract_covariates(the_data,
                                       cov_list,
                                       id_col              = "a_rowid",
                                       return_only_results = TRUE,
                                       centre_scale        = FALSE,
                                       period_var          = 'year',
                                       period_map          = adult_prev_period_map)
    
    the_stackers <- the_stackers[, c("a_rowid", child_model_names, child_model_names), with = F]
    
    stacker_names <- c(paste0(child_model_names, "_full_pred"), paste0(child_model_names, "_cv_pred"))
    setnames(the_stackers, c("a_rowid", stacker_names))
    
    the_data <- merge(the_data, the_stackers, by = "a_rowid")
  }
  
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~ Final Pre-MBG Processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  ## set the fixed effects to use in INLA based on config args
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
  if(as.logical(use_stacking_covs) & as.logical(use_gbd_covariate)){
    all_fixed_effects <- paste(stacked_fixed_effects, 'gbd_est', sep = " + ")
  }
  
  ## copy things back over to df
  df <- copy(the_data)
  ## remove the covariate columns so that there are no name conflicts when they get added back in
  df <- df[,paste0(the_covs) := rep(NULL, length(the_covs))]
  
  ## Double-check that gaul codes get dropped before extracting in save_mbg_input()
  df <- df[, grep('gaul_code_*', names(df), value = T) := rep(NULL, length(grep('gaul_code_*', names(df), value = T)))]
  
  ## create a full raster list to carry though to the shiny/next steps
  if(as.logical(use_stacking_covs) & !as.logical(use_adult_prev_stackers)){
    cov_list    <- c(unlist(stacked_rasters),unlist(all_cov_layers))
    child_mod_ras <- unlist(stacked_rasters) 
  }else if(!as.logical(use_stacking_covs)){
    cov_list <- unlist(all_cov_layers)
    child_model_names <- ''
  }
  
  toc(log = T) ## End stacking master timer
  
  ## make sure this inla patch is implemented if running on geos
  if(grepl("geos", Sys.info()[4])) INLA:::inla.dynload.workaround()
  
  ## Build spatial mesh over modeling area--for interacting space terms
  mesh_int <- build_space_mesh(d           = df,
                               simple      = simple_polygon,
                               max_edge    = mesh_s_max_edge,
                               mesh_offset = mesh_s_offset,
                               s2mesh = as.logical(use_s2_mesh),
                               s2params = s2_mesh_params_int)
  
  ##build spatial mesh for space-only GP
  mesh_s   <- build_space_mesh(d           = df,
                               simple      = simple_polygon,
                               max_edge    = mesh_s_max_edge,
                               mesh_offset = mesh_s_offset,
                               s2mesh = as.logical(use_s2_mesh),
                               s2params = s2_mesh_params_s)
  
  
  ## Build temporal mesh (standard for now)
  if (length(unique(year_list)) == 1) {
    mesh_t <- NULL
  } else {
    mesh_t <- build_time_mesh(periods = eval(parse(text = mesh_t_knots)))
  }
  
  
  ## ## For raw covs, don't want to center-scale (as that will happen in `build_mbg_data_stack()`)
  ##
  ## ## When stacking covs are used the oos-stackers (used in `fit_mbg()`)
  ## ## do not get center scaled in save_mbg_input() - this just harmonizes the measures.  If this
  ## ## step isn't done, then the covs get double-center-scaled with odd results.
  ##
  ## ## For predict_mbg, non-center-scaled covs are pulled from cov_list (either stackes or raw) and
  ## ## center-scaled within the function.  So both fit and predict take place on center-scaled covs
  ##
  ## ## But for now, we just do this
  if (as.logical(use_raw_covs) == TRUE) {
    centre_scale_covs <- FALSE
  } else {
    centre_scale_covs <- TRUE
  }
  
  ## Save all inputs for MBG model into correct location on /share
  save_mbg_input(indicator         = indicator,
                 indicator_group   = indicator_group,
                 df                = df,
                 simple_raster     = simple_raster,
                 mesh_s            = mesh_s,
                 mesh_t            = NULL,
                 mesh_int          = mesh_int,
                 cov_list          = cov_list,
                 pathaddin         = pathaddin,
                 run_date          = run_date,
                 child_model_names = child_model_names,
                 all_fixed_effects = all_fixed_effects,
                 period_map        = period_map,
                 centre_scale      = centre_scale_covs,
                 reformat_df_covs  = F)
  
} else { ## END !SKIPTOINLA
  message(paste0('You have chosen to skip directly to INLA. Picking up objects from run_date ',skiptoinla_from_rundate))
  message('Now copying saved MBG inputs from that chosen run_date.')
  
  file.copy(from = paste0(<<<< FLIEPATH REDACTED >>>> '/model_image_history/', skiptoinla_from_rundate, pathaddin, '.RData'),
            to = paste0(<<<< FLIEPATH REDACTED >>>> '/model_image_history/', run_date, pathaddin, '.RData'))
  
  file.copy(from = paste0(<<<< FLIEPATH REDACTED >>>> '/output/', skiptoinla_from_rundate, '/input_data', pathaddin, '.csv'),
            to = paste0(<<<< FLIEPATH REDACTED >>>> '/output/', run_date, '/input_data', pathaddin, '.csv'))
  
  file.copy(from = paste0(<<<< FLIEPATH REDACTED >>>> '/output/', skiptoinla_from_rundate, 
                          '/processed_input_data_', reg, '_holdout_', holdout, '_agebin_', age, '.rda'),
            to = paste0(<<<< FLIEPATH REDACTED >>>>'/output/', run_date,
                        '/processed_input_data_', reg, '_holdout_', holdout, '_agebin_', age, '.rda'))
}

## reload data an prepare for MBG
load(paste0(<<<< FLIEPATH REDACTED >>>>'/model_image_history/', run_date, pathaddin, '.RData'))

# Bound GBM to 0-1 if desired
if (exists("gbm_bounded_0_1")) {
  if (as.logical(gbm_bounded_0_1) == T) {
    message("Truncating GBM values > 1 to 0.999")
    values(cov_list[["gbm"]])[values(cov_list[["gbm"]]) >= 1 & !is.na(values(cov_list[["gbm"]]))] <- 0.999
    gbm_cols <- grep(paste0("(gbm)(.*_pred)"), names(df), value=T)
    replace_one <- function(x) {
      x[x>=1 & !is.na(x)] <- 0.999
      return(x)
    }
    df[, (gbm_cols) := lapply(.SD, replace_one), .SDcols = gbm_cols]
  }
}

## convert stackers to transform space, if desired
## NOTE: we do this here to ensure that the stacker rasters are saved in prevalence/untransformed space
## this is useful for diagnostics and other code that was built expecting the untransformed rasters
if (as.logical(stackers_in_transform_space) & indicator_family == 'binomial' & as.logical(use_stacking_covs)){
  message('Converting stackers to logit space')
  
  if(as.logical(agesex_in_stacking)) stacked_rasters <- cov_list[1:(length(z_list)*length(child_model_names)*length(sex_list))] else stacked_rasters <- cov_list
  
  ## transform the rasters
  for (ii in 1:length(stacked_rasters)) {      
    
    ## Preserve variable names in the raster first
    tmp_rastvar <- names(stacked_rasters[[ii]])
    
    ## Logit
    stacked_rasters[[ii]] <- logit(stacked_rasters[[ii]])
    
    ## Reassign names
    names(stacked_rasters[[ii]]) <- tmp_rastvar
    rm(tmp_rastvar)
  }
  
  cov_list <- stacked_rasters
  remove(stacked_rasters)
  
  ## transform the stacker values that are in df
  stacker_cols <- grep(paste0("(", paste(child_model_names, collapse="|"), ")(.*_pred)"), names(df), value=T)
  
  if(use_gbd_covariate) stacker_cols <- c(stacker_cols, 'gbd_est')
  df[, (stacker_cols) := lapply(.SD, logit), .SDcols = stacker_cols]
  
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Run MBG ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


tic("MBG - all") ## Start MBG master timer

## for stacking, overwrite the columns matching the model_names so that we can trick inla into being our stacker
if(as.logical(use_stacking_covs)){
  df[, paste0(child_model_names) := lapply(child_model_names, function(x) get(paste0(x,'_cv_pred')))]
}

## Generate MBG formula for INLA call (will run but not used by TMB)
if(!as.logical(fit_with_tmb)){
  mbg_formula <- build_mbg_formula_with_priors(fixed_effects               = all_fixed_effects,
                                               add_nugget                  = use_inla_nugget,
                                               nugget_prior                = nugget_prior,
                                               add_ctry_res                = use_country_res,
                                               ctry_re_prior               = ctry_re_prior,
                                               temporal_model_theta1_prior = rho_prior,
                                               no_gp                       = !as.logical(use_gp),
                                               use_space_only_gp           = as.logical(use_space_only_gp),
                                               use_time_only_gmrf          = as.logical(use_time_only_gmrf),
                                               time_only_gmrf_type         = time_only_gmrf_type,
                                               coefs.sum1                  = coefs_sum1,
                                               subnat_RE                   = use_subnat_res,
                                               subnat_country_to_get       = subnat_country_to_get,
                                               subnat_re_prior             = subnat_re_prior,
                                               adm0_list                    = gaul_list)
  
  ## Add ANC offset
  if (anc_correction == 'iid') {
    df[, anc := as.numeric(type == "ANC")]
    df[, site := as.character(site)]
    df[!is.na(site), site := paste(country, site)]
    df[, site := factor(site)]
    mbg_formula <- as.formula(paste("covered ~", mbg_formula[3], '+ anc + f(site, model = "iid")'))
    
  } else if (anc_correction == 'gp') {
    df[, anc := as.numeric(type == "ANC")]
    mbg_formula <- as.formula(paste("covered ~", mbg_formula[3], '+ anc + f(site, model = spde)'))
    
  }
}
## For INLA we need to add data for missing time points to ensure we get predictions
##  for all relevant time points. The 0 observations do not contribute to the 
##  model fitting but they prevent INLA from auto-removing 
##  random effects that (conditionally) have no data impacting their fit
if (!as.logical(fit_with_tmb)) {
  if(use_timebyctry_res) {
    ## If we are using a time only effect by country then we need to make sure 
    ##  all year effects are estimated for each country.
    df$adm0code <- gaul_convert(df$country)
    for(adm0_code in gaul_list) {
      dfsub <- df[df$adm0code == adm0_code, ]
      missing_years <- setdiff(year_list, dfsub$year)
      if (length(missing_years) > 0) {
        fake_data <- dfsub[1:length(missing_years), ]
        fake_data[, year := missing_years]
        fake_data[, c(indicator, 'N', 'weight') := 0]
        fake_data[, period := NULL]
        fake_data <- merge(fake_data, period_map)
        df <- rbind(df, fake_data)
      }
    }
  } else {
    ## If not, we only need to make sure we have an observation for each missing
    ##  year (country irrelevant)
    missing_years <- setdiff(year_list, df$year)
    if (length(missing_years) > 0) {
      fake_data <- df[1:length(missing_years), ]
      fake_data[, year := missing_years]
      fake_data[, c(indicator, 'N', 'weight') := 0]
      fake_data[, period := NULL]
      fake_data <- merge(fake_data, period_map)
      df <- rbind(df, fake_data)
    }
  }
}

## Dont rescale year as covariate
if ("year_cov" %in% all_fixed_effects) {
  exclude_cs <- c(child_model_names, 'year_cov')
} else {
  exclude_cs <- child_model_names
}

############ SET UP SBH TO BE MODELLED AS GAUSSIAN
if(testgauss==TRUE){
  df[,sd_i:=0]
  df[data_type == 'sbh', sd_i := sqrt(var_logit_lss_w)] 
  df[data_type == 'sbh', died := qlogis(haz)]
  df[,lik_fam_binom := ifelse(data_type=='sbh',0,1)]
  df[,lik_fam_gauss := ifelse(data_type=='sbh',1,0)]
  scale_gaussian_variance_N<- FALSE
}

if(use_gbd_covariate) child_model_names <- c(child_model_names, 'gbd_est')

## Create SPDE INLA stack
input_data <- build_mbg_data_stack(df            = df,
                                   fixed_effects = all_fixed_effects,
                                   mesh_s        = mesh_s,
                                   mesh_int = mesh_int,
                                   mesh_t        = mesh_t,
                                   spde_prior    = spde_prior,
                                   use_ctry_res  = use_country_res,
                                   use_subnat_res  = use_subnat_res, 
                                   use_gp = as.logical(use_gp),
                                   st_gp_int_zero = as.logical(st_gp_int_zero), 
                                   use_space_only_gp = as.logical(use_space_only_gp),
                                   s_gp_int_zero = as.logical(s_gp_int_zero),
                                   use_time_only_gmrf = as.logical(use_time_only_gmrf),
                                   use_age_only_gmrf = as.logical(use_age_only_gmrf),
                                   use_sex_only_gmrf = as.logical(use_sex_only_gmrf),
                                   use_nugget    = use_inla_nugget,
                                   exclude_cs    = child_model_names, # raw covs will get center scaled here though (see notes above)
                                   coefs.sum1    = coefs_sum1, 
                                   tmb           = fit_with_tmb,
                                   scale_gaussian_variance_N = scale_gaussian_variance_N,
                                   shapefile_version = modeling_shapefile_version, 
                                   zl            = z_list, 
                                   zcol          = zcol,
                                   use_sz_gp = as.logical(use_sz_gp),
                                   use_sx_gp = as.logical(use_sx_gp),
                                   use_tx_gp = as.logical(use_tx_gp),
                                   use_zx_gp = as.logical(use_zx_gp),
                                   use_cre_z_gp = as.logical(use_cre_z_gp),
                                   use_cre_x_gp = as.logical(use_cre_x_gp),
                                   use_anc_corrections = as.logical(use_anc_corrections),
                                   use_error_iid_re = as.logical(use_error_iid_re),
                                   use_observation_level_error_iid = as.logical(use_observation_level_error_iid),
                                   use_cyzx_error_iid = as.logical(use_cyzx_error_iid),
                                   xcol = sex_col,
                                   xl = sex_list)

## combine all the inputs, other than cs_df these are not used if you are using TMB
stacked_input  <- input_data[[1]]
spde_int           <- input_data[[2]] ## used for space-time gp
cs_df          <- input_data[[3]]
#spde.sp        <- input_data[[4]] ## used for space only (time stationary) gp

## Generate other inputs necessary
outcome <- df[[indicator]] # N+_i - event obs in cluster
N       <- df$N                  # N_i - total obs in cluster
weights <- df$weight

## catch in case there is no weight column
if(is.null(weights)){
  weights = rep(1,nrow(df))
}

tic("MBG - fit model") ## Start MBG - model fit timer

## Set the number of cores to be equal to input;
## If missing, then revert to cores_to_use value
if(Sys.getenv("OMP_NUM_THREADS") != "") {
  setompthreads(Sys.getenv("OMP_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input OpenMP threads.")
  setompthreads(cores_to_use)
}

if(Sys.getenv("MKL_NUM_THREADS") != "") {
  setmklthreads(Sys.getenv("MKL_NUM_THREADS"))
} else {
  print("Threading information not found; setting cores_to_use as the input MKL threads.")
  setmklthreads(cores_to_use)
}

# save RDS file of input data for replication
saveRDS(object = input_data, ## save this here in case predict dies
        file = sprintf('<<<< FLIEPATH REDACTED >>>>/output/%s/%s_data_input_list_%s_holdout_%i_agebin_%i.RDS',
                     run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))

input_data<-readRDS( ## save this here in case predict dies
  file = sprintf('<<<< FLIEPATH REDACTED >>>>/output/%s/%s_data_input_list_%s_holdout_%i_agebin_%i.RDS',
                 run_date, ifelse(fit_with_tmb,'tmb','inla'), reg, holdout, age))
#########################################################
##If specified in config, launch fitting in a separate job
    memory <- c(cssa = 200, 'essa_sdn-COM' = 550, sssa = 600, 'wssa-CPV-STP-MRT' = 850)[reg]
  
    project = 'proj_geo_nodes'
    que   = 'geospatial.q'
  
  qsub <- make_qsub_share(code_path      = sprintf("%s/mbg_central/share_scripts/parallel_model_fitting.R", core_repo),
                          addl_job_name  ='fit',
                          age            = age,
                          reg            = reg,
                          holdout        = holdout,
                          test           = testing,
                          indic          = indicator,
                          saveimage      = FALSE,
                          memory         = memory,
                          cores          = 5,
                          proj           = project,
                          geo_nodes      = FALSE,
                          use_c2_nodes   = FALSE,
                          singularity    = '<<<< FILEPATH REDACTED >>>>/lbd_full_20200128.simg',
                          singularity_opts = list(SET_OMP_THREADS=1, SET_MKL_THREADS=1),
                          queue          = que,
                          run_time       = '05:00:00:00')
  
  system(qsub)


message('Reached model fitting, which has been launched. Done with this job!')
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ FIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

