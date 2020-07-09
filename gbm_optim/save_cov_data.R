#####################################################################
## First part of generic parallel script for running MBG models    ##
#####################################################################


## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~ SETUP ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## grab arguments from qsub
## note this requires a shell script with "<$1 --no-save $@", because its starting at 4
reg                      <- as.character(commandArgs()[4])
user                     <- as.numeric(commandArgs()[5])
core_repo                <- as.character(commandArgs()[6])
indicator_group          <- as.character(commandArgs()[7])
indicator                <- as.character(commandArgs()[8])
config_par               <- as.character(commandArgs()[9])
config_file              <- as.character(commandArgs()[10])
cov_par                  <- as.character(commandArgs()[11])
cov_file                 <- as.character(commandArgs()[12])
individual_countries     <- as.character(commandArgs()[13])
file_addin               <- as.character(commandArgs()[14])


## set other arguments
test <- 0
age <- 0
holdout <- 0

# pathaddin that gets used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

# print out session info so we have it on record
sessionInfo()

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv',header=FALSE)))
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Read config file
config <- set_up_config(repo            = core_repo,
                        indicator_group = '',
                        indicator       = '',
                        config_name     = paste0(config_file, config_par),
                        covs_name       = paste0(cov_file, cov_par))
                      
# Create a run date
run_date <- 'gbm_optim_data'

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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Prep MBG inputs/Load Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PID <- Sys.getpid()
tic("Entire script") # Start master timer

## Set seed for reproducibility
message('Setting seed 98112 for reproducibility')
set.seed(98112)


## some set up
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (class(z_list)    == "character") z_list    <- eval(parse(text=z_list))

## Load simple polygon template to model over
gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

## Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']]

## Load input data based on stratification and holdout, OR pull in data as normal and run with the whole dataset if holdout == 0.
## For holdouts, we have depreciated val level, so each val level must be recorded in a different run date
if(holdout!=0) {
  message(paste0('Holdout != 0 so loading holdout data only from holdout ',holdout))
  message('Please be sure you have a list object called stratum_ho in your environment.')
  ## if stratifies by age then make sure loads correctly
  if(age!=0) df <- as.data.table(stratum_ho[[paste('region',reg,'_age',age,sep='__')]])
  if(age==0) df <- as.data.table(stratum_ho[[paste('region',reg,sep='__')]])
  df <- df[fold != holdout, ]
}
if(holdout==0) {
  message('Holdout == 0 so loading in full dataset using load_input_data()')
  df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                        simple      = simple_polygon,
                        agebin      = age,
                        pathaddin   = pathaddin,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        use_share   = as.logical(use_share),
                        yl          = year_list)
}

## Remove any data outside of region for ORS
reg_list <- fread('<<<< FILEPATH REDACTED >>>>/dia_region_iso.csv')
reg_list[, gadm_adm0_code := get_adm0_codes(iso3, shapefile_version = modeling_shapefile_version), by = iso3]
iso3_list <- filter(reg_list, gadm_adm0_code %in% gaul_list)
iso3_list <- unique(iso3_list$iso3)
df <- filter(df, country %in% iso3_list)
df <- as.data.table(df)

## if testing, we only keep 1000 or so observations
if(test == 1){
  test_pct <- as.numeric(test_pct)
  
  message(paste0('Test option was set on and the test_pct argument was found at ',test_pct,'% \n
                 ... keeping only ', round(nrow(df)*(test_pct/100),0),' random rows of data.'))
  df <- df[sample(nrow(df),round(nrow(df)*(test_pct/100),0)),]
  
  message('Also, making it so we only take 100 draws')
  samples <- 100
}

## Some built in data checks that cause known problems later on
if(indicator_family=='binomial' & any(df[,get(indicator)]/df$N > 1))
  stop('You have binomial data where k > N. Check your data before proceeding')
if(any(df[['weight']] %in% c(Inf,-Inf) | any(is.na(df[['weight']] ))))
  stop('You have illegal weights (NA,Inf,-Inf). Check your data before proceeding')

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~ Pull Covariates ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
  
  effects <- trim(strsplit(fixed_effects, "\\+")[[1]])
  measures <- trim(strsplit(fixed_effects_measures, "\\+")[[1]])
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
all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")

## Make stacker-specific formulas where applicable
all_fixed_effects_brt <- all_fixed_effects

## Set Up Country Fixed Effects
if(use_child_country_fes == TRUE | use_inla_country_fes == TRUE) {
  message('Setting up country fixed effects')
  fe_gaul_list <- unique(c(gaul_convert(unique(df[, country]), shapefile_version = modeling_shapefile_version), gaul_list))
  fe_template  <- cov_layers[[1]][[1]]
  simple_polygon_list <- load_simple_polygon(gaul_list   = fe_gaul_list,
                                             buffer      = 0.4,
                                             subset_only = TRUE,
                                             shapefile_version = modeling_shapefile_version)
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
## ~~~~~~~~~~~~~~~~~~~~~~~~ Stacking prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
                              period_var          = 'year',
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

the_data <- merge(the_data, cs_covs[[1]], by = "a_rowid", all.x = F, all.y = F)

## store the centre scaling mapping
covs_cs_df  <-  cs_covs[[2]]

## this will drop rows with NA covariate values
the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))

## stop if this na omit demolished the whole dataset
if(nrow(the_data) == 0) stop('You have an empty df, make sure one of your covariates was not NA everywhere.')

# save the data with covariates cropped
write.csv(the_data, file = paste0('<<<< FILEPATH REDACTED >>>>', indicator, '_', reg, 
                                  ifelse(file_addin == FALSE, '', paste0('_', file_addin)), '.csv'))
