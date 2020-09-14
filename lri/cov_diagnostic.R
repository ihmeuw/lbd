# ---------------------------------------------------------------------------------------------------
# Purpose: Plot covariates by pixel in countries or admins with weird stacker patterns
# 
# 1. Load covariate and stacker rasters (in the future could save these in parallel model)
# 2. Extract pixel in country or admin that's causing problems
# 3. Plot a pixel  to investigate which covariates follow these weird patterns
# 4. Save results
# ---------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------
# 1. Setup

#pull arguments from qsub call
run_date <- commandArgs()[4]
reg <- commandArgs()[5]
country <- commandArgs()[6]
lat <- as.numeric(commandArgs()[7])
long <- as.numeric(commandArgs()[8])
stacker_names <- unlist(strsplit(commandArgs()[9], '~'))
test <- commandArgs()[10]
holdout <- as.numeric(commandArgs()[11])
age <- as.numeric(commandArgs()[12])
pixel_id <- as.numeric(commandArgs()[13])

stacker_names

#set other arguments
indicator                <- 'has_lri'
indicator_group          <- 'lri'

# make a pathaddin that gets used widely
pathaddin <- paste0('_bin',age,'_',reg,'_',holdout)

# set output directory
outputdir <- '<<<< FILEPATH REDACTED >>>>'

# load an image of the main environment
load('<<<< FILEPATH REDACTED >>>>')

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

#source any functions with custom edits from lbd_custom folder
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/stacking_functions.R')

# Throw a check for things that are going to be needed later
message('Looking for things in the config that will be needed for this script to run properly')
check_config()

# Make sure year list is correct format
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (class(z_list)    == "character") z_list    <- eval(parse(text=z_list))
# ---------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------
# 2. Load covariate and stacker rasters (in the future could save these in parallel model)

# set seed for reproducibility
set.seed(98112)

# Load simple polygon template to model over
gaul_list           <- get_adm0_codes(reg, shapefile_version = modeling_shapefile_version)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

# Load list of raster inputs (pop and simple)
raster_list        <- build_simple_raster_pop(subset_shape)
simple_raster      <- raster_list[['simple_raster']]
pop_raster         <- raster_list[['pop_raster']]

# Define modeling space. In years only for now.
if(yearload=='annual') period_map <-
  make_period_map(modeling_periods = c(min(year_list):max(year_list)))
if(yearload=='five-year') period_map <-
  make_period_map(modeling_periods = seq(min(year_list),max(year_list),by=5))

# load the data
df <- load_input_data(indicator   = gsub(paste0('_age',age),'',indicator),
                      simple      = simple_polygon,
                      agebin      = age,
                      removeyemen = TRUE,
                      pathaddin   = pathaddin,
                      years       = yearload,
                      withtag     = as.logical(withtag),
                      datatag     = datatag,
                      use_share   = as.logical(use_share),
                      yl          = year_list)

## Remove any data outside of region for ORS
if (indicator != 'had_diarrhea') {
  reg_list <- fread('<<<< FILEPATH REDACTED >>>>')
  reg_list[, gadm_adm0_code := get_adm0_codes(iso3, shapefile_version = modeling_shapefile_version), by = V1]
  iso3_list <- filter(reg_list, gadm_adm0_code %in% gaul_list)
  iso3_list <- unique(iso3_list$iso3)
  df <- filter(df, country %in% iso3_list)
  df <- as.data.table(df)
}

# Make placeholders for covariates
cov_layers <- gbd_cov_layers <- NULL

# Pull all covariate bricks/layers
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

# Pull country level gbd covariates
gbd_fixed_effects <- ''
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

# Combine all covariates
all_cov_layers <- c(cov_layers, gbd_cov_layers)

# regenerate all fixed effects equation from the cov layers
all_fixed_effects <- paste(names(all_cov_layers), collapse = " + ")

## Make stacker-specific formulas where applicable
all_fixed_effects_brt <- all_fixed_effects

# Format covariates
the_covs <- format_covariates(all_fixed_effects)

# copy data
the_data <- copy(df)

# add a row id column
the_data[, a_rowid := seq(1:nrow(the_data))]

# Figure out which models we're going to use
child_model_names <- stacker_names

the_covs <- format_covariates(all_fixed_effects)

# copy the dataset to avoid unintended namespace conflicts
the_data <- copy(df)

# shuffle the data into six folds
the_data <- the_data[sample(nrow(the_data)),]
the_data[,fold_id := cut(seq(1,nrow(the_data)),breaks=as.numeric(n_stack_folds),labels=FALSE)]

# add a row id column
the_data[, a_rowid := seq(1:nrow(the_data))]

# extract covariates to the points and subset data where its missing covariate values
cs_covs <- extract_covariates(the_data,
                              all_cov_layers,
                              id_col              = "a_rowid",
                              return_only_results = TRUE,
                              centre_scale        = TRUE,
                              period_var          = 'year',
                              period_map          = period_map)

# Check for data where covariate extraction failed
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

# store the centre scaling mapping
covs_cs_df  <-  cs_covs[[2]]

# this will drop rows with NA covariate values
the_data    <- na.omit(the_data, c(indicator, 'N', the_covs))

# stop if this na omit demolished the whole dataset
if(nrow(the_data) == 0) stop('You have an empty df, make sure one of your covariates was not NA everywhere.')

# sourcing this script will run the child stackers:
source(paste0('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/run_child_stackers.R'))

# combine the children models
the_data  <- cbind(the_data, do.call(cbind, lapply(lapply(child_model_names, 'get'), function(x) x[[1]])))
child_model_objs <- setNames(lapply(lapply(child_model_names, 'get'), function(x) x[[2]]), child_model_names)

# fit GAM stacker
stacked_results <- gam_stacker(the_data,
                               model_names      = child_model_names,
                               indicator        = indicator,
                               indicator_family = indicator_family)

# return the stacked rasters
stacked_rasters <- make_stack_rasters(covariate_layers = all_cov_layers, #raster layers and bricks
                                      period           = min(period_map[, period_id]):max(period_map[, period_id]),
                                      child_models     = child_model_objs,
                                      stacker_model    = stacked_results[[2]],
                                      indicator_family = indicator_family,
                                      return_children  = TRUE,
                                      centre_scale_df  = covs_cs_df)

# save original raster stacks in case we need them later (for example, if you're plotting multiple countries/region)
stacked_rasters_region <- copy(stacked_rasters)
all_cov_layers_region <- copy(all_cov_layers)
# ---------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------
# 2. Extract pixel in country or admin that's causing problems

# Load simple polygon for country that's causing porblems
gaul_list           <- get_adm0_codes(country, shapefile_version = modeling_shapefile_version)
simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
subset_shape        <- simple_polygon_list[[1]]
simple_polygon      <- simple_polygon_list[[2]]

# crop covariate rasters
covs <- names(all_cov_layers)
for (cov in covs) {
  all_cov_layers[[cov]] <- crop(all_cov_layers[[cov]], extent(simple_polygon))
}

# crop stacker rasters
stacks <- names(stacked_rasters)
for (stack in stacks) {
  stacked_rasters[[stack]] <- crop(stacked_rasters[[stack]], extent(simple_polygon))
}

# plot to double-check
plot(all_cov_layers[[1]])

# exclude non-time varying covariates
for (cov in covs) {
  if(nlayers(all_cov_layers[[cov]]) == 1) all_cov_layers[[cov]] <- NULL
}
covs <- names(all_cov_layers)

# exclude stacked results
stacked_rasters[['stacked_results']] <- NULL
# ---------------------------------------------------------------------------------------------------



# ---------------------------------------------------------------------------------------------------
# 3. Plot a pixel to investigate which covariates follow these weird patterns

# set up data table 
dt <- data.table(year = 2000:2017)

# convenience
stackers <- stacker_names

# extract values from rasters for a pixel by lat-long
for (cov in covs) {
  mat <- extract(all_cov_layers[[cov]], SpatialPoints(cbind(long, lat)))
  dt[, (cov) := mat[1, ]]
}
for (stack in stackers) {
  mat <- extract(stacked_rasters[[stack]], SpatialPoints(cbind(long, lat)))
  dt[, (stack) := mat[1, ]]
}

# reshape for group plotting in ggplot
dt <- melt(dt, id.vars = 'year', measure.vars = c(stackers, covs))

# ---------------------------------------------------------------------------------------------------
# 4. Save results

# create directory to save in
savedir <- paste0(outputdir, 'covariate_diagnostics/')
dir.create(savedir)

# save stacked rasters and covariates
save(list = c('stacked_rasters', 'all_cov_layers'), file = paste0(savedir, 'cov_stacker_rasters_', country, '_', pixel_id, '.RData'))

# save plot
pdf(paste0(savedir, 'cov_stacker_plots_by_pixel_', country, '_', pixel_id, '.pdf'), width = 12, height = 7)

#try a different plotting scheme
# plot stackers over time
p1 <- ggplot(data = dt[variable == 'gam' | variable == 'enet' | variable == 'xgboost' | variable == 'gbm'], 
             aes(x = year, y = value, group = variable)) + 
  geom_line(aes(color = variable)) + 
  geom_point(aes(color = variable)) +
  ggtitle(paste0('Stackers for ', toupper(indicator), ' in a pixel at ', lat, ', ', long, ' in ', country))
plot(p1)

for (cov in covs){
  p2 <- ggplot(data = dt[variable == cov], 
               aes(x = year, y = value, group = variable)) + 
    geom_line(aes(color = variable)) + 
    geom_point(aes(color = variable)) +
    ggtitle(paste0(cov, ' at ', lat, ', ', long, ' in ', country))
  plot(p2)
}

dev.off()
# ---------------------------------------------------------------------------------------------------
