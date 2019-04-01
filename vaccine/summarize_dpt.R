###############################################################################
###############################################################################
## MBG Summarization Script
##
## Purpose: Summarize draw-level objects into means, upper and lower CIs, 
##          and other summary statistics
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- '<<<< FILEPATH REDACTED >>>>'
indic_repo         <- '<<<< FILEPATH REDACTED >>>>'

## sort some directory stuff
commondir      <- '<<<< FILEPATH REDACTED >>>>'
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific code begins here ##########################################

# Load from qsub
load_from_parallelize()

## Read config file (from sharedir) and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,           
                      run_date        = run_date)

## Ensure you have defined all necessary settings in your config
check_config()

## Create a few objects from options above
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
gaul_list <- get_gaul_codes(Regions)
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
test <- as.logical(test)

## Convenience
reg <- region
rundir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", indicator, "/output/", run_date, "/")

## Load gaul list
gaul_list <- get_gaul_codes(region)

#######################################################################################
## Load Objects
#######################################################################################

message("\n######################################################################")
message(paste0("Summarizing ", indicator, " | ", region, " | ", raked))
message("\n######################################################################")

## Prepare rasters, etc. --------------------------------------------------------------

message('Loading Map Objects...')
# Get aggregated estimates for all admin0. Aggregate to level you rake to
simple_polygon_list <- load_simple_polygon(gaul_list = get_gaul_codes(reg), buffer = 0.4, subset_only = FALSE)
subset_shape   <- simple_polygon_list[[1]]
simple_polygon <- simple_polygon_list[[2]]
raster_list    <- build_simple_raster_pop(subset_shape)
simple_raster  <- raster_list[['simple_raster']]
pop_raster     <- raster_list[['pop_raster']]

## Pull 2000-2016 annual population brick using new covariates function
pop_raster_annual <- load_and_crop_covariates_annual(covs = 'worldpop',                
                                                     measures = pop_measure,          
                                                     simple_polygon = simple_polygon,
                                                     start_year  = min(year_list),
                                                     end_year    = max(year_list),
                                                     interval_mo = as.numeric(interval_mo),
                                                     agebin=1)
pop_raster_annual  <- pop_raster_annual[[1]]
pop_raster_annual  <- crop(pop_raster_annual, extent(simple_raster))
pop_raster_annual  <- setExtent(pop_raster_annual, simple_raster)
pop_raster_annual  <- mask(pop_raster_annual, simple_raster)

## Create population weights using the annual brick and feed custom year argument to aggregation function
pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                        simple_raster = simple_raster,
                                        pop_raster    = pop_raster_annual,
                                        gaul_list     = get_gaul_codes(reg))

#######################################################################################
## RF table (if needed)
#######################################################################################

## Make a raking factors table --------------------------------------------------------

if (raked == T) {
  ## Raked indicators only
  message("Making a raking factors table...")

  # Get GBD for raking factor csvs
  message("  Loading GBD estimates...")
  if (grepl("3_cov", indicator)) {
    
    gbd <- load_newest_gbd_vax(vaccine = paste0(vaccine, "3"),
                               gaul_list = gaul_list,
                               years = year_list,
                               gbd_date = gbd_date)

  } else if (grepl("1_cov", indicator)) {

    gbd <- load_newest_gbd_vax(vaccine = paste0(vaccine, "1"),
                               gaul_list = gaul_list,
                               years = year_list,
                               gbd_date = gbd_date)

  } else if (grepl("1_3_abs_dropout", indicator)) {

    gbd <- load_dropout_for_raking(vaccine = vaccine,
                                   gaul_list,
                                   years = year_list,
                                   abs_rel = "absolute",
                                   gbd_date = gbd_date)

  } else if (grepl("1_3_rel_dropout", indicator)) {

    gbd <- load_dropout_for_raking(vaccine = vaccine,
                                 gaul_list,
                                 years = year_list,
                                 abs_rel = "relative",
                                 gbd_date = gbd_date)
  }

  # Load the unraked cell pred for condsim
  message("  Loading unraked cell_pred estimates...")
  load(paste0(rundir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))

  # Run condsim
  message("  Run make_condSim()...")
  condsim <- make_condSim(pop_wts_object = pop_wts_adm0,
                          gaul_list      = get_gaul_codes(reg),
                          admin_level    = 0,
                          cell_pred      = cell_pred,
                          summarize      = TRUE,
                          years          = year_list)

  rm(cell_pred)

  ## Get rf into standard format - get aggregated values to append to rf
  message("  Making and saving raking factor table...")
  aggs <- get_aggs(agg_geo_est = condsim,
                   gaul_list   = get_gaul_codes(reg),
                   rake_to     = gbd)

  aggs$name <- as.character(aggs$name)
  aggs$year <- as.character(aggs$year)

  ## Additional format (no raking factors apply)
  aggs[,raking_factor := NA]
  aggs[,country_year := paste0(name, "_", year)]
  rf <- aggs

  ## save RF
  save_post_est(rf, "csv", paste0(reg, "_rf"))
}

## Create summaries -------------------------------------------------------------------
## Load cell pred

if (raked == F) load(paste0(rundir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))
if (raked == T) load(paste0(rundir, indicator, '_raked_cell_draws_eb_bin0_', reg, '_0.RData'))

message(paste0("Making summary outputs for ", indicator, " | ", region , "..."))

if (raked == T) rake_list <- "raked"
if (raked == F) rake_list <- "unraked"
summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list, stringsAsFactors = F)

message("  Saving summary statistic rasters...")

for (i in 1:nrow(summ_list)) {
  summstat <- as.character(summ_list[i, 1])
  raked <- as.character(summ_list[i, 2])

  message(paste0('Making summary raster for: ',summstat, " (", raked, ")"))
  if (raked=="unraked") cpred <- "cell_pred"
  if (raked=="raked")   cpred <- "raked_cell_pred"
 
  ras <- make_cell_pred_summary(
          draw_level_cell_pred = get(cpred),
          mask                 = simple_raster,
          return_as_raster     = TRUE,
          summary_stat         = summstat)

  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ""),'_',summstat,'_raster'))
  rm(ras)

}

## Can't pass additional params in the above framework, so will code by hand here
## Run only for the first and third dose indicators

for (r in rake_list) {
  if (("p_above" %in% summstats) & (indicator %in% c(paste0(vaccine, "1_cov"), paste0(vaccine, "3_cov")))) {
    save_cell_pred_summary(summstat = "p_above",
                           raked = r,
                           value = 0.8,
                           equal_to = T)
  }
}


for (r in rake_list) {
  if (("p_below" %in% summstats) & (indicator %in% c(paste0(vaccine, "1_cov"), paste0(vaccine, "3_cov")))) {
    save_cell_pred_summary(summstat = "p_below",
                           raked = r,
                           value = 0.8,
                           equal_to = F)
  }
}

# Submit a job to create model diagnostics (modeled indicators only)
if (indicator %in% c(paste0(vaccine, "3_cov"), paste0(vaccine, "1_cond"), paste0(vaccine, "2_cond")) & raked == F) {
  make_model_diagnostics(indic      = indicator,
                         ig         = indicator_group,
                         rd         = run_date,
                         geo_nodes  = TRUE,
                         cores      = 10)
}
