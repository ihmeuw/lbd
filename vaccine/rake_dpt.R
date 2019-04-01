###############################################################################
###############################################################################
## MBG Raking Script
##
## Purpose: This script calibrates ("rakes") modeled outputs to national-level 
##          estimates of vaccine coverage from the Global Burden of Disease 
##          Study, such that the population-weighted average of the modeled
##          pixel-level estimates is equal to the GBD estimate 
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
indicator_group <- 'vaccine'

# Load objects from qsub
load_from_parallelize()

## Read config file (from sharedir) and save all parameters in memory
config <- load_config(repo            = indic_repo,
                      indicator_group = indicator_group,
                      indicator       = paste0(vaccine, "3_cov"),
                      post_est_only   = TRUE,           
                      run_date        = run_date)

## Ensure you have defined all necessary settings in your config
check_config()

## Create a few objects from options above
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
gaul_list <- get_gaul_codes(region)
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
test <- as.logical(test)

###############################################################################
## Load GBD Estimates
###############################################################################

## Check to see if input_data objects have already been saved; save if not
message("Checking to see if input_data csv files present & creating if not...")
for (ind in c(paste0(vaccine, "3_cov"), paste0(vaccine, "1_cov"), paste0(vaccine, "1_3_abs_dropout"))) {
  
  check_input_data(indicator = ind,
                   indicator_group = indicator_group,
                   age = 0,
                   reg = region,
                   run_date = run_date,
                   use_share = F) 
}

## Load GBD Estimates for this indicator which will be used in raking

dose_3 <- list()
dose_12_cond <- list()

if (vaccine == "dpt") {

  # Load 3rd dose (3_cov)
  gbd3_cov <- load_newest_gbd_vax(vaccine = paste0(vaccine, "3"),
                             gaul_list = gaul_list,
                             years = year_list,
                             gbd_date = gbd_date)

  setnames(gbd3_cov, "mean", "gbd3_cov")
  
  # Load 1+ dose (1_cov)
  gbd1_cov <- load_newest_gbd_vax(vaccine = paste0(vaccine, "1"),
                                  gaul_list = gaul_list,
                                  years = year_list,
                                  gbd_date = gbd_date)

  setnames(gbd1_cov, "mean", "gbd1_cov")

  gbd <- merge(gbd3_cov, gbd1_cov) %>% as.data.table

  gbd_unrealistic <- subset(gbd, gbd1_cov < gbd3_cov)
  if (nrow(gbd_unrealistic) > 0) {
    warning("You have some GBD values where 1_cov < 3_cov!!")
    message("Replacing 1_cov < 3_cov with 1_cov = 3_cov+(1-3_cov)*0.1")
    message("Affected countries & years:")
    print(gbd_unrealistic)
    gbd[gbd1_cov < gbd3_cov, gbd1_cov := (gbd3_cov + ((1-gbd3_cov) * 0.1))]
  }

  # Calculate gbd_12_cond (dose = 1|2 / dose = 0|1|2)
  gbd[, gbd_12_cond := (gbd1_cov - gbd3_cov)/(1-(gbd3_cov))]

  # Return doses and put in lists
  format_gbd <- function(colname) {
    d <- subset(gbd, select = c("name", "year", colname))
    setnames(d, colname, "mean")
    return(d)
  }

  dose_3[["gbd"]] <- format_gbd("gbd3_cov")
  dose_12_cond[["gbd"]] <- format_gbd("gbd_12_cond")

} else {

  stop("No GBD estimates currently available for this indicator...")

}

###############################################################################
## Load Rasters, Etc.
###############################################################################

reg <- region # convenience

## LOAD CELL PREDS ------------------------------------------------------------
message('Loading data for dose 3...')
ind <- paste0(vaccine, "3_cov")
main_dir <- paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', ind, '/output/', run_date, '/')
load(paste0(main_dir, ind, '_cell_draws_eb_bin0_', reg, '_0.RData'))
dose_3_cell_pred <- cell_pred
rm(cell_pred)

message('Loading data for dose 0...')
ind <- paste0(vaccine, "0_dose")
main_dir <- paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', ind, '/output/', run_date, '/')
load(paste0(main_dir, ind, '_cell_draws_eb_bin0_', reg, '_0.RData'))
dose_0_cell_pred <- cell_pred
rm(cell_pred)

# Calculate dose_12_cond (dose = 1|2 / dose = 0|1|2)
dose_12_cond_cell_pred <- (1-(dose_0_cell_pred + dose_3_cell_pred))/(1-dose_3_cell_pred)
dose_12_cond[["cell_pred"]] <- dose_12_cond_cell_pred
rm(dose_12_cond_cell_pred)

# Store dose_3 cell pred
dose_3[["cell_pred"]] <- dose_3_cell_pred
rm(dose_3_cell_pred)

# Clean up
rm(dose_0_cell_pred)

## PREPARE RASTERS, ETC. ------------------------------------------------------

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

## POPULATION WEIGHTS ---------------------------------------------------------

## Create population weights using the annual brick and feed custom year argument to aggregation function
pop_wts_adm0 <- make_population_weights(admin_level   = 0,
                                        simple_raster = simple_raster,
                                        pop_raster    = pop_raster_annual,
                                        gaul_list     = get_gaul_codes(reg))


# Run make_condSim() on all the indicators
dose_3[["cond_sim_raw_adm0"]] <- make_condSim(pop_wts_object = pop_wts_adm0,
                                              gaul_list      = get_gaul_codes(reg),
                                              admin_level    = 0,
                                              cell_pred      = dose_3[["cell_pred"]],
                                              summarize      = TRUE,
                                              years          = year_list)

dose_12_cond[["cond_sim_raw_adm0"]] <- make_condSim(pop_wts_object = pop_wts_adm0,
                                                    gaul_list      = get_gaul_codes(reg),
                                                    admin_level    = 0,
                                                    cell_pred      = dose_12_cond[["cell_pred"]],
                                                    summarize      = TRUE,
                                                    years          = year_list)

###############################################################################
## Raking
###############################################################################

## Get raking factors
dose_3[["rf"]] <- calc_logit_raking_factors(pop_wts_object = pop_wts_adm0,
                                            gaul_list      = get_gaul_codes(reg),
                                            admin_level    = 0,
                                            cell_pred      = dose_3[["cell_pred"]],
                                            rake_to        = dose_3[["gbd"]],
                                            years          = year_list,
                                            cores_to_use   = 4)

dose_12_cond[["rf"]] <- calc_logit_raking_factors(pop_wts_object = pop_wts_adm0,
                                                  gaul_list      = get_gaul_codes(reg),
                                                  admin_level    = 0,
                                                  cell_pred      = dose_12_cond[["cell_pred"]],
                                                  rake_to        = dose_12_cond[["gbd"]],
                                                  years          = year_list,
                                                  cores_to_use   = 4)
 
## Make raked cell pred
dose_3[["raked_cell_pred"]] <- rake_predictions(raking_factors = dose_3[["rf"]],
                                                pop_wts_object = pop_wts_adm0,
                                                cell_pred      = dose_3[["cell_pred"]],
                                                logit_rake     = TRUE)

dose_12_cond[["raked_cell_pred"]] <- rake_predictions(raking_factors = dose_12_cond[["rf"]],
                                                      pop_wts_object = pop_wts_adm0,
                                                      cell_pred      = dose_12_cond[["cell_pred"]],
                                                      logit_rake     = TRUE)

## Clean out unraked cell preds for memory reasons
dose_3[["cell_pred"]] <- NULL
dose_12_cond[["cell_pred"]] <- NULL

## Calculate other cell preds of interest
save_rcp <- function(raked_cell_pred, ind, ig = indicator_group, rd = run_date) {
  message(paste0("Saving raked cell pred object for ", ind, "..."))
  save(raked_cell_pred, file = paste0("<<<< FILEPATH REDACTED >>>>/", ig, "/", 
                                    ind, "/output/", rd, "/", 
                                    ind, "_raked_cell_draws_eb_bin0_", reg, "_0.RData" ))
}

### Vaccine3_cov (p >= 3 doses)
save_rcp(dose_3[["raked_cell_pred"]], paste0(vaccine, "3_cov"))

### Absolute dropout
### Math: (p(dose != 3)) * (p(dose = 1 or 2 | dose != 3))
abs_dropout_rcp <- dose_12_cond[["raked_cell_pred"]]*(1-dose_3[["raked_cell_pred"]])
save_rcp(abs_dropout_rcp, paste0(vaccine, "1_3_abs_dropout"))
rm(abs_dropout_rcp)

### Relative dropout
### Math: abs_dropout / (abs_dropout + dose 3)
rel_dropout_rcp <- (dose_12_cond[["raked_cell_pred"]]*(1-dose_3[["raked_cell_pred"]]))/((dose_12_cond[["raked_cell_pred"]]*(1-dose_3[["raked_cell_pred"]]))+dose_3[["raked_cell_pred"]])
save_rcp(rel_dropout_rcp, paste0(vaccine, "1_3_rel_dropout"))
rm(rel_dropout_rcp)

### Vaccine1_cov (p>= 1 dose)
### Math: abs dropout + dose 3
vax1_cov_rcp <- (dose_12_cond[["raked_cell_pred"]]*(1-dose_3[["raked_cell_pred"]])) + dose_3[["raked_cell_pred"]]
save_rcp(vax1_cov_rcp, paste0(vaccine, "1_cov"))
rm(vax1_cov_rcp)


