# Parallel script for post-estimation
#
# Note that prep_postest is used to pass objects to the parallel script
# (saved in a temporarly RData file).  If you want to use logit raking,
# pass "rake_transform" as a string that equals "logit".  You can set this
# up easily in the config file with a row for rake_transform with value "logit"
#
# Sample use:
#
# # Prepare for parallel post-estimation - save file with objs to re-load
# prep_postest(indicator = indicator,
#              indicator_group = indicator_group,
#              run_date = run_date,
#              save_objs = c("core_repo", "gbd", "year_list", "summstats", "rake_transform", "pop_measure", "pop_release"))
#
# ## Parallelized post-estimation over region
# postest_script <- "postest_script"
#
# for (s in strata) {
#  qsub <- make_qsub_postest(code = postest_script,
#                            stratum = s,
#                            log_location = 'sharedir',
#                            memory        = mem,
#                            cores         = cores)
#  system(qsub)
# }
#
# ## check to make sure post-est done before continuing
# waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

## SETUP ################################################################################

stratum=as.character(commandArgs()[4])
run_date=as.character(commandArgs()[5])
indicator=as.character(commandArgs()[6])
indicator_group=as.character(commandArgs()[7])
geos_node = as.logical(commandArgs()[8])

interval_mo <- 12

# Define directories
main_dir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date, '/')
temp_dir <- paste0(main_dir, "temp_post_est/")

# Load objects from convenience temp file
load(paste0(temp_dir, "post_est_temp_objs.RData"))

# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(stratum)
message(pop_measure)
message(pop_release)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
reg <- stratum

# do we really need this again?
sharedir       <- sprintf('/share/geospatial/mbg/%s/%s',indicator_group,indicator)
commondir      <- sprintf('/share/geospatial/mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Get the necessary variables out from the config object into global env
rake_subnational <- eval(parse(text = config[V1 == "subnational_raking", V2]))
modeling_shapefile_version <- config[V1 == "modeling_shapefile_version", V2]
raking_shapefile_version <- config[V1 == "raking_shapefile_version", V2]
countries_not_to_rake <- config[V1 == "countries_not_to_rake", V2]
countries_not_to_subnat_rake <- config[V1 == "countries_not_to_subnat_rake", V2]
rake_countries <- eval(parse(text = config[V1 == "rake_countries", V2]))

## PREPARE RASTERS, ETC. ################################################################

# Load cell draws
message('Loading Data...')
load(paste0(main_dir, indicator, '_cell_draws_eb_bin0_', reg, '_0.RData'))

# Check if load was successful; stop if not
if (!exists("cell_pred")){
  message(filename_rds)
  stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")
}

# Rake estimates 
if (rake_countries) {
  if(is.null(gbd)){
    stop("rake_countries was specified as T in config, gbd raking targets must be provided.")
  }

  ## determine if a crosswalk is needed
  if (modeling_shapefile_version == raking_shapefile_version) crosswalk <- F else crosswalk <- T
  
  # Assume linear raking unless specified as logit
  if (rake_transform == "logit") rake_method <- "logit" else rake_method <- "linear" 
  
  #if countries not to rake are specified, sets "name" in gbd (raking_targets) to -1. 
  #the if_no_gbd parameter in rake_cell_pred will return unraked estimates if no targets are found for that country
  if(!is.null(countries_not_to_rake)) {
    connector <- get_gbd_locs(
      rake_subnational = F,
      reg = reg,
      shapefile_version = raking_shapefile_version
    )
    not_to_rake_adm_codes <- get_adm0_codes(countries_not_to_rake, shapefile_version=raking_shapefile_version)
    not_to_rake_gbd_loc_ids <- connector[ADM_CODE %in% not_to_rake_adm_codes, location_id]
    if(length(not_to_rake_gbd_loc_ids) > 0){
      setDT(gbd)
      gbd[name %in% not_to_rake_gbd_loc_ids, name := -1]
    }
  }
  
  #make custom raking shapefile if countries_not_to_subnat_rake is set in config
  if (rake_subnational & !is.null("countries_not_to_subnat_rake")){
    custom_raking_shapefile <- make_custom_raking_shapefile(countries_not_to_subnat_rake, raking_shapefile_version)
  } else {
    custom_raking_shapefile <- NULL
  }
  
  #Rake cell pred
  outputs <-
    rake_cell_pred(cell_pred = cell_pred,
                   rake_to = gbd,
                   reg = reg,
                   year_list = year_list,
                   pop_measure = pop_measure,
                   rake_method = rake_method,
                   crosswalk = crosswalk,
                   shapefile_path = get_admin_shapefile(admin_level = 0,
                                                        raking = rake_subnational, ## if true, subnat raking shp is grabbed. otherwise, admin0 is grabbed
                                                        version = raking_shapefile_version), 
                   rake_subnational = rake_subnational,
                   modeling_shapefile_version = modeling_shapefile_version,
                   raking_shapefile_version = raking_shapefile_version, 
                   if_no_gbd = "return_unraked",
                   custom_raking_shapefile = custom_raking_shapefile,
                   countries_not_to_subnat_rake = countries_not_to_subnat_rake)
  
  # Extract raking factor and raked cell pred
  rf <- outputs[["raking_factors"]]
  raked_cell_pred <- outputs[["raked_cell_pred"]]
  
  # Define raking mask to be used with raked_cell_pred (needed if crosswalk occured)
  raked_simple_raster <- outputs[["new_simple_raster"]]
  # we keep the original simple_raster for unraked cell preds
  simple_raster <- outputs[["simple_raster"]] 
  
} else {
  rf <- NULL
  raked_cell_pred <- NULL
  
  # Define simple raster for mask
  simple_polygon_list <- load_simple_polygon(gaul_list = get_adm0_codes(reg,
                                                                        shapefile_version = modeling_shapefile_version),
                                             buffer = 0.4, subset_only = FALSE,
                                             shapefile_version = modeling_shapefile_version)
  subset_shape   <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list    <- build_simple_raster_pop(subset_shape)
  simple_raster  <- raster_list[['simple_raster']]
  rm(simple_polygon_list)
}

## SAVE THE RESULTS #####################################################################
message('Saving results...')

## save RF
save_post_est(rf, "csv", paste0(reg, "_rf"))

## save raked cell preds
save(raked_cell_pred, file = paste0(sharedir, "/output/", run_date, "/",
                                    indicator, "_raked_cell_draws_eb_bin0_", reg, "_0.RData" ))

# make and save summaries

save_cell_pred_summary <- function(summstat, raked, ...) {
  message(paste0('Making summmary raster for: ',summstat, " (", raked, ")"))
  if (raked=="unraked"){
    cpred <- "cell_pred"
    mask_raster <- 'simple_raster'
  }
  if (raked=="raked"){
    cpred <- "raked_cell_pred"
    mask_raster <- 'raked_simple_raster'
  }
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask                 = get(mask_raster), 
    return_as_raster     = TRUE,
    summary_stat         = summstat,
    ...)
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ""),'_',summstat,'_raster'))
}

# Do this as lapply to not fill up memory in global env with big obs
if (is.null(gbd))  rake_list <- c("unraked")
if (!is.null(gbd)) rake_list <- c("unraked", "raked")

summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)

lapply(1:nrow(summ_list), function(i) {
  summstat <- as.character(summ_list[i, 1])
  raked <- as.character(summ_list[i, 2])
  save_cell_pred_summary(summstat, raked)
})

## Can't pass additional params in the above framework, so will code by hand here
for (r in rake_list) {
  if ("p_below" %in% summstats) {
    save_cell_pred_summary(summstat = "p_below",
                           raked = r,
                           value = 0.8,
                           equal_to = F)
  }
}

# Write a file to mark done
output_dir <- paste0('/share/geospatial/mbg/', indicator_group, '/', indicator, '/output/', run_date)
pathaddin <- paste0('_bin0_',reg,'_0') # To allow us to use waitformodelstofinish()
write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

# All done
message(paste0("Done with post-estimation for ", stratum))
