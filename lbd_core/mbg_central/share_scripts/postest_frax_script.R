# Parallel script for post-estimation and aggregation
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
#              save_objs = c("core_repo", "gbd", "year_list", "summstats", "rake_transform", "config"))
#
# ## Parallelized post-estimation over region
# postest_script <- "postest_frac_script"
#
# for (s in strata) {
#   qsub <- make_qsub_postest(
#     code = postest_script,
#     stratum = s,
#     log_location = "sharedir",
#     memory = 9.,
#     subnat_raking = subnational_raking,
#     modeling_shapefile_version = modeling_shapefile_version,
#     raking_shapefile_version = raking_shapefile_version,
#     singularity = "/share/singularity-images/lbd/testing_INLA_builds/lbd_rpkgs3.5.3gcc8mklrstudioserver1.1.463_v3.simg",
#     queue = "geospatial.q",
#     run_time = "1:00:00",
#     proj = "proj_geo_nodes",
#     cores = 2
#   )
#   system(qsub)
# }

#
# ## check to make sure post-est done before continuing
# waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

## SETUP ################################################################################

stratum <- as.character(commandArgs()[4])
run_date <- as.character(commandArgs()[5])
indicator <- as.character(commandArgs()[6])
indicator_group <- as.character(commandArgs()[7])

interval_mo <- 12

# Define directories
main_dir <- paste0("/share/geospatial/mbg/", indicator_group, "/", indicator, "/output/", run_date, "/")
temp_dir <- paste0(main_dir, "temp_post_est/")

# Load objects from convenience temp file
load(paste0(temp_dir, "post_est_temp_objs.RData"))


# do we really need this again?
sharedir <- sprintf("/share/geospatial/mbg/%s/%s", indicator_group, indicator)
commondir <- sprintf("/share/geospatial/mbg/common_inputs")
package_list <- c(t(read.csv(sprintf("%s/package_list.csv", commondir), header = FALSE)))

message("Loading in required R packages and MBG functions")
source(paste0(core_repo, "/mbg_central/setup.R"))
mbg_setup(package_list = package_list, repos = core_repo)


# Get the necessary variables out from the config object into global env
rake_countries <- eval(parse(text = config[V1 == "rake_countries", V2]))
rake_subnational <- eval(parse(text = config[V1 == "subnational_raking", V2]))
modeling_shapefile_version <- config[V1 == "modeling_shapefile_version", V2]
raking_shapefile_version <- config[V1 == "raking_shapefile_version", V2]
countries_not_to_rake <- config[V1 == "countries_not_to_rake", V2]
countries_not_to_subnat_rake <- config[V1 == "countries_not_to_subnat_rake", V2]
year_list <- eval(parse(text = config[V1 == "year_list", V2]))
metric_space <- config[V1 == "metric_space", V2]



# Print some settings to console
message(indicator)
message(indicator_group)
message(run_date)
message(stratum)
message(pop_measure)
message(pop_release)
message(paste0("Summary stats: ", paste0(summstats, collapse = ", ")))

# Print raking info
print(paste0("Metric Space                       : ", metric_space))
print(paste0("Subnational raking                 : ", rake_subnational))
print(paste0("Countries not to rake at all       : ", countries_not_to_rake))
print(paste0("Countries not to rake subnationally: ", countries_not_to_subnat_rake))



# For now just assign stratum to reg (will need to modify the below for strata beyond reg)
reg <- stratum
age <- 0
holdout <- 0


## PREPARE RASTERS, ETC. ################################################################

# Load cell draws
message("Loading Data...")
load(paste0(main_dir, indicator, "_cell_draws_eb_bin0_", reg, "_0.RData"))

# Check if load was successful; stop if not
if (!exists("cell_pred")) {
  message(filename_rds)
  stop("Unable to load cell_pred object! Check to make sure that the relevant object exists.")
}

# Rake estimates
if (rake_countries) {
  if (!exists("gbd")) {
    stop("rake_countries was specified as T in config, gbd raking targets must be provided.")
  }

  ## determine if a crosswalk is needed
  if (modeling_shapefile_version == raking_shapefile_version) crosswalk <- F else crosswalk <- T

  # Assume linear raking unless specified as logit
  if (rake_transform == "logit") rake_method <- "logit" else rake_method <- "linear"


  ##### Prep input data into raking:

  ## Get the simple and new_simple rasters prepped up for us
  print("Getting simple and prepped rasters")
  raster_outputs <- prep_shapes_for_raking(
    reg = reg,
    modeling_shapefile_version = modeling_shapefile_version,
    raking_shapefile_version = raking_shapefile_version,
    field = "loc_id"
  )

  ## Take out the objects from the list that actually matters to us:
  simple_raster <- raster_outputs[["simple_raster"]]
  new_simple_raster <- raster_outputs[["new_simple_raster"]]

  simple_polygon <- raster_outputs[["simple_polygon"]]
  new_simple_polygon <- raster_outputs[["new_simple_polygon"]]

  pixel_id <- raster_outputs[["pixel_id"]]



  ##### Using fractional raking #####

  if (metric_space == "rates") {
    print("Get GBD populations")
    gbd_pops <- prep_gbd_pops_for_fraxrake(pop_measure = pop_measure,
                                           reg = reg,
                                           year_list = year_list,
                                           gbd_round_id = 6,
                                           decomp_step = "step4")

    print("Using the rates raking and aggregation functions:")

    ## First, create all the fractional rake factors
    fractional_rake_rates(
      cell_pred = cell_pred,
      simple_raster = simple_raster,
      simple_polygon = simple_polygon,
      pixel_id = pixel_id,
      shapefile_version = raking_shapefile_version,
      reg = reg,
      pop_measure = pop_measure,
      year_list = year_list,
      interval_mo = interval_mo,
      rake_subnational = rake_subnational,
      age_group = age_group,
      sex_id = sex_id,
      sharedir = sharedir,
      run_date = run_date,
      indicator = indicator,
      gbd = gbd,
      rake_method = rake_method,
      gbd_pops = gbd_pops,
      countries_not_to_rake = countries_not_to_rake,
      countries_not_to_subnat_rake = countries_not_to_subnat_rake
    )

    ## Now, create the raked cell pred files!
    outputs <- fractional_agg_rates(
      cell_pred = cell_pred,
      simple_raster = simple_raster,
      simple_polygon = simple_polygon,
      pixel_id = pixel_id,
      shapefile_version = raking_shapefile_version,
      reg = reg,
      pop_measure = pop_measure,
      year_list = year_list,
      interval_mo = interval_mo,
      rake_subnational = rake_subnational,
      sharedir = sharedir,
      run_date = run_date,
      indicator = indicator,
      main_dir = main_dir,
      rake_method = rake_method,
      age = age,
      holdout = holdout,
      countries_not_to_subnat_rake = countries_not_to_subnat_rake,
      return_objects = TRUE
    )

    ## Get the necessary outputs and rename the columns
    rf <- data.table(outputs[["rf"]])[, .(loc = location_id, year, start_point = mbg_prev, target = gbd_prev, raking_factor = rf)]
    raked_cell_pred <- outputs[["raked_cell_pred"]]


    ## Raked simple raster has been made above
    raked_simple_raster <- new_simple_raster

  } else if (metric_space == "counts") {
    print("Using the counts raking and aggregation functions:")

    ## Rake counts
    outputs <- fractionally_rake_counts(
      count_cell_pred = data.table(cell_pred),
      rake_to = gbd,
      reg = reg,
      year_list = year_list,
      rake_subnational = rake_subnational,
      countries_not_to_subnat_rake = countries_not_to_subnat_rake,
      countries_not_to_rake = countries_not_to_rake,
      simple_raster = simple_raster,
      modeling_shapefile_version = modeling_shapefile_version,
      raking_shapefile_version = raking_shapefile_version
    )




    ## Get the necessary outputs
    rf <- outputs$raking_factors
    raked_cell_pred <- outputs$raked_cell_pred
    raked_simple_raster <- new_simple_raster

    ## Save out the aggregate files
    raked_frax_counts_save(output_list = outputs, sharedir = sharedir, indicator = indicator, age = age, reg = reg, holdout = holdout)
  }
} else {
  rf <- NULL
  raked_cell_pred <- NULL

  # Define simple raster for mask
  simple_polygon_list <- load_simple_polygon(
    gaul_list = get_adm0_codes(reg,
      shapefile_version = modeling_shapefile_version
    ),
    buffer = 0.4, subset_only = FALSE,
    shapefile_version = modeling_shapefile_version
  )
  subset_shape <- simple_polygon_list[[1]]
  simple_polygon <- simple_polygon_list[[2]]
  raster_list <- build_simple_raster_pop(subset_shape)
  simple_raster <- raster_list[["simple_raster"]]
  rm(simple_polygon_list)
}

## SAVE THE RESULTS #####################################################################
message("Saving results...")

## save RF
save_post_est(rf, "csv", paste0(reg, "_rf"))

## save raked cell preds
save(raked_cell_pred, file = paste0(
  sharedir, "/output/", run_date, "/",
  indicator, "_raked_cell_draws_eb_bin0_", reg, "_0.RData"
))

# make and save summaries

save_cell_pred_summary <- function(summstat, raked, ...) {
  message(paste0("Making summmary raster for: ", summstat, " (", raked, ")"))
  
  if (raked == "unraked") {
    cpred <- "cell_pred"
    mask_raster <- "simple_raster"
  }
  if (raked == "raked") {
    cpred <- "raked_cell_pred"
    mask_raster <- "raked_simple_raster"
  }
  if (raked == "raked_c") {
    cpred <- "raked_cell_pred_c"
    load(paste0(sharedir, "/output/", run_date, "/", indicator, "_raked_c_cell_draws_eb_bin0_", reg, "_0.RData" ))
    mask_raster <- "raked_simple_raster"
  }
  ras <- make_cell_pred_summary(
    draw_level_cell_pred = get(cpred),
    mask = get(mask_raster),
    return_as_raster = TRUE,
    summary_stat = summstat,
    ...
  )
  save_post_est(ras,'raster',paste0(reg, ifelse(raked == "raked", "_raked", ifelse(raked == 'raked_c', '_raked_c', '')), '_', summstat, '_raster'))
}

# Do this as lapply to not fill up memory in global env with big obs
if (is.null(gbd)) {
    rake_list <- c("unraked")
} else {
    rake_list <- c("unraked", "raked")
}


summ_list <- expand.grid(summstats[summstats != "p_below"], rake_list)

lapply(1:nrow(summ_list), function(i) {
  summstat <- as.character(summ_list[i, 1])
  raked <- as.character(summ_list[i, 2])
  save_cell_pred_summary(summstat, raked)
})

## Can't pass additional params in the above framework, so will code by hand here
for (r in rake_list) {
  if ("p_below" %in% summstats) {
    save_cell_pred_summary(
      summstat = "p_below",
      raked = r,
      value = 0.8,
      equal_to = F
    )
  }
}

# Write a file to mark done
output_dir <- paste0("/share/geospatial/mbg/", indicator_group, "/", indicator, "/output/", run_date)
pathaddin <- paste0("_bin0_", reg, "_0") # To allow us to use waitformodelstofinish()
write(NULL, file = paste0(output_dir, "/fin_", pathaddin))

# All done
message(paste0("Done with post-estimation and aggregation for ", stratum))
