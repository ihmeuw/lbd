###############################################################################
###############################################################################
## MBG Main Launch Script
##
## This is the main launch script.  Will run individual models for each needed
## vaccine or conditional vaccine data set for ordinal regression using a
## continuation-ratio model, then combine these arithmatically in order to
## the desired vaccine coverage metrics. 
##
###############################################################################
###############################################################################

###############################################################################
## SETUP
###############################################################################

## clear environment
rm(list=ls())

## Set repo location and indicator group
user               <- Sys.info()['user']
core_repo          <- sprintf("FILEPATH")
indic_repo         <- sprintf("FILEPATH")
remote             <- 'origin'
branch             <- 'develop'
pullgit            <- FALSE

## sort some directory stuff
commondir      <- sprintf("FILEPATH")
package_list <- c(t(read.csv(sprintf("FILEPATH"),header=FALSE)))

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

# Custom load indicator-specific functions
source(paste0(indic_repo,'functions/misc_vaccine_functions.R'))

## Script-specific stuff begins here ##########################################
indicator_group <- "vaccine"
coverage <- TRUE # if false, modeling ratio! i.e. the timeliness ratio
vacc <- 'mcv'

if(coverage == TRUE){
  vaccine <- vacc
  set_up_indicators(stem = vaccine, 
                  doses = 1, # if modeling mcv, this should be 1
                  single_dose = T, # if modeling mcv, this should be TRUE
                  save_doses = F,
                  save_2_cov = F)
}


if(vacc == "timeliness"){
  vaccine <- paste0("dpt3_timeliness_ratio")
  model_indicators <- vaccine
  dose_indicators <- vaccine
  rake_indicators <- vaccine
  postest_indicators <- vaccine
  all_indicators <- vaccine
}

 
slots              <- 4

## Create run date in correct format
run_date <- make_time_stamp(TRUE)

# Define a log directory and clean out any files that haven't been touched in the last week
log_dir <- paste0("FILEPATH")
dir.create(log_dir, recursive = T, showWarnings = F)
system(paste0("FILEPATH"), intern=T)

  config <- load_config(repo            = indic_repo,
                    indicator_group = "",
                    indicator       = NULL,
                    config_name     = paste0("config_", vaccine), 
                    covs_name       = paste0(vaccine, "_covs"),
                    run_test        = FALSE) 

## Ensure you have defined all necessary settings in your config
check_config()

# distribute config to all indicator directories
distribute_config(cfg = config, indicators = all_indicators)

# parse formatting
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(summstats) == "character" & length(summstats) == 1) summstats <- eval(parse(text=summstats))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

###############################################################################
## Make Holdouts
###############################################################################
if(as.logical(makeholdouts)){

  if(as.logical(load_other_holdouts)){
    message(paste0("Loading holdouts from ", holdout_rundate, "..."))

    holdout_rundate <- '2020_07_19_19_07_38'

    stratum_filename <- paste0("FILEPATH") 
    
    stratum_ho <- readRDS(stratum_filename)

    saveRDS(stratum_ho,
            file = paste0("FILEPATH"))

  } else {
    df <- load_input_data(indicator   = paste0("mcv1_cov"),
                          #simple      = NULL,
                          removeyemen = FALSE,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list)


    # load the full input data for last dose of vaccine in series
    df <- load_input_data_special_cohorts(indicator   = paste0("mcv1_cov"),
                          #simple      = NULL,
                          removeyemen = FALSE,
                          years       = yearload,
                          withtag     = as.logical(withtag),
                          datatag     = datatag,
                          use_share   = as.logical(use_share),
                          yl          = year_list)

    # add in location information
    df <- merge_with_ihme_loc(df, shapefile_version = modeling_shapefile_version)




    # add in location information
    df <- merge_with_ihme_loc(df, shapefile_version = modeling_shapefile_version)

    indicator <- paste0(vaccine, doses, "_cov") # referenced internally in make_folds

  # Load simple polygon template by region
    gaul_list           <- get_adm0_codes(Regions, shapefile_version = modeling_shapefile_version)
    simple_polygon_list <- load_simple_polygon(gaul_list = gaul_list, buffer = 1, tolerance = 0.4, shapefile_version = modeling_shapefile_version)
    subset_shape        <- simple_polygon_list[[1]]
    # Load simple raster by region
    raster_list        <- build_simple_raster_pop(subset_shape)
    simple_raster      <- raster_list[['simple_raster']]

    # load admin 2 shapefile and crop by region
    shapefile_admin <- shapefile(get_admin_shapefile(admin_level = 2, version = modeling_shapefile_version))
    shapefile_admin <- shapefile_admin[which(shapefile_admin$ADM0_CODE %in% gaul_list),]
    # load mask raster and crop by region
    raster_mask <- raster("FILEPATH")
    raster_mask <- crop(raster_mask, extent(simple_raster))
    # make a list of dfs for each region, with 5 qt folds identified in each    
    lat_col = 'latitude'
    long_col = 'longitude';
    ss_col = 'N'
    yr_col = 'year'
    admin_raster <- simple_raster
    admin_shps=shapefile_admin
shape_ident="ADM1_CODE"

    stratum_ho <- make_folds(data = df, n_folds = 5, spat_strat = 'poly',
                              temp_strat = NULL, strat_cols = 'region',
                              #admin_raster=simple_raster,
                              shape_ident="ADM1_CODE",
                              ts         = as.numeric(ho_ts),
                              mb         = as.numeric(ho_mb),
                              admin_shps=shapefile_admin,
                              admin_raster=simple_raster,
                              mask_shape=subset_shape,
                              mask_raster=raster_mask,
                              lat_col = 'latitude', long_col = 'longitude',
                              ss_col = 'N', yr_col = 'year', seed = 98112)
    # Recreate & distribute holdouts
  } 

  indicators <- indicator
  invisible(lapply((unique(c(indicators, postest_indicators))), function(ind) {

    sharedir <- sprintf("FILEPATH")

    ## Create holdouts if not already present
    if(as.logical(makeholdouts) & !(file.exists(paste0(sharedir,'/output/',run_date,"/stratum.rds")))){

      message(paste0("Recreating holdouts for ", ind, "..."))
      # load the full input data
      df <- load_input_data(indicator   = ind,
                            simple      = NULL,
                            removeyemen = FALSE,
                            years       = yearload,
                            withtag     = as.logical(withtag),
                            datatag     = datatag,
                            use_share   = as.logical(use_share),
                            yl          = year_list)

      # add in location information
      df <- merge_with_ihme_loc(df, shapefile_version = modeling_shapefile_version)

      # make a list of dfs for each region, with 5 qt folds identified in each

      # New way: use the loaded data and the dpt3_cov holdouts to assign each
      #          row in the data to the same data across modeled indicators
      stratum_ho <- recreate_holdouts(data = df,
                                      row_id_col = "row_id",
                                      load_from_indic = paste0(vaccine, "1_cov"),
                                      rd = run_date,
                                      ig = indicator_group)

      # Save stratum_ho object for reference
      save_stratum_ho(indic       = ind,
                      ig          = indicator_group,
                      rd          = run_date,
                      stratum_obj = stratum_ho)

    } # if(as.logical(makeholdouts) & !(file.exists(...
  })) # invisible(lapply...
} #if(as.logical(makeholdouts))

###############################################################################
## Launch parallel models
###############################################################################

launch_lv = list(indicator = model_indicators,
                 use_gn = TRUE)

launch_para_output <- parallelize(script = "launch",
                                  log_location = paste0(log_dir, "launch/"),
                                  expand_vars = launch_lv,
                                  save_objs = c("run_date", "indicator_group", "vaccine", "core_repo"),
                                  prefix = "launch",
                                  slots = 2,
                                  memory = 4,
                                  script_dir = indic_repo,
                                  geo_nodes = TRUE,
                                  singularity = 'default')

monitor_jobs(launch_para_output, notification = "pushover")


###############################################################################
## Rake indicators
###############################################################################

# Skip some steps if validation metrics only
if (!as.logical(validation_metrics_only)) {

  pushover_notify("Master script: starting raking", 
                  title = "Raking")

  rake_lv <- list(region = Regions1,
                  holdout = 0)

  rake_para_output <- parallelize(script = "rake_fractional",
                                  log_location = paste0(log_dir, "rake_fractional/"),
                                  expand_vars = rake_lv,
                                  save_objs = c("indicator_group", "run_date", "vaccine", 
                                                "rake_indicators", "doses"),
                                  prefix = "rake_fractional",
                                  slots = 8, # geo_nodes
                                  #slots = 3,
                                  use_c2_nodes = FALSE,
                                  script_dir = indic_repo, 
                                  memory = 300,
                                  geo_nodes = T,
                                  singularity = 'default')


} 


#######################################################################################
## Post-estimation
#######################################################################################
sharedir <- sprintf("FILEPATH")

indicator <- rake_indicators

strata <- Regions

#######################################################################################
## Merge!! For raked
#######################################################################################
  rr <- c("raked")
  rf_tab <- T

  run_summary <- T
  summstats <- c('mean')

if (("p_below" %in% summstats) & !(indicator %in% c(paste0(vaccine, "1_cov"), paste0(vaccine, "3_cov")))) {
  # Only run this for 1st and 3rd doses
  summstats <- summstats[summstats != "p_below"]
}

post_load_combine_save(regions    = Regions,
                       summstats  = summstats,
                       raked      = rr,
                       rf_table   = rf_tab,
                       run_summ   = run_summary,
                       indic      = indicator,
                       ig         = indicator_group,
                       sdir       = sharedir)

# Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)


#######################################################################################
## Merge!! Now for unraked
#######################################################################################
  rr <- c("unraked")
  rf_tab <- F
  run_summary <- T

if (("p_below" %in% summstats) & !(indicator %in% c(paste0(vaccine, "1_cov"), paste0(vaccine, "3_cov")))) {
  # Only run this for 1st and 3rd doses
  summstats <- summstats[summstats != "p_below"]
}

post_load_combine_save(regions    = Regions,
                       summstats  = summstats,
                       raked      = rr,
                       rf_table   = rf_tab,
                       run_summ   = run_summary,
                       indic      = indicator,
                       ig         = indicator_group,
                       sdir       = sharedir)

# Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)


#######################################################################################
## Combine aggreagations
#######################################################################################
holdouts <- 0
combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages = 0, 
                    regions = strata,
                    holdouts = holdouts,
                    raked = c(T, F))


#######################################################################################
## Summarize
#######################################################################################

summarize_admins(summstats = c("mean", "lower", "upper", "cirange", "cfb"), 
                 ad_levels = c(0,1,2), 
                 raked = c(T,F))

message("Summarize p_above")
summarize_admins(indicator, indicator_group,
                 summstats = c("p_above"),
                 raked = c(T,F),
                 ad_levels = c(0,1,2),
                 file_addin = "p_0.8_or_better",
                 value = 0.8,
                 equal_to = T)

message("Summarize p_above")
summarize_admins(indicator, indicator_group,
                 summstats = c("p_above"),
                 raked = c(T,F),
                 ad_levels = c(0,1,2),
                 file_addin = "p_0.95_or_better",
                 value = 0.95,
                 equal_to = T)


#######################################################################################
## Save everything for mapping directory
#######################################################################################

postest_indicators <- indicator
#summstats <- summstats[-5]

  tifs_to_copy <- data.table(expand.grid(summstats, c(T), postest_indicators, stringsAsFactors = F))
  names(tifs_to_copy) <- c("measure", "raked", "indicator") 

  for (i in 1:nrow(tifs_to_copy)) {
    copy_tif_to_map_input_dir(ind = tifs_to_copy[i, "indicator"],
                              ig = indicator_group,
                              measure = tifs_to_copy[i, "measure"],
                              rd = run_date,
                              raked = tifs_to_copy[i, "raked"],
                              yl = year_list)
  }

  # Save ads
  admin_lvs <- data.table(expand.grid(summstats, c(T), c(0,1,2), postest_indicators, stringsAsFactors = F))
  admin_lvs <- rbind(admin_lvs,
                     data.table(expand.grid("psup80", 
                                            c(T,F), 
                                            c(0,1,2),
                                            postest_indicators, 
                                stringsAsFactors = F)))

  admin_lvs <- rbind(admin_lvs,
                     data.table(expand.grid("psup95", 
                                            c(T,F), 
                                            c(0,1,2),
                                            postest_indicators, 
                                stringsAsFactors = F)))

  names(admin_lvs) <- c("measure", "raked", "ad_level", "indicator") 

  for (i in 1:nrow(admin_lvs)) {
    copy_admins_to_map_input_dir(ind = admin_lvs[i, "indicator"],
                                 ig = indicator_group, 
                                 measure = as.character(admin_lvs[i, "measure"]), 
                                 rd = run_date, 
                                 raked = admin_lvs[i, "raked"], 
                                 yl = year_list, 
                                 ad_level = admin_lvs[i, "ad_level"]) 
  }


#######################################################################################
## In sample metric validation
#######################################################################################

# Combine csv files only if none present
csvs <- list.files(paste0(sharedir, '/output/', run_date, '/'), 
                   pattern = "input_data(.*).csv", 
                   full.names = T)

if (sum(grepl("input_data.csv", csvs)) == 0) {
  csv_master <- lapply(csvs, read.csv, stringsAsFactors = F) %>% 
                rbindlist %>%
                subset(., select = !(names(.) %in% c("X.2", "X.1"))) %>%
                as.data.table
  write.csv(csv_master, file=paste0(sharedir, '/output/', run_date, '/input_data.csv'))
}

# Get in and out of sample draws
run_in_oos <- get_is_oos_draws(ind_gp = indicator_group,
                               ind = indicator,
                               rd = run_date,
                               ind_fm = 'binomial',
                               model_domain = Regions,
                               age = 0,
                               nperiod = length(year_list),
                               yrs = year_list,
                               get.oos = as.logical(makeholdouts),
                               write.to.file = TRUE,
                               year_col = "year",
                               shapefile_version = modeling_shapefile_version) # uses unraked

## set out_dir
out_dir <- paste0(sharedir, "/output/", run_date, "/summary_metrics/")
dir.create(out_dir, recursive = T, showWarnings = F)

## set up titles
if (indicator == "mcv1_cov") plot_title <- "MCV1 Coverage"

draws.df <- fread(sprintf("FILEPATH"))

# By region
pvtable.reg <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            result_agg_over = c("year", "oos", "region"),
                            coverage_probs = seq(from=5, to=95, by=10),
                            aggregate_on=c("country", "ad1", "ad2"),
                            draws = sum(grepl("draw[0-9]+",names(draws.df))),
                            out.dir = out_dir,
                            plot = TRUE,
                            plot_by = "region",
                            plot_by_title = "Region",
                            plot_ci = TRUE,
                            point_alpha = 0.5,
                            point_color = "black",
                            ci_color = "gray",
                            plot_title = plot_title)

# All regions
pvtable.all <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            result_agg_over = c("year", "oos"),
                            coverage_probs = seq(from=5, to=95, by=10),
                            aggregate_on=c("country", "ad1", "ad2"),
                            draws = sum(grepl("draw[0-9]+",names(draws.df))),
                            out.dir = out_dir,
                            plot = TRUE,
                            plot_ci = TRUE,
                            point_alpha = 0.1,
                            point_color = "black",
                            ci_color = "gray",
                            plot_title = plot_title)

