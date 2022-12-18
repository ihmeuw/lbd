###########################################################################################
###########################################################################################
## Launch all models
##
###########################################################################################
###########################################################################################

###########################################################################################
## SETUP
###########################################################################################

## Clear environment
rm(list=ls())

## Write message for clarity
message(commandArgs()[4])
message(commandArgs()[5])
message(commandArgs()[6])
message(commandArgs()[7])

## Set indicator
indicator_group <- 'hiv'
indicator       <- commandArgs()[4]
prior           <- commandArgs()[5]
population      <- commandArgs()[6]
data_tag        <- commandArgs()[7]
run_date        <- commandArgs()[8]

### Set repo
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
if (!dir.exists(core_repo)) core_repo <- "<<<< FILEPATH REDACTED >>>>"
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")

setwd(core_repo)

## Load libraries and  MBG project functions.
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "<<<< FILEPATH REDACTED >>>>"))
mbg_setup(package_list = package_list, repos = core_repo)

## Record git status info
outputdir <- paste("<<<< FILEPATH REDACTED >>>>", indicator_group, indicator, 'output', run_date, '', sep='/')
record_git_status(core_repo, indic_repo, show_diff = T, file = paste0(outputdir, "git_at_launch.txt"))

# Load and check config
config <- set_up_config(repo            = indic_repo,
                        core_repo       = core_repo,
                        indicator       = "",
                        indicator_group = "",
                        config_name     = paste0('mbg/covariates/2_modeling/config_', prior),
                        covs_name       = 'mbg/covariates/2_modeling/cov_list')


message("Done with config")

## Create a few objects from the config file loaded above
if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text = Regions))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
if (length(summstats) == 1 & grepl(",", summstats)) summstats <- eval(parse(text=summstats))

## Make sure pc prior options are correct class
#if (is.character(p_above_value)) p_above_value <- as.numeric(p_above_value)
if (is.character(use_inla_nugget)) use_inla_nugget <- as.logical(use_inla_nugget)
if (is.character(use_inla_country_res)) use_inla_country_res <- as.logical(use_inla_country_res)

# Correct from range in earth radii to degrees, rough approximation (6371 km in radii, each degree is ~ 100km)
#if (use_s2_mesh == F & pc_prior == T) range_lower_tail <- range_lower_tail*637.1

## Make sure correct population is specified
pop_measure <- population

# Make sure p_above_value is correctly specified
if (is.character(p_above_value)) p_above_value <- as.numeric(p_above_value)

if (data_tag == "none") {
  withtag <- FALSE
  datatag <- ""
} else {
  withtag <- TRUE
  datatag <- data_tag
}

##Set cluster environment for fair cluster
Sys.setenv(SGE_ENV=SGE_ENV)

## Get cluster project-related arguments, default right now is geospatial queue 
#queue <- 'long.q'
stopifnot(queue %in% c('long.q', 'geospatial.q'))
if (queue == 'geospatial.q') {
  project <- 'proj_geo_nodes_hiv'
  use_geos_nodes <- T
}
if (queue == 'long.q') {
  project <- 'proj_geospatial_hiv'
  use_geos_nodes <- F
}

###############################################################################
## Make Holdouts
###############################################################################
if(makeholdouts) {

  #set seed?

  # load the full input data
  df <- load_input_data(indicator   = indicator,
                        removeyemen = TRUE,
                        years       = yearload,
                        withtag     = as.logical(withtag),
                        datatag     = datatag,
                        yl          = year_list,
                        use_share   = as.logical(use_share))

  # add in location information
  df <- merge_with_ihme_loc(df, shapefile_version = shapefile_version)

  if (holdout_strategy == "qt") {
    stratum_ho <- make_folds(data       = df,
                             n_folds    = as.numeric(n_ho_folds),
                             spat_strat = 'qt',
                             temp_strat = 'prop_comb',
                             strat_cols = 'region',
                             ts         = as.numeric(ho_ts),
                             mb         = as.numeric(ho_mb))

  } else if (holdout_strategy == "nids") {
    stratum_ho <- make_folds(data       = df,
                             n_folds    = as.numeric(n_ho_folds),
                             spte_strat = 'nids',
                             strat_cols = 'region',
                             ss_col     = ss_col,
                             yr_col     = yr_col)

  } else if (holdout_strategy == "poly") {
    stratum_ho <- make_folds(data         = df,
                             n_folds      = as.numeric(n_ho_folds),
                             spat_strat   = 'poly',
                             temp_strat   = 'prop_comb',
                             strat_cols   = 'region',
                             shape_ident  = 'gaul_code',
                             admin_shps   = "<<<< FILEPATH REDACTED >>>>/shapefiles/ad2_raster.grd",
                             admin_raster = "<<<< FILEPATH REDACTED >>>>/shapefiles/africa_ad2.shp",
                             mask_shape   = "<<<< FILEPATH REDACTED >>>>/shapefiles/africa_simple.shp",
                             mask_raster  = "<<<< FILEPATH REDACTED >>>>/shapefiles/ad0_raster.grd",
                             lat_col      = lat_col,
                             long_col     = long_col,
                             ss_col       = ss_col,
                             yr_col       = yr_col,
                             ts           = as.numeric(ho_ts),
                             mb           = as.numeric(ho_mb))
  } else {
    stop("holdout strategy not recognized")
  }
}

###############################################################################
## Launch Parallel Script
###############################################################################

### Make loopvars aka strata grid (format = regions, ages, holdouts)
if (makeholdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)



#Loop over them, save images and submit qsubs
for (i in 1:nrow(loopvars)) {

  message(paste(loopvars[i, 2], as.character(loopvars[i, 1]), loopvars[i, 3]))

  # set memory
  memory <- c("cssa-GNQ" = 100, "essa-DJI-SOM-SSD-YEM" = 190, "sssa" = 150, "wssa-CPV-GMB-MRT-STP" = 190)[as.character(loopvars[i, 1])]
  memory <- c("cssa-GNQ" = 100, "essa" = 190, "sssa" = 150, "wssa" = 190)[as.character(loopvars[i, 1])]
  if (is.na(memory)) memory <- 200

  # Restrict using country random effects for larger regions
  if (as.character(loopvars[i, 1]) %in% c("IND", "KHM+VNM")) use_country_res <- FALSE else use_country_res <- TRUE

  # For male circumcision do not use year fixed effects in region with only one year data
  if (indicator == "male_circumcision" & as.character(loopvars[i, 1]) == "KHM+VNM") fixed_effects <- ""

  # make a qsub string
  qsub <- make_qsub_share(age            = loopvars[i, 2],
                          reg            = as.character(loopvars[i, 1]),
                          holdout        = loopvars[i, 3],
                          test           = as.logical(test),
                          indic          = indicator,
                          saveimage      = TRUE,
                          memory         = memory,
                          cores          = 5,
                          proj           = project,
                          geo_nodes      = use_geos_nodes,
                          singularity    = '"<<<< FILEPATH REDACTED >>>>"',
                          singularity_opts = list(SET_OMP_THREADS=3, SET_MKL_THREADS=1),
                          queue          = queue,
                          run_time       = run_time)
  system(qsub)
}


## check to make sure models are done before continuing
waitformodelstofinish(lv = cbind(as.character(loopvars[,1]),loopvars[,3]),sleeptime=60)

# plot hyperparameters ##################################################
source(paste0(indic_repo, "<<<< FILEPATH REDACTED >>>>"))
try(plot_hyperparameters(indicator, "hiv", run_date, 0, 0))

###############################################################################
## Make summary metrics
###############################################################################

## Combine csv files
csvs <- list.files(outputdir, pattern = "input_data_(.*).csv", full.names = T)
csv_master <- rbindlist(lapply(csvs, fread))
csv_master[, V1 := NULL]
write.csv(csv_master, file=paste0(outputdir, '/input_data.csv'))

qsub <- paste('qsub -e', paste0(outputdir, 'errors'), '-o', paste0(outputdir, 'output'),
              '-P', project, paste0('-q ', queue), '-l fthread=10 -l m_mem_free=100G',
              "-l h_rt=0:05:00:00", '-N', paste0('pv_', run_date), '-l archive=TRUE',
              '-v sing_image=default -v SET_OMP_THREADS=1 -v SET_MKL_THREADS=1',
              paste0("<<<< FILEPATH REDACTED >>>>", '/lbd_core/mbg_central/share_scripts/shell_sing.sh'),
              paste0("<<<< FILEPATH REDACTED >>>>", Sys.info()['user'], '/lbd_hiv/mbg/functions/predictive_validity.r'),
              indicator, indicator_group, run_date, core_repo, shapefile_version)
system(qsub)

###############################################################################
## Post-Estimation
###############################################################################

## Save strata for Shiny to use in producing aggregated fit statistics
strata <- unique(as.character(loopvars[, 1]))
dir.create(paste0(outputdir, '/fit_stats'))
save(strata, file = paste0(outputdir, '/fit_stats/strata.RData'))

## Load GBD Estimates for this indicator which will be used in raking
age_group <- 24
sex_id <- 1
connector <- get_gbd_locs(rake_subnational = rake_subnational, reg = Regions, shapefile_version = shapefile_version)

# calculate the complement (can be added to the config file)
calc_complement = TRUE

source("<<<< FILEPATH REDACTED >>>>/current/r/get_outputs.R")
gbd <- get_outputs(topic = "cause",
                   version = "best", 
                   gbd_round_id = 5,
                   cause_id = 298,
                   measure_id = 5,
                   metric_id = 3,
                   age_group_id = age_group,
                   location_id = connector$location_id,
                   year_id = year_list)
gbd <- gbd[, list(name = location_id, year = year_id, mean = val)]

## Prepare for parallel post-estimation - save file with objs to re-load in child processes
prep_postest(indicator        = indicator,
             indicator_group  = indicator_group,
             run_date         = run_date,
             save_objs        = c("core_repo", "gbd", "year_list", 
                                  "summstats", "pop_measure", 
                                  "rake_subnational", "shapefile_version", 
                                  "config", "rake_transform", 
                                  "calc_complement", "pop_release",
                                  "p_above_value")) # can add gbd and rake_transform

#p_above_value should be added for male circumcision model runs

## Parallelized post-estimation over region
for (s in strata) {
  
  memory <- c("cssa-GNQ" = 450, "essa-DJI-SOM-SSD-YEM" = 650, "sssa" = 400, "wssa-CPV-GMB-MRT-STP" = 700)[s]
  memory <- c("cssa" = 450, "essa_sdn" = 650, "sssa" = 400, "wssa" = 700)[s]
  if (is.na(memory)) memory <- 600

  # make qsub string
  qsub <- make_qsub_postest(code           = "postest_frax_script",
                            stratum        = s,
                            log_location   = 'sharedir',
                            memory         = memory,
                            cores          = 5,
                            run_time       = '01:00:00:00',
                            proj           = project,
                            geo_nodes      = use_geos_nodes,
                            singularity    = "<<<< FILEPATH REDACTED >>>>",
                            modeling_shapefile_version = shapefile_version,
                            raking_shapefile_version = shapefile_version)
  system(qsub)
}


## check to make sure post-est done before continuing
waitformodelstofinish(lv = cbind(strata, 0), sleeptime = 60)

# Source accomodate counts function 
source("<<<< FILEPATH REDACTED >>>>/male_circumcision/accomodating_counts_functions.R")

## Combine post est stuff across regions and save needed outputs

post_load_combine_save(regions    = strata,
                       summstats  = setdiff(summstats, "p_above"),
                       raked      = c("unraked", "raked"),
                       rf_table   = TRUE,
                       run_summ   = TRUE,
                       indic      = indicator,
                       ig         = indicator_group)

# Combine p_above into a raster brick
if ("p_above" %in% summstats) {
  post_load_combine_save_counts(regions    = strata,
                                summstats  = "p_above",
                                raked      = "raked",
                                rf_table   = F,
                                run_summ   = F,
                                indic      = indicator,
                                ig         = indicator_group)
}

## Clean up / delete unnecessary files
clean_after_postest(indicator             = indicator,
                    indicator_group       = indicator_group,
                    run_date              = run_date,
                    strata                = strata,
                    delete_region_rasters = F)

combine_aggregation(rd = run_date,
                    indic = indicator,
                    ig = indicator_group,
                    ages = 0,
                    regions = strata,
                    holdouts = 0,
                    raked = c("unraked", "raked"), #c("unraked", "raked", "raked_c")
                    dir_to_search = NULL,
                    delete_region_files = F,
                    merge_hierarchy_list = F)

summarize_admins(ind = indicator,
                 ig = indicator_group,
                 summstats = setdiff(summstats, "p_above"),
                 raked    = c("unraked", "raked"), #c("unraked", "raked", "raked_c")
                 ad_levels = c(0, 1, 2),
                 file_addin = NULL)

# Create uncircumcised estimates for admin aggregations----------
summarize_admins_counts(ind = indicator,
                        ig = indicator_group,
                        summstats = setdiff(summstats, "p_above"),
                        raked    = "raked_c",
                        ad_levels = c(0, 1, 2),
                        file_addin = NULL,
                        calc_complement = TRUE)

# Make admin uncertainty plots over time used in Male circumcision paper ----------------------------------------------------
maps_path <- paste0("<<<< FILEPATH REDACTED >>>>", run_date)
 
# source(""<<<< FILEPATH REDACTED >>>>"")
make_admin_change_uncertainty_tables(run_date = run_date,
                                     indicator = "male_circumcision",
                                     indicator_group = "hiv",
                                     regions = Regions, #c('sssa','cssa-GNQ','wssa-CPV-GMB-MRT','essa-DJI-SOM-SSD-YEM'),
                                     outdir = maps_path,
                                     raked = T,
                                     years = c(2000, 2008))

make_admin_change_uncertainty_tables(run_date = run_date,
                                     indicator = "male_circumcision",
                                     indicator_group = "hiv",
                                     regions = Regions, # c('sssa','cssa-GNQ','wssa-CPV-GMB-MRT','essa-DJI-SOM-SSD-YEM'),
                                     outdir = maps_path,
                                     raked = T,
                                     years = c(2008, 2017))

## admin estimates of prevalence over 80 ----------------------------------------------------------------
country_map <-
  fread(paste0(maps_path, "/pred_derivatives/admin_summaries/", indicator, "_admin_2_raked_summary.csv")) %>%
  dplyr::select(ADM0_CODE, ADM0_NAME) %>%
  unique()

admin1_map <-
  fread(paste0(maps_path, "/pred_derivatives/admin_summaries/", indicator, "_admin_2_raked_summary.csv")) %>%
  dplyr::select(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME) %>%
  unique()

admin2_map <-
  fread(paste0(maps_path, "/pred_derivatives/admin_summaries/", indicator, "_admin_2_raked_summary.csv")) %>%
  dplyr::select(ADM0_CODE, ADM0_NAME, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME) %>%
  unique()

load(paste0(maps_path,"/male_circumcision_raked_admin_draws_eb_bin0_0.RData"))
library(tidyr)
#
# # Create coverage of 80% estimates
dir.create(paste0(maps_path, "/pct_80/"))
adm0_coverage <-
  admin_0 %>%
  filter(year == 2017) %>%
  dplyr::select(-year, -pop) %>%
  gather(draws, value, -ADM0_CODE) %>%
  group_by(ADM0_CODE) %>%
  dplyr::summarize(n_draws = n(),
            n_over_80 = sum(value > 0.8),
            pct_80 = 100*n_over_80 / n_draws) %>%
  ungroup() %>%
  left_join(country_map, by = "ADM0_CODE") %>%
  dplyr::select(ADM0_NAME, ADM0_CODE, pct_80) %>%
  arrange(ADM0_CODE) %>%
  data.table() %>%
  fwrite(paste0(maps_path, "/pct_80/admin_0.csv"))

adm1_coverage <-
  admin_1 %>%
  filter(year == 2017) %>%
  dplyr::select(-year, -pop) %>%
  gather(draws, value, -ADM1_CODE) %>%
  group_by(ADM1_CODE) %>%
  dplyr::summarize(n_draws = n(),
            n_over_80 = sum(value > 0.8),
            pct_80 = 100*n_over_80 / n_draws) %>%
  ungroup() %>%
  left_join(admin1_map, by = "ADM1_CODE") %>%
  dplyr::select(ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME, pct_80) %>%
  arrange(ADM0_CODE, ADM1_CODE) %>%
  data.table() %>%
  fwrite(paste0(maps_path, "/pct_80/admin_1.csv"))

adm2_coverage <-
  admin_2 %>%
  filter(year == 2017) %>%
  dplyr::select(-year, -pop) %>%
  gather(draws, value, -ADM2_CODE) %>%
  group_by(ADM2_CODE) %>%
  dplyr::summarize(n_draws = n(),
            n_over_80 = sum(value > 0.8),
            pct_80 = 100*n_over_80 / n_draws) %>%
  ungroup() %>%
  left_join(admin2_map, by = "ADM2_CODE") %>%
  dplyr::select(ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME, pct_80) %>%
  arrange(ADM0_CODE, ADM1_CODE, ADM2_CODE) %>%
  data.table() %>%
  fwrite(paste0(maps_path, "/pct_80/admin_2.csv"))

rm(admin_0)
rm(admin_1)
rm(admin_2)

###############################################################################
## Make maps and plots
###############################################################################

#Levels and differences maps
# source(paste0(indic_repo, "<<<< FILEPATH REDACTED >>>>"))
try(map_model_results(indicator, indicator_group, run_date, type = "mean", raked = 'raked', plot_by_year = F, limits_type = "quantile", lvl_limits = c(0.01, 0.99)))
try(map_model_results(indicator, indicator_group, run_date, type = "cirange", raked = 'raked', plot_by_year = F, include_diff = F, limits_type = "quantile", lvl_limits = c(0.01, 0.99)))
try(map_model_results(indicator, indicator_group, run_date, type = 'mean', limits_type = "quantile",
                             lvl_limits = c(0.01, 0.99), raked = 'raked_c', plot_by_year = F, geo_levels = "raster"))
try(map_model_results(indicator, indicator_group, run_date, type = 'mean', limits_type = "quantile",
                      lvl_limits = c(0.01, 0.99), raked = 'raked_c', plot_by_year = F, geo_levels = "admin1"))
try(map_model_results(indicator, indicator_group, run_date, type = 'mean', limits_type = "quantile",
                             lvl_limits = c(0.01, 0.99), raked = 'raked_c', plot_by_year = F, geo_levels = "admin2"))
#
# source(paste0(indic_repo, "<<<< FILEPATH REDACTED >>>>"))
try(data_and_estimates_maps_simplified(indicator, indicator_group, run_date, correct_anc = F))


################################################################################
# Prepare for visualization function
share_dir <- paste0("<<<< FILEPATH REDACTED >>>>", indicator_group, "/",
                    indicator, "/output/", run_date, "/")
in_dir <- paste0(share_dir, "<<<< FILEPATH REDACTED >>>>")
out_dir <- paste0(share_dir, "admin_plots/")

in_file_ad0 <- paste0(in_dir, indicator, "_admin_0_raked_summary.csv")
in_file_ad1 <- paste0(in_dir, indicator, "_admin_1_raked_summary.csv")
in_file_ad2 <- paste0(in_dir, indicator, "_admin_2_raked_summary.csv")

# Prepare inputs #######################################################

ad0_df <- fread(in_file_ad0)
ad1_df <- fread(in_file_ad1)
ad2_df <- fread(in_file_ad2)

# Drop Ma'tan al-Sarra if present
ad0_df <- subset(ad0_df, ADM0_CODE != 40762)
ad1_df <- subset(ad1_df, ADM0_CODE != 40762)
ad2_df <- subset(ad2_df, ADM0_CODE != 40762)

# Load input_data
admin_data <- 
  input_aggregate_admin(indicator = indicator,
                        indicator_group = indicator_group,
                        run_date = run_date,
                        regions = Regions)

ad0_data <- admin_data$ad0
ad1_data <- admin_data$ad1
ad2_data <- admin_data$ad2

# Make sure ZAF survey has conbined point and polygon for these surveys 
ad0_data <- 
  ad0_data %>% 
  filter(svy_id %in% c(228102, 313076)) %>% 
  group_by(svy_id, source, ADM0_NAME, ADM0_CODE, year) %>%
  summarize(point = 1, 
            outcome = weighted.mean(outcome, N),
            N = sum(N)) %>% 
  ungroup() %>% 
  rbind(filter(ad0_data, !svy_id %in% c(228102, 313076))) %>% 
  data.table()


ad1_data <- 
  ad1_data %>% 
  filter(svy_id %in% c(228102, 313076)) %>% 
  group_by(svy_id, source, ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, year) %>%
  summarize(point = 1, 
            outcome = weighted.mean(outcome, N),
            N = sum(N)) %>% 
  ungroup() %>% 
  rbind(filter(ad1_data, !svy_id %in% c(228102, 313076))) %>% 
  data.table()


ad2_data <- 
  ad2_data %>% 
  filter(svy_id %in% c(228102, 313076)) %>% 
  group_by(svy_id, source, ADM0_NAME, ADM0_CODE, ADM1_NAME, ADM1_CODE, ADM2_NAME, ADM2_CODE, year) %>%
  summarize(point = 1, 
            outcome = weighted.mean(outcome, N),
            N = sum(N)) %>% 
  ungroup() %>% 
  rbind(filter(ad2_data, !svy_id %in% c(228102, 313076))) %>% 
  data.table()


# Run the plotting code ################################################
source(paste0(core_repo, "mbg_central/visualization_functions.R"))
subnational_ts_plots(ad0_df = ad0_df,
                     ad1_df = ad1_df,
                     ad2_df = ad2_df,
                     ind_title = indicator,
                     out_dir = out_dir,
                     highisbad = F,
                     val_range = c(0,1),
                     ad0_map_regions = Regions, # c("cssa", "essa_sdn", "sssa", "wssa"),
                     ad0_map_region_titles = Regions, #c("Central Sub-Saharan Africa",
                                               # "Eastern Sub-Saharan Africa",
                                               # "Southern Sub-Saharan Africa",
                                               # "Western Sub-Saharan Africa"),
                     verbose = T,
                     plot_data = T,
                     ad0_data = ad0_data,
                     ad1_data = ad1_data,
                     ad2_data = ad2_data,
                     plot_levels = c("ad0", "ad1", "ad2"),
                     counts = F)


message('Done!')
