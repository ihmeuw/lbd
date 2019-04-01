###############################################################################
###############################################################################
## MBG Post-prediction script
##
## Purpose: This script takes modeled outputs and performs additional analyses
##          including aggregation, assessment of progress towards goals, and
##          analyses relating to model validation
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

## Load from qsub
load_from_parallelize() #indicator, vaccine, run_date

indicator_group <- 'vaccine'
sharedir <- sprintf('<<<< FILEPATH REDACTED >>>>/mbg/%s/%s',indicator_group, indicator)

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
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))
test <- as.logical(test)

gaul_list <- get_gaul_codes(Regions)

# Stop if individual countries (no need for post-est)
if(as.logical(individual_countries) == F) {

###############################################################################
## Post-Estimation
###############################################################################

## Make loopvars aka strata grid (format = regions, ages, holdouts)
if(makeholdouts) loopvars <- expand.grid(Regions, 0, 0:n_ho_folds) else loopvars <- expand.grid(Regions, 0, 0)

## Check to see if input_data objects have already been saved; save if not
message("Checking to see if input_data csv files present & creating if not...")
for (rr in Regions) {
  
  check_input_data(indicator = indicator,
                   indicator_group = indicator_group,
                   age = 0,
                   reg = rr,
                   run_date = run_date,
                   use_share = F) 
}

## Load holdouts
if(as.logical(makeholdouts)) {
  stratum_ho <- readRDS(paste0(sharedir,'/output/',run_date,"/stratum.rds"))
}

## Save strata
strata <- unique(as.character(loopvars[,1]))

# Skip the following if all you want is validation metrics
if (!as.logical(validation_metrics_only)) {

###############################################################################
## Aggregate to admin2, admin1, and national levels
###############################################################################
holdouts <- ifelse(makeholdouts, 0:n_ho_folds, 0)

log_dir <- paste0(sharedir, '/output/', run_date, '/logs/agg_logs/')
dir.create(log_dir, recursive = T)

# Define variables to loop over
agg_lv <- list(region = Regions,
               holdout = holdouts,
               age = 0,
               raked = c(TRUE, FALSE),
               overwrite = TRUE)

agg_qsub_output <- parallelize(script = "aggregate_results",
                               expand_vars = agg_lv,
                               save_objs = c("indicator", "indicator_group", "run_date", "pop_measure", "year_list"),
                               prefix = paste0("agg_", indicator),
                               slots = 10,
                               script_dir = paste0(core_repo, "/mbg_central/share_scripts/"), 
                               memory = 20,
                               geo_nodes = TRUE,
                               singularity = 'default')

monitor_jobs(agg_qsub_output)

combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages = 0, 
                    regions = strata,
                    holdouts = holdouts,
                    raked = c(T, F))

summarize_admins(summstats = c("mean", "lower", "upper", "cirange", "cfb"), 
                 ad_levels = c(0,1,2), 
                 raked = c(T,F))

summarize_admins(indicator, indicator_group,
                 summstats = c("p_above"),
                 raked = c(T,F),
                 ad_levels = c(0,1,2),
                 file_addin = "p_0.8_or_better",
                 value = 0.8,
                 equal_to = T)

summarize_admins(indicator, indicator_group,
                 summstats = c("p_below"),
                 raked = c(T,F),
                 ad_levels = c(0,1,2),
                 file_addin = "p_under_0.8",
                 value = 0.8,
                 equal_to = T)

###############################################################################
# Create AROC objects & do projections
###############################################################################

make_aroc(ind_gp = indicator_group, 
          ind = indicator,
          rd = run_date,
          matrix_pred_name = NULL,
          type = c("cell", "admin"),
          measure = "prevalence",
          year_list = year_list,
          uselogit = TRUE,
          raked = TRUE)

make_proj(ind_gp = indicator_group, 
          ind = indicator, 
          rd = run_date, 
          type = c("cell", "admin"),
          proj_years = c(2020, 2025, 2030),
          measure = "prevalence",
          skip_cols = NULL,
          year_list = year_list, 
          uselogit = TRUE,
          raked = TRUE)

###############################################################################
# Look at performance against goals
###############################################################################

# Define goals: start by initializing goal object
goals <- add_goal(target_year = 2030, 
                  target = 0.8,
                  target_type = "greater",
                  abs_rel = "absolute",
                  pred_type = c("cell", "admin"))

# Add goals to existing goal object by specifying goal_obj
goals <- add_goal(goal_obj = goals,
                  target_year = 2020, 
                  target = 0.8,
                  target_type = "greater",
                  abs_rel = "absolute",
                  pred_type = c("cell", "admin"))

goals <- add_goal(goal_obj = goals,
                  target_year = 2020, 
                  target = 0.9,
                  target_type = "greater",
                  abs_rel = "absolute",
                  pred_type = c("cell", "admin"))

compare_to_target(ind_gp = indicator_group, 
                  ind = indicator, 
                  rd = run_date, 
                  goal_obj = goals,
                  measure = "prevalence", 
                  year_list = year_list,
                  uselogit = TRUE,
                  raked = TRUE)

## Compare against all-admin2 > 80% target for GVAP
prob_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", indicator, 
                   "/output/", run_date, "/pred_derivatives/target_probs/")
dir.create(prob_dir, showWarnings = F, recursive = T)

gvap_table <- compare_all_admin_target(ind_gp = indicator_group, 
                                       ind = indicator, 
                                       rd = run_date,
                                       measure = "prevalence",
                                       target_year = 2020,
                                       target = 0.8,
                                       target_type = "greater",
                                       admin_level = 2,
                                       uselogit = T,
                                       proj = T)

write.csv(gvap_table, file = paste0(prob_dir, indicator, 
                                    "_all_admins_2020_absolute_greater_0.8_admin_target_probs.csv"))

## Compare against all-admin2 > 80% target for GVAP, but using 2015
gvap_2016  <- compare_all_admin_target(ind_gp = indicator_group, 
                                       ind = indicator, 
                                       rd = run_date,
                                       measure = "prevalence",
                                       target_year = max(year_list),
                                       target = 0.8,
                                       target_type = "greater",
                                       admin_level = 2,
                                       uselogit = T,
                                       proj = F) 

write.csv(gvap_2016, file = paste0(prob_dir, indicator, 
                                    "_all_admins_2016_absolute_greater_0.8_admin_target_probs.csv"))

## Compare against all-admin2 > 80% target for GVAP, but using 2000
gvap_2000  <- compare_all_admin_target(ind_gp = indicator_group, 
                                       ind = indicator,
                                       rd = run_date,
                                       measure = "prevalence",
                                       target_year = min(year_list),
                                       target = 0.8,
                                       target_type = "greater",
                                       admin_level = 2,
                                       uselogit = T,
                                       proj = F) 

write.csv(gvap_2000, file = paste0(prob_dir, indicator, 
                                    "_all_admins_2000_absolute_greater_0.8_admin_target_probs.csv"))

# Now look at the inverse - probability that all are *under* 80%
gvap_under <- compare_all_admin_target(ind_gp = indicator_group, 
                                       ind = indicator, 
                                       rd = run_date,
                                       measure = "prevalence",
                                       target_year = 2020,
                                       target = 0.8,
                                       target_type = "less",
                                       admin_level = 2,
                                       uselogit = T,
                                       proj = T)

write.csv(gvap_under, file = paste0(prob_dir, indicator, 
                                    "_all_admins_2020_absolute_less_0.8_admin_target_probs.csv"))

## Compare against all-admin2 > 80% target for GVAP, but using 2016
gvap_2016_under  <- compare_all_admin_target(ind_gp = indicator_group, 
                                             ind = indicator, 
                                             rd = run_date,
                                             measure = "prevalence",
                                             target_year = max(year_list),
                                             target = 0.8,
                                             target_type = "less",
                                             admin_level = 2,
                                             uselogit = T,
                                             proj = F) 

write.csv(gvap_2016_under, file = paste0(prob_dir, indicator, 
                                         "_all_admins_2016_absolute_less_0.8_admin_target_probs.csv"))

###############################################################################
## Look at (first year - last year) differences
###############################################################################
ad_summary_dir <- paste0("<<<< FILEPATH REDACTED >>>>/", indicator_group, "/", indicator, 
                         "/output/", run_date, "/pred_derivatives/admin_summaries/")
dir.create(ad_summary_dir, showWarnings = F, recursive = T)

for (ad in c(0,1,2)) {

  difference_df <- compare_admin_years(ind_gp = indicator_group,
                                       ind = indicator,
                                       rd = run_date,
                                       measure = "prevalence",
                                       start_year = min(year_list),
                                       end_year = max(year_list),
                                       admin_level = ad,
                                       uselogit = T,
                                       year_list = year_list,
                                       summstats = c("mean", "upper", "lower", "cirange"))

  write.csv(difference_df, file = paste0(ad_summary_dir, indicator, "_admin_", ad, "_raked_diff_",
                                    min(year_list), "-", max(year_list), ".csv"))

}


# ##############################################################################
# ## Make summary metrics
# ##############################################################################

# Combine csv files only if none present
csvs <- list.files(paste0(sharedir, '/output/', run_date, '/'), 
                   pattern = "input_data(.*).csv", 
                   full.names = T)

if (sum(grepl("input_data.csv", csvs)) == 0) {
  csv_master <- lapply(csvs, fread) %>% 
                rbindlist %>%
                subset(., select = names(.) != "V1")
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
                               year_col = "year")

## set out_dir
out_dir <- paste0(sharedir, "/output/", run_date, "/summary_metrics/")
dir.create(out_dir, recursive = T, showWarnings = F)

## set up titles
if (indicator == "dpt3_cov") plot_title <- "DPT3 Coverage"
if (indicator == "dpt1_cov") plot_title <- "DPT1 Coverage"
if (indicator == "dpt1_3_abs_dropout") plot_title <- "DPT1-3 Absolute Dropout"
if (indicator == "dpt1_3_rel_dropout") plot_title <- "DPT1-3 Relative Dropout"

## for admin0
draws.df <- fread(sprintf("<<<< FILEPATH REDACTED >>>>/%s/%s/output/%s/output_draws_data.csv",
                          indicator_group, indicator, run_date))

country.pvtable <- get_pv_table(d = draws.df,
                                indicator_group = indicator_group,
                                rd = run_date,
                                indicator=indicator,
                                coverage_probs = c(95),
                                aggregate_on='country',
                                draws = as.numeric(samples),
                                out.dir = out_dir,
                                plot_ci = TRUE, 
                                point_alpha = 0.5,
                                point_color = "black",
                                ci_color = "gray",
                                plot_title = plot_title)

# By region
country.pvtable.reg <- get_pv_table(d = draws.df,
                                    indicator_group = indicator_group,
                                    rd = run_date,
                                    indicator=indicator,
                                    result_agg_over = c("year", "oos", "region"),
                                    coverage_probs = c(95),
                                    aggregate_on='country',
                                    draws = as.numeric(samples),
                                    out.dir = out_dir,
                                    plot_ci = TRUE, 
                                    point_alpha = 0.5,
                                    point_color = "black",
                                    ci_color = "gray",
                                    plot_title = plot_title)

write.csv(country.pvtable,
          file = sprintf("<<<< FILEPATH REDACTED >>>>/%s/%s/output/%s/summary_metrics/country_metrics.csv",
                         indicator_group, indicator, run_date))

ad1.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            coverage_probs = c(95),
                            aggregate_on='ad1',
                            draws = as.numeric(samples),
                            out.dir = out_dir,
                            plot_ci = TRUE,
                            point_alpha = 0.2,
                            point_color = "black",
                            ci_color = "gray",
                            plot_title = plot_title)

# By region

ad1.pvtable.reg <- get_pv_table(d = draws.df,
                                indicator_group = indicator_group,
                                rd = run_date,
                                indicator=indicator,
                                result_agg_over = c("year", "oos", "region"),
                                coverage_probs = c(95),
                                aggregate_on='ad1',
                                draws = as.numeric(samples),
                                out.dir = out_dir,
                                plot_ci = TRUE,
                                point_alpha = 0.2,
                                point_color = "black",
                                ci_color = "gray",
                                plot_title = plot_title)

write.csv(ad1.pvtable,
          file = sprintf("<<<< FILEPATH REDACTED >>>>/%s/%s/output/%s/summary_metrics/ad1_metrics.csv",
                           indicator_group, indicator, run_date))

ad2.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            indicator=indicator,
                            coverage_probs = c(95),
                            aggregate_on='ad2',
                            draws = as.numeric(samples),
                            out.dir = out_dir,
                            plot_ci = TRUE, 
                            point_alpha = 0.1,
                            point_color = "black",
                            ci_color = "gray",
                            plot_title = plot_title)

# By region

ad2.pvtable.reg <- get_pv_table(d = draws.df,
                                indicator_group = indicator_group,
                                rd = run_date,
                                indicator=indicator,
                                result_agg_over = c("year", "oos", "region"),
                                coverage_probs = c(95),
                                aggregate_on='ad2',
                                draws = as.numeric(samples),
                                out.dir = out_dir,
                                plot_ci = TRUE, 
                                point_alpha = 0.1,
                                point_color = "black",
                                ci_color = "gray",
                                plot_title = plot_title)

write.csv(ad2.pvtable,
          file = sprintf("<<<< FILEPATH REDACTED >>>>/%s/%s/output/%s/summary_metrics/ad2_metrics.csv",
                           indicator_group, indicator, run_date))

# End individual country loop
}