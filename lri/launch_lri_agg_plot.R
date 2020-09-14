#####################################################################

# LRI aggregation and diagnostic plots
# 1) save strata
# 2) launch model diagnostics shiny
# 3) aggregate results
# 4) produce csvs
# 5) make line plots

####################################################################

# Setup -------------------------------------------------------------
#libraries
library(dplyr)
library(data.table)

# pull arguments from command_args
regions <- unlist(strsplit(commandArgs()[4], '~'))
holdout <- commandArgs()[5]
run_date <- commandArgs()[6]
age <- commandArgs()[7]

indicator <- 'has_lri'
indicator_group <- 'lri'

# load model image history
reg_ref <- regions[1]
pathaddin <- paste0('_bin',age,'_',reg_ref,'_',holdout)
load('<<<< FILEPATH REDACTED >>>>')

#custom functions
source('<<<< FILEPATH REDACTED >>>>/lbd_core_custom/misc_functions.R')
source('<<<< FILEPATH REDACTED >>>>/lri/stacker_admin_time_series.R')

#--------------------------------------------------------------------
# (1) save strata
#--------------------------------------------------------------------

strata <- unique(as.character(loopvars[,1]))
dir.create(paste0(sharedir, '/fit_stats'))
save(strata, file = paste0(sharedir, '/fit_stats/strata.RData'))

#--------------------------------------------------------------------
# (2) launch model diagnostics shiny
#--------------------------------------------------------------------
make_model_diagnostics(indic      = indicator,
                       ig         = indicator_group,
                       rd         = run_date,
                       geo_nodes  = TRUE,
                       cores      = 10,
                       corerepo = core_repo)

#--------------------------------------------------------------------
# (3) aggregate results
#--------------------------------------------------------------------

#aggregate raked
submit_aggregation_script (indicator       = indicator,
                           indicator_group = indicator_group,
                           run_date        = run_date,
                           raked           = T,
                           pop_measure     = pop_measure,
                           overwrite       = T,
                           ages            = 0, # Note: can take vector of ages
                           holdouts        = 0,
                           regions         = strata,
                           corerepo        = '<<<< FILEPATH REDACTED >>>>',
                           log_dir         = paste0(sharedir, "/output/", run_date),
                           geo_nodes       = T,
                           singularity     = "default",
                           slots           = 10,
                           modeling_shapefile_version = modeling_shapefile_version,
                           raking_shapefile_version = raking_shapefile_version,
                           run_time        = '10:00:00',
                           memory          = 100,
                           measures        = agg_measures)


#aggregate unraked
submit_aggregation_script (indicator       = indicator,
                           indicator_group = indicator_group,
                           run_date        = run_date,
                           raked           = F,
                           pop_measure     = pop_measure,
                           overwrite       = T,
                           ages            = 0, # Note: can take vector of ages
                           holdouts        = 0,
                           regions         = strata,
                           corerepo        = '<<<< FILEPATH REDACTED >>>>',
                           log_dir         = paste0(sharedir, "/output/", run_date),
                           geo_nodes       = T,
                           singularity     = "default",
                           slots           = 10,
                           modeling_shapefile_version = modeling_shapefile_version,
                           raking_shapefile_version = raking_shapefile_version,
                           run_time        = '10:00:00',
                           memory          = 100,
                           measures        = 'prevalence')

#check that all raked files exist
waitforaggregation(rd = run_date, indic = indicator, ig = indicator_group,
                   ages     = 0,
                   regions  = strata,
                   holdouts = 0,
                   raked    = TRUE,
                   measure  = agg_measures)

#check that all unraked files exist
waitforaggregation(rd = run_date, indic = indicator, ig = indicator_group,
                   ages     = 0,
                   regions  = strata,
                   holdouts = 0,
                   raked    = FALSE,
                   measure  = 'prevalence')

#--------------------------------------------------------------------
# (4) make csvs
#--------------------------------------------------------------------

combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages     = 0,
                    Regions  = strata,
                    holdouts = 0,
                    raked    = F,
                    measures = 'prevalence',
                    check_for_dupes = F)

summarize_admins(summstats = c("mean", "upper", "lower", "cirange"),
                 ad_levels = c(0,1,2),
                 raked     = F,
                 measures  = 'prevalence')

#raked combine and summarize
combine_aggregation(rd = run_date, indic = indicator, ig = indicator_group,
                    ages     = 0,
                    Regions  = strata,
                    holdouts = 0,
                    raked    = T,
                    measures = agg_measures,
                    check_for_dupes = F)

summarize_admins(summstats = c("mean", "upper", "lower", "cirange"),
                 ad_levels = c(0,1,2),
                 raked     = T,
                 measures  = agg_measures)


#--------------------------------------------------------------------
# (5) make line plots
#--------------------------------------------------------------------

plot_stackers_by_adm01(indicator = indicator,
                       indicator_group = indicator_group,
                       run_date = run_date,
                       regions = strata,
                       draws = T,
                       raked = T,
                       credible_interval = 0.95,
                       N_breaks = c(0, 10, 50, 100, 500, 1000, 2000, 4000),
                       admin_data = NULL,
                       data_tag = 'has_lri_stackerplot',
                       measure <- 'prevalence',
                       rm_yemen <- T)
