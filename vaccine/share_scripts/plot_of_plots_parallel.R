###############################################################################
###############################################################################
## Create plots comparing model specifications (parallel / child script)
##
## Purpose: Create plots comparing OOS validation metrics for various 
##          combinations of stacking vs raw covariates with and without 
##          Gaussian process. Called by `plot_of_plots_master.R`
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

## Set up some required objects and directories

load_from_parallelize()
#   Varying: indicator, run_date
#   Static:indicator_group, samples

ind <- indicator
rd <- run_date

mod.dir <- sprintf('<<<< FILEPATH REDACTED >>>>/%s/%s/output/%s/',
                   indicator_group, indicator, run_date)
si.fig.dir <- paste0(mod.dir, 'si_figs/')
dir.create(si.fig.dir, showWarnings = F)

## For admin0
draws.df <- fread(sprintf("<<<< FILEPATH REDACTED >>>>/%s/%s/output/%s/output_draws_data.csv",
                          indicator_group, indicator, run_date))

samples <- length(grep('draw', colnames(draws.df)))

## admin 0
message("Admin 0")
ad0.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            plot=F,
                            indicator=indicator,
                            aggregate_on='country',
                            result_agg_over="oos",
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            out.dir = si.fig.dir,
                            save_csv = F)
write.csv(ad0.pvtable,
          file = sprintf("%s/ad0_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)

## admin 1
message("Admin 1")
ad1.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            plot=F,
                            indicator=indicator,
                            aggregate_on='ad1',
                            result_agg_over="oos",
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            out.dir = si.fig.dir,
                            save_csv = F)
write.csv(ad1.pvtable,
          file = sprintf("%s/ad1_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)

## admin 2
message("Admin 2")
ad2.pvtable <- get_pv_table(d = draws.df,
                            indicator_group = indicator_group,
                            rd = run_date,
                            plot=F,
                            indicator=indicator,
                            aggregate_on='ad2',
                            result_agg_over="oos",
                            draws = as.numeric(samples),
                            plot_ci = TRUE,
                            out.dir = si.fig.dir,
                            save_csv = F)

write.csv(ad2.pvtable,
          file = sprintf("%s/ad2_metrics.csv",
                         si.fig.dir),
          row.names = FALSE)