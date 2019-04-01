###############################################################################
###############################################################################
## Plot covariate importance
##
## Purpose: Produce plots of covariate influence from stacking and INLA
###############################################################################
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

## Script-specific code begins here #######################################################

## Options ---------------------------------------

run_date <- NULL # Change to model run date
indicator <- "dpt3_cov"
indicator_group <- "vaccine"
regions <- c("cssa", "essa", "name", "sssa", "wssa")

output_plot_dir <- paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', indicator, '/output/', run_date, "/cov_importance_plots/")
dir.create(output_plot_dir)

## Make covariate importance heat map 
if(length(list.files(paste0('<<<< FILEPATH REDACTED >>>>/', indicator_group, '/', indicator, '/output/', run_date), pattern = 'child_model_list')) != 0) {

  for (rr in regions) {

    message(reg)
    
    cov_importance <- get.cov.wts(rd = run_date,
                                  ind = indicator, 
                                  ind_gp = indicator_group, 
                                  reg = rr,
                                  age = 0,
                                  holdout = 0)

    cov_gg <- plot.cov.wts(rd = run_date,
                           ind = indicator, 
                           ind_gp = indicator_group, 
                           reg = rr,
                           age = 0,
                           holdout = 0,
                           plot.inla.col = F)

    png(file = paste0(output_plot_dir, indicator, "_rel_cov_imp_", reg, ".png"),
        height = 4,
        width = 8,
        units = "in", 
        res = 300)

    print(cov_gg)
    dev.off()
  }
  
}
