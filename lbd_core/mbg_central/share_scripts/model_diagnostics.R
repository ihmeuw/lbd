##########################################################################
# HEADER #################################################################
# Author: Jonathan Mosser (jmosser@uw.edu)
# Date: 2017-11-08
# Project: Geospatial - all
# Purpose: Model diagnostics (stacker plots, etc) for MBG
# Details:
##########################################################################
##########################################################################

##########################################################################
## SETUP
##########################################################################

# clear memory
rm(list=ls())

## Load from qsub
indicator       <- as.character(commandArgs()[4]); message(paste0("indicator: ", indicator))
indicator_group <- as.character(commandArgs()[5]); message(paste0("indicator_group: ", indicator_group))
run_date        <- as.character(commandArgs()[6]); message(paste0("run_date: ", run_date))
n_cores         <- as.numeric(commandArgs()[7]); message(paste0("n_cores: ", n_cores))
core_repo       <- as.character(commandArgs()[8]); message(paste0("core_repo: ", core_repo))

## Set repo location and indicator group
user            <- Sys.info()['user']

## Folders & drive locations
commondir       <- '/share/geospatial/mbg/common_inputs'
sharedir        <- paste0("/share/geospatial/mbg/", indicator_group, "/", indicator,
                         "/output/", run_date, "/")
out_dir         <- paste0(sharedir, "shiny_diagnostics/")

## Load libraries and  MBG project functions.
package_list    <- c(t(read.csv(sprintf('%s/package_list.csv', commondir), header = FALSE)))
package_list    <- c(package_list, "gridExtra")

# Load MBG packages and functions
message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

## Read config file (from sharedir) and save all parameters in memory
config <- load_config(repo            = core_repo, ## since this is post_est_only, the repo doesn't matter. config is loaded from run_date dir
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,
                      run_date        = run_date)

## Create a few objects from options above

# Check in case region_list is used instead of Regions
if (exists("region_list")) { 
  Regions <-  unlist(strsplit(region_list, " \\+ "))
}

# If all else fails, try to extract from the files in the directory
if (!exists("Regions")) {
  regions <- get_output_regions(in_dir = paste0("/share/geospatial/mbg/",indicator_group,"/",indicator,"/output/",run_date,"/"))
}

if (class(Regions) == "character" & length(Regions) == 1) Regions <- eval(parse(text=Regions))
if (class(year_list) == "character") year_list <- eval(parse(text=year_list))

##########################################################################
## MAIN CODE
##########################################################################

# Make sure that out_dir exists & create if not
dir.create(out_dir, recursive = T, showWarnings = F)

##########################################################################
## Plot model parameters #################################################

# Load model results
if (!file.exists(paste0(sharedir, indicator, "_model_results_table.csv"))) {
  clean_model_results_table() # create a model results table
}
model_results <- read.csv(paste0(sharedir, indicator, "_model_results_table.csv"),
                          stringsAsFactors = F) %>% as.data.table
names(model_results) <- c("region", "age", "parameter", "q0.025", "q0.5", "q0.975")

# Grab list of stackers
stackers <- model_results$parameter %>% unique
stackers <- stackers[!(stackers %in% c("int", "Nominal Range", "Nominal Variance",
                                       "GPRandom rho for time", "Precision for IID.ID",
                                       "Precision for CTRY.ID"))]

# Plot stacker betas -----------------------------------------------------
#   Output of this part: plot of stacker betas & 95% CIs
#   Saved as [out_dir]/stacker_betas.png

gg_betas <- plot_stacker_betas(model_results = model_results,
                               stackers = stackers,
                               xaxis = "stacker")

png(filename = paste0(out_dir, "stacker_betas.png"),
    width = 10,
    height = 5,
    units = "in",
    type = "cairo-png",
    res = 200,
    pointsize = 10)
print(gg_betas)
dev.off()

gg_betas <- plot_stacker_betas(model_results = model_results,
                               stackers = stackers,
                               xaxis = "region")

png(filename = paste0(out_dir, "stacker_betas_by_region.png"),
    width = 10, 
    height = 5, 
    units = "in",
    type = "cairo-png",
    res = 200,
    pointsize = 10)
print(gg_betas)
dev.off()

# Grab list of other parameters
other_params <- unique(model_results$parameter[!(model_results$parameter %in% c(stackers))])

# Plot other parameters --------------------------------------------------
#   Output of this part: plot of rho, int, variance, etc. from INLA
#   Saved as [out_dir]/other_inla_params.png

gg_other_params <- plot_other_params(model_results = model_results,
                                     other_params = other_params)

png(filename = paste0(out_dir, "other_inla_params.png"),
    width = 10,
    height = 5,
    units = "in",
    type = "cairo-png",
    res = 200,
    pointsize = 10)
print(gg_other_params)
dev.off()

##########################################################################
## Plot maps of stackers #################################################

# Plot stackers vs true data for each region -----------------------------
#   Output of this part: maps of each stacker along with the mean raster
#   (unraked) and a plot of the data (with weight = alpha for resampled
#   polygons --> pseudoclusters), one for each year, for each region &
#   country.
#
#   Saved as [out_dir]/stacker_maps/[ctry]/stacker_map_[ctry]_YYYY.png
#   (or equivalently [reg] instead of [ctry] where applicable)

# Set multithreading to serial for `mclapply()`:
set_serial_threads()

mclapply(Regions, function(r) {
  message(paste0("Making stacker maps for region: ", r))
  plot_stackers(reg = r, highisbad = F, individual_countries = T, shapefile_version = modeling_shapefile_version)
}, mc.cores = n_cores)

# Return to multithreading (if any):
set_original_threads()

##########################################################################
## Plot diagnostics for stackers

# Get a list of child models for each region

message("Loading child models for each region")
child_model_list <- lapply(Regions, function(reg) {
  child_model_file <- paste0(sharedir, "child_model_list_", reg, "_0.RData")
  load(child_model_file, verbose = F)
  return(child_models)
})

names(child_model_list) <- Regions

# This list is now structured such that child_model_list[["cssa"]][["gam"]]
# will return the appropriate gam child object for the cssa region

# Plot GAM models ------------------------------------------------------
#   Output of this part: Plots of the component smooth functions that make
#   up the fitted gam model from stacking for each region's model.
#   Uses mgcv::plot.gam(). Limited to 12 plots per page, and will produce
#   multiple pages if needed to ensure that all covariates are graphed.
#
#   Saved as [out_dir]/gam/[reg]/gam_plot_[reg]_page_[i].png
#   where i is {1,2...n} pages needed to fit all covariates

if ("gam" %in% stackers) {
  plot_gam_models(child_model_list = child_model_list,
                  regions = Regions,
                  o_dir = out_dir)
}