# Project: Geospatial - all
# Purpose: Model diagnostics (stacker plots, etc) for MBG

##########################################################################
## SETUP
##########################################################################

# clear memory
rm(list=ls())

## Load from qsub
indicator       <-  as.character(commandArgs()[3]); message(paste0("indicator: ", indicator))
indicator_group <-  as.character(commandArgs()[4]); message(paste0("indicator_group: ", indicator_group))
run_date        <-  as.character(commandArgs()[5]); message(paste0("run_date: ", run_date))
n_cores         <-  as.numeric(commandArgs()[6]); message(paste0("n_cores: ", n_cores))

## Set repo location and indicator group
user               <- Sys.info()['user']
repo               <- sprintf('<<<< FILEPATH REDACTED >>>>>',user)

## Folders & drive locations
root           <- ifelse(Sys.info()[1]=='Windows', '<<<< FILEPATH REDACTED >>>>>', '<<<< FILEPATH REDACTED >>>>>')
package_lib    <- ifelse(grepl('geos', Sys.info()[4]),
                         paste0(root,'<<<< FILEPATH REDACTED >>>>>'),
                         paste0('<<<< FILEPATH REDACTED >>>>>'))
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>>')
sharedir       <- paste0("<<<< FILEPATH REDACTED >>>>>", indicator_group, "/", indicator, 
                         "/output/", run_date, "/") 
out_dir        <- paste0(sharedir, "shiny_diagnostics/")

## Load libraries and  MBG project functions.
.libPaths(package_lib)
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir),header=FALSE)))
package_list <- c(package_list, "gridExtra")

if(Sys.info()[1] == 'Windows'){
  stop('STOP! you will overwrite these packages if you run from windows\n
        STOP! also, lots of this functions wont work so get on the cluster!')
} else {
  for(package in package_list)
    require(package, lib.loc = package_lib, character.only=TRUE)
  for(funk in list.files(recursive=TRUE,pattern='functions')){
    if(length(grep('central',funk))!=0){
      message(paste0('loading ',funk))
      source(funk)
    }
  }
}

# Source some functions
setwd(repo)
source('mbg_central/prep_functions.R')
source('mbg_central/diagnostic_functions.R')

## Read config file (from sharedir) and save all parameters in memory
config <- load_config(repo            = repo,
                      indicator_group = indicator_group,
                      indicator       = indicator,
                      post_est_only   = TRUE,           
                      run_date        = run_date)

## Create a few objects from options above
if (exists("region_list")) {
  # Check in case region_list is used instead of Regions
  Regions <-  unlist(strsplit(region_list, " \\+ "))
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
                               stackers = stackers)

png(filename = paste0(out_dir, "stacker_betas.png"),
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


print(gg_other_params)

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

mclapply(Regions, function(r) {
  message(paste0("Making stacker maps for region: ", r))
  plot_stackers(reg = r, highisbad = F, individual_countries = T)
}, mc.cores = n_cores)

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

plot_gam_models(child_model_list = child_model_list,
                regions = Regions,
                o_dir = out_dir)
