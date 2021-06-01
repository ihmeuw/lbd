################################################################################
## Purpose: This is the script that takes in the data objects for used to run a
##          draw of an admin 2 model and preforms the IMIS fitting and results
##          simulation.  It does not do data prep at all which is the main
##          difference between it and the main_geo script which does the data
##          prep for EPP runs.
## Run instructions: This gets launched from the main_geo.R script as part of
##                   an array job.  Each draw of the higher level model then
##                   launches an array that is 1 job for each admin2 in the
##                   that falls into the higher level geography.  This script
##                   does not require archive access.
################################################################################
## LBD Base Set Up
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>/mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

### GBD addititions Setup


windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0(indic_repo, "epp/2_modeling/epp_geo/")

## Packages
library(data.table); library(mvtnorm); library(survey)
new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>", "/r_packages/")
dir.create(new_pkg_lib)
.libPaths(new_pkg_lib)
test_pkg_list <- c('slackr', "assertable")
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(slackr)
library(assertable)
## Arguments

run.name <- as.character(commandArgs()[4])
message(paste0("the run name is ",run.name))
loc <- as.character(commandArgs()[5])
message(paste0("the loc is ",loc))
proj.end <- as.integer(commandArgs()[6])
message(paste0("the proj.end is ",proj.end))
no.anc <- as.logical(as.character(commandArgs()[7]))
message(paste0("the no.anc setting is ",no.anc))
id <- as.integer(Sys.getenv("SGE_TASK_ID"))
message(paste0("the ad2 id is ", id))
ad2_list <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
subnat <- as.character(ad2_list[[id]])
message(paste0("the subnat setting is ",subnat))
i <- as.integer(commandArgs()[8])
message(paste0("the mort draw id is ", i))
### read in the config file

config <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))

### Arguments
cluster.project <- as.character(config$Value[which(config$Setting == "cluster.project")])
cluster.queue <- as.character(config$Value[which(config$Setting == "cluster.queue")])
start.year <- as.numeric(as.character(config$Value[which(config$Setting == "start.year")]))
message(paste0("the demo start year is ", start.year))
stop.year <- as.numeric(as.character(config$Value[which(config$Setting == "stop.year")]))
message(paste0("the estimation stop year is ", stop.year))
trans.params <- as.logical(config$Value[which(config$Setting == "trans.params")])
message(paste0("the trans.params setting is ", trans.params))
stop.collapse <- as.logical(config$Value[which(config$Setting == "stop.collapse")])
message(paste0("the stop.collapse setting is ", stop.collapse))
gbd.pop <- as.logical(config$Value[which(config$Setting == "gbd.pop")])
message(paste0("the gbd.pop setting is ", gbd.pop))
if (no.anc == FALSE) {
  anc.prior <- TRUE
} else {anc.prior <- FALSE}
message(paste0("the anc.prior setting is ", anc.prior))
art.sub <- as.logical(config$Value[which(config$Setting == "art.sub")])
message(paste0("the art.sub setting is ", art.sub))
eq.prior <- as.logical(config$Value[which(config$Setting == "eq.prior")])
message(paste0("the ep.prior setting is ", eq.prior))
anc.backcast <- as.logical(config$Value[which(config$Setting == "anc.backcast")])
message(paste0("the anc.backcast setting is ", anc.backcast))
num.knots <- as.numeric(as.character(config$Value[which(config$Setting == "num.knots")]))
message(paste0("there will be ", num.knots, " knots in the r-spline model"))
run.admin2 <- as.logical(config$Value[which(config$Setting == "run.admin2")])
message(paste0("the run.admin2 setting is ", run.admin2))
shapefile_version <- (config$Value[which(config$Setting == "shapefile_version")])
message(paste0("the shapefile version is ", shapefile_version))
prev_date <- (config$Value[which(config$Setting == "prev_date")])
message(paste0("the prev_date version is ", prev_date))
migration_run  <- (config$Value[which(config$Setting == "migration_run")])
message(paste0("the migration_run version is ", migration_run))
migration_type <- (config$Value[which(config$Setting == "migration_type")])
message(paste0("the migration_type version is ", migration_type))
popadjust <- as.logical(config$Value[which(config$Setting == "popadjust")])
message(paste0("the popadjust setting is ", popadjust))
debug <- as.logical(config$Value[which(config$Setting == "debug")])
message(paste0("the debug setting is ", debug))
gbd_mort_run <- as.character(config$Value[which(config$Setting == "gbd_mort_run")])
eppmod_setting <- as.character(config$Value[which(config$Setting == "eppmod_setting")])
likelihood_type <- as.character(config$Value[which(config$Setting == "likelihood")])
prev_type <- as.character(config$Value[which(config$Setting == "prev_transform")])

### Paths
input.dir <- paste0()
out.dir <- paste0("<<<< FILEPATH REDACTED >>>>")
out.path <- paste0(out.dir, "/results", i, ".RData")
pdf.path <- paste0(out.dir, "/test_results", i, ".pdf")

### Functions
## GBD
source(paste0(code.dir,"gbd/prep_data_geo.R"))
source(paste0(code.dir,"gbd/prep_output_geo.R"))
source(paste0(code.dir,"gbd/data_sub_geo.R"))
source(paste0(code.dir,"gbd/plot_fit_geo.R"))
source(paste0(code.dir,"gbd/ind_data_prep.R"))
source(paste0(code.dir,"gbd/data_sub_subnat_geo.R"))
source(paste0(code.dir,"gbd/prep_data_admin2_geo.R"))


## EPP
source(paste0(code.dir,"R/epp_geo.R"))
source(paste0(code.dir,"R/fit-model_geo.R"))
source(paste0(code.dir,"R/generics.R"))
source(paste0(code.dir,"R/IMIS_geo.R"))
source(paste0(code.dir,"R/likelihood_geo.R"))
source(paste0(code.dir,"R/read-epp-files_geo.R"))

## Model in C
if (trans.params) {
  # Load C version
  dyn.load(paste0(code.dir, "src/fnEPP_popadjust_mort_v2", .Platform$dynlib.ext))  # Load C version with time series for transition parameters
} else {
  dyn.load(paste0(code.dir, "src/epp", .Platform$dynlib.ext))
}

### Tables
loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))


dt <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))

## Fit national model
fit <- list()
for (subpop in names(dt)) {
  print(subpop)
  no.anc <<- no.anc

  if (no.anc == TRUE) {attr(dt[[subpop]], "eppd")$anc.used <- FALSE}

  attr(dt[[subpop]], "eppfp")$artelig.idx.ts <- as.integer(attr(dt[[subpop]], "eppfp")$artelig.idx.ts)


  if (anc.prior) {
    set.anc.prior(loc, subpop)
  }


  if (length(names(dt)) == 1) {
    fit[[subpop]] <- fitmod(dt[[subpop]], equil.rprior = eq.prior, B0 = 2e5, B = 1e3, B.re = 1e3, number_k = 500, D = 4, opt_iter = 3, no.anc = no.anc, likelihood_type = likelihood_type, prev_type = prev_type)
  } else {
    fit[[subpop]] <- fitmod(dt[[subpop]], equil.rprior = eq.prior, B0 = 2e5, B = 1e3, B.re = 1e3, no.anc = no.anc, likelihood_type = likelihood_type, prev_type = prev_type)
  }
}

## Prepare output
result <- c()
for (subpop in names(fit)) {
result[[subpop]] <- simfit.gbd(fit = fit[[subpop]], random_walk = F, no.anc = no.anc)
}

save(result, file = paste0("<<<< FILEPATH REDACTED >>>>"))

## write prevalence and incidence draws
years <- unique(floor(result[[1]]$fp$proj.steps))
nat.data <- nat.draws(result)
var_names <- sapply(1:ncol(nat.data$prev), function(a) {paste0('draw',a)})
out_data <- lapply(nat.data, data.frame)
for (n in c('prev', 'incid', 'art', 'art_num', 'pop')) {
  names(out_data[[n]]) <- var_names
  out_data[[n]]$year <- years
  col_idx <- grep("year", names(out_data[[n]]))
  out_data[[n]] <- out_data[[n]][, c(col_idx, (1:ncol(out_data[[n]]))[-col_idx])]
  write.csv(out_data[[n]], paste0("<<<< FILEPATH REDACTED >>>>"), row.names = F)
}


## Plot results
plot.fit(result, pdf.path, nat.data, no.anc = no.anc)


#if this isn't a debugging run delete the input object for all but the first draw.
if (debug == FALSE & i != 1) {
  file.remove(paste0("<<<< FILEPATH REDACTED >>>>"))
}

### END
