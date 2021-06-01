################################################################################
## Purpose: This is the main EPP model control script.  It creates the data
##          objects needed to run EPP both at the high level and at the admin 2
##          level.  It then spins up separate jobs to run the actual IMIS
##          simulations for the admin 2 models and finishes the high level model.
##          It gets launched once per draw per gbd location that we are modeling.
##          So if we are doing 1000 draws for 5 countries and 9 provinces in ZAF
##          this script would get launched 14,000 times.
## Run instructions: This is launched as part of array jobs from the launch EPP
##                   script.  This particular process accesses data on the J drive
##                   and thus requires archive access.  The child processes do not.
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
i <- as.integer(Sys.getenv("SGE_TASK_ID"))
message(paste0("the mort draw id is ", i))
ds <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
mort_draw <- ds[[i]]

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
shapefile_version <- as.character(config$Value[which(config$Setting == "shapefile_version")])
message(paste0("the shapefile version is ", shapefile_version))
prev_date <- as.character(config$Value[which(config$Setting == "prev_date")])
message(paste0("the prev_date version is ", prev_date))
migration_run  <- as.character(config$Value[which(config$Setting == "migration_run")])
message(paste0("the migration_run version is ", migration_run))
migration_type <- as.character(config$Value[which(config$Setting == "migration_type")])
message(paste0("the migration_type version is ", migration_type))
gbd.mort <- (config$Value[which(config$Setting == "gbd.mort")])
message(paste0("the gbd.mort version is ", gbd.mort))
mig_set <- (config$Value[which(config$Setting == "mig")])
message(paste0("the mig setting to use gbd migration is ", mig_set))
art_type <- (config$Value[which(config$Setting == "art_type")])
message(paste0("the admin2 art_type setting is ", art_type))
popadjust <- as.logical(config$Value[which(config$Setting == "popadjust")])
message(paste0("the popadjust setting is ", popadjust))
debug <- as.logical(config$Value[which(config$Setting == "debug")])
message(paste0("the debug setting is ", debug))
gbd_mort_run <- as.character(config$Value[which(config$Setting == "gbd_mort_run")])
eppmod_setting <- as.character(config$Value[which(config$Setting == "eppmod_setting")])
likelihood_type <- as.character(config$Value[which(config$Setting == "likelihood")])
prev_type <- as.character(config$Value[which(config$Setting == "prev_transform")])
mean_test <- as.logical(config$Value[which(config$Setting == "mean_test")])
tirr <- as.logical(config$Value[which(config$Setting == "tirr")])



### get randomized mort

### GBD input directories
gbd_mig_run <- as.character(config$Value[which(config$Setting == "gbd_mig_run")])
nat_mig_path <- paste0("<<<< FILEPATH REDACTED >>>>")
art.dir <- as.character(config$Value[which(config$Setting == "art.dir")])
ad0_prev_path <- paste0("<<<< FILEPATH REDACTED >>>>")
mortnoart_path <- as.character(config$Value[which(config$Setting == "mortnoart_path")])
mortart_path <- as.character(config$Value[which(config$Setting == "mortart_path")])
age_irr_path <- as.character(config$Value[which(config$Setting == "age_irr_path")])
sex_irr_path <- as.character(config$Value[which(config$Setting == "sex_irr_path")])
progdata_path <- as.character(config$Value[which(config$Setting == "progdata_path")])




### Paths
input.dir <- paste0()
out.dir <- paste0('<<<< FILEPATH REDACTED >>>>', run.name, "/", loc)
dir.create(out.dir)
out.path <- paste0(out.dir, "/results", i, ".RData")
pdf.path <- paste0(out.dir, "/test_results", i, ".pdf")

### Functions
## GBD
source(paste0(code.dir,"gbd/prep_data_geo.R"))#
source(paste0(code.dir,"gbd/prep_output_geo.R"))# ish
source(paste0(code.dir,"gbd/data_sub_geo.R"))#
source(paste0(code.dir,"gbd/plot_fit_geo.R"))#
source(paste0(code.dir,"gbd/ind_data_prep.R"))
source(paste0(code.dir,"gbd/data_sub_subnat_geo.R"))#
source(paste0(code.dir,"gbd/prep_data_admin2_geo.R"))#


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

library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")

### Tables
loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

loc.table$ad2_time <- 60


### Code
## Prep data and collapse location subpopulations
dt <- prep_epp_data(loc, proj.end = (stop.year + 0.5), stop_collapse = stop.collapse, gbd_pop = gbd.pop, art_sub = art.sub, num_knots = num.knots, no.anc = no.anc, run.name = run.name, gbd.mort = gbd.mort, mig = mig_set, shapefile_version = shapefile_version, i = i, gbd_mort_run = gbd_mort_run, mort_draw = mort_draw)
if (popadjust == TRUE) {attributes(dt[[1]])$eppfp$popadjust <- as.integer(1)} else {attributes(dt[[1]])$eppfp$popadjust <- as.integer(0)}
attributes(dt[[1]])$eppfp$likelihood_type <- likelihood_type
if (eppmod_setting == "rhybrid") {
  fp <- attributes(dt[[1]])$eppfp
  fp <- prepare_rhybrid(fp, tsEpidemicStart = fp$tsEpidemicStart, rw_start = 2003, rw_trans = 5, rw_dk = 5)
  attributes(dt[[1]])$eppfp <- fp
  attributes(dt[[1]])$eppfp$eppmod <- eppmod_setting
}

## Substitute IHME data
# Prevalence surveys
if ((!stop.collapse & length(dt) == 1) | grepl("IND_", loc)) {
  print("Substituting prevalence surveys")
  dt <- sub.prev(loc, dt, likelihood_type, prev_date, prev_type)
}

## ANC Back cast data
if (anc.backcast) {
  dt <- sub.anc(loc, dt)
}

# Transition parameters
if (trans.params) {
  dt <- extend.trans.params(dt, start.year, stop.year)
  dt <- sub.off.art(dt, loc, mort_draw, mean_test = mean_test, tirr = tirr)
  dt <- sub.on.art(dt, loc, mort_draw, mean_test = mean_test, tirr = tirr)
  dt <- sub.cd4.prog(dt, loc, mort_draw, mean_test = mean_test, tirr = tirr)
}




## Prep Subnational data
if (run.admin2) {
  admins <- read.dbf(get_admin_shapefile(admin_level = 2, suffix = ".dbf", version = shapefile_version))
  ad0_code <- get_adm0_codes(loc)
  admins <- admins[which(admins$ADM0_CODE == ad0_code), ]
  ad2_list <- admins$ADM2_CODE
  subs_dir <- paste0("<<<< FILEPATH REDACTED >>>>",run.name,"/", loc, "/admin2_models/")
  dir.create(subs_dir)
  if (i == 1 ) {saveRDS(ad2_list, paste0(subs_dir,"ad2_list.rds"))}
  for (subnat in ad2_list) {

    message(paste0("this is admin 2 ", subnat))
    subnat_dir <- paste0(subs_dir,subnat,"/")
    dir.create(subnat_dir)
    dt_subnat <- 1

    try(dt_subnat <- prep_epp_data_adm2(loc = loc, proj.end = 2019.5, subnat = subnat, run.name = run.name, i = i, stop_collapse = FALSE, gbd_pop = TRUE, art_sub = TRUE, num_knots = 7, no.anc = no.anc, prev_date = prev_date, migration_run = migration_run, migration_type = migration_type, stop.year = stop.year, gbd.mort = gbd.mort, art_type = art_type, shapefile_version = shapefile_version, gbd_mort_run = gbd_mort_run))

    if (is.numeric(dt_subnat[[1]])) {
      print(paste0("skipping ", subnat, " for draw", i, " because subnational data object failed to be made"))
      next
    }

    if (popadjust == TRUE) {attributes(dt_subnat[[1]])$eppfp$popadjust <- as.integer(1)} else {attributes(dt_subnat[[1]])$eppfp$popadjust <- as.integer(0)}
    attributes(dt_subnat[[1]])$eppfp$likelihood_type <- likelihood_type
    if (eppmod_setting == "rhybrid") {
      fp <- attributes(dt_subnat[[1]])$eppfp
      fp <- prepare_rhybrid(fp, tsEpidemicStart = fp$tsEpidemicStart, rw_start = 2003, rw_trans = 5, rw_dk = 5)
      attributes(dt_subnat[[1]])$eppfp <- fp
      attributes(dt_subnat[[1]])$eppfp$eppmod <- eppmod_setting
    }


    dt_subnat <- extend.trans.params(dt_subnat, start.year, stop.year)
    print("extended params")
    dt_subnat <- sub.off.art.adm2(dt_subnat, loc, mort_draw, subnat, shapefile_version, mean_test = mean_test, tirr = tirr)
    print("subbed off art mort")
    dt_subnat <- sub.on.art.adm2(dt_subnat, loc, mort_draw, subnat, shapefile_version, mean_test = mean_test, tirr = tirr)
    print("subbed on art mort")
    dt_subnat <- sub.cd4.prog.adm2(dt_subnat, loc, mort_draw, subnat, shapefile_version, mean_test = mean_test, tirr = tirr)
    print("subbed pregression params")

    if (likelihood_type == "multi_norm") {
      load(paste0("<<<< FILEPATH REDACTED >>>>"))
      attributes(dt_subnat[[1]])$likdat$means <- means
      attributes(dt_subnat[[1]])$likdat$cov_mat <- mat

    }

    saveRDS(dt_subnat, file = paste0(subnat_dir,subnat,"_", i,"_input.rds"))


  }
  epp.string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=1.5G -l h_rt=01:",loc.table$ad2_time[which(loc.table$ihme_loc_id == loc)],":00 -p -1023 ",
                       "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/",loc,"/ ",
                       "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/",loc,"/ ",
                       "-N ",loc,"_adm2_epp ",
                       "-t 1:", length(ad2_list), " ",
                       "-v sing_image=<<<< FILEPATH REDACTED >>>>singularity-images/lbd/releases/lbd_full_20200128.simg ",
                       core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                       code.dir, "gbd/main_subnat_geo.R ",
                       run.name, " ", loc, " ", proj.end," ", no.anc, " ", i)
  print(epp.string)
  system(epp.string)
}

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

  saveRDS(dt, file = paste0("<<<< FILEPATH REDACTED >>>>"))
  if (length(names(dt)) == 1) {
    fit[[subpop]] <- fitmod(dt[[subpop]], equil.rprior = eq.prior, B0 = 2e5, B = 1e3, B.re = 1e3, number_k = 500, D = 4, opt_iter = 3, no.anc = no.anc, likelihood_type = likelihood_type, prev_type = prev_type)
  } else {# B0 should be 2e5, and B.re should be 1e3, k should be 500 currently testing behavior with the multi norm
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
