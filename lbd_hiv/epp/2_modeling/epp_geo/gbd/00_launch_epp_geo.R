################################################################################
## Purpose: This is to launch the geospatial epp model

################################################################################

## LBD Base Set Up #############################################################
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))

## Set repo
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>/mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_hiv/")

setwd(core_repo)

### GBD addititional setup ##################################################

windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0(indic_repo, "epp/2_modeling/epp_geo/")
date <- substr(gsub("-","",Sys.Date()),3,8)

## Packages
library(data.table)


new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>", "/r_packages/")
if (!dir.exists(new_pkg_lib)) {dir.create(new_pkg_lib)}

.libPaths(new_pkg_lib)
test_pkg_list <- c('slackr', "RMariaDB", "DBI", "RMySQL", "vctrs")
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}




message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

library(slackr, lib = new_pkg_lib)
library(sf)
source("<<<< FILEPATH REDACTED >>>>/get_outputs.R")

### load the config csv from the geo code directory

config <- read.csv(paste0(indic_repo,"epp/2_modeling/epp_config.csv"))



######################################################################################################################
run.name <- as.character(config$Value[which(config$Setting == "run.name")])
######################################################################################################################
proj.end <- as.numeric(as.character(config$Value[which(config$Setting == "proj.end")]))
n.draws <- as.numeric(as.character(config$Value[which(config$Setting == "n.draws")]))
cluster.project <- as.character(config$Value[which(config$Setting == "cluster.project")])
cluster.queue <- as.character(config$Value[which(config$Setting == "cluster.queue")])
no.anc <- as.logical(config$Value[which(config$Setting == "no.anc")])

backcast <- as.logical(config$Value[which(config$Setting == "backcast")])
plot <- as.logical(config$Value[which(config$Setting == "plot")])
run.admin2 <- as.logical(config$Value[which(config$Setting == "run.admin2")])
comp_run <- as.character(config$Value[which(config$Setting == "comp_run")])
rake_inc_internal <- as.logical(config$Value[which(config$Setting == "rake_inc_internal")])
rake_inc_gbd <- as.logical(config$Value[which(config$Setting == "rake_inc_gbd")])
use_gbd_subnats <- as.character(config$Value[which(config$Setting == "use_gbd_subnats")])
use_gbd_subnats <- strsplit(use_gbd_subnats, "+", fixed = T)
loc.list <- as.character(config$Value[which(config$Setting == "loc.list")])
loc.list <- strsplit(loc.list, "+", fixed = T)
gbd_rake_target <- as.character(config$Value[which(config$Setting == "gbd_rake_target")])
rake_mort_internal <- as.logical(config$Value[which(config$Setting == "rake_mort_internal")])
rake_mort_gbd <- as.logical(config$Value[which(config$Setting == "rake_mort_gbd")])
art_type <- as.character(config$Value[which(config$Setting == "art_type")])
gbd_round <- as.numeric(as.character(config$Value[which(config$Setting == "gbd_round")]))
decomp_step <- as.character(config$Value[which(config$Setting == "decomp_step")])
lt <- as.character(config$Value[which(config$Setting == "likelihood")])


message(paste0("run name is ",run.name))
message(paste0("proj end is ",proj.end))
message(paste0("n.draws is ",n.draws))
message(paste0("cluster_project is ",cluster.project))
message(paste0("cluster_queue is ",cluster.queue))
message(paste0("no.anc is ",no.anc))
message(paste0("backcast is ",backcast))
message(paste0("plot is ",plot))
message(paste0("run.admin2 is ", run.admin2))
message(paste0("comp_run is ",comp_run))
message(paste0("art_type is ",art_type))


### Paths
dir <- paste0("<<<< FILEPATH REDACTED >>>>", run.name, "/")
dir.create(dir)

### Save the config and record the relative git statuses
write.csv(config, file = paste0(dir, "config_used.csv"))

record_git_status(core_repo, indic_repo, show_diff = F, file = paste0(dir, "git_at_launch.txt"), check_core_repo = FALSE)


### generate random draws from the mortality and transition parameters

draws <- sample(1:1000, n.draws)
saveRDS(draws, file = paste0(dir, "mort_draws.RDS"))

### Functions
library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")
source(paste0(code.dir, "cache_population.R"))

### Tables
loc.table <- fread("<<<< FILEPATH REDACTED >>>>")
KEN_list <- grep("KEN", loc.table$ihme_loc_id, value = TRUE)
loc.table$epp[which(loc.table$ihme_loc_id %in% KEN_list & loc.table$unaids_recent == 2018)] <- 1 #switch to 2018 KEN file
loc.table$epp[which(loc.table$ihme_loc_id %in% KEN_list & loc.table$unaids_recent != 2018)] <- 0 #switch to 2018 KEN file
loc.table <- loc.table[which(loc.table$ihme_loc_id != "KEN_35646"), ] #remove duplicate Nairobi data
loc.table$epp[which(loc.table$ihme_loc_id == "MRT")] <-

write.csv(loc.table, file = paste0(dir, "loc_table_used.csv"))

### Code
epp.list <- sort(loc.table[epp == 1, ihme_loc_id])

if (length(loc.list[[1]]) == 0) {
loc.list <- epp.list  # if this is set to a list that is shorter than the total possibility then it is for testing
loc.list <- c(loc.list, "ZAF", "NGA", "KEN", "ETH")
loc.list <- as.list(loc.list)
} else {loc.list <- loc.list[[1]]}

if (length(use_gbd_subnats[[1]]) == 0) {
  loc.list <- loc.list[which(nchar(loc.list) <= 3)] #only want to do this for countries not all of the GBD subnationals
} else {
    gbd_subnats <- c()
    for (p in use_gbd_subnats[[1]]) {
      gbd_subnats[[p]] <- grep(paste0(as.character(p), "_"), loc.list, value = TRUE)
    }
    loc.list <- loc.list[which(nchar(loc.list) <= 3)]
    for (p in use_gbd_subnats[[1]]) {
      loc.list <- c(loc.list, gbd_subnats[[p]])
    }
    gbd_subnats <- NULL
}


# Cache populations
if (!file.exists(paste0(dir, "populations/"))) {
  print("Caching populations")
  ind.ur.list <- loc.table[grepl("IND", ihme_loc_id) & level == 5, location_id]
  cache_population(location_id = c(loc.table[ihme_loc_id %in% epp.list, location_id], ind.ur.list, 179,180,196,214), age_group_id = 8:15, dir = paste0(dir, "populations/"), gbd_round = gbd_round, decomp_step = decomp_step)
}


#Check Inputs
source(paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_hiv/epp/2_modeling/epp_geo/gbd/check_epp_inputs_geo.r"))
check_epp_inputs(config)


## Launch EPP
for (loc in loc.list) {

  errors.dir <- paste0("<<<< FILEPATH REDACTED >>>>", run.name,"/errors/")
  output.dir <- paste0("<<<< FILEPATH REDACTED >>>>", run.name,"/output/")
  dir.create(errors.dir)
  dir.create(output.dir)

  errors.dir <- paste0("<<<< FILEPATH REDACTED >>>>", run.name,"/errors/", loc, "/")
  output.dir <- paste0("<<<< FILEPATH REDACTED >>>>", run.name,"/output/", loc, "/")
  dir.create(errors.dir)
  dir.create(output.dir)

  epp.string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=3G -l h_rt=7:00:00 -l archive=TRUE -p -1023 ",
                       "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/",loc,"/ ",
                       "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/",loc,"/ ",
                       "-N ", loc, "_epp_test ",
                       "-t 1:", n.draws, " ",
                       "-v sing_image=<<<< FILEPATH REDACTED >>>>",
                       core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                       code.dir, "gbd/main_geo.R ",
                       run.name, " ", loc, " ", proj.end," ", no.anc)
  print(epp.string)
  system(epp.string)
}




### End
