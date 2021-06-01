################################################################################
## Purpose: This script searches for EPP model draws that should exist but don't
##          and it releaunches those draws so that the EPP model can be used
##          to create a full posterior.
## Run instructions: Launched from the launch script, once per high level
##                   geography.  This script is also the main clean up script
##                   used to minimze the size of data on disk.  So, once you run
##                   this you loose a bunch of intermediary files.
################################################################################

# Set up
rm(list = ls())
'%!in%' <- function(x,y)!('%in%'(x,y))
gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
windows <- Sys.info()[1] == "Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>")

## Packages
library(data.table); library(parallel)

## Arguments
core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0(indic_repo, "epp/2_modeling/epp_geo/")
date <- substr(gsub("-","",Sys.Date()),3,8)

run_name <- as.character(commandArgs()[5])
message(paste0("the run name is ",run_name))
cluster.project <- "proj_geo_nodes_hiv"
cluster.queue <- "geospatial.q"

### read in the config file

config <- read.csv("<<<< FILEPATH REDACTED >>>>")

proj.end <- as.numeric(as.character(config$Value[which(config$Setting == "proj.end")]))
no.anc <- as.logical(config$Value[which(config$Setting == "no.anc")])
n.draws <- as.numeric(as.character(config$Value[which(config$Setting == "n.draws")]))
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


n.runs <- as.numeric(as.character(config$Value[which(config$Setting == "n.draws")]))

debug <- as.logical(config$Value[which(config$Setting == "debug")])
message(paste0("the debug setting is ", debug))

loc.list <- as.character(config$Value[which(config$Setting == "loc.list")])
loc.list <- strsplit(loc.list, "+", fixed = T)
use_gbd_subnats <- as.character(config$Value[which(config$Setting == "use_gbd_subnats")])
use_gbd_subnats <- strsplit(use_gbd_subnats, "+", fixed = T)

loc.table <- fread("<<<< FILEPATH REDACTED >>>>")

epp.list <- sort(loc.table[epp == 1, ihme_loc_id])

if (length(loc.list[[1]]) == 0) {
  loc.list <- epp.list  # if this is set to a list that is shorter than the total possibility then it is for testing
  loc.list <- c(loc.list, "ZAF", "NGA", "KEN", "ETH", "IND", "MAR", "VNM")
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

loc.list <- loc.list[which(loc.list %!in% c("IND", "MAR", "VNM", "MMR", "KHM", "PNG", "HTI", "DOM"))] # removing countries that I currently expect to fail


for (loc in loc.list) {
  print(loc)
  ### Paths
  in.dir <- paste0('<<<< FILEPATH REDACTED >>>>')
  inc_dir <- paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/")
  inc_path <- paste0(inc_dir, loc, "_SPU_inc_draws.csv")
  ### Code
  missing.runs <- c()
  j = 1
  for (i in c(1:n.runs)) {
    file <- paste0(loc,'_results_incid',i,'.csv')
    if (file.exists(paste0(in.dir, file)) & (file.info(paste0(in.dir, file))$size > 0)) {
      #print(paste("Found run",i))
    }  else {
      #print(paste("Missing run",i))
      missing.runs <- c(missing.runs, i)
    }
  }
  run.admin2 <- as.logical(config$Value[which(config$Setting == "run.admin2")])
  message(paste0("the run.admin2 setting is ", run.admin2))
  stash_dir <- paste0("<<<< FILEPATH REDACTED >>>>",run_name,"/missing_runs_for_rerun/")
  dir.create(stash_dir)
  saveRDS(missing.runs , file = paste0(stash_dir, loc,"_missing_runs.rds"))
  if (run.admin2 == TRUE) {
    ad2_dir <- paste0(in.dir, "admin2_models/")
    admin_dirs <- list.dirs(ad2_dir, full.names = FALSE)
    for (f in 2:length(admin_dirs)) {
      subnat <- admin_dirs[[f]]
      print(subnat)
      subnat_dir <- paste0(ad2_dir, subnat, "/")
      inc_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_inc_draws.csv")

      ### Code
      ## Incidence


      missing.runs <- c()

      j = 1


      for (i in c(1:n.runs)) {
        file <- paste0(subnat,'_results_incid',i,'.csv')
        if ((file.exists(paste0(subnat_dir, file))) & (file.info(paste0(subnat_dir, file))$size > 0)) {
          #print(paste("Found run",i))
        }  else {
          #print(paste("Missing run",i))
          missing.runs <- c(missing.runs, i)
        }
      }
      saveRDS(missing.runs , file = paste0(stash_dir, loc,"_",subnat,"_missing_runs.rds"))
    }
  }
}


missing <- list.files(stash_dir)
todo <- list()
for (id in c(1:length(missing))) {
  list <- readRDS(paste0(stash_dir,missing[[id]]))
  if (length(list) > 0) {
    todo <- c(todo, missing[[id]])
  }
}


nat_todo <- todo[which(is.na(as.numeric(substr(todo,5,5))))]
subnat_todo <- todo[which(!is.na(as.numeric(substr(todo,5,5))))]



for (f in nat_todo) {
  list <- readRDS(paste0(stash_dir,f))
  loc <- substr(f,0,3)

  for (i in list) {
  epp.string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=7G -l h_rt=12:00:00 -l archive=TRUE -p -1023 ",
                       "-e <<<< FILEPATH REDACTED >>>>", run_name,"/errors/",loc,"/ ",
                       "-o <<<< FILEPATH REDACTED >>>>", run_name,"/output/",loc,"/ ",
                       "-N ", loc, "_epp_test ",
                       "-t ",i,":",i, " ",
                       "-v sing_image=<<<< FILEPATH REDACTED >>>>singularity-images/lbd/releases/lbd_full_20200128.simg ",
                       core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                       code.dir, "gbd/main_geo.R ",
                       run_name, " ", loc, " ", proj.end," ", no.anc)
  print(epp.string)
  system(epp.string)
  }

  nat_covered_subnats <- grep(loc, subnat_todo, value = T)
  if (length(nat_covered_subnats) > 0) {
  for (fi in nat_covered_subnats) {
    list_subnat <- readRDS(paste0(stash_dir,fi))
    subnat <- strsplit(fi,"_")
    subnat <- subnat[[1]][[2]]
    for (id in list_subnat[which(list_subnat %!in% list)]) {
      epp.string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=1.5G -l h_rt=00:60:00 -p -1023 ",
                           "-e <<<< FILEPATH REDACTED >>>>", run_name,"/errors/",loc,"/ ",
                           "-o <<<< FILEPATH REDACTED >>>>", run_name,"/output/",loc,"/ ",
                           "-N ",loc,"_adm2_epp ",
                           "-v sing_image=<<<< FILEPATH REDACTED >>>>singularity-images/lbd/releases/lbd_full_20200128.simg ",
                           core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                           code.dir, "gbd/main_subnat_geo_relaunch.R ",
                           run_name, " ", loc, " ", proj.end," ", no.anc, " ", id, " ", subnat)
      print(epp.string)
      system(epp.string)
    }
  }
  }
}

if (!exists("nat_covered_subnats")) { nat_covered_subnats <- list()}

for (fi in subnat_todo[which(subnat_todo %!in% nat_covered_subnats)]) {
  list_subnat <- readRDS(paste0(stash_dir,fi))
  subnat <- strsplit(fi,"_")
  subnat <- subnat[[1]][[2]]
  loc <- substr(fi,0,3)
  for (id in list_subnat) {
    epp.string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=1.5G -l h_rt=00:60:00 -p -1023 ",
                         "-e <<<< FILEPATH REDACTED >>>>", run_name,"/errors/",loc,"/ ",
                         "-o <<<< FILEPATH REDACTED >>>>", run_name,"/output/",loc,"/ ",
                         "-N ",loc,"_adm2_epp ",
                         "-v sing_image=<<<< FILEPATH REDACTED >>>>singularity-images/lbd/releases/lbd_full_20200128.simg ",
                         core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                         code.dir, "gbd/main_subnat_geo_relaunch.R ",
                         run_name, " ", loc, " ", proj.end," ", no.anc, " ", id, " ", subnat)
    print(epp.string)
    system(epp.string)
  }
}


### End



