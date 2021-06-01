################################################################################
## Purpose: This is to launch the clean up phase of the geospatial epp model
################################################################################

## LBD Base Set Up #############################################################
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>/"mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)


### load the config csv from the geo code directory
######################################################################################################################
run.name <- "2020_05_18_art"
######################################################################################################################
config <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))

### GBD addititional setup ##################################################

windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0(indic_repo, "epp/2_modeling/epp_geo/")
date <- substr(gsub("-","",Sys.Date()),3,8)

## Packages
library(data.table)

new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>", "/r_packages/")
if (!dir.exists(new_pkg_lib)) {dir.create(new_pkg_lib)}

.libPaths(new_pkg_lib)
test_pkg_list <- c('slackr', "RMariaDB", "DBI", "RMySQL", "emdbook", "mvtnorm")
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(slackr, lib = new_pkg_lib)



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
compile_art <- TRUE
rake_prev_internal <- TRUE

### Paths
dir <- paste0("<<<< FILEPATH REDACTED >>>>", run.name, "/")


### Functions
library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")


### Tables
loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

### Code
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



# Draw compilation


for (loc in loc.list) {
   draw.string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=2G -l h_rt=24:00:00 -l archive=TRUE ",
                        "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                        "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                        "-N ", loc, "_save_draws_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                        core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                        code.dir, "gbd/save_paired_draws_geo_new.R ",
                        loc, " ", run.name, " ", n.draws)
  print(draw.string)
  system(draw.string)
}

combine.holds <- paste(paste0(loc.list, "_save_draws_",run.name), collapse = ",")

## Raking
if (rake_inc_internal == TRUE) {
  for (loc in loc.list) {

    internal_inc_rake_string <-  paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=5G -l h_rt=03:20:00 -l archive=TRUE ",
                                        "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                                        "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                                        "-N ", loc, "_internal_inc_raking_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                                        "-hold_jid ", combine.holds, " ",
                                        core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                                        code.dir, "gbd/internal_rake_inc_geo.R ",
                                        loc, " ", run.name)
    print(internal_inc_rake_string)
    system(internal_inc_rake_string)
  }
}


if (rake_mort_internal == TRUE) {
  for (loc in loc.list) {

    internal_mort_rake_string <-  paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=5G -l h_rt=03:10:00 -l archive=TRUE ",
                                         "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                                         "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                                         "-N ", loc, "_internal_mort_raking_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                                         "-hold_jid ", combine.holds, " ",
                                         core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                                         code.dir, "gbd/internal_rake_mort_geo.R ",
                                         loc, " ", run.name)
    print(internal_mort_rake_string)
    system(internal_mort_rake_string)
  }
}


if (rake_prev_internal == TRUE) {
  for (loc in loc.list) {

    internal_prev_rake_string <-  paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=5G -l h_rt=03:10:00 -l archive=TRUE ",
                                         "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                                         "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                                         "-N ", loc, "_internal_prev_raking_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                                         "-hold_jid ", combine.holds, " ",
                                         core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                                         code.dir, "gbd/internal_rake_prev_geo.R ",
                                         loc, " ", run.name)
    print(internal_prev_rake_string)
    system(internal_prev_rake_string)
  }
}


if (compile_art == TRUE) {
  for (loc in loc.list) {

    compile_art_string <-  paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=2G -l h_rt=01:10:00 -l archive=TRUE ",
                                         "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                                         "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                                         "-N ", loc, "_compile_art_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                                         "-hold_jid ", combine.holds, " ",
                                         core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                                         code.dir, "gbd/compile_art_geo.R ",
                                         loc, " ", run.name)
    print(compile_art_string)
    system(compile_art_string)
  }
}






in.inc.holds <- paste(paste0(loc.list, "_internal_inc_raking_",run.name), collapse = ",")
in.mort.holds <- paste(paste0(loc.list, "_internal_mort_raking_",run.name), collapse = ",")
in.prev.holds <- paste(paste0(loc.list, "_internal_prev_raking_",run.name), collapse = ",")

if (rake_inc_gbd == TRUE) {
  for (loc in loc.list) {

    gbd_inc_rake_string <-  paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=4G -l h_rt=01:40:00 -l archive=TRUE ",
                                   "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                                   "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                                   "-N ", loc, "_gbd_inc_raking_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                                   "-hold_jid ", in.inc.holds, " ",
                                   core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                                   code.dir, "gbd/gbd_rake_inc_geo.R ",
                                   loc, " ", run.name, " ", gbd_rake_target)
    print(gbd_inc_rake_string)
    system(gbd_inc_rake_string)
  }
}

if (rake_mort_gbd == TRUE) {
  for (loc in loc.list) {

    gbd_mort_rake_string <-  paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=4G -l h_rt=01:20:00 -l archive=TRUE ",
                                    "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                                    "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                                    "-N ", loc, "_gbd_mort_raking_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                                    "-hold_jid ", in.mort.holds, " ",
                                    core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                                    code.dir, "gbd/gbd_rake_mort_geo.R ",
                                    loc, " ", run.name, " GBD2019")
    print(gbd_mort_rake_string)
    system(gbd_mort_rake_string)
  }
}

rake_prev_gbd <- TRUE
if (rake_prev_gbd == TRUE) {
  for (loc in loc.list) {

    gbd_prev_rake_string <-  paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=1 -l m_mem_free=4G -l h_rt=01:20:00 -l archive=TRUE ",
                                    "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                                    "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                                    "-N ", loc, "_gbd_prev_raking_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                                    "-hold_jid ", in.prev.holds, " ",
                                    core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                                    code.dir, "gbd/gbd_rake_prev_geo.R ",
                                    loc, " ", run.name, " GBD2019")
    print(gbd_prev_rake_string)
    system(gbd_prev_rake_string)
  }
}






gbd.inc.holds <- paste(paste0(loc.list, "_gbd_inc_raking_",run.name), collapse = ",")
gbd.mort.holds <- paste(paste0(loc.list, "_gbd_mort_raking_",run.name), collapse = ",")
plot.holds <- c(in.inc.holds, in.mort.holds)

## Plot
for (loc in loc.list) {
  # Plot individual locations
  plot.string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=2 -l m_mem_free=15G -l h_rt=18:00:00 -l archive=TRUE ",
                        "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                        "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                        "-N ", loc, "_plots_",run.name," -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                        "-hold_jid ", plot.holds, " ",
                        core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                        code.dir, "gbd/plot_epp_geo.R ",
                        loc, " ", run.name, " ", comp_run)
  print(plot.string)
  system(plot.string)
}


rfs.holds <- paste(paste0(loc.list, "_plots_",run.name), collapse = ",")
# Plot rfs

plot_inc_rfs_string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=2 -l m_mem_free=15G -l h_rt=10:00:00 -l archive=TRUE ",
                              "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                              "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                              "-N inc_rfs_plots -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                              "-hold_jid ", plot.holds, " ",
                              core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                              code.dir, "gbd/plot_epp_inc_rfs_geo.r ",
                              " ", run.name, " GBD2019")
print(plot_inc_rfs_string)
system(plot_inc_rfs_string)

plot_mort_rfs_string <- paste0("qsub -P ", cluster.project, " -q ",cluster.queue," -l fthread=2 -l m_mem_free=15G -l h_rt=10:00:00 -l archive=TRUE ",
                               "-e <<<< FILEPATH REDACTED >>>>", run.name,"/errors/ ",
                               "-o <<<< FILEPATH REDACTED >>>>", run.name,"/output/ ",
                               "-N mort_rfs_plots -v sing_image=<<<< FILEPATH REDACTED >>>> ",
                               "-hold_jid ", plot.holds, " ",
                               core_repo, "mbg_central/share_scripts/shell_sing.sh ",
                               code.dir, "gbd/plot_epp_mort_rfs_geo.r ",
                               " ", run.name, " GBD2019")
print(plot_mort_rfs_string)
system(plot_mort_rfs_string)

### End
