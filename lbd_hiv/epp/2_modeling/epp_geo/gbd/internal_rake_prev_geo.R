################################################################################
## Purpose: This script rakes the prevalence curves from the constituent admin 2
##          models to be equal to the prevalence curve of the high level geography.
##          This process is useful both as a diagnostic and ad prep to create a
##          standard set as a jumping off point for GBD raking.
## Run instructions: Launched form the general launch script once for each country.
##                   This script compiles the admin 2 data into one location and
##                   performs the initial internal raking.
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
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>")

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

run.name <- as.character(commandArgs()[5])
message(paste0("the run name is ",run.name))
loc <- as.character(commandArgs()[4])
message(paste0("the loc is ",loc))


### read in the config file
config <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))


### Paths
input.dir <- paste0('<<<< FILEPATH REDACTED >>>>', run.name, "/", loc)
adm2.dir <- paste0(input.dir,"/admin2_models")
subnats <- list.dirs(adm2.dir, full.names = FALSE)
subnats <- as.list(subnats)
subnats <- subnats[which(nchar(subnats) > 0)]

nat_prev <- fread(paste0(input.dir,"/",loc,"_SPU_prev_draws.csv"))
ndraws <- NCOL(nat_prev) - 1
overs <- paste0("V", 1:ndraws)
names(nat_prev) <- c("year", overs)
nat_prev[,mean := rowMeans(.SD) , .SD = overs]


for (s in subnats) {
  s_prev <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_prev_draws.csv"))
  s_pop <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_pop_draws.csv"))
  s_prev$subnat <- s
  s_pop$subnat <- s
  if (s == subnats[[1]]) {
    adm2_prev <- s_prev
    adm2_pop <- s_pop} else {
      adm2_prev <- rbind(adm2_prev, s_prev, fill = TRUE)
      adm2_pop <- rbind(adm2_pop, s_pop, fill = TRUE)
    }
}


prev_run <- config$Value[which(config$Setting == "prev_date")]
ex <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
ex <- as.character(unlist(ex))


adm2_prev[which(adm2_prev$subnat %in% ex), c(2:1001)] <- NA # remove the estimates for the supressed admin 2 units


adm2_pop$pop <- adm2_pop$`1`
adm2_pop <- adm2_pop[ ,c("year", "pop", "subnat")]
adm2_prev <- merge(adm2_prev, adm2_pop, by = c("year", "subnat"))
names(adm2_prev) <- c("year", "subnat", overs, "pop")
adm2_prev[,(overs) := lapply(overs, function(x) get(x) * pop)]

comp_prev <- adm2_prev[, lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("year")]
names(comp_prev) <- c("year", overs, "pop")
comp_prev[,(overs) := lapply(overs, function(x) get(x) / pop)]
comp_prev[,mean := rowMeans(.SD) , .SD = overs]

rfs_prev_internal <- merge(nat_prev[ ,c("year", "mean")], comp_prev[ ,c("year", "mean")], by = c("year"))
names(rfs_prev_internal) <- c("year","nat_model", "subnat_model")
rfs_prev_internal$rf <- rfs_prev_internal$nat_model / rfs_prev_internal$subnat_model
rfs_prev_internal$rf[which(is.na(rfs_prev_internal$rf))] <- 1

adm2_prev <- merge(adm2_prev, rfs_prev_internal, by = c("year"))
adm2_prev[,(overs) := lapply(overs, function(x) get(x) / pop)]

ur <- adm2_prev
ur$nat_model <- NULL
ur$subnat_model <- NULL
ur$rf <- NULL
write.csv(ur, paste0(input.dir, "/", loc,"_unraked_adm2_prev.csv"), row.names = FALSE)

adm2_prev[,(overs) := lapply(overs, function(x) get(x) * rf)]
adm2_prev$nat_model <- NULL
adm2_prev$subnat_model <- NULL
adm2_prev$rf <- NULL


write.csv(rfs_prev_internal, paste0(input.dir, "/internal_prev_rfs.csv"), row.names = FALSE)
write.csv(adm2_prev, paste0(input.dir, "/", loc,"_internal_raked_adm2_prev.csv"), row.names = FALSE)

