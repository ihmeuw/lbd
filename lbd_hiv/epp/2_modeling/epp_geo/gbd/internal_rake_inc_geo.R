################################################################################
## Purpose: This script rakes the incidence curves from the constituent admin 2
##          models to be equal to the incidence curve of the high level geography.
##          This process is useful both as a diagnostic and ad prep to create a
##          standard set as a jumping off point for GBD raking.
################################################################################
## LBD Base Set Up
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

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

### GBD addititions Setup


windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_hiv/epp/2_modeling/hiv_gbd2019/decomp2019/epp-feature-reset/")

## Packages
library(data.table); library(mvtnorm); library(survey)
new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>/", "/r_packages/")
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


nat_inc <- fread(paste0(input.dir,"/",loc,"_SPU_inc_draws.csv"))
ndraws <- NCOL(nat_inc) - 1
overs <- paste0("V", 1:ndraws)
names(nat_inc) <- c("year", overs)
nat_inc[,mean := rowMeans(.SD) , .SD = overs]


for (s in subnats) {
  s_inc <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_inc_draws.csv"))
  s_pop <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_pop_draws.csv"))
  s_inc$subnat <- s
  s_pop$subnat <- s
  if (s == subnats[[1]]) {
    adm2_inc <- s_inc
    adm2_pop <- s_pop} else {
      adm2_inc <- rbind(adm2_inc, s_inc, fill = TRUE)
      adm2_pop <- rbind(adm2_pop, s_pop, fill = TRUE)
    }
}

prev_run <- config$Value[which(config$Setting == "prev_date")]
ex <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
ex <- as.character(unlist(ex))

adm2_inc[which(adm2_inc$subnat %in% ex), c(2:1001)] <- NA # remove the estimates for the supressed admin 2 units



adm2_pop$pop <- adm2_pop$`1`
adm2_pop <- adm2_pop[ ,c("year", "pop", "subnat")]
adm2_inc <- merge(adm2_inc, adm2_pop, by = c("year", "subnat"))
names(adm2_inc) <- c("year", "subnat", overs, "pop")


adm2_inc[,(overs) := lapply(overs, function(x) get(x) * pop)]

comp_inc <- adm2_inc[, lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("year")]
names(comp_inc) <- c("year", overs, "pop")
comp_inc[,(overs) := lapply(overs, function(x) get(x) / pop)]
comp_inc[,mean := rowMeans(.SD) , .SD = overs]

rfs_inc_internal <- merge(nat_inc[ ,c("year", "mean")], comp_inc[ ,c("year", "mean")], by = c("year"))
names(rfs_inc_internal) <- c("year","nat_model", "subnat_model")
rfs_inc_internal$rf <- rfs_inc_internal$nat_model / rfs_inc_internal$subnat_model
rfs_inc_internal$rf[which(is.na(rfs_inc_internal$rf))] <- 1

adm2_inc <- merge(adm2_inc, rfs_inc_internal, by = c("year"))
adm2_inc[,(overs) := lapply(overs, function(x) get(x) / pop)]

ur <- adm2_inc
ur$nat_model <- NULL
ur$subnat_model <- NULL
ur$rf <- NULL
write.csv(ur, paste0(input.dir, "/", loc,"_unraked_adm2_inc.csv"), row.names = FALSE)

adm2_inc[,(overs) := lapply(overs, function(x) get(x) * rf)]
adm2_inc$nat_model <- NULL
adm2_inc$subnat_model <- NULL
adm2_inc$rf <- NULL


write.csv(rfs_inc_internal, paste0(input.dir, "/internal_inc_rfs.csv"), row.names = FALSE)
write.csv(adm2_inc, paste0(input.dir, "/", loc,"_internal_raked_adm2_inc.csv"), row.names = FALSE)

