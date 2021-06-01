################################################################################
## Purpose: This script rakes the HIV+ deaths comming out of the admin 2 EPP
##          models to the HIV+ mortality coming out ouf the higher level EPP
##          model.  It also compiles all of the mortaltiy numbers from the admin
##          2 models ad saves them in one csv.
## Run instructions: This gets launched from the main EPP launch script once for
##                   each high level geography.  Make sure that the correct GBD
##                   raking targets are set.  The final production run should
##                   get set to 'GBD2019' but there are a bunch of other GBD
##                   raking options so that you can track the incremental
##                   changes in an estimate as it progresses through the GBD.
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

nat_mort <- fread(paste0(input.dir,"/",loc,"_SPU_mort_draws.csv"))
ndraws <- NCOL(nat_mort) - 1
overs <- paste0("V", 1:ndraws)
names(nat_mort) <- c("year", overs)
nat_mort[,mean := rowMeans(.SD) , .SD = overs]


for (s in subnats) {
  s_mort <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_mort_draws.csv"))
  s_pop <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_pop_draws.csv"))
  s_mort$subnat <- s
  s_pop$subnat <- s
  #message(s)
  if (s == subnats[[1]]) {
    adm2_mort <- s_mort
    adm2_pop <- s_pop} else {
      adm2_mort <- rbind(adm2_mort, s_mort, fill = TRUE)
      adm2_pop <- rbind(adm2_pop, s_pop, fill = TRUE)
    }
}


prev_run <- config$Value[which(config$Setting == "prev_date")]
ex <- readRDS(paste0("<<<< FILEPATH REDACTED >>>>"))
ex <- as.character(unlist(ex))


adm2_mort[which(adm2_mort$subnat %in% ex), c(2:1001)] <- NA # remove the estimates for the supressed admin 2 units



adm2_pop$pop <- adm2_pop$`1`
adm2_pop <- adm2_pop[ ,c("year", "pop", "subnat")]
adm2_mort <- merge(adm2_mort, adm2_pop, by = c("year", "subnat"))
names(adm2_mort) <- c("year", "subnat", overs, "pop")


comp_mort <- adm2_mort[, lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("year")]
comp_mort[,mean := rowMeans(.SD) , .SD = overs]

rfs_mort_internal <- merge(nat_mort[ ,c("year", "mean")], comp_mort[ ,c("year", "mean")], by = c("year"))
names(rfs_mort_internal) <- c("year","nat_model", "subnat_model")
rfs_mort_internal$rf <- rfs_mort_internal$nat_model / rfs_mort_internal$subnat_model
rfs_mort_internal$rf[which(is.na(rfs_mort_internal$rf))] <- 1

adm2_mort <- merge(adm2_mort, rfs_mort_internal, by = c("year"))


ur <- adm2_mort
ur$nat_model <- NULL
ur$subnat_model <- NULL
ur$rf <- NULL
write.csv(ur, paste0(input.dir, "/", loc,"_unraked_adm2_mort.csv"), row.names = FALSE)

adm2_mort[,(overs) := lapply(overs, function(x) get(x) * rf)]
adm2_mort$nat_model <- NULL
adm2_mort$subnat_model <- NULL
adm2_mort$rf <- NULL


write.csv(rfs_mort_internal, paste0(input.dir, "/internal_mort_rfs.csv"), row.names = FALSE)
write.csv(adm2_mort, paste0(input.dir, "/", loc,"_internal_raked_adm2_mort.csv"), row.names = FALSE)

