################################################################################
## Purpose: This script compiles the modeled art coverage.
################################################################################
## LBD Base Set Up
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir <- sprintf("<<<< FILEPATH REDACTED>>>>/mbg/common_inputs")
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED>>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED>>>>", "/lbd_hiv/")

setwd(core_repo)ellfldkdeeznz

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

### GBD addititions Setup


windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED>>>>", "/lbd_hiv/epp/2_modeling/hiv_gbd2019/decomp2019/epp-feature-reset/")

## Packages
library(data.table); library(mvtnorm); library(survey)
new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED>>>>", "/r_packages/")
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
config <- read.csv(paste0("<<<< FILEPATH REDACTED>>>>"))


### Paths
input.dir <- paste0("<<<< FILEPATH REDACTED>>>>", run.name, "/", loc)
adm2.dir <- paste0(input.dir,"/admin2_models")
subnats <- list.dirs(adm2.dir, full.names = FALSE)
subnats <- as.list(subnats)
subnats <- subnats[which(nchar(subnats) > 0)]

nat_art <- fread(paste0(input.dir,"/",loc,"_SPU_art_draws.csv"))
ndraws <- NCOL(nat_art) - 1
overs <- paste0("V", 1:ndraws)
names(nat_art) <- c("year", overs)
nat_art[,mean := rowMeans(.SD) , .SD = overs]


for (s in subnats) {
  s_art <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_art_draws.csv"))
  s_pop <- fread(paste0(adm2.dir,"/", s,"/", loc,"_",s,"_SPU_pop_draws.csv"))
  s_art$subnat <- s
  s_pop$subnat <- s
  if (s == subnats[[1]]) {
    adm2_art <- s_art
    adm2_pop <- s_pop} else {
      adm2_art <- rbind(adm2_art, s_art, fill = TRUE)
      adm2_pop <- rbind(adm2_pop, s_pop, fill = TRUE)
    }
}

adm2_pop$pop <- adm2_pop$`1`
adm2_pop <- adm2_pop[ ,c("year", "pop", "subnat")]
adm2_art <- merge(adm2_art, adm2_pop, by = c("year", "subnat"))
names(adm2_art) <- c("year", "subnat", overs, "pop")

adm2_art[,mean := rowMeans(.SD) , .SD = overs]


ur <- adm2_art

write.csv(ur, paste0(input.dir, "/", loc,"_unraked_adm2_art.csv"), row.names = FALSE)


