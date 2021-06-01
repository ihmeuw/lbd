## LBD Base Set Up
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('"<<<< FILEPATH REDACTED >>>>//mbg/common_inputs')
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


start.year <- 1970
stop.year <- 2019


loc.table <- fread("<<<< FILEPATH REDACTED >>>>")

locs <- list.dirs("<<<< FILEPATH REDACTED >>>>", recursive = F, full.names = F)
locs <- locs[which(locs %in% loc.table$ihme_loc_id)]


gbd_flat <- read.csv("<<<< FILEPATH REDACTED >>>>")

gbd_flat <- gbd_flat[which(gbd_flat$measure == "Scaled_Inc" & gbd_flat$metric == "Rate"), ]
gbd_flat <- gbd_flat[which(gbd_flat$loc %in% locs), ]
gbd_flat <- gbd_flat[which(gbd_flat$age_group_id %in% c(8:14)), ]


gbd_flat$base <- NA
gbd_flat$IRR <- NA

for (l in locs) {
  loc <- l
  ####### using the GBD flat file ###############
  for (y in unique(gbd_flat$year_id)) {
    val <- gbd_flat$mean[which(gbd_flat$sex_id == 1 & gbd_flat$age_group_id == 10 & gbd_flat$year_id == y & gbd_flat$loc == loc)]
    gbd_flat$base[which(gbd_flat$year_id == y & gbd_flat$loc == loc)] <- val

  }

}

gbd_flat$IRR <- gbd_flat$mean / gbd_flat$base

write.csv(gbd_flat[ ,c("loc", "year_id", "sex_id", "age_group_id", "measure", "mean", "IRR")], file = "<<<< FILEPATH REDACTED >>>>")
