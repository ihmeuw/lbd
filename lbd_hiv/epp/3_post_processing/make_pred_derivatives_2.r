################################################################################
## Purpose: Create Summary Estimates at the various geographic levels.

################################################################################



## LBD Base Set Up
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('/mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

new_pkg_lib   <- paste0("/<<<< FILEPATH REDACTED >>>>/", "/r_packages/")
dir.create(new_pkg_lib)
.libPaths(new_pkg_lib)
test_pkg_list <- c('haven', 'readstata13', 'stringdist', 'Metrics', 'tidyverse')
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(dplyr)
library(data.table)
library(rgeos)
library(rgdal)
library(maptools)
library(raster)
library(sf)


### Setup
windows <- Sys.info()[1] == "Windows"
code.dir <- paste0("/<<<< FILEPATH REDACTED >>>>/", "lbd_hiv/epp/2_modeling/epp_geo/gbd/")

## Packages
library(data.table); library(ggplot2)
library(parallel); library(mortdb, lib = "<<<< FILEPATH REDACTED")

## Arguments
runs <- c("2020_05_18_art")


for (r in runs) {
run.name <- r
message(r)
config <- fread(paste0("<<<< FILEPATH REDACTED >>>>"), header = TRUE)

shapefile_version <- config$Value[which(config$Setting == "shapefile_version")]

sp_h <- st_read(paste0("<<<< FILEPATH REDACTED >>>>"))
sp_h <- sp_h[ , c("ADM2_CODE", "ADM2_NAME", "ADM1_CODE", "ADM1_NAME", "ADM0_CODE", "ADM0_NAME")]
sp_h$geometry <- NULL

dir <- paste0("<<<< FILEPATH REDACTED >>>>",run.name,"/admin_summaries")


levels <- c(0:2)
indic <- c("inc", "mort")

for (ind in indic) {
  message(ind)
  for (l in levels) {
    message(l)

    data_c <- fread(paste0(dir,"/admin",l,"_",ind,"_c_GBD2019_raked_draws.csv"))
    data_r <- fread(paste0(dir,"/admin",l,"_",ind,"_GBD2019_raked_draws.csv"))

    data_c <- data_c[which(data_c$year %in% c(2010, 2018)), ]
    data_r <- data_r[which(data_r$year %in% c(2010, 2018)), ]

    ids <- names(data_c)
    ids <- grep("V", ids, value = T, invert = T)
    data_c_long <- melt(data_c, id.vars = ids)
    dcast.vars <- names(data_c_long)
    dcast.vars <- dcast.vars[which(dcast.vars %!in% c("year", "value", "pop", ".id", "gbd_rf", "double_raked_mean"))]
    f <- paste(dcast.vars, collapse = "+")
    f <- paste0(f,"~year")
    f <- as.formula(f)
    data_c_long <- dcast(data_c_long, formula = f , value.var = c("value"))
    data_c_long$diff <- (data_c_long$`2018` - data_c_long$`2010`)
    data_c_long$inc_change_p <- data_c_long$diff / (-1 * data_c_long$`2010`)
    dcast.vars2 <- grep("ADM", names(data_c_long), value = T)
    f2 <- paste(dcast.vars2, collapse = "+")
    f2 <- paste0(f2,"~variable")
    f2 <- as.formula(f2)
    data_c_long <- dcast(data_c_long, formula = f2, value.var = c("inc_change_p"))
    data_c_long$year <- 2018
    data_c_long$pop <- NA
    data_c_long_summary <- make_admin_pred_summary(data_c_long, sp_h, summary_stats = c("mean", "lower", "upper"))

    write.csv(data_c_long_summary, paste0(dir, "/admin",l,"_",ind,"_c_change_2010-2018.csv"))

    ids <- names(data_r)
    ids <- grep("V", ids, value = T, invert = T)
    data_r_long <- melt(data_r, id.vars = ids)
    dcast.vars <- names(data_r_long)
    dcast.vars <- dcast.vars[which(dcast.vars %!in% c("year", "value", "pop", ".id", "gbd_rf", "double_raked_mean"))]
    f <- paste(dcast.vars, collapse = "+")
    f <- paste0(f,"~year")
    f <- as.formula(f)
    data_r_long <- dcast(data_r_long, formula = f , value.var = c("value"))
    data_r_long$diff <- (data_r_long$`2018` - data_r_long$`2010`)
    data_r_long$inc_change_p <- data_r_long$diff / (-1 * data_r_long$`2010`)
    dcast.vars2 <- grep("ADM", names(data_r_long), value = T)
    f2 <- paste(dcast.vars2, collapse = "+")
    f2 <- paste0(f2,"~variable")
    f2 <- as.formula(f2)
    data_r_long <- dcast(data_r_long, formula = f2, value.var = c("inc_change_p"))
    data_r_long$year <- 2018
    data_r_long$pop <- NA
    data_r_long_summary <- make_admin_pred_summary(data_r_long, sp_h, summary_stats = c("mean", "lower", "upper"))

    write.csv(data_r_long_summary, paste0(dir, "/admin",l,"_",ind,"_r_change_2010-2018.csv"))
  }

}



for (l in levels) {
  message(l)
  inc_c <- fread(paste0(dir,"/admin",l,"_inc_c_GBD2019_raked_draws.csv"))
  prev_c <- fread(paste0(dir,"/admin",l,"_prev_c_GBD2019_raked_draws.csv"))
  inc_c$year <- as.numeric(as.character(inc_c$year))
  inc_c$year <- as.numeric(as.character(inc_c$year))

  ids <- names(inc_c)
  ids <- grep("V", ids, value = T, invert = T)
  inc_c_long <- melt(inc_c, id.vars = ids)


  ids <- names(prev_c)
  ids <- grep("V", ids, value = T, invert = T)
  prev_c_long <- melt(prev_c, id.vars = ids)

  ids <- ids[which(ids %!in% c("gbd_rf", "double_raked_mean", ".id"))]

  ids <- c(ids, "variable")
  comb <- merge(inc_c_long, prev_c_long, by = ids)
  comb$ipr <- comb$value.x / comb$value.y

  dcast.vars <- names(comb)
  dcast.vars <- dcast.vars[which(dcast.vars %!in% c("value.x", "value.y", "variable", "ipr", ".id.x", ".id.y", "gbd_rf.x", "gbd_rf.y", "double_raked_mean.x", "double_raked_mean.y"))]
  f <- paste0(dcast.vars, collapse = "+")
  f <- paste0(f, "~variable")
  f <- as.formula(f)

  ipr <- dcast(comb, formula = f, value.var = c("ipr"))

  ipr_summary <- make_admin_pred_summary(ipr, sp_h, summary_stats = c("mean", "lower", "upper"))

  write.csv(ipr_summary, paste0(dir, "/admin",l,"_ipr_summary.csv"))

}

}
