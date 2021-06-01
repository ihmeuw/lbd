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

windows <- Sys.info()[1][["sysname"]] == "Windows"

## Packages
library(data.table); library(mvtnorm); library(survey)
new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>", "/r_packages/")
.libPaths(new_pkg_lib)
test_pkg_list <- c('slackr', "assertable", "RMariaDB", "rlang")
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(slackr)
library(assertable)


message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

### GBD addititions Setup


#library(RMariaDB, lib = new_pkg_lib)
library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>/02_mortality/shared/r")

## Arguments

run.name <- as.character(commandArgs()[5])
message(paste0("the run name is ",run.name))
loc <- as.character(commandArgs()[4])
message(paste0("the loc is ",loc))
gbd_rake_target <- as.character(commandArgs()[6])
if (length(grep("+", gbd_rake_target, fixed = TRUE)) > 0 ) {
  gbd_rake_target <- strsplit(gbd_rake_target, "+", fixed = T)
  gbd_rake_target <- gbd_rake_target[[1]]
}
gbd_rake_target <- "GBD2019"
### read in the config file
config <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))
shp_version <- config$Value[which(config$Setting == "shapefile_version")]

### Paths
input.dir <- paste0('<<<< FILEPATH REDACTED >>>>', run.name, "/", loc)


loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

for (t in gbd_rake_target) {
  t <- "GBD2019"
  if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
    nat_prev <- fread(paste0(input.dir,"/",loc,"_SPU_prev_draws.csv"))
    ndraws <- NCOL(nat_prev) - 1
    overs <- paste0("V", 1:ndraws)
    names(nat_prev) <- c("year", overs)
    nat_prev[,mean := rowMeans(.SD) , .SD = overs]
  } else {
    nat_prev <- fread(paste0(input.dir,"/",loc,"_unraked_adm2_prev.csv"))
    sp_h <- read.dbf(paste0("<<<< FILEPATH REDACTED >>>>"))
    nat_prev <- as.data.table(merge(sp_h, nat_prev, by.x = c("ADM2_CODE"), by.y = c("subnat")))
    overs <- names(nat_prev)
    overs <- overs[which(overs %!in% c("ADM2_CODE","NAME_0", "NAME_1", "NAME_2", "geo_id", "ad2_id", "ad0_parent", "ad1_parent", "ADM2_NAME","ADM1_CODE", "ADM1_NAME", "ADM0_CODE", "ADM0_NAME", "year", "pop"))]
    nat_prev[, (overs) := lapply(overs, function(x) get(x) * pop )]
    nat_prev <- nat_prev[,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM1_CODE","ADM1_NAME", "year")]
    names(nat_prev) <- c("ADM1_CODE","ADM1_NAME", "year", overs, "pop")
    nat_prev[, (overs) := lapply(overs, function(x) get(x) / pop )]
    nat_prev[,mean := rowMeans(.SD) , .SD = overs]
    nat_prev <- nat_prev[,c("ADM1_CODE","ADM1_NAME", "year", "mean", "pop")]

  }


  if (t == "GBD2017") {
    source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
    comp_prev <- get_outputs(topic = "cause",
                            version = "best", # using the latest version instead of the best version on Sept 18, 2018. Revisit as other versions become available
                            gbd_round_id = 5,
                            cause_id = 298,
                            measure_id = 5,
                            metric_id = 3,
                            age_group_id = 24,
                            location_id = loc.table$location_id[which(loc.table$ihme_loc_id == loc)],
                            year_id = c(1990:2017))
    comp_prev <- comp_prev[,c("year_id", "val")]
    names(comp_prev) <- c("year", "mean")
    rfs_prev_gbd <- merge(comp_prev, nat_prev[ ,c("year", "mean")], by = c("year"))
    names(rfs_prev_gbd) <- c("year","gbd_model", "lbd_model")
    rfs_prev_gbd$gbd_rf <- (rfs_prev_gbd$gbd_model * 100) / rfs_prev_gbd$lbd_model # the scalar of 100 is applied to the output of get output because EPP reports out previdence in cases per 100 (or at least I hope so)
    rfs_prev_gbd$gbd_rf[which(is.na(rfs_prev_gbd$gbd_rf))] <- 1
  }






  ##################################################################################################################################


  if (t == "GBD2019") {
    source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')

    if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
      comp_prev <- get_outputs(topic = "cause",
                              version = 'latest',
                              gbd_round_id = 6,
                              cause_id = 298,
                              measure_id = 5,
                              metric_id = 3,
                              age_group_id = 24,
                              location_id = loc.table$location_id[which(loc.table$ihme_loc_id == loc)],
                              year_id = c(1990:2019),
                              decomp_step = "step5")
      comp_prev <- comp_prev[,c("year_id", "val")]
      names(comp_prev) <- c("year", "mean")
      rfs_prev_gbd <- merge(comp_prev, nat_prev[ ,c("year", "mean")], by = c("year"))
      names(rfs_prev_gbd) <- c("year","gbd_model", "lbd_model")
      rfs_prev_gbd$gbd_rf <- (rfs_prev_gbd$gbd_model * 100) / rfs_prev_gbd$lbd_model # the scalar of 100 is applied to the output of get output because EPP reports out previdence in cases per 100 (or at least I hope so)
      rfs_prev_gbd$gbd_rf[which(is.na(rfs_prev_gbd$gbd_rf))] <- 1
    } else {

      connector <- get_gbd_locs(rake_subnational = T,
                                reg = loc,
                                shapefile_version = shp_version)
      nat_prev <- merge(nat_prev, connector, by = c("ADM1_CODE"))
      comp_prev <- get_outputs(topic = "cause",
                              version = 'latest',
                              gbd_round_id = 6,
                              cause_id = 298,
                              measure_id = 5,
                              metric_id = 3,
                              age_group_id = 24,
                              location_id = unique(nat_prev$location_id),
                              year_id = c(1990:2019),
                              decomp_step = "step5")
      comp_prev <- comp_prev[,c("location_id","year_id", "val")]
      names(comp_prev) <- c("location_id","year", "mean")
      rfs_prev_gbd <- merge(comp_prev, nat_prev[ ,c("location_id", "ADM1_CODE","year", "mean")], by = c("location_id", "year"))
      names(rfs_prev_gbd) <- c("location_id","year","gbd_model","ADM1_CODE","lbd_model")
      rfs_prev_gbd$gbd_rf <- (rfs_prev_gbd$gbd_model * 100) / rfs_prev_gbd$lbd_model # the scalar of 100 is applied to the output of get output because EPP reports out previdence in cases per 100 (or at least I hope so)
      rfs_prev_gbd$gbd_rf[which(is.na(rfs_prev_gbd$gbd_rf))] <- 1
    }



  }


  ####################################################################################################################################

  if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
    nat_prev <- merge(nat_prev, rfs_prev_gbd, by = c("year"))
    nat_prev[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
    nat_prev$unraked_mean <- nat_prev$mean
    nat_prev$mean <- NULL
    nat_prev[,raked_mean := rowMeans(.SD) , .SD = overs]
    nat_prev$lbd_model <- NULL


    adm2_prev <- fread(paste0(input.dir, "/", loc,"_internal_raked_adm2_prev.csv"))
    adm2_prev <- merge(adm2_prev, rfs_prev_gbd, by = c("year"))
    adm2_prev[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
    adm2_prev$lbd_model <- NULL
    adm2_prev$gbd_model <- NULL
    adm2_prev[,double_raked_mean := rowMeans(.SD) , .SD = overs]


    valid_prev <- adm2_prev
    valid_prev[,(overs) := lapply(overs, function(x) get(x) * pop)]

    comp_prev <- valid_prev[, lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("year")]
    names(comp_prev) <- c("year", overs, "pop")
    comp_prev[,(overs) := lapply(overs, function(x) get(x) / pop)]
    comp_prev[,mean := rowMeans(.SD) , .SD = overs]
    valid_prev <- merge(nat_prev[ ,c("year", "raked_mean")], comp_prev[ ,c("year", "mean")], by = c("year"))
    valid_prev$test <- valid_prev$raked_mean / valid_prev$mean
    valid_prev$test[which(is.na(valid_prev$test))] <- 1
    if (max(valid_prev$test) == Inf) {
      message( "there is a /0 causing issues, troubleshooting needed.")}
    valid_prev$test[which(valid_prev$test == Inf)] <- 1

    adm2_prev[,(overs) := lapply(overs, function(x) get(x) / pop)]

    if (as.integer(sum(valid_prev$test)) != as.integer(NROW(valid_prev))) {
      message(paste0("error in ",t," raking")) } else {
        write.csv(rfs_prev_gbd, paste0(input.dir, "/",t,"_prev_rfs.csv"), row.names = FALSE)
        write.csv(adm2_prev, paste0(input.dir, "/", loc,"_",t,"_raked_adm2_prev.csv"), row.names = FALSE)
        write.csv(nat_prev, paste0(input.dir, "/", loc,"_",t,"_raked_nat_prev.csv"), row.names = FALSE)
      }

  } else {

    adm2_prev <- fread(paste0(input.dir,"/",loc,"_unraked_adm2_prev.csv"))
    adm2_prev <- merge(adm2_prev, sp_h, by.x = c("subnat"), by.y = c("ADM2_CODE"))
    adm2_prev <- merge(adm2_prev, rfs_prev_gbd, by = c("year", "ADM1_CODE"))
    adm2_prev[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
    adm2_prev$lbd_model <- NULL
    adm2_prev$gbd_model <- NULL
    adm2_prev[,double_raked_mean := rowMeans(.SD) , .SD = overs]
    adm2_prev <- as.data.frame(adm2_prev)
    adm2_prev <- adm2_prev[,c("year","subnat",overs,"pop","gbd_rf","double_raked_mean")]

    write.csv(rfs_prev_gbd, paste0(input.dir, "/",t,"_prev_rfs.csv"), row.names = FALSE)
    write.csv(adm2_prev, paste0(input.dir, "/", loc,"_",t,"_raked_adm2_prev.csv"), row.names = FALSE)
  }

}
