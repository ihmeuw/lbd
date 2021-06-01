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
library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")

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
  if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
    nat_mort <- fread(paste0(input.dir,"/",loc,"_SPU_mort_draws.csv"))
    ndraws <- NCOL(nat_mort) - 1
    overs <- paste0("V", 1:ndraws)
    names(nat_mort) <- c("year", overs)
    nat_mort[,mean := rowMeans(.SD) , .SD = overs]
  } else {
    nat_mort <- fread(paste0(input.dir,"/",loc,"_unraked_adm2_mort.csv"))
    sp_h <- read.dbf(paste0("<<<< FILEPATH REDACTED >>>>"))
    nat_mort <- as.data.table(merge(sp_h, nat_mort, by.x = c("ADM2_CODE"), by.y = c("subnat")))
    overs <- names(nat_mort)
    overs <- overs[which(overs %!in% c("ADM2_CODE","NAME_0", "NAME_1", "NAME_2", "geo_id", "ad2_id", "ad0_parent", "ad1_parent", "ADM2_NAME","ADM1_CODE", "ADM1_NAME", "ADM0_CODE", "ADM0_NAME", "year", "pop"))]
    nat_mort <- nat_mort[,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM1_CODE","ADM1_NAME", "year")]
    names(nat_mort) <- c("ADM1_CODE","ADM1_NAME", "year", overs, "pop")
    nat_mort[,mean := rowMeans(.SD) , .SD = overs]
    nat_mort <- nat_mort[,c("ADM1_CODE","ADM1_NAME", "year", "mean", "pop")]

  }


  if (t == "GBD2017") {
    source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
    comp_mort <- get_outputs(topic = "cause",
                            version = "best", # using the latest version instead of the best version on Sept 18, 2018. Revisit as other versions become available
                            gbd_round_id = 5,
                            cause_id = 298,
                            measure_id = 1,
                            metric_id = 1,
                            age_group_id = 24,
                            location_id = loc.table$location_id[which(loc.table$ihme_loc_id == loc)],
                            year_id = c(1990:2017))
    comp_mort <- comp_mort[,c("year_id", "val")]
    names(comp_mort) <- c("year", "mean")
    rfs_mort_gbd <- merge(comp_mort, nat_mort[ ,c("year", "mean")], by = c("year"))
    names(rfs_mort_gbd) <- c("year","gbd_model", "lbd_model")
    rfs_mort_gbd$gbd_rf <- (rfs_mort_gbd$gbd_model) / rfs_mort_gbd$lbd_model
    rfs_mort_gbd$gbd_rf[which(is.na(rfs_mort_gbd$gbd_rf))] <- 1
  }






  ##################################################################################################################################


  if (t == "GBD2019") {
    source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')

    if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
      comp_mort <- get_outputs(topic = "cause",
                              version = 'latest',
                              gbd_round_id = 6,
                              cause_id = 298,
                              measure_id = 1,
                              metric_id = 1,
                              age_group_id = 24,
                              location_id = loc.table$location_id[which(loc.table$ihme_loc_id == loc)],
                              year_id = c(1990:2019),
                              decomp_step = "step5")
      comp_mort <- comp_mort[,c("year_id", "val")]
      names(comp_mort) <- c("year", "mean")
      rfs_mort_gbd <- merge(comp_mort, nat_mort[ ,c("year", "mean")], by = c("year"))
      names(rfs_mort_gbd) <- c("year","gbd_model", "lbd_model")
      rfs_mort_gbd$gbd_rf <- (rfs_mort_gbd$gbd_model) / rfs_mort_gbd$lbd_model
      rfs_mort_gbd$gbd_rf[which(is.na(rfs_mort_gbd$gbd_rf))] <- 1
    } else {

      connector <- get_gbd_locs(rake_subnational = T,
                                reg = loc,
                                shapefile_version = shp_version)
      nat_mort <- merge(nat_mort, connector, by = c("ADM1_CODE"))
      comp_mort <- get_outputs(topic = "cause",
                              version = 'latest',
                              gbd_round_id = 6,
                              cause_id = 298,
                              measure_id = 1,
                              metric_id = 1,
                              age_group_id = 24,
                              location_id = unique(nat_mort$location_id),
                              year_id = c(1990:2019),
                              decomp_step = "step5")
      comp_mort <- comp_mort[,c("location_id","year_id", "val")]
      names(comp_mort) <- c("location_id","year", "mean")
      rfs_mort_gbd <- merge(comp_mort, nat_mort[ ,c("location_id", "ADM1_CODE","year", "mean")], by = c("location_id", "year"))
      names(rfs_mort_gbd) <- c("location_id","year","gbd_model","ADM1_CODE","lbd_model")
      rfs_mort_gbd$gbd_rf <- (rfs_mort_gbd$gbd_model) / rfs_mort_gbd$lbd_model
      rfs_mort_gbd$gbd_rf[which(is.na(rfs_mort_gbd$gbd_rf))] <- 1
    }



  }


  ####################################################################################################################################

  if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
    nat_mort <- merge(nat_mort, rfs_mort_gbd, by = c("year"))
    nat_mort[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
    nat_mort$unraked_mean <- nat_mort$mean
    nat_mort$mean <- NULL
    nat_mort[,raked_mean := rowMeans(.SD) , .SD = overs]
    nat_mort$lbd_model <- NULL


    adm2_mort <- fread(paste0(input.dir, "/", loc,"_internal_raked_adm2_mort.csv"))
    adm2_mort <- merge(adm2_mort, rfs_mort_gbd, by = c("year"))
    adm2_mort[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
    adm2_mort$lbd_model <- NULL
    adm2_mort$gbd_model <- NULL
    adm2_mort[,double_raked_mean := rowMeans(.SD) , .SD = overs]


    valid_mort <- adm2_mort

    comp_mort <- valid_mort[, lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("year")]
    names(comp_mort) <- c("year", overs, "pop")
    comp_mort[,mean := rowMeans(.SD) , .SD = overs]
    valid_mort <- merge(nat_mort[ ,c("year", "raked_mean")], comp_mort[ ,c("year", "mean")], by = c("year"))
    valid_mort$test <- valid_mort$raked_mean / valid_mort$mean
    valid_mort$test[which(is.na(valid_mort$test))] <- 1
    if (max(valid_mort$test) == Inf) {
      message( "there is a /0 somewhere, please troubleshoot.")}
    valid_mort$test[which(valid_mort$test == Inf)] <- 1


    if (as.integer(sum(valid_mort$test)) != as.integer(NROW(valid_mort))) {
      message(paste0("error in ",t," raking")) } else {
        write.csv(rfs_mort_gbd, paste0(input.dir, "/",t,"_mort_rfs.csv"), row.names = FALSE)
        write.csv(adm2_mort, paste0(input.dir, "/", loc,"_",t,"_raked_adm2_mort.csv"), row.names = FALSE)
        write.csv(nat_mort, paste0(input.dir, "/", loc,"_",t,"_raked_nat_mort.csv"), row.names = FALSE)
      }

  } else {

    adm2_mort <- fread(paste0(input.dir,"/",loc,"_unraked_adm2_mort.csv"))
    adm2_mort <- merge(adm2_mort, sp_h, by.x = c("subnat"), by.y = c("ADM2_CODE"))
    adm2_mort <- merge(adm2_mort, rfs_mort_gbd, by = c("year", "ADM1_CODE"))
    adm2_mort[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
    adm2_mort$lbd_model <- NULL
    adm2_mort$gbd_model <- NULL
    adm2_mort[,double_raked_mean := rowMeans(.SD) , .SD = overs]
    adm2_mort <- as.data.frame(adm2_mort)
    adm2_mort <- adm2_mort[,c("year","subnat",overs,"pop","gbd_rf","double_raked_mean")]

    write.csv(rfs_mort_gbd, paste0(input.dir, "/",t,"_mort_rfs.csv"), row.names = FALSE)
    write.csv(adm2_mort, paste0(input.dir, "/", loc,"_",t,"_raked_adm2_mort.csv"), row.names = FALSE)
  }

}

