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
config <- read.csv(paste0("'<<<< FILEPATH REDACTED >>>>"))
shp_version <- config$Value[which(config$Setting == "shapefile_version")]

### Paths
input.dir <- paste0('<<<< FILEPATH REDACTED >>>>')


loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

for (t in gbd_rake_target) {
   if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
     nat_inc <- fread(paste0(input.dir,"/",loc,"_SPU_inc_draws.csv"))
     ndraws <- NCOL(nat_inc) - 1
     overs <- paste0("V", 1:ndraws)
     names(nat_inc) <- c("year", overs)
     nat_inc[,mean := rowMeans(.SD) , .SD = overs]
   } else {
     nat_inc <- fread(paste0(input.dir,"/",loc,"_unraked_adm2_inc.csv"))
     sp_h <- read.dbf(paste0("<<<< FILEPATH REDACTED >>>>"))
     nat_inc <- as.data.table(merge(sp_h, nat_inc, by.x = c("ADM2_CODE"), by.y = c("subnat")))
     overs <- names(nat_inc)
     overs <- overs[which(overs %!in% c("ADM2_CODE","NAME_0", "NAME_1", "NAME_2", "geo_id", "ad2_id", "ad0_parent", "ad1_parent", "ADM2_NAME","ADM1_CODE", "ADM1_NAME", "ADM0_CODE", "ADM0_NAME", "year", "pop"))]
     nat_inc[, (overs) := lapply(overs, function(x) get(x) * pop )]
     nat_inc <- nat_inc[,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM1_CODE","ADM1_NAME", "year")]
     names(nat_inc) <- c("ADM1_CODE","ADM1_NAME", "year", overs, "pop")
     nat_inc[, (overs) := lapply(overs, function(x) get(x) / pop )]
     nat_inc[,mean := rowMeans(.SD) , .SD = overs]
     nat_inc <- nat_inc[,c("ADM1_CODE","ADM1_NAME", "year", "mean", "pop")]

   }


if (t == "GBD2017") {
  source('/<<<< FILEPATH REDACTED >>>>/get_outputs.R')
  comp_inc <- get_outputs(topic = "cause",
                          version = "best", # using the latest version instead of the best version on Sept 18, 2018. Revisit as other versions become available
                          gbd_round_id = 5,
                          cause_id = 298,
                          measure_id = 6,
                          metric_id = 3,
                          age_group_id = 24,
                          location_id = loc.table$location_id[which(loc.table$<<<< FILEPATH REDACTED >>>>_loc_id == loc)],
                          year_id = c(1990:2017))
  comp_inc <- comp_inc[,c("year_id", "val")]
  names(comp_inc) <- c("year", "mean")
  rfs_inc_gbd <- merge(comp_inc, nat_inc[ ,c("year", "mean")], by = c("year"))
  names(rfs_inc_gbd) <- c("year","gbd_model", "lbd_model")
  rfs_inc_gbd$gbd_rf <- (rfs_inc_gbd$gbd_model * 100) / rfs_inc_gbd$lbd_model # the scalar of 100 is applied to the output of get output because EPP reports out incidence in cases per 100 (or at least I hope so)
  rfs_inc_gbd$gbd_rf[which(is.na(rfs_inc_gbd$gbd_rf))] <- 1
}






##################################################################################################################################


  if (t == "GBD2019") {
    source('/<<<< FILEPATH REDACTED >>>>/get_outputs.R')

    if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
    comp_inc <- get_outputs(topic = "cause",
                            version = 'latest',
                            gbd_round_id = 6,
                            cause_id = 298,
                            measure_id = 6,
                            metric_id = 3,
                            age_group_id = 24,
                            location_id = loc.table$location_id[which(loc.table$<<<< FILEPATH REDACTED >>>>_loc_id == loc)],
                            year_id = c(1990:2019),
                            decomp_step = "step5")
    comp_inc <- comp_inc[,c("year_id", "val")]
    names(comp_inc) <- c("year", "mean")
    rfs_inc_gbd <- merge(comp_inc, nat_inc[ ,c("year", "mean")], by = c("year"))
    names(rfs_inc_gbd) <- c("year","gbd_model", "lbd_model")
    rfs_inc_gbd$gbd_rf <- (rfs_inc_gbd$gbd_model * 100) / rfs_inc_gbd$lbd_model # the scalar of 100 is applied to the output of get output because EPP reports out incidence in cases per 100 (or at least I hope so)
    rfs_inc_gbd$gbd_rf[which(is.na(rfs_inc_gbd$gbd_rf))] <- 1
    } else {

      connector <- get_gbd_locs(rake_subnational = T,
                                reg = loc,
                                shapefile_version = shp_version)
      nat_inc <- merge(nat_inc, connector, by = c("ADM1_CODE"))
      comp_inc <- get_outputs(topic = "cause",
                              version = 'latest',
                              gbd_round_id = 6,
                              cause_id = 298,
                              measure_id = 6,
                              metric_id = 3,
                              age_group_id = 24,
                              location_id = unique(nat_inc$location_id),
                              year_id = c(1990:2019),
                              decomp_step = "step5")
      comp_inc <- comp_inc[,c("location_id","year_id", "val")]
      names(comp_inc) <- c("location_id","year", "mean")
      rfs_inc_gbd <- merge(comp_inc, nat_inc[ ,c("location_id", "ADM1_CODE","year", "mean")], by = c("location_id", "year"))
      names(rfs_inc_gbd) <- c("location_id","year","gbd_model","ADM1_CODE","lbd_model")
      rfs_inc_gbd$gbd_rf <- (rfs_inc_gbd$gbd_model * 100) / rfs_inc_gbd$lbd_model # the scalar of 100 is applied to the output of get output because EPP reports out incidence in cases per 100 (or at least I hope so)
      rfs_inc_gbd$gbd_rf[which(is.na(rfs_inc_gbd$gbd_rf))] <- 1
    }



  }


####################################################################################################################################

  if (loc %!in% c("NGA", "ETH", "KEN", "ZAF")) {
nat_inc <- merge(nat_inc, rfs_inc_gbd, by = c("year"))
nat_inc[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
nat_inc$unraked_mean <- nat_inc$mean
nat_inc$mean <- NULL
nat_inc[,raked_mean := rowMeans(.SD) , .SD = overs]
nat_inc$lbd_model <- NULL


adm2_inc <- fread(paste0(input.dir, "/", loc,"_internal_raked_adm2_inc.csv"))
adm2_inc <- merge(adm2_inc, rfs_inc_gbd, by = c("year"))
adm2_inc[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
adm2_inc$lbd_model <- NULL
adm2_inc$gbd_model <- NULL
adm2_inc[,double_raked_mean := rowMeans(.SD) , .SD = overs]


valid_inc <- adm2_inc
valid_inc[,(overs) := lapply(overs, function(x) get(x) * pop)]

comp_inc <- valid_inc[, lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("year")]
names(comp_inc) <- c("year", overs, "pop")
comp_inc[,(overs) := lapply(overs, function(x) get(x) / pop)]
comp_inc[,mean := rowMeans(.SD) , .SD = overs]
valid_inc <- merge(nat_inc[ ,c("year", "raked_mean")], comp_inc[ ,c("year", "mean")], by = c("year"))
valid_inc$test <- valid_inc$raked_mean / valid_inc$mean
valid_inc$test[which(is.na(valid_inc$test))] <- 1
if (max(valid_inc$test) == Inf) {
  message( "there is a /0 causing issues somewhere, troubleshoot needed.")}
valid_inc$test[which(valid_inc$test == Inf)] <- 1

adm2_inc[,(overs) := lapply(overs, function(x) get(x) / pop)]

if (as.integer(sum(valid_inc$test)) != as.integer(NROW(valid_inc))) {
  message(paste0("error in ",t," raking")) } else {
  write.csv(rfs_inc_gbd, paste0(input.dir, "/",t,"_inc_rfs.csv"), row.names = FALSE)
  write.csv(adm2_inc, paste0(input.dir, "/", loc,"_",t,"_raked_adm2_inc.csv"), row.names = FALSE)
  write.csv(nat_inc, paste0(input.dir, "/", loc,"_",t,"_raked_nat_inc.csv"), row.names = FALSE)
}

  } else {

    adm2_inc <- fread(paste0(input.dir,"/",loc,"_unraked_adm2_inc.csv"))
    adm2_inc <- merge(adm2_inc, sp_h, by.x = c("subnat"), by.y = c("ADM2_CODE"))
    adm2_inc <- merge(adm2_inc, rfs_inc_gbd, by = c("year", "ADM1_CODE"))
    adm2_inc[,(overs) := lapply(overs, function(x) get(x) * gbd_rf)]
    adm2_inc$lbd_model <- NULL
    adm2_inc$gbd_model <- NULL
    adm2_inc[,double_raked_mean := rowMeans(.SD) , .SD = overs]
    adm2_inc <- as.data.frame(adm2_inc)
    adm2_inc <- adm2_inc[,c("year","subnat",overs,"pop","gbd_rf","double_raked_mean")]

    write.csv(rfs_inc_gbd, paste0(input.dir, "/",t,"_inc_rfs.csv"), row.names = FALSE)
    write.csv(adm2_inc, paste0(input.dir, "/", loc,"_",t,"_raked_adm2_inc.csv"), row.names = FALSE)
}

}
