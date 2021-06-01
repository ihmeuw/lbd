################################################################################
## Purpose: Script to cache the gbd migration and mortality estimates needed to
##          get EPP to run.
## Run instructions: Run this in the collapse phase of EPP to prepare the input
##                   data needed to run EPP
################################################################################


## LBD Base Set Up #############################################################
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>","/r_packages/")
if (!dir.exists(new_pkg_lib)) {dir.create(new_pkg_lib)}


.libPaths(new_pkg_lib)
test_pkg_list <- c('rlang', "DBI","RMariaDB", "RMySQL")
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}

library(rlang, lib = new_pkg_lib)



## Set repo
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>/mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
#mbg_setup(package_list = package_list, repos = core_repo)


.libPaths(new_pkg_lib)
test_pkg_list <- c('assertable')
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}


library(assertable, lib = new_pkg_lib)



library(mortcore, lib = "mortality/shared/r")
library(mortdb, lib = "mortality/shared/r")


#Settings:
age_groups <- c(8:15)
gbd_year <- 2019
life_table_parameter_id <- 1 #this is the mortaltiy rate per 1
run_name <- "gbd_2019_april23_m"
dir.create(paste0("<<<< FILEPATH REDACTED >>>>",run_name,"/"))
dir.create(paste0("<<<< FILEPATH REDACTED >>>>",run_name,"/"))

gbd_q <- get_mort_outputs("no shock life table", "estimate", run_id = "best", life_table_parameter_id = life_table_parameter_id, age_group_id = age_groups, gbd_year = gbd_year)

saveRDS(gbd_q, paste0("<<<< FILEPATH REDACTED >>>>"))



gbd_m <- get_mort_outputs(model_name = "migration", model_type = "estimate", run_id = "best", gbd_year = gbd_year, age_group_id = age_groups)
gbd_m_2018 <- gbd_m[which(gbd_m$year_id == 2018), ]
gbd_m_2018$year_id <- 2019
gbd_m <- rbind(gbd_m, gbd_m_2018)


data <- gbd_m[which(gbd_m$year_id %in% c(1970:2019)), ]
data <- data[which(data$sex_id == 3), ]
data <- data[which(data$age_group_id %in% c(8:14)), ]



col <- data[,lapply("mean", function(x) sum(get(x))), by = c("year_id", "ihme_loc_id")]


write.csv(gbd_m, paste0("<<<< FILEPATH REDACTED >>>>"))
write.csv(col, paste0("<<<< FILEPATH REDACTED >>>>"))

