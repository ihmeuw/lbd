################################################################################
## Purpose: Create Summary Estimates at the various geographic levels.
################################################################################



## LBD Base Set Up
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>mbg/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("/<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("/<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

new_pkg_lib   <- paste0("/<<<< FILEPATH REDACTED >>>>", "/r_packages/")
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
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>", "lbd_hiv/epp/2_modeling/epp_geo/gbd/")

## Packages
library(data.table); library(ggplot2)
library(parallel); library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")

## Arguments
runs <- c("2020_05_18_art")

for (r in runs) {

  run.name <- r

config <- fread(paste0("<<<< FILEPATH REDACTED >>>>"), header = TRUE)

shapefile_version <- config$Value[which(config$Setting == "shapefile_version")]


locs <- list.dirs(paste0("<<<< FILEPATH REDACTED >>>>"), recursive = FALSE, full.names = FALSE)
locs <- locs[which(nchar(locs) == 3)]

data_list <- list()

for (loc in locs) {
  if (file.exists(paste0("<<<< FILEPATH REDACTED >>>>"))) {
    data <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    data_list[[loc]] <- data
  }
}

tot_inc <- ldply(data_list, data.frame)

for (loc in locs) {
  if (file.exists(paste0("<<<< FILEPATH REDACTED >>>>"))) {
    data <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    data_list[[loc]] <- data
  }
}

tot_mort <- ldply(data_list, data.frame)


for (loc in locs) {
  if (file.exists(paste0("<<<< FILEPATH REDACTED >>>>"))) {
    data <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    data_list[[loc]] <- data
  }
}

tot_prev <- ldply(data_list, data.frame)

for (loc in locs) {
  if (file.exists(paste0("<<<< FILEPATH REDACTED >>>>"))) {
    data <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    data_list[[loc]] <- data
  }
}

tot_art <- ldply(data_list, data.frame)




vars <- names(tot_inc)
overs <- grep("V", vars, value = TRUE)

sp_h <- st_read(paste0("/<<<< FILEPATH REDACTED >>>>"))
sp_h <- sp_h[ , c("ADM2_CODE", "ADM2_NAME", "ADM1_CODE", "ADM1_NAME", "ADM0_CODE", "ADM0_NAME")]
sp_h$geometry <- NULL

tot_mort <- as.data.table(merge(sp_h, tot_mort, by.x = c("ADM2_CODE"), by.y = c("subnat"), all.y = TRUE))
tot_inc <- as.data.table(merge(sp_h, tot_inc, by.x = c("ADM2_CODE"), by.y = c("subnat"), all.y = TRUE))
tot_prev <- as.data.table(merge(sp_h, tot_prev, by.x = c("ADM2_CODE"), by.y = c("subnat"), all.y = TRUE))
tot_art <- as.data.table(merge(sp_h, tot_art, by.x = c("ADM2_CODE"), by.y = c("subnat"), all.y = TRUE))

admin2_pop <- tot_mort[, c("ADM2_CODE", "ADM2_NAME", "ADM1_CODE", "ADM1_NAME", "ADM0_CODE", "ADM0_NAME", "year", "pop")]

admin2_inc <- tot_inc
admin2_inc$ADM0_CODE <- NULL
admin2_inc$ADM1_CODE <- NULL
admin2_mort <- tot_mort
admin2_mort$ADM0_CODE <- NULL
admin2_mort$ADM1_CODE <- NULL
admin2_prev <- tot_prev
admin2_prev$ADM0_CODE <- NULL
admin2_prev$ADM1_CODE <- NULL
admin2_art <- tot_art
admin2_art$ADM0_CODE <- NULL
admin2_art$ADM1_CODE <- NULL


admin2_inc <- admin2_inc[, (overs) := lapply(overs, function(x) ((get(x))*1000)) ]
admin2_mort <- admin2_mort[, (overs) := lapply(overs, function(x) ((get(x) / pop)*100000)) ]

admin2_inc_c <- as.data.table(as.data.frame(admin2_inc)) #there was a scoping issue needed to ensure an actuall copy
admin2_mort_c <- as.data.table(as.data.frame(admin2_mort)) #there was a scoping issue needed to ensure an actuall copy
admin2_prev_c <- as.data.table(as.data.frame(admin2_prev)) #there was a scoping issue needed to ensure an actuall copy

admin2_inc_c <- admin2_inc_c[, (overs) := lapply(overs, function(x) ((get(x)/100000)*pop)) ]
admin2_mort_c <- admin2_mort_c[, (overs) := lapply(overs, function(x) ((get(x)/100000)*pop)) ]
admin2_prev_c <- admin2_prev_c[, (overs) := lapply(overs, function(x) ((get(x)/100)*pop)) ]

tot_inc <- tot_inc[, (overs) := lapply(overs, function(x) get(x) * pop) ]
tot_prev <- tot_prev[, (overs) := lapply(overs, function(x) get(x) * pop) ]
tot_art <- tot_art[, (overs) := lapply(overs, function(x) get(x) * pop) ]


admin1_inc <- tot_inc[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM1_CODE", "year")]
admin1_prev <- tot_prev[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM1_CODE", "year")]
#admin1_art <- tot_art[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM1_CODE", "year")]
admin1_mort <- tot_mort[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM1_CODE", "year")]

names(admin1_inc) <- c("ADM1_CODE", "year", overs, "pop")
names(admin1_mort) <- c("ADM1_CODE", "year", overs, "pop")
names(admin1_prev) <- c("ADM1_CODE", "year", overs, "pop")
#names(admin1_art) <- c("ADM1_CODE", "year", overs, "pop")

admin1_inc <- admin1_inc[, (overs) := lapply(overs, function(x) ((get(x) / pop)*1000)) ]
admin1_mort <- admin1_mort[, (overs) := lapply(overs, function(x) ((get(x) / pop)*100000)) ]
admin1_prev <- admin1_prev[, (overs) := lapply(overs, function(x) ((get(x) / pop))) ]
#admin1_art <- admin1_art[, (overs) := lapply(overs, function(x) ((get(x) / pop))) ]


admin1_inc_c <- as.data.table(as.data.frame(admin1_inc)) #there was a scoping issue needed to ensure an actuall copy
admin1_mort_c <- as.data.table(as.data.frame(admin1_mort)) #there was a scoping issue needed to ensure an actuall copy
admin1_prev_c <- as.data.table(as.data.frame(admin1_prev)) #there was a scoping issue needed to ensure an actuall copy

admin1_inc_c <- admin1_inc_c[, (overs) := lapply(overs, function(x) ((get(x)/100000)*pop)) ]
admin1_mort_c <- admin1_mort_c[, (overs) := lapply(overs, function(x) ((get(x)/100000)*pop)) ]
admin1_prev_c <- admin1_prev_c[, (overs) := lapply(overs, function(x) ((get(x)/100)*pop)) ]

admin0_inc <- tot_inc[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM0_CODE", "year")]
admin0_mort <- tot_mort[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM0_CODE", "year")]
admin0_prev <- tot_prev[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM0_CODE", "year")]
#admin0_art <- tot_art[ ,lapply(c(overs, "pop"), function(x) sum(get(x), na.rm = T)), by = c("ADM0_CODE", "year")]

names(admin0_inc) <- c("ADM0_CODE", "year", overs, "pop")
names(admin0_mort) <- c("ADM0_CODE", "year", overs, "pop")
names(admin0_prev) <- c("ADM0_CODE", "year", overs, "pop")
#names(admin0_art) <- c("ADM0_CODE", "year", overs, "pop")

admin0_inc <- admin0_inc[, (overs) := lapply(overs, function(x) ((get(x) / pop)*1000)) ]
admin0_mort <- admin0_mort[, (overs) := lapply(overs, function(x) ((get(x) / pop)*100000)) ]
admin0_prev <- admin0_prev[, (overs) := lapply(overs, function(x) ((get(x) / pop))) ]
#admin0_art <- admin0_art[, (overs) := lapply(overs, function(x) ((get(x) / pop))) ]

admin0_inc_c <- as.data.table(as.data.frame(admin0_inc)) #there was a scoping issue needed to ensure an actuall copy
admin0_mort_c <- as.data.table(as.data.frame(admin0_mort)) #there was a scoping issue needed to ensure an actuall copy
admin0_prev_c <- as.data.table(as.data.frame(admin0_prev)) #there was a scoping issue needed to ensure an actuall copy


admin0_inc_c <- admin0_inc_c[, (overs) := lapply(overs, function(x) ((get(x)/100000)*pop)) ]
admin0_mort_c <- admin0_mort_c[, (overs) := lapply(overs, function(x) ((get(x)/100000)*pop)) ]
admin0_prev_c <- admin0_prev_c[, (overs) := lapply(overs, function(x) ((get(x)/100)*pop)) ]


admin0_inc_summary <- make_admin_pred_summary(admin0_inc[which(admin0_inc$year > 1999 & admin0_inc$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin0_mort_summary <- make_admin_pred_summary(admin0_mort[which(admin0_mort$year > 1999 & admin0_mort$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin0_prev_summary <- make_admin_pred_summary(admin0_prev[which(admin0_prev$year > 1999 & admin0_prev$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))


admin1_inc_summary <- make_admin_pred_summary(admin1_inc[which(admin1_inc$year > 1999 & admin1_inc$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin1_mort_summary <- make_admin_pred_summary(admin1_mort[which(admin1_mort$year > 1999 & admin1_mort$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin1_prev_summary <- make_admin_pred_summary(admin1_prev[which(admin1_prev$year > 1999 & admin1_prev$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))


admin2_inc_summary <- make_admin_pred_summary(admin2_inc[which(admin2_inc$year > 1999 & admin2_inc$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin2_mort_summary <- make_admin_pred_summary(admin2_mort[which(admin2_mort$year > 1999 & admin2_mort$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin2_prev_summary <- make_admin_pred_summary(admin2_prev[which(admin2_prev$year > 1999 & admin2_prev$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin2_art_summary <- make_admin_pred_summary(admin2_art[which(admin2_art$year > 1999 & admin2_art$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))



admin0_inc_c_summary <- make_admin_pred_summary(admin0_inc_c[which(admin0_inc_c$year > 1999 & admin0_inc_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin0_mort_c_summary <- make_admin_pred_summary(admin0_mort_c[which(admin0_mort_c$year > 1999 & admin0_mort_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin0_prev_c_summary <- make_admin_pred_summary(admin0_prev_c[which(admin0_prev_c$year > 1999 & admin0_prev_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))


admin1_inc_c_summary <- make_admin_pred_summary(admin1_inc_c[which(admin1_inc_c$year > 1999 & admin1_inc_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin1_mort_c_summary <- make_admin_pred_summary(admin1_mort_c[which(admin1_mort_c$year > 1999 & admin1_mort_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin1_prev_c_summary <- make_admin_pred_summary(admin1_prev_c[which(admin1_prev_c$year > 1999 & admin1_prev_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))


admin2_inc_c_summary <- make_admin_pred_summary(admin2_inc_c[which(admin2_inc_c$year > 1999 & admin2_inc_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin2_prev_c_summary <- make_admin_pred_summary(admin2_prev_c[which(admin2_prev_c$year > 1999 & admin2_prev_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))
admin2_mort_c_summary <- make_admin_pred_summary(admin2_mort_c[which(admin2_mort_c$year > 1999 & admin2_mort_c$year < 2019), ], sp_h, summary_stats = c("mean", "lower", "upper"))



dir <- paste0("<<<< FILEPATH REDACTED >>>>",run.name,"/admin_summaries")
dir.create(dir)

write.csv(admin2_pop, file = paste0(dir,"/admin2_pop.csv"), row.names = F)

write.csv(admin0_inc_summary, file = paste0(dir,"/admin0_inc_GBD2019_raked.csv"), row.names = F)
write.csv(admin0_mort_summary, file = paste0(dir,"/admin0_mort_GBD2019_raked.csv"), row.names = F)
write.csv(admin0_prev_summary, file = paste0(dir,"/admin0_prev_GBD2019_raked.csv"), row.names = F)
write.csv(admin1_inc_summary, file = paste0(dir,"/admin1_inc_GBD2019_raked.csv"), row.names = F)
write.csv(admin1_mort_summary, file = paste0(dir,"/admin1_mort_GBD2019_raked.csv"), row.names = F)
write.csv(admin1_prev_summary, file = paste0(dir,"/admin1_prev_GBD2019_raked.csv"), row.names = F)
write.csv(admin2_inc_summary, file = paste0(dir,"/admin2_inc_GBD2019_raked.csv"), row.names = F)
write.csv(admin2_mort_summary, file = paste0(dir,"/admin2_mort_GBD2019_raked.csv"), row.names = F)
write.csv(admin2_prev_summary, file = paste0(dir,"/admin2_prev_GBD2019_raked.csv"), row.names = F)
write.csv(admin2_art_summary, file = paste0(dir,"/admin2_art.csv"), row.names = F)



write.csv(admin0_inc, file = paste0(dir,"/admin0_inc_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin0_mort, file = paste0(dir,"/admin0_mort_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin0_prev, file = paste0(dir,"/admin0_prev_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin1_inc, file = paste0(dir,"/admin1_inc_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin1_mort, file = paste0(dir,"/admin1_mort_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin1_prev, file = paste0(dir,"/admin1_prev_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin2_inc, file = paste0(dir,"/admin2_inc_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin2_mort, file = paste0(dir,"/admin2_mort_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin2_prev, file = paste0(dir,"/admin2_prev_GBD2019_raked_draws.csv"), row.names = F)




write.csv(admin0_inc_c_summary, file = paste0(dir,"/admin0_inc_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin0_mort_c_summary, file = paste0(dir,"/admin0_mort_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin0_prev_c_summary, file = paste0(dir,"/admin0_prev_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin1_inc_c_summary, file = paste0(dir,"/admin1_inc_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin1_mort_c_summary, file = paste0(dir,"/admin1_mort_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin1_prev_c_summary, file = paste0(dir,"/admin1_prev_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin2_inc_c_summary, file = paste0(dir,"/admin2_inc_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin2_mort_c_summary, file = paste0(dir,"/admin2_mort_c_GBD2019_raked.csv"), row.names = F)
write.csv(admin2_prev_c_summary, file = paste0(dir,"/admin2_prev_c_GBD2019_raked.csv"), row.names = F)




write.csv(admin0_inc_c, file = paste0(dir,"/admin0_inc_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin0_mort_c, file = paste0(dir,"/admin0_mort_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin0_prev_c, file = paste0(dir,"/admin0_prev_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin1_inc_c, file = paste0(dir,"/admin1_inc_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin1_mort_c, file = paste0(dir,"/admin1_mort_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin1_prev_c, file = paste0(dir,"/admin1_prev_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin2_inc_c, file = paste0(dir,"/admin2_inc_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin2_mort_c, file = paste0(dir,"/admin2_mort_c_GBD2019_raked_draws.csv"), row.names = F)
write.csv(admin2_prev_c, file = paste0(dir,"/admin2_prev_c_GBD2019_raked_draws.csv"), row.names = F)

}


