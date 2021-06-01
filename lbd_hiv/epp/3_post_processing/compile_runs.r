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

new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>", "/r_packages/")
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


ecov_run <- "2020_05_13_ecov"
art_run <- "2020_05_18_art"


runs <- c(ecov_run,art_run)
levels <- c(0,1,2)
types <- c(".csv", "_draws.csv")

output_dir <- "<<<< FILEPATH REDACTED >>>>"
dir.create(output_dir)

## Summaries of Rates

for (l in c(1:3)) {

ad2_1_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

ad2_1_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

ad2_1_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

if (levels[[l]] == 2) {
ad2_1_art <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
ad2_art <- rbind(ad2_1_art)
write.csv(ad2_art, paste0(output_dir, "/admin",levels[[l]],"_art_ecov.csv"), row.names = F)
}



ad2_inc <- rbind(ad2_1_inc)
ad2_mort <- rbind(ad2_1_mort)
ad2_prev <- rbind(ad2_1_prev)

write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_raked_ecov.csv"), row.names = F)
write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_raked_ecov.csv"), row.names = F)
write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_raked_ecov.csv"), row.names = F)


}

for (l in c(1:3)) {

  ad2_1_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_1_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_1_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))


  if (levels[[l]] == 2) {
  ad2_1_art <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_art <- rbind(ad2_1_art)
  write.csv(ad2_art, paste0(output_dir, "/admin",levels[[l]],"_art_art.csv"), row.names = F)
  }

  ad2_inc <- rbind(ad2_1_inc)
  ad2_mort <- rbind(ad2_1_mort)
  ad2_prev <- rbind(ad2_1_prev)


  write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_raked_art.csv"), row.names = F)
  write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_raked_art.csv"), row.names = F)
  write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_raked_art.csv"), row.names = F)

}

for (l in c(1:3)) {

  ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM0_NAME %in% ad2_4_inc$ADM0_NAME), ]

  ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM0_NAME %in% ad2_4_mort$ADM0_NAME), ]

  ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM0_NAME %in% ad2_4_prev$ADM0_NAME), ]


  if (levels[[l]] == 2) {
  ad2_2_art <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_art <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_art <- ad2_2_art[!which(ad2_2_art$ADM0_NAME %in% ad2_4_art$ADM0_NAME), ]
  ad2_art <- rbind(ad2_2_art, ad2_4_art)
  write.csv(ad2_art, paste0(output_dir, "/admin",levels[[l]],"_art_final.csv"), row.names = F)
  }

  ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
  ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
  ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

  write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_raked_final.csv"), row.names = F)
  write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_raked_final.csv"), row.names = F)
  write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_raked_final.csv"), row.names = F)

}

## Summaries of counts

for (l in c(1:3)) {

  ad2_1_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_1_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_1_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_inc <- rbind(ad2_1_inc)
  ad2_mort <- rbind(ad2_1_mort)
  ad2_prev <- rbind(ad2_1_prev)

  write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_c_raked_ecov.csv"), row.names = F)
  write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_c_raked_ecov.csv"), row.names = F)
  write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_c_raked_ecov.csv"), row.names = F)


}

for (l in c(1:3)) {

  ad2_1_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_1_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_1_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

  ad2_inc <- rbind(ad2_1_inc)
  ad2_mort <- rbind(ad2_1_mort)
  ad2_prev <- rbind(ad2_1_prev)


  write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_c_raked_art.csv"), row.names = F)
  write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_c_raked_art.csv"), row.names = F)
  write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_c_raked_art.csv"), row.names = F)

}

for (l in c(1:3)) {
  ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM0_NAME %in% ad2_4_inc$ADM0_NAME), ]

  ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM0_NAME %in% ad2_4_mort$ADM0_NAME), ]

  ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM0_NAME %in% ad2_4_prev$ADM0_NAME), ]

  ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
  ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
  ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

  write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_c_raked_final.csv"), row.names = F)
  write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_c_raked_final.csv"), row.names = F)
  write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_c_raked_final.csv"), row.names = F)

}

for (l in c(1:3)) {
  ad2_2_inc_c_c <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_inc_c_c <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_inc_c_c <- ad2_2_inc_c_c[!which(ad2_2_inc_c_c$ADM0_NAME %in% ad2_4_inc_c_c$ADM0_NAME), ]

  ad2_2_inc_c_r <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_inc_c_r <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_inc_c_r <- ad2_2_inc_c_r[!which(ad2_2_inc_c_r$ADM0_NAME %in% ad2_4_inc_c_r$ADM0_NAME), ]

  ad2_2_mort_c_c <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_mort_c_c <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_mort_c_c <- ad2_2_mort_c_c[!which(ad2_2_mort_c_c$ADM0_NAME %in% ad2_4_mort_c_c$ADM0_NAME), ]

  ad2_2_mort_c_r <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_mort_c_r <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_mort_c_r <- ad2_2_mort_c_r[!which(ad2_2_mort_c_r$ADM0_NAME %in% ad2_4_mort_c_r$ADM0_NAME), ]

  ad2_2_ipr <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_ipr <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_ipr <- ad2_2_ipr[!which(ad2_2_ipr$ADM0_NAME %in% ad2_4_ipr$ADM0_NAME), ]

  ad2_inc_c_c <- rbind(ad2_2_inc_c_c, ad2_4_inc_c_c)
  ad2_inc_c_r <- rbind(ad2_2_inc_c_r, ad2_4_inc_c_r)
  ad2_mort_c_c <- rbind(ad2_2_mort_c_c, ad2_4_mort_c_c)
  ad2_mort_c_r <- rbind(ad2_2_mort_c_r, ad2_4_mort_c_r)
  ad2_ipr <- rbind(ad2_2_ipr, ad2_4_ipr)

  write.csv(ad2_inc_c_c, paste0(output_dir, "/admin",levels[[l]],"_inc_change_c_final.csv"), row.names = F)
  write.csv(ad2_inc_c_r, paste0(output_dir, "/admin",levels[[l]],"_inc_change_r_final.csv"), row.names = F)
  write.csv(ad2_mort_c_c, paste0(output_dir, "/admin",levels[[l]],"_mort_change_c_final.csv"), row.names = F)
  write.csv(ad2_mort_c_r, paste0(output_dir, "/admin",levels[[l]],"_mort_change_r_final.csv"), row.names = F)
  write.csv(ad2_ipr, paste0(output_dir, "/admin",levels[[l]],"_ipr_final.csv"), row.names = F)

}

######### Draw level rates

for (l in c(1:3)) {
  if (l == 1) {
  ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM0_CODE %in% ad2_4_inc$ADM0_CODE), ]

  ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM0_CODE %in% ad2_4_mort$ADM0_CODE), ]

  ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
  ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM0_CODE %in% ad2_4_prev$ADM0_CODE), ]

  ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
  ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
  ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

  write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_raked_final_draws.csv"), row.names = F)
  write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_raked_final_draws.csv"), row.names = F)
  write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_raked_final_draws.csv"), row.names = F)
  } else {if (l == 2) {
    ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM1_CODE %in% ad2_4_inc$ADM1_CODE), ]

    ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM1_CODE %in% ad2_4_mort$ADM1_CODE), ]

    ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM1_CODE %in% ad2_4_prev$ADM1_CODE), ]

    ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
    ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
    ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

    write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_raked_final_draws.csv"), row.names = F)
  } else {
    l <- 3
    ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM2_CODE %in% ad2_4_inc$ADM2_CODE), ]

    ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM2_CODE %in% ad2_4_mort$ADM2_CODE), ]

    ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM2_CODE %in% ad2_4_prev$ADM2_CODE), ]

    ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
    ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
    ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

    write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_raked_final_draws.csv"), row.names = F)


    }}

}


######### Draw level counts

for (l in c(1:3)) {
  if (l == 1) {
    ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM0_CODE %in% ad2_4_inc$ADM0_CODE), ]

    ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM0_CODE %in% ad2_4_mort$ADM0_CODE), ]

    ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM0_CODE %in% ad2_4_prev$ADM0_CODE), ]

    ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
    ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
    ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

    write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_c_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_c_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_c_raked_final_draws.csv"), row.names = F)
  } else {if (l == 2) {
    ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM1_CODE %in% ad2_4_inc$ADM1_CODE), ]

    ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM1_CODE %in% ad2_4_mort$ADM1_CODE), ]

    ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM1_CODE %in% ad2_4_prev$ADM1_CODE), ]

    ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
    ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
    ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

    write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_c_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_c_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_c_raked_final_draws.csv"), row.names = F)
  } else {
    l <- 3
    ad2_2_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_inc <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_inc <- ad2_2_inc[!which(ad2_2_inc$ADM2_CODE %in% ad2_4_inc$ADM2_CODE), ]

    ad2_2_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_mort <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_mort <- ad2_2_mort[!which(ad2_2_mort$ADM2_CODE %in% ad2_4_mort$ADM2_CODE), ]

    ad2_2_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_4_prev <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
    ad2_2_prev <- ad2_2_prev[!which(ad2_2_prev$ADM2_CODE %in% ad2_4_prev$ADM2_CODE), ]

    ad2_inc <- rbind(ad2_2_inc, ad2_4_inc)
    ad2_mort <- rbind(ad2_2_mort, ad2_4_mort)
    ad2_prev <- rbind(ad2_2_prev, ad2_4_prev)

    write.csv(ad2_inc, paste0(output_dir, "/admin",levels[[l]],"_inc_c_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_mort, paste0(output_dir, "/admin",levels[[l]],"_mort_c_raked_final_draws.csv"), row.names = F)
    write.csv(ad2_prev, paste0(output_dir, "/admin",levels[[l]],"_prev_c_raked_final_draws.csv"), row.names = F)


  }}

}


############ compile populations.

ad2_pops <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))


write.csv(ad2_pops, paste0(output_dir, "/admin2_pops.csv"), row.names = F)








