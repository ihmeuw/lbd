
#### Prepping UNAIDS data for later use



rm(list = ls())

#########
## Set Up
#########
setwd("~")
rm(list = ls())
# Set up
'%!in%' <- function(x,y)!('%in%'(x,y))
## Set repo
commondir      <- sprintf('<<<< FILEPATH REDACTED >>>>/common_inputs')
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>/", "/lbd_hiv/")

setwd(core_repo)

message('Loading in required R packages and MBG functions')
source(paste0(core_repo, '/mbg_central/setup.R'))
mbg_setup(package_list = package_list, repos = core_repo)

library(dplyr)
library(data.table)
library(lme4)
library(sf)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtools)
require(nnet)
library(spdep)
library(INLA)
library(stringr)

new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>", "/r_packages/")
.libPaths(new_pkg_lib)
test_pkg_list <- c('compositions', 'dplyr')
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(compositions)


'%!in%' <- function(x,y)!('%in%'(x,y))



#####
## Determine what countries we have ART data for
#####

countries <- list.dirs(paste0("<<<< FILEPATH REDACTED >>>>/"), full.names = F, recursive = F)

countries <- countries[which(countries %!in% c("HTI", "diag_plots"))]
#######
## set the HIV prevalence extraction to use for modeling ART coverage.
#######
prev_date <- "2020_04_10"
wp_version <- "wp1"

prev_info <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
mbg_run <- as.character(prev_info$Value[which(prev_info$Setting == "mbg_run")])

shp_date <- as.character(prev_info$Value[which(prev_info$Setting == "shapefile_version")])

plhiv <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))
plhiv$ADM1_CODE <- as.character(plhiv$ADM1_CODE)
##_read in geographic hierarchy
ad1 <- st_read(paste0("<<<< FILEPATH REDACTED >>>>"))
ad1$geometry <- NULL
ad1 <- ad1[ ,c("ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE")]
ad1$ADM1_CODE <- as.character(ad1$ADM1_CODE)




#####
## Identify the appropriate file to read in for each country
#####
files <- list()
plhiv_list <- list()

for (loc in countries) {
  print(loc)

  fs <- list.files(paste0("<<<< FILEPATH REDACTED >>>>/"), full.names = TRUE, recursive = FALSE)
  fs <- grep("province_unaids.csv", fs, value = TRUE)

  if (loc == "ZAF" | loc == "CMR") {
    fs <- list.files(paste0("<<<< FILEPATH REDACTED >>>>/"), full.names = FALSE, recursive = FALSE)
    fs <- grep(paste0("modeled_admin1_art_",wp_version), fs, value = TRUE)
    dates <- list()
    for (id in c(1:length(fs))) {
      d <- as.Date(str_replace_all(substr(fs[[id]], 1, 10), "-", "/"))
      l2 <- list()
      l2[["date"]] <- d
      l2[["file"]] <- fs[id]
      dates[[id]] <- l2
    }
    dates <- bind_rows(dates)
    fs <- paste0("<<<< FILEPATH REDACTED >>>>/")

    }








  files[[loc]] <- fread(fs)
  if (loc == "CMR") {
    files[[loc]] <- files[[loc]][,c("year", "ADM1_CODE", "modeled_prov_tot", "plhiv")]
    names(files[[loc]]) <- c("year", "name", "prov_tot", "mean")
    files[[loc]]$name <- as.character(files[[loc]]$name)

  } else{if (loc == "ZAF") {

    files[[loc]] <- files[[loc]][,c("year", "ADM1_CODE", "modeled_prov_tot", "plhiv")]
    files[[loc]]$ADM1_CODE <- as.character(files[[loc]]$ADM1_CODE)
    names(files[[loc]]) <- c("year", "name", "prov_tot", "mean")

  } else {
    lu <- fread(paste0("<<<< FILEPATH REDACTED >>>>/"))

    plhiv_lu <- merge(plhiv, lu, by = c("ADM1_NAME"))

    plhiv_lu <- plhiv_lu[,lapply("mean", function(x) sum(get(x), na.rm = TRUE)), by = c("year", "UNAIDS_NAME")]
    plhiv_lu$mean <- plhiv_lu$V1
    plhiv_lu$V1 <- NULL
    plhiv_list[[loc]] <- plhiv_lu

    files[[loc]] <- merge(files[[loc]], plhiv_lu, by.x = c("year", "name"), by.y = c("year", "UNAIDS_NAME"))
    files[[loc]]$name <- as.character(files[[loc]]$name)
    files[[loc]] <- files[[loc]][,c("year", "name", "prov_tot", "mean")]

  }}

}


#######
## create a look up table to connect location code to country name
#######
loc.table <- read.csv("<<<< FILEPATH REDACTED >>>>")

loc.table <- loc.table[which(loc.table$ihme_loc_id %in% countries), ]

c_lookup <- loc.table[ ,c("ihme_loc_id", "location_name", "location_id")]
c_lookup$loc <- c_lookup$ihme_loc_id
c_lookup$ADM0_NAME <- as.character(c_lookup$location_name)





###########
## combine extracted data into one data frame
###########

data_s <- bind_rows(files, .id = "ISO")

data_s$full_name <- paste0(data_s$ISO, "_", data_s$name)

data_s <- data_s[which(data_s$year < 2019), ] # 2019 data is not trusted at the moment

data_s$prov_cov <- data_s$prov_tot / data_s$mean


dir <- paste0("<<<< FILEPATH REDACTED >>>>")
dir.create(dir)


pdf(paste0(dir,Sys.Date(),"_admin1_ts.pdf"), width = 11, height = 8)
for (id in unique(data_s$full_name)) {

  df2 <- data_s[which(data_s$full_name == id), ]
  gg_count <- ggplot() +
    geom_point(data = df2, aes(x = year, y = prov_tot, color = "data"), size = 3) +
    scale_color_manual(values = c("#000000","#FF0000","#FFA500", "#0000FF", "#008000")) +
    xlab("year") +
    ylab("People on ART") +
    labs(title = paste0("# ART in ", id)) +
    theme(legend.position = "none", plot.title = element_text(size = 12)) +
    coord_cartesian(xlim = c(1995, 2018), ylim = c(0, max(df2$prov_tot)))


  gg_cov <- ggplot() +
    geom_point(data = df2, aes(x = year, y = prov_cov, color = "data"), size = 3) +
    scale_color_manual(values = c("#000000","#FF0000","#FFA500", "#0000FF", "#008000")) +
    xlab("year") +
    ylab("ART Coverage") +
    labs(title = "ART coverage") +
    theme(legend.position = c(-0.5,-0.20), legend.justification = c(0, 0), legend.direction = "horizontal", plot.title = element_text(size = 12)) +
    coord_cartesian(xlim = c(1995, 2018), ylim = c(0,1))

  grid.newpage()
  vp_1 <- viewport(width = 0.5, height = 0.85, x = 0.25, y = 0.575)
  vp_2 <- viewport(width = 0.5, height = 0.85, x = 0.75, y = 0.575)
  print(gg_count, vp = vp_1)
  print(gg_cov, vp = vp_2)
  print(id)
}
dev.off()





## Export

write.csv(data_s, paste0("<<<< FILEPATH REDACTED >>>>"))

