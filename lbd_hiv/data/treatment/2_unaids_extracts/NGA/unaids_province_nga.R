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
commondir      <- sprintf("<<<< FILEPATH REDACTED >>>>")
package_list <- c(t(read.csv(sprintf('%s/package_list.csv',commondir), header = FALSE)))

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_core/")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/")

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
## set country
#####

countries <- c("NGA")
#######
## set the HIV prevalence extraction to use for modeling ART coverage.
#######
prev_date <- "2020_04_10" #string that indicates where to save the desired outputs and diagnostic plots.

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

AD1_NGA <- ad1[which(ad1$ADM0_NAME == "Nigeria"), ]



#######
## create a look up table to connect location code to country name
#######
loc.table <- read.csv("<<<< FILEPATH REDACTED >>>>")

locs_list <- grep("NGA", loc.table$ihme_loc_id, value = T)

loc.table <- loc.table[which(loc.table$ihme_loc_id %in% locs_list), ]

c_lookup <- loc.table[ ,c("ihme_loc_id", "location_name", "location_id")]
c_lookup$loc <- c_lookup$ihme_loc_id
c_lookup$ADM0_NAME <- as.character(c_lookup$location_name)





###########
## read in the preliminary data sheet from UNAIDS that has compleete counts time series.
###########

UNAIDS <- readxl::read_excel("<<<< FILEPATH REDACTED >>>>", sheet = 4)

UNAIDS_NGA <- UNAIDS[which(UNAIDS$ISO3 %in% grep("NGA",UNAIDS$ISO3, value = T)), ]
UNAIDS_NGA <- UNAIDS_NGA[which(UNAIDS_NGA$...1 != "J- Total number receiving ART (0-14) - (Dec 31) Male+Female - (Dec 31)"), ]
UNAIDS_NGA <- UNAIDS_NGA[which(UNAIDS_NGA$ISO3 != "NGA"), ]

UNAIDS_NGA$sex_id <- NA

UNAIDS_NGA$sex_id[which(UNAIDS_NGA$...1 == "J- Total number receiving ART (15+) - (Dec 31) Male - (Dec 31)")] <- 1
UNAIDS_NGA$sex_id[which(UNAIDS_NGA$...1 == "J- Total number receiving ART (15+) - (Dec 31) Female - (Dec 31)")] <- 2

UNAIDS_NGA$...1 <- NULL
UNAIDS_NGA$name <- UNAIDS_NGA$...3
UNAIDS_NGA$...3 <- NULL
UNAIDS_NGA$ISO3 <- NULL

data_long <- reshape2::melt(UNAIDS_NGA, id.vars = c("name", "sex_id"))
data_long$name <- substr(data_long$name, 17, nchar(data_long$name))

data_long$name[which(data_long$name == "Akwa-Ibom")] <- "Akwa Ibom"
data_long$name[which(data_long$name == "Cross-River")] <- "Cross River"
data_long$name[which(data_long$name == "FCT")] <- "Federal Capital Territory"
data_long$name[which(data_long$name == "Nasarawa")] <- "Nassarawa"
data_long$name[which(data_long$name == "Tarabe")] <- "Taraba"

data_long <- as.data.table(data_long)

data_long <- data_long[,lapply("value", function(x) sum(get(x), na.rm = TRUE)), by = c("variable", "name")]
names(data_long) <- c("year", "name", "prov_tot")




## Export

write.csv(data_long, paste0("<<<< FILEPATH REDACTED >>>>"))
