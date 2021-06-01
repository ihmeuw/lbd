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

core_repo  <- paste0("<<<< FILEPATH REDACTED >>>>")
indic_repo <- paste0("<<<< FILEPATH REDACTED >>>>")

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

new_pkg_lib   <- paste0("<<<< FILEPATH REDACTED >>>>","/r_packages/")
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

countries <- c("CIV")

###########
## read in the preliminary data sheet from UNAIDS that has compleete counts time series.
###########

UNAIDS <- readxl::read_excel("<<<< FILEPATH REDACTED >>>>", sheet = 4)

UNAIDS_CIV <- UNAIDS[which(UNAIDS$ISO3 %in% grep("CIV",UNAIDS$ISO3, value = T)), ]
UNAIDS_CIV <- UNAIDS_CIV[which(UNAIDS_CIV$...1 != "J- Total number receiving ART (0-14) - (Dec 31) Male+Female - (Dec 31)"), ]
UNAIDS_CIV <- UNAIDS_CIV[which(UNAIDS_CIV$ISO3 != "CIV"), ]

UNAIDS_CIV$sex_id <- NA

UNAIDS_CIV$sex_id[which(UNAIDS_CIV$...1 == "J- Total number receiving ART (15+) - (Dec 31) Male - (Dec 31)")] <- 1
UNAIDS_CIV$sex_id[which(UNAIDS_CIV$...1 == "J- Total number receiving ART (15+) - (Dec 31) Female - (Dec 31)")] <- 2

UNAIDS_CIV$...1 <- NULL
UNAIDS_CIV$name <- UNAIDS_CIV$...3
UNAIDS_CIV$...3 <- NULL
UNAIDS_CIV$ISO3 <- NULL

data_long <- reshape2::melt(UNAIDS_CIV, id.vars = c("name", "sex_id"))


data_long$name <- substr(data_long$name, 22, nchar(data_long$name))


data_long <- as.data.table(data_long)

data_long <- data_long[,lapply("value", function(x) sum(get(x), na.rm = TRUE)), by = c("variable", "name")]
names(data_long) <- c("year", "name", "prov_tot")




## Export

write.csv(data_long, paste0("<<<< FILEPATH REDACTED >>>>"))
