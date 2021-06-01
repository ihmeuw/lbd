## LBD Base Set Up
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

### GBD addititions Setup


windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>/lbd_hiv/epp/2_modeling/epp_geo/")

## Packages
library(data.table); library(mvtnorm); library(survey)
new_pkg_lib   <- "<<<< FILEPATH REDACTED >>>>"
.libPaths(new_pkg_lib)
test_pkg_list <- c('slackr', "assertable")
for (pkg in test_pkg_list) {
  if (!pkg %in% as.vector(installed.packages(lib.loc = new_pkg_lib)[, "Package"])) {
    install.packages(pkg, lib = new_pkg_lib)
  }
}
library(slackr)
library(assertable)

## Set UNAIDS yr to search and loc

loc <- "ZAF"

y <- "2017"


### Functions
## GBD
source(paste0(code.dir,"gbd/prep_data_geo.R"))
source(paste0(code.dir,"gbd/prep_output_geo.R"))
source(paste0(code.dir,"gbd/data_sub_geo.R"))
source(paste0(code.dir,"gbd/plot_fit_geo.R"))
source(paste0(code.dir,"gbd/ind_data_prep.R"))
source(paste0(code.dir,"gbd/data_sub_subnat_geo.R"))
source(paste0(code.dir,"gbd/prep_data_admin2_geo.R"))


## EPP
source(paste0(code.dir,"R/epp_geo.R"))
source(paste0(code.dir,"R/fit-model_geo.R"))
source(paste0(code.dir,"R/generics.R"))
source(paste0(code.dir,"R/IMIS_geo.R"))
source(paste0(code.dir,"R/likelihood_geo.R"))
source(paste0(code.dir,"R/read-epp-files_geo.R"))

paste0("<<<< FILEPATH REDACTED >>>>",y,"/",loc,"/")

files <- list.files(paste0("<<<< FILEPATH REDACTED >>>>",y,"/",loc,"/"))
pjnz_files <- grep("PJNZ", files, value = TRUE)
input_list <- list()

c.year <- y

if (c.year == 2016 | c.year == 2017 | c.year == 2018) {
  dir <- paste0("<<<< FILEPATH REDACTED >>>>", c.year, "/ZAF/")
} else {
  dir <- paste0("<<<< FILEPATH REDACTED >>>>", c.year, "/", temp.loc, "/")
}


pjnz <- paste0(dir,pjnz_files)

eppd <- read_epp_data(pjnz, c.year = 2017)
epp.subp <- read_epp_subpops(pjnz, no.anc = no.anc)
epp.input <- read_epp_input(pjnz)

epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd, no.anc = no.anc)
dataset <- c()
for (s in names(epp.subp.input)) {
  df <- epp.subp.input[[s]]$epp.art
  dataset[[s]] <- df
}

dataset <- bind_rows(dataset, .id = "name")

dataset$prov_tot <- dataset$m.val + dataset$f.val
dataset <- dataset[which(dataset$f.isperc == "N"), ]

df <- dataset[ ,c("year", "name", "prov_tot")]

write.csv(df, paste0("<<<< FILEPATH REDACTED >>>>"), row.names = FALSE)
