## LBD Base Set Up
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

### GBD addititions Setup


windows <- Sys.info()[1][["sysname"]] == "Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>/epp/2_modeling/epp_geo/")

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

loc <- "BEN"

y <- "2019"


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

for (i in 1:length(pjnz_files)) {
  f <- pjnz_files[[i]]
  breaks <- strsplit(f, "_")[[1]]
  pjnz <- paste0("<<<< FILEPATH REDACTED >>>>",y,"/",loc,"/", f)
  input <- read_epp_input(pjnz)
  art <- input$epp.art
  art$name <- breaks[[2]] ## this is specific to the naming convention of the UNAIDS files and needs to be set manually
  input_list[[f]] <- art
}


dataset <- do.call(rbind, input_list)


dataset$prov_tot <- dataset$m.val + dataset$f.val
dataset <- dataset[which(dataset$f.isperc == "N"), ]

df <- dataset[ ,c("year", "name", "prov_tot")]

write.csv(df, paste0("<<<< FILEPATH REDACTED >>>>"), row.names = FALSE)

