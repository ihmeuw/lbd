################################################################################
## Purpose:
## Date created:
## Date modified:

################################################################################

### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
code.dir <- "<<<< FILEPATH REDACTED >>>>"

## Packages
library(data.table)

# Arguments
 args <- commandArgs(trailingOnly = TRUE)
 if(length(args) > 0) {
  gbd.year <- args[1]
 } else {
  gbd.year <- 2019
 }

### Paths
epp.dt.dir <- paste0("<<<< FILEPATH REDACTED >>>>")
out.path <- paste0("<<<< FILEPATH REDACTED >>>>")

### Functions
library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")

### Tables
loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

### Code
file.list <- list.files(epp.dt.dir)
dt <- rbindlist(lapply(file.list, function(file) {
	path <- paste0(epp.dt.dir, file)
	temp.dt <- fread(path)
}))

write.csv(dt, out.path, row.names = F)

### End
