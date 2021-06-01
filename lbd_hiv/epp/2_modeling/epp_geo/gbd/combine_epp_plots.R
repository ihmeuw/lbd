################################################################################
## Purpose: Combine location specific PDF's

################################################################################

### Setup
rm(list=ls())
windows <- Sys.info()[1]=="Windows"
code.dir <- paste0(ifelse(windows, "<<<< FILEPATH REDACTED>>>>"), "/HIV/")

## Packages
library(data.table)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
	run.name <- args[1]
} else {
	run.name <- "gbd17/test"
}
run.name <- paste0("gbd17/", run.name)
### Paths

### Functions
# source(paste0(code.dir, "shared_functions/get_locations.R"))

### Tables
# loc.table <- fread(paste0("<<<< FILEPATH REDACTED>>>>"))

### Code
if(grepl("/", run.name)) {
	file.name <- strsplit(run.name, "/")[[1]][2]
} else {
	file.name <- run.name
}

plot.dir <- paste0("<<<< FILEPATH REDACTED>>>>")
setwd(paste0(plot.dir, "locations/"))
# Combine location-specific plots
system(paste0("/usr/bin/ghostscript -dBATCH -dSAFER -DNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=", file.name, "_epp_plots.pdf -f *"))
# Move to parent directory
system(paste0("mv ", plot.dir, "locations/", file.name, "_epp_plots.pdf ", plot.dir, file.name, "_epp_plots.pdf"))
# Delete location specific plots
system(paste0("rm -r -f ", plot.dir, "locations/"))

### End
