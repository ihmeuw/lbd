################################################################################
## Purpose: Split India states into U/R using NFHS survey values

################################################################################

### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>", "/hiv_gbd2019/decomp_2019/")

## Packages
library(data.table)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
	run.name <- args[1]
} else {
	run.name <- "190124_decomp"
}


### Paths
dir.list <- list(inc = paste0("<<<< FILEPATH REDACTED >>>>"),
                 prev = paste0("<<<< FILEPATH REDACTED >>>>"))
suffix.list <- list(inc = "_SPU_inc_draws.csv", prev = "_SPU_prev_draws.csv")

### Functions
library(mortdb, lib = "<<<< FILEPATH REDACTED >>>>")

### Tables
loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

### Code
## India
file.list <- list.files(dir.list[[2]], "IND")
locs <- gsub(suffix.list[[2]], "", file.list)
ind.locs <- loc.table[grepl("IND", ihme_loc_id) & level == 4, ihme_loc_id]
missing.locs <- setdiff(ind.locs, locs)
for(measure in names(dir.list)) {
	# measure <- "prev"
	print(measure)
	bound.dt <- rbindlist(lapply(locs, function(ind.loc) {
		path <- paste0(dir.list[[measure]], ind.loc, suffix.list[[measure]])
		dt <- fread(path)
	}))
	mean.dt <- bound.dt[, lapply(.SD, mean), by = "year"]
	for(loc in missing.locs) {
		out.path <- paste0(dir.list[[measure]], loc, suffix.list[[measure]])
		write.csv(mean.dt, out.path, row.names = F)
	}
}

## Kenya
file.list <- list.files(dir.list[[2]], "KEN")
locs <- gsub(suffix.list[[2]], "", file.list)
ken.locs <- loc.table[grepl("KEN", ihme_loc_id) & level == 5, ihme_loc_id]
missing.locs <- setdiff(ken.locs, locs)
for(measure in names(dir.list)) {
	# measure <- "prev"
	print(measure)
	bound.dt <- rbindlist(lapply(locs, function(ken.loc) {
		path <- paste0(dir.list[[measure]], ken.loc, suffix.list[[measure]])
		dt <- fread(path)
	}))
	mean.dt <- bound.dt[, lapply(.SD, mean), by = "year"]
	for(loc in missing.locs) {
		out.path <- paste0(dir.list[[measure]], loc, suffix.list[[measure]])
		write.csv(mean.dt, out.path, row.names = F)
	}
}
