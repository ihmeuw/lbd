################################################################################
## Purpose: This caches the populations used in the large geography models of
##          epp.  It gets them form the get_population() central function and
##          stores them in a run specific folder.
## Run instructions: Run from the main epp lauch script.
################################################################################

## Packages
if (!"data.table" %in% (.packages())) {
	library(data.table)
}

### Tables
if (!"loc.table" %in% ls()) {
	windows <- Sys.info()[1][["sysname"]] == "Windows"
	source(paste0("functions/get_locations.r"))
	loc.table <- fread("<<<< FILEPATH REDACTED >>>>")
}

### Functions
if (!"get_population" %in% ls()) {
	source(paste0("<<<< FILEPATH REDACTED >>>>/get_population.R"))
}

### Code
cache_population <- function(location_id  = 1, age_group_id = 21, dir = NULL, gbd_round = 6, decomp_step = "step4") {
	if (is.null(dir)) {
		stop("directory must be specified")
	}
	dir.create(dir, recursive = T, showWarnings = F)
	in.pop <- get_population(age_group_id = age_group_id,
	                         location_id = location_id,
	                         gbd_round_id = gbd_round,
	                         year_id = -1,
	                         sex_id = 1:2,
	                         location_set_id = 21,
	                         status = "best",
	                         decomp_step = decomp_step)

	invisible(lapply(location_id, function(loc.id) {
		out.pop <- copy(in.pop[location_id == loc.id])
		write.csv(out.pop, paste0(dir, loc.id, ".csv"), row.names = F)
	}))
}

### End
