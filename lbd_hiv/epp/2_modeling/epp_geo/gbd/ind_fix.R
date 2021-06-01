################################################################################
## Purpose:
## Date created:
## Date modified:
## Run instructions:
## Notes:
################################################################################

### Setup
rm(list=ls())
windows <- Sys.info()[1][["sysname"]]=="Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>", "/HIV/epp-feature-reset/")

## Packages
library(data.table); library(mvtnorm); library(survey)

## Arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args) > 0) {
	run.name <- args[1]
	loc <- args[2]
	proj.end <- as.integer(args[3])
	i <- Sys.getenv("SGE_TASK_ID")
} else {
	# run.name <- paste0(substr(gsub("-","",Sys.Date()),3,8), "_ken_county")
	run.name <- paste0("171219_ken_ind")
	loc <- "IND_4844"
	proj.end <- 2017
	i <- 1
}


### Temp Arguments
start.year <- 1970
stop.year <- 2017
trans.params <- TRUE
stop.collapse <- FALSE
anc.prior <- TRUE
gbd.pop <- TRUE
no.anc <- FALSE

### Paths
input.dir <- paste0('<<<< FILEPATH REDACTED >>>>', run.name, "/", loc)
out.dir <- paste0('<<<< FILEPATH REDACTED >>>>')
dir.create(out.dir, recursive = T, showWarnings = F)
out.path <- paste0(out.dir, "/results", i, ".RData")
pdf.path <- paste0(out.dir, "/test_results", i, ".pdf")

### Functions
## GBD
source(paste0(code.dir,"gbd/prep_data.R"))
source(paste0(code.dir,"gbd/prep_output.R"))
source(paste0(code.dir,"gbd/data_sub.R"))
source(paste0(code.dir,"gbd/plot_fit.R"))
source(paste0(code.dir,"gbd/ind_data_prep.R"))


## EPP
source(paste0(code.dir,"R/epp.R"))
source(paste0(code.dir,"R/fit-model.R"))
source(paste0(code.dir,"R/generics.R"))
source(paste0(code.dir,"R/IMIS.R"))
source(paste0(code.dir,"R/likelihood.R"))
source(paste0(code.dir,"R/read-epp-files.R"))

## Model in C
if(trans.params) {
 # Load C version
	dyn.load(paste0(code.dir, "src/fnEPP", .Platform$dynlib.ext))  # Load C version with time series for transition parameters
} else {
	dyn.load(paste0(code.dir, "src/epp", .Platform$dynlib.ext))
}

source("<<<<FILEPATH REDACTED>>>>>/get_locations.r")

### Tables
loc.table <- fread(paste0("<<<< FILEPATH REDACTED >>>>"))

### Code
## Prep data and collapse location subpopulations
## Get location of file

unaids.year <- loc.table[ihme_loc_id == loc, unaids_recent]
if(unaids.year == 2017) {
    dir <- paste0("<<<< FILEPATH REDACTED >>>>")
} else {
    dir <- paste0("<<<< FILEPATH REDACTED >>>>")
}

filepath <- paste0(dir, loc)
ind_totals <- ind.prepare.epp.fit(filepath, proj.end = proj.end + 0.5)
eppd <- ind_totals$eppd.ind
epp.subp <- ind_totals$epp.subp.ind
epp.input <- ind_totals$epp.input.ind
## Fill in missing data (epp_totals vs. ind_totals)
epp.input$epidemic.start <- 1985
pop.dim <- length(epp.input$epp.pop$year)
# epp.pop
epp.input$epp.pop$cd4median <- rep(0, times = pop.dim)
epp.input$epp.pop$hivp15yr <- rep(0, times = pop.dim) # can pull from spectrum estimates (prevelance among 15 year olds)

# Expand progression parameters to matrix
prog.names <- c("cd4stage.dur", "cd4mort", "artmort.less6mos", "artmort.6to12mos", "artmort.after1yr")
for(param in prog.names) {
    # param <- prog.names[1]
    epp.input[[param]] <- matrix(rep(epp.input[[param]], times = 8), nrow = 8, byrow = T)
}

# cd4initperc
std.cd4initperc <- matrix(data = c(64.3, 35.7, 0, 0, 0, 0, 0, 60.7, 39.3, 0, 0, 0, 0, 0, 58.5, 41.5, 0, 0, 0, 0, 0, 55.2, 44.8, 0, 0, 0, 0, 0, 64.3, 35.7, 0, 0, 0, 0, 0, 60.7, 39.3, 0, 0, 0, 0, 0, 58.5, 41.5, 0, 0, 0, 0, 0, 55.2, 44.8, 0, 0, 0, 0, 0),
                            nrow = 8, ncol = 7, byrow = T)
row.names(std.cd4initperc) <- c("NEWINFECTSCD4_M_15_24", "NEWINFECTSCD4_M_25_34", "NEWINFECTSCD4_M_35_44", "NEWINFECTSCD4_M_45_54", "NEWINFECTSCD4_F_15_24", "NEWINFECTSCD4_F_25_34", "NEWINFECTSCD4_F_35_44", "NEWINFECTSCD4_F_45_54")
colnames(std.cd4initperc) <- paste0("V", 2:8)
epp.input$cd4initperc <- std.cd4initperc

# epp.art
art.dim <- nrow(epp.input$epp.art)
epp.input$epp.art$'1stto2ndline'  <- rep(0, times = art.dim)
epp.input$epp.art$art15yr  <- rep(0, times = art.dim)
epp.input$epp.art$artdropout  <- rep(0, times = art.dim)

# art.specpop
specpop.sub <- data.frame(specpop = c("PW", "TBHIV", "DC", "FSW", "MSM"),
                          percelig = c(0, 0, 0, 0, 0),
                          yearelig = c(2015, 2015, 2015, 2015, 2015))
epp.input$art.specpop <- specpop.sub

# hivp15yr.cd4dist
epp.input$hivp15yr.cd4dist <- c(0.056, 0.112, 0.112, 0.070, 0.140, 0.230, 0.280)

# art15yr.cd4dist
epp.input$art15yr.cd4dist <- c(0.00, 0.00, 0.11, 0.23, 0.23, 0.14, 0.29)

epp.subp.input <- fnCreateEPPSubpops(epp.input, epp.subp, eppd)
val <- setNames(vector("list", length(eppd)), names(eppd))
set.list.attr <- function(obj, attrib, value.lst) mapply(function(set, value) {
    attributes(set)[[attrib]] <- value
    set
}, obj, value.lst)
val <- set.list.attr(val, "eppd", eppd)
val <- set.list.attr(val, "likdat", lapply(eppd, fnCreateLikDat, anchor.year=epp.input$start.year))
val <- set.list.attr(val, "eppfp", lapply(epp.subp.input, fnCreateEPPFixPar, proj.end = proj.end + 0.5))
val <- set.list.attr(val, "country", attr(eppd, "country"))
val <- set.list.attr(val, "region", names(eppd))


load(paste0("<<<< FILEPATH REDACTED >>>>"))
names(result)

## Reset subpopulation proportions
for(subpop in names(result)) {
	# subpop <- names(result)[1]
	print(subpop)
	new.pop <- attr(val[[subpop]], "eppfp")$epp.pop.ts
	print(new.pop)
	# result[[subpop]]$fp$epp.pop.ts <- new.pop
}
result <- lapply(result, simfit.gbd)

## Aggregate subpopulations to national and write prevalence and incidence draws
years <- unique(floor(result[[1]]$fp$proj.steps))
nat.data <- nat.draws(result)
var_names <- sapply(1:ncol(nat.data$prev), function(a) {paste0('draw',a)})
out_data <- lapply(nat.data, data.frame)
for (n in c('prev', 'incid', 'art', 'art_num', 'pop')) {
  names(out_data[[n]]) <- var_names
  out_data[[n]]$year <- years
  col_idx <- grep("year", names(out_data[[n]]))
  out_data[[n]] <- out_data[[n]][, c(col_idx, (1:ncol(out_data[[n]]))[-col_idx])]
  write.csv(out_data[[n]], paste0(out.dir, "/results_", n , i, ".csv"), row.names=F)
}

## Plot results
plot.fit(result, pdf.path, nat.data)
