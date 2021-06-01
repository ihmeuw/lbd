################################################################################
## Purpose: Compile draws of incidence, prevalence, mortality, art coverage, and
##          population from EPP outputs so that we have the aproximated posterior
##          distribution to use as estimates.
## Run instructions: Launched from the launch script, once per high level
##                   geography.  This script is also the main clean up script
##                   used to minimze the size of data on disk.  So, once you run
##                   this you loose a bunch of intermediary files.
################################################################################

### Setup
rm(list = ls())
gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
windows <- Sys.info()[1] == "Windows"
code.dir <- paste0("<<<< FILEPATH REDACTED >>>>", "/lbd_hiv/epp/2_modeling/epp_geo/")

## Packages
library(data.table); library(parallel)

## Arguments

run_name <- as.character(commandArgs()[5])
message(paste0("the run name is ",run_name))
loc <- as.character(commandArgs()[4])
message(paste0("the loc is ",loc))

### read in the config file

config <- read.csv(paste0("<<<< FILEPATH REDACTED >>>>"))
mean_test <- as.logical(config$Value[which(config$Setting == "mean_test")])


if (mean_test == FALSE) {

n.runs <- as.numeric(as.character(config$Value[which(config$Setting == "n.draws")]))

debug <- as.logical(config$Value[which(config$Setting == "debug")])
message(paste0("the debug setting is ", debug))

### Paths
in.dir <- paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/")
prev_dir <- paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/")
inc_dir <- paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/")
dir.create(prev_dir, showWarnings = F)
dir.create(inc_dir, showWarnings = F)

prev_path <- paste0(prev_dir, loc, "_SPU_prev_draws.csv")
inc_path <- paste0(inc_dir, loc, "_SPU_inc_draws.csv")
art_path <- paste0(prev_dir, loc, "_SPU_art_draws.csv")
pop_path <- paste0(inc_dir, loc, "_SPU_pop_draws.csv")
r_path <- paste0(inc_dir, loc, "_SPU_r_draws.csv")
mort_path <- paste0(inc_dir, loc, "_SPU_mort_draws.csv")
pops_path <- paste0(inc_dir, loc, "_SPU_pop_s_draws.csv")

### Code
## Incidence
dt1 <- NULL

missing.runs <- c()
kept.draws <- c()
replacement.draws <- c()
j = 1

nums_dt <- c(1:1000)
names_dt <- c("year", nums_dt)




for (i in c(1:n.runs)) {
  file <- paste0(loc,'_results_incid',i,'.csv')
  if (file.exists(paste0(in.dir, file)) & (file.info(paste0(in.dir, file))$size > 0)) {
    print(paste("Found run",i))
    dt <- fread(paste0(in.dir, file))
    tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
    keep.draw <- sample(1:tmp.n.draws, ifelse(n.runs == 1, 1000, (1000/n.runs)))
    k <- 1
    while (any(is.na(dt[,c('year', paste0('draw',keep.draw)), with = F])) & k != 50) {
      keep.draw <- sample(1:tmp.n.draws, ifelse(n.runs == 1000, 1, (1000/n.runs)))
      k <- k + 1
    }
    if (k == 50) {
      print(paste("Missing draw",i))
      missing.runs <- c(missing.runs, i)
      kept.draws <- c(kept.draws, rep.int(0, (1000/n.runs)))
      next
    }
    kept.draws <- c(kept.draws, keep.draw)
    dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "incid")
    dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
    dt1 <- rbindlist(list(dt1,dt))
  }  else {
    print(paste("Missing run",i))
    missing.runs <- c(missing.runs, i)
    kept.draws <- c(kept.draws, rep.int(0, (1000/n.runs)))
  }
}

if (length(missing.runs) > 0) {
  replace.with <- sample(dt1[,unique(draw)], length(missing.runs), replace = TRUE)
  for (i in 1:length(missing.runs)) {
    print(paste("replacement run",i))
    file <- paste0(loc,'_results_incid',replace.with[[i]],'.csv')
    dt <- fread(paste0(in.dir, file))
    tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
    keep.draw <- sample(1:tmp.n.draws, ifelse(n.runs == 1, 1000, (1000/n.runs)))
    replacement.draws <- c(replacement.draws, keep.draw)
    dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "incid")
    dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
    dt$run <- paste0(dt$run, "_r")
    dt1 <- rbindlist(list(dt1,dt))
  }
}

dt1[,draw2 := paste0(draw,"_", run)]
dt1[,draw := draw2]
dt1[,draw2 := NULL]
dt1[,run := NULL]
dt1[,incid := incid*100]
setnames(dt1, c('draw', 'incid'), c('run', 'draw'))
reshaped.dt <- dcast.data.table(dt1,year~run, value.var = "draw")
col_order <- names(reshaped.dt)
names(reshaped.dt) <- names_dt
reshaped.dt <- reshaped.dt[order(year),]

write.csv(reshaped.dt, file = inc_path, row.names = F)


## Prevalence
prev_dt <- NULL

j = 1
for (i in setdiff(c(1:n.runs), missing.runs)) {
  file <- paste0(loc,'_results_prev',i,'.csv')
  if (file.exists(paste0(in.dir, file))) {
    print(paste('prev',i))
    dt <- fread(paste0(in.dir, file))
    dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "prev")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    prev_dt <- rbindlist(list(prev_dt,dt))
    j <- j + 1
  }
}

if (length(missing.runs) > 0) {
  for (i in 1:length(missing.runs)) {
    print(paste("replacement run",i))
    file <- paste0(loc,'_results_prev',replace.with[[i]],'.csv')
    dt <- fread(paste0(in.dir, file))
    tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
    keep.draw <- replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))]
    dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "prev")
    dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
    dt$run <- paste0(dt$run, "_r")
    prev_dt <- rbindlist(list(prev_dt,dt))
  }
}

prev_dt[,draw2 := paste0(draw,"_", run)]
prev_dt[,draw := draw2]
prev_dt[,draw2 := NULL]
prev_dt[,run := NULL]
prev_dt[,prev := prev*100]
setnames(prev_dt, c('draw', 'prev'), c('run', 'draw'))
out.prev <- dcast.data.table(prev_dt, year~run, value.var = "draw")
setcolorder(out.prev, col_order)
names(out.prev) <- names_dt
out.prev <- out.prev[order(year),]

write.csv(out.prev, file = prev_path, row.names = F)


## ART coverage
art_dt <- NULL

j = 1
for (i in setdiff(c(1:n.runs), missing.runs)) {
  file <- paste0(loc,'_results_art',i,'.csv')
  if (file.exists(paste0(in.dir, file))) {
    print(paste('art',i))
    dt <- fread(paste0(in.dir, file))
    dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "art")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    art_dt <- rbindlist(list(art_dt,dt))
    j <- j + 1
  }
}

if (length(missing.runs) > 0) {
  for (i in 1:length(missing.runs)) {
    print(paste("replacement run",i))
    file <- paste0(loc,'_results_art',replace.with[[i]],'.csv')
    dt <- fread(paste0(in.dir, file))
    tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
    keep.draw <- replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))]
    dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "art")
    dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
    dt$run <- paste0(dt$run, "_r")
    art_dt <- rbindlist(list(art_dt,dt))
  }
}

art_dt[,draw2 := paste0(draw,"_", run)]
art_dt[,draw := draw2]
art_dt[,draw2 := NULL]
art_dt[,run := NULL]
art_dt[,art := art*100]
setnames(art_dt, c('draw', 'art'), c('run', 'draw'))
out.art <- dcast.data.table(art_dt, year~run, value.var = "draw")
setcolorder(out.art, col_order)
names(out.art) <- names_dt
out.art <- out.art[order(year),]

write.csv(out.art, file = art_path, row.names = F)


## Pop
pop_dt <- NULL

j = 1
for (i in setdiff(c(1:n.runs), missing.runs)) {
  file <- paste0(loc,'_results_pop',i,'.csv')
  if (file.exists(paste0(in.dir, file))) {
    print(paste('pop',i))
    dt <- fread(paste0(in.dir, file))
    dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    pop_dt <- rbindlist(list(pop_dt,dt))
    j <- j + 1
  }
}

if (length(missing.runs) > 0) {
  for (i in 1:length(missing.runs)) {
    print(paste("replacement run",i))
    file <- paste0(loc,'_results_pop',replace.with[[i]],'.csv')
    dt <- fread(paste0(in.dir, file))
    tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
    keep.draw <- replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))]
    dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop")
    dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
    dt$run <- paste0(dt$run, "_r")
    pop_dt <- rbindlist(list(pop_dt,dt))
  }
}

pop_dt[,draw2 := paste0(draw,"_", run)]
pop_dt[,draw := draw2]
pop_dt[,draw2 := NULL]
pop_dt[,run := NULL]
setnames(pop_dt, c('draw', 'pop'), c('run', 'draw'))
out.pop <- dcast.data.table(pop_dt, year~run, value.var = "draw")
setcolorder(out.pop, col_order)
names(out.pop) <- names_dt
out.pop <- out.pop[order(year),]

write.csv(out.pop, file = pop_path, row.names = F)


## force of infection
rvec_dt <- NULL

j = 1
for (i in setdiff(c(1:n.runs), missing.runs)) {
  file <- paste0(loc,'_',i,'_results.RData')
  if (file.exists(paste0(in.dir, file))) {
    print(paste('rvec',i))
    load(paste0(in.dir, file))
    dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$rvec))
    overs <- paste0("draw", 1:1000)
    names(dt) <- c("year", overs)
    yrs <- as.list(out.pop$year)
    dt <- dt[which(dt$year %in% yrs), ]
    dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "r")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    rvec_dt <- rbindlist(list(rvec_dt,dt))
    j <- j + 1
  }
}

if (length(missing.runs) > 0) {
  for (i in 1:length(missing.runs)) {
    print(paste("replacement run",i))
    file <- paste0(loc,'_',replace.with[[i]],'_results.RData')
    load(paste0(in.dir, file))
    dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$rvec))
    overs <- paste0("draw", 1:1000)
    names(dt) <- c("year", overs)
    yrs <- as.list(out.pop$year)
    dt <- dt[which(dt$year %in% yrs), ]
    dt <- dt[,c('year', paste0('draw',replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "r")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    dt$run <- paste0(dt$run, "_r")
    rvec_dt <- rbindlist(list(rvec_dt,dt))
  }
}

rvec_dt[,draw2 := paste0(draw,"_", run)]
rvec_dt[,draw := draw2]
rvec_dt[,draw2 := NULL]
rvec_dt[,run := NULL]
setnames(rvec_dt, c('draw', 'r'), c('run', 'draw'))
out.r <- dcast.data.table(rvec_dt, year~run, value.var = "draw")
setcolorder(out.r, col_order)
names(out.r) <- names_dt
out.r <- out.r[order(year),]

write.csv(out.r, file = r_path, row.names = F)


## mortality
mvec_dt <- NULL

j = 1
for (i in setdiff(c(1:n.runs), missing.runs)) {
  file <- paste0(loc,'_',i,'_results.RData')
  if (file.exists(paste0(in.dir, file))) {
    print(paste('mvec',i))
    load(paste0(in.dir, file))
    dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$mvec))
    overs <- paste0("draw", 1:1000)
    names(dt) <- c("year", overs)
    dt$year <- as.numeric(substr(as.character(dt$year),1,4))
    collapsed_mort <- dt[,lapply(overs, function(x) sum((get(x)*0.1))), by = c("year")] #note the 0.1 is a relic of the location in the EPP code where we extracted these values.  it is intentinal and appropriate.
    names(collapsed_mort) <- c("year", overs)
    dt <- collapsed_mort
    dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "mort")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    mvec_dt <- rbindlist(list(mvec_dt,dt))
    j <- j + 1
  }
}

if (length(missing.runs) > 0) {
  for (i in 1:length(missing.runs)) {
    print(paste("replacement run",i))
    file <- paste0(loc,'_',replace.with[[i]],'_results.RData')
    load(paste0(in.dir, file))
    dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$mvec))
    overs <- paste0("draw", 1:1000)
    names(dt) <- c("year", overs)
    dt$year <- as.numeric(substr(as.character(dt$year),1,4))
    collapsed_mort <- dt[,lapply(overs, function(x) sum((get(x)*0.1))), by = c("year")] #note the 0.1 is a relic of the location in the EPP code where we extracted these values.  it is intentinal and appropriate.
    names(collapsed_mort) <- c("year", overs)
    dt <- collapsed_mort
    dt <- dt[,c('year', paste0('draw',replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "mort")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    dt$run <- paste0(dt$run, "_r")
    mvec_dt <- rbindlist(list(mvec_dt,dt))
  }
}

mvec_dt[,draw2 := paste0(draw,"_", run)]
mvec_dt[,draw := draw2]
mvec_dt[,draw2 := NULL]
mvec_dt[,run := NULL]
setnames(mvec_dt, c('draw', 'mort'), c('run', 'draw'))
out.mort <- dcast.data.table(mvec_dt, year~run, value.var = "draw")
setcolorder(out.mort, col_order)
names(out.mort) <- names_dt
out.mort <- out.mort[order(year),]

write.csv(out.mort, file = mort_path, row.names = F)



## pop Scalar
svec_dt <- NULL

j = 1
for (i in setdiff(c(1:n.runs), missing.runs)) {
  file <- paste0(loc,'_',i,'_results.RData')
  if (file.exists(paste0(in.dir, file))) {
    print(paste('svec',i))
    load(paste0(in.dir, file))
    dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$svec))
    overs <- paste0("draw", 1:1000)
    names(dt) <- c("year", overs)
    yrs <- as.list(out.pop$year)
    dt <- dt[which(dt$year %in% yrs), ]
    dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop_scalar")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    svec_dt <- rbindlist(list(svec_dt,dt))
    j <- j + 1
  }
}

if (length(missing.runs) > 0) {
  for (i in 1:length(missing.runs)) {
    print(paste("replacement run",i))
    file <- paste0(loc,'_',replace.with[[i]],'_results.RData')
    load(paste0(in.dir, file))
    dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$svec))
    overs <- paste0("draw", 1:1000)
    names(dt) <- c("year", overs)
    yrs <- as.list(out.pop$year)
    dt <- dt[which(dt$year %in% yrs), ]
    dt <- dt[,c('year', paste0('draw',replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
    dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop_scalar")
    dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
    dt$run <- paste0(dt$run, "_r")
    svec_dt <- rbindlist(list(svec_dt,dt))
  }
}

svec_dt[,draw2 := paste0(draw,"_", run)]
svec_dt[,draw := draw2]
svec_dt[,draw2 := NULL]
svec_dt[,run := NULL]
setnames(svec_dt, c('draw', "pop_scalar"), c('run', 'draw'))
out.s <- dcast.data.table(svec_dt, year~run, value.var = "draw")
setcolorder(out.s, col_order)
names(out.s) <- names_dt
out.s <- out.s[order(year),]

write.csv(out.s, file = pops_path, row.names = F)





run.admin2 <- as.logical(config$Value[which(config$Setting == "run.admin2")])
message(paste0("the run.admin2 setting is ", run.admin2))

stash_dir <- paste0("<<<< FILEPATH REDACTED >>>>",run_name,"/missing_runs/")
dir.create(stash_dir)

saveRDS(missing.runs , file = paste0(stash_dir, loc,"_missing_runs.rds"))


# clean up files

#if (debug == TRUE) {
#  files <- list.files(paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/"), full.names = TRUE)
#  files_to_delete <- grep(paste0(loc,"_results"), files, value = TRUE)
#  for (fi in files_to_delete) {
#    file.remove(fi)
#  }
#}

#if (debug == FALSE) {
#  files <- list.files(paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/"), full.names = TRUE)
#  files_to_delete <- grep(paste0("_results"), files, value = TRUE)
#  for (fi in files_to_delete) {
#    file.remove(fi)
#  }
#}



if (run.admin2 == TRUE) {
  ad2_dir <- paste0(in.dir, "admin2_models/")
  admin_dirs <- list.dirs(ad2_dir, full.names = FALSE)
  for (f in 2:length(admin_dirs)) {

    subnat <- admin_dirs[[f]]
    subnat_dir <- paste0(ad2_dir, subnat, "/")
    prev_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_prev_draws.csv")
    inc_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_inc_draws.csv")
    art_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_art_draws.csv")
    pop_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_pop_draws.csv")
    r_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_r_draws.csv")
    mort_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_mort_draws.csv")
    pops_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_pop_s_draws.csv")

    ### Code
    ## Incidence
    dt1 <- NULL

    missing.runs <- c()
    kept.draws <- c()
    replacement.draws <- c()
    j = 1


    for (i in c(1:n.runs)) {
      file <- paste0(subnat,'_results_incid',i,'.csv')
      if ((file.exists(paste0(subnat_dir, file))) & (file.info(paste0(subnat_dir, file))$size > 0)) {
        #print(paste("Found run",i))
        dt <- fread(paste0(subnat_dir, file))
        tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
        keep.draw <- sample(1:tmp.n.draws, ifelse(n.runs == 1, 1000, (1000/n.runs)))
        k <- 1
        while (any(is.na(dt[,c('year', paste0('draw',keep.draw)), with = F])) & k != 50) {
          keep.draw <- sample(1:tmp.n.draws, ifelse(n.runs == 1000, 1, 1000))
          k <- k + 1
        }
        if (k == 50) {
          print(paste("Missing draw",i))
          missing.runs <- c(missing.runs, i)
          kept.draws <- c(kept.draws, rep.int(0, (1000/n.runs)))

          next
        }
        kept.draws <- c(kept.draws, keep.draw)
        dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "incid")
        dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
        dt1 <- rbindlist(list(dt1,dt))
      }  else {
        print(paste("Missing run",i))
        missing.runs <- c(missing.runs, i)
        kept.draws <- c(kept.draws, rep.int(0, (1000/n.runs)))
      }
    }
    if (length(missing.runs) == n.runs) { next } #skip to the next admin 2 if there is no info for this admin 2
    saveRDS(missing.runs , file = paste0(stash_dir, loc,"_",subnat,"_missing_runs.rds"))
    if (length(missing.runs) > 0) {
      replace.with <- sample(dt1[,unique(draw)], length(missing.runs), replace = TRUE)
      for (i in 1:length(missing.runs)) {
        print(paste("replacement run",i))
        file <- paste0(subnat,'_results_incid',replace.with[[i]],'.csv')
        dt <- fread(paste0(subnat_dir, file))
        tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
        keep.draw <- sample(1:tmp.n.draws, ifelse(n.runs == 1, 1000, (1000/n.runs)))
        replacement.draws <- c(replacement.draws, keep.draw)
        dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "incid")
        dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
        dt$run <- paste0(dt$run, "_r")
        dt1 <- rbindlist(list(dt1,dt))
      }
    }

    dt1[,draw2 := paste0(draw,"_", run)]
    dt1[,draw := draw2]
    dt1[,draw2 := NULL]
    dt1[,run := NULL]
    dt1[,incid := incid*100]
    setnames(dt1, c('draw', 'incid'), c('run', 'draw'))
    reshaped.dt <- dcast.data.table(dt1,year~run, value.var = "draw")
    col_order <- names(reshaped.dt)
    names(reshaped.dt) <- names_dt
    reshaped.dt <- reshaped.dt[order(year),]

    write.csv(reshaped.dt, file = inc_path, row.names = F)


    ## Prevalence
    prev_dt <- NULL

    j = 1
    for (i in setdiff(c(1:n.runs), missing.runs)) {
      file <- paste0(subnat,'_results_prev',i,'.csv')
      if (file.exists(paste0(subnat_dir, file)) & (file.info(paste0(subnat_dir, file))$size > 0)) {
        #print(paste('prev',i))
        dt <- fread(paste0(subnat_dir, file))
        dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "prev")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        prev_dt <- rbindlist(list(prev_dt,dt))
        j <- j + 1
      }
    }

    if (length(missing.runs) > 0) {
      for (i in 1:length(missing.runs)) {
        print(paste("replacement run",i))
        file <- paste0(subnat,'_results_prev',replace.with[[i]],'.csv')
        dt <- fread(paste0(subnat_dir, file))
        tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
        keep.draw <- replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))]
        dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "prev")
        dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
        dt$run <- paste0(dt$run, "_r")
        prev_dt <- rbindlist(list(prev_dt,dt))
      }
    }
    prev_dt <- as.data.table(prev_dt)
    prev_dt[,draw2 := paste0(draw,"_", run)]
    prev_dt[,draw := draw2]
    prev_dt[,draw2 := NULL]
    prev_dt[,run := NULL]
    prev_dt[,prev := prev*100]
    setnames(prev_dt, c('draw', 'prev'), c('run', 'draw'))
    out.prev <- dcast.data.table(prev_dt, year~run, value.var = "draw")
    setcolorder(out.prev, col_order)
    names(out.prev) <- names_dt
    out.prev <- out.prev[order(year),]


    write.csv(out.prev, file = prev_path, row.names = F)

    ## art
    art_dt <- NULL

    j = 1
    for (i in setdiff(c(1:n.runs), missing.runs)) {
      file <- paste0(subnat,'_results_art',i,'.csv')
      if (file.exists(paste0(subnat_dir, file)) & (file.info(paste0(subnat_dir, file))$size > 0)) {
        #print(paste('art ',i))
        dt <- fread(paste0(subnat_dir, file))
        dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "art")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        art_dt <- rbindlist(list(art_dt,dt))
        j <- j + 1
      }
    }

    if (length(missing.runs) > 0) {
      for (i in 1:length(missing.runs)) {
        print(paste("replacement run",i))
        file <- paste0(subnat,'_results_art',replace.with[[i]],'.csv')
        dt <- fread(paste0(subnat_dir, file))
        tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
        keep.draw <- replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))]
        dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "art")
        dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
        dt$run <- paste0(dt$run, "_r")
        art_dt <- rbindlist(list(art_dt,dt))
      }
    }
    art_dt <- as.data.table(art_dt)
    art_dt[,draw2 := paste0(draw,"_", run)]
    art_dt[,draw := draw2]
    art_dt[,draw2 := NULL]
    art_dt[,run := NULL]
    art_dt[,art := art*100]
    setnames(art_dt, c('draw', 'art'), c('run', 'draw'))
    out.art <- dcast.data.table(art_dt, year~run, value.var = "draw")
    setcolorder(out.art, col_order)
    names(out.art) <- names_dt
    out.art <- out.art[order(year),]

    write.csv(out.art, file = art_path, row.names = F)


    ## Pop
    pop_dt <- NULL

    j = 1
    for (i in setdiff(c(1:n.runs), missing.runs)) {
      file <- paste0(subnat,'_results_pop',i,'.csv')
      if (file.exists(paste0(subnat_dir, file)) & (file.info(paste0(subnat_dir, file))$size > 0)) {
        #print(paste('pop ',i))
        dt <- fread(paste0(subnat_dir, file))
        dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        pop_dt <- rbindlist(list(pop_dt,dt))
        j <- j + 1
      }
    }

    if (length(missing.runs) > 0) {
      for (i in 1:length(missing.runs)) {
        print(paste("replacement run",i))
        file <- paste0(subnat,'_results_pop',replace.with[[i]],'.csv')
        dt <- fread(paste0(subnat_dir, file))
        tmp.n.draws <- length(names(dt)[names(dt) != 'year'])
        keep.draw <- replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))]
        dt <- dt[,c('year', paste0('draw',keep.draw)), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop")
        dt <- dt[,draw := ifelse(n.runs != 1, i, as.integer(sub("draw", "", run))), by = "run"]
        dt$run <- paste0(dt$run, "_r")
        pop_dt <- rbindlist(list(pop_dt,dt))
      }
    }
    pop_dt <- as.data.table(pop_dt)
    pop_dt[,draw2 := paste0(draw,"_", run)]
    pop_dt[,draw := draw2]
    pop_dt[,draw2 := NULL]
    pop_dt[,run := NULL]
    setnames(pop_dt, c('draw', 'pop'), c('run', 'draw'))
    out.pop <- dcast.data.table(pop_dt, year~run, value.var = "draw")
    setcolorder(out.pop, col_order)
    names(out.pop) <- names_dt
    out.pop <- out.pop[order(year),]

    write.csv(out.pop, file = pop_path, row.names = F)

    ## force of infection
    rvec_dt <- NULL

    j = 1
    for (i in setdiff(c(1:n.runs), missing.runs)) {
      file <- paste0(subnat,'_',i,'_results.RData')
      if (file.exists(paste0(subnat_dir, file))) {
        #print(paste('rvec ',i))
        load(paste0(subnat_dir, file))
        dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$rvec))
        overs <- paste0("draw", 1:1000)
        names(dt) <- c("year", overs)
        yrs <- as.list(out.pop$year)
        dt <- dt[which(dt$year %in% yrs), ]
        dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "r")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        rvec_dt <- rbindlist(list(rvec_dt,dt))
        j <- j + 1
      }
    }

    if (length(missing.runs) > 0) {
      for (i in 1:length(missing.runs)) {
        print(paste("replacement run",i))
        file <- paste0(subnat,'_',replace.with[[i]],'_results.RData')
        load(paste0(subnat_dir, file))
        dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$rvec))
        overs <- paste0("draw", 1:1000)
        names(dt) <- c("year", overs)
        yrs <- as.list(out.pop$year)
        dt <- dt[which(dt$year %in% yrs), ]
        dt <- dt[,c('year', paste0('draw',replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "r")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        dt$run <- paste0(dt$run, "_r")
        rvec_dt <- rbindlist(list(rvec_dt,dt))
      }
    }

    rvec_dt[,draw2 := paste0(draw,"_", run)]
    rvec_dt[,draw := draw2]
    rvec_dt[,draw2 := NULL]
    rvec_dt[,run := NULL]
    setnames(rvec_dt, c('draw', 'r'), c('run', 'draw'))
    out.r <- dcast.data.table(rvec_dt, year~run, value.var = "draw")
    setcolorder(out.r, col_order)
    names(out.r) <- names_dt
    out.r <- out.r[order(year),]

    write.csv(out.r, file = r_path, row.names = F)

    ## mortality
    mvec_dt <- NULL

    j = 1
    for (i in setdiff(c(1:n.runs), missing.runs)) {
      file <- paste0(subnat,'_',i,'_results.RData')
      if (file.exists(paste0(subnat_dir, file))) {
        #print(paste('mvec',i))
        load(paste0(subnat_dir, file))
        dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$mvec))
        overs <- paste0("draw", 1:1000)
        names(dt) <- c("year", overs)
        dt$year <- as.numeric(substr(as.character(dt$year),1,4))
        collapsed_mort <- dt[,lapply(overs, function(x) sum((get(x)*0.1))), by = c("year")] #note the 0.1 is a relic of the location in the EPP code where we extracted these values.  it is intentinal and appropriate.
        names(collapsed_mort) <- c("year", overs)
        dt <- collapsed_mort
        dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "mort")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        mvec_dt <- rbindlist(list(mvec_dt,dt))
        j <- j + 1
      }
    }

    if (length(missing.runs) > 0) {
      for (i in 1:length(missing.runs)) {
        print(paste("replacement run",i))
        file <- paste0(subnat,'_',replace.with[[i]],'_results.RData')
        load(paste0(subnat_dir, file))
        dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$mvec))
        overs <- paste0("draw", 1:1000)
        names(dt) <- c("year", overs)
        dt$year <- as.numeric(substr(as.character(dt$year),1,4))
        collapsed_mort <- dt[,lapply(overs, function(x) sum((get(x)*0.1))), by = c("year")] #note the 0.1 is a relic of the location in the EPP code where we extracted these values.  it is intentinal and appropriate.
        names(collapsed_mort) <- c("year", overs)
        dt <- collapsed_mort
        dt <- dt[,c('year', paste0('draw',replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "mort")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        dt$run <- paste0(dt$run, "_r")
        mvec_dt <- rbindlist(list(mvec_dt,dt))
      }
    }

    mvec_dt[,draw2 := paste0(draw,"_", run)]
    mvec_dt[,draw := draw2]
    mvec_dt[,draw2 := NULL]
    mvec_dt[,run := NULL]
    setnames(mvec_dt, c('draw', 'mort'), c('run', 'draw'))
    out.mort <- dcast.data.table(mvec_dt, year~run, value.var = "draw")
    setcolorder(out.mort, col_order)
    names(out.mort) <- names_dt
    out.mort <- out.mort[order(year),]

    write.csv(out.mort, file = mort_path, row.names = F)

    ## pop scalar
    svec_dt <- NULL

    j = 1
    for (i in setdiff(c(1:n.runs), missing.runs)) {
      file <- paste0(subnat,'_',i,'_results.RData')
      if (file.exists(paste0(subnat_dir, file))) {
        #print(paste('svec ',i))
        load(paste0(subnat_dir, file))
        dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$svec))
        overs <- paste0("draw", 1:1000)
        names(dt) <- c("year", overs)
        yrs <- as.list(out.pop$year)
        dt <- dt[which(dt$year %in% yrs), ]
        dt <- dt[,c('year', paste0('draw',kept.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop_scalar")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        svec_dt <- rbindlist(list(svec_dt,dt))
        j <- j + 1
      }
    }

    if (length(missing.runs) > 0) {
      for (i in 1:length(missing.runs)) {
        print(paste("replacement run",i))
        file <- paste0(subnat,'_',replace.with[[i]],'_results.RData')
        load(paste0(subnat_dir, file))
        dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$svec))
        overs <- paste0("draw", 1:1000)
        names(dt) <- c("year", overs)
        yrs <- as.list(out.pop$year)
        dt <- dt[which(dt$year %in% yrs), ]
        dt <- dt[,c('year', paste0('draw',replacement.draws[(((i - 1)*(1000/n.runs)) + 1):((i*(1000/n.runs)))])), with = F]
        dt <- data.table::melt(dt, id.vars = "year", variable.name = "run", value.name = "pop_scalar")
        dt <- dt[,draw := ifelse(n.runs == 1, as.integer(sub("draw", "", run)), i), by = "run"]
        dt$run <- paste0(dt$run, "_r")
        svec_dt <- rbindlist(list(svec_dt,dt))
      }
    }

    svec_dt[,draw2 := paste0(draw,"_", run)]
    svec_dt[,draw := draw2]
    svec_dt[,draw2 := NULL]
    svec_dt[,run := NULL]
    setnames(svec_dt, c('draw', "pop_scalar"), c('run', 'draw'))
    out.s <- dcast.data.table(svec_dt, year~run, value.var = "draw")
    setcolorder(out.s, col_order)
    names(out.s) <- names_dt
    out.s <- out.s[order(year),]

    write.csv(out.s, file = pops_path, row.names = F)



    print(f)
    print(subnat)

    # clean up files

    #if (debug == TRUE) {
    #  files <- list.files(subnat_dir, full.names = TRUE)
    #  files_to_delete <- grep(paste0(subnat,"_results"), files, value = TRUE)
    #  for (fi in files_to_delete) {
    #    file.remove(fi)
    #  }
    #}

    #if (debug == FALSE) {
    #  files <- list.files(subnat_dir, full.names = TRUE)
    #  files_to_delete <- grep(paste0("_results"), files, value = TRUE)
    #  for (fi in files_to_delete) {
    #    file.remove(fi)
    #  }
    #}
  }
}

} else if (mean_test == TRUE) {

  debug <- as.logical(config$Value[which(config$Setting == "debug")])
  message(paste0("the debug setting is ", debug))

  ### Paths
  in.dir <- paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/")
  prev_dir <- paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/")
  inc_dir <- paste0('<<<< FILEPATH REDACTED >>>>', run_name, "/", loc, "/")
  dir.create(prev_dir, showWarnings = F)
  dir.create(inc_dir, showWarnings = F)

  prev_path <- paste0(prev_dir, loc, "_SPU_prev_draws.csv")
  inc_path <- paste0(inc_dir, loc, "_SPU_inc_draws.csv")
  art_path <- paste0(prev_dir, loc, "_SPU_art_draws.csv")
  pop_path <- paste0(inc_dir, loc, "_SPU_pop_draws.csv")
  r_path <- paste0(inc_dir, loc, "_SPU_r_draws.csv")
  mort_path <- paste0(inc_dir, loc, "_SPU_mort_draws.csv")


  nums_dt <- c(1:1000)
  names_dt <- c("year", nums_dt)



  inc <- fread(paste0(in.dir, loc,"_results_incid1.csv"))
  overs <- paste0("draw", 1:1000)
  inc <- inc[,(overs) := lapply(overs, function(x) get(x) * 100) ]
  names(inc) <- names_dt
  write.csv(inc, inc_path, row.names = F)

  prev <- fread(paste0(in.dir, loc,"_results_prev1.csv"))
  overs <- paste0("draw", 1:1000)
  prev <- prev[,(overs) := lapply(overs, function(x) get(x) * 100) ]
  names(prev) <- names_dt
  write.csv(prev, prev_path, row.names = F)

  pop <- fread(paste0(in.dir, loc,"_results_pop1.csv"))
  names(pop) <- names_dt
  write.csv(pop, pop_path, row.names = F)

  art <- fread(paste0(in.dir, loc,"_results_art1.csv"))
  overs <- paste0("draw", 1:1000)
  art <- art[,(overs) := lapply(overs, function(x) get(x) * 100) ]
  names(art) <- names_dt
  write.csv(art, art_path, row.names = F)

  file <- paste0(loc,'_1_results.RData')
  load(paste0(in.dir, file))
  dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$mvec))
  overs <- paste0("draw", 1:1000)
  names(dt) <- c("year", overs)
  dt$year <- as.numeric(substr(as.character(dt$year),1,4))
  collapsed_mort <- dt[,lapply(overs, function(x) sum((get(x)*0.1))), by = c("year")] #note the 0.1 is a relic of the location in the EPP code where we extracted these values.  it is intentinal and appropriate.
  names(collapsed_mort) <- c("year", nums_dt)
  dt <- collapsed_mort
  write.csv(dt, mort_path, row.names = F)


  file <- paste0(loc,'_1_results.RData')
  load(paste0(in.dir, file))
  dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$rvec))
  overs <- paste0("draw", 1:1000)
  names(dt) <- c("year", nums_dt)
  yrs <- as.list(inc$year)
  dt <- dt[which(dt$year %in% yrs), ]
  write.csv(dt, r_path, row.names = F)



  run.admin2 <- as.logical(config$Value[which(config$Setting == "run.admin2")])
  message(paste0("the run.admin2 setting is ", run.admin2))

  stash_dir <- paste0("<<<< FILEPATH REDACTED >>>>")
  dir.create(stash_dir)

  missing.runs <- list()

  saveRDS(missing.runs , file = paste0(stash_dir, loc,"_missing_runs.rds"))

  if (run.admin2 == TRUE) {
    ad2_dir <- paste0(in.dir, "admin2_models/")
    admin_dirs <- list.dirs(ad2_dir, full.names = FALSE)
    for (f in 2:length(admin_dirs)) {
      subnat <- admin_dirs[[f]]
      subnat_dir <- paste0(ad2_dir, subnat, "/")
      prev_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_prev_draws.csv")
      inc_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_inc_draws.csv")
      art_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_art_draws.csv")
      pop_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_pop_draws.csv")
      r_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_r_draws.csv")
      mort_path <- paste0(subnat_dir, loc,"_",subnat,"_SPU_mort_draws.csv")

      inc <- fread(paste0(subnat_dir, subnat,"_results_incid1.csv"))
      overs <- paste0("draw", 1:1000)
      inc <- inc[,(overs) := lapply(overs, function(x) get(x) * 100) ]
      names(inc) <- names_dt
      write.csv(inc, inc_path, row.names = F)

      prev <- fread(paste0(subnat_dir, subnat,"_results_prev1.csv"))
      overs <- paste0("draw", 1:1000)
      prev <- prev[,(overs) := lapply(overs, function(x) get(x) * 100) ]
      names(prev) <- names_dt
      write.csv(prev, prev_path, row.names = F)

      pop <- fread(paste0(subnat_dir, subnat,"_results_pop1.csv"))
      names(pop) <- names_dt
      write.csv(pop, pop_path, row.names = F)

      art <- fread(paste0(subnat_dir, subnat,"_results_art1.csv"))
      overs <- paste0("draw", 1:1000)
      art <- art[,(overs) := lapply(overs, function(x) get(x) * 100) ]
      names(art) <- names_dt
      write.csv(art, art_path, row.names = F)

      file <- paste0(subnat,'_1_results.RData')
      load(paste0(subnat_dir, file))
      dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$mvec))
      overs <- paste0("draw", 1:1000)
      names(dt) <- c("year", overs)
      dt$year <- as.numeric(substr(as.character(dt$year),1,4))
      collapsed_mort <- dt[,lapply(overs, function(x) sum((get(x)*0.1))), by = c("year")] #note the 0.1 is a relic of the location in the EPP code where we extracted these values.  it is intentinal and appropriate.
      names(collapsed_mort) <- c("year", nums_dt)
      dt <- collapsed_mort
      write.csv(dt, mort_path, row.names = F)


      file <- paste0(subnat,'_1_results.RData')
      load(paste0(subnat_dir, file))
      dt <- as.data.table(cbind(result[[1]]$fp$proj.steps, result[[1]]$rvec))
      overs <- paste0("draw", 1:1000)
      names(dt) <- c("year", nums_dt)
      yrs <- as.list(inc$year)
      dt <- dt[which(dt$year %in% yrs), ]
      write.csv(dt, r_path, row.names = F)

      saveRDS(missing.runs , file = paste0(stash_dir, loc,"_",subnat,"_missing_runs.rds"))
    }

  }

}




### End
