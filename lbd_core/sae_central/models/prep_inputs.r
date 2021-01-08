####################################################################################################
## Description: Prep inputs for small area models for area-year-age-sex all-cause or cause-specific
##              mortality rates. Specifically, collapse and merge events and population data and
##              then combine with covariates data. Recode age groups and years to run sequentially
##              from zero. Also, prep CAR structure matrices for spatial, temporal, and age
##              effects.
##
## Passed args: main_dir [character] -- home directory for settings and final output
##
## Requires:    events data (events_file)
##              population data (pop_file)
##              covariates data, if used (covars, covars_as, covars_trans, covar_file, covar_as_file)
##              neighborhood adjacency matrix (adjmat_file)
##
## Outputs:     prepped data file ('data.rdata' in temp_dir)
##
## Run from within 'sae_central' directory!!
####################################################################################################

library(Matrix)
library(data.table)

rm(list=ls())

## Get settings ------------------------------------------------------------------------------------
main_dir <- commandArgs()[4]

source("settings.r")
get_settings(main_dir)

## Load and prep inputs ----------------------------------------------------------------------------
# combine events and population
message("Combining events and population")
load(events_file)
setnames(events, area_var, "area")
events <- events[year %in% years & age %in% ages, list(events = sum(events)), by='area,year,sex,age'] # removes any stratifying variables (e.g., race)

load(pop_file)
setnames(pop, area_var, "area")
pop <- pop[year %in% years & age %in% ages, list(pop = sum(pop)), by='area,year,sex,age']

data <- merge(pop, events, by=c("area", "year", "sex", "age"), all=T)
data[is.na(events), events := 0]
rm(pop, events); gc()

# adjust population data, if necessary, so that events <= pop
data[events > pop, pop := events]

# load covariates
message("Loading covariates")
if (!is.null(covars) | !is.null(covars_as)) {
  if (is.null(covars_as)) {
    load(covar_file)
    covar <- covar[year %in% years, c(area_var, "year", covars), with=F]
  } else if (is.null(covars)) {
    load(covar_as_file)
    covar <- covar[year %in% years, c(area_var, "year", "sex", "age", covars_as), with=F]
  } else {
    load(covar_file)
    covar1 <- covar[year %in% years, c(area_var, "year", covars), with=F]
    load(covar_as_file)
    covar2 <- covar[year %in% years, c(area_var, "year", "sex", "age", covars_as), with=F]
    covar <- merge(covar1, covar2, by=c(area_var, "year"))
    rm(covar1, covar2); gc()
  }

# transform and standardize covariates
  if (!is.null(covars_trans)) {
    for (var in names(covars_trans)) {
      fun <- eval(parse(text = paste("function(x)", covars_trans[var])))
      covar[[var]] <- fun(covar[[var]])
    }
  }

  for (var in c(covars, covars_as)) {
    if (covar[, uniqueN(get(var))] == 2) next # don't standardize indicator covariates
    covar[[var]] <- (covar[[var]] - mean(covar[[var]])) / sd(covar[[var]])
  }
  setnames(covar, area_var, "area")

# merge covariates onto the dataset
  if (is.null(covars_as)) data <- merge(data, covar, by=c("area", "year"))
  else data <- merge(data, covar, by=c("area", "year", "sex", "age"))
  rm(covar); gc()
}

# recode year from 0
data[, year := as.integer(year - years[1])]

# recode age from 0
data[, age := as.integer(factor(age, levels=ages)) - 1L]

# create structure matrices for space, time, and age effects
load(adjmat_file)
graph_j <- diag(apply(adjmat, 1, sum)) - adjmat
graph_j <- as(graph_j, "dgTMatrix")
rm(adjmat)

num_t <- max(data$year) + 1
graph_t <- matrix(rep(0, num_t^2), nrow=num_t, ncol=num_t)
graph_t[abs(row(graph_t) - col(graph_t)) == 1] <- 1
graph_t <- diag(apply(graph_t, 1, sum)) - graph_t
graph_t <- as(graph_t, "dgTMatrix")

if (length(ages) > 1) {
  num_a <- max(data$age) + 1
  graph_a <- matrix(rep(0, num_a^2), nrow=num_a, ncol=num_a)
  graph_a[abs(row(graph_a) - col(graph_a)) == 1] <- 1
  graph_a <- diag(apply(graph_a, 1, sum)) - graph_a
  graph_a <- as(graph_a, "dgTMatrix")
} else {
  graph_a <- NULL
}

# If using completeness, save distributions here
if (grepl("_c", model)) {
  completeness_prior <- data.table(readRDS(completeness_prior_file))[, year := as.integer(year - years[1])][order(age_group, area, year)]
  completeness_prior <- completeness_prior[, C := 1:nrow(completeness_prior) - 1L]
  # MAke sure no SD is 0 
  completeness_prior[sd == 0, sd := 0.01]
  map  <- data[, .(area, year, age, sex)][, age_group := ifelse(age < 3, 0, 1)]
  
  if (!admin1_completeness) {
    
    #Completeness is at the area level
    map  <- merge(map, completeness_prior, by = c("area", "year", "age_group"), allow.cartesian = T)[, .(area, year, age, sex, C)]
    data <- merge(data, map, by = c("area", "year", "sex", "age"))
  } else {
    
    # Completeness is at the admin1 level
    setnames(completeness_prior, "area", "admin1")
    load(geoagg_files[["admin1"]])
    area_admin1 <- unique(weights[, .(area, admin1)])
    map <- merge(map, area_admin1, by = "area")
    map <- merge(map, completeness_prior, by = c("admin1", "year", "age_group"), allow.cartesian = T)[, .(area, year, age, sex, C)]
    data <- merge(data, map, by = c("area", "year", "sex", "age"))
  }
  
  # Save prepped data
  setkeyv(data, c("area", "year", "sex", "age"))
  save(data, completeness_prior, graph_j, graph_t, graph_a, file=paste0(temp_dir, "/data.rdata"))
} else {
  # save prepped data
  setkeyv(data, c("area", "year", "sex", "age"))
  save(data, graph_j, graph_t, graph_a, file=paste0(temp_dir, "/data.rdata")) 
}
