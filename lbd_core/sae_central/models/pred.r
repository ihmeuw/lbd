####################################################################################################
## Description: Generate draws from a fitted SAE model for all areas/ages/years for a given sex.
##              Note that all draws from the model are generated simultaneously so that the
##              correlation structure is appropriately reflected. Also generates point estimates,
##              confidence intervals, and standard errors for all areas/ages/years.
##
## Passed args: main_dir [character] -- home directory for settings and final output
##              sex [integer] -- sex to generate predictions for
##
## Requires:    fitted model object ('model_fit_[sex].rdata' in temp_dir)
##              prepped data file ('data.rdata' in temp_dir)
##              populations (pop_file)
##              age standard file (age_std_file)
##
## Outputs:     draws for all areas/ages/years, saved in separate files by year
##                ('draws_[area_var]_[year]_[sex].rdata' in temp_dir)
##
##              point estimates, confidence intervals, and standard errors, saved in
##                separate files by year
##                ('est_[area_var]_[year]_[sex].rdata' in temp_dir)
##
## Run from within 'sae_central' directory!!
####################################################################################################

library(Matrix)
library(data.table)

rm(list=ls())
set.seed(98121)

## Get settings and functions ----------------------------------------------------------------------
main_dir <- commandArgs()[4]
sex <- as.integer(commandArgs()[5])

source("settings.r")
get_settings(main_dir)

source("post_estimation/calc_all_ages.r")
source("post_estimation/collapse.r")

## Load data and model fit -------------------------------------------------------------------------
load(paste0(temp_dir, "/data.rdata"))
data <- data[sex == get("sex", .GlobalEnv),]
rm(list=grep("graph", ls(), value=T))

load(paste0(temp_dir, "/model_fit_", sex, ".rdata"))

## Generate draws ----------------------------------------------------------------------------------
simvars <- paste0("V", 1:n.sims)

# extract mean and precision matrix for all model parameters
mean <- out$par.random
ii <- grep("^B|^re", rownames(out$jointPrecision))
prec <- out$jointPrecision[ii,ii]
rm(out, ii); gc()

# draw from a multivariate normal distribution given the mean and precision matrix for all parameters
gen_sims <- function(mu, prec, n.sims) {
  z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  L_inv = Cholesky(prec)
  mu + solve(as(L_inv, "pMatrix"), solve(t(as(L_inv, "Matrix")), z))
}
sims <- gen_sims(mu=mean, prec=prec, n.sims=n.sims)
rm(prec, gen_sims); gc()

# calculate draws for the fixed portion of the model (same for all models)
B <- sims[grepl("^B", names(mean)),]
if (class(B) == "numeric") B <- matrix(B, nrow=1) else B <- as.matrix(B)
data[, int := 1]
fe <- as.matrix(data[, c("int", covars, covars_as), with=F]) %*% B
rm(B); gc()

# calculate draws for the random portion of the model (differs by model)
if (model == "1") {
  # pull out draws for the age-year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := rep(1:num_t, each = num_a) - 1L]
  re1[, age := rep(1:num_a, num_t) - 1L]
  setkeyv(re1, c("age", "year"))

  # pull out draws for the area-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := 1:num_j - 1]
  setkeyv(re2, "area")

  # pull out draws for the area-level random slope on year
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1]
  setkeyv(re3, "area")

  # pull out draws for the area-level random slope on age
  re4 <- data.table(as.matrix(sims[names(mean) == "re4",]))
  re4[, area := 1:num_j - 1]
  setkeyv(re4, "area")

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[data[, list(age, year)], simvars, with=F]) +
    as.matrix(re2[J(data$area), simvars, with=F]) +
    data$year * as.matrix(re3[J(data$area), simvars, with=F]) +
    data$age * as.matrix(re4[J(data$area), simvars, with=F])
  rm(re1, re2, re3, re4); gc()
}

if (model == "1b") {
  # pull out draws for the year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := 1:num_t - 1]
  setkeyv(re1, "year")

  # pull out draws for the area-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := 1:num_j - 1]
  setkeyv(re2, "area")

  # pull out draws for the area-level random slope on year
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1]
  setkeyv(re3, "area")

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[J(data$year), simvars, with=F]) +
    as.matrix(re2[J(data$area), simvars, with=F]) +
    data$year * as.matrix(re3[J(data$area), simvars, with=F])
  rm(re1, re2, re3); gc()
}

if (model == "2") {
  # pull out draws for the age-year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := rep(1:num_t, each = num_a) - 1L]
  re1[, age := rep(1:num_a, num_t) - 1L]
  setkeyv(re1, c("age", "year"))

  # pull out draws for the area-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := 1:num_j - 1]
  setkeyv(re2, "area")

  # pull out draws for the area-level random slope on year
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1]
  setkeyv(re3, "area")

  # pull out draws for the area-level random slope on age
  re4 <- data.table(as.matrix(sims[names(mean) == "re4",]))
  re4[, area := 1:num_j - 1]
  setkeyv(re4, "area")

  # pull out draws from the area-age-year-level random intercept
  re5 <- data.table(as.matrix(sims[names(mean) == "re5",]))
  re5[, area := rep(1:num_j, num_a * num_t) - 1]
  re5[, age := rep(rep(1:num_a, each = num_j), num_t) - 1]
  re5[, year := rep(1:num_t, each = num_j * num_a) - 1]
  setkeyv(re5, c("area", "age", "year"))

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[data[, list(age, year)], simvars, with=F]) +
    as.matrix(re2[J(data$area), simvars, with=F]) +
    data$year * as.matrix(re3[J(data$area), simvars, with=F]) +
    data$age * as.matrix(re4[J(data$area), simvars, with=F]) +
    as.matrix(re5[data[, list(area, age, year)], simvars, with=F])
  rm(re1, re2, re3, re4, re5); gc()
}

if (model == "2b") {
  # pull out draws for the year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := 1:num_t - 1]
  setkeyv(re1, "year")

  # pull out draws for the area-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := 1:num_j - 1]
  setkeyv(re2, "area")

  # pull out draws for the area-level random slope on year
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1]
  setkeyv(re3, "area")

  # pull out draws from the area-year-level random intercept
  re5 <- data.table(as.matrix(sims[names(mean) == "re5",]))
  re5[, area := rep(1:num_j, num_t) - 1]
  re5[, year := rep(1:num_t, each = num_j) - 1]
  setkeyv(re5, c("area", "year"))

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[J(data$year), simvars, with=F]) +
    as.matrix(re2[J(data$area), simvars, with=F]) +
    data$year * as.matrix(re3[J(data$area), simvars, with=F]) +
    as.matrix(re5[data[, list(area, year)], simvars, with=F])
  rm(re1, re2, re3, re5); gc()
}

if (model == "3") {
  # pull out draws for the age-year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := rep(1:num_t, each = num_a) - 1L]
  re1[, age := rep(1:num_a, num_t) - 1L]
  setkeyv(re1, c("age", "year"))

  # pull out draws for the area-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := 1:num_j - 1]
  setkeyv(re2, "area")

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[data[, list(age, year)], simvars, with=F]) +
    as.matrix(re2[J(data$area), simvars, with=F])
  rm(re1, re2); gc()
}

if (model == "4") {
  # pull out draws for the age-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, age := 1:num_a - 1]
  setkeyv(re1, "age")

  # pull out draws for the year-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, year := 1:num_t - 1]
  setkeyv(re2, "year")

  # pull out draws for the area-level random intercept
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1]
  setkeyv(re3, "area")

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[J(data$age), simvars, with=F]) +
    as.matrix(re2[J(data$year), simvars, with=F]) +
    as.matrix(re3[J(data$area), simvars, with=F])
  rm(re1, re2, re3); gc()
}

if (model == "5") {
  # pull out draws for the age-year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := rep(1:num_t, each = num_a) - 1L]
  re1[, age := rep(1:num_a, num_t) - 1L]
  setkeyv(re1, c("age", "year"))

  # pull out draws for the area-year-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := rep(1:num_j, num_t) - 1]
  re2[, year := rep(1:num_t, each = num_j) - 1]
  setkeyv(re2, c("area", "year"))

  # pull out draws for the area-level random slope on age
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1]
  setkeyv(re3, "area")

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[data[, list(age, year)], simvars, with=F]) +
    as.matrix(re2[data[, list(area, year)], simvars, with=F]) +
    data$age * as.matrix(re3[J(data$area), simvars, with=F])
  rm(re1, re2, re3); gc()
}

if (model == "6") {
  # pull out draws for the age-year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := rep(1:num_t, each = num_a) - 1L]
  re1[, age := rep(1:num_a, num_t) - 1L]
  setkeyv(re1, c("age", "year"))

  # pull out draws for the area-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := 1:num_j - 1L]
  setkeyv(re2, "area")

  # pull out draws for the area-level random slope on year
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1L]
  setkeyv(re3, "area")

  # pull out draws for the area-level random slope on age
  re4 <- data.table(as.matrix(sims[names(mean) == "re4",]))
  re4[, area := 1:num_j - 1L]
  setkeyv(re4, "area")

  # pull out draws from the area-year-level random intercept
  re5 <- data.table(as.matrix(sims[names(mean) == "re5",]))
  re5[, area := rep(1:num_j, num_t) - 1L]
  re5[, year := rep(1:num_t, each = num_j) - 1L]
  setkeyv(re5, c("area", "year"))

  # pull out draws from the area-age-level random intercept
  re6 <- data.table(as.matrix(sims[names(mean) == "re6",]))
  re6[, area := rep(1:num_j, num_a) - 1L]
  re6[, age := rep(1:num_a, each = num_j) - 1L]
  setkeyv(re6, c("area", "age"))

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[data[, list(age, year)], simvars, with=F]) +
    as.matrix(re2[J(data$area), simvars, with=F]) +
    data$year * as.matrix(re3[J(data$area), simvars, with=F]) +
    data$age * as.matrix(re4[J(data$area), simvars, with=F]) +
    as.matrix(re5[data[, list(area, year)], simvars, with=F]) +
    as.matrix(re6[data[, list(area, age)], simvars, with=F])
  rm(re1, re2, re3, re4, re5, re6); gc()
}

if (model == "7") {
  # pull out draws for the age-year-level random intercept
  re1 <- data.table(as.matrix(sims[names(mean) == "re1",]))
  re1[, year := rep(1:num_t, each = num_a) - 1L]
  re1[, age := rep(1:num_a, num_t) - 1L]
  setkeyv(re1, c("age", "year"))

  # pull out draws for the area-level random intercept
  re2 <- data.table(as.matrix(sims[names(mean) == "re2",]))
  re2[, area := 1:num_j - 1L]
  setkeyv(re2, "area")

  # pull out draws for the area-level random slope on year
  re3 <- data.table(as.matrix(sims[names(mean) == "re3",]))
  re3[, area := 1:num_j - 1L]
  setkeyv(re3, "area")

  # pull out draws from the area-year-level random intercept
  re5 <- data.table(as.matrix(sims[names(mean) == "re5",]))
  re5[, area := rep(1:num_j, num_t) - 1L]
  re5[, year := rep(1:num_t, each = num_j) - 1L]
  setkeyv(re5, c("area", "year"))

  # combine all random effects
  rm(sims, mean); gc()
  re <- as.matrix(re1[data[, list(age, year)], simvars, with=F]) +
    as.matrix(re2[J(data$area), simvars, with=F]) +
    data$year * as.matrix(re3[J(data$area), simvars, with=F]) +
    as.matrix(re5[data[, list(area, year)], simvars, with=F])
  rm(re1, re2, re3, re5); gc()
}

# calculate the rate and combine with identifying information
draws <- fe + re
rm(fe, re); gc()

draws <- data.table(exp(draws))
draws[, area := as.integer(data$area)]
draws[, year := as.integer(data$year + years[1])]
draws[, sex := as.integer(data$sex)]
draws[, age := as.integer(as.character(factor(data$age, labels=ages)))]
rm(data); gc()

## Split out draws by year, calculate all-ages and age-standardized draws, collapse, and save ------
# load and subset the population file, then merge onto draws
load(pop_file)
setnames(pop, area_var, "area")
pop <- pop[sex == get("sex", .GlobalEnv) & year %in% years,
           list(pop = as.numeric(sum(pop))), keyby='area,year,sex,age']
draws <- merge(pop, draws, by=c("area", "year", "sex", "age"), all=T)
rm(pop); gc()

# load the age standard
std_wt <- fread(age_std_file)[, list(age, wt = wt/sum(wt))]

# split up data and process separately by year
all_draws <- draws
rm(draws)

all_draws <- lapply(years, function(this_year) all_draws[year == this_year,])
for (this_year in years) {
  cat(paste(Sys.time(), this_year, "\n"))
  draws <- all_draws[[which(years == this_year)]]
  draws[, level := area_var]

  # reshape long for easier processing
  draws <- melt(draws, id.vars=c("level", "area", "year", "sex", "age", "pop"), variable.name="sim")
  draws[, sim := as.integer(sim)]
  draws[is.na(value), value := 0]

  # add all-ages rates
  draws <- calc_all_ages(draws, std_wt)

  # calculate point estimates and CIs
  est <- collapse_draws(draws, c("level", "area", "year", "sex", "age"))

  # save draws and estimates
  save(draws, file=paste0(temp_dir, "/draws_", area_var, "_", this_year, "_", sex, ".rdata"))
  save(est, file=paste0(temp_dir, "/est_", area_var, "_", this_year, "_", sex, ".rdata"))

  # clear some memory
  all_draws[[which(years == this_year)]] <- NA
  rm(draws, est); gc()
}
