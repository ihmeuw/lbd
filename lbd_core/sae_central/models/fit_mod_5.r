####################################################################################################
## Description: Fit the following model in TMB:
##
##              D_{j,a,t} ~ Poisson(m_{j,a,t} * P_{j,a,t})
##              log(m_{j,a,t}) = B0 + B1*X_{j,t} + u1_{a,t} + u2_{j,t} + u3_{j}*a
##              u1_{a,t} ~ LCAR:LCAR(rho_1a, rho_1t, sigma_1)
##              u2_{j,t} ~ LCAR(rho_2j, rho_2t, sigma_2)
##              u3_{j} ~ LCAR(rho_3, sigma_3)
##              where j = area, a = age group, t = year, D = events, P = population, m = event rate
##
## Passed args: main_dir [character] -- home directory for settings and final output
##              sex [integer] -- sex to run model for
##
## Requires:    prepped data file ('data.rdata' in temp_dir)
##
## Outputs:     fitted model object ('model_fit_[sex].rdata' in temp_dir)
##              model-fitting log ('model_fitting_[sex].txt' in temp_dir)
##
## Run from within 'sae_central' directory!!
####################################################################################################

library(TMB)
library(optimx)
library(data.table)

rm(list=ls())
set.seed(98121)

## Get settings ------------------------------------------------------------------------------------
main_dir <- commandArgs()[4]
sex <- as.integer(commandArgs()[5])

source("settings.r")
get_settings(main_dir)

sink(file=paste0(temp_dir, "/model_fitting_", sex, ".txt"), type="output", split=T)

## Format data for TMB -----------------------------------------------------------------------------
cat("\n\n***** Format data\n"); flush.console()
load(paste0(temp_dir, "/data.rdata"))
data <- data[sex == get("sex", .GlobalEnv),]
data[, int := 1L]

num_j <- max(data$area) + 1
num_t <- max(data$year) + 1
num_a <- max(data$age) + 1

ii <- which(data$pop > 0)
tmb_data <- list(
  Y = data$events[ii],
  N = data$pop[ii],
  J = data$area[ii],
  A = data$age[ii],
  T = data$year[ii],
  X = as.matrix(data[ii, c("int", covars, covars_as), with=F]),
  graph_j = graph_j,
  graph_t = graph_t,
  graph_a = graph_a)

## Set parameters for TMB --------------------------------------------------------------------------
cat("\n\n***** Set parameters\n"); flush.console()
tmb_par <- list(
  B = rep(0, length(covars) + length(covars_as) + 1), # intercept & covariate effects

  re1=array(rep(0, num_a * num_t), dim=c(num_a, num_t)), # age-year-level random intercept
  log_sigma2_1=-5,
  logit_rho_1a=0,
  logit_rho_1t=0,

  re2=array(rep(0, num_j * num_t), dim=c(num_j, num_t)), # area-year-level random intercept
  log_sigma2_2=-5,
  logit_rho_2j=0,
  logit_rho_2t=0,

  re3=rep(0, num_j), # area-level random slope on age
  log_sigma2_3=-5,
  logit_rho_3=0)

## Fit model ---------------------------------------------------------------------------------------
# compile CPP code for objective function
openmp(16)
TMB::compile("models/mod_5.cpp")
dyn.load(dynlib("models/mod_5"))
config(tape.parallel=0, DLL="mod_5")

# make objective function
cat("\n\n***** Make objective function\n"); flush.console()
obj <- MakeADFun(tmb_data, tmb_par, random=c("B", paste0("re", 1:3)), DLL="mod_5")

# optimize objective function
cat("\n\n***** Optimize objective function\n"); flush.console()
opt_time <- proc.time()
for (method in c("nlminb", "L-BFGS-B", "BFGS", "Nelder-Meade", "CG")) {
  cat(paste("\n  *** Method: ", method, "\n"))
  opt <- optimx(par=obj$par, fn=function(x) as.numeric(obj$fn(x)), gr=obj$gr,
                method=method, control=list(iter.max=500, eval.max=500))
  print(opt)
  if (opt$convcod == 0) break
}
(opt_time <- proc.time() - opt_time)
if (opt$convcod != 0) stop("Model did not converge")

# get standard errors
cat("\n\n***** Extract standard errors\n"); flush.console()
se_time <- proc.time()
out <- sdreport(obj, getJointPrecision=T)
(se_time <- proc.time() - se_time)

# save model output
cat("\n\n***** Save model output\n"); flush.console()
save(out, num_j, num_t, num_a, opt_time, se_time, file=paste0(temp_dir, "/model_fit_", sex, ".rdata"))
sink()
