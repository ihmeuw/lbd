
root <- "<<<< FILEPATH REDACTED >>>>"

# Load packages
package_lib <- "<<<< FILEPATH REDACTED >>>>" # Library for all MBG versioned packages. Ensures that none of this code is dependent on the machine where the user runs the code.
.libPaths(package_lib) # Ensures packages look for dependencies here when called with library(). Necessary for seeg libraries.
library(INLA, lib.loc = package_lib) # The fmesher executables in the INLA library can't be run off the /share drives or J:/WORK. However, J:/temp works fine... some sort of permissions problem.
library(raster, lib.loc = package_lib)
library(seegMBG, lib.loc = package_lib)
library(data.table, lib.loc = package_lib)
library(rgdal, lib.loc = package_lib)

print(commandArgs())
indicator <- commandArgs()[3]
repo <- commandArgs()[4]
indicator_group <- commandArgs()[5]
run_date <- commandArgs()[6]

if(is.na(indicator)) {
  indicator <- 'edu_0'
  repo <- "<<<< FILEPATH REDACTED >>>>"
  indicator_group <- 'education'
  run_date <- '2016_09_11'
}

print(indicator)

setwd(repo)
source('mbg_central/functions.R')

# Load model input image
load("<<<< FILEPATH REDACTED >>>>")
# Process parameters from config file
config <- fread(paste0(indicator_group, '/config.csv'), header=FALSE)
for(param in config[, V1]) {
  assign(param, config[V1==param, V2])
}

# Set up model equation
f_null  <- formula('covered~-1+int')
f_lin <- reformulate(fixed_effects)

f_space <- ~ f(
  space,
  model = spde,
  group = space.group,
  control.group = list(model = 'ar1'))


f_mbg = f_null + f_lin + f_space

# construct an SPDE model with a Matern kernel
spde <- inla.spde2.matern(mesh = mesh_s,
                          alpha = 2)


# Projector Matrix
A <- inla.spde.make.A(
  mesh = mesh_s,
  loc = as.matrix(df[, c('longitude', 'latitude'),with=F]),
  group = df$period,
  group.mesh = mesh_t
)
space = inla.spde.make.index("space",
                             n.spde = spde$n.spde,
                             n.group = mesh_t$m)


#find cov indices
covs_indices <- unique(c(match(all.vars(f_lin), colnames(df))))

# make design matrix, center the scaling
design_matrix <- data.frame(int = 1,
                            df[,covs_indices,with=F])


cs_df <- getCentreScale(design_matrix,
                        exclude = c('int'))


design_matrix <- centreScale(design_matrix,
                             df = cs_df)

# construct a 'stack' object for observed data
cov=df[[indicator]] # N+_i
N=df$N                 # N_i

stack.obs <- inla.stack(
  data = list(covered = cov),
  A = list(A, 1),
  effects = list(space,
                 design_matrix),
  tag = 'est'
)



# get prior mean for the intercept
int_prior_mn <- 0

# ~~~~~~~~~~~~~~~~
# fit the model
# enable weights
inla.setOption("enable.inla.argument.weights", TRUE)

message('Fitting INLA model')

# set a prior variance of 1.96 on the intercept as this is
# roughly as flat as possible on logit scale without >1 inflection

# code just to fit the model (not predict)
if(keep_inla_files==FALSE) inla_working_dir <- "<<<< FILEPATH REDACTED >>>>"
if(keep_inla_files==TRUE) inla_working_dir <- "<<<< FILEPATH REDACTED >>>>"
system.time(
  res_fit <- inla(f_mbg,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           link = 1,
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       prec.intercept = 1,
                                       mean.intercept = int_prior_mn),
                  control.compute = list(dic = TRUE,
                                         cpo = TRUE,
                                         config = TRUE),
                  control.inla = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                  family = 'binomial',
                  num.threads = getOption('cores'),
                  Ntrials = N,
                  verbose = TRUE,
                  working.directory = inla_working_dir,
                  keep = TRUE)
)
# Clean up INLA intermediate directories unless user has specified to keep them.
inla_working_dirs <- list.dirs("<<<< FILEPATH REDACTED >>>>", recursive = FALSE)
inla_working_dirs <- inla_working_dirs[grepl("inla_scratch", inla_working_dirs)]
for(inla_dir in inla_working_dirs) {
  unlink(inla_dir, recursive = TRUE)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~
# predictions

message('Making predictions')

# number of samples
n_draws <- samples

# dummy raster
template <- raster("<<<< FILEPATH REDACTED >>>>")

cell_idx <- seegSDM:::notMissingIdx(template)

# sample from posterior over latents
suppressWarnings(draws <- inla.posterior.sample(n_draws, res_fit))

# get parameter names
par_names <- rownames(draws[[1]]$latent)

# index to spatial field and linear coefficient samples
s_idx <- grep('^s.*', par_names)
l_idx <- match(sprintf('%s.1', res_fit$names.fixed),
               par_names)


# get samples as matrices
pred_s <- sapply(draws, function (x)
  x$latent[s_idx])
pred_l <- sapply(draws, function (x)
  x$latent[l_idx])
if(length(l_idx)==1) pred_l=t(as.matrix(pred_l))
rownames(pred_l) <- res_fit$names.fixed

# get coordinates of cells to predict to
coords <- xyFromCell(template, seegSDM:::notMissingIdx(template))

# ~~~~~~~~~~~~~~
# project spatial effect


# make predictions for all periods

# get predictor matrix between spatial nodes and prediction locations
nperiod <- mesh_t$n

# replicate coordinates and years
coords_periods <- do.call(rbind,
                          replicate(nperiod,
                                    coords,
                                    simplify = FALSE))

groups_periods <- rep(1:nperiod,
                      each = nrow(coords))

# Projector matrix
A.pred <- inla.spde.make.A(
  mesh = mesh_s,
  loc = coords_periods,
  group = groups_periods,
  group.mesh = mesh_t
)


# get samples of s for all cells
cell_s <- A.pred %*% pred_s
cell_s <- as.matrix(cell_s)



# get sub-raster of covariates (remove int and otehrs not in there)
pars <- res_fit$names.fixed
# remove out temporally varying covariates for now, will deal with them later
tvnames=pars[!(pars %in% c('int',names(covs)))]
pars <- pars[(pars %in% names(covs))]

# extract
covs <- brick("<<<< FILEPATH REDACTED >>>>") # non-varying
message('Time-varying covariate rasters')
for(c in c('lights_new','evi','LST_day','total_pop','edu_0','edu_mean')) # time-varying
  if(c %in% fixed_effects) assign(paste0('tv_',c),brick(paste0(cov_dir, '/time_varying_covs_transformed_',c,'.grd')))

covs_sub <- covs[[pars]]
covs_sub = resample(covs_sub,template)
covs_sub = mask(covs_sub,template)

# extract cell values  non-temporally varying covariates
vals <- extract(covs_sub, coords)
if(dim(covs_sub)[3]==1)
  vals= data.frame("irrigation"=vals)


# add intercept term and dummy variables for other terms
vals <- data.frame(cbind(int = 1, vals))
for(n in tvnames)   vals[,n]=0


# apply centreing/scaling
vals_cs <- centreScale(vals, cs_df)
vals_cs[,tvnames]=0
# put names back to where they were

# get order of predictors
vals_cs = as.matrix(vals_cs[,c(rownames(pred_l))])

# get predictions (Without time varying covariates for now..)
cell_l <- vals_cs %*% pred_l


# replicate over periods
cell_l <- do.call(rbind,
                  replicate(nperiod,
                            cell_l,
                            simplify = FALSE))



#Adding in temporally varying covariate effects
for(tv in ls()[grep('tv_',ls())]){

  message(tv)

  #print(tv)
  r <- get(tv)

  #rates_all <- extract(rates_ras, cell_idx)
  r = resample(r,template)
  r = mask(r,template)

  d <- extract(r, coords)

  # turn to vector
  d_vals = as.vector(d)
  n = gsub('tv_','',tv)

  # scale them
  tmp=data.frame(d_vals)
  colnames(tmp)=n
  d_vals <- centreScale(tmp,
                        cs_df)[,n]

  # multiply by parameter draws
  d_effect <- d_vals %*% t(pred_l[n, ])

  # add to linear component
  cell_l <- cell_l + d_effect

  rm(d_effect)

}


# ~~~~~~~~~~~~~~
# project model uncertainty from node-level sd
cell_l_sd <-  apply(cell_l, 1, sd)



node_s_sd <- apply(pred_s, 1, sd)
cell_s_sd <- as.matrix(A.pred %*% node_s_sd)
# combine to get cell-level sd

# project to cells
cell_sd <- sqrt(cell_l_sd ^ 2 + cell_s_sd ^ 2)

# ~~~~~~~~~~~~~~
# combine and summarise them
cell_all <- cell_l + cell_s


# get predictive draws on probability scale
cell_pred <- plogis(as.matrix(cell_all))


# get prediction mean (integrated probability)
pred_mean <- rowMeans(cell_pred)


# create multi-band rasters for each of these metrics
# each band giving the values for a different period
mean_ras <- insertRaster(template,
                         matrix(pred_mean,
                                ncol = nperiod))

sd_ras <- insertRaster(template,
                       matrix(cell_sd,
                              ncol = nperiod))

names(mean_ras) <-
names(sd_ras) <-
  paste0('period_', 1:nperiod)






# ~~~~~~~~~~~~
# save the outputs
if(time_stamp==TRUE) output_dir <- "<<<< FILEPATH REDACTED >>>>"
if(time_stamp==FALSE) output_dir <- "<<<< FILEPATH REDACTED >>>>"
dir.create(output_dir)

# Save log of config file
fwrite(config, paste0(output_dir,'/config.csv'))

writeRaster(
  mean_ras,
  file = (paste0(output_dir, '/', indicator,'_prediction_eb')),
  overwrite = TRUE
)
# latent sd
writeRaster(
  sd_ras,
  file = (paste0(output_dir, '/', indicator,'_sd_eb')),
  overwrite = TRUE
)
# save model
save(res_fit,
     file = (paste0(output_dir, '/', indicator,'_model_eb.RData')))
# save draws (with compression) to recombine later
save(
  cell_pred,
  file = (paste0(output_dir, '/', indicator,'_cell_draws_eb.RData')),
  compress = TRUE
)
# save training data
write.csv(
  df,
  file = (paste0(output_dir, '/', indicator,'_trainingdata.RData')),
  row.names = FALSE
)


