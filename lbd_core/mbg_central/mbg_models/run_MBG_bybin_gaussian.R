
root <- "<<<< FILEPATH REDACTED >>>>"
hroot <- "<<<< FILEPATH REDACTED >>>>"
setwd("<<<< FILEPATH REDACTED >>>>")
source('functions.R')
inla_working_dir <- "<<<< FILEPATH REDACTED >>>>"

defaultOptions(resolution = 5,        # raster resolution
               location = 'seattle',  # location for final run
               cores = 30,            # number of cores to use
               start = Sys.time())    # start time


# reload packages, just in case
library(INLA, lib.loc = "<<<< FILEPATH REDACTED >>>>")
library(raster, lib.loc = "<<<< FILEPATH REDACTED >>>>")
library(seegMBG, lib.loc = "<<<< FILEPATH REDACTED >>>>")
library(data.table)
library(rgdal)

## Logit functions
logit <- function(x) {
  log(x/(1-x))
}
invlogit <- function(x) {
  exp(x)/(1+exp(x))
}

print(commandArgs())
indicator <- commandArgs()[4]
print(indicator)

f_null  <- formula('covered~-1+int')
f_lin <- reformulate(paste0("irrigation"," + ","lights_new"))

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
unlink(inla_working_dir, recursive = T)
dir.create(inla_working_dir)
system.time(
  res_fit <- inla(f_mbg,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       prec.intercept = 1,
                                       mean.intercept = int_prior_mn),
                  control.compute = list(dic = TRUE,
                                         cpo = TRUE,
                                         config = TRUE),
                  control.inla = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                  family = 'Gaussian',
                  num.threads = getOption('cores'),
                  verbose = TRUE,
                  working.directory = inla_working_dir,
                  keep = TRUE)
)
unlink(inla_working_dir)

# ~~~~~~~~~~~~~~~~~~~~~~~~
# predictions

message('Making predictions')

# number of samples
n_draws <- 100

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
#covs <- brick(paste0(root,'/WORK/01_covariates/02_inputs/education/update_2017/data/geospatial_data/data/',target_country,'/clean/covs_transformed.grd'))
covs <- brick("<<<< FILEPATH REDACTED >>>>") # non-varying
message('Time-varying covariate rasters')
# for(c in c('total_pop','lights_new'))
#   assign(paste0('tv_',c),brick(paste0(root,'/WORK/01_covariates/02_inputs/education/update_2017/data/geospatial_data/data/',target_country,'/clean/time_varying_covs_transformed_',c,'.grd')))
#for(c in c('lights_new','total_pop')) # time-varying
for(c in c('lights_new')) # time-varying
  assign(paste0('tv_',c),brick("<<<< FILEPATH REDACTED >>>>"))

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
  r <- get(tv)

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
cell_pred <- (invlogit(as.matrix(cell_all)))*18

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

## Test taking probability threshold of candidate maps (summary of full posterior)
# Probability that over 50% of women 15-54 have 0 years of education.
cell_pred.dt <- as.data.table(cell_pred)
cols <- names(cell_pred.dt)
cell_pred.dt[, (cols) := lapply(.SD, function(x){ifelse((x)>6,1,0)}), .SDcols=cols]
over6_mean_edu <- rowMeans(cell_pred.dt)

over6_mean_edu_ras <- insertRaster(template,
                                   matrix(over6_mean_edu,
                                          ncol = nperiod))

names(over6_mean_edu_ras) <-
  names(mean_ras) <-
  names(sd_ras) <-

  paste0('period_', 1:nperiod)






# ~~~~~~~~~~~~
# save the outputs


# output the prediction samples as a GeoTiff
writeRaster(
  mean_ras,
  file = ("<<<< FILEPATH REDACTED >>>>"),
  overwrite = TRUE
)
# probability threshold raster
writeRaster(
  over6_mean_edu_ras,
  file = "<<<< FILEPATH REDACTED >>>>",
  overwrite = TRUE
)
# latent sd
writeRaster(
  sd_ras,
  file = "<<<< FILEPATH REDACTED >>>>",
  overwrite = TRUE
)

# save model
save(res_fit,
     file = "<<<< FILEPATH REDACTED >>>>")


# save draws (with compression) to recombine later
save(
  cell_pred,
  file = "<<<< FILEPATH REDACTED >>>>",
  compress = TRUE
)


# save training data
write.csv(
  df,
  file = "<<<< FILEPATH REDACTED >>>>",
  row.names = FALSE
)


