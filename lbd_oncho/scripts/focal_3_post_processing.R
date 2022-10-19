########################################################################################################
### Miscellaenous post-processing code for Focal 3 models
########################################################################################################

################################################ Setup #################################################
library(data.table)
library(INLA)
library(ggplot2)


######################## 1. Plot random effects for oncho suitability random walk model ################

run_date <- "2020_08_08_15_01_40"
model_fit <- readRDS(<<<< FILEPATH REDACTED >>>>)
ggplot(data = model_fit$summary.random$`inla.group(oncho_suitability, n = 25)`) + theme_classic() + geom_line(aes(x = ID, y = mean)) + geom_ribbon(aes(x = ID, ymin = `0.025quant`, ymax = `0.975quant`), alpha = 0.2) + geom_abline(aes(intercept = 0, slope = 0), linetype = "dotted") + labs(x = "Normalized Suitability Index", y = "Random Effect")


######################## 2. Output INLA priors and spatial hyperparameter estimates ####################

run_date <- "2020_08_08_15_01_40"
model_fit <- readRDS(<<<< FILEPATH REDACTED >>>>)
load(<<<< FILEPATH REDACTED >>>>)

build_spde_prior <- function(spde_prior, mesh_s, st_gp_int_zero) {
  if(spde_prior$type=="pc"){ # PC prior
    if(is.null(spde_prior$prior$sigma)) {
      spde_prior$prior$sigma <- c(3, 0.05) # P(sigma > 3) = 0.05
    }
    if(is.null(spde_prior$prior$range)) {
      mesh.range <- max(c(diff(range(mesh_s$loc[, 1])), 
                          diff(range(mesh_s$loc[, 2])), 
                          diff(range(mesh_s$loc[, 3]))))
      spde_prior$prior$range <- c(mesh.range*0.05, 0.05) # P(range < 5% max extent of mesh) = 0.05
    }
    message(paste("Building spde with pc prior,",
                  spde_prior$prior$range[2]*100,
                  "% probability that the range is lower than",
                  spde_prior$prior$range[1], 
                  "and a",
                  spde_prior$prior$sigma[2]*100,
                  "% probability that sigma is greater than",
                  spde_prior$prior$sigma[1]))
    spde <- inla.spde2.pcmatern(mesh = mesh_s, 
                                alpha = 2,
                                prior.range = spde_prior$prior$range,
                                prior.sigma = spde_prior$prior$sigma,
                                constr = st_gp_int_zero)
  } else { # Non PC prior
    if(is.null(spde_prior$prior$variance.nominal)) {
      spde_prior$prior$variance.nominal <- 1
    }
    spde <- inla.spde2.matern(mesh = mesh_s,  alpha = 2, constr = st_gp_int_zero,
                              prior.range.nominal = spde_prior$prior$range.nominal,
                              prior.variance.nominal = spde_prior$prior$variance.nominal)
  }
  return(list(spde=spde, spde_prior=spde_prior))
}

spde_prior <- "list(type='pc')"
spde_prior <- eval(parse(text=spde_prior)) # convert from string to list
st_gp_int_zero <- FALSE
spde_list <- build_spde_prior(spde_prior, mesh_s, st_gp_int_zero)
spde_prior <- spde_list$spde_prior
spde_prior
spde <- spde_list$spde
spde$param.inla$theta.initial

# res.field <- inla.spde2.result(model_fit, 'space', spde, do.transf=TRUE)
res.field <- inla.spde2.result(model_fit, 'sp.no.t', spde, do.transf=TRUE)

## nominal range at 0.025, 0.5, 0.975 quantiles
range   <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.range.nominal[[1]]) * 6371
nom.var <- inla.qmarginal(c(0.025, 0.5, 0.975), res.field$marginals.variance.nominal[[1]])
spat.hyps <- rbind(range, nom.var)
rownames(spat.hyps) <- c('Nominal Range', 'Nominal Variance')
spat.hyps


######################## 3. Save stacker rasters ####################
library(data.table)
library(raster)

run_date <- "2020_05_20_21_52_36"

load(<<<< FILEPATH REDACTED >>>>)

gam_2016 <- cov_list$gam[[27]]
gbm_2016 <- cov_list$gbm[[27]]
lasso_2016 <- cov_list$lasso[[27]]
values(gbm_2016)[values(gbm_2016) >= 1 & !is.na(values(gbm_2016))] <- 0.999

writeRaster(gam_2016, <<<< FILEPATH REDACTED >>>>, format = "GTiff")
writeRaster(gbm_2016, <<<< FILEPATH REDACTED >>>>, format = "GTiff")
writeRaster(lasso_2016, <<<< FILEPATH REDACTED >>>>, format = "GTiff")


######################## 4. Output partial dependence plot and lasso coefficients ####################
library(INLA)
library(data.table)
library(ggplot2)
library(dismo)
library(gbm)
library(glmnet)

load(<<<< FILEPATH REDACTED >>>>)
gbm <- child_models[[2]]
summary(gbm)
gbm.plot(gbm)

covlist <- c("crutswet", "mswep", "onchomda", "crutstmx", "crutsdtr", "evi_v6", "tcb_v6", "sgbdrlog", "worldpop", "river_size", "lfmda", "sgcrfvol")
results <- data.table("x"=NA, "y"=NA, "covariate"=NA)
results <- results[-1,]
for (a in covlist) {
  output <- plot.gbm(gbm, i.var=a, return.grid=TRUE, type="link", continuous.resolution=1000)
  output$covariate <- a
  colnames(output)[1] <- "x"
  results <- rbind(results, output)
}

ggplot(data=transform(results, covariate=factor(covariate, levels=covlist)), aes(x=x, y=y)) + theme_classic() + geom_line(color="#666666") + facet_wrap(covariate ~ ., nrow = 4, ncol = 3, scales = "free") + labs(title="GBM Partial Dependence\n", x="\nStandardized covariate value", y="Partial dependence\n") + theme(plot.title=element_text(hjust=0.5, face="bold"))# + theme(axis.label.y=element_text(size=10, margin=margin(t=150)), panel.spacing=unit(1, "lines"))

lasso <- child_models[[3]]
plot(lasso)
l <- lasso$cv_1se_lambda
l.idx <- which(lasso$lambda == l)
lasso$beta[, l.idx]

######################## 5. Determine number of data rows used in model ####################
library(data.table)

df <- fread(<<<< FILEPATH REDACTED >>>>)

points <- nrow(unique(df[point == 1, c("Master_UID", "point", "latitude", "longitude")]))
polys <- nrow(unique(df[point == 0, c("Master_UID", "point", "shapefile", "location_code")]))
all <- points + polys
c(all, points, polys)

table(unique(df[point == 1, c("Master_UID", "point", "latitude", "longitude", "diagnostic")])$diagnostic) + table(unique(df[point == 0, c("Master_UID", "point", "shapefile", "location_code", "diagnostic")])$diagnostic)

by_country <- merge(as.data.table(table(unique(df[point == 1, c("Master_UID", "point", "latitude", "longitude", "country")])$country)), as.data.table(table(unique(df[point == 0, c("Master_UID", "point", "shapefile", "location_code", "country")])$country)), by = "V1", all = TRUE)
by_country[is.na(N.x), N.x := 0]
by_country[is.na(N.y), N.y := 0]
by_country$total <- by_country$N.x + by_country$N.y
sum(by_country$total)

by_year <- merge(as.data.table(table(unique(df[point == 1, c("Master_UID", "point", "latitude", "longitude", "year")])$year)), as.data.table(table(unique(df[point == 0, c("Master_UID", "point", "shapefile", "location_code", "year")])$year)), by = "V1", all = TRUE)
by_year[is.na(N.x), N.x := 0]
by_year[is.na(N.y), N.y := 0]
by_year$total <- by_year$N.x + by_year$N.y
sum(by_year$total)
