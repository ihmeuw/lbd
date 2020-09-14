########################################################################################################################################
###### Plot prior and posterior distributions of spatial hyperparameters (Theta1 and Theta2 for space) from MBG INLA fit.
###### Note: This code currently assumes that Theta1 and Theta2 are the first two hyperparameters in the INLA object; adjust indices
###### if this is not the case, or change code to explicitly call these hyperparameters by name.
########################################################################################################################################

indicator <- 'has_lri'
indicator_group <- 'lri'
run_date <- '2019_09_16_16_39_00'
region <- 'dia_name-ESH'
holdout <- 0
age <- 0

### Turning Chris's code into a function
plot_spatial_priors <- function(indicator,
                                indicator_group,
                                run_date,
                                region,
                                holdout = 0,
                                age = 0) {
  
  
  ### Load required libraries
  library(INLA)
  library(MASS)
  library(ggplot2)
  library(actuar, lib.loc = '<<<< FILEPATH REDACTED >>>>')
  library(mvtnorm)
  library(data.table)
  
  
  ### Load a fitted INLA model object
  fit <- readRDS('<<<< FILEPATH REDACTED >>>>')
  
  ### Plot priors and posteriors using INLA, for comparison
  plot(fit, plot.fixed.effects=FALSE, plot.lincomb=FALSE, plot.random.effects=FALSE, plot.predictor=FALSE, plot.q=FALSE, plot.cpo=FALSE, plot.prior=TRUE)
  
  ### Retrieve index of "space" ----------------------------------------------------------------------------------------------------------
  index_space <- NA
  for (a in 1:length(fit$all.hyper$random)) {
    if (fit$all.hyper$random[[a]][["hyperid"]] == "space") {
      index_space <- a
      break
    }
  }
  
  ### Retrieve prior means and precisions for Theta1 and Theta2 (provided in param for theta1)
  means_space <- fit$all.hyper$random[[index_space]]$hyper$theta1$param[1:2]
  covar_space <- 1/fit$all.hyper$random[[index_space]]$hyper$theta1$param[3:6] ## Convert precision to variance
  covar_space[covar_space == Inf] <- 0
  
  ### Sample from multivariate normal distribution for Theta1 and Theta2
  priors_space <- as.data.frame(mvrnorm(n=10000, means_space, matrix(covar_space, 2, 2)))
  colnames(priors_space) <- c("Theta1.Prior", "Theta2.Prior")
  
  ### Sample from posterior distributions of all hyperparameters
  posteriors <- data.frame(inla.hyperpar.sample(10000, fit, improve.marginals=TRUE))

  ### Combine prior and posterior samples for Theta1, and plot
  Theta1 <- as.data.frame(rbind(cbind(priors_space[, 1], Distribution="Prior"), cbind(posteriors[, "Theta1.for.space"], Distribution="Posterior")))
  colnames(Theta1)[1] <- "Value"
  Theta1$Value <- as.numeric(levels(Theta1$Value))[Theta1$Value]
  plot_th1 <- ggplot() + theme_classic() + stat_density(data=Theta1, aes(x=Value, group=Distribution, color=Distribution), position='identity', geom='line') + labs(title=paste0(toupper(indicator), ' theta1 [log(Tau)] for space in ', region)) + theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(expand=TRUE)
  
  ### Combine prior and posterior samples for Theta2, and plot
  Theta2 <- as.data.frame(rbind(cbind(priors_space[, 2], Distribution="Prior"), cbind(posteriors[, "Theta2.for.space"], Distribution="Posterior")))
  colnames(Theta2)[1] <- "Value"
  Theta2$Value <- as.numeric(levels(Theta2$Value))[Theta2$Value]
  plot_th2 <- ggplot() + theme_classic() + stat_density(data=Theta2, aes(x=Value, group=Distribution, color=Distribution), position='identity', geom='line') + labs(title=paste0(toupper(indicator), ' theta2 [log(Kappa)] for space in ', region)) + theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(expand=TRUE)
  
  
  ### Retrieve index of "rho" -------------------------------------------------------------------------------------------------------------------------
  index_rho <- NA
  for (a in 1:length(fit$all.hyper$random)) {
    if (fit$all.hyper$random[[a]]$group.hyper$theta$hyperid == 41001) {
      index_rho <- a
      break
    }
  }
  
  ### Retrieve prior mean and precision for rho
  mean_rho <- fit$all.hyper$random[[index_rho]]$group.hyper$theta$param[1]
  precision_rho <- fit$all.hyper$random[[index_rho]]$group.hyper$theta$param[2]
  
  ### Sample from normal distribution for rho
  prior_rho <- as.data.frame(rnorm(n=100000, mean=mean_rho, sd=sqrt(1/precision_rho)))
  colnames(prior_rho) <- c("Rho.Prior")
  
  ### Combine prior and posterior samples for Rho
  Rho <- as.data.frame(rbind(cbind(prior_rho[, 1], Distribution="Prior"), cbind(posteriors[, "GroupRho.for.space"], Distribution="Posterior")))
  colnames(Rho)[1] <- "Value"
  Rho$Value <- as.numeric(levels(Rho$Value))[Rho$Value]
  plot_rho <- ggplot() + theme_classic() + stat_density(data=Rho, aes(x=Value, group=Distribution, color=Distribution), position="identity", geom="line") + labs(title="Rho for space-time AR1 model")+ theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(expand=TRUE)
  
  ### Transform priors from theta-space to rho-space; from.theta: 2*exp(x)/(1+exp(x))-1
  Rho_mod <- as.data.table(Rho)
  Rho_mod[Distribution == "Posterior", "ModValue"] <- Rho_mod[Distribution == "Posterior", "Value"]
  Rho_mod[Distribution == "Prior", "ModValue"] <- 2*exp(Rho_mod[Distribution == "Prior", "Value"])/(1 + exp(Rho_mod[Distribution == "Prior", "Value"])) - 1
  plot_rho <- ggplot() + theme_classic() + stat_density(data=Rho_mod, aes(x=ModValue, group=Distribution, color=Distribution), position="identity", geom="line") +  labs(title=paste0(toupper(indicator), ' Rho for space-time AR1 model in ', region))  + theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(expand=TRUE)
  
  ### NOTE: While the shape of the prior as plotted above appears to be correct, the magnitudes of the peaks differ from those in INLA's plots. This needs to be investigated further.
  
  
  
  ### Retrieve index of "IID.ID" --------------------------------------------------------------------------------------------------------------------------------
  index_IID <- NA
  for (a in 1:length(fit$all.hyper$random)) {
    if (fit$all.hyper$random[[a]]$hyperid == "IID.ID") {
      index_IID <- a
      break
    }
  }
  
  ### Retrieve prior shape and inverse scale for IID.ID
  shape_IID <- fit$all.hyper$random[[index_IID]]$hyper$theta$param[1]
  inverse_scale_IID <- fit$all.hyper$random[[index_IID]]$hyper$theta$param[2]
  
  ### Derive loggamma density distribution for IID.ID
  prior.density <- as.data.table(cbind(seq(0, 20, 0.1), as.numeric(dlgamma(seq(0, 20, 0.1), shape_IID, inverse_scale_IID, log=FALSE)), "Prior"))
  
  ### Specify posterior density distribution for IID.ID
  post.density <- density(posteriors[, "Precision.for.IID.ID"])
  post.density <- as.data.table(cbind(as.numeric(post.density$x), as.numeric(post.density$y), "Posterior"))
  
  ### Combine prior and posterior densities for IID.ID
  IID <- as.data.frame(rbind(post.density, prior.density))
  colnames(IID) <- c("x", "y", "Distribution")
  IID$x <- as.numeric(IID$x)
  IID$y <- as.numeric(IID$y)
  
  ### Plot density distributions for IID.ID
  plot_iid <- ggplot() + theme_classic() + geom_line(data=IID, aes(x=x, y=y, group=Distribution, color=Distribution)) +  labs(title=paste0(toupper(indicator), ' precision for IID.ID (nugget) ', region)) + ylab("Density") + xlab("Value") + theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(expand=TRUE)
  
  
  
  ### Retrieve index of "CTRY.ID" -------------------------------------------------------------------------------------------------------------------------------------
  index_CTRY <- NA
  for (a in 1:length(fit$all.hyper$random)) {
    if (fit$all.hyper$random[[a]]$hyperid == "CTRY.ID") {
      index_CTRY <- a
      break
    }
  }
  
  ### Retrieve prior shape and inverse scale for CTRY.ID
  shape_CTRY <- fit$all.hyper$random[[index_CTRY]]$hyper$theta$param[1]
  inverse_scale_CTRY <- fit$all.hyper$random[[index_CTRY]]$hyper$theta$param[2]
  
  ### Derive loggamma density distribution for CTRY.ID
  prior.density <- as.data.table(cbind(seq(0, 20, 0.1), as.numeric(dlgamma(seq(0, 20, 0.1), shape_CTRY, inverse_scale_CTRY, log=FALSE)), "Prior"))
  
  ### Specify posterior density distribution for CTRY.ID
  post.density <- density(posteriors[, "Precision.for.CTRY.ID"])
  post.density <- as.data.table(cbind(as.numeric(post.density$x), as.numeric(post.density$y), "Posterior"))
  
  ### Combine prior and posterior densities for CTRY.ID
  CTRY <- as.data.frame(rbind(post.density, prior.density))
  colnames(CTRY) <- c("x", "y", "Distribution")
  CTRY$x <- as.numeric(CTRY$x)
  CTRY$y <- as.numeric(CTRY$y)
  
  ### Plot density distributions for CTRY.ID
  plot_ctry <- ggplot() + theme_classic() + geom_line(data=CTRY, aes(x=x, y=y, group=Distribution, color=Distribution)) +  labs(title=paste0(toupper(indicator), ' precision for CTRY.ID (country random effects) in ', region)) + ylab("Density") + xlab("Value") + theme(plot.title = element_text(hjust = 0.5)) + coord_cartesian(expand=TRUE)
  
  ### NOTE: While the shape of the IID.ID and CTRY.ID priors as plotted above appear to be correct, the magnitude differs a bit from those in INLA's plots. This needs to be investigated further.
  
  
  ### Get summary table ------------------------------------------------------------------------------------------------------
  fs <- summary(fit)
  hypers <- tableGrob(fs$hyperpar[1:6])
  
  ### Save plots
  savedir <- '<<<< FILEPATH REDACTED >>>>'
  pdf(paste0(savedir, 'posterior_prior_plots_', region, '.pdf'), width = 12, height = 7)
  grid.arrange(hypers, plot_th1, plot_th2,
               layout_matrix = rbind(c(1, 1),
                                     c(2, 3),
                                     c(2, 3)))
  grid.arrange(plot_rho, plot_iid, plot_ctry, 
               layout_matrix = rbind(c(1, 1),
                                     c(2, 3)))
  dev.off()
  
  ### End function
  return('Plots saved!')
}
