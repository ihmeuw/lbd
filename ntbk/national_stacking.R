setwd("<<<< FILEPATH REDACTED >>>>")

library(data.table)
library(dplyr)
library(mgcv)
library(ggplot2)
library(glmnet)
indi <- 'w_piped'
iso <- 'afg'

cov <- read.csv(paste0(iso, '.csv'), stringsAsFactors = FALSE)
mydat <- read.csv("<<<< FILEPATH REDACTED >>>>", stringsAsFactors = FALSE)
mydat <- as.data.frame(mydat)
mydat$indicator <- mydat[,indi]
mydat <- filter(mydat, country == toupper(iso)) %>%
	group_by(nid, year, subnat) %>%
	filter(subnat == 0) %>%
	summarize(mean = weighted.mean(x = indicator/N, w = N*weight),
		ss = sum(weight*N))

mod_dat <- left_join(mydat, cov)

nids <- unique(mod_dat$nid)
if (length(nids) > 2) {
	if (length(nids) <= 5) {
		mod_dat$fold <- 1:nrow(mod_dat)
	} else {
		mod_dat$fold <- sample(1:5, nrow(mod_dat), replace = TRUE)
	}

} else {
	mod_dat$fold <- 1
}

# gam
if (length(nids) > 3) {

	gam_k <- if (length(nids) >= 5) {4} else {2}
	spline_args = paste0("bs = 'ts', k = ", as.numeric(gam_k))
	covariates <- 'year'
	f_bod = paste(paste0('s(',covariates,', ',spline_args,')'),collapse = " + ")
	gam_formula = paste0("mean ~ 1 + ", f_bod)
	gam_formula = as.formula(gam_formula)
	gam = mgcv::gam(gam_formula, data = mod_dat, family = 'binomial', weights = mod_dat[['ss']])
	a0_mod_gam <- data.frame(year = 2000:2017,
		gam = predict(gam, newdata = data.frame(year = 2000:2017), type = 'response'))

	results <- list()
	for (ii in unique(mod_dat$fold)) {
		pred_dat <- filter(mod_dat, fold == ii)
		fit_dat <- filter(mod_dat, fold != ii)
		gam = mgcv::gam(gam_formula, data = fit_dat,
			family = 'binomial', weights = fit_dat[['ss']])
		pred_dat$gam <- predict(gam, newdata = pred_dat, type = 'response')
		pred_dat <- as.data.frame(pred_dat)
		results[[1+length(results)]] <- as.data.frame(pred_dat)
	}

	if (length(results) == 1) {
		gam_fit_df <- results[[1]]
	} else {gam_fit_df <- data.frame(do.call(rbind, results))}


	# glm
	glm <- glm(mean ~ 1 + year, data = mod_dat, weights = mod_dat$ss, family = binomial(link = 'logit'))
	a0_mod_glm <- data.frame(year = 2000:2017,
		glm = predict(glm, newdata = data.frame(year = 2000:2017), type = 'response'))

	results <- list()
	for (ii in unique(mod_dat$fold)) {
		pred_dat <- filter(mod_dat, fold == ii)
		fit_dat <- filter(mod_dat, fold != ii)
		glm = glm(mean ~ 1 + year, data = fit_dat, weights = fit_dat$ss,
			family = binomial(link = 'logit'))
		 pred_dat$glm <- predict(glm, newdata = pred_dat, type = 'response')
		 pred_dat <- as.data.frame(pred_dat)
		 results[[1+length(results)]] <- pred_dat
	}

	if (length(results) == 1) {
		glm_fit_df <- results[[1]]
	} else {glm_fit_df <- data.frame(do.call(rbind, results))}

	ensemble_cov <- left_join(glm_fit_df, gam_fit_df)
	ensemble_dat <- left_join(mod_dat, ensemble_cov)
	ensemble <- glm(mean ~ 1 + gam + glm, data = ensemble_dat, weights = ensemble_dat$ss, family = binomial(link = 'logit'))
	a0_mod_ensemble <- data.frame(year = 2000:2017,
		ensemble = predict(ensemble, newdata = left_join(a0_mod_glm, a0_mod_gam), type = 'response'))

	results <- left_join(left_join(a0_mod_glm, a0_mod_gam), a0_mod_ensemble)
	results$mean <- results$ensemble

} else {
	glm <- glm(mean ~ 1 + year, data = mod_dat, weights = mod_dat$ss, family = binomial(link = 'logit'))
	results <- data.frame(year = 2000:2017,
		glm = predict(glm, newdata = data.frame(year = 2000:2017), type = 'response'))
	results$mean <- results$glm
}

ggplot() +
	geom_point(data = mydat, aes(x = year, y = mean, size = ss)) +
	geom_line(data = results, aes(x = year, y = mean)) +
	xlim(2000, 2017) + ylim(0, 1)
