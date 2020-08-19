setwd("<<<< FILEPATH REDACTED >>>>")

rm(list = ls())
library(dplyr)
library(mgcv)
library(ggplot2)
library(glmnet)
library(scam)
library(data.table)

indicators <- c('w_piped')

regions <- 'ind'

source('<<<< FILEPATH REDACTED >>>>/mbg/mbg_central/prep_functions.R')
source('<<<< FILEPATH REDACTED >>>>/mbg/mbg_central/shapefile_functions.R')
source('<<<< FILEPATH REDACTED >>>>/mbg/mbg_central/post_estimation_functions.R')
gbd <- get_gbd_locs('ind', TRUE, 'current')

for (indi in indicators) {
all_mod <- list()
	for (iso in regions) {
		message(paste(indi, iso))
		mydat <- read.csv("<<<< FILEPATH REDACTED >>>>",
		 stringsAsFactors = FALSE)
		mydat <- as.data.frame(mydat)
		mydat$indicator <- mydat[,indi]

		mydat <- mydat %>%
			group_by(nid, year, subnat, ADM1_CODE, ADM1_NAME) %>%
			filter(subnat == 0, !is.na(ADM1_CODE), year %in% c(2000:2020)) %>%
			summarize(mean = weighted.mean(x = indicator/N, w = N*weight),
				ss = sum(weight*N))

		for (a1 in unique(mydat$ADM1_CODE)) {
			mod_dat <- filter(mydat, ADM1_CODE == a1)

			if (a1 %in% setdiff(unique(mydat$ADM1_CODE), c())) {
				if (a1 %in% c(32105)) {
					gam_k <- 4
				} else {
					gam_k <- 6
				}
				spline_args = paste0("bs = 'mpi', k = ", as.numeric(gam_k))
				covariates <- 'year'
				f_bod = paste(paste0('s(',covariates,', ',spline_args,')'),collapse = " + ")
				gam_formula = paste0("mean ~ 1 + ", f_bod)
				gam_formula = as.formula(gam_formula)

				gam = try(scam(gam_formula, data = mod_dat, family = 'binomial', weights = mod_dat[['ss']]))

				if (class(gam) == 'try-error') {
					results <- data.frame(year = 2000:2020,
						model = NA, type = 'response')
				} else {
					results <- data.frame(year = 2000:2020,
					model = predict(gam, newdata = data.frame(year = 2000:2020), type = 'response'))
				}
				results$mean <- results$model

			} else {
				if (a1 %in% c()) {
					# linear model
					glm <- glm(mean ~ 1 + year, data = mod_dat, weights = mod_dat$ss, family = binomial(link = 'logit'))
					results <- data.frame(year = 2000:2020,
						model = predict(glm, newdata = data.frame(year = 2000:2020), type = 'response'))
					results$mean <- results$model

				} else {
					# gam
					gam_k <- 3

					if (a1 %in% c()) {
						gam_k <- 3
					}

					if (a1 %in% c()) {
						gam_k <- 4
					}

					spline_args = paste0("bs = 'ts', k = ", as.numeric(gam_k))
					covariates <- 'year'
					f_bod = paste(paste0('s(',covariates,', ',spline_args,')'),collapse = " + ")
					gam_formula = paste0("mean ~ 1 + ", f_bod)
					gam_formula = as.formula(gam_formula)
					gam = mgcv::gam(gam_formula, data = mod_dat, family = 'binomial', weights = mod_dat[['ss']])
					results <- data.frame(year = 2000:2020,
						model = predict(gam, newdata = data.frame(year = 2000:2020), type = 'response'))
					results$mean <- results$model
				}
			}

			results <- dplyr::select(results, year, mean)
			results$indi <- indi
			results$reg <- a1
			results$ADM1_NAME <- unique(filter(mydat, ADM1_CODE == a1)$ADM1_NAME)
			all_mod[[1+length(all_mod)]] <- results
		}
	}
	all_mod <- do.call(rbind, all_mod)
	all_mod$ADM1_CODE <- all_mod$reg
	all_mod <- left_join(all_mod, dplyr::select(gbd, ADM1_CODE, location_id))
	all_mod$name <- all_mod$location_id
	write.csv(all_mod, "<<<< FILEPATH REDACTED >>>>", row.names = FALSE)
}

mydat$reg <- mydat$ADM1_CODE
pdf("<<<< FILEPATH REDACTED >>>>",
	height = 8.5, width = 11)
print(
ggplot() +
	geom_point(data = mydat, aes(x = year, y = mean, size = ss)) +
	geom_line(data = all_mod, aes(x = year, y = mean)) +
	xlim(2000, 2020) + ylim(0, 1) +
	facet_wrap(. ~ reg)

)
dev.off()

pdf("<<<< FILEPATH REDACTED >>>>",
	height = 8.5, width = 11)
print(
ggplot() +
	geom_point(data = mydat, aes(x = year, y = mean, size = ss)) +
	geom_line(data = all_mod, aes(x = year, y = mean)) +
	xlim(2000, 2020) + ylim(0, 1) +
	facet_wrap(. ~ ADM1_NAME)

)
dev.off()
