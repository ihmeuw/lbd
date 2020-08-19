setwd("<<<< FILEPATH REDACTED >>>>")

library(data.table)
rm(list = ls())
library(dplyr)
library(mgcv)
library(ggplot2)
library(glmnet)

indicators <- c('s_piped')

region_list <- 'chn + mng + khm + lao + mmr + tha + vnm + idn + phl + tls + png + ind + bgd + btn + lka + npl + pak + irq + jor + syr + afg + kgz + tjk + tkm + uzb + dza + egy + lby + mar + tun + ben + bfa + civ + cmr + gha + gin + gmb + gnb + lbr + mli + mrt + ner + nga + sen + sle + tcd + tgo + eth + sdn + som + ssd + yem + eri + ken + uga + tza + bdi + com + lso + mdg + moz + mwi + rwa + swz + zmb + zwe + bwa + nam + zaf + ago + caf + cog + gab + cod + bol + per + ecu + col + bra + guy + pry + mex + dom + gtm + hnd + hti + nic + slv + pan + cri'
regions <- region_list       %>%
                  gsub(" ", "", .)          %>%
                  strsplit(., "+", fixed=T) %>%
                  unlist

lin_rake <- read.csv("<<<< FILEPATH REDACTED >>>>",
	stringsAsFactors = FALSE)
for (indi in indicators) {
all_mod <- list()
	for (iso in regions) {
		message(paste(indi, iso))
		cov <- read.csv(paste0(iso, '.csv'), stringsAsFactors = FALSE)
		mydat <- read.csv("<<<< FILEPATH REDACTED >>>>", stringsAsFactors = FALSE)
		mydat <- as.data.frame(mydat)
		mydat$indicator <- mydat[,indi]
		mydat <- filter(mydat, country == toupper(iso)) %>%
			group_by(nid, year, subnat) %>%
			filter(subnat == 0, year %in% c(2000:2017)) %>%
			summarize(mean = weighted.mean(x = indicator/N, w = N*weight),
				ss = sum(weight*N))

		mod_dat <- left_join(mydat, cov)

		nids <- unique(mod_dat$nid)

		if (iso %in% c(lin_rake$iso, 'uzb')) {
			if (indi %in% c(filter(lin_rake, iso == iso)$indicator, 's_piped')) {
				if (length(nids) < 2) {
					nids <- nids} else {nids <- rep(2, 2)}
			}
		}

		# gam
		if (length(nids) > 3) {

			gam_k <- if (length(nids) >= 5) {
				if (iso == 'ind' & indi %in% c('s_piped', 's_imp_cr', 's_unimp_cr')) {
					6
				} else {4}
				} else {2}
			spline_args = paste0("bs = 'ts', k = ", as.numeric(gam_k))
			covariates <- 'year'
			f_bod = paste(paste0('s(',covariates,', ',spline_args,')'),collapse = " + ")
			gam_formula = paste0("mean ~ 1 + ", f_bod)
			gam_formula = as.formula(gam_formula)
			gam = mgcv::gam(gam_formula, data = mod_dat, family = 'binomial', weights = mod_dat[['ss']])
			results <- data.frame(year = 2000:2017,
				model = predict(gam, newdata = data.frame(year = 2000:2017), type = 'response'))
			results$mean <- results$model

		} else {
			if (length(nids) >= 2) {
				glm <- glm(mean ~ 1 + year, data = mod_dat, weights = mod_dat$ss, family = binomial(link = 'logit'))
				results <- data.frame(year = 2000:2017,
					model = predict(glm, newdata = data.frame(year = 2000:2017), type = 'response'))
				results$mean <- results$model
			} else {
				glm <- glm(mean ~ 1, data = mod_dat, weights = mod_dat$ss, family = binomial(link = 'logit'))
				results <- data.frame(year = 2000:2017,
					model = predict(glm, newdata = data.frame(year = 2000:2017), type = 'response'))
				results$mean <- results$model
			}
		}

		results <- select(results, year, mean)
		results$indi <- indi
		results$reg <- iso
		all_mod[[1+length(all_mod)]] <- results
	}
	all_mod <- do.call(rbind, all_mod)
	write.csv(all_mod, "<<<< FILEPATH REDACTED >>>>", row.names = FALSE)
}
