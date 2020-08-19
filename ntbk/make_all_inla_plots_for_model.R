rm(list = ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

indicators <- c('s_piped', 's_imp_cr', 's_unimp_cr',
	'w_piped', 'w_imp_cr', 'w_unimp_cr')

region_list <- 'chn + mng + khm + lao + mmr + tha + vnm + idn + phl + tls + png + ind + bgd + btn + lka + npl + pak + irq + jor + syr + afg + kgz + tjk + tkm + uzb + dza + egy + lby + mar + tun + ben + bfa + civ + cmr + gha + gin + gmb + gnb + lbr + mli + mrt + ner + nga + sen + sle + tcd + tgo + eth + sdn + som + ssd + yem + eri + ken + uga + tza + bdi + com + lso + mdg + moz + mwi + rwa + swz + zmb + zwe + bwa + nam + zaf + ago + caf + cog + gab + cod + bol + per + ecu + col + bra + guy + pry + mex + dom + gtm + hnd + hti + nic + slv + pan + cri'
regions <- region_list       %>%
                  gsub(" ", "", .)          %>%
                  strsplit(., "+", fixed=T) %>%
                  unlist

region_list <- regions
pdf("<<<< FILEPATH REDACTED >>>>", width = 11, height = 8.5)
for (iso  in region_list) {
	plots <- list()
	for (indi in indicators) {
		message(paste(indi, iso))
		mydat <- read.csv("<<<< FILEPATH REDACTED >>>>",
			stringsAsFactors = FALSE)
		mydat <- as.data.frame(mydat)
		mydat$indicator <- mydat[,indi]
		mydat <- filter(mydat, country == toupper(iso), year %in% 2000:2017) %>%
			group_by(nid, year, subnat) %>%
			summarize(mean = weighted.mean(x = indicator/N, w = N*weight),
				ss = sum(weight*N))
		if (length(mydat$subnat) == 0) {mydat$subnat <- 0}
		raking <- read.csv("<<<< FILEPATH REDACTED >>>>",
			stringsAsFactors = FALSE)
		raking <- filter(raking, reg == iso)

		setwd("<<<< FILEPATH REDACTED >>>>")
		files <- list.files(pattern = '_raked_admin_draws')
		ff <- files[grep(iso, files)]

		if (length(ff) != 0) {
			load(ff)
			draws <- as.data.frame(admin_0)[,grep('V', names(admin_0))]
			mean <- apply(draws, 1, mean)
			if (!is.na(mean)) {
				lci <- apply(draws, 1, quantile, probs = 0.025)
				uci <- apply(draws, 1, quantile, probs = 0.975)
				proc <- data.frame(year = admin_0$year, ADM0_CODE = admin_0$ADM0_CODE,
					mean = mean, lci = lci, uci = uci)

				plots[[1+length(plots)]] <- ggplot() +
					geom_point(data = mydat, aes(x = year, y = mean, size = ss, col = as.factor(subnat))) +
					geom_line(data = proc, aes(x = year, y = mean)) +
					geom_ribbon(data = proc, aes(x = year, ymin = lci, ymax = uci), alpha = 0.3) +
					geom_line(data = raking, aes(x = year, y = mean), col = 'blue') +
					xlim(2000, 2017) + ylim(0, 1) +
					theme(legend.position = 'none') +
					ggtitle(paste(iso, indi))
			}
		}

	}
	if (length(plots) > 0) {
		print(ggarrange(plotlist = plots, ncol = 3, nrow = 2))
	}
}

dev.off()
