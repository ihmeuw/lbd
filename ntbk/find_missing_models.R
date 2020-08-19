rm(list = ls())
library(dplyr)
library(data.table)


indicators <- c('w_piped', 'w_network_cr', 'w_imp_cr',
  'w_unimp_cr', 's_piped', 's_imp_cr', 's_unimp_cr')

region_list <- 'chn + mng + khm + lao + mmr + tha + vnm + idn + phl + tls + png + ind + bgd + btn + lka + npl + pak + irq + jor + syr + afg + kgz + tjk + tkm + uzb + dza + egy + lby + mar + tun + ben + bfa + civ + cmr + gha + gin + gmb + gnb + lbr + mli + mrt + ner + nga + sen + sle + tcd + tgo + eth + sdn + som + ssd + yem + eri + ken + uga + tza + bdi + com + lso + mdg + moz + mwi + rwa + swz + zmb + zwe + bwa + nam + zaf + ago + caf + cog + gab + cod + bol + per + ecu + col + bra + guy + pry + mex + dom + gtm + hnd + hti + nic + slv + pan + cri'

region_list <- region_list       %>%
                  gsub(" ", "", .)          %>%
                  strsplit(., "+", fixed=T) %>%
                  unlist

miss <- list()
for (ii in indicators) {
	for (jj in region_list) {
		setwd("<<<< FILEPATH REDACTED >>>>")

		bricks <- list.files(pattern = '.grd')
		map <- bricks[grep(jj, bricks)]

		if (length(map) != 1) {
			miss[[1+length(miss)]] <- c(ii, jj)
		}
	}
}

miss_df <- as.data.frame(do.call(rbind, miss))
names(miss_df) <- c('indicator', 'iso')

arrange(count(miss_df, iso), -n)

filter(miss_df, indicator != 'w_network_cr')
