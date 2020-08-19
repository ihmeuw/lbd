library(dplyr)


 indis <- c('s_piped', 's_imp_cr', 's_unimp_cr', 's_network_cr',
  	's_imp', 's_imp_other', 's_unimp', 's_od', 's_network',
  'w_piped', 'w_imp_cr', 'w_unimp_cr', 'w_network_cr', 'w_imp',
  'w_imp_other', 'w_unimp', 'w_surface', 'w_network')

check_cond <- 'mod'
	indis <- c('s_piped', 's_imp', 's_unimp', 's_od', 's_imp_other',
		  'w_piped', 'w_imp', 'w_unimp', 'w_imp_other', 'w_surface')

region_list <- 'chn + mng + khm + lao + mmr + tha + vnm + idn + phl + tls + png + ind + bgd + btn + lka + npl + pak + irq + jor + syr + afg + kgz + tjk + tkm + uzb + dza + egy + lby + mar + tun + ben + bfa + civ + cmr + gha + gin + gmb + gnb + lbr + mli + mrt + ner + nga + sen + sle + tcd + tgo + eth + sdn + som + ssd + yem + eri + ken + uga + tza + bdi + com + lso + mdg + moz + mwi + rwa + swz + zmb + zwe + bwa + nam + zaf + ago + caf + cog + gab + cod + bol + per + ecu + col + bra + guy + pry + mex + dom + gtm + hnd + hti + nic + slv + pan + cri'
regions <- region_list       %>%
                  gsub(" ", "", .)          %>%
                  strsplit(., "+", fixed=T) %>%
                  unlist

results <- data.frame(matrix(0, nrow = length(regions), ncol = length(indis)))
names(results) <- indis
results$regions <- regions


for (ii in indis) {
	for (rr in regions) {
		setwd("<<<< FILEPATH REDACTED >>>>")

		file <- list.files(pattern = paste0('unraked_admin_draws_eb_bin0_', rr))
		if (length(file) == 1) {
			results[which(results$regions == rr),
			ii] <- 1
		}
	}
}

results$total <- apply(select(results, -regions), 1, sum)
miss_iso <- as.character(select(filter(results, total != length(indis)), regions)$regions)
miss_df <- filter(results, regions %in% miss_iso)

setdiff(miss_df$regions, mod_miss)

mod_miss <- miss_df$regions

