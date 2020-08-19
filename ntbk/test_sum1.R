library(dplyr)
region_list <- 'chn + mng + khm + lao + mmr + tha + vnm + idn + phl + tls + png + ind + bgd + btn + lka + npl + pak + irq + jor + syr + afg + kgz + tjk + tkm + uzb + dza + egy + lby + mar + tun + ben + bfa + civ + cmr + gha + gin + gmb + gnb + lbr + mli + mrt + ner + nga + sen + sle + tcd + tgo + eth + sdn + som + ssd + yem + eri + ken + uga + tza + bdi + com + lso + mdg + moz + mwi + rwa + swz + zmb + zwe + bwa + nam + zaf + ago + caf + cog + gab + cod + bol + per + ecu + col + bra + guy + pry + mex + dom + gtm + hnd + hti + nic + slv + pan + cri'
regions <- region_list       %>%
                  gsub(" ", "", .)          %>%
                  strsplit(., "+", fixed=T) %>%
                  unlist

regions <- 'chn'
draws <- round(runif(5, 3, 252))
results <- c()
i <- 0
for (dd in draws) {
	print("####################################")
	print("####################################")
	for (ii in regions) {
		# UNIMP
		setwd("<<<< FILEPATH REDACTED >>>>")
		load(list.files(pattern = "<<<< FILEPATH REDACTED >>>>"))
		a <- ((admin_2[,dd]))

		# PIPED
		setwd("<<<< FILEPATH REDACTED >>>>")
		load(list.files(pattern = "<<<< FILEPATH REDACTED >>>>"))
		b <- ((admin_2[,dd]))

		# IMP
		setwd("<<<< FILEPATH REDACTED >>>>")
		load(list.files(pattern = "<<<< FILEPATH REDACTED >>>>"))
		c <- ((admin_2[,dd]))

		# IMP OTHER
		setwd("<<<< FILEPATH REDACTED >>>>")
		load(list.files(pattern = "<<<< FILEPATH REDACTED >>>>"))
		d <- ((admin_2[,dd]))

		# OD
		setwd("<<<< FILEPATH REDACTED >>>>")
		load(list.files(pattern = "<<<< FILEPATH REDACTED >>>>"))
		e <- ((admin_2[,dd]))

		test1x <- max(a + b + d + e, na.rm = TRUE)
		test1n <- min(a + b + d + e, na.rm = TRUE)

		test2x <- max(a + c + e, na.rm = TRUE)
		test2n <- min(a + c + e, na.rm = TRUE)


		if (
			round(test1x, 3) != 1 | round(test2x, 3) != 1 |
			round(test1n, 3) != 1 | round(test2n, 3) != 1
			) {
			print(paste(ii, 'test1:', test1x, test1n,
			 'test2:', test2x, test2n))
			print(filter(sp_hierarchy_list,
			 ADM1_CODE %in% admin_1[which((a + c + e) > 1),'ADM1_CODE'])$ADM1_NAME)
		}
	}
	i <- i + 1
}

