library(dplyr)
library(ggplot2)
library(ggrepel)

repo <- "<<<< FILEPATH REDACTED >>>>"
reg_list <- readRDS(paste0(repo, 'wash/00_reg_list.rds'))
reg_list[['dia_sessa']] <- setdiff(reg_list[['dia_sessa']], 'SYC')
reg_list[['dia_ncaca']] <- setdiff(reg_list[['dia_ncaca']], 'DMA')
reg_list[['dia_ncaca']] <- setdiff(reg_list[['dia_ncaca']], c('GRD', 'VCT', 'LCA'))

regs <- c('dia_chn_mng', 'dia_se_asia', 'idn', 'phl', 'tls', 'png', 'ind', 'dia_south_asia-ind',
    'dia_mid_east-afg', 'afg', 'dia_central_asia', 'dia_name', 'dia_wssa',
    'dia_afr_horn-dji-eri', 'dji_eri',
    'dia_nessa', 'dia_sessa',
    'dia_sssa', 'dia_cssa-cod', 'cod',
    'dia_sam_andes', 'dia_sam_trop',
    'mex', 'dia_ncaca', 'dia_scaca')


mydat <- read.csv("<<<< FILEPATH REDACTED >>>>",
	stringsAsFactors = FALSE)

pdf("<<<< FILEPATH REDACTED >>>>", width = 11, height = 11)
for (i in regs) {
	message(i)
	for (j in reg_list[[i]]) {
		stack_dat <- list()

		for (k in paste0('2020_01_01_00_00_6', 3:5)) {

			file.dir <- "<<<< FILEPATH REDACTED >>>>"
			setwd(file.dir)
			a0_stack <- list.files(pattern = 'agg_stacker')[grep('0.rds', list.files(pattern = 'agg_stacker'))]

			stackfile <- a0_stack[grep(paste0('_', i), a0_stack)]
			if (length(stackfile) > 0) {
			    result <- do.call(rbind, readRDS(stackfile))

				vars <- c('mean', 'iso3', 'stacker',
			    	names(result)[grep('ADM', names(result))])
				result <- result[,vars]

				result$year <- sapply(result$stacker, function(x) {
			    	       as.numeric(strsplit(x, '.', fixed = TRUE)[[1]][2]) + 1999
			        	 })
				result$stacker <- sapply(result$stacker, function(x) {
			              strsplit(x, '.', fixed = TRUE)[[1]][1]
			            })
				result <- filter(result, iso3 == j)
				result$rd <- k

				stack_dat[[1+length(stack_dat)]] <- result

			}
		}

		stack_dat <- do.call(rbind, stack_dat)

		if (!is.null(stack_dat)) {
			if (nrow(stack_dat) > 0) {
				subset <- filter(mydat, country == j) %>%
					group_by(year, nid, subnat) %>%
					summarize(prev = weighted.mean(x = s_imp_cr/N, w = weight*N),
						ss = sum(weight*N))
				print(
					ggplot(stack_dat) +
						geom_line(aes(x = year, y = mean)) +
						geom_point(data = subset, aes(x = year, y = prev, col = as.factor(subnat),
							size = log(ss))) +
						geom_text_repel(data = subset, aes(x = year, y = prev, label = nid)) +
						facet_wrap(rd ~ stacker) +
						xlim(2000, 2017) + ylim(0, 1) + ggtitle(j)
				)
			}
		}
	}
}
dev.off()
