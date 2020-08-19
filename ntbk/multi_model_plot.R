library(dplyr)
library(ggplot2)
library(ggpubr)
library(raster)

indi <- c('w_piped', 'w_imp_cr', 's_piped', 's_imp_cr')
rundates <- paste0('2020_01_01_00_00_', 11:62)
regs <- c('dia_mid_east-afg')
cntry <- c('SYR')
ind_subnats <- c(23219, 129905, 227139, 228081, 303171, 133397, 129783, 59275, 174154,
234353, 93804, 129770, 163066, 165390, 238494, 409011, 285980, 349833,  349843, 414781, 398721,
414793)

geotbl <- read.csv("<<<< FILEPATH REDACTED >>>>", stringsAsFactors = FALSE)
a0_code <- unique(filter(geotbl, iso3 == cntry)$ADM0_CODE)

for (jj in indi) {
	filedir <- "<<<< FILEPATH REDACTED >>>>"
	setwd(filedir)
	stackagg <- readRDS(paste0('agg_stacker_admin_', regs, '_0.rds'))
	stack_dat <- data.frame(gam = sapply(1:19, function(w, x, y, z) {
		return(filter(y[w, x][[1]], geoid == z)$mean)
		},
		w = 1, y = stackagg, z = a0_code),
		gbm = sapply(1:19, function(w, x, y, z) {
		return(filter(y[w, x][[1]], geoid == z)$mean)
		},
		w = 1, y = stackagg, z = a0_code),
		lasso = sapply(1:19, function(w, x, y, z) {
		return(filter(y[w, x][[1]], geoid == z)$mean)
		},
		w = 1, y = stackagg, z = a0_code),
		year = 2000:2018)

	idat <- read.csv(paste0('input_data_bin0_', regs, '_0.csv'),
	 stringsAsFactors = FALSE)
	idat$indicator <- idat[[jj]]
	pidat <- idat %>%
		filter(!(nid %in% ind_subnats), country == cntry) %>%
		group_by(year, nid) %>%
		summarize(mean = weighted.mean(x = indicator/N, w = N*weight),
			N = sum(weight*N))

	plist <- list()
	plist[[1]] <- ggplot() +
		geom_line(data = stack_dat, aes(x = year, y = gam)) +
		ylim(0, 1) +
		xlim(2000, NA) +
		theme_bw() +
		theme(legend.position = 'none') +
		geom_point(data = pidat, aes(x = year, y = mean, size = log(N)),
			alpha = 0.75, col = 'red') +
		ggtitle('GAM')

	plist[[2]] <- ggplot() +
		geom_line(data = stack_dat, aes(x = year, y = gam)) +
		ylim(0, 1) +
		xlim(2000, NA) +
		theme_bw() +
		theme(legend.position = 'none') +
		geom_point(data = pidat, aes(x = year, y = mean, size = log(N)),
			alpha = 0.75, col = 'red') +
		ggtitle('GBM')

	plist[[3]] <- ggplot() +
		geom_line(data = stack_dat, aes(x = year, y = gam)) +
		ylim(0, 1) +
		xlim(2000, NA) +
		theme_bw() +
		geom_point(data = pidat, aes(x = year, y = mean, size = log(N)),
			alpha = 0.75, col = 'red') +
		theme(legend.position = 'none') +
		ggtitle('Lasso')

	png("<<<< FILEPATH REDACTED >>>>", width = 2560, height = 1440)
	print(
		ggarrange(plotlist = plist)
		)
	dev.off()

}
