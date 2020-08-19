rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)

setwd("<<<< FILEPATH REDACTED >>>>")

indi_group <- 'sani'

read_dat <- function(ig) {
	files <- c('a1_idat_agg', 'a1_mod_agg', 'stack_aggs')

	sets <- (lapply(files, function(x, group = ig) {
		as.data.frame(fread(paste0(group, '_', x, '.csv')))
	}))

}

datasets <- read_dat(indi_group)
oname <- c('idat', 'mod', 'stack')
for (i in 1:3) {
	assign(oname[i], datasets[[i]])
}

stack2 <- stack %>%
	filter(ADM0_NAME == 'India') %>%
	select(year, ADM1_NAME, ADM1_CODE, mod, prev = imp)
mod$mod <- 'inla'
all_mod <- bind_rows(mod, stack2)

idat <- filter(idat, !is.na(ADM1_NAME))

pdf('sani_vetting.pdf', width = 11, height = 8)

for (i in unique(all_mod$ADM1_NAME)) {
	print(
	ggplot() +
	geom_point(data = filter(idat, ADM1_NAME == i),
		aes(x = year, y = prev, size = ss, col = survey_series), alpha  =0.7) +
	geom_line(data = filter(all_mod, ADM1_NAME == i),
		aes(x = year, y = prev)) +
	facet_wrap(.~mod) +
	xlim(2000, 2017) +
	ylim(0, 1) +
	theme_bw() +
	ggtitle(paste0('Improved Sanitation ', i))
	)
}

dev.off()
