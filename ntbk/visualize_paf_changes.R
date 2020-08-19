rm(list = ls())
# Load libraries
libs <- c('dplyr', 'purrr', 'tidyr', 'data.table', 'ggplot2')
lapply(libs, library, character.only = TRUE)

# function to reformat draw obj into long data.frame
reformat_agg_draws <- function(indicator, run_date, region, ad_lvl) {
	setwd("<<<< FILEPATH REDACTED >>>>")
	########################################################
	all_files <- list.files(pattern = '_unraked_admin_draws_')
	########################################################
	load(all_files[grep(paste0('_', region, '_0'), all_files)])
	x <- get(paste0('admin_', ad_lvl))
	x <- dplyr::select(x, -pop)
	x <- mutate(gather(x, 'draw', 'prev', 3:252), agg_level = ad_lvl)
	names(x)[grep('ADM', names(x))] <- 'code'
	########################################################
	return(x)
}

# Read in a indicator-regions data.frame and make a master agg_draw frame
get_draws <- function(ind, rd, reg) {
	message(paste("Reformating all admin levels via reformat_agg_draws for", ind, reg))
	agg_draw <- do.call(bind_rows, lapply(0:2,
		reformat_agg_draws, indicator = ind, run_date = rd, region = reg))
	return(agg_draw)
}

# function to load multiple indicators across a specified region list & run date
load_agg_levels <- function(rd, levels, reg_list) {
	results <- list()
	for (ii in levels) {
		message("Loading all regions via get_draws")
		xx <- do.call(bind_rows, lapply(reg_list,
			get_draws, ind = ii, rd = rd))
		names(xx)[grep('prev', names(xx))] <- ii
		results[[length(results)+1]] <- data.table(xx)
	}
	return(reduce(results, merge))
}

# function to load sp hierarchy list
get_sp_tbl <- function(ind, rd, reg) {
	setwd("<<<< FILEPATH REDACTED >>>>")
	all_files <- list.files(pattern = '_unraked_admin_draws_')
	load(all_files[grep(paste0('_', reg, '_0'), all_files)])
	sp_tbl <- as.data.frame(sp_hierarchy_list)
	ad_codes <- grep('CODE', names(sp_tbl))
	ad_name <- grep('NAME', names(sp_tbl))
	for (cc in ad_codes){
		sp_tbl[,cc] <- as.numeric(as.character(sp_tbl[,cc]))
	}
	for (cc in ad_name){
		sp_tbl[,cc] <- as.character(sp_tbl[,cc])
	}
	return(sp_tbl)
}

# Load Data
levels <- c('w_piped', 'w_imp_other', 'w_unimp', 'w_surface')
indi_group <- 'water'
## WASH has to use two run dates since regions are split into two model runs
run_date <- '2019_09_29_00_00_01'
regions <- c('dia_wssa-civ')
all_draws <- load_agg_levels(run_date, levels, regions)
sp_tbl <- do.call(bind_rows, lapply(regions, get_sp_tbl, ind = levels[1], rd = run_date))

# WASH specific PAF math
if (indi_group  == 'water') {
	intmd <- mutate(all_draws,
	piped = w_piped,
	imp = w_imp_other,
	unimp = w_unimp + w_surface
	) %>%
	select(year, code, draw, agg_level, piped, imp, unimp)

	# Add ADM0_NAME to draw obj
	intmd_0 <- filter(intmd, agg_level == 0) %>%
		mutate(ADM0_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE)))
	intmd_1 <- filter(intmd, agg_level == 1) %>%
		mutate(ADM1_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE, ADM1_CODE))) %>%
		select(-ADM1_CODE)
	intmd_2 <- filter(intmd, agg_level == 2) %>%
		mutate(ADM2_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE, ADM2_CODE))) %>%
		select(-ADM2_CODE)
	intmd <- bind_rows(intmd_0, intmd_1, intmd_2)


	# Merge on HWT information
	w_treat <- read.csv("<<<< FILEPATH REDACTED >>>>")
	intmd <- left_join(intmd,
		distinct(dplyr::select(w_treat, year, ADM0_CODE, no_treat, filtered, boil))) %>%
		mutate(p_f = piped*filtered,
			p_b = piped*boil,
			p_n = piped*no_treat,
			i_f = imp*filtered,
			i_b = imp*boil,
			i_n = imp*no_treat,
			u_f = unimp*filtered,
			u_b = unimp*boil,
			u_n = unimp*no_treat)

	# Calculate PAF
	intmd <- mutate(intmd,
		rrsum = p_b*1 + p_f*1.652 + p_n*2.401 +
		        i_b*1.118 + i_f*1.848 + i_n*2.685 +
		        u_b*1.364 + u_f*2.254 + u_n*3.276) %>%
		mutate(paf = (rrsum - 1)/rrsum)

	intmd$draw <- gsub('V', '', intmd$draw)

	cntry_dat <- intmd %>%
		filter(agg_level == 0, ADM0_NAME == 'Nigeria') %>%
		group_by(year, code) %>%
		summarize(
			piped = mean(piped),
			imp = mean(imp),
			unimp = mean(unimp),
			no_treat = mean(no_treat),
			filtered = mean(filtered),
			boil = mean(boil),
			p_f = mean(p_f),
			p_b = mean(p_b),
			p_n = mean(p_n),
			i_f = mean(i_f),
			i_b = mean(i_b),
			i_n = mean(i_n),
			u_f = mean(u_f),
			u_b = mean(u_b),
			u_n = mean(u_n),
			paf = mean(paf)
			) %>%
		gather("level", "prev",
			piped, imp, unimp, no_treat, filtered, boil, p_f, p_b, p_n, i_f, i_b, i_n,
			u_f, u_b, u_n)
}

if (indi_group == 'sani') {
	# Calculate PAF
	intmd <- mutate(all_draws,
		piped = s_piped,
		imp = s_imp_other,
		unimp = s_unimp + s_od
		) %>%
		mutate(rrsum = (piped*1) + (imp*2.595) + (unimp*3.242),
			paf = (rrsum -1)/rrsum)

	# Add ADM0_NAME to draw obj
	intmd_0 <- filter(intmd, agg_level == 0) %>%
		mutate(ADM0_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE)))
	intmd_1 <- filter(intmd, agg_level == 1) %>%
		mutate(ADM1_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE, ADM1_CODE))) %>%
		select(-ADM1_CODE)
	intmd_2 <- filter(intmd, agg_level == 2) %>%
		mutate(ADM2_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE, ADM2_CODE))) %>%
		select(-ADM2_CODE)
	intmd <- bind_rows(intmd_0, intmd_1, intmd_2)

	intmd$draw <- gsub('V', '', intmd$draw)
	cntry_dat <- intmd %>%
		filter(agg_level == 0, ADM0_NAME == 'Nigeria') %>%
		group_by(year, code) %>%
		summarize(
			piped = mean(piped),
			imp = mean(imp),
			unimp = mean(unimp),
			paf = mean(paf)
			) %>%
		gather("level", "prev", piped, imp, unimp)
}

# Add mortality to country data
dia_mort <- read.csv("<<<< FILEPATH REDACTED >>>>")

dia_mort2 <- filter(dia_mort, ADM0_NAME == 'Nigeria', year %in% c(2000, 2017)) %>%
	rename(code = ADM2_CODE, dia_deaths = mean) %>%
	select(code, year, dia_deaths)
cntry_dat <- left_join(cntry_dat, dia_mort2) %>%
	mutate(paf_deaths = dia_deaths * paf)
# Trim to variables of interest & format variables
if (indi_group == 'sani') {
	setwd("<<<< FILEPATH REDACTED >>>>")
	pdf('Nigeria_a2_sani_dx.pdf')
		cntry_dat$level_f <- factor(cntry_dat$level, levels = c('piped', 'imp', 'unimp'))
		print(
			ggplot(cntry_dat) +
				geom_histogram(aes(prev, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 0.01) +
				facet_grid(level_f ~ .) +
				xlim(0, 1) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Levels for Sanitation')
		)

		cntry_dat3 <- cntry_dat %>%
			ungroup() %>%
			select(year, code, paf, dia_deaths, paf_deaths) %>%
			distinct()

		print(
			ggplot(cntry_dat3) +
				geom_histogram(aes(paf, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 0.01) +
				xlim(0, 1) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 PAF for Sanitation')
		)

		print(
			ggplot(cntry_dat3) +
				geom_histogram(aes(dia_deaths, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 10) +
				xlim(0, NA) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Diarrhea Deaths')
		)

		print(
			ggplot(cntry_dat3) +
				geom_histogram(aes(paf_deaths, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 10) +
				xlim(0, NA) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Diarrhea Deaths for Sanitation')
		)

		cntry_dat4 <- cntry_dat3 %>%
			mutate(year = paste0('y', year)) %>%
			select(year, code, paf_deaths) %>%
			spread(year, paf_deaths)

		print(
			ggplot(cntry_dat4) +
				geom_point(aes(x = y2000, y = y2017)) +
				geom_abline(slope = 1, sintercept = 0, col = 'red') +
				xlim(0, 2000) + ylim(0, 2000) +
				theme_bw() +
				ggtitle('Nigeria: Change in A2 Diarrhea Deaths for Sanitation')
		)

	dev.off()
}

if (indi_group == 'water') {
	setwd("<<<< FILEPATH REDACTED >>>>")
	pdf('Nigeria_a2_water_dx.pdf')

		sub_lvl <- c('piped', 'imp', 'unimp')
		cntry_dat2 <- filter(cntry_dat, level %in% sub_lvl)
		cntry_dat2$level_f <- factor(cntry_dat2$level, levels = sub_lvl)
		print(
			ggplot(cntry_dat2) +
				geom_histogram(aes(prev, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 0.01) +
				facet_grid(level_f ~ .) +
				xlim(0, 1) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Levels for Water')
		)

		sub_lvl <- c('boil', 'filtered', 'no_treat')
		cntry_dat2 <- filter(cntry_dat, level %in% sub_lvl)
		cntry_dat2$level_f <- factor(cntry_dat2$level, levels = sub_lvl)
		print(
			ggplot(cntry_dat2) +
				geom_histogram(aes(prev, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 0.01) +
				facet_grid(level_f ~ .) +
				xlim(0, 1) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Levels for Water Treatment')
		)

		sub_lvl <- c('p_b', 'i_b', 'u_b', 'p_f', 'i_f', 'u_f', 'p_n', 'i_n', 'u_n')
		cntry_dat2 <- filter(cntry_dat, level %in% sub_lvl)
		cntry_dat2$level_f <- factor(cntry_dat2$level, levels = sub_lvl)
		print(
			ggplot(cntry_dat2) +
				geom_histogram(aes(prev, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 0.01) +
				facet_grid(level_f ~ .) +
				xlim(0, 1) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Levels for Water Exposures')
		)


		cntry_dat3 <- cntry_dat %>%
			ungroup() %>%
			select(year, code, paf, dia_deaths, paf_deaths) %>%
			distinct()

		print(
			ggplot(cntry_dat3) +
				geom_histogram(aes(paf, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 0.01) +
				xlim(0, 1) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 PAF for Water')
		)

		print(
			ggplot(cntry_dat3) +
				geom_histogram(aes(dia_deaths, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 10) +
				xlim(0, NA) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Diarrhea Deaths')
		)

		print(
			ggplot(cntry_dat3) +
				geom_histogram(aes(paf_deaths, col = as.factor(year), fill = as.factor(year)), position = 'identity', alpha = 0.3, binwidth = 10) +
				xlim(0, NA) +
				theme_bw() +
				ggtitle('Nigeria: Distribution of A2 Diarrhea Deaths for Water')
		)

		cntry_dat4 <- cntry_dat3 %>%
			mutate(year = paste0('y', year)) %>%
			select(year, code, paf_deaths) %>%
			spread(year, paf_deaths)

		print(
			ggplot(cntry_dat4) +
				geom_point(aes(x = y2000, y = y2017)) +
				geom_abline(slope = 1, sintercept = 0, col = 'red') +
				xlim(0, 2000) + ylim(0, 2000) +
				theme_bw() +
				ggtitle('Nigeria: Change in A2 Diarrhea Deaths for Water')
		)
	dev.off()
}
