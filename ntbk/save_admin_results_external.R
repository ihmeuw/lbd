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
levels <- c('s_piped', 's_imp_other', 's_unimp', 's_od')
indi_group <- 'sani'
## WASH has to use two run dates since regions are split into two model runs
run_date <- '2019_09_29_00_00_01'
regions <- c('dia_wssa-civ')
all_draws <- load_agg_levels(run_date, levels, regions)
sp_tbl <- do.call(bind_rows, lapply(regions, get_sp_tbl, ind = levels[1], rd = run_date))

# Add ADM0_NAME to draw obj
intmd_0 <- filter(all_draws, agg_level == 0) %>%
	mutate(ADM0_CODE = code) %>%
	left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE)))
intmd_1 <- filter(all_draws, agg_level == 1) %>%
	mutate(ADM1_CODE = code) %>%
	left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME)))
intmd_2 <- filter(all_draws, agg_level == 2) %>%
	mutate(ADM2_CODE = code) %>%
	left_join(distinct(select(sp_tbl, ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME,
		ADM2_CODE, ADM2_NAME)))

bmgf_0 <- filter(intmd_0, ADM0_NAME == 'Nigeria') %>%
	mutate(piped = w_piped,
		imp = w_piped + w_imp_other) %>%
	group_by(ADM0_NAME, ADM0_CODE, year) %>%
	summarize(
		piped_mean = mean(piped),
		piped_lci = quantile(piped, probs = 0.025, na.rm = TRUE),
		piped_uci = quantile(piped, probs = 0.975, na.rm = TRUE),
		imp_mean = mean(imp),
		imp_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		imp_uci = quantile(imp, probs = 0.975, na.rm = TRUE))

bmgf_1 <- filter(intmd_1, ADM0_NAME == 'Nigeria') %>%
	mutate(piped = w_piped,
		imp = w_piped + w_imp_other) %>%
	group_by(ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME, year) %>%
	summarize(
		piped_mean = mean(piped),
		piped_lci = quantile(piped, probs = 0.025, na.rm = TRUE),
		piped_uci = quantile(piped, probs = 0.975, na.rm = TRUE),
		imp_mean = mean(imp),
		imp_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		imp_uci = quantile(imp, probs = 0.975, na.rm = TRUE))

bmgf_2 <- filter(intmd_2, ADM0_NAME == 'Nigeria') %>%
	mutate(piped = w_piped,
		imp = w_piped + w_imp_other) %>%
	group_by(ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME, year) %>%
	summarize(
		piped_mean = mean(piped),
		piped_lci = quantile(piped, probs = 0.025, na.rm = TRUE),
		piped_uci = quantile(piped, probs = 0.975, na.rm = TRUE),
		imp_mean = mean(imp),
		imp_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		imp_uci = quantile(imp, probs = 0.975, na.rm = TRUE))

setwd("<<<< FILEPATH REDACTED >>>>")
write.csv(bmgf_0, 'nga_a0_water_estimates.csv', row.names = FALSE)
write.csv(bmgf_1, 'nga_a1_water_estimates.csv', row.names = FALSE)
write.csv(bmgf_2, 'nga_a2_water_estimates.csv', row.names = FALSE)

bmgf_0 <- filter(intmd_0, ADM0_NAME == 'Nigeria') %>%
	mutate(sewer = s_piped,
		imp = s_piped + s_imp_other,
		od = s_od) %>%
	select(-s_piped, -s_imp_other, -s_od, -s_unimp) %>%
	group_by(ADM0_NAME, ADM0_CODE, year) %>%
	summarize(
		sewer_mean = mean(sewer),
		sewer_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		sewer_uci = quantile(imp, probs = 0.975, na.rm = TRUE),
		imp_mean = mean(imp),
		imp_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		imp_uci = quantile(imp, probs = 0.975, na.rm = TRUE),
		od_mean = mean(od),
		od_lci = quantile(od, probs = 0.025, na.rm = TRUE),
		od_uci = quantile(od, probs = 0.975, na.rm = TRUE))

bmgf_1 <- filter(intmd_1, ADM0_NAME == 'Nigeria') %>%
	mutate(sewer = s_piped,
		imp = s_piped + s_imp_other,
		od = s_od) %>%
	select(-s_piped, -s_imp_other, -s_od, -s_unimp) %>%
	group_by(ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME, year) %>%
	summarize(
		sewer_mean = mean(sewer),
		sewer_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		sewer_uci = quantile(imp, probs = 0.975, na.rm = TRUE),
		imp_mean = mean(imp),
		imp_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		imp_uci = quantile(imp, probs = 0.975, na.rm = TRUE),
		od_mean = mean(od),
		od_lci = quantile(od, probs = 0.025, na.rm = TRUE),
		od_uci = quantile(od, probs = 0.975, na.rm = TRUE))

bmgf_2 <- filter(intmd_2, ADM0_NAME == 'Nigeria') %>%
	mutate(sewer = s_piped,
		imp = s_piped + s_imp_other,
		od = s_od) %>%
	select(-s_piped, -s_imp_other, -s_od, -s_unimp) %>%
	group_by(ADM0_NAME, ADM0_CODE, ADM1_CODE, ADM1_NAME, ADM2_CODE, ADM2_NAME, year) %>%
	summarize(
		sewer_mean = mean(sewer),
		sewer_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		sewer_uci = quantile(imp, probs = 0.975, na.rm = TRUE),
		imp_mean = mean(imp),
		imp_lci = quantile(imp, probs = 0.025, na.rm = TRUE),
		imp_uci = quantile(imp, probs = 0.975, na.rm = TRUE),
		od_mean = mean(od),
		od_lci = quantile(od, probs = 0.025, na.rm = TRUE),
		od_uci = quantile(od, probs = 0.975, na.rm = TRUE))


setwd("<<<< FILEPATH REDACTED >>>>")
write.csv(bmgf_0, 'nga_a0_sani_estimates.csv', row.names = FALSE)
write.csv(bmgf_1, 'nga_a1_sani_estimates.csv', row.names = FALSE)
write.csv(bmgf_2, 'nga_a2_sani_estimates.csv', row.names = FALSE)
