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
	x <- filter(x, year %in% c(2000, 2017))
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
levels <- c('s_piped', 's_imp_other', 's_unimp', 's_od')
indi_group <- 'sani'
## WASH has to use two run dates since regions are split into two model runs
run_date <- '2019_09_29_00_00_01'
regions <- c('dia_essa-ken', 'ken', 'dia_afr_horn', 'dia_cssa-cod', 'cod', 'dia_wssa-civ', 'civ', 'dia_name', 'dia_sssa',
	'dia_mcaca-mex-cri', 'mex', 'cri', 'per', 'dia_central_asia', 'dia_chn_mng', 'dia_se_asia', 'dia_malay-idn-png',
	'idn', 'dia_south_asia', 'dia_mid_east')
all_draws <- load_agg_levels(run_date, levels, regions)
sp_tbl <- do.call(bind_rows, lapply(regions, get_sp_tbl, ind = levels[1], rd = run_date))

run_date <- '2019_05_20_00_00_02'
regions <- c('png', 'dia_s_america-per')
all_draws_0 <- load_agg_levels(run_date, levels, regions)
sp_tbl_0 <- do.call(bind_rows, lapply(regions, get_sp_tbl, ind = levels[1], rd = run_date))

all_draws <- bind_rows(all_draws, all_draws_0)
sp_tbl <- bind_rows(sp_tbl, sp_tbl_0)

#
if (indi_group  == 'water') {
	intmd <- mutate(all_draws,
	piped = w_piped,
	imp = w_imp_other,
	unimp = w_unimp + w_surface
	) %>%
	select(year, code, draw, agg_level, piped, imp, unimp)

	intmd_0 <- filter(intmd, agg_level == 0) %>%
		mutate(ADM0_CODE = code)
	intmd_1 <- filter(intmd, agg_level == 1) %>%
		mutate(ADM1_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_CODE, ADM1_CODE))) %>%
		select(-ADM1_CODE)
	intmd_2 <- filter(intmd, agg_level == 2) %>%
		mutate(ADM2_CODE = code) %>%
		left_join(distinct(select(sp_tbl, ADM0_CODE, ADM2_CODE))) %>%
		select(-ADM2_CODE)
	intmd <- bind_rows(intmd_0, intmd_1, intmd_2)

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

	intmd <- mutate(intmd,
		rrsum = p_b*1 + p_f*1.652 + p_n*2.401 +
		        i_b*1.118 + i_f*1.848 + i_n*2.685 +
		        u_b*1.364 + u_f*2.254 + u_n*3.276) %>%
		mutate(paf = (rrsum - 1)/rrsum)
}

if (indi_group == 'sani') {
	intmd <- mutate(all_draws,
		piped = s_piped,
		imp = s_imp_other,
		unimp = s_unimp + s_od
		) %>%
		mutate(rrsum = (piped*1) + (imp*2.595) + (unimp*3.242),
			paf = (rrsum -1)/rrsum)
}

intmd <- select(intmd, year, code, draw, agg_level, paf)
intmd$draw <- gsub('V', '', intmd$draw)
intmd <- spread(intmd, year, paf)
names(intmd)[4:5] <- c(paste0(indi_group, '_paf_2000'), paste0(indi_group, '_paf_2017'))

setwd("<<<< FILEPATH REDACTED >>>>")
write.csv(intmd, paste0(indi_group,'_paf_all_countries.csv'), row.names = FALSE)
