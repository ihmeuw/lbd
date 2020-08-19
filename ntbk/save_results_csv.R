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
	# x <- filter(x, year %in% c(2000, 2017))
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
regions <- c('dia_afr_horn', 'dia_essa-ken', 'ken',
	'dia_sssa', 'dia_cssa-cod', 'cod',
	'civ', 'dia_wssa-civ')
all_draws <- load_agg_levels(run_date, levels, regions)
sp_tbl <- do.call(bind_rows, lapply(regions, get_sp_tbl, ind = levels[1], rd = run_date))

all_draws2 <- filter(all_draws, year %in% c(2000, 2017)) %>%
	mutate(w_unimp = w_unimp + w_surface) %>%
	rename(w_imp = w_imp_other) %>%
	group_by(year, code, agg_level) %>%
	summarize(w_imp = mean(w_imp),
		w_piped = mean(w_piped),
		w_unimp = mean(w_unimp))

setwd("<<<< FILEPATH REDACTED >>>>")
write.csv(all_draws2, 'water_00_17.csv')
