rm(list = ls())
libs <- c('data.table', 'dplyr')
lapply(libs, library, character.only = TRUE)

# function to clean output csvs


# save out main indicators at admin levels for sanitation
setwd("<<<< FILEPATH REDACTED >>>>")
## admin0
for (i in c('piped', 'imp_total', 'unimp', 'no_facility')) {
	message(i)
	a0 <- fread('a0_results.csv') %>%
		select_('year', 'ADM0_CODE', i) %>%
		rename(value = i)
	if (i == 'no_facility') {i <- 'od'}
	if (i == 'imp_total') {i <- 'imp'}
	dir.create("<<<< FILEPATH REDACTED >>>>", recursive = TRUE)
	write.csv(
		a0, "<<<< FILEPATH REDACTED >>>>",
		row.names = FALSE
		)
}

## admin1
for (i in c('piped', 'imp_total', 'unimp', 'no_facility')) {
	message(i)
	a0 <- fread('a1_results.csv') %>%
		select_('year', 'ADM1_CODE', i) %>%
		rename(value = i)
	if (i == 'no_facility') {i <- 'od'}
	if (i == 'imp_total') {i <- 'imp'}
	dir.create("<<<< FILEPATH REDACTED >>>>", recursive = TRUE)
	write.csv(
		a0, "<<<< FILEPATH REDACTED >>>>",
		row.names = FALSE
		)
}

## admin2
for (i in c('piped', 'imp_total', 'unimp', 'no_facility')) {
	message(i)
	a0 <- fread('a2_results.csv') %>%
		select_('year', 'ADM2_CODE', i) %>%
		rename(value = i)
	if (i == 'no_facility') {i <- 'od'}
	if (i == 'imp_total') {i <- 'imp'}
	dir.create("<<<< FILEPATH REDACTED >>>>", recursive = TRUE)
	write.csv(
		a0, "<<<< FILEPATH REDACTED >>>>",
		row.names = FALSE
		)
}

# save out main indicators at admin levels for water
setwd("<<<< FILEPATH REDACTED >>>>")
## admin0
for (i in c('piped', 'imp_total', 'unimp', 'no_facility')) {
	message(i)
	a0 <- fread('a0_results.csv') %>%
		select_('year', 'ADM0_CODE', i) %>%
		rename(value = i)
	if (i == 'no_facility') {i <- 'surface'}
	if (i == 'imp_total') {i <- 'imp'}
	dir.create("<<<< FILEPATH REDACTED >>>>", recursive = TRUE)
	write.csv(
		a0, "<<<< FILEPATH REDACTED >>>>",
		row.names = FALSE
		)
}

## admin1
for (i in c('piped', 'imp_total', 'unimp', 'no_facility')) {
	message(i)
	a0 <- fread('a1_results.csv') %>%
		select_('year', 'ADM1_CODE', i) %>%
		rename(value = i)
	if (i == 'no_facility') {i <- 'surface'}
	if (i == 'imp_total') {i <- 'imp'}
	dir.create("<<<< FILEPATH REDACTED >>>>", recursive = TRUE)
	write.csv(
		a0, "<<<< FILEPATH REDACTED >>>>",
		row.names = FALSE
		)
}

## admin2
for (i in c('piped', 'imp_total', 'unimp', 'no_facility')) {
	message(i)
	a0 <- fread('a2_results.csv') %>%
		select_('year', 'ADM2_CODE', i) %>%
		rename(value = i)
	if (i == 'no_facility') {i <- 'surface'}
	if (i == 'imp_total') {i <- 'imp'}
	dir.create("<<<< FILEPATH REDACTED >>>>", recursive = TRUE)
	write.csv(
		a0, "<<<< FILEPATH REDACTED >>>>",
		row.names = FALSE
		)
}
