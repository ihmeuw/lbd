libs <- c('data.table', 'dplyr')
lapply(libs, library, character.only = TRUE)

# function to clean output csvs


# save out main indicators at admin levels for sanitation
setwd("<<<< FILEPATH REDACTED >>>>")
## admin0
for (i in c('piped', 'imp_total', 'unimp', 'no_facility')) {
	message(i)
	a0 <- as.data.frame(fread('a0_results.csv'))
	a0 <- a0[,c('year', 'ADM0_CODE', paste0(i, '_uci'),
		paste0(i, '_lci'))]
	a0$value <- a0[,paste0(i, '_uci')] - a0[,paste0(i, '_lci')]
	a0 <- select(a0, value, year, ADM0_CODE)
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
	a0 <- as.data.frame(fread('a1_results.csv'))
	a0 <- a0[,c('year', 'ADM1_CODE', paste0(i, '_uci'),
		paste0(i, '_lci'))]
	a0$value <- a0[,paste0(i, '_uci')] - a0[,paste0(i, '_lci')]
	a0 <- select(a0, value, year, ADM1_CODE)
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
	a0 <- as.data.frame(fread('a2_results.csv'))
	a0 <- a0[,c('year', 'ADM2_CODE', paste0(i, '_uci'),
		paste0(i, '_lci'))]
	a0$value <- a0[,paste0(i, '_uci')] - a0[,paste0(i, '_lci')]
	a0 <- select(a0, value, year, ADM2_CODE)
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
	a0 <- as.data.frame(fread('a0_results.csv'))
	a0 <- a0[,c('year', 'ADM0_CODE', paste0(i, '_uci'),
		paste0(i, '_lci'))]
	a0$value <- a0[,paste0(i, '_uci')] - a0[,paste0(i, '_lci')]
	a0 <- select(a0, value, year, ADM0_CODE)
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
	a0 <- as.data.frame(fread('a1_results.csv'))
	a0 <- a0[,c('year', 'ADM1_CODE', paste0(i, '_uci'),
		paste0(i, '_lci'))]
	a0$value <- a0[,paste0(i, '_uci')] - a0[,paste0(i, '_lci')]
	a0 <- select(a0, value, year, ADM1_CODE)
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
	a0 <- as.data.frame(fread('a2_results.csv'))
	a0 <- a0[,c('year', 'ADM2_CODE', paste0(i, '_uci'),
		paste0(i, '_lci'))]
	a0$value <- a0[,paste0(i, '_uci')] - a0[,paste0(i, '_lci')]
	a0 <- select(a0, value, year, ADM2_CODE)
	if (i == 'no_facility') {i <- 'surface'}
	if (i == 'imp_total') {i <- 'imp'}
	dir.create("<<<< FILEPATH REDACTED >>>>", recursive = TRUE)
	write.csv(
		a0, "<<<< FILEPATH REDACTED >>>>",
		row.names = FALSE
		)
}
