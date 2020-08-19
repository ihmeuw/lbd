rm(list = ls())
file_path <- "<<<< FILEPATH REDACTED >>>>"

setwd(file_path)

library(data.table)
library(dplyr)
library(ggplot2)

process_data <- function(lvl) {
	a0 <- fread(paste0(file_path,'/sani/','a', lvl, '_results.csv'))

	a0_sani <- select(a0, year, paste0('ADM', lvl, '_NAME'),
		paste0('ADM', lvl, '_CODE'), ADM0_NAME, ADM0_CODE, piped, imp) %>%
		mutate(imp_sani = piped + imp) %>%
		select(-piped, -imp)

	a0 <- fread(paste0(file_path,'/water/','a', lvl, '_results.csv'))

	a0_water <- select(a0, year, paste0('ADM', lvl, '_NAME'),
		paste0('ADM', lvl, '_CODE'), ADM0_NAME, ADM0_CODE, piped, imp) %>%
		mutate(imp_water = piped + imp) %>%
		select(-piped, -imp)

	a0_results <- left_join(a0_sani, a0_water)

	if (lvl %in% 1:2) {

	}

	print(paste(nrow(a0_sani), nrow(a0_results)))
	return(a0_results)
}

a0_results <- process_data(0)
a1_results <- process_data(1)
a2_results <- process_data(2)

print_miss <- function(mydat) {
	mydat <- filter(mydat, is.na(imp_water)|is.na(imp_sani))
	print(unique(mydat$ADM0_NAME))

}

print_miss(a0_results)
print_miss(a1_results)
print_miss(a2_results)

a0_results <- filter(a0_results, year == 2017 & ADM0_NAME != 'India')
a1_results <- filter(a1_results, year == 2017 & ADM0_NAME != 'India')
a2_results <- filter(a2_results, year == 2017 & ADM0_NAME != 'India')
setwd("<<<< FILEPATH REDACTED >>>>")

write.csv(a0_results, 'a0_results.csv', row.names = FALSE)
write.csv(a1_results, 'a1_results.csv', row.names = FALSE)
write.csv(a2_results, 'a2_results.csv', row.names = FALSE)
