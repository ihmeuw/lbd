rm(list = ls())
library(dplyr)

gbd <- read.csv("<<<< FILEPATH REDACTED >>>>",
	stringsAsFactors = FALSE)

lbd <- read.csv("<<<< FILEPATH REDACTED >>>>",
	stringsAsFactors = FALSE)

w_typology <- gbd %>%
	filter(indicator %in% c('w_source_drink', 'variable_name')) %>%
	select(string, sdg) %>%
	rename(gbd = sdg) %>%
	distinct() %>%
	left_join(select(rename(lbd, lbd = sdg), string, lbd)) %>%
	distinct() %>%
	filter((gbd != lbd)|(is.na(gbd))|(is.na(lbd))) %>%
	select(lbd, gbd) %>%
	distinct() %>%
	arrange(lbd)
w_mismatch <- w_typology <- gbd %>%
	filter(indicator %in% c('w_source_drink', 'variable_name')) %>%
	select(string, sdg) %>%
	rename(gbd = sdg) %>%
	distinct() %>%
	left_join(select(rename(lbd, lbd = sdg), string, lbd)) %>%
	distinct() %>%
	filter((gbd != lbd)|(is.na(gbd))|(is.na(lbd)))

lbd <- read.csv("<<<< FILEPATH REDACTED >>>>",
	stringsAsFactors = FALSE)

s_typology <- gbd %>%
	filter(indicator %in% c('t_type', 'variable_name')) %>%
	select(string, sdg) %>%
	rename(gbd = sdg) %>%
	distinct() %>%
	left_join(select(rename(lbd, lbd = sdg_2), string, lbd)) %>%
	distinct() %>%
	filter((gbd != lbd)|(is.na(gbd))|(is.na(lbd))) %>%
	select(lbd, gbd) %>%
	distinct() %>%
	arrange(lbd)
s_mismatch <- w_typology <- gbd %>%
	filter(indicator %in% c('t_type', 'variable_name')) %>%
	select(string, sdg) %>%
	rename(gbd = sdg) %>%
	distinct() %>%
	left_join(select(rename(lbd, lbd = sdg_2), string, lbd)) %>%
	distinct() %>%
	filter((gbd != lbd)|(is.na(gbd))|(is.na(lbd)))

setwd("<<<< FILEPATH REDACTED >>>>")

write.csv(w_typology, 'w_typology.csv', row.names = FALSE)
write.csv(w_mismatch, 'w_mismatch.csv', row.names = FALSE)
write.csv(s_typology, 's_typology.csv', row.names = FALSE)
write.csv(s_mismatch, 's_mismatch.csv', row.names = FALSE)
