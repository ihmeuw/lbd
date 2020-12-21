# -------------------------------------------------------------
# Script to get numbers for manuscript text
# -------------------------------------------------------------


# ------------------------------------------------
# Set up

# clear space
rm(list=ls())

# set repo
repo <- '<<<< FILEPATH REDACTED >>>>'

# set model arguments
run_date          <- '2019_12_13_12_00_46'
shapefile_version <- '2019_09_10'

# set countries
countries <- c('Sierra Leone', 'Mali', 'Senegal')

# directories
dir <- paste0('<<<< FILEPATH REDACTED >>>>')
# ------------------------------------------------


# --------------------------------------------------------------------------------
# Packages and functions

# load required packages
libs <- list('raster', 'sf', 'sp', 'dplyr', 'fasterize', 'ggplot2', 'scales',
             'RColorBrewer', 'rgeos', 'rgdal', 'data.table', 'viridis')
lapply(libs, library, character.only = TRUE)
library(scico, lib.loc = '<<<< FILEPATH REDACTED >>>>')

# load custom functions
funcs <- paste0(repo, 'post_estimation/',
                list('squeeze_outputs', 'area_plots', 'anchor_maps',
                     'treatment_change_analysis', 'efficacy_analysis'),
                '.R')
lapply(funcs, source)

# get year function
get_years <- function(c) {
  y <- config[name == c, grep('date', names(config)), with = FALSE]
  y <- c(y[,1][[1]], y[,2][[1]], y[,3][[1]], y[,4][[1]])
  return(y)
}
# --------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------------------------------------
# Load all data summaries

# load config
config <- fread(paste0(repo, 'post_estimation/key_dates.csv'))

# load any ORS coverage estimates
ors0 <- fread(paste0('<<<< FILEPATH REDACTED >>>>'))
ors2 <- fread(paste0('<<<< FILEPATH REDACTED >>>>'))

# load only RHF coverage estimates
rhf0 <- fread(paste0('<<<< FILEPATH REDACTED >>>>'))
rhf2 <- fread(paste0('<<<< FILEPATH REDACTED >>>>'))

# load no ORT coverage estimates
nort0 <- fread(paste0('<<<< FILEPATH REDACTED >>>>'))
nort2 <- fread(paste0('<<<< FILEPATH REDACTED >>>>'))

# load summary tables
dt1 <- fread(paste0(dir, 'cleaned_treatment_change_summary_table.csv'))
dt2 <- fread(paste0(dir, 'cleaned_rhf_replacement_summary_table_v2.csv'))
dt3 <- fread(paste0(dir, 'summary_data_for_figs.csv'))

# create subnational groups
groups_master <- list(
  list(
    'N' = c('Tonkolili', 'Port Loko', 'Koinadugu', 'Kambia', 'Bombali'),
    'C' = c('Western Rural', 'Western Urban'),
    'S' = c('Bo', 'Bonthe', 'Moyamba', 'Pujehun', 'Kailahun', 'Kenema', 'Kono')
  ),
  groups <- list(
    'N' = unique(ors2[ADM1_NAME %in% c('Timbuktu', 'Kidal', 'Gao', 'Mopti'), ADM2_NAME]),
    'S' = unique(ors2[ADM1_NAME %in% c('Kayes', 'Koulikoro', 'Sikasso', 'Ségou'), ADM2_NAME]),
    'C' = unique(ors2[ADM2_NAME %in% c('Bamako'), ADM2_NAME])
  ),
  groups <- list(
    'N' = unique(ors2[ADM1_NAME %in% c('Saint-Louis', 'Louga', 'Matam'), ADM2_NAME]),
    'C' = unique(ors2[ADM1_NAME %in% c('Thiès', 'Dakar', 'Kaolack', 'Fatick', 'Kaffrine', 'Thiès', 'Dakar', 'Diourbel', 'Tambacounda', 'Kédougou', 'Kolda'), ADM2_NAME]),
    'S' = unique(ors2[ADM1_NAME %in% c('Sédhiou', 'Ziguinchor'), ADM2_NAME])
  )
)
names(groups_master) <- countries
# ------------------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# National-level changes

message('any ORS')
for (i in countries) {
  message(i)
  print(ors0[ADM0_NAME == i & year == config[name == i, date_start],
             c('year', 'mean', 'lower', 'upper')])
  print(ors0[ADM0_NAME == i & year == config[name == i, date_end],
             c('year', 'mean', 'lower', 'upper')])
}

message('no ORT')
for (i in countries) {
  message(i)
  print(nort0[ADM0_NAME == i & year == config[name == i, date_start],
              c('year', 'mean', 'lower', 'upper')])
  print(nort0[ADM0_NAME == i & year == config[name == i, date_end],
              c('year', 'mean', 'lower', 'upper')])
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Changes within countries (ranges)

message('any ORS')
for (c in countries) {
  message(c)
  
  # get years
  y <- get_years(c)
  
  # get groups
  groups <- groups_master[[c]]
  
  # clean data
  dat <- copy(ors2)
  dat <- dat[ADM0_NAME == c & year %in% y]
  dat[, year := as.factor(year)]
  if (c != 'Senegal') {
    for (i in names(groups)) dat[ADM2_NAME %in% groups[[i]], group := i]
  } else {
    dat[, group := 'all']
  }
  dat[, c('min_mean', 'min_lower', 'min_upper') := lapply(.SD, min), 
      by = c('group', 'year'), .SDcols = c('mean', 'lower', 'upper')]
  dat[, c('max_mean', 'max_lower', 'max_upper') := lapply(.SD, max), 
      by = c('group', 'year'), .SDcols = c('mean', 'lower', 'upper')]
  dat <- unique(dat[, c('group', 'year', 'min_mean', 'min_lower', 'min_upper',
                        'max_mean', 'max_lower', 'max_upper')])
  print(dat)
}

message('only RHF')
for (c in countries) {
  message(c)
  
  # get years
  y <- get_years(c)
  
  # get groups
  groups <- groups_master[[c]]
  
  # clean data
  dat <- copy(rhf2)
  dat <- dat[ADM0_NAME == c & year %in% y]
  dat[, year := as.factor(year)]
  if (c != 'Senegal') {
    for (i in names(groups)) dat[ADM2_NAME %in% groups[[i]], group := i]
  } else {
    dat[, group := 'all']
  }
  dat[, c('min_mean', 'min_lower', 'min_upper') := lapply(.SD, min), 
      by = c('group', 'year'), .SDcols = c('mean', 'lower', 'upper')]
  dat[, c('max_mean', 'max_lower', 'max_upper') := lapply(.SD, max), 
      by = c('group', 'year'), .SDcols = c('mean', 'lower', 'upper')]
  dat <- unique(dat[, c('group', 'year', 'min_mean', 'min_lower', 'min_upper',
                        'max_mean', 'max_lower', 'max_upper')])
  print(dat)
}
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Changes within admin 2 units (specific districts)

message('ORS decline in Freetown during Ebola')
print(ors2[ADM2_NAME == 'Western Urban' & year %in% c(2013, 2017)])

message('ORS change in areas of intervention')
print(ors2[ADM2_NAME == 'Bougouni' & year %in% c(2001, 2004)])
print(ors2[ADM2_NAME == 'Sikasso' & year %in% c(2001, 2004)])
print(ors2[ADM2_NAME == 'Bla' & year %in% c(2001, 2004)])

message('Max coverage in Mali')
print(setorderv(ors2[ADM0_NAME == 'Mali' & year == 2018], 'mean', order = -1)[1,])
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Posterior probabilities > 95%

# change within time periods for ORS
message('any ORS')
for (c in countries) {
  message('\n\n', c)
  
  tmp <- dt1[ADM0_NAME == c]
  groups <- groups_master[[c]]
  for (i in names(groups)) tmp[ADM2_NAME %in% groups[[i]], group := i]
  
  message('\nincreases')
  for (i in c('N', 'C', 'S')) {
    message(i)
    message('Number of ', i, ' districts: ', length(unique(tmp[group == i, ADM2_NAME])))
    message('Number with significant increases before policy change: ', 
            length(unique(tmp[group == i & period == 'pre' & lower_any_ors > 0, ADM2_NAME])))
    message('Number with significant increases during policy change: ', 
            length(unique(tmp[group == i & period == 'during' & lower_any_ors > 0, ADM2_NAME])))
    message('Number with significant increases after policy change: ', 
            length(unique(tmp[group == i & period == 'after' & lower_any_ors > 0, ADM2_NAME])))
  }
  
  message('\ndecreases')
  for (i in c('N', 'C', 'S')) {
    message(i)
    message('Number of ', i, ' districts: ', length(unique(tmp[group == i, ADM2_NAME])))
    message('Number with significant decreases before policy change: ', 
            length(unique(tmp[group == i & period == 'pre' & upper_any_ors < 0, ADM2_NAME])))
    message('Number with significant decreases during policy change: ', 
            length(unique(tmp[group == i & period == 'during' & upper_any_ors < 0, ADM2_NAME])))
    message('Number with significant decreases after policy change: ', 
            length(unique(tmp[group == i & period == 'after' & upper_any_ors < 0, ADM2_NAME])))
  }
  
  # check which districts in Mali showed significant increases
  message('Significant increases in ORS in N Mali after war')
  if (c == 'Mali') print(tmp[group == 'N' & period == 'after' & lower_any_ors > 0])
  
  # check which districts that showed or didn't show change during Ebola
  message('No change in North during Ebola in N')
  if (c == 'Sierra Leone') print(tmp[group == 'N' & period == 'after' & upper_any_ors > 0])
  message('Change in South during Ebola')
  if (c == 'Sierra Leone') print(tmp[group == 'S' & period == 'after' & upper_any_ors < 0])
}

# change within time periods for RHF
message('RHF only')
for (c in countries) {
  message('\n\n', c)
  
  tmp <- dt1[ADM0_NAME == c]
  groups <- groups_master[[c]]
  for (i in names(groups)) tmp[ADM2_NAME %in% groups[[i]], group := i]
  
  message('\nincreases')
  for (i in c('N', 'C', 'S')) {
    message(i)
    message('Number of ', i, ' districts: ', length(unique(tmp[group == i, ADM2_NAME])))
    message('Number with significant increases before policy change: ', 
            length(unique(tmp[group == i & period == 'pre' & lower_rhf_only > 0, ADM2_NAME])))
    message('Number with significant increases during policy change: ', 
            length(unique(tmp[group == i & period == 'during' & lower_rhf_only > 0, ADM2_NAME])))
    message('Number with significant increases after policy change: ', 
            length(unique(tmp[group == i & period == 'after' & lower_rhf_only > 0, ADM2_NAME])))
  }
  
  message('\ndecreases')
  for (i in c('N', 'C', 'S')) {
    message(i)
    message('Number of ', i, ' districts: ', length(unique(tmp[group == i, ADM2_NAME])))
    message('Number with significant decreases before policy change: ', 
            length(unique(tmp[group == i & period == 'pre' & upper_rhf_only < 0, ADM2_NAME])))
    message('Number with significant decreases during policy change: ', 
            length(unique(tmp[group == i & period == 'during' & upper_rhf_only < 0, ADM2_NAME])))
    message('Number with significant decreases after policy change: ', 
            length(unique(tmp[group == i & period == 'after' & upper_rhf_only < 0, ADM2_NAME])))
  }
  
  # check which districts in Mali showed significant changes
  message('Significant declines in RHF in S Mali during interventions')
  if (c == 'Mali') print(tmp[group == 'S' & period == 'pre' & upper_rhf_only < 0])
  
}

# check to see if significant changes in RHF in Mali intervention districts
print(dt1[ADM2_NAME %in% c('Bla', 'Sikasso', 'Bougouni') & period == 'pre'])
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Numbers untreated

# load and prep ad0 data
dia <- fread('<<<< FILEPATH REDACTED >>>>/had_diarrhea_admin_0_raked_prevalence_summary.csv')
dia <- dia[year == 2017 & (ADM0_NAME == 'Senegal' | ADM0_NAME == 'Mali')]
dia[, num_kids := mean*pop]
dt5 <- merge(ors0, dia[, c('ADM0_CODE', 'year', 'num_kids')], by = 'ADM0_CODE', allow.cartesian = T)

no_ort_sub <- nort0[((year == 2000 | year == 2017) & ADM0_NAME == 'Senegal') | 
                      (ADM0_NAME == 'Mali' & (year == 2001 | year == 2018))]
no_ort_sub[, year := ifelse(year < 2005, 'start', 'end')]
no_ort_sub <- dcast(no_ort_sub, formula = ADM0_CODE ~ year, value.var = c('mean', 'upper', 'lower'))

dt5 <- merge(dt5, no_ort_sub, by = 'ADM0_CODE')

dt5[, no_ors_mean := (1 - mean) * num_kids, by = .I]
dt5[, no_ors_lower := (1 - upper) * num_kids, by = .I]
dt5[, no_ors_upper := (1 - lower) * num_kids, by = .I]

dt5[, would_have_ort_mean := no_ors_mean*(1-mean_start)]
dt5[, would_have_ort_upper := no_ors_upper*(1-lower_start)]
dt5[, would_have_ort_lower := no_ors_lower*(1-upper_start)]

# print summaries for admin 0
message('did not receive any ORS')
for (i in c('Mali', 'Senegal')) {
  message(i)
  print(dt5[ADM0_NAME == i & year.x == config[name == i, date_end],
             c('year.x', 'no_ors_mean', 'no_ors_lower', 'no_ors_upper')])
}

message('would have received ORT at start of study')
for (i in c('Mali', 'Senegal')) {
  message(i)
  print(dt5[ADM0_NAME == i & year.x == config[name == i, date_end],
            c('year.x', 'would_have_ort_mean', 'would_have_ort_lower', 'would_have_ort_upper')])
}

# load and prep ad2 data
dia <- fread('<<<< FILEPATH REDACTED >>>>/had_diarrhea_admin_2_raked_prevalence_summary.csv')
dia <- dia[year == 2017 & (ADM0_NAME == 'Senegal' | ADM0_NAME == 'Mali')]
dia[, num_kids := mean*pop]
dt4 <- merge(ors2, dia[, c('ADM2_CODE', 'year', 'num_kids')], by = 'ADM2_CODE', allow.cartesian = T)

no_ort_sub <- nort2[((year == 2000 | year == 2017) & ADM0_NAME == 'Senegal') | 
                      (ADM0_NAME == 'Mali' & (year == 2001 | year == 2018))]
no_ort_sub[, year := ifelse(year < 2005, 'start', 'end')]
no_ort_sub <- dcast(no_ort_sub, formula = ADM2_CODE ~ year, value.var = c('mean', 'upper', 'lower'))

dt4 <- merge(dt4, no_ort_sub, by = 'ADM2_CODE')

dt4[, no_ors_mean := (1 - mean) * num_kids, by = .I]
dt4[, no_ors_lower := (1 - upper) * num_kids, by = .I]
dt4[, no_ors_upper := (1 - lower) * num_kids, by = .I]

dt4[, would_have_ort_mean := no_ors_mean*(1-mean_start)]
dt4[, would_have_ort_upper := no_ors_upper*(1-lower_start)]
dt4[, would_have_ort_lower := no_ors_lower*(1-upper_start)]

# print summaries for admin 2 units where the most lives would have been saved
message('would have received ORT at start of study')
for (i in c('Mali', 'Senegal')) {
  message(i)
  tmp <- dt4[ADM0_NAME == i & year.x == config[name == i, date_end]]
  setorderv(tmp, 'would_have_ort_mean', order = -1)
  print(tmp[1:3, c('ADM1_NAME', 'ADM2_NAME', 'would_have_ort_mean', 'would_have_ort_lower', 'would_have_ort_upper')])
}
# ------------------------------------------------------------------------------












