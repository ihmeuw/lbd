# -------------------------------------------------------------
# Script to make additional tables for manuscript
# -------------------------------------------------------------

# set up --------------------------------------------------

library(data.table)

repo <- '<<<< FILEPATH REDACTED >>>>'
indics <- c('any_ors', 'rhf_only', 'no_ort')
isos <- c('MLI', 'SEN', 'SLE')
countries <- c('Sierra Leone', 'Mali', 'Senegal')

config <- fread(paste0(repo, 'anchor/key_dates.csv'))[, iso3 := NULL]

# get year function
get_years <- function(c) {
  y <- config[name == c, grep('date', names(config)), with = FALSE]
  y <- c(y[,1][[1]], y[,2][[1]], y[,3][[1]], y[,4][[1]])
  return(y)
}

# covariate tables ------------------------------------------------

setwd('<<<< FILEPATH REDACTED >>>>')

for (i in indics) {

  cov <- rbindlist(lapply(isos, 
                          function(c) {
                            fread(paste0('covs_', i, '/covs_ort_', c, '.csv'))[, country := c]
                            }
                          ))

  cov <- dcast(cov, covariate ~ country, value.var = 'include')
  
  write.csv(cov, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_8_', i, '.csv'))

}

# input data table ----------------------------------------------------------

dt <- fread('<<<< FILEPATH REDACTED >>>>/no_ort.csv')
dt <- unique(dt[country %in% isos, c('country', 'nid', 'year', 'source')])
write.csv(dt, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_2.csv'))

# admin 0 coverage estimates ----------------------------------------------------

dt <- rbindlist(lapply(indics,
                       function(i) {
                         fread(paste0('<<<< FILEPATH REDACTED >>>>', 
                                      i, '_admin_0_squeezed_summary.csv'))[, indicator := i]
                       }
)
)

dia <- fread('<<<< FILEPATH REDACTED >>>>/had_diarrhea_admin_0_raked_prevalence_summary.csv')
dia[, num_kids := mean*pop]
dt <- merge(dt, dia[, c('ADM0_CODE', 'year', 'num_kids')], by = c('ADM0_CODE', 'year'), 
            allow.cartesian = T, all.x = T)

kids17mli <- dia[ADM0_NAME == 'Mali' & year == 2017, num_kids]
dt[is.na(num_kids), num_kids := kids17mli]

dt[, coverage_ui := paste0(round(mean*100, 1), ' (', 
                           round(lower*100, 1), ' to ', 
                           round(upper*100, 1), ')'), by = .I]

dt[indicator == 'no_ort', untreated_ui := paste0(round((mean)*num_kids/1000, 1), ' (',
                                                 round((lower)*num_kids/1000, 1), ' to ',
                                                 round((upper)*num_kids/1000, 1), ')'), by = .I] # in thousands

dt[indicator != 'no_ort', untreated_ui := paste0(round((1-mean)*num_kids/1000, 1), ' (',
                                                 round((1-upper)*num_kids/1000, 1), ' to ',
                                                 round((1-lower)*num_kids/1000, 1), ')'), by = .I] # in thousands

dt2 <- dcast(dt,  ADM0_NAME + year ~ indicator, value.var = 'coverage_ui')
dt3 <- dcast(dt,  ADM0_NAME + year ~ indicator, value.var = 'untreated_ui')[, rhf_only := NULL]
dt <- merge(dt2, dt3)

write.csv(dt, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_9.csv'))

# admin 2 coverage estimates ------------------------------------------------

dt <- rbindlist(lapply(indics,
                       function(i) {
                         fread(paste0('<<<< FILEPATH REDACTED >>>>', 
                                      i, '_admin_2_squeezed_summary.csv'))[, indicator := i]
                       }
                       )
)

groups_master <- list(
  list(
    'B' = c('Tonkolili', 'Port Loko', 'Koinadugu', 'Kambia', 'Bombali'),
    'A' = c('Western Rural', 'Western Urban'),
    'C' = c('Bo', 'Bonthe', 'Moyamba', 'Pujehun', 'Kailahun', 'Kenema', 'Kono')
  ),
  groups <- list(
    'B' = unique(dt[ADM1_NAME %in% c('Timbuktu', 'Kidal', 'Gao', 'Mopti'), ADM2_NAME]),
    'C' = unique(dt[ADM1_NAME %in% c('Kayes', 'Koulikoro', 'Sikasso', 'Ségou'), ADM2_NAME]),
    'A' = unique(dt[ADM2_NAME %in% c('Bamako'), ADM2_NAME])
  ),
  groups <- list(
    'A' = unique(dt[ADM1_NAME %in% c('Saint-Louis', 'Louga', 'Matam'), ADM2_NAME]),
    'B' = unique(dt[ADM1_NAME %in% c('Thiès', 'Dakar', 'Kaolack', 'Fatick', 'Kaffrine', 'Thiès', 'Dakar', 'Diourbel', 'Tambacounda', 'Kédougou', 'Kolda'), ADM2_NAME]),
    'C' = unique(dt[ADM1_NAME %in% c('Sédhiou', 'Ziguinchor'), ADM2_NAME])
  )
)
names(groups_master) <- countries

dt[, mean_ui := paste0(round(mean*100, 1), ' (', round(lower*100, 1), ' to ', round(upper*100, 1), ')'), by = .I]

clean_data <- function(c) {
  groups <- groups_master[[c]]
  tmp <- dt[ADM0_NAME == c & year %in% get_years(c)]
  for (i in names(groups)) tmp[ADM2_NAME %in% groups[[i]], group := i]
  tmp <- dcast(tmp, group + ADM1_NAME + ADM2_NAME + indicator ~ year, value.var = 'mean_ui')
  tmp$indicator <- as.factor(tmp$indicator)
  levels(tmp$indicator) <- indics
  setorderv(tmp, c('group', 'ADM1_NAME', 'ADM2_NAME', 'indicator'))
  tmp[indicator == 'any_ors', indicator := 'Any ORS']
  tmp[indicator == 'rhf_only', indicator := 'Only RHF']
  tmp[indicator == 'no_ort', indicator := 'No ORT']
  return(tmp)
}

mli <- clean_data('Mali')
sle <- clean_data('Sierra Leone')
sen <- clean_data('Senegal')

write.csv(sle, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_10A.csv'))
write.csv(mli, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_10B.csv'))
write.csv(sen, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_10C.csv'))

# posterior probability of change between each time period -------------------------------------------

dt2 <- fread('<<<< FILEPATH REDACTED >>>>/cleaned_treatment_change_summary_table.csv')
dt2 <- merge(dt2, dt[, c('ADM2_NAME', 'ADM1_NAME')], by = 'ADM2_NAME', allow.cartesian = T)

for (i in indics) {
  
  dt2[get(paste0('lower_', i)) > 0, (paste0(i, '_result')) := 'Increase']
  dt2[get(paste0('upper_', i)) < 0, (paste0(i, '_result')) := 'Decrease']
  dt2[, (i) := paste0(round(get(paste0('mean_', i))*100, 1), ' (', 
                      round(get(paste0('lower_', i))*100, 1), ' to ', 
                      round(get(paste0('upper_', i))*100, 1), ')'), 
                         
      by = .I]
}

dt3 <- melt(dt2, id.vars = c('ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME', 'period'),
            measure.vars = c('any_ors', 'rhf_only', 'no_ort'), variable.name = 'indicator', value.name = 'mean_ui')
dt4 <- melt(dt2, id.vars = c('ADM0_NAME', 'ADM1_NAME', 'ADM2_NAME', 'period'),
            measure.vars = c('any_ors_result', 'rhf_only_result', 'no_ort_result'), variable.name = 'indicator', value.name = 'result')
dt4[is.na(result), result := 'Uncertain']
dt4[, indicator := gsub('_result', '', indicator)]

clean_data <- function(c) {
  groups <- groups_master[[c]]
  tmp1 <- unique(dt3[ADM0_NAME == c])
  for (i in names(groups)) tmp1[ADM2_NAME %in% groups[[i]], group := i]
  tmp2 <- unique(dt4[ADM0_NAME == c])
  for (i in names(groups)) tmp2[ADM2_NAME %in% groups[[i]], group := i]
  tmp1 <- dcast(tmp1, group + ADM1_NAME + ADM2_NAME + indicator ~ period, value.var = 'mean_ui')
  tmp2 <- dcast(tmp2, group + ADM1_NAME + ADM2_NAME + indicator ~ period, value.var = 'result')
  setnames(tmp2, c('after', 'during', 'pre'), c('after_pp', 'during_pp', 'pre_pp'))
  tmp <- merge(tmp1, tmp2)
  tmp$indicator <- as.factor(tmp$indicator)
  levels(tmp$indicator) <- indics
  setorderv(tmp, c('group', 'ADM1_NAME', 'ADM2_NAME', 'indicator'))
  tmp[indicator == 'any_ors', indicator := 'Any ORS']
  tmp[indicator == 'rhf_only', indicator := 'Only RHF']
  tmp[indicator == 'no_ort', indicator := 'No ORT']
  setcolorder(tmp, c('group', 'ADM1_NAME', 'ADM2_NAME', 'indicator', 'pre', 'pre_pp', 'during', 'during_pp', 'after', 'after_pp'))
  return(tmp)
}

mli <- clean_data('Mali')
sle <- clean_data('Sierra Leone')
sen <- clean_data('Senegal')

write.csv(sle, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_11A.csv'))
write.csv(mli, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_11B.csv'))
write.csv(sen, paste0('<<<< FILEPATH REDACTED >>>>/SI_table_11C.csv'))

# admin 2 csvs for additional file 1 ArcGIS number untreated maps ----------------------------------------------------

dt <- rbindlist(lapply(indics,
                       function(i) {
                         fread(paste0('<<<< FILEPATH REDACTED >>>>', 
                                      i, '_admin_2_squeezed_summary.csv'))[, indicator := i]
                       }
)
)

dia <- fread('<<<< FILEPATH REDACTED >>>>/had_diarrhea_admin_2_raked_prevalence_summary.csv')
dia[, num_kids := mean*pop]
dt <- merge(dt, dia[, c('ADM2_CODE', 'year', 'num_kids')], by = c('ADM2_CODE', 'year'), 
            allow.cartesian = T, all.x = T)

kids17mli <- dia[ADM0_NAME == 'Mali' & year == 2017, c('ADM2_CODE', 'num_kids')]
setnames(kids17mli, 'num_kids', 'num_kids17')
dt <- merge(dt, kids17mli, by = 'ADM2_CODE', all.x = T, allow.cartesian = T)
dt[is.na(num_kids), num_kids := num_kids17]

dt[indicator == 'no_ort', num_untreated := mean*num_kids]
dt[indicator != 'no_ort', num_untreated := (1-mean)*num_kids]

dt <- rbindlist(lapply(countries,
                       function(c) {
                         return(dt[ADM0_NAME == c & year %in% get_years(c)])
                         }
                       )
)

dt1 <- dt[, c('ADM2_CODE', 'year', 'mean', 'indicator')]
dt2 <- dt[, c('ADM2_CODE', 'year', 'num_untreated', 'indicator')]
setnames(dt1, 'mean', 'value')
setnames(dt2, 'num_untreated', 'value')

map_dir <- '<<<< FILEPATH REDACTED >>>>'

for (i in indics) {
  
  tmp1 <- dt1[indicator == i]
  tmp1$indicator <- NULL
  write.csv(tmp1, paste0(map_dir, i, '/2019_12_13_12_00_46/', i, '_proportion_mean_ad2.csv'))
  
  tmp2 <- dt2[indicator == i]
  tmp2$indicator <- NULL
  write.csv(tmp2, paste0(map_dir, i, '/2019_12_13_12_00_46/', i, '_count_mean_ad2.csv'))
  
}
