#####################################################################
# Combine aggregated DALYs into csvs formatted for analysis
#####################################################################

#(1) Setup ---------------------------------------------------------
rm(list = ls())

#set user arguments
causes <- c('diarrhea','lri','malaria')
admins <- c(0,1,2)
check_csvs <- TRUE #if true, will check number of countries in each csv
dia_run_date <- '2019_09_17_14_12_53'
lri_run_date <- '2019_09_16_16_39_00'
mal_run_date <- '2019_10_28'
tridaly_run_date <- '2019_10_28'
lri_dir <- '<<< FILEPATH REDACTED >>>'
dia_dir <- '<<< FILEPATH REDACTED >>>'
mal_dir <- '<<< FILEPATH REDACTED >>>'
tridaly_dir <- '<<< FILEPATH REDACTED >>>'
exclude_countries <- c('Western Sahara', 'Yemen')

#load libraries and packages
library(dplyr)
library(data.table)

for (cause in causes){
  for (admin_lvl in admins){
    
    if (cause == 'diarrhea') {
      cause_dir <- dia_dir
      indicator_group <- 'ort'
      indicator <- 'had_diarrhea'
      run_date <- dia_run_date
      }
    
    if (cause == 'lri') {
      cause_dir <- lri_dir
      indicator_group <- 'lri'
      indicator <- 'has_lri'
      run_date <- lri_run_date
    }
    
    if (cause == 'malaria') {
      cause_dir <- mal_dir
      indicator_group <- 'malaria'
      indicator <- 'had_malaria'
      run_date <- mal_run_date
    }
    
    if (cause %in% c('lri', 'diarrhea')){
      
      #load in dalys, ylls, yld from tridalys dir
      
      ylds <- fread(paste0(tridaly_dir, indicator, '_admin_', admin_lvl, '_raked_yld_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylds_rate_mean = mean, ylds_rate_upper = upper, ylds_rate_lower = lower)
      
      ylls <- fread(paste0(tridaly_dir, indicator, '_admin_', admin_lvl, '_raked_yll_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylls_rate_mean = mean, ylls_rate_upper = upper, ylls_rate_lower = lower)
      
      dalys <- fread(paste0(tridaly_dir, indicator, '_admin_', admin_lvl, '_raked_daly_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(dalys_rate_mean = mean, dalys_rate_upper = upper, dalys_rate_lower = lower)
      
      ylds_c <- fread(paste0(tridaly_dir, indicator, '_c_admin_', admin_lvl, '_raked_yld_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylds_count_mean = mean, ylds_count_upper = upper, ylds_count_lower = lower)
      
      ylls_c <- fread(paste0(tridaly_dir, indicator, '_c_admin_', admin_lvl, '_raked_yll_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylls_count_mean = mean, ylls_count_upper = upper, ylls_count_lower = lower)
      
      dalys_c <- fread(paste0(tridaly_dir, indicator, '_c_admin_', admin_lvl, '_raked_daly_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(dalys_count_mean = mean, dalys_count_upper = upper, dalys_count_lower = lower)
      
      #merge
      loc_info <- grep(names(dalys), pattern = 'ADM', value = T)
      
      admin_dalys <- merge(ylds, ylls, by = c('year', loc_info)) %>%
        merge(dalys, by = c('year', loc_info)) %>%
        merge(ylds_c, by = c('year', loc_info)) %>%
        merge(ylls_c, by = c('year', loc_info)) %>%
        merge(dalys_c, by = c('year', loc_info))
      
      #(3) Read in aggregated mortality and incidence rates, counts -----
      
      if (cause == 'lri'){
        load('<<< FILEPATH REDACTED >>>')
        pop <- get(paste0('admin_', admin_lvl)) %>%
          select(c('year',paste0('ADM', admin_lvl, '_CODE'), 'pop'))
        
        mort <- fread(paste0(cause_dir, indicator, '_admin_', admin_lvl, '_raked_mortality_summary.csv')) %>%
          dplyr::select(-region, -cirange) %>%
          rename(mort_rate_mean = mean, mort_rate_upper = upper, mort_rate_lower = lower) %>%
          merge(pop, by = c('year',paste0('ADM', admin_lvl, '_CODE')))
        
        inc <- fread(paste0(cause_dir, indicator, '_admin_', admin_lvl, '_raked_mortality_summary.csv')) %>%
          dplyr::select(-region, -cirange) %>%
          rename(inc_rate_mean = mean, inc_rate_upper = upper, inc_rate_lower = lower)
        
        mort_c <- fread(paste0(cause_dir, indicator, '_c_admin_', admin_lvl, '_raked_mortality_summary.csv')) %>%
          dplyr::select(-region, -cirange) %>%
          rename(mort_count_mean = mean, mort_count_upper = upper, mort_count_lower = lower)
        
        inc_c <- fread(paste0(cause_dir, indicator, '_c_admin_', admin_lvl, '_raked_incidence_summary.csv')) %>%
          dplyr::select(-region, -cirange) %>%
          rename(inc_count_mean = mean, inc_count_upper = upper, inc_count_lower = lower)
        
      }
      
      if (cause == 'diarrhea'){
        
        inc <- fread(paste0(cause_dir, 'had_diarrhea_incidence_estimate_ad', admin_lvl, '.csv')) %>%
          dplyr::select(c(paste0('ADM', admin_lvl, '_CODE'),'year','mean','upper','lower','pop')) %>%
          rename(incidence_rate_mean = mean, incidence_rate_upper = upper, incidence_rate_lower = lower)
        
        mort <- fread(paste0(cause_dir, 'had_diarrhea_deaths_estimate_ad', admin_lvl, '.csv')) %>%
          dplyr::select(c(paste0('ADM', admin_lvl, '_CODE'),'year','mean','upper','lower')) %>%
          rename(mortality_rate_mean = mean, mortality_rate_upper = upper, mortality_rate_lower = lower)
        
        inc_c <- fread(paste0(cause_dir, 'had_diarrhea_incidence_counts_ad', admin_lvl, '.csv')) %>%
          dplyr::select(c(paste0('ADM', admin_lvl, '_CODE'),'year','mean','upper','lower')) %>%
          rename(incidence_count_mean = mean, incidence_count_upper = upper, incidence_count_lower = lower)
        
        mort_c <- fread(paste0(cause_dir, 'had_diarrhea_deaths_counts_ad', admin_lvl, '.csv')) %>%
          dplyr::select(c(paste0('ADM', admin_lvl, '_CODE'),'year','mean','upper','lower')) %>%
          rename(mortality_count_mean = mean, mortality_count_upper = upper, mortality_count_lower = lower)

      }
      
      #merge onto dalys
      loc_info <- grep(names(inc), pattern = 'ADM', value = T)
      
      admin_dalys <- merge(admin_dalys, inc, by = c('year', loc_info)) %>%
        merge(mort, by = c('year', loc_info)) %>%
        merge(inc_c, by = c('year', loc_info)) %>%
        merge(mort_c, by = c('year', loc_info)) %>%
        filter(!(ADM0_NAME %in% exclude_countries)) %>%
        as.data.table()
      
      #add indicator column
      admin_dalys$indicator <- cause
      
      
    } else if (cause == 'malaria'){
      load('<<< FILEPATH REDACTED >>>')
      pop <- get(paste0('admin_', admin_lvl)) %>%
        select(c('year',paste0('ADM', admin_lvl, '_CODE'), 'pop'))
      
      ylds <- fread(paste0(cause_dir, 'had_malaria_admin_', admin_lvl, '_raked_yld_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylds_rate_mean = mean, ylds_rate_upper = upper, ylds_rate_lower = lower)
      
      ylls <- fread(paste0(cause_dir, 'had_malaria_admin_', admin_lvl, '_raked_yll_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylls_rate_mean = mean, ylls_rate_upper = upper, ylls_rate_lower = lower)
      
      dalys <- fread(paste0(cause_dir, 'had_malaria_admin_', admin_lvl, '_raked_daly_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(dalys_rate_mean = mean, dalys_rate_upper = upper, dalys_rate_lower = lower)
      
      mort <- fread(paste0(cause_dir, 'had_malaria_admin_', admin_lvl, '_raked_mortality_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(mort_rate_mean = mean, mort_rate_upper = upper, mort_rate_lower = lower)
      
      inc <- fread(paste0(cause_dir, 'had_malaria_admin_', admin_lvl, '_raked_incidence_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(inc_rate_mean = mean, inc_rate_upper = upper, inc_rate_lower = lower)
      
      ylds_c <- fread(paste0(cause_dir, 'had_malaria_c_admin_', admin_lvl, '_raked_yld_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylds_count_mean = mean, ylds_count_upper = upper, ylds_count_lower = lower)
      
      ylls_c <- fread(paste0(cause_dir, 'had_malaria_c_admin_', admin_lvl, '_raked_yll_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(ylls_count_mean = mean, ylls_count_upper = upper, ylls_count_lower = lower)
      
      dalys_c <- fread(paste0(cause_dir, 'had_malaria_c_admin_', admin_lvl, '_raked_daly_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(dalys_count_mean = mean, dalys_count_upper = upper, dalys_count_lower = lower)
      
      mort_c <- fread(paste0(cause_dir, 'had_malaria_c_admin_', admin_lvl, '_raked_mortality_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(mort_count_mean = mean, mort_count_upper = upper, mort_count_lower = lower)
      
      inc_c <- fread(paste0(cause_dir, 'had_malaria_c_admin_', admin_lvl, '_raked_incidence_summary.csv')) %>%
        dplyr::select(-region, -cirange) %>%
        rename(inc_count_mean = mean, inc_count_upper = upper, inc_count_lower = lower)
      
      #merge onto dalys
      loc_info <- grep(names(dalys), pattern = 'ADM', value = T)
      
      admin_dalys <- merge(ylds, ylls, by = c('year', loc_info)) %>%
        merge(dalys, by = c('year', loc_info)) %>%
        merge(mort, by = c('year', loc_info)) %>%
        merge(inc, by = c('year', loc_info)) %>%
        merge(ylds_c, by = c('year', loc_info)) %>%
        merge(ylls_c, by = c('year', loc_info)) %>%
        merge(dalys_c, by = c('year', loc_info)) %>%
        merge(inc_c, by = c('year', loc_info)) %>%
        merge(mort_c, by = c('year', loc_info)) %>%
        merge(pop, by = c('year', paste0('ADM', admin_lvl, '_CODE'))) %>%
        filter(!(ADM0_NAME %in% exclude_countries)) %>%
        as.data.table()
      
      admin_dalys$indicator <- 'malaria'
      
    }
    
    #(6) Reshape all csvs long  -------------------------------------
    loc_info <- grep(names(admin_dalys), pattern = 'ADM', value = T)
    
    measure_names <- names(admin_dalys)
    measure_names <- measure_names[!measure_names %in% c('year', loc_info,'pop','indicator')]
    admin_dalys <- melt(admin_dalys, id.vars = c('year', loc_info, 'pop','indicator'), measure.vars = measure_names)
    
    admin_dalys[, c('measure', 'metric', 'stat') := tstrsplit(variable, "_", fixed = TRUE)]
    admin_dalys <- dplyr::select(admin_dalys, -variable)
    save <- admin_dalys
    
    admin_dalys <- save
    
    if (admin_lvl == 0) admin_dalys <- dcast(admin_dalys, year + ADM0_CODE + pop + indicator + ADM0_NAME + measure ~ stat + metric, value.var = 'value')
    if (admin_lvl == 1) admin_dalys <- dcast(admin_dalys, year + ADM0_CODE + pop + indicator + ADM0_NAME + measure + ADM1_NAME + ADM1_CODE ~ stat + metric, value.var = 'value')
    if (admin_lvl == 2) admin_dalys <- dcast(admin_dalys, year + ADM0_CODE + pop + indicator + ADM0_NAME + measure + ADM1_NAME + ADM1_CODE + ADM2_CODE + ADM2_NAME ~ stat + metric, value.var = 'value')
    
    #(5) Write csv to outputdir --------------------------------------
    outdir <- '<<< FILEPATH REDACTED >>>'
    dir.create(outdir)
    write.csv(admin_dalys, paste0(outdir, cause, '_dalys_admin', admin_lvl, '.csv'))
  }
}

#(7) Check all csvs  --------------------------------------
if (check_csvs){
  ind_list <- c('diarrhea', 'lri','malaria')
  admin_list <- c(0:2)
  
  for (ind in ind_list){
    for (admin in admin_list){
      results <- fread('<<< FILEPATH REDACTED >>>')
      countries <- length(unique(results$ADM0_NAME))
      
      print(paste0(ind, ' ADM', admin, ' has ', countries, ' unique countries.'))
      
    }
  }
}

ind <- 'diarrhea'
admin <- 2
results <- fread('<<< FILEPATH REDACTED >>>')
message('Countries with NA values for any value:')
unique(filter(results, is.na(lower_count) | is.na (upper_count) | is.na(mean_count) |
                is.na(lower_rate) | is.na(upper_rate) | is.na(mean_rate))$ADM0_NAME)

missing <- filter(results, is.na(lower_count) | is.na (upper_count) | is.na(mean_count) |
                is.na(lower_rate) | is.na(upper_rate) | is.na(mean_rate))

message('Population distribution for these NA-valued admins:')
summary(missing$pop)

message('NA admins with pop != 0:')
disp <- filter(missing, pop!=0)
disp <- disp[,names(disp) %like% "ADM"] %>%
  unique.data.frame()
disp

message('Countries with value == 0 for any value:')
unique(filter(results, lower_count == 0 |upper_count == 0 | mean_count == 0 |
                lower_rate == 0 | upper_rate == 0  | mean_rate == 0)$ADM0_NAME)

message('0 admins with pop != 0:')
disp <- filter(results, lower_count == 0 |upper_count == 0 | mean_count == 0 |
         lower_rate == 0 | upper_rate == 0  | mean_rate == 0)
disp <- filter(missing, pop!=0)
disp <- disp[,names(disp) %like% "ADM"] %>%
  unique.data.frame()
disp

message('Table of admins with value == 0 for any value, by country:')
table(filter(results, lower_count == 0 |upper_count == 0 | mean_count == 0 |
         lower_rate == 0 | upper_rate == 0  | mean_rate == 0)$ADM0_NAME)

if (admin == 0) cast <- dcast(results, year + ADM0_CODE + pop + indicator + ADM0_NAME ~ measure, value.var = c('lower_count', 'lower_rate', 'upper_count','upper_rate','mean_count', 'mean_rate'))
if (admin == 1) cast <- dcast(results, year + ADM0_CODE + pop + indicator + ADM0_NAME + ADM1_NAME + ADM1_CODE ~ measure, value.var = c('lower_count', 'lower_rate', 'upper_count','upper_rate','mean_count', 'mean_rate'))
if (admin == 2) cast <- dcast(results, year + ADM0_CODE + pop + indicator + ADM0_NAME + ADM0_NAME + ADM1_NAME + ADM1_CODE + ADM2_NAME + ADM2_CODE ~ measure, value.var = c('lower_count', 'lower_rate', 'upper_count','upper_rate','mean_count', 'mean_rate'))

cast[,`:=` (yll_yld_mean_rate = mean_rate_ylls + mean_rate_ylds,
            yll_yld_upper_rate = upper_rate_ylls + upper_rate_ylds,
            yll_yld_lower_rate = lower_rate_ylls + lower_rate_ylds,
            yll_yld_mean_count = mean_count_ylls + mean_count_ylds,
            yll_yld_upper_count = upper_count_ylls + upper_count_ylds,
            yll_yld_lower_count = lower_count_ylls + lower_count_ylds)]

plot <- ggplot(data = cast, aes(x = yll_yld_mean_rate, y = mean_rate_dalys)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red') + ggtitle(paste0('ADM', admin, ' mean_rate'))
plot(plot)

plot <- ggplot(data = cast, aes(x = yll_yld_lower_rate, y = lower_rate_dalys)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red') + ggtitle(paste0('ADM', admin, ' lower_rate'))
plot(plot)

plot <- ggplot(data = cast, aes(x = yll_yld_upper_rate, y = upper_rate_dalys)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red') + ggtitle(paste0('ADM', admin, ' upper_rate'))
plot(plot)

plot <- ggplot(data = cast, aes(x = yll_yld_mean_count, y = mean_count_dalys)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red') + ggtitle(paste0('ADM', admin, ' mean_count'))
plot(plot)

plot <- ggplot(data = cast, aes(x = yll_yld_lower_count, y = lower_count_dalys)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red') + ggtitle(paste0('ADM', admin, ' lower_count'))
plot(plot)

plot <- ggplot(data = cast, aes(x = yll_yld_upper_count, y = upper_count_dalys)) + geom_point() + geom_abline(slope = 1, intercept = 0, color = 'red') + ggtitle(paste0('ADM', admin, ' upper_count'))
plot(plot)
