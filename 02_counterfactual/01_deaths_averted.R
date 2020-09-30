##############################################################################
# Calculated PAFs and deaths averted for LRI, diarrhea due to stunting/wasting
# and water/sanitation
############################################################################

# (1) Setup ---------------------------------------------------------

rm(list = ls())
library(stringr)

cause <- 'lri' #lri or dia
risk <- 'child_growth_failure' #wash or child_growth_failure
admin_levels <- c(0:2)

custom_wash_filepath <- TRUE

#user arguments
lri_run_date <- '2019_09_16_16_39_00'
stunt_wast_run_date <- '2019_09_30_11_45_54'
underweight_run_date <- '2019_09_30_11_45_54'
wash_run_date <- '2019_10_31' #'2019_05_20_00_00_02' 
dia_run_date <- '2019_09_17_14_12_53'

run_date <- '2019_10_28' #for tridalys results

if (risk == 'child_growth_failure'){
  group_list <- c('wasting', 'stunting','underweight')
} else if (risk == 'wash'){
  group_list <- c('water', 'sani')
}

#specify RR vectors
rr_lri_stunting <- c(1.125,1.318,2.355) # mild, mod, sev  https://www.thelancet.com/cms/10.1016/S0140-6736(17)32366-8/attachment/090de11b-0f50-4019-bb51-e10a2c782e35/mmc1.pdf
rr_lri_wasting <- c(5.941,20.455,47.67) # mild, mod, sev  https://www.thelancet.com/cms/10.1016/S0140-6736(17)32366-8/attachment/090de11b-0f50-4019-bb51-e10a2c782e35/mmc1.pdf
rr_lri_underweight <- c(1.145,1.365,2.593) #mild, mod, sev https://www.thelancet.com/cms/10.1016/S0140-6736(17)32366-8/attachment/090de11b-0f50-4019-bb51-e10a2c782e35/mmc1.pdf

rr_dia_sani <- c(1,1.249) # imp, unimp
rr_dia_water <- c(1,1.118,1.364) # piped, imp, unimp
rr_dia_stunting <- c(1.111,1.222,1.851) # mild, mod, sev https://www.thelancet.com/cms/10.1016/S0140-6736(17)32366-8/attachment/090de11b-0f50-4019-bb51-e10a2c782e35/mmc1.pdf
rr_dia_wasting <- c(6.601,23.261,105.759) # mild, mod, sev https://www.thelancet.com/cms/10.1016/S0140-6736(17)32366-8/attachment/090de11b-0f50-4019-bb51-e10a2c782e35/mmc1.pdf
rr_dia_underweight <- c(1.088,1.23,2.332) #mild, mod, sev https://www.thelancet.com/cms/10.1016/S0140-6736(17)32366-8/attachment/090de11b-0f50-4019-bb51-e10a2c782e35/mmc1.pdf

libs <- c('data.table', 'raster', 'dplyr', 'ggplot2', 'sf', 'viridis')
lapply(libs, library, character.only = TRUE)

#define function to load in the admin draws for wash and cgf
load_dat <- function(indi, rd, risk, admin_level) {
  setwd('<<< FILEPATH REDACTED >>>')
  
  if (admin_level == 2) admin_code <- 'ADM2_CODE'
  if (admin_level == 1) admin_code <- 'ADM1_CODE'
  if (admin_level == 0) admin_code <- 'ADM0_CODE'
  
  if (risk == 'wash'){
    rake_string <- '_unraked'
    mydat <- do.call(rbind, lapply(list.files(pattern = 
                                                paste0(indi, rake_string, '_admin_draws_eb_bin0_.*_0.RData$')), function(x) {
                                                  load(x)
                                                  if (admin_level == 2) return(admin_2)
                                                  if (admin_level == 1) return(admin_1)
                                                  if (admin_level == 0) return(admin_0)
                                                }))
  } else if (risk == 'child_growth_failure'){
    rake_string <- '_raked'
    mydat <- do.call(rbind, lapply(list.files(pattern = 
                                                paste0(indi, rake_string, '_admin_draws_eb_bin0_0.RData')), function(x) {
                                                  load(x)
                                                  if (admin_level == 2) return(admin_2)
                                                  if (admin_level == 1) return(admin_1)
                                                  if (admin_level == 0) return(admin_0)
                                                }))
  }
  
  mydat <- as.data.frame(mydat)
  draws <- mydat[, grep('V', names(mydat))]
  mydat$mean <- apply(draws, 1, mean, na.rm = TRUE)
  mydat <- dplyr::select(mydat, mean, pop, admin_code, year)
  return(mydat)
}

for (admin_level in admin_levels){
  
  #set identifiers based on admin level
  if (admin_level == 2){
    admin_code <- 'ADM2_CODE'
    admin_obj <- 'admin_2'
  }
  
  if (admin_level == 1) {
    admin_code <- 'ADM1_CODE'
    admin_obj <- 'admin_1'
  }
  
  if (admin_level == 0) {
    admin_code <- 'ADM0_CODE'
    admin_obj <- 'admin_0'
  }
  
  # (2) Read in deaths data  ---------------------------------------------------------
  
  #read in cause data and merge on aggregated population values
  if (cause == 'lri'){
    cause_data_lbd <- fread('<<< FILEPATH REDACTED >>>')
    load('<<< FILEPATH REDACTED >>>')
    admin_pop <- dplyr::select(get(admin_obj), pop, admin_code, year)
    rm(admin_2, admin_1, admin_0)
    cause_data_lbd <- merge(cause_data_lbd, admin_pop, by = c(admin_code, 'year'))
    rm(admin_pop)
  } else if (cause == 'dia'){
    cause_data_lbd <- fread('<<< FILEPATH REDACTED >>>')
  }
  
  #read in risk data 
  for (group in group_list){
    if (group == 'stunting'){
      indis <- c('stunting_mil_c', 'stunting_mod_c', 'stunting_sev_c')
      cgf_run_date <- stunt_wast_run_date
    } else if (group == 'wasting'){
      indis <- c('wasting_mil_c', 'wasting_mod_c', 'wasting_sev_c')
      cgf_run_date <- stunt_wast_run_date
    } else if (group == 'underweight'){
      indis <- c('underweight_mil_c', 'underweight_mod_c', 'underweight_sev_c')
      cgf_run_date <- underweight_run_date 
    } else if (group == 'sani') {
      indis <- c('s_piped', 's_imp_cr', 's_unimp_cr')
    } else if (group == 'water'){
      indis <- c('w_piped', 'w_imp_cr', 'w_unimp_cr')
    }
    
    #change indis if using custom wash file path
    if (custom_wash_filepath == TRUE & group == 'water'){
      indis <- c('w_piped', 'w_imp', 'w_unimp')
    }
    
    if (custom_wash_filepath == TRUE & group == 'sani'){
      indis <- c('s_piped', 's_imp', 's_unimp')
    }
    
    if (group %in% c('stunting', 'wasting','underweight')){
      # Load data for indicators
      mydat <- load_dat(indis[1], cgf_run_date, risk, admin_level)
      mild <- mydat %>% rename(mild = mean)
      
      mydat <- load_dat(indis[2], cgf_run_date, risk, admin_level)
      mod <- mydat %>% rename(mod = mean)
      
      mydat <- load_dat(indis[3], cgf_run_date, risk, admin_level)
      sev <- mydat %>% rename(sev = mean)
      
      # Combine indicator data across all levels
      tern_df <- left_join(mild, mod)
      tern_df <- left_join(tern_df, sev)
      tern_df <- tern_df[complete.cases(tern_df),]
      
      #make sure levels don't add to greater than 1
      #if sum exceeds one, normalize all levels to sum of levels to remedy
      tern_df <- as.data.table(tern_df)
      tern_df[, level_sum := mild + mod + sev]
      tern_df[level_sum > 1, `:=`(mild = mild/level_sum, mod = mod/level_sum, sev = sev/level_sum)]
      tern_df[, level_sum := mild + mod + sev]
      tern_df[, other := 1 - level_sum]
      tern_df$level_sum <- NULL
      
      # Subset to first and last year
      y00 <- filter(tern_df, year == 2000) %>%
        select(-pop, -year) %>%
        rename(mild00 = mild, mod00 = mod, sev00 = sev, other00 = other)
      y17 <- filter(tern_df, year == 2017) %>%
        select(-pop, -year) %>%
        rename(mild17 = mild, mod17 = mod, sev17 = sev, other17 = other)
      y_all <- left_join(y17, y00)
      
      # Add on cause data
      cause_deaths_17 <- cause_data_lbd %>%
        as.data.frame() %>%
        dplyr::filter(year == 2017) %>%
        dplyr::select(mean, admin_code) %>%
        rename(cause_deaths_17 = mean)
      cause_deaths_00 <- cause_data_lbd %>%
        as.data.frame() %>%
        dplyr::filter(year == 2000) %>%
        dplyr::select(mean, admin_code) %>%
        rename(cause_deaths_00 = mean)
      cause_deaths <- left_join(cause_deaths_17, cause_deaths_00)
      y_all <- left_join(y_all, cause_deaths)
      y_all <- y_all[complete.cases(y_all),]
      
      # calculate deaths averted
      if (cause == 'lri' & group == 'wasting'){
        rr <- rr_lri_wasting
      } else if (cause == 'dia' & group == 'wasting'){
        rr <- rr_dia_wasting
      } else if (cause == 'lri' & group == 'stunting'){
        rr <- rr_lri_stunting
      } else if (cause == 'dia' & group == 'stunting'){
        rr <- rr_dia_stunting
      } else if (cause == 'lri' & group == 'underweight'){
        rr <- rr_lri_underweight
      } else if (cause == 'dia' & group == 'underweight'){
        rr <- rr_dia_underweight
      }
      y_all <- mutate(y_all,
                      paf17 = (((mild17 * rr[1]) + (mod17*rr[2]) + (sev17*rr[3]) + other17) - 1)/
                        ((mild17 * rr[1]) + (mod17*rr[2]) + (sev17*rr[3]) + other17),
                      paf00 = (((mild00 * rr[1]) + (mod00*rr[2]) + (sev00*rr[3]) + other00) - 1)/
                        ((mild00 * rr[1]) + (mod00*rr[2]) + (sev00*rr[3]) + other00)) %>%
        mutate(deaths_averted = cause_deaths_17*((1-paf17)/(1-paf00)*paf00-paf17))
      
      save_df <- dplyr::select(y_all, admin_code, paf17, paf00, deaths_averted, cause_deaths_17, cause_deaths_00)
      setwd('<<< FILEPATH REDACTED >>>')
      
      write.csv(save_df, paste0(cause, '_', group, '_paf_ad', admin_level,'.csv'))
    } else if (group %in% c('water', 'sani')){
      
      # Load data for indicators
      
      #standard file path
      if (custom_wash_filepath == FALSE){
        mydat <- load_dat(indis[1], wash_run_date, risk, admin_level)
        piped <- mydat %>% rename(piped = mean)
        
        mydat <- load_dat(indis[2], wash_run_date, risk, admin_level)
        imp <- mydat %>% rename(imp = mean)
        
        mydat <- load_dat(indis[3], wash_run_date, risk, admin_level)
        unimp <- mydat %>% rename(unimp = mean)
        
        # Combine indicator data across all levels
        tern_df <- left_join(piped, imp)
        tern_df <- left_join(tern_df, unimp)
        tern_df <- tern_df[complete.cases(tern_df),] %>%
          select(-pop)
      } else if (custom_wash_filepath == TRUE){
        tern_df <- fread('<<< FILEPATH REDACTED >>>') %>%
          select(-V1) %>%
          rename_all(funs(str_replace(., indis[1], "piped"))) %>%
          rename_all(funs(str_replace(., indis[2], "imp"))) %>%
          rename_all(funs(str_replace(., indis[3], "unimp"))) %>%
          select(-agg_level) %>%
          plyr::rename(replace = c('code' = admin_code))
        }
      
      #make sure levels don't add to greater than 1
      tern_df <- tern_df %>%
        mutate(imp = (1 - piped) * imp) %>%
        mutate(unimp = 1 - piped - imp)
      
      # Subset to first and last year
      y00 <- filter(tern_df, year == 2000) %>%
        select(-year) %>%
        rename(piped00 = piped, imp00 = imp, unimp00 = unimp)
      y17 <- filter(tern_df, year == 2017) %>%
        select(-year) %>%
        rename(piped17 = piped, imp17 = imp, unimp17 = unimp)
      y_all <- left_join(y17, y00)
      
      # Add on cause data
      cause_deaths_17 <- cause_data_lbd %>%
        as.data.frame() %>%
        dplyr::filter(year == 2017) %>%
        dplyr::select(mean, admin_code) %>%
        rename(cause_deaths_17 = mean)
      cause_deaths_00 <- cause_data_lbd %>%
        as.data.frame() %>%
        dplyr::filter(year == 2000) %>%
        dplyr::select(mean, admin_code) %>%
        rename(cause_deaths_00 = mean)
      cause_deaths <- left_join(cause_deaths_17, cause_deaths_00)
      y_all <- left_join(y_all, cause_deaths)
      y_all <- y_all[complete.cases(y_all),]
      
      # calculate deaths averted
      if (group == 'sani') {
        if (cause == 'lri'){
          stop('no sanitation PAF for LRI')
        } else if (cause == 'dia'){
          rr <- rr_dia_sani
        }
        
        #transform data to accomodate fact that we only have "unimproved" and "improved" categories
        y_all <- mutate(y_all, imp00 = imp00 + piped00,
                        imp17 = imp17 + piped17) %>%
          dplyr::select(-piped00, -piped17)
        
        y_all <- mutate(y_all,
                        paf17 = (((imp17*rr[1]) + (unimp17*rr[2])) - 1)/
                          ((imp17*rr[1]) + (unimp17*rr[2])),
                        paf00 = (((imp00*rr[1]) + (unimp00*rr[2])) - 1)/
                          ((imp00*rr[1]) + (unimp00*rr[2]))) %>%
          mutate(deaths_averted = cause_deaths_17*((1-paf17)/(1-paf00)*paf00-paf17))
        
      } else if (group == 'water') {
        if (cause == 'lri'){
          stop('no water PAF for LRI')
          
        } else if (cause == 'dia'){
          rr <- rr_dia_water
        }
        y_all <- mutate(y_all,
                        paf17 = (((piped17 * rr[1]) + (imp17*rr[2]) + (unimp17*rr[3])) - 1)/
                          ((piped17 * rr[1]) + (imp17*rr[2]) + (unimp17*rr[3])),
                        paf00 = (((piped00 * rr[1]) + (imp00*rr[2]) + (unimp00*rr[3])) - 1)/
                          ((piped00 * rr[1]) + (imp00*rr[2]) + (unimp00*rr[3]))) %>%
          mutate(deaths_averted = cause_deaths_17*((1-paf17)/(1-paf00)*paf00-paf17))
      }
      
      save_df <- dplyr::select(y_all, admin_code, paf17, paf00, deaths_averted, cause_deaths_00, cause_deaths_17)
      setwd('<<< FILEPATH REDACTED >>>')
      
      write.csv(save_df, paste0(cause, '_', group, '_paf_ad', admin_level,'.csv'))
    }
  }
}

