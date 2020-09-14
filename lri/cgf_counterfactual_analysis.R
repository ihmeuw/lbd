#################################################################################################################################################################
# Calculated PAFs and deaths averted for LRI due to stunting, wasting, underweight
##################################################################################################################################################################

# (1) Setup -----------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls())

#user arguments
lri_run_date <- commandArgs()[4]
stunting_run_date <- commandArgs()[5]
wasting_run_date <- commandArgs()[6]
underweight_run_date <- commandArgs()[7]

#additional arguments for this risk
risk <- 'child_growth_failure'
group_list <- c('wasting', 'stunting','underweight')
admin_levels <- c(0:2)

#specify RR vectors
#citation: Supplement to GBD 2016 Risk Factors Collaborators. Global, regional, and national comparative risk assessment of 84 behavioural, environmental and 
#occupational,and metabolic risks or clusters of risks, 1990–2016: a systematic analysis for the Global Burden of Disease Study 2016. Lancet 2017; 390: 1345–422.

rr_lri_stunting <- c(1.125,1.318,2.355) # mild, mod, sev
rr_lri_wasting <- c(5.941,20.455,47.67) # mild, mod, sev
rr_lri_underweight <- c(1.145,1.365,2.593) #mild, mod, sev

#libraries and packages
libs <- c('data.table', 'raster', 'dplyr', 'ggplot2', 'sf', 'viridis')
lapply(libs, library, character.only = TRUE)

#define and create savedir
save_dir <- '<<<< FILEPATH REDACTED >>>>'
dir.create(save_dir)

#define function to load in the admin draws for cgf
load_dat <- function(indi, rd, risk, admin_level) {
  setwd('<<<< FILEPATH REDACTED >>>>')
  
  if (admin_level == 2) admin_code <- 'ADM2_CODE'
  if (admin_level == 1) admin_code <- 'ADM1_CODE'
  if (admin_level == 0) admin_code <- 'ADM0_CODE'
  
  rake_string <- '_raked'
  mydat <- do.call(rbind, lapply(list.files(pattern = 
                                              paste0(indi, rake_string, '_admin_draws_eb_bin0_0.RData')), function(x) {
                                                load(x)
                                                if (admin_level == 2) return(admin_2)
                                                if (admin_level == 1) return(admin_1)
                                                if (admin_level == 0) return(admin_0)
                                              }))
  
  mydat <- as.data.frame(mydat)
  draws <- mydat[, grep('V', names(mydat))]
  mydat$mean <- apply(draws, 1, mean, na.rm = TRUE)
  mydat <- dplyr::select(mydat, mean, pop, admin_code, year)
  return(mydat)
}

#(2) Loop over admin levels  -----------------------------------------------------------------------------------------------------------------------------------------------
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
  
  # (A) Read in deaths data  -----------------------------------------------------------------------------------------------------------------------------------------------
  
  #read in cause data and merge on aggregated population values
  cause_data_lbd <- fread('<<<< FILEPATH REDACTED >>>>')
  load('<<<< FILEPATH REDACTED >>>>')
  admin_pop <- dplyr::select(get(admin_obj), pop, admin_code, year)
  rm(admin_2, admin_1, admin_0)
  cause_data_lbd <- merge(cause_data_lbd, admin_pop, by = c(admin_code, 'year'))
  rm(admin_pop)
  
  # (B) Loop over groups in risk to calcuate PAFs and deaths averted ----------------------------------------------------------------------------------------------------------
  for (group in group_list){

    #define indicators, run date, and rr for group
    if (group == 'stunting'){
      indis <- c('stunting_mil_c', 'stunting_mod_c', 'stunting_sev_c')
      risk_run_date <- stunting_run_date
      rr <- rr_lri_stunting
    } else if (group == 'wasting'){
      indis <- c('wasting_mil_c', 'wasting_mod_c', 'wasting_sev_c')
      risk_run_date <- wasting_run_date
      rr <- rr_lri_wasting
    } else if (group == 'underweight'){
      indis <- c('underweight_mil_c', 'underweight_mod_c', 'underweight_sev_c')
      risk_run_date <- underweight_run_date
      rr <- rr_lri_underweight
    }
    
    # Load data for indicators
    mydat <- load_dat(indis[1], risk_run_date, risk, admin_level)
    mild <- mydat %>% rename(mild = mean)
    
    mydat <- load_dat(indis[2], risk_run_date, risk, admin_level)
    mod <- mydat %>% rename(mod = mean)
    
    mydat <- load_dat(indis[3], risk_run_date, risk, admin_level)
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
      dplyr::select(-pop, -year) %>%
      rename(mild00 = mild, mod00 = mod, sev00 = sev, other00 = other)
    y17 <- filter(tern_df, year == 2017) %>%
      dplyr::select(-pop, -year) %>%
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
    y_all <- mutate(y_all,
                    paf17 = (((mild17 * rr[1]) + (mod17*rr[2]) + (sev17*rr[3]) + other17) - 1)/
                      ((mild17 * rr[1]) + (mod17*rr[2]) + (sev17*rr[3]) + other17),
                    paf00 = (((mild00 * rr[1]) + (mod00*rr[2]) + (sev00*rr[3]) + other00) - 1)/
                      ((mild00 * rr[1]) + (mod00*rr[2]) + (sev00*rr[3]) + other00)) %>%
      mutate(deaths_averted = cause_deaths_17*((1-paf17)/(1-paf00)*paf00-paf17))
    
    save_df <- dplyr::select(y_all, admin_code, paf17, paf00, deaths_averted, cause_deaths_17, cause_deaths_00)
    
    write.csv(save_df, paste0(save_dir, 'has_lri_', group, '_paf_ad', admin_level,'.csv'))
  }
}

