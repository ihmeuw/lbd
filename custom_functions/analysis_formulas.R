# ---------------------------------------------------------------------------------------------
# Functions useful for exploratory analysis of ORS, RHF, ORT, and diarrhea
#
# Written by Kirsten Wiens
# 2019-01-25
# ---------------------------------------------------------------------------------------------


# -------------------------------------------------------------------
# Read in and clean data

get_data <- function(indi, run_date, indicator_group = 'ort', ad = 'admin_0') {
  
  # load data
  mydat <- fread(paste0('<<<< FILEPATH REDACTED >>>>', 
                        indi, '_', ad, '_unraked_summary.csv'))
  
  # clean data
  mydat <- mydat[!is.na(mean)]
  mydat[, cirange := NULL]
  mydat[, indicator := indi]
  
  # end function
  return(mydat)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Calculate fold difference between max and min admin unit

min_max_diff <- function(x) {
  
  # calculate fold difference
  diff <- max(x)/min(x)
  
  # end function
  return(diff)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Calculate relative inequity

rel_dev <- function(x, y) {
  
  # calculate fold difference
  ineq <- (x - y) / y
  
  # end function
  return(ineq)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Calculate absolute inequity

abs_dev <- function(x, y) {
  
  # calculate fold difference
  ineq <- x - y
  
  # end function
  return(ineq)
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------
# Correlation coefficient function

# get mean spearman statistic
get_cor_coef <- function(x, y) {
  if (!is.na(mean(x)) | !is.na(mean(y))) { 
    return(cor.test(x, y, method = 'spearman')[['estimate']][['rho']])
  } else {
    return(NA)
  }
}
# -------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# General deaths averted for comparative risk assessment

# get deaths averted
get_deaths_averted <- function(observed_rate, counterfactual_rate, observed_deaths, rr) {
  
  paf_obs <- (observed_rate*(rr - 1))/(observed_rate*(rr - 1) + 1)
  
  paf_ctf <- (counterfactual_rate*(rr - 1))/(counterfactual_rate*(rr - 1) + 1)
  
  averted <- observed_deaths * (((1 - paf_obs)/(1 - paf_ctf)) * paf_ctf - paf_obs)
  
  return(averted)
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Get PAF

get_paf <- function(rate, rr) {
  
  paf <- (rate*(rr - 1))/(rate*(rr - 1) + 1)
  
  return(paf)
  
}

# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Deaths averted across four indicators

# get total deaths averted
get_total_deaths_averted <- function(observed_deaths, 
                                     paf_O1, paf_O2, paf_O3, paf_O4,
                                     paf_C1, paf_C2, paf_C3, paf_C4,
                                     only_keep_pos = FALSE) {
  
  paf_obs <- 1 - (1 - paf_O1) * (1 - paf_O2) * (1 - paf_O3) * (1 - paf_O4)
  
  paf_ctf <- 1 - (1 - paf_C1) * (1 - paf_C2) * (1 - paf_C3) * (1 - paf_C4)
  
  averted <- observed_deaths * paf_ctf - observed_deaths * paf_obs
  
  if (only_keep_pos) averted <- max(0, averted)
  
  return(averted)
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# Deaths averted across two indicators

# get total deaths averted
get_half_deaths_averted <- function(observed_deaths, 
                                    paf_O1, paf_O2,
                                    paf_C1, paf_C2,
                                    only_keep_pos = FALSE) {
  
  paf_obs <- 1 - (1 - paf_O1) * (1 - paf_O2)
  
  paf_ctf <- 1 - (1 - paf_C1) * (1 - paf_C2)
  
  averted <- observed_deaths * paf_ctf - observed_deaths * paf_obs
  
  if (only_keep_pos) averted <- max(0, averted)
  
  return(averted)
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# General number of deaths averted by to a given risk factor per 1,000 (rate)

# get ratio averted
get_rate_averted <- function(deaths_averted, population) {
  
  rate <- deaths_averted/population*1000
  
  return(rate)
  
}
# -------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------
# General ratio of change in deaths attributable to a given risk factor

# get ratio averted
get_ratio_averted <- function(deaths_averted, change_in_deaths) {
  
  ratio <- deaths_averted/change_in_deaths
  
  return(ratio)
  
}
# -------------------------------------------------------------------------------------------