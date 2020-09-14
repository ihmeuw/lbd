# ---------------------------------------------------------------------------------------------
# get_gbd_group_sf()
#
# Function that pulls gbd results by cause, measure, age, and sex group and returns a data 
# frame with scaling factors (sf) by each group
#
# Inputs:
#   region - region or or vector of regions to pull results for
#   modeling_shapefile_version - shapefile version used in your current MBG model
#   year_id - integer vector of years to pull estimates for
#   cause_id - GBD cause ID for your indicator
#   age_group_id - GBD age group ids that you want to create groups for
#   sex_id - GBD sex ids that you want to create groups for
#   measure_id - GBD measure you want to pull (e.g. prevalence, incidence, mortality)
#   gbd_round_id - GBD round that you want to pull
# 
#
# Outputs:
# Data table with the following columns:
#   year - modeling year
#   group - group identifier named a#_s#, with the first number being the GBD age_group_id
#           and the second number being the GBD sex_id
#   ADM0_CODE - LBD admin 0 code
#   scaling_factor - GBD scaling factor by location and year that can be multiplied by 
#                    LBD rate estimates to get groups by age and/or sex
# ---------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------
# Start function
get_gbd_group_sf <- function(region,
                             modeling_shapefile_version,
                             year_id = 2000:2019,
                             cause_id = 302,
                             age_group_id = c(2,3,4,5),
                             sex_id = c(1,2),
                             measure_id = 5,
                             gbd_round_id = 6) {
  # -----------------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------------------
  # Set-up
  
  # Report which type of groups we're creating
  message('Performing age/sex splits for:')
  message('cause_id = ', cause_id)
  message('age_group_id = ', paste0(age_group_id, ' '))
  message('sex_id = ', paste0(sex_id, ' '))
  message('This will generate ', length(age_group_id)*length(sex_id), ' total groups.\n')
  
  # GBD get_outputs function
  source('<<< FILEPATH REDACTED >>>/get_outputs.R')
  # ----------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------------------------------
  # Pull and clean GBD estimates
  
  # get GBD age and sex split diarrhea data
  message('Pulling GBD estimates')
  results <- get_outputs(topic = 'cause', 
                         cause_id = cause_id,
                         location_id = 'all', 
                         year_id = year_id,
                         age_group_id = age_group_id,
                         sex_id = sex_id,
                         metric_id = c(1,3),
                         measure_id = measure_id,
                         gbd_round_id = gbd_round_id,
                         decomp_step='step4',
                         version='latest')
  
  # create grouping variables based on inputs
  group_vars <- c()
  if (length(age_group_id > 1)) {
    group_vars <- c(group_vars, 'age_group_id')
    results[, age_group_id := paste0('a', age_group_id)]
  }
  if (length(sex_id > 1)) {
    group_vars <- c(group_vars, 'sex_id')
    results[, sex_id := paste0('s', sex_id)]
  }
  
  # clean up and get GBD population
  results <- dcast(results, 
                   as.formula(paste(paste(c(group_vars, 'location_id', 'year_id'), collapse = '+'), ' ~ metric_name')), 
                   value.var = 'val')
  results[, pop_gbd := Number/Rate]
  setnames(results, c('Number', 'year_id'), c('val_gbd', 'year'))
  
  # create group identifier
  results[, group := Reduce(function(...) paste(..., sep = '_'), .SD), 
          .SDcols = group_vars]
  
  # remove unneeded columns
  results[, c(group_vars, 'Rate') := NULL]
  # ---------------------------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------------------
  # Calculate scaling factors
  
  # get total GBD prevalence and population for a given year/location
  message('Calculating scaling factors')
  results[, total_val := sum(val_gbd), by = c('location_id', 'year')]
  results[, total_pop := sum(pop_gbd), by = c('location_id', 'year')]
  
  # get ADM0_CODES by gbd location
  connector <-
    data.table(get_gbd_locs(rake_subnational = F ,
                            reg = "all",
                            shapefile_version = modeling_shapefile_version))
  setnames(connector, 'ADM_CODE', 'ADM0_CODE')
  results <- merge(results, connector, by = 'location_id')
  
  # get scaling factor by location and clean up
  results[, scaling_factor := get_scaling_factor(val_gbd, total_val, pop_gbd, total_pop), 
          by = c('location_id', 'year', 'group')]
  results[, c('location_id', 'val_gbd', 'total_val', 'pop_gbd', 'total_pop') := NULL]
  # --------------------------------------------------------------------------------------
  
  
  # -------------------------------------------------------------------------------------
  # End function
  return(results)
  
}
# -----------------------------------------------------------------------------------




# ---------------------------------------------------------------------------------------------
# get_scaling_factor()
#
# Function that calculates scaling factor by GBD age and/or sex group
#
# Inputs:
#   cases_by_group - number of cases by GBD group, location, and year
#   cases_total - number of cases by GBD location, and year
#   population_by_group - population by GBD group, location, and year
#   population_total - population by GBD location, and year
#
# Outputs:
#   scaling_factor - constant to multiple pixel or admin rates by to get age splits
# ---------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------
# Start function
get_scaling_factor <- function(cases_by_group, cases_total,
                               population_by_group, population_total) {
  
  # Get scaling factor
  scaling_factor <- (cases_by_group/cases_total) * (population_total/population_by_group)
  
  # End function
  return(scaling_factor)
  
}
# -------------------------------------------------------------------------




# ---------------------------------------------------------------------------------------------
# get_gbd_measure_id()
#
# Function that determines GBD measure ID based on MBG measure
#
# Inputs:
#   measure - measure you'd like to pull the ID for (currently not compatible with etiologies)
#
# Outputs:
#   measure_id - GBD measure ID that can be input into get_outputs shared function
# ---------------------------------------------------------------------------------------------

# ----------------------------------------
# Start function
get_gbd_measure_id <- function(measure) {
  
  # Define measures
  m <- list(
    'deaths' = 1,
    'mortality' = 1,
    'daly' = 2,
    'yld' = 3,
    'yll' = 4,
    'prevalence' = 5,
    'incidence' = 6
  )
  
  # Get id based on measure
  id <- m[[measure]]
  
  # End function
  return(id)
  
}
# ----------------------------------------