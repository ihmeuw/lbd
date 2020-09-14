##############################################################################
## Combine results from BRT optimization
##############################################################################


## Setup -------------------------------------------------------------------------

# clear workspace
rm(list=ls())
library(data.table)

# region tested
region <- c('dia_se_asia')

# indicator tested
indicator <- c('had_diarrhea')

# space bounds version
experiment_version <- 'test18_2018-12-18'

# set directories
dir_train <- '<<<< FILEPATH REDACTED >>>>'
dir_test <- '<<<< FILEPATH REDACTED >>>>'
dir_model <- '<<<< FILEPATH REDACTED >>>>'
dir_params <- '<<<< FILEPATH REDACTED >>>>'


## Combine results from BRT optimization runs -------------------------------------------------------------------------

# function to read in, clean, and save model outputs
combine_outputs <- function(dir, order_col = 'exp_num', save = T, reg = region, ind = indicator) {
  
  # read files into a list
  results <- lapply(paste0(dir, list.files(dir, pattern = '*exp*')), fread)
  
  # rbind into one file and order by experiment
  dt <- setorderv(rbindlist(results), cols = order_col)
  
  # save, if desired
  if (save == T) write.csv(dt, paste0(dir, 'summary_', region, '_', indicator, '.csv'))
}


# combine in-sample stats
combine_outputs(dir_train)

# combine out-of-sample stats
combine_outputs(dir_test)

# combine overall model stats
combine_outputs(dir_model)

# combine parameter files
combine_outputs(dir_params, order_col = c('cv', 'bag.fraction'))
