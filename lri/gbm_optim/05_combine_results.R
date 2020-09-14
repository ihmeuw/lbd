##############################################################################
## Combine results from BRT optimization
##############################################################################


## Setup -------------------------------------------------------------------------

# clear workspace
rm(list=ls())
library(data.table)

# indicator tested
indicators <- c('lri')

# space bounds version
experiment_version <- '<<<< FILEPATH REDACTED >>>>'

# set directories
dir_train <- '<<<< FILEPATH REDACTED >>>>'
dir_test <- '<<<< FILEPATH REDACTED >>>>'
dir_model <- '<<<< FILEPATH REDACTED >>>>'
dir_params <- '<<<< FILEPATH REDACTED >>>>'

for (indicator in indicators) {
  message(indicator)
  
  # get completed regions
  regions <- list.files(dir_params, pattern = paste0('best_pars_', indicator))
  regions <- gsub(paste0('best_pars_', indicator, '_'), '', regions)
  regions <- gsub('_exp1.csv', '', regions)
  if (indicator == 'lri') regions <- gsub('or_rhf_', '', regions)
  regions <- unique(regions)
  
  
  ## Combine results from BRT optimization runs -------------------------------------------------------------------------
  
  # function to read in, clean, and save model outputs
  combine_outputs <- function(dir, save = T, regs = regions, ind = indicator) {
    
    # initiate a list
    all_results <- as.data.table(lapply(paste0(dir, list.files(dir, pattern = paste0(ind, '_', regs[[1]]))), fread))
    all_results[, region := NA]
    all_results <- all_results[-1, ]
    
    for (r in regs) {
      message(r)
    
      # read files into a list
      results <- as.data.table(lapply(paste0(dir, list.files(dir, pattern = paste0(ind, '_', r, '_exp'))), fread))
      results[, region := r]
      
      # rbind into one file and order by experiment
      all_results <- rbind(all_results, results, fill = T)
    }
      
    # save, if desired
    if (save == T) write.csv(all_results, paste0(dir, 'summary_', indicator, '.csv'))
  }
  
  # combine parameter files
  combine_outputs(dir_params)
  
  # combine in-sample stats
  combine_outputs(dir_train)
  
  # combine out-of-sample stats
  combine_outputs(dir_test)
  
  # combine overall model stats
  combine_outputs(dir_model)

}