# --------------------------------------------------------------------------------------------------
# squeeze_outputs()
#
# Function that squeezes rhf_only and any_ors estimates to sum to 1 - no_ort
#
# Inputs:
# run_date - run date for current model
# shapefile_version - shapefile version corresponding to run date
# repo - parent modeling repo
# holdout - which holdout we're running this for (0 is the full model)
#
# Outputs (saved in main repo and /pred_derivatives/admin_summaries/, respectively):
# - Matrices containing squeezed draw-level results
# - CSVs containing mean, upper, and lower estimates
# --------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# Start function
squeeze_outputs <- function(run_date,
                            shapefile_version,
                            repo = '<<<< FILEPATH REDACTED >>>>',
                            holdout = 0) {
  # -------------------------------------------------------------------------


  # --------------------------------------------------------------------------------------------------------
  # Setup
  
  # set indicator group and indicators
  indicator_group <-  'ort'
  indicators <- c('no_ort', 'rhf_only', 'any_ors')
  
  # Load MBG packages and functions
  package_list <- c(t(read.csv('<<<< FILEPATH REDACTED >>>>/package_list.csv', header=FALSE)))
  source(paste0(repo, '/mbg_central/setup.R'))
  mbg_setup(package_list = package_list, repos = repo)
  
  # Load custom functions
  source(paste0(repo, '/custom_functions/format_draw_objects.R'))
  
  # Define share directories
  share_dirs <- paste0('<<<< FILEPATH REDACTED >>>>')
  names(share_dirs) <- indicators
  # --------------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------------------------------------------
  # Load and format outputs
  
  # format admin draws
  dt <- format_admin_results(ind_gp = rep(indicator_group, 3),
                             ind = indicators,
                             rd = rep(run_date, 3),
                             measure = rep('', 3),
                             suffix = rep('_eb_bin0_0', 3),
                             rk = rep(FALSE, 3))
  
  # format hierarchy list
  sp_list <- list()
  sp_list[[1]] <- unique(sp_hierarchy_list[, grep('ADM0', names(sp_hierarchy_list)), with = FALSE])
  sp_list[[2]] <- unique(sp_hierarchy_list[, grep('ADM1|ADM0', names(sp_hierarchy_list)), with = FALSE])
  sp_list[[3]] <- unique(sp_hierarchy_list[, grep('ADM0|ADM1|ADM2', names(sp_hierarchy_list)), with = FALSE])
  
  # add region name to dt
  dt[substring(code, nchar(code)-2) == 146, region := 'MLI']
  dt[substring(code, nchar(code)-2) == 194, region := 'SEN']
  dt[substring(code, nchar(code)-2) == 200, region := 'SLE']
  
  # correct the years for Mali, if Mali is being modeled
  dt$year <- as.numeric(dt$year)
  dt[region == 'MLI', year := year+1, by = .I]
  # ---------------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------------------------------------------
  # Squeeze               
  
  # squeeze
  dt[, any_ors_sqz := any_ors/(any_ors + rhf_only) * (1 - no_ort), by = .I] 
  dt[, rhf_only_sqz := rhf_only/(any_ors + rhf_only) * (1 - no_ort), by = .I]
  dt[, no_ort_sqz := 1 - (any_ors_sqz + rhf_only_sqz)]
  
  # save a copy
  for (i in indicators) save(dt, file = paste0(share_dirs[[i]], 'formatted_combined_admin_results.RData'))
  # ---------------------------------------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------------------------
  # Save formatted squeezed results
  
  for (a in 0:2) {
    message(a)
  
    for (i in indicators) {
      message(i)
    
      # subset by admin
      dt_ad <- dt[agg_level == paste0('ADM', a)]
      
      # get mean, upper, and lower
      dt_ad[, mean := lapply(.SD, mean, na.rm = T), 
            .SDcols = paste0(i, '_sqz'), by = c('code', 'year')]
      dt_ad[, upper := lapply(.SD, upper), 
            .SDcols = paste0(i, '_sqz'), by = c('code', 'year')]
      dt_ad[, lower := lapply(.SD, lower), 
            .SDcols = paste0(i, '_sqz'), by = c('code', 'year')]
      
      # save summary and free up space
      ad_summary <- merge(sp_list[[a+1]], 
                          unique(dt_ad[, c('code', 'year', 'mean', 'upper', 'lower')]), 
                          by.x = paste0('ADM', a, '_CODE'), by.y = 'code',
                          allow.cartesian = T)
      write.csv(ad_summary, paste0(share_dirs[[i]], 'pred_derivatives/admin_summaries/',
                                   i, '_admin_', a, '_squeezed_summary.csv'))
    
      # reshape draw object wide
      dt_ad <- dcast(dt_ad, year + code ~ draw, value.var = paste0(i, '_sqz'))
      setnames(dt_ad, 'code', paste0('ADM', a, '_CODE'))
      
      # check for NAs
      message('TESTING: Percent of NA rows per column is: ', mean(is.na(dt_ad[, V1])), '%')
      
      # save draw objects
      outpath <- paste0(share_dirs[[i]], i, '_squeezed_', 
                        'admin_', a, '_draws_eb_bin0_', holdout, '.RData')
      save(dt_ad, file = outpath)
    
    }
  
  }
  # --------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------
  # End function
  return('Squeezing complete. Results saved.')
  
}
# -----------------------------------------------

