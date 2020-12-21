# --------------------------------------------------------------------------------------------------
# analyze_treatment_changes()
#
# Function that calculate changes in treatment at the admin 0, 1, and 2 draw level and save outputs
#
# Inputs:
# run_date - run date for current model
# year_start - start date for analysis (start of modeling period, pre/post intervention, etc.)
# year_end - end date for analysis (end of modeling period, pre/post intervention, etc.)
# holdout - which holdout we're running this for (0 is the full model)
# country - iso3 for single country to run analysis for
#
# Outputs (saved in /pred_derivatives/treatment_change/):
# - Matrices of draw level differences in treatment between year_start and year_end
# - CSVs of summarized mean, upper, and lower values
# --------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Start function
analyze_treatment_changes <- function(run_date,
                                      year_start,
                                      year_end,
                                      country,
                                      holdout = 0) {
  # -----------------------------------------------------------------------------------


  # ---------------------------------------------------------------------------------------------------------
  # Setup
  
  # set indicator group and indicators
  indicator_group <-  'ort'
  indicators <- c('no_ort', 'rhf_only', 'any_ors')
  
  # define share directories
  share_dirs <- paste0('<<<< FILEPATH REDACTED >>>>')
  share_paths <- paste0(share_dirs, indicators, '_squeezed_admin_2_draws_eb_bin0_', holdout, '.RData')
  names(share_dirs) <- names(share_paths) <- indicators
  
  # create directories
  for (i in indicators) dir.create(paste0(share_dirs[[i]], 'pred_derivatives/treatment_change/'), showWarnings = F)
  # ---------------------------------------------------------------------------------------------------------
  
  
  # ----------------------------------------------------------------------------------------------------------
  # Calculate differences by draw
  
  # load hierarchy list
  load(paste0(share_dirs[[1]], indicators[[1]], '_unraked_admin_draws_eb_bin0_0.RData'))
  sp_list <- list()
  sp_list[[1]] <- unique(sp_hierarchy_list[, grep('ADM0', names(sp_hierarchy_list)), with = FALSE])
  sp_list[[2]] <- unique(sp_hierarchy_list[, grep('ADM1|ADM0', names(sp_hierarchy_list)), with = FALSE])
  sp_list[[3]] <- unique(sp_hierarchy_list[, grep('ADM0|ADM1|ADM2', names(sp_hierarchy_list)), with = FALSE])
  rm(admin_0, admin_1, admin_2, sp_hierarchy_list)
  
  # load formatted outputs
  load(paste0(share_dirs[[1]], 'formatted_combined_admin_results.RData'))
  
  # clean up columns
  dt[, (indicators) := NULL]
  setnames(dt, paste0(indicators, '_sqz'), indicators)
  
  # subset by years and region in analysis
  dt <- dt[region == country & (year == year_start | year == year_end)]
  
  # reshape
  dt <- melt(dt, id.vars = c('code', 'year', 'name', 'draw', 'agg_level'), measure.vars = indicators)
  dt <- dcast(dt, code + name + agg_level + draw ~ year + variable, value.var = 'value')
  
  # calculate aroc
  for (i in indicators) {
    dt[, (paste0(i, '_change')) := (get(paste0(year_end, '_', i)) - get(paste0(year_start, '_', i)))/(as.numeric(year_end) - as.numeric(year_start))]
  }
  
  # clean up
  dt[, grep('20', names(dt)) := NULL]
  # ----------------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------------------------------------------
  # Save formatted squeezed results
  
  # loop over admins
  for (a in 0:2) {
    message(a)
  
    # loop over indicators
    for (i in indicators) {
      message(i)
      
      # subset by admin
      dt_ad <- dt[agg_level == paste0('ADM', a)]
      
      # get mean, upper, and lower
      dt_ad[, mean := lapply(.SD, mean, na.rm = T), 
            .SDcols = paste0(i, '_change'), by = 'code']
      dt_ad[, upper := lapply(.SD, quantile, probs = 0.975, na.rm = T), 
            .SDcols = paste0(i, '_change'), by = 'code']
      dt_ad[, lower := lapply(.SD, quantile, probs = 0.025, na.rm = T), 
            .SDcols = paste0(i, '_change'), by = 'code']
      
      # save summary and free up space
      ad_summary <- merge(sp_list[[a+1]], 
                          unique(dt_ad[, c('code', 'mean', 'upper', 'lower')]), 
                          by.x = paste0('ADM', a, '_CODE'), by.y = 'code',
                          allow.cartesian = T)
      write.csv(ad_summary, paste0(share_dirs[[i]], 'pred_derivatives/treatment_change/',
                                   i, '_', year_start, '_', year_end, '_', country, '_change_summary_ad', a, '.csv'))
      
      # reshape draw object wide
      dt_ad <- dcast(dt_ad, code ~ draw, value.var = paste0(i, '_change'))
      setnames(dt_ad, 'code', paste0('ADM', a, '_CODE'))
      
      # check for NAs
      message('TESTING: Percent of NA rows per column is: ', mean(is.na(dt_ad[, V1])), '%')
      
      # save draw objects
      outpath <- paste0(paste0(share_dirs[[i]], 'pred_derivatives/treatment_change/',
                               i, '_', year_start, '_', year_end, '_', country, '_change_ad', a, '.RData'))
      save(dt_ad, file = outpath)
      
    }
    
  }
  # ---------------------------------------------------------------------------------------------------------
  

  # ------------------------------------------------------
  # End function
  return(paste0('Treatment analysis complete for ', country, '. Results saved.'))
  
}
# --------------------------------------------------------