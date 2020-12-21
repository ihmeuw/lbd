# ---------------------------------------------------------------------------------------------------------
# analyze_efficacy()
#
# Function that calculates how effective RHF would have to be (hypothetically) to have not lost lives in
# areas where RHF has decreased to a greater extend than ORS increased
#
# NB: We decided to NOT include this analysis in the manuscript, but kept the code as some of the outputs
#     are used in later steps of the full analysis pipeline
#
# Inputs:
# run_date - run date for current model
# year_start - start date for analysis (start of modeling period, pre/post intervention, etc.)
# year_end - end date for analysis (end of modeling period, pre/post intervention, etc.)
# nreps - number of ORS:RHF efficacy ratios to test
# ors_ef - ORS efficacy (percent of deaths that could be prevented with perfect coverage)
# max_ef_ratio - the maximum ratio to test before capping it
# holdout - which holdout we're running this for (0 is the full model)
# country - iso3 for single country to run analysis for
#
# Outputs (saved in /pred_derivatives/ors_rhf_efficacy/ in the no_ort folder):
# - Matrices of draw level ORS:RHF efficacy between year_start and year_end
# - CSVs of summarized mean, upper, and lower values
# - Draw summary file including both the number and proportion of draws that showed ORS decline,
#   RHF replaced, and ORS:RHF efficacy ratio needed of greater than 100
# ---------------------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
# Start function
analyze_efficacy <- function(run_date,
                             year_start,
                             year_end,
                             country,
                             nreps = 10000,
                             ors_ef = 0.69,
                             max_ef_ratio = 100,
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
  dir.create(paste0(share_dirs[[1]], 'pred_derivatives/ors_rhf_efficacy/'), showWarnings = F)
  # ---------------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------------------------------------------
  # Load and clean modeled results
  
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
  dt[, c(indicators) := NULL]
  setnames(dt, paste0(indicators, '_sqz'), indicators)
  
  # subset by years and region in analysis
  dt <- dt[region == country & (year == year_start | year == year_end)]
  
  # reshape
  dt <- melt(dt, id.vars = c('code', 'year', 'name', 'draw', 'agg_level'), measure.vars = indicators)
  dt <- dcast(dt, code + name + agg_level + draw ~ year + variable, value.var = 'value')
  names(dt) <- c('code', 'name', 'agg_level', 'draw', 'c1', 'b1', 'a1', 'c2', 'b2', 'a2')
  
  # create column to indicate whether or not RHF was replaced and whether ORS declined
  dt[, ors_declined := a2 < a1]
  dt[, rhf_replaced := c2 < c1]
  dt[, c('c1', 'c2') := NULL]
  # ---------------------------------------------------------------------------------------------------------
  
  
  # ------------------------------------------------------------------------------------------------
  # Calculate ORS:RHF efficacy by draw
  
  # number of draws and admins we're using
  ndraws <-  length(unique(dt$draw))
  nads <- length(unique(dt$code))
  
  # get theoretical RHF efficacies (with ORS max 10 times as effective as RHF)
  rel_ef <- seq(1, max_ef_ratio, length.out = nreps)
  rhf_ef <- ors_ef/rel_ef
  
  # create RHF efficacy table
  rhf_ef <- rep(rhf_ef, times = nrow(dt))
  draws <- rep(sort(rep(unique(dt$draw), times = nreps)), times = nads)
  ads <- sort(rep(unique(dt$code), times = ndraws*nreps))
  rhf_ef_ads <- data.table(code = ads, draw = draws, rhf_efficacy = rhf_ef)
  
  # add to data table
  dt <- merge(dt, rhf_ef_ads, by = c('code', 'draw'), allow.cartesian = T)
  
  # calculate
  dt[, lives_saved := (a2*ors_ef + b2*rhf_efficacy)/(a1*ors_ef + b1*rhf_efficacy),
     by = c('code', 'draw')]
  dt[, efficacy_ratio := 0.69/rhf_efficacy]
  
  # clean up
  dt[, c('rhf_efficacy', 'a1', 'a2', 'b1', 'b2') := NULL]
  # ------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------------------------------------------
  # Create indicator for ORS:RHF efficacy needed to have saved lives
  
  # set ratio to the point at which the efficacy ratio is closest to 1
  get_ratio <- function(x, y) {
    idx <- min(which(x >= 1))
    return(y[idx])
  }
  dt[, ratio_needed := lapply(.SD, get_ratio, y = efficacy_ratio), 
     .SDcols = 'lives_saved', by = c('code', 'draw')]
  
  # if ORS declined, set ratio to 999
  dt[ors_declined == TRUE, ratio_needed := 999]
  
  # if RHF was replaced, set ratio to 0
  dt[rhf_replaced == TRUE, ratio_needed := 0]
  
  # if any NAs remain, set to 100 and will report that efficacy ratio needs to be > 100
  dt[is.na(ratio_needed), ratio_needed := 100]
  
  # cleanup
  dt[, c('rhf_replaced', 'lives_saved', 'efficacy_ratio', 'ors_declined') := NULL]
  dt <- unique(dt)
  # ---------------------------------------------------------------------------------------------------------
  
  
  # --------------------------------------------------------------------------------------------------------------------------
  ## Summarize and subset places with declines in ORT but not in ORS
  
  # get number of draws that show ORS decline
  dt[, ors_declined_num:= sum(ratio_needed == 999),
     by = 'code']
  
  # get number of draws that show ORT increased
  dt[, rhf_replaced_num := sum(ratio_needed == 0),
     by = 'code']
  
  # get number of draws that show max ORS:RHF efficacy greater than maximum of 100
  dt[, efficacy_100_num := sum(ratio_needed == 100),
     by = 'code']
  
  # summarize and save
  dt_summary <- unique(dt[, c('code', 'name', 'agg_level', 'ors_declined_num', 'rhf_replaced_num', 'efficacy_100_num')])
  dt_summary[, ors_declined_prop := ors_declined_num/ndraws]
  dt_summary[, rhf_replaced_prop := rhf_replaced_num/ndraws]
  dt_summary[, efficacy_100_prop := efficacy_100_num/ndraws]
  write.csv(dt_summary, file = paste0(share_dirs[[1]], 'pred_derivatives/ors_rhf_efficacy/draw_summary_', year_start, '_', year_end, '_', country, '.csv'))
  
  # subset
  dt <- dt[, c('code', 'draw', 'name', 'agg_level', 'ratio_needed')]
  dt <- dt[ratio_needed == 999 | ratio_needed == 0, ratio_needed := NA]
  # --------------------------------------------------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------------------------------------
  ## Save formatted squeezed results
  
  for (a in 0:2) {
    message(a)
      
    # subset by admin
    dt_ad <- dt[agg_level == paste0('ADM', a)]
    
    # get mean, upper, and lower
    dt_ad[, mean := lapply(.SD, mean, na.rm = T), 
          .SDcols = 'ratio_needed', by = 'code']
    dt_ad[, upper := lapply(.SD, quantile, probs = 0.975, na.rm = T), 
          .SDcols = 'ratio_needed', by = 'code']
    dt_ad[, lower := lapply(.SD, quantile, probs = 0.025, na.rm = T), 
          .SDcols = 'ratio_needed', by = 'code']
    
    # save summary and free up space
    ad_summary <- merge(sp_list[[a+1]],
                        unique(dt_ad[, c('code', 'mean', 'upper', 'lower')]), 
                        by.x = paste0('ADM', a, '_CODE'), by.y = 'code',
                        allow.cartesian = T)
    write.csv(ad_summary, paste0(share_dirs[[1]], 'pred_derivatives/ors_rhf_efficacy/',
                                 'efficacy_ratio_needed_', year_start, '_', year_end, '_', country, '_ad', a, '.csv'))
    
    # reshape draw object wide
    dt_ad <- dcast(dt_ad, code ~ draw, value.var = 'ratio_needed')
    setnames(dt_ad, 'code', paste0('ADM', a, '_CODE'))
    
    # check for NAs
    message('TESTING: Percent of NA rows per column is: ', mean(is.na(dt_ad[, V1])), '%')
    
    # save draw objects
    outpath <- paste0(paste0(share_dirs[[1]], 'pred_derivatives/ors_rhf_efficacy/',
                             'efficacy_ratio_needed_', year_start, '_', year_end, '_', country, '_ad', a, '.RData'))
    save(dt_ad, file = outpath)
    
  }
  # ---------------------------------------------------------------------------------------------------
  
  
  # ------------------------------------------------------
  # End function
  return(paste0('Efficacy analysis complete for ', country, '. Results saved.'))
  
}
# --------------------------------------------------------
