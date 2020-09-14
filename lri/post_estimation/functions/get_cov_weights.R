####################################################################################################
## Description:   Get covariate "importance" scores.
####################################################################################################
get_cov_weights <- function(indicator, indicator_group, run_date, Regions, outputdir) {
  
  # use mbg_central functions to calculate covariate weights and make plots
  all_plots <- lapply(Regions, function(r) {
    
    # calculate weights (this auto-saves the output)
    tree_error <- try(get.cov.wts(rd = run_date,
                                  ind = indicator,
                                  ind_gp = indicator_group,
                                  reg = r,
                                  age = 0,
                                  holdout = 0), silent = 1)
    if (inherits(tree_error, 'try-error')) {
      if (grepl('Non-tree model detected!', tree_error, fixed = TRUE)) {
        message('Non-tree model detected error: will not generate covariate importance plot for xgboost.')
      }
    } else {
      cov.wts <- get.cov.wts(rd = run_date,
                             ind = indicator,
                             ind_gp = indicator_group,
                             reg = r,
                             age = 0,
                             holdout = 0)
      
      # make plots
      cov.plots <- plot.cov.wts(rd = run_date,
                                ind = indicator,
                                ind_gp = indicator_group,
                                reg = r,
                                age = 0,
                                holdout = 0)
      cov.plots <- cov.plots + labs(x='', y='', title=r)
      return(cov.plots)
    }
  })
  names(all_plots) <- Regions
  
  # save plots (combined, and individually)
  pdf(paste0(outputdir, '/diagnostic_plots/cov_wts_all_regions.pdf'), width=10, height=7)
  for (r in Regions) print(all_plots[[r]])
  dev.off()
  
  return('Plots saved!')
}