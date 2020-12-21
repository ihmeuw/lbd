# --------------------------------------------------------------------------------------------------
# make_area_plots()
#
# Function that makes area plots for changes over time in any_ors, rhf_only, and no_ort
#
# Inputs:
# run_date - run date for current model
# holdout - which holdout we're running this for (0 is the full model)
# regions - country or vector of countries (order doesn't matter)
# region_names - full names of the countries (same order as regions)
#
# Outputs (saved in /results_maps/):
# - PDFs containing the area plots for all countries, each saved separately
# --------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------
# Start function
make_area_plots <- function(run_date,
                            holdout = 0,
                            regions = c('MLI', 'SEN', 'SLE'),
                            region_names = c('Mali', 'Senegal', 'Sierra Leone')) {
  # --------------------------------------------------------------------------------


  # -------------------------------------------------------------------------------
  # Set-up
  
  # set indicator group and indicators
  indicator_group <-  'ort'
  indicators <- c('no_ort', 'rhf_only', 'any_ors')
  
  # Define share directories
  share_dirs <- paste0('<<<< FILEPATH REDACTED >>>>/')
  names(share_dirs) <- indicators
  outdir <- paste0('<<<< FILEPATH REDACTED >>>>')
  # -------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------------------------
  # Load and clean data
  
  # load admin 2 data
  filepaths <- list()
  for (i in indicators) filepaths[[i]] <- list.files(share_dirs[[i]], pattern = 'admin_2_squeezed')
  dat <- lapply(paste0(share_dirs, filepaths), fread)
  names(dat) <- indicators
  
  # clean admin 2
  for (i in indicators) dat[[i]][, indicator := i]
  dat <- rbindlist(dat)
  
  # load admin 0 data
  filepaths <- list()
  for (i in indicators) filepaths[[i]] <- list.files(share_dirs[[i]], pattern = 'admin_0_squeezed')
  dat0 <- lapply(paste0(share_dirs, filepaths), fread)
  names(dat0) <- indicators
  
  # clean admin 0
  for (i in indicators) dat0[[i]][, indicator := i]
  dat0 <- rbindlist(dat0)
  # -----------------------------------------------------------------------------------------------
  
  
  # -----------------------------------------------------------------------------------------
  # Make and save the area plots
  
  # plotting function
  plot_treatment <- function(ad_dt, ad_name, years) {
    ad_dt$indicator <- gsub('any_ors', 'C', ad_dt$indicator)
    ad_dt$indicator <- gsub('rhf_only', 'B', ad_dt$indicator)
    ad_dt$indicator <- gsub('no_ort', 'A', ad_dt$indicator)
    year_labels <- years
    year_labels[c(TRUE, FALSE)] <- ''
    ggplot(ad_dt, aes(x = year, y = mean, fill = indicator)) +
      geom_area(position = 'stack', colour = 'black') + 
      ggtitle(ad_name) + 
      scale_x_continuous(expand = c(0, 0), breaks = years, labels = year_labels) + 
      scale_y_continuous(expand = c(0, 0)) + ylab('Coverage') + xlab('Year') + 
      scale_fill_manual(values = c('#BEBADA', '#FFFFB3', '#8DD3C7'), 
                        name = 'Treatment', 
                        labels = c('No ORT', 'RHF only', 'Any ORS')) +
      theme(text = element_text(size=36), axis.text.x = element_text(angle=45, hjust=1))
  }
  
  # loop over countries to plot for admin 2
  for (i in 1:length(region_names)) {
    message(region_names[[i]])
    
    # subset
    cty_dat <- dat[ADM0_NAME == region_names[[i]]]
    
    pdf(paste0(outdir, 'ad2_area_plots_', region_names[[i]], '.pdf'), width = 11, height = 8)
    
    # plot proportional distribution over time at admin 2 level
    for (ad in unique(cty_dat$ADM2_CODE)) {
      ad_dt <- cty_dat[ADM2_CODE == ad]
      print(plot_treatment(ad_dt, 
                           ad_name = paste0(region_names[[i]], ': ', ad_dt$ADM2_NAME[1]),
                           years = unique(ad_dt$year)))
    }
    
    dev.off()
    
  }
  
  # loop over countries to plot for admin 0
  for (i in 1:length(region_names)) {
    message(region_names[[i]])
    
    # subset
    cty_dat <- dat0[ADM0_NAME == region_names[[i]]]
    
    pdf(paste0(outdir, 'ad0_area_plots_', region_names[[i]], '.pdf'), width = 11, height = 8)
    
    # plot proportional distribution over time at admin 2 level
    for (ad in unique(cty_dat$ADM0_CODE)) {
      ad_dt <- cty_dat[ADM0_CODE == ad]
      print(plot_treatment(ad_dt, 
                           ad_name = paste0(ad_dt$ADM0_NAME[1]),
                           years = unique(ad_dt$year)))
    }
    
    dev.off()
    
  }
  # -----------------------------------------------------------------------------------------
  

  # ----------------------------------------------------
  # End function
  return('Making area plots complete. Plots saved.')
  
}
# ------------------------------------------------------
