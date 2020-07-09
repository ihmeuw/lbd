####################################################################################################
## Description:   Make plots of predictive validity metrics.
##
## Inputs:        Predictive validity metrics for all specified model run_dates (ie, pv_metrics.csv)
##
## Output:        PDF of plots in the specified save location.
####################################################################################################

require(data.table)
require(ggplot2)

plot_pv_results <- function(indicator, indicator_group, run_dates, labels = NULL, save_file) {
  
  # load and combine PV metrics
  pv <- rbindlist(lapply(run_dates, function(rd) {
    cbind(model = rd, fread(paste0('<<<< FILEPATH REDACTED >>>>/pv_metrics.csv')))
  }))
  pv <- melt(pv, id.vars = c('model', 'Level', 'Year', 'Region', 'OOS', 'Mean Pred.', 'Mean Obs.', 'Median SS'))
  
  includes_oos <- pv[OOS == T, .N] > 0
  
  # format for plotting
  pv[, model := factor(model, levels = run_dates, if (!is.null(labels)) labels = labels else labels = run_dates)]
  pv[, Level := factor(Level, levels = c('country', 'ad1', 'ad2'), labels = c('Country', 'Admin. 1', 'Admin. 2'))]
  pv[, Year := factor(ifelse(is.na(Year), 'All Years', Year), levels = c('All Years', sort(unique(Year))))]
  pv[, Region := factor(ifelse(is.na(Region), 'All Regions', Region), levels = c('All Regions', sort(unique(Region))))]
  pv[, OOS := factor(OOS, levels = c(F, T), labels = c('In sample', 'Out of sample'))]
  pv[, variable := factor(variable, levels = c('Corr.', 'Mean Err.', 'RMSE', '95% Cov.'), labels = c('Correlation', 'Mean Error', 'RMSE', 'Coverage'))]
  pv[variable == 'Coverage', value := 100*value]
  pv <- pv[!is.na(variable)]
  
  # create file
  pdf(save_file, width = 14, height = 8)
  
  # plot results for all years, all regions combined
  p <- ggplot(pv[Region == 'All Regions' & Year == 'All Years',],
              aes(x = model, y = value, color = OOS, group = OOS)) +
    facet_grid(variable ~ Level, scales = 'free_y') +
    geom_point(size = 3) + geom_line() +
    labs(title = 'Model Performance') +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  
  # plot in sample results by year
  p <- ggplot(pv[Region == 'All Regions' & OOS == 'In sample',],
              aes(x = Year, y = value, color = model, group = model)) +
    facet_grid(variable ~ Level, scales = 'free_y') +
    geom_point() + geom_line() +
    labs(title = 'Model Performance (In sample, by year)') +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  
  # plot out of sample results by year
  if (includes_oos) {
    p <- ggplot(pv[Region == 'All Regions' & OOS == 'Out of sample',],
                aes(x = Year, y = value, color = model, group = model)) +
      facet_grid(variable ~ Level, scales = 'free_y') +
      geom_point() + geom_line() +
      labs(title = 'Model Performance (Out of sample, by year)') +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  }
  
  # plot in sample results by region
  p <- ggplot(pv[Year == 'All Years' & OOS == 'In sample',],
              aes(x = Region, y = value, color = model, group = model)) +
    facet_grid(variable ~ Level, scales = 'free_y') +
    geom_point() + geom_line() +
    labs(title = 'Model Performance (In sample, by region)') +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  
  # plot out of sample results by year
  if (includes_oos) {
    p <- ggplot(pv[Year == 'All Years' & OOS == 'Out of sample',],
                aes(x = Region, y = value, color = model, group = model)) +
      facet_grid(variable ~ Level, scales = 'free_y') +
      geom_point() + geom_line() +
      labs(title = 'Model Performance (Out of sample, by region)') +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  }
  
  dev.off()
  return(paste('Plots saved:', save_file))
}
