####################################################################################################
## Description:   Make plots of predictive validity metrics.
##
## Inputs:        Predictive validity metrics for all specified model run_dates (ie, pv_metrics.csv)
##
## Output:        PDF of plots in the specified save location.
####################################################################################################

require(data.table)
require(ggplot2)

plot_pv_results <- function(indicator, indicator_group, run_dates, labels = NULL, save_file, raked = NULL) {

  # load and combine PV metrics
  pv <- rbindlist(lapply(run_dates, function(rd) {
    cbind(model = rd, fread(paste("<<<< FILEPATH REDACTED >>>>", sep = "/")))
  }))
  if (!is.null(raked)) pv <- pv[Raked == raked, ]
  if (!"Year" %in% names(pv)) pv$Year <- NA
  if (!"Region" %in% names(pv)) pv$Region <- NA
  if (!"Point" %in% names(pv)) pv$Point <- NA
  pv <- pv[, c("model", "Level", "Year", "Region", "Point", "OOS", "Mean Pred.", "Mean Obs.",
               "Median SS", "Corr.", "Mean Err.", "Mean Abs. Err.", "RMSE", "95% Cov.")]
  pv <- melt(pv, id.vars = c("model", "Level", "Year", "Region", "Point", "OOS", "Mean Pred.", "Mean Obs.", "Median SS"))
  includes_oos <- pv[OOS == T, .N] > 0

  # format for plotting
  pv[, model := factor(model, levels = run_dates, if (!is.null(labels)) labels = labels else labels = run_dates)]
  pv[, Level := factor(Level, levels = c("country", "ad1", "ad2"), labels = c("Country", "Admin. 1", "Admin. 2"))]
  pv[, Year := factor(Year, levels = sort(unique(Year)))]
  pv[, Region := factor(Region, levels = sort(unique(Region)))]
  pv[, Point := factor(Point, levels = c(0, 1), labels = c('Polygons', 'Points'))]
  pv[, OOS := factor(OOS, levels = c(F, T), labels = c("In sample", "Out of sample"))]
  pv[, variable := factor(variable, levels = c("Corr.", "Mean Err.", "Mean Abs. Err.", "RMSE", "95% Cov."),
                          labels = c("Correlation", "Mean Error", "Mean Absolute Error", "RMSE", "Coverage"))]
  pv[variable == "Coverage", value := 100*value]

  # create file
  pdf(save_file, width = 14, height = 8)

  # plot overall results
  p <- ggplot(pv[is.na(Year) & is.na(Region) & is.na(Point), ], aes(x = model, y = value, color = OOS, group = OOS)) +
    facet_grid(variable ~ Level, scales = "free_y") +
    geom_point(size = 3) + geom_line() +
    labs(title = "Model Performance") +
    theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)

  # plot results by year
  if (pv[!is.na(Year), .N] > 0) {
    p <- ggplot(pv[!is.na(Year) & OOS == "In sample",], aes(x = Year, y = value, color = model, group = model)) +
      facet_grid(variable ~ Level, scales = "free_y") +
      geom_point() + geom_line() +
      labs(title = "Model Performance (In sample, by year)") +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)

    if (includes_oos) {
      p <- ggplot(pv[!is.na(Year) & OOS == "Out of sample",], aes(x = Year, y = value, color = model, group = model)) +
        facet_grid(variable ~ Level, scales = "free_y") +
        geom_point() + geom_line() +
        labs(title = "Model Performance (Out of sample, by year)") +
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
    }
  }

  # plot results by region
  if (pv[!is.na(Region), .N] > 0) {
    p <- ggplot(pv[!is.na(Region) & OOS == "In sample",], aes(x = Region, y = value, color = model, group = model)) +
      facet_grid(variable ~ Level, scales = "free_y") +
      geom_point() + geom_line() +
      labs(title = "Model Performance (In sample, by region)") +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)

    if (includes_oos) {
      p <- ggplot(pv[!is.na(Region) & OOS == "Out of sample",], aes(x = Region, y = value, color = model, group = model)) +
        facet_grid(variable ~ Level, scales = "free_y") +
        geom_point() + geom_line() +
        labs(title = "Model Performance (Out of sample, by region)") +
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
    }
  }

  # plot results by data type
  if (pv[!is.na(Point), .N] > 0) {
    p <- ggplot(pv[!is.na(Point) & OOS == "In sample",], aes(x = Point, y = value, color = model, group = model)) +
      facet_grid(variable ~ Level, scales = "free_y") +
      geom_point() + geom_line() +
      labs(title = "Model Performance (In sample, by data type)") +
      theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)

    if (includes_oos) {
      p <- ggplot(pv[!is.na(Point) & OOS == "Out of sample",], aes(x = Point, y = value, color = model, group = model)) +
        facet_grid(variable ~ Level, scales = "free_y") +
        geom_point() + geom_line() +
        labs(title = "Model Performance (Out of sample, by data type)") +
        theme(axis.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
      print(p)
    }
  }

  dev.off()
  return(paste("Plots saved:", save_file))
}
