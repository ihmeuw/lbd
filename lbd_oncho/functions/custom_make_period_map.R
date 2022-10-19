## Make map of period indices to run any time periods in your data (like annual instead of 5-year)
make_period_map <- function(modeling_periods) {
  data_period <- sort(modeling_periods)
  if (length(modeling_periods) == 1) {
    period_ids <- 1
  } else {
    period_ids <- seq(data_period)
  }
  period_map <- as.data.table(data_period)
  period_map <- period_map[, period_id := period_ids]
  return(period_map)
}