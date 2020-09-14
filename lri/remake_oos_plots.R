library('data.table')
library('ggplot2')

indicator = 'has_lri'
indicator_group = 'lri'
run_date = c('2018_05_10_15_58_20_1000_draws')
samples = 1000
ho_strat = data.frame(run_date = run_date)
ho_strat[,'Spatial Holdout'] = c('ho_id')
period = data.frame(ppp = 0:4, Period = c('2000 - 2016', '2000 - 2003', '2004 - 2007', '2008 - 2011', '2012 - 2017'))
source('<<<< FILEPATH REDACTED >>>>/mbg_central/validation_functions.R')


diagfolder = '<<<< FILEPATH REDACTED >>>>'
dir.create(diagfolder)
draws.df = fread('<<<< FILEPATH REDACTED >>>>')

if(!any(names(draws.df) %in% 'ho_id')){
  draws.df[, ho_id:=nid]
}

dim(draws.df)
draws.df = merge(draws.df, period, by = 'ppp')
dim(draws.df)
drawcols <- grep("draw", names(draws.df), value = T)
data <- dplyr::select(draws.df, c(-drawcols))

#get_pv_table inputs
d = draws.df
indicator_group = indicator_group
rd = run_date
indicator=indicator
aggregate_on='ho_id'
result_agg_over = c('oos')
draws = as.numeric(samples)
out.dir = paste0(diagfolder,'/')
plot_ci = T
plot_ncol = 2

draws = 1000
coverage_probs = c(95)
result_agg_over = c("year", "oos")
weighted = TRUE
family = "binomial"
plot = TRUE
plot_by = NULL
plot_by_title = NULL
plot_ci_level = 95
ci_color = "grey"
point_alpha = 1
point_color = "black"
plot_title = indicator
save_csv = T

str_match <- stringr::str_match

d <- data.table(d)

if (!is.null(plot_by)) {
  if (!(plot_by %in% result_agg_over)) {
    stop("If you specify `plot_by`, you must also include that item in `result_agg_over`")
  }
}

## Get binomial predictions
message("Simulating predictive draws")
if (family == "binomial") {
  x <- sapply(1:draws, function(i) rbinom(nrow(d), round(d$N), d[[paste0("draw", i)]]))
  d[, Y := get(indicator) * round(N) / N] # adjust observed number of cases for effect of rounding N
}
if (family == "gaussian") {
  message("NICK: THIS SHOULD WORK IF YOU HAVE DRAWS OF PRECISION NAMED tau_1,...tau_1000 BUT THIS IS UNTESTED")
  x <- sapply(1:draws, function(i) rnorm(nrow(d), sqrt(d[[paste0("tau_", i)]]), d[[paste0("draw", i)]]))
  d[, Y := get(indicator)]
}

message(paste0("...NOTE: missing predictions for ", sum(is.na(x[, 1])), " of ", nrow(x), " (", round(100 * mean(is.na(x[, 1]))), "%) data points."))

## Get coverage for each data point
message("Calculate coverage")
one_side <- (1 - coverage_probs / 100) / 2
ui <- t(apply(x, 1, quantile, c(one_side, 1 - one_side), na.rm = T))
for (c in coverage_probs) {
  c_id <- which(c == coverage_probs)
  d[, paste0("clusters_covered_", c) := Y %between% list(ui[, c_id], ui[, c_id + length(coverage_probs)])]
}

## Collapse data and draws spatially
message("Collapse data and draws spatially")
d[, oos := (fold != 0)]
d[, p := get(indicator) / N]
d[, exposure := N * weight]

# Create list of pv metrics for all levels of aggregation in `aggregate_on`
message("Outputting predictive validity metrics for your models at each level of aggregation")
pv_table_list <- lapply(aggregate_on, function(agg_on) {
  message(paste0("  ", agg_on))
  by_vars <- unique(c("oos", "year", result_agg_over, agg_on))
  collapse_vars <- c("p", paste0("clusters_covered_", coverage_probs), paste0("draw", 1:draws))
  
  res <- d[!is.na(draw1),
           c(
             list(total_clusters = .N, exposure = sum(exposure)),
             lapply(.SD, function(x) weighted.mean(x, exposure, na.rm = T))
           ),
           keyby = by_vars, .SDcols = collapse_vars
           ]
  res[, mean_draw := rowMeans(.SD), .SDcols = paste0("draw", 1:draws)]
  res[, error := p - mean_draw]
  res[, abs_error := abs(error)]
  
  ## Collapse to calculate predictive validity metrics
  weighted.rmse <- function(error, w) {
    sqrt(sum((w / sum(w)) * ((error) ^ 2)))
  }
  if (weighted) res$weight <- res$exposure else res$weight <- 1
  res2 <- 
    res[, c(lapply(.SD, function(x) weighted.mean(x, weight)),
            rmse = weighted.rmse(error, weight),
            median_SS = median(exposure),
            cor = boot::corr(cbind(p, mean_draw), weight)),
        by = result_agg_over,
        .SDcols = c("error", "abs_error", "mean_draw", "p", paste0("clusters_covered_", coverage_probs))]
  
  setnames(res2, c("error", "abs_error", "p", paste0("clusters_covered_", coverage_probs)),
           c("me", "mae", "mean_p", paste0("coverage_", coverage_probs)))
  
  return(list(res = res, res2 = res2))
})

names(pv_table_list) <- aggregate_on

message("Making plots of aggregated data and estimates")

# Get unique levels of `plot_by` and set up a plot_by title if needed
if (!is.null(plot_by)) {
  plot_by_levels <- unique(pv_table_list[[1]][["res"]][, get(plot_by)])
} else {
  plot_by_levels <- NULL
}

if (is.null(plot_by_title)) plot_by_title <- plot_by

# Create a table of things to plot
plot_table <- CJ(
  aggregate_on = aggregate_on,
  oos = unique(pv_table_list[[1]][["res"]]$oos),
  plot_by_value = if (is.null(plot_by_levels)) NA else plot_by_levels
)

message("...saving plots here: ", out.dir)
# Loop over plots
i <- 2 #OOS plot
# Grab items from plot table
agg_on <- plot_table[i, aggregate_on]
oosindic <- plot_table[i, oos]
pb_val <- plot_table[i, plot_by_value]

# Set up titles
if (agg_on == "country") agg_title <- "Country"
if (agg_on == "ad1") agg_title <- "Admin 1"
if (agg_on == "ad2") agg_title <- "Admin 2"

res <- pv_table_list[[agg_on]][["res"]]
res2 <- pv_table_list[[agg_on]][["res2"]]

# Make a validation plot -----------------------------------------------------

# Set up filename and file
plot_filename <- paste0(
  indicator, "_validation_plot_",
  paste(c(
    as.character(agg_on),
    setdiff(result_agg_over, "oos")
  ),
  collapse = "_"
  ), "_",
  ifelse(oosindic, "OOS", "IS"),
  ifelse(is.na(pb_val), "", paste0("_", pb_val)),
  ".png"
)
message(paste("    ", plot_filename))
#png(paste0(out.dir, plot_filename), width = 12, height = 12, units = "in", res = 350)

# Subset data
fdata <- res[oos == oosindic, ]
if (!is.na(pb_val)) {
  setnames(fdata, plot_by, "plot_by_column") # convenience
  fdata <- fdata[plot_by_column == pb_val, ]
}

plot_data <- dplyr::select(fdata, c(-drawcols))
plot_data[order(-p)]

#worst nid fits
nids <- c(7721, 200617, 155335, 20796)
surveys <- filter(data, nid %in% nids) %>%
  dplyr::select(c('source', 'country','nid','year')) %>%
  unique.data.frame()
