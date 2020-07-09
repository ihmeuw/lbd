# --------------------------------------------------------------------------------
# Purpose: Functions used more than once to present ORT data scoping results
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Packages

library('data.table')
library('knitr')
library('magrittr')
library('ggplot2')
library('kableExtra')
source('<<<< FILEPATH REDACTED >>>>/GBD_WITH_INSETS_MAPPING_FUNCTION.R')

# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Plotting functions

# make a pretty table
pretty_table <- function(data, num_row=4) {
  kable(data) %>%
    # simple striped, bordered table
    kable_styling(bootstrap_options = c('striped', 'bordered'), position = 'center') %>%
    # rows are black
    row_spec(1:num_row, color='black') %>%
    # return table
    return()
}

# make a pretty graph
pretty_graph <- function(data) {
  # increase font size
  theme_set(theme_gray(base_size = 14))
  # stacked bar chart
  ggplot(data=data, aes(x=Question, y=Number, fill=Survey)) +
    geom_bar(stat='identity') + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(~Group, scales = 'free_x', space = 'free_x')
}

# national data coverage plots using GBD mapping function
coverage_map <- function(inc_data, year_end) {
  # subset by year to plot
  map_data <- inc_data[year_group == year_end]
  # make dataframes
  map_data <- as.data.frame(map_data)
  # Define deciles for bucket cutoffs
  quantiles <- c(1L, 1L, 2L, 3L, 4L, 6L, 8L, 30L)
  # Define labels for bucket cutoffs
  labels <- c("1", "2", "3", "4", "5-6", "7-8", "9+")
  # map and output using the mapping function
  # National
  gbd_map(data = map_data,
          limits = quantiles,
          labels=labels,
          col= 'PiYG', col.reverse= TRUE,
          na.color= '#999999',
          title = 'Number of surveys available for ORT in Stage 1 and 2 countries',
          legend.title = 'Surveys per location/year',
          legend.cex = 1.1,
          legend.columns = 2,
          legend.shift = c(-5,-10))
}

# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# File saving functions (optional)

# save plot object as a PDF
save_plot <- function(plot, filepath) {
  pdf(filepath)
  print(plot)
  dev.off()
}

# --------------------------------------------------------------------------------