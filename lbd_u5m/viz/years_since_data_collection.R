## Visualize input data by years since collection

library(data.table)
library(ggplot2)
library(scales)

## Define filepaths

in_fp <- '<<<< FILEPATH REDACTED >>>>'
out_dir   <- paste0(
  '<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/'
)
dir.create(out_dir, showWarnings = FALSE) 
out_fp <- paste0(out_dir, 'years_since_data_collection.png')

## ~~~ Read and prep data ~~~
dt <- readRDS(in_fp)

# Keep only surveys pre-2018 and observations since 2000
dt <- dt[(year >= 2000) & (year <= 2017) & (svyyr <= 2017), ]

# Get time (in years) from life-year experienced to survey date
dt[, yr_diff := svyyr - year ]
# (this should never be negative)
dt <- dt[ yr_diff >= 0, ]
dt[yr_diff <= 1, yr_diff := 1]

# Aggregate by years since data collection
dt_agg <- dt[, .(N = sum(N, na.rm=T)), by=yr_diff ][order(yr_diff)]
dt_agg[, prop_n := N/sum(N)]

# Make a plot
plt <- ggplot(data=dt_agg, aes(x=yr_diff, y=prop_n)) + 
    geom_bar(stat='identity') + 
    theme_minimal() + 
    scale_y_continuous(
        breaks=seq(0, 0.12, by=0.02),
        labels=scales::percent_format()
    ) + 
    labs(
        x = 'Years since data collection',
        y = 'Proportion of data',
        title = ''
    )

png(out_fp, height=600, width=1200)
print(plt)
dev.off()
