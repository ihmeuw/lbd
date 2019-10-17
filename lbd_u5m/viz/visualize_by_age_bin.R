##
## Visualize distributions of prepped dataset by age groups
## Purpose: Checking the distribution of age groups across two datasets
##

library(data.table)
library(ggplot2)


## Set input and output filepaths

data_title <- 'RTR Data (June 7 2019)'
data_file <- 'died_global_7ab_annual_weightedSS_20190606_rtr.RDS'
in_fp <- paste0('<<<< FILEPATH REDACTED >>>>',data_file)

out_dir   <- paste0(
  '<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/'
)
dir.create(out_dir, showWarnings = FALSE) 
out_fp <- paste0(out_dir, '<<<< FILEPATH REDACTED >>>>',gsub('.RDS','',data_file),'.png')


##  Load data and format by age bin

dt <- readRDS(in_fp)[(year >= 2000) & (year <= 2017),]

# Get proportion of deaths and sample size in each age bin
dt_agg <- dt[, .(N=sum(N,na.rm=T), died=sum(died, na.rm=T)), by=.(ab, age)][order(age)]
dt_agg[, N := N/sum(N)]
dt_agg[, died := died/sum(died)]
dt_agg$ab <- factor(x=dt_agg$ab, levels=dt_agg$ab)

dt_agg <- melt(
    dt_agg, 
    id.vars=c('ab','age'),
    value.vars=c('died','N')
)
dt_agg[variable =='N', variable:="Life years"]
dt_agg[variable =='died', variable:="Deaths"]

# Visualize as ggplot

plt <- ggplot(data=dt_agg, aes(x=ab, y=value)) + 
    geom_bar(stat='identity') + 
    facet_wrap('variable', ncol=2, scales='free_y') +
    theme_bw() + 
    labs(
        x = 'Age bin',
        y = 'Proportion of input data',
        title = paste0('Distribution of data by age bin: ', data_title)
    )

png(out_fp, height=600, width=1200, units='px', res=150)
print(plt)
dev.off()
