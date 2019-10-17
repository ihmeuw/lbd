## #############################################################################
## 
## Make India diagnostic plots
## 
## Purpose: Plot the following for India:
##  - Range and IQR of mean district observations by year
##  - Data availability (number of individuals observed) by year
## 
## #############################################################################

## SET INPUTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

run_date <- '<<<< REDACTED >>>>'
shp_version <- '<<<< REDACTED >>>>'
plot_yrs <- c(2000, 2005, 2010, 2015, 2017)

mbg_input_data_file <- '<<<< FILEPATH REDACTED >>>>'

## Load libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
user            <- Sys.info()['user']                                  # user running this code
core_repo       <- sprintf('<<<< FILEPATH REDACTED >>>>',user) # main MBG repo
ig_repo         <- sprintf('<<<< FILEPATH REDACTED >>>>',user)      # u5m specific repo
indicator_group <- 'u5m'                                               # indicator group is u5m, this shouldn't change
indicator       <- 'died'                                              # indicator is died, this shouldn't change 
commondir <- '<<<< FILEPATH REDACTED >>>>'

viz_dir <- paste0('<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()))
dir.create(viz_dir, showWarnings=FALSE)

# run some setup scripts
source(sprintf('<<<< FILEPATH REDACTED >>>>',ig_repo)) 

source(sprintf('<<<< FILEPATH REDACTED >>>>',core_repo))
source(sprintf('<<<< FILEPATH REDACTED >>>>',core_repo))
# Load all needed packages and MBG functions
package_list <- c(t(read.csv(sprintf('<<<< FILEPATH REDACTED >>>>',commondir),header=FALSE)))
mbg_setup(package_list=package_list, repos=core_repo)


## ** INDIA RANGE AND UI PLOTS ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Prep data

# Get admin2 data for India only and match on admin1
ad2_raw <- fread(sprintf(
  '<<<< FILEPATH REDACTED >>>>',
  run_date
))
ad2_meta <- read.dbf(
  get_admin_shapefile(admin_level=2, suffix='.dbf', version=shp_version)
) %>% as.data.table
ind_meta <- ad2_meta[NAME_0=='India', .(NAME_1, ADM2_CODE)]
setnames(ind_meta, 'NAME_1', 'State')
ind_full <- merge(
  x = ad2_raw,
  y = ind_meta,
  by = c('ADM2_CODE')
)
# Keep only plotting years
ind_full <- ind_full[ year %in% plot_yrs, ]
# Get range and IQR of mean district observations by year
ind_agg <- ind_full[, 
  .(IQR = range(mean, na.rm=T),
    q_range = quantile(mean, .75, na.rm=T) - quantile(mean, .25, na.rm=T)
  ), 
  by=.(year, State)
]
ind_agg[, State := as.character(State) ]

## PLOT
# IQR plot
iqr_fig <- ggplot(data=ind_agg, aes(x=year, y=IQR, label=State)) + 
  geom_text(col='blue') + 
  theme_minimal()
pdf(paste0(viz_dir,'<<<< FILEPATH REDACTED >>>>'), height=7, width=15)
plot(iqr_fig)
dev.off()

# Range plot
range_fig <- ggplot(data=ind_agg, aes(x=year, y=q_range, label=State)) + 
  geom_text(col='blue') + 
  theme_minimal()
pdf(paste0(viz_dir,'<<<< FILEPATH REDACTED >>>>'), height=7, width=15)
plot(iqr_fig)
dev.off()


## ** INDIA INPUT DATA PLOTS ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Prep data
in_data <- readRDS(mbg_input_data_file)
ind_only <- in_data[country=='IND' & year >= 2000 & year <= 2017,]
# Sum observations by year
ind_data_by_yr <- ind_only[, .(N = sum(N, na.rm=T)), by=year]

## Plot it
fig2 <- ggplot(data=ind_data_by_yr, aes(x=year, y=N)) + 
  geom_line() + 
  theme_minimal()
pdf(paste0(viz_dir,'<<<< FILEPATH REDACTED >>>>'), height=7, width=10)
plot(fig2)
dev.off()
