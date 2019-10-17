## #############################################################################
## Compare India UIs between LBD and GBD
## Purpose: Discuss the behavior of LBD uncertainty intervals in a meeting
##   with Lalit (on behalf of PHFI collaborators)
## #############################################################################

# Source libraries
library(data.table)
library(foreign)
library(ggplot2)
# Source GBD functions
library(mortdb, lib='<<<< FILEPATH REDACTED >>>>') # Must be in Rstudio
source('<<<< FILEPATH REDACTED >>>>')

## Set inputs
run_date <- '<<<< REDACTED >>>>'
shp_version <- '<<<< REDACTED >>>>'
out_dir <- paste0('<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/')
dir.create(out_dir, showWarnings = FALSE)

## Run script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Pull GBD inputs
gbd_locs <- get_location_metadata(
  location_set_id=21,
  gbd_round_id=5
)
ind_ad1 <- gbd_locs[parent_id==163,]

gbd_5q0_by_state <- mortdb::get_mort_outputs(
  model_name = '5q0',
  model_type = 'estimate',
  run_id = 'best',
  location_ids = ind_ad1$location_id,
  year_ids = 2000:2017,
  gbd_year = 2017,
  sex_ids = 3
)
gbd_5q0_by_state <- as.data.table( gbd_5q0_by_state )[
  estimate_stage_name=='GPR', .(location_id, year_id, mean, upper, lower)
]
gbd_5q0_by_state <- merge(
  x=gbd_5q0_by_state,
  y=ind_ad1[, .(location_id, location_name, map_id)],
  by='location_id',
  all.x=TRUE
)
gbd_5q0_by_state[, source := 'GBD 2017 Estimates']
gbd_5q0_by_state[, location_id := NULL ]


## Pull LBD inputs
in_dir <- sprintf('<<<< FILEPATH REDACTED >>>>',run_date)
lbd_ad1_est <- fread(paste0(in_dir,'<<<< FILEPATH REDACTED >>>>'))
lbd_ad1_est[, V1 := NULL]

ad1_meta <- foreign::read.dbf(paste0(
  '<<<< FILEPATH REDACTED >>>>',shp_version,
  '<<<< FILEPATH REDACTED >>>>'
))
ad1_meta <- as.data.table(ad1_meta)[ADM0_NAME=='India',.(ADM1_CODE,ADM1_NAME)]
# Update to match GBD names
ad1_meta[, ADM1_NAME := as.character(ADM1_NAME)]
ad1_meta[
  ADM1_NAME=='The Six Minor Territories',
  ADM1_NAME:='Union Territories other than Delhi'
]
ad1_meta[ADM1_NAME=='NCT of Delhi', ADM1_NAME:='Delhi']

ind_lbd_est <- merge(
  x=lbd_ad1_est,
  y=ad1_meta,
  by=c('ADM1_CODE'),
  all.y=TRUE
)
setnames(ind_lbd_est,c('ADM1_NAME','year'),c('location_name','year_id'))
ind_lbd_est[, ADM1_CODE := NULL ]
ind_lbd_est[, source := "LBD Estimates" ]
ind_lbd_est <- merge(
  x = ind_lbd_est,
  y = ind_ad1[, .(location_name, map_id)],
  by = 'location_name',
  all.x = TRUE
)

## Combine results
full_data <- rbindlist(
  list(ind_lbd_est, gbd_5q0_by_state),
  use.names = TRUE
)
full_data[
  location_name=='Union Territories other than Delhi',
  location_name:='Union Territories'
]

## Create plots by year
full_data[, lower_diff := lower - mean ]
full_data[, upper_diff := upper - mean ]
full_data[, mean_diff := 0 ]

for(this_year in 2000:2017){
  # Make plot
  year_plot <- ggplot(data=full_data[year_id==this_year,]) +
    geom_crossbar(
      data=function(x) x[source=='GBD 2017 Estimates'],
      aes(x=location_name, ymin=lower_diff, ymax=upper_diff, y=mean_diff, color=source),
      size=2.2, width=0, position = position_nudge(x = -.075)
    ) +
    geom_crossbar(
      data=function(x) x[source=='LBD Estimates'],
      aes(x=location_name, ymin=lower_diff, ymax=upper_diff, y=mean_diff, color=source),
      size=2.2, width=0, position = position_nudge(x = .075)
    ) +
    labs(
      title=paste0("Width of 95% UIs in ",this_year),
      x='State of India',
      y="Uncertainty Intervals (Normalized to Mean)",
      color='Source'
    ) +
    geom_hline(yintercept = 0, color='#666666', linetype=3) +
    scale_y_continuous(
      limits=c(-0.015, 0.015),
      breaks=seq(-0.015, 0.015, by=0.005)
    ) +
    theme(
      text = element_text(size=20),
      axis.text.x = element_text(size=12, angle = 90, hjust=1, vjust=0.5)
    )
  # Save to file
  pdf(sprintf('<<<< FILEPATH REDACTED >>>>',out_dir,this_year), width=15, height=10)
  print(year_plot)
  dev.off()
}
message("All done!")
