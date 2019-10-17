## #############################################################################
## 
## MAP PROPORTION OF U5 DEATHS DUE TO SHOCKS, 2000-2017
## 
## Purpose: Check the proportion of under-5 deaths due to conflict during the 
##   study period.
## 
## #############################################################################

rm(list=ls())

## SET INPUTS

run_date <- '<<<< REDACTED >>>>'


## Imports and function sourcing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(data.table)

# GBD mapping function
source(paste0("<<<< FILEPATH REDACTED >>>>",
              "<<<< FILEPATH REDACTED >>>>"))
# GBD shared functions
source('<<<< FILEPATH REDACTED >>>>')
source('<<<< FILEPATH REDACTED >>>>')

# Set up output directory
out_dir   <- paste0(
  '<<<< FILEPATH REDACTED >>>>',gsub('-','_',Sys.Date()),'/'
)
dir.create(out_dir, showWarnings = FALSE)


## Get locations and causes to map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Locations
lookup_table <- fread('<<<< FILEPATH REDACTED >>>>')
adm0_codes <- unique(fread(paste0(
  '<<<< FILEPATH REDACTED >>>>',run_date,
  '<<<< FILEPATH REDACTED >>>>'
))$ADM0_CODE)
pull_locs <- lookup_table[ gadm_geoid %in% adm0_codes, .(loc_id, iso3) ]
setnames(pull_locs, 'loc_id', 'location_id')


# Pull all-cause deaths for LMICs in the under-5 age groups ~~~~~~~~~~~~~~~~~~~~
u5_age_groups <- 2:5

deaths_no_shock <- get_envelope(
  age_group_id=u5_age_groups,
  sex_id=3,
  year_id=2000:2017,
  location_id=pull_locs$location_id,
  gbd_round_id=5, # GBD 2017
  with_shock=0,
  with_hiv=1
)
setnames(deaths_no_shock, 'mean', 'val_noshock')
deaths_no_shock <- deaths_no_shock[,
  .(val_noshock = sum(val_noshock)),
  by=.(location_id, year_id)
]

deaths_no_hiv <- get_envelope(
  age_group_id=u5_age_groups,
  sex_id=3,
  year_id=2000:2017,
  location_id=pull_locs$location_id,
  gbd_round_id=5, # GBD 2017
  with_shock=0,
  with_hiv=0
)
setnames(deaths_no_hiv, 'mean', 'val_nohiv')
deaths_no_hiv <- deaths_no_hiv[,
  .(val_nohiv = sum(val_nohiv)),
  by=.(location_id, year_id)
]


deaths_with_shock <- get_envelope(
  age_group_id=u5_age_groups,
  sex_id=3,
  year_id=2000:2017,
  location_id=pull_locs$location_id,
  gbd_round_id=5, # GBD 2017
  with_shock=1,
  with_hiv=1
)
setnames(deaths_with_shock, 'mean', 'val_withshock')
deaths_with_shock <- deaths_with_shock[,
  .(val_withshock = sum(val_withshock)),
  by=.(location_id, year_id)
]

# Combine into a single group
deaths_combined <- merge(
  x = deaths_no_shock,
  y = deaths_with_shock,
  by = c('location_id','year_id')
)
deaths_combined <- merge(
  x = deaths_combined,
  y = deaths_no_hiv,
  by = c('location_id','year_id')
)
deaths_combined <- merge(
  x = deaths_combined,
  y = pull_locs,
  by = 'location_id'
)

# Get number and proportion of deaths from shocks BY COUNTRY AND YEAR
deaths_combined[, hiv_deaths := val_noshock - val_nohiv]
deaths_combined[, shock_deaths := val_withshock - val_noshock]
deaths_combined[, pct_from_shocks := shock_deaths/val_withshock]
deaths_combined[, pct_from_hiv := hiv_deaths/val_withshock]

all_years_by_country <- deaths_combined[,
  .(val_noshock = sum(val_noshock), 
    val_withshock=sum(val_withshock), 
    shock_deaths=sum(shock_deaths),
    hiv_deaths=sum(hiv_deaths)),
  by=.(iso3)
]
all_years_by_country[, pct_from_shocks := shock_deaths/val_withshock]
all_years_by_country[, pct_from_hiv := hiv_deaths/val_withshock]

lmics_by_year <- deaths_combined[,
  .(val_noshock = sum(val_noshock), 
    val_withshock=sum(val_withshock), 
    shock_deaths=sum(shock_deaths)),
  by=.(year_id)
]
lmics_by_year[, pct_from_shocks := shock_deaths/deaths_with_shock]


## NUMBER PLUG
plug_sentence_1 <- paste0(
  "The country with the highest proportion of U5 shock deaths in a single year ",
  "was %s in %s, with %s child deaths from shocks (%s%% of total)."
)
highest_deaths <- deaths_combined[ pct_from_shocks == max(pct_from_shocks),]
highest_5 <- head(deaths_combined[order(-pct_from_shocks)], 5)
haiti_deaths <- deaths_combined[iso3=='HTI' & year_id==2010,]

syria_deaths <- deaths_combined[iso3=='SYR' & year_id >= 2012,]
syria_deaths <- syria_deaths[, .(val_withshock = sum(val_withshock), shock_deaths=sum(shock_deaths)), ]
syria_deaths[, pct_from_shocks := shock_deaths/val_withshock]

iqr_data <- quantile(deaths_combined[, pct_from_shocks], c(0.025, .25, .5, .75, 0.975))*100

print(sprintf(
  plug_sentence_1,
  highest_deaths[,iso3],
  highest_deaths[,year_id],
  round(highest_deaths[,shock_deaths],1),
  round(highest_deaths[,pct_from_shocks]*100,1)
))

plug_sentence_2 <- paste0(
  "Out of the 123 million under-5 deaths which we map, %s million (%s%%) ",
  "are attributable to fatal discontinuities."
)
print(sprintf(
  plug_sentence_2, 
  round(sum(deaths_combined$shock_deaths)/1E6,1),
  round(sum(deaths_combined$shock_deaths)/sum(deaths_combined$val_withshock)*100,2)
))

print(sprintf(
  "Of those, %s (%s%% of total shock deaths) are in countries with subnational raking.",
  round(sum(deaths_combined[iso3 %in% c('IND', 'IDN'), shock_deaths]), 0),
  round(sum(deaths_combined[iso3 %in% c('IND', 'IDN'), shock_deaths])/sum(deaths_combined$shock_deaths)*100,1)
))

## MAKE MAP
all_years_by_country[, ihme_loc_id := iso3 ]
all_years_by_country[, mapvar := pct_from_shocks ]

limits <- c(0:5/100,999)
labels <- c("0-1%","1-2%","2-3%","3-4%","4-5%","Over 5%")

gbd_map(
  data           = all_years_by_country,
  limits         = limits,
  labels         = labels,
  col            = 'YlGnBu',
  col.reverse    = FALSE,
  na.color       = "#999999",
  title          = "Proportion of under-5 deaths from conflict, 2000-2017",
  legend.title   = "Proportion\nfrom conflict",
  legend.cex     = 1,
  legend.columns = 1,
  legend.shift   = c(-5,-20),
  fname          = paste0(out_dir, '/proportion_of_shocks_deaths.pdf')
)


## MAKE HIV MAP
all_years_by_country[, ihme_loc_id := iso3 ]
all_years_by_country[, mapvar := pct_from_hiv ]

limits <- c(0:5/50,999)
labels <- c("0-2%","2-4%","4-6%","6-8","8-10%","Over 10%")

gbd_map(
  data           = all_years_by_country,
  limits         = limits,
  labels         = labels,
  col            = 'YlGnBu',
  col.reverse    = FALSE,
  na.color       = "#999999",
  title          = "Proportion of under-5 deaths from HIV, 2000-2017",
  legend.title   = "Proportion\nfrom HIV",
  legend.cex     = 1,
  legend.columns = 1,
  legend.shift   = c(-5,-20),
  fname          = paste0(out_dir, '/proportion_of_hiv_deaths.pdf')
)
