## ###########################################################################
## 
## MAP TOTAL NUMBER OF CBH/SBH SOURCES BY COUNTRY
## 
## Purpose: Map the total number of CBH and SBH sources by country, 1998-2017,
##    using the standard GBD world map
## 
## ###########################################################################

rm(list=ls())

## Imports and function sourcing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(data.table)

# GBD mapping function
source(paste0("<<<< FILEPATH REDACTED >>>>",
              "<<<< FILEPATH REDACTED >>>>"))

## Prep data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# YEAR RANGE
start_year <- 2000
end_year   <- 2017

# INPUT FILEPATHS
in_file <- '<<<< FILEPATH REDACTED >>>>'
stages_fp <- '<<<< FILEPATH REDACTED >>>>'

# OUTPUT FILEPATHS
main_out_dir <- '<<<< FILEPATH REDACTED >>>>'
sub_out_dir <- paste0(main_out_dir, gsub('-', '_', Sys.Date()) )
dir.create(sub_out_dir, showWarnings = FALSE)
out_file_base <- paste0(sub_out_dir, '<<<< FILEPATH REDACTED >>>>')
in_data <- fread(in_file)
# Keep only years in year range
in_data <- in_data[(svyyr >= start_year) & (svyyr <= end_year), ]
# Subset to only needed columns and drop duplicates
in_data <- unique(in_data[,.(nid, data_type, country)])
in_data <- in_data[, .(mapvar = .N), by=.(country, data_type)]
setnames(in_data, 'country', 'ihme_loc_id')
sbh_data <- in_data[data_type =='sbh',]
cbh_data <- in_data[data_type =='cbh',]

# Add all countries in Stage 2, and drop countries in Stage 3
stage2 <- fread(stages_fp)[Stage != '3', .(iso3)]
setnames(stage2, 'iso3', 'ihme_loc_id')

sbh_data <- merge(
  x     = stage2,
  y     = sbh_data,
  by    = c('ihme_loc_id'),
  all.x = TRUE
)
sbh_data[is.na(mapvar), mapvar := 0]

cbh_data <- merge(
  x     = stage2,
  y     = cbh_data,
  by    = c('ihme_loc_id'),
  all.x = TRUE
)
cbh_data[is.na(mapvar), mapvar := 0]

# Add MCHS data to CBH data
cbh_data[ihme_loc_id=="CHN",mapvar := mapvar + (2012 - start_year + 1)]


## Create maps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

limits <- c(0, .5, 1.5, 2.5, 3.5, 4.5, 999)
labels <- c("None","1","2","3","4","5 or more")

## CBH Map
gbd_map(
  data           = cbh_data,
  limits         = limits,
  labels         = labels,
  col            = 'YlGnBu',
  col.reverse    = FALSE,
  na.color       = "#999999",
  title          = sprintf("Number of CBH data sources by country, %s-%s",
                           start_year,end_year),
  legend.title   = "CBH Sources",
  legend.cex     = 1,
  legend.columns = 1,
  legend.shift   = c(-5,-20),
  fname          = sprintf(out_file_base, 'cbh')
)

## SBH Map
gbd_map(
  data           = sbh_data,
  limits         = limits,
  labels         = labels,
  col            = 'YlGnBu',
  col.reverse    = FALSE,
  na.color       = "#999999",
  title          = sprintf("Number of SBH data sources by country, %s-%s",
                           start_year,end_year),
  legend.title   = "SBH Sources",
  legend.cex     = 1,
  legend.columns = 1,
  legend.shift   = c(-5,-20),
  fname          = sprintf(out_file_base, 'sbh')
)
