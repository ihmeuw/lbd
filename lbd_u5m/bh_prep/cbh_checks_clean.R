## ##############################################################
## Inputs: Child-Level CBH data prepared from Ubcov
## Outputs: CBH checks report pdf
## Description: This script creates some diagnostic variables
##   then runs the cbh checks markdown script to generate
##   a pdf - simplified version of cbh_collapse.R
##   Dependencies: data.table, Biograph, survival
## ##############################################################

## ##############################################################
## USER DEFINED GLOBALS FOR AGE BIN AND PERIOD TABULATION

## set file sytem paths
# where to pull in the ubcov'd CBH data in from
dat_dir  <- '<<<< FILEPATH REDACTED >>>>'

## ##############

## ##############################################################
## SET UP
## Loading Packages
library(data.table)
library(Biograph)
library(survival)

## Loading in child-level CBH data
d <- setDT(readRDS(sprintf("%<<<< FILEPATH REDACTED >>>>",dat_dir)))

# limiting to 1998 for now..
d <- subset(d, year >= 1998)

## ##############################################################
## CLEAN DATA

## Identify the smallest possible geography
d[, location_code := as.character(location_code)]
d$location_code[d$location_code=='NA']=''
d[, latnum        := as.numeric(latnum)]
d[, longnum       := as.numeric(longnum)]
d$point[is.na(d$point) & (d$location_code != "" & d$shapefile != "")] = 0
d$point[is.na(d$point) & (!is.na(d$latnum) & !is.na(d$longnum))] = 1
d[point==0, geo_unit := paste0(nid,location_code,shapefile)]
d[point==1 , geo_unit := paste0(nid,cluster_number)]

## identify rows with no geography
d[, no_geog := (is.na(latnum) | is.na(longnum)) & (shapefile == '' | location_code == '')]

## Standardize yearborn definition and check for missingness
d[, yrborn  := cmc_as_year(interview_date_cmc - birthtointerview_cmc)]
d[, yrborn2 := cmc_as_year(child_dob_cmc)]
d[is.na(yrborn), yrborn := yrborn2]

## Standardize year of interview definition
d[, yrint := cmc_as_year(interview_date_cmc)]

## Standardize aod definition
## Also, define variable 'alive', will get used as a censoring indicator
setnames(d,old='child_age_at_death_months',new='aod')
d[,alive := aod == 6000]  # if alive at time of survey aod is coded as 6000, so consider them alive
d[,alive := aod >= 60]   # Consider all those who died after 60 months as 'alive' for our purposes of studying under-5 mortality. Equivalent to being right-censored from this perspective
d$alive[is.na(d$child_alive)] <- NA # if child_alive is unknown in the ubcov data, keep as NA
d[,died  := alive == FALSE] # if not alive, then they have died

## Standardize child age definition
d[,childage := round((yrint - yrborn)*12)] # if alive, record their age
d$childage[d$died==TRUE & !is.na(d$died)] <- d$aod[d$died==TRUE & !is.na(d$died)] # if died, record their age at death
d[,childage := childage + 0.0001] # add a tiny bit to help with discrete age category definitions

## Point data will not be weighted, so we convert at point row weights to 1
d$weight[d$point == TRUE] <- 1

## if an entire survey is missing weights, then we give them a 1 weight.
#   This follows on the explicit choice to include surveys without weights as unweighted
#   The only non-weighted rows remaining after this will be polygon rows in surveys with
#   weights, theyll be dropped later on
d[nid==26661, weight := 1] # for this eth GF survey for some reason the weights are all 0.. 
nowgt <- d[,.(pct_unweighted = sum(is.na(weight))/.N), by = nid]
nowgt <- subset(nowgt, pct_unweighted == 1)$nid
if(length(nowgt) > 0) {
  d$weight[d$nid %in% nowgt] <- 1
}

## Generate variable for rows to be dropped
d[, must_drop := is.na(weight) | is.na(childage) | childage > 41*12 | is.na(point) |
    childage < 0 | no_geog == TRUE | is.na(aod) | is.na(alive) | is.na(yrborn)]

#generate cbh checks report
rmarkdown::render("<<<< FILEPATH REDACTED >>>>", params = list(
  dt = d,
  outlier_path = "<<<< FILEPATH REDACTED >>>>"))
