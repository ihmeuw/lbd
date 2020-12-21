####################################################################################################
## Collapse script for oral rehydration therapy (ORT) indicators
####################################################################################################
## Description: 
##  -Setup
##  -Load & clean data
##  -Data checks for DOM and PNG surveys
##  -Drop observations without diarrhea or ORT information
##  -Map string variables
##  -Data checks for Niger 9522 survey
##  -Data checks for IRQ and MEX surveys
##  -Exclude surveys identified as problematic during data vetting
##  -Generate all binary indicator variables
##  -Collapse to points & polygons
##  -Perform RHF definition crosswalk
##  -Format data & save
####################################################################################################
## Mutually-exclusive and collectively exhaustive (sum to 1) indicator list:
# 1) any_ors
# 2) rhf_only
# 3) no_ort
####################################################################################################


## Setup -------------------------------------------------------------------------
message('Setting up environment...')
rm(list=ls())

# set arguments
indics <- c('any_ors', 'rhf_only', 'no_ort')
sing <- T # running on RStudio singularity image? T/F
repo <- '<<<< FILEPATH REDACTED >>>>'
indicator_repo <- '<<<< FILEPATH REDACTED >>>>'
setwd(repo)
prev_version <- as.Date('', '%m-%d-%Y') # previous output data version (i.e., last time this code was run)

# get input version from most recently modified data file
files <- file.info(list.files('<<<< FILEPATH REDACTED >>>>', pattern = '*.feather'))
files <- rownames(files)
files <- gsub('.feather', '', files)
input_version <- max(files)

# load libraries
if (Sys.info()[1] == 'Windows') stop('This must be run on the cluster')
package_list <- readLines('<<<< FILEPATH REDACTED >>>>/package_list.csv')

# load functions
message('Loading in required R packages and MBG functions')
source(paste0(repo, '/setup.R'))
mbg_setup(package_list = package_list, repos = repo)
library('feather')


## Load & clean data -------------------------------------------------------------------------

# load data
message('Loading input dataset...')
data <- read_feather(paste0('<<<< FILEPATH REDACTED >>>>', input_version, '.feather'))
data <- data.table(data)
message(paste0('Input data starts with ', nrow(data), ' rows.'))

# clean data
message('Cleaning up dataset...')

# rename variables
setnames(data, c('iso3', 'survey_series', 'lat', 'long', 'year_start'), c('country', 'source', 'latitude', 'longitude', 'year'))

# fix variable classes
num_cols <- c('year', 'point', 'latitude', 'longitude')
data[, (num_cols) := lapply(.SD, function(x){as.numeric(x)}), .SDcols=num_cols]

# drop surveys started prior to 2000
message(paste0('Dropping ', nrow(data[year < 2000]), ' rows where the survey started prior to the year 2000.'))
data <- data[year >= 2000]

# drop rows missing age, or with age over 5
message(paste0('Dropping ', nrow(data[is.na(age_year) & is.na(age_month)]), ' rows with no age info.'))
data <- data[!is.na(age_year) | !is.na(age_month),]
message(paste0('Dropping ', nrow(data[age_year >= 5 & (age_month >= 60 | is.na(age_month))]), ' rows with neither age_year/age_month under 5y/60mo.'))
data <- data[age_year < 5 | age_month < 60,]
message(paste0('Dropping ', nrow(data[age_year >= 5 & age_month < 60]), ' rows with age_month under 60 but age_year over 5.'))
data <- data[!(age_year >=5 & age_month < 60),]
message(paste0('Dropping ', nrow(data[age_month >= 60 & age_year < 5]), ' rows with age_year under 5 but age_month over 60.'))
data <- data[(age_month < 60 | is.na(age_month)) & age_year < 5]

# drop observations with missing survey weights
if ('hhweight' %in% names(data)) {
  message(paste0('Dropping ', nrow(data[is.na(pweight) & is.na(hhweight)]), ' rows with no survey weight info.'))
  data <- data[!is.na(pweight) | !is.na(hhweight),]
  # assign hhweight to pweight if pweight is missing
  message(paste0('Setting ', nrow(data[is.na(pweight)]), ' rows with no pweight to hhweight.'))
  data[is.na(pweight), pweight := hhweight]
} else {
  message(paste0('Dropping ', nrow(data[is.na(pweight)]), ' rows with no survey weight info.'))
  data <- data[!is.na(pweight),]
}

# create new point column because the one in geocodebooks is unreliable
data[, point := NULL]
message(paste0('Assigning ', nrow(data[!is.na(latitude) & !is.na(longitude) & latitude != '' & longitude != '']), ' rows with lat/longs to point data.'))
data <- data[(!is.na(latitude) & !is.na(longitude)) & latitude != '' & longitude != '', point := 1]
message(paste0('Assigning ', nrow(data[is.na(latitude) & is.na(longitude) & shapefile != '' & location_code != '' & !is.na(shapefile) & !is.na(location_code)]), ' rows with shapefile/location code to polygon data.'))
data <- data[is.na(latitude) & is.na(longitude) & shapefile != '' & location_code != '' & !is.na(shapefile) & !is.na(location_code), point := 0]

# drop rows missing geographic information
message(paste0('Dropping ', nrow(data[point != 1 & point != 0]), ' rows not yet geo-matched.'))
data <- data[point == 1 | point == 0]


## Check Dominican Republic survey NID 77819 and NID 165645 to see whether points close together report similar diarrhea prevalences ---------------------------------------
if (nrow(data[nid == 165645]) > 0) {
  dom <- data[nid == 77819 | nid == 165645]
  dom <- dom[!is.na(had_diarrhea)]
  dom[, cluster_n := .N, by = geospatial_id]
  dom[, cluster_dia := sum(had_diarrhea), by = geospatial_id]
  dom[, dia_prop := cluster_dia/cluster_n]
  dom1 <- dom[latitude > 18.3 & latitude < 18.4 & longitude > -71.2 & longitude < -71.1]
  dom2 <- dom[latitude > 19.7 & latitude < 19.8 & longitude > -70.7 & longitude < -70.6]
  dom3 <- dom[latitude > 18.5 & latitude < 18.6 & longitude > -69.4 & longitude < -69.3]
  print(unique(dom1[, c('nid', 'geospatial_id', 'dia_prop', 'cluster_n', 'cluster_dia')]))
  print(unique(dom2[, c('nid', 'geospatial_id', 'dia_prop', 'cluster_n', 'cluster_dia')]))
  print(unique(dom3[, c('nid', 'geospatial_id', 'dia_prop', 'cluster_n', 'cluster_dia')]))
  message('No evidence to exclude survey 165645 as adjascent clusters show similar variation and ranges to survey 77819.')
  rm(dom1); rm(dom2); rm(dom3); rm(dom)
}


## Deal with weird PNG survey that has NA's instead of 0's where child did not recieve ORS or RHF ---------------------------------------
if (nrow(data[nid == 44870]) > 0) {
  data[nid == 44870 & diarrhea_tx == 1 & is.na(diarrhea_ors), diarrhea_ors := 0]
  data[nid == 44870 & diarrhea_tx == 1 & is.na(diarrhea_rhf), diarrhea_rhf := 0]
}


## Drop observations without diarrhea and missing any ORT information ---------------------------------------------------------------------------

# drop rows where child did not have diarrhea
message(paste0('Dropping ', nrow(data[had_diarrhea != 1]), ' rows where the child did not have diarrhea.'))
data <- data[had_diarrhea == 1]

# check rows with diarrhea, treatment, and ort information
message(paste0('There are ', nrow(data), ' observations that had diarrhea.'))
message(paste0('There are ', nrow(data[!is.na(diarrhea_tx)]), ' observations that asked about treatment and ', 
               nrow(data[is.na(diarrhea_tx)]), ' observations that are missing treatment info.'))
message(paste0('There are ', nrow(data[!is.na(diarrhea_ors) | !is.na(diarrhea_rhf) | !is.na(diarrhea_ort) |!is.na(diarrhea_ort_string_mapped) | !is.na(diarrhea_zinc) | !is.na(diarrhea_amount_eat_mapped)]), ' observations that asked about ORT-related info and ', 
               nrow(data[is.na(diarrhea_ors) & is.na(diarrhea_rhf) & is.na(diarrhea_ort) & is.na(diarrhea_ort_string_mapped) & is.na(diarrhea_zinc) & is.na(diarrhea_amount_eat_mapped)]), ' observations that are missing ORT-related info.'))

# drop rows missing ort information
message(paste0('Dropping ', nrow(data[is.na(diarrhea_ors) & is.na(diarrhea_rhf) & is.na(diarrhea_ort) & is.na(diarrhea_ort_string_mapped) & is.na(diarrhea_zinc) & is.na(diarrhea_amount_eat_mapped)]), ' rows with no ORT-related info.'))
data <- data[!is.na(diarrhea_ors) | !is.na(diarrhea_rhf) | !is.na(diarrhea_ort) | !is.na(diarrhea_ort_string_mapped) | !is.na(diarrhea_zinc) | !is.na(diarrhea_amount_eat_mapped)]


## Map all ORT string variables to binary variables -------------------------------------------------------------

# if any ORT string variables are extracted
if (!is.na(match('diarrhea_ort_string_mapped', names(data)))) {
  
  # create binary mapped variables
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_ors_mapped := as.numeric(diarrhea_ort_string_mapped %like% 'ors')]
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_rhf_mapped := as.numeric(diarrhea_ort_string_mapped %like% 'rhf')]
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_zinc_mapped := as.numeric(diarrhea_ort_string_mapped %like% 'zinc')]
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_ort_mapped := as.numeric(diarrhea_ort_string_mapped %like% 'ort')]
  
  # check the types of variables extracted per study (ideally we don't extract both string and binary variables)
  message(paste0('There are ', nrow(data[!is.na(diarrhea_ors) & !is.na(diarrhea_ors_mapped)]), ' rows with both binary and string variables extracted for ors. Make sure only one type is extracted per survey.'))
  message(paste0('There are ', nrow(data[!is.na(diarrhea_rhf) & !is.na(diarrhea_rhf_mapped)]), ' rows with both binary and string variables extracted for rhf. Make sure only one type is extracted per survey.'))
  message(paste0('There are ', nrow(data[!is.na(diarrhea_zinc) & !is.na(diarrhea_zinc_mapped)]), ' rows with both binary and string variables extracted for zinc. Make sure only one type is extracted per survey.'))
  message(paste0('There are ', nrow(data[!is.na(diarrhea_ort) & !is.na(diarrhea_ort_mapped)]), ' rows with both binary and string variables extracted for ort. Make sure only one type is extracted per survey.'))
  
  # combine into one variable
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_ors := as.numeric(diarrhea_ors_mapped == 1)]
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_rhf := as.numeric(diarrhea_rhf_mapped == 1)]
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_zinc := as.numeric(diarrhea_zinc_mapped == 1)]
  data[!is.na(diarrhea_ort_string_mapped), diarrhea_ort := as.numeric(diarrhea_ort_mapped == 1)]
}


## Check whether the Niger survey NID 9522 ORT question encompasses ORS and RHF completely ---------------------------------------
if (nrow(data[nid == 9522]) > 0) {
  # subset data by study and ort variable
  niger <- data[nid == 9522 & diarrhea_ort == 1]
  # determine number of ort that aren't designated as ORS or RHF
  num_unclear_ort <- nrow(niger[diarrhea_ors == 0 & diarrhea_rhf == 0])
  # notify if study will be excluded because ORT can't be separated into ORS and RHF confidently for each child
  if (is.na(unique(niger$diarrhea_ors))) {
    message('Niger survey 9522 has not yet been properly extracted. It will be excluded for now.')
    data <- data[nid != 9522]
    } else {
      if (num_unclear_ort > 0) {
        message(paste0(num_unclear_ort, ' children with diarrhea in Niger survey 9522 received "ORT", but their treatment was not designated as either "ORS" or "RHF".'))
        message('We will set all treated with either ORS or RHF as ORT=1, set ORS=NA and RHF=NA, and only use this survey to model ORT (ORS or RHF).')
        data[nid == 9522 & diarrhea_ors == 1, diarrhea_ort := 1]
        data[nid == 9522 & diarrhea_rhf == 1, diarrhea_ort := 1]
        data[nid == 9522, c('diarrhea_ors', 'diarrhea_rhf') := NA]
      } else {message('For each child with diarrhea in Niger survey 9522 who received "ORT", their treatment was designated as either "ORS" or "RHF".') 
              message('This survey will be included and the "ORT" variable can be ignored.')}
    }
  rm(niger)
}


## Check Mexico survey NID 8618 and Iraq survey NID 23565 with multiple binary variables for ORT ---------------------------------------
if (nrow(data[nid == 8618]) > 0) {
  mex <- data[nid == 8618]
  message('2005 Mexico survey 8618 with multiple ORT variables reported:')
  message(paste0(nrow(mex[had_diarrhea == 1]), ' children with diarrhea'))
  message(paste0(round(nrow(mex[diarrhea_ors == 1])/nrow(mex)*100,0), '% of children with diarrhea received ORS'))
  message(paste0(round(nrow(mex[diarrhea_rhf == 1])/nrow(mex)*100,0), '% of children with diarrhea received RHF'))
  message(paste0(round(nrow(mex[diarrhea_ors == 1 | diarrhea_rhf == 1])/nrow(mex)*100,0), '% of children with diarrhea received either ORS or RHF'))
  message('These numbers are similar to numbers reported in 2011. We will include survey 8618 for now.')
  rm(mex)
}
if (nrow(data[nid == 23565]) > 0) {
  irq <- data[nid == 23565]
  message('2004 Iraq survey 23565 with multiple ORT variables reported:')
  message(paste0(nrow(irq), ' children with diarrhea'))
  message(paste0(round(nrow(irq[diarrhea_ors == 1])/nrow(irq)*100,0), '% of children with diarrhea received ORS'))
  message(paste0(round(nrow(irq[diarrhea_rhf == 1])/nrow(irq)*100,0), '% of children with diarrhea received RHF'))
  message(paste0(round(nrow(irq[diarrhea_ors == 1 | diarrhea_rhf == 1])/nrow(irq)*100,0), '% of children with diarrhea received either ORS or RHF'))
  message('These numbers are ~8-fold lower than those from MICS surveys in 2000 and 2006. We will exclude survey 23565 for now.')
  data <- data[nid != 23565]
  rm(irq)
}


## Exclude additional surveys with problematic, inaccurate, or duplicate data --------------------------------------------------------------------------------

# Identified for exclustion through systematic data vetting
# Documented here: https://docs.google.com/spreadsheets/d/1z44LQiUzrDHtWIIBJraARgSXvtamJBBgVoiSSyfKwl4/edit?ts=5b3feb4f#gid=1571196334)

# Exclude Dominican Republic NID 21198 because it is a survey that specifically sampled Batey villages and we only have polygon, not point, data
if (nrow(data[nid == 21198]) > 0) data <- data[nid != 21198]

# Exclude Dominican Republic NID 627 (also named as 265014) because it is a duplicate of NID 200697
if (nrow(data[nid == 627]) > 0) data <- data[nid != 627]
if (nrow(data[nid == 265014]) > 0) data <- data[nid != 265014]

# Exclude Bangladesh NID 942 because it has problematic case definition for diarrhea
if (nrow(data[nid == 942]) > 0) data <- data[nid != 942]

# Exclude Benin NID 79839 because it has known quality issues
if (nrow(data[nid == 79839]) > 0) data <- data[nid != 79839]

# Exclude longitudinal surveys that sample the same children more than once
if (nrow(data[nid == 286657]) > 0) data <- data[nid != 286657]
if (nrow(data[nid == 93848]) > 0) data <- data[nid != 93848]
if (nrow(data[nid == 58419]) > 0) data <- data[nid != 58419]
if (nrow(data[nid == 142934]) > 0) data <- data[nid != 142934]
if (nrow(data[nid == 160781]) > 0) data <- data[nid != 160781]
if (nrow(data[nid == 81005]) > 0) data <- data[nid != 81005]
if (nrow(data[nid == 81004]) > 0) data <- data[nid != 81004]
if (nrow(data[nid == 224096]) > 0) data <- data[nid != 224096]
if (nrow(data[nid == 142935]) > 0) data <- data[nid != 142935]
if (nrow(data[nid == 249499]) > 0) data <- data[nid != 249499]

# Exclude Niger NID 94140 that souble counts children in two seasons
if (nrow(data[nid == 94140]) > 0) data <- data[nid != 94140]

# Exclude the Angola NID 687 because it is known to have serious methodological issues and UNICEF doesn't use it in its ORS estimates
if (nrow(data[nid == 687]) > 0) data <- data[nid != 687]

# Exclude Guatemala NID 352625 because it only allows one treatment to be listed and vastly underestimates ORS
if (nrow(data[nid == 352625]) > 0) data <- data[nid != 352625]

# Exlude Angola 2011 survey because we have no denominator
if (nrow(data[nid == 151568]) > 0) data <- data[nid != 151568]

# Exclude Burkina Faso 2005 survey 22950 and 2017 survey 37530 because age groups are very non-representative
if (nrow(data[nid == 22950]) > 0) data <- data[nid != 22950]
if (nrow(data[nid == 375030]) > 0) data <- data[nid != 375030]

# Exclude Kenya 2017 survey 375035 because age group looks very non-representative
if (nrow(data[nid == 375035]) > 0) data <- data[nid != 375035]

# Exclude Malawi surveys 224223 and 93806 because they are longitudinal and data look very off from DHS and MICS surveys in those years
if (nrow(data[nid == 224223]) > 0) data <- data[nid != 224223]
if (nrow(data[nid == 93806]) > 0) data <- data[nid != 93806]

# Exclude South Africa survey 20798 because we don't have permission to use it yet
if (nrow(data[nid == 20798]) > 0) data <- data[nid != 20798]

# Exclude Brazil data because it is not below admin 1
if(nrow(data[nid == 141948]) > 0) data <- data[nid != 141948]

# Exclude Syria and Jamaica surveys that have unreasonably low and high ORS estimates, respectively, and are also missing a lot of ORS information
if(nrow(data[nid == 126911]) > 0) data <- data[nid != 126911]
if(nrow(data[nid == 7140]) > 0) data <- data[nid != 7140]


## Generate mutually-exclusive and collectively exhaustive indicator variables -----------------------------------------------------------------

# Any ORS
if('any_ors' %in% indics) {
  # remove missing values
  any_ors <- data[!is.na(diarrhea_ors) & !is.na(diarrhea_rhf)]
  # create variable
  any_ors <- any_ors[, value := ifelse(diarrhea_ors == 1, 1, 0)]
  # store list of NIDs for comparison at the end
  any_ors_nids <- unique(any_ors$nid)
}
# RHF only
if('rhf_only' %in% indics) {
  # remove missing values
  rhf_only <- data[!is.na(diarrhea_ors) & !is.na(diarrhea_rhf)]
  # create variable
  rhf_only <- rhf_only[, value := ifelse(diarrhea_ors == 0 & diarrhea_rhf == 1, 1, 0)]
  # store list of NIDs for comparison at the end
  rhf_only_nids <- unique(rhf_only$nid)
}
# no treatment
if('no_ort' %in% indics) {
  # remove missing values
  no_ort <- data[!is.na(diarrhea_ors) & !is.na(diarrhea_rhf)]
  # create variable
  no_ort <- no_ort[, value := ifelse(diarrhea_ors == 0 & diarrhea_rhf == 0, 1, 0)]
  # store list of NIDs for comparison at the end
  no_ort_nids <- unique(no_ort$nid)
}


## Collapse to points & polygons --------------------------------------------------------------------
collapse <- function(i) {
  message('Collapsing ', i)
  dt <- get(i)
  # collapse
  byvars <- 'nid,country,source,year,point,shapefile,location_code,latitude,longitude'
  dt <- dt[, list(int_year = floor(median(int_year, na.rm=T)), 
                  mean = weighted.mean(value, pweight),
                  N_obs = .N,
                  N = sum(pweight)^2/sum(pweight^2),
                  sum_of_sample_weights = sum(pweight)),
           by = byvars]
  # set year to be the weighted mean of median interview year by cluster
  dt[, year := NULL]
  dt[, year := round(weighted.mean(int_year, N, na.rm = T)), by = 'nid']
  dt[, int_year := NULL]
  # create cluster id
  dt[, cluster_id := .I]
  # return data table
  return(dt)
}
alldata <- lapply(indics, collapse)
names(alldata) <- indics


## Set years manually for surveys with NA year after collapse --------------------------------------------------

for (i in indics) {
  
  # manually assign years
  if (nrow(alldata[[i]][nid == 286783]) > 0) alldata[[i]][nid == 286783, year := 2018]
  if (nrow(alldata[[i]][nid == 399853]) > 0) alldata[[i]][nid == 399853, year := 2018]
  if (nrow(alldata[[i]][nid == 18815]) > 0) alldata[[i]][nid == 18815, year := 2007]
  
  # check for other NAs
  if (nrow(alldata[[i]][is.na(year)]) > 0) stop('NA years are present in the data file after collapse.\nFix this before proceeding.')
}


## Perform RHF definition cross-walk --------------------------------------------------------------------

# loop over indicators that use RHF survey questions
for (i in c('rhf_only', 'no_ort')) {

  # perform crosswalk for each indicator if present
  if (i %in% indics) {
    message('\nCrosswalking ', i)

    # get copy of indicator data table
    coverage_data <- copy(alldata[[i]])

    # perform definition crosswalk
    source(paste0(indicator_repo, 'ors/2_data_prep/rhf_definition_crosswalk.R'))
    alldata[[i]] <- copy(coverage_data)

  }

}


## Save data for data coverage plots -----------------------------------------------------
save(alldata, file = paste0('<<<< FILEPATH REDACTED >>>>/alldata_', input_version, '.RData'))


## Convert to counts and save data for polygon resampling -----------------------------------------------------
countdata <- lapply(indics, function(i) {dt <- alldata[[i]][, mean := mean*N]; return(dt)})
names(countdata) <- indics
save(countdata, file = paste0('<<<< FILEPATH REDACTED >>>>/countdata_', input_version, '.RData'))