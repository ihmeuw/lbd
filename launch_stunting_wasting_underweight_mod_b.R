###########################################################################################
###########################################################################################
## Run MBG model for proportion of moderate and severe stunting children using binomial cluster data
## written by AOZ
##
###########################################################################################
###########################################################################################

## Set repo location and indicator group
repo            <- <<<< FILEPATH REDACTED >>>>>
indicator_group <- 'child_growth_failure'
indicator       <- 'stunting_mod_b'
Regions         <- c('cssa','essa','name','sssa','wssa') ## could be put in config

## set repo and source the launch script
setwd(repo)
source('./child_growth_failure/launch.R')


## Set repo location and indicator group
repo            <- <<<< FILEPATH REDACTED >>>>>
indicator_group <- 'child_growth_failure'
indicator       <- 'wasting_mod_b'
Regions         <- c('cssa','essa','name','sssa','wssa') ## could be put in config

## set repo and source the launch script
setwd(repo)
source('./child_growth_failure/launch.R')



## Set repo location and indicator group
repo            <- <<<< FILEPATH REDACTED >>>>>
indicator_group <- 'child_growth_failure'
indicator       <- 'underweight_mod_b'
Regions         <- c('cssa','essa','name','sssa','wssa') ## could be put in config

## set repo and source the launch script
setwd(repo)
source('./child_growth_failure/launch.R')

