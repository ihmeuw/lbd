## CMR Shell ##
rm(list = ls())
libs <- c("data.table", "raster", "foreign", "rgdal", "geosphere", "colorspace",
          "dplyr", "rgeos", "car","plyr", "ggplot2", "scales", "sf", "gridExtra",
          "tidyr", "readxl", "reticulate", "stringr")
sapply(libs, require, character.only = T)

#######################
#required things
iso3 <- "CMR"
country <- "Cameroon"
input_data <- "<<<< FILEPATH REDACTED >>>>"
username <- "USER"

#kwargs for things
## for plot_time_series function
unaids_compare = T
adm2_compare = T
first_yr = 2000
last_yr = 2019

#######################

source("<<<< FILEPATH REDACTED >>>>/dcp_functions.R")
data <- standardize_file(input_file=input_data)
final_dt <- art_rate_1(iso3, username, data)
# you can just pull UNAIDS numbers without graphing pretty easily
UNAIDS <- pull_UNAIDS(iso3) #unaids year needed
# Otherwise it'll get calculated in the plotting function to be used for the comparison plot.
p_timeseries <- plot_time_series(country, iso3, final_dt, unaids_compare = T,
                                 adm2_compare = T, first_yr = 2000, last_yr = 2019)

test <- perc_unit_covered(country, iso3, final_dt)
launch_dcps(username, final_dt, test, iso3, plot_type = "alive_ART")
launch_dcps(username, final_dt, test, iso3, plot_type = "art_rate")
