#######################
rm(list = ls())
#required things
iso3 <- "ETH"
country <- "Ethiopia"
#input the treat_extract fi
input_data <- paste0("<<<< FILEPATH REDACTED >>>>")
username <- "USER"
source(paste0("<<<< FILEPATH REDACTED >>>>/dcp_functions.R"))
#kwargs for things
## for plot_time_series function
unaids_compare = T
adm2_compare = T
first_yr = 2000
last_yr = 2019

#######################

data <- standardize_file(input_file=input_data)
final_dt <- art_rate_1(iso3, username, data)
# you can just pull UNAIDS numbers without graphing pretty easily
UNAIDS <- pull_UNAIDS(iso3) #unaids year needed
# Otherwise it'll get calculated in the plotting function to be used for the comparison plot.
p_timeseries <- plot_time_series(country, iso3, final_dt, unaids_compare = T,
                                 adm2_compare = T, first_yr = 2000, last_yr = 2019)

test <- perc_unit_covered(iso3, data)
launch_dcps(username, final_dt, test, iso3, plot_type = "alive_ART")
launch_dcps(username, final_dt, test, iso3, plot_type = "art_rate")

