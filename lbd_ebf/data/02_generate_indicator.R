####################################################################################################
## Script for generating breastfeeding indicators
##  Part 1: Setup, load & clean data
##  Part 2: Generate indicator variables
####################################################################################################

message("Setting up environment...")
rm(list = ls())

#set arguments 

input_version  <- "<<<< FILEPATH REDACTED >>>>"
core_repo      <- paste0("<<<< FILEPATH REDACTED >>>>") 
indics         <- "ebf"
commondir      <- "<<<< FILEPATH REDACTED >>>>"
indic_repo     <- "<<<< FILEPATH REDACTED >>>>"
stg <- "1"
stage <- ifelse("3" %in% stg, "stage3", ifelse((("2a" %in% stg) | ("2b" %in% stg)), "stage2", ""))
start_year <- 1998 
keepyr <- FALSE
message("Loading packages...")
package_list <- readLines("<<<< FILEPATH REDACTED >>>>")

source("<<<< FILEPATH REDACTED >>>>")
mbg_setup(package_list = package_list, repos = c(indic_repo, core_repo))

## Loading and cleaning data--------------------------------------------------------------------------

data <- fread("<<<< FILEPATH REDACTED >>>>", format.Date(input_version, "%Y_%m_%d"), ".csv")
data <- data[, .(nid, iso3, survey_series, survey_name, year_start, strata, psu,
                 pweight, hhweight, point, shapefile, location_code, lat, long,
                 int_month, int_year, sex_id, age_year, age_month, bf_ever, bf_dur, cur_bf,
                 bf_still, food_24h, liq_24h, food_24h_list, liq_24h_list, bf_excl_dur, 
                 bf_excl, bf_init, ebf_first3, water, other_milk, other_foods)] 

for (ii in 1:ncol(data)) attributes(data[[ii]]) <- NULL

setnames(data, c("iso3", "survey_series", "lat", "long", "year_start"), c("country", "source", "latitude", "longitude", "year"))
num_cols <- c("year", "point", "latitude", "longitude", "location_code")
data[, (num_cols) := lapply(.SD, function(x){as.numeric(x)}), .SDcols = num_cols]
n_surv <- length(unique(data[!is.na(nid), nid]))

locs <- fread("<<<< FILEPATH REDACTED >>>>")
stg_list <- locs[Stage %in% as.character(stg), iso3]
data[, country := substr(country, 1, 3)] 
data <- data[country %in% stg_list,]
nid_list_all <- unique(data$nid)

data <- data[!is.na(age_year) | !is.na(age_month),]
data <- data[!is.na(pweight) | !is.na(hhweight),]
data[!is.na(hhweight) & is.na(pweight), pweight := hhweight]
data[, point := as.integer(!is.na(latitude) & !is.na(longitude))]
data <- data[(!is.na(shapefile) & !is.na(location_code)) | point == 1,]
nid_list_all_indic <- unique(data$nid)
data <- data[ year >= start_year, ]

message(paste(c("Generating indicators:", indics), collapse=" "))

if("ebf" %in% indics) {
  message("Generating indicator ebf")
  ebf <- if(keepyr) data[age_year < 0.5 ] else data[age_month < 6,] 
  
  ebf[, askstill := ifelse((length(na.omit(bf_still)) > 0) | (length(na.omit(cur_bf)) > 0), 1, 0), by = c("nid")]
  ebf <- ebf[askstill == 1,]
  ebf[, hasvars := as.numeric((!is.na(bf_still) | !is.na(cur_bf)) & !is.na(food_24h) & !is.na(liq_24h)), by = c("nid")] 
  ebf <- ebf[hasvars == 1,]
  
  ebf[, value := 0]
  ebf[bf_still == 1 & food_24h == "0" & liq_24h == "0", value := 1]
  
  write.csv(ebf, file = paste0("<<<< FILEPATH REDACTED >>>>"), row.names=F)
}
