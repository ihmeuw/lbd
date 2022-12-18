####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, ANC data
####################################################################################################

## Setup -------------------------------------------------------------------------------------------
rm(list = ls())

# set arguments for this run and user.
core_repo <- paste0(<<<< FILEPATH REDACTED >>>> '/lbd_core/')
indic_repo <- paste0(<<<< FILEPATH REDACTED >>>>'/lbd_hiv/')
cores <- 10
indicator <- 'hiv_prev_disagg'
shapefile_version <- 'current'
modeling_shapefile_version<-shapefile_version

# find most recent versions
most_recent_date <- function(dir, date_format = "%Y_%m_%d", out_format = "%Y_%m_%d", file_pattern = NULL) {
  date_pattern <- gsub("%y|%m|%d", "[[:digit:]]{2}", gsub("%Y", "[[:digit:]]{4}", date_format))
  dates <- dir(dir, pattern = date_pattern)
  if (!is.null(file_pattern)) dates <- grep(file_pattern, dates, value = T)
  dates <- gsub(paste0("(.*)(", date_pattern, ")(.*)"), "\\2", dates)
  dates <- as.Date(dates, date_format)
  format(max(dates), out_format)
}

geomatch_date <- most_recent_date(<<<< FILEPATH REDACTED >>>>)
anc_date <- most_recent_date(<<<< FILEPATH REDACTED >>>>)
prev_collapse_date <- most_recent_date(<<<< FILEPATH REDACTED >>>>, file_pattern = paste0(indicator, "_"))
collapse_date <- format(Sys.time(), "%Y_%m_%d")
report_commit <- system(paste0('cd ', indic_repo, ';git log -1 --format="%h"'), intern = T)

# load libraries & functions
source(paste0(core_repo, '/mbg_central/setup.R'))
package_list <- readLines(paste0(core_repo, "/mbg_central/share_scripts/common_inputs/package_list.csv"))
mbg_setup(package_list = c(package_list, 'openxlsx'), repos = c(core_repo, indic_repo))

setompthreads(cores)

## Load and subset data ----------------------------------------------------------------------------
# load geomatched data
data <- readRDS(paste0(<<<< FILEPATH REDACTED >>>>,
                       format.Date(as.Date(geomatch_date, "%Y_%m_%d"), "%Y_%m_%d"), "/",
                       geomatch_date, ".RDS"))
data <- data.table(data)

# drop unneeded variables and remove annoying stata attributes
data <- data[, list(nid, country, survey_series, survey_name, year,
                    strata, psu, point, shapefile, location_code, latitude, longitude,
                    sex_id, age_year, int_year, hiv_test, hiv_weight)]
for (ii in 1:ncol(data)) attributes(data[[ ii]]) <- NULL

# rename variables
setnames(data, c("survey_series"), c("source"), skip_absent=TRUE)

# fix variable class
data[, point := as.numeric(point)]
data[, latitude := as.numeric(latitude)]
data[, longitude := as.numeric(longitude)]

# subset to ages 15-59 
data <- data[between(age_year, 15, 59),]
data <- data[!is.na(age_year),]

# drop observations with missing hiv test result or hiv weight
data <- data[!is.na(hiv_test) & !is.na(hiv_weight),]

# drop PSUs where all hiv weights are 0 (otherwise this causes NaNs)
drop <- data[, sum(hiv_weight), by = 'nid,psu'][V1 == 0, list(nid, psu)]
data <- data[!drop, on = c('nid', 'psu')]

# drop point data with missing latitude/longitude
data <- data[(!is.na(latitude) & !is.na(longitude)) | point == 0,]

# drop polygon data with missing shapefile/location code
data <- data[(!is.na(shapefile) & !is.na(location_code)) | point == 1,]

# drop point data with missing latitude/longitude
drop_pt_polys <- data[(!is.na(latitude) & !is.na(longitude) & point == 0),]
if(nrow(drop_pt_polys) > 0) message('Warning: some rows are listed at polygons but contain longitudes & latitudes. Please inspect')
data <- data[!(!is.na(latitude) & !is.na(longitude) & point == 0),]

# drop countries outside Africa--Will want to change when non-african countries are included
loc <- source("<<<< FILEPATH REDACTED >>>>/get_location_metadata.R")
loc <- data.table(get_location_metadata(location_set_id = 2, gbd_round_id = 4))
loc1 <- loc[grepl("Africa", region_name) & location_type == "admin0", ihme_loc_id]

#include select Non-african locations
loc2 <-   loc[ihme_loc_id== 'MAR' | ihme_loc_id== 'KHM' |
                ihme_loc_id== 'IND' | ihme_loc_id== 'HTI' |
                ihme_loc_id== 'DOM' | ihme_loc_id== 'VNM' |  
                ihme_loc_id== 'PNG', ihme_loc_id]

data <- data[country %in% loc1 | country %in%  loc2,]

# Drop North Africa Data
north_africa <- c(
  "ESH", "LBY","EGY","DZA","TUN" )
data <- subset(data, !(country %in% north_africa))

#Function to round age_year to 5-year age group
mfloor <- function(age,base=5){ 
  base*floor(age/base) 
} 
data[,round_age:=age_year]
data<-data[, round_age:=mfloor(age_year)]

#Convert round_age to categorized bins 1:9
data[,agebin:=as.factor(round_age)]
levels(data$agebin)<-c(1:9)
data$agebin<-as.numeric(data$agebin)

# store NIDs and sample size to compare at the end
nid_list <- unique(data$nid)
start_n <- data[, list(start = .N), by = 'country,nid']

## Collapse all data (polygon and point) -----------------------------------------------------------
# collapse data by survey and location, age-sex disaggregated
data_as_disagg <- data[, list(int_year = floor(median(int_year, na.rm = T)),
                              hiv_test = weighted.mean(hiv_test, hiv_weight),
                              N = sum(hiv_weight) ^ 2 / sum(hiv_weight ^ 2),
                              N_obs = .N,
                              sum_of_sample_weights = sum(hiv_weight)),
                       by = 'nid,country,source,year,point,shapefile,location_code,latitude,longitude,agebin,sex_id']

#----------------------------------------------------------------------------------------------------------
#Read in report data for inclusion in disagg_ag_data
report <- fread(<<<<FILEPATH REDACTED >>>>)
report <- report[, list(country, nid, location_type, survey_series, survey_name, int_year,
                        start_age, end_age, sex_id, point, shapefile, location_code, latitude,
                        longitude, hiv_test, hiv_test_lb, hiv_test_ub, N)]

report[, id := paste(nid, country, int_year)]
#For now take out the sex-aggregated reports
report <- report[sex_id < 3]
# drop national estimates, countries outside Africa, and NIDs we have microdata for
report <- report[location_type != "National",]
report <- report[country %in% loc1 | country %in% loc2]
report <- report[!nid %in% data$nid]
report <- report[nid!=414568]
report <- report[nid!=409506]
report <- report[, hiv_test := as.numeric(hiv_test)]
report <- report[, N := as.numeric(N)]
# change hiv_test from a percentage to a proportion
report[, `:=` (hiv_test = hiv_test/100, hiv_test_lb = hiv_test_lb/100, hiv_test_ub = hiv_test_ub/100)]

# estimate N if only CIs are available (using a normal approximation, assuming 95% CIs)
report[is.na(N) & !is.na(hiv_test_lb) & !is.na(hiv_test_ub),
       N := (hiv_test * (1 - hiv_test)) / (((hiv_test_ub - hiv_test_lb) / (2 * 1.96))^2)]
ci0 <- report[is.na(N) & hiv_test_lb == 0 & hiv_test_ub == 0, .N]
if (ci0 > 0) warning(paste("Report data dropped because no sample size was provided and lower and upper bound are both 0 for", ci0, "source-country-year-locations"))
report <- report[!(is.na(N) & hiv_test_lb == 0 & hiv_test_ub == 0),]
report[, c("hiv_test_lb", "hiv_test_ub") := NULL]

# drop alternate ages if 15 to 49 is available
report <- report[start_age %in% c(15, 20, 25, 30, 35, 40, 45, 50, 55) & 
                   end_age %in% c(19, 24, 29, 34, 39, 44, 49, 54, 59)]

report[, type := "Survey report"]
setnames(report, 'int_year', 'year')

report <- report[,list(country, nid, point, shapefile, location_code, latitude, longitude, year,   
                       sex_id, start_age, end_age, hiv_test,N, survey_name, survey_series, type )]

ag_data <- copy(report)
ag_data[, N := as.numeric(N)]
setnames(ag_data, 'survey_series', 'source', skip_absent=TRUE)
rm(report)


## Add ANC data ------------------------------------------------------------------------------------
# load anc data
anc <- readRDS(<<<< FILEPATH REDACTED >>>>)
anc <- data.table(mutate_if(anc, is.factor, as.character))
setnames(anc, "anc_hiv_test", "hiv_test")

anc <- anc[country != 'VNM']

# remove potentially problematic non-standard characters in site names
for (i in 1:nrow(anc)) {
  if (class(try(nchar(anc[i, site]), silent = T)) == 'try-error') {
    prob <- which(sapply(1:50, function(x) class(try(substr(anc[i, site], x, x), silent = T))) == 'try-error')
    anc[i, site := paste0(substr(site, 1, prob - 1), 'XX')]
  }
}

# combine with survey data, assign other values
anc[, sum_of_sample_weights := N]
anc[, type := "ANC"]
anc[, sex_id := 2]
anc[, agebin := NA]


##--------------------------------------------------------------------------
###Drop outlier ANC data--determined following model fitting. Review when adding new data 
##----------------------------------------------------------------------
print(paste0('before dropping outliers, nrow ANC==' , nrow(anc)))
anc <- anc[!(country=='BDI' & year==2006)]
anc <- anc[!(country=='ERI' & year %in% c(2003, 2005) & site=='Asseb')]
anc <- anc[!(country=='KEN' & year==2000 & site %in% c('Mbale', 'Meru', 'Thika'))]
anc <- anc[!(country=='KEN' & year==2001 & site=='Fatima')]
anc <- anc[!(country=='UGA' & year==2011 & site=='Lira Hosp')]
anc <- anc[!(country=='BEN' & year==2003 & site %in% c('Angaradébou', 'Kotopounga'))]
anc <- anc[!(country=='BEN' & year==2005 & site=='Comè')]
anc <- anc[!(country=='BEN' & year==2008 & site=='Abomey−Calavi / Godomey')]
anc <- anc[!(country=='CMR' & year==2009 & site=='Sangmelima')]
anc <- anc[!(country=='GHA' & year %in% c(2003, 2011) & site=='Cape Coast')]
anc <- anc[!(country=='GIN' & year==2004 & site=='GOUECKE')]
anc <- anc[!(country=='NGA' & year==2008 & site %in% c('Sokoto (Dogon Daji)', 'Yola', 'Abuja Bwari'))]

anc <- anc[!(country=='LSO' & year==2018)]


anc <- anc[!(country=='ETH' & year %in% c(2000, 2001) & site=='Adama HC')]
anc <- anc[!(country=='KEN' & year==2000 & site %in% c('Chulaimbo', 'Kisumu', 'UON4 - Baba Dogo'))]
anc <- anc[!(country=='MOZ' & year==2000 & site %in% c('9 - Ponta-Gea', '11 - Mondlane', '12 - No. 3'))]
anc <- anc[!(country=='TZA' & year==2000 & site %in% c('Mbeya (Itete)', 'Mbeya (Kiwanjampaka)', 'Mbeya (Kyela)', 'Mbeya(Ruanda)'))]
anc <- anc[!(country=='TZA' & year==2007 & site %in% c('Iringa (Ngome)', 'Vwawa Hospital'))]
anc <- anc[!(country=='SLE' & year==2006 & site %in% c('Makeni', 'Pujehun'))]
anc <- anc[!(country=='SLE' & year==2007 & site %in% c('PCMH', 'Kissy Urban centre', 'Koidu'))]
anc <- anc[!(country=='SLE' & year==2008 & site %in% c('Kambia', 'Mattru UBC'))]
anc <- anc[!(country=='LBR' & year==2006 & site %in% c('Firestone Hospital', 'Phebe Hospital', 'Liberia Govt Hosp−Buchanan', 'Martha Tubman Hospital'))]
anc <- anc[!(country=='LBR' & year==2007 & site %in% c('JJ Dossen Hospital', 'Senji Health Center', 'Karnplay Health Center', 'Voinjama Health Center', 'Fish Town Health Center'))]
anc <- anc[!(country=='GMB' & year %in% c(2006, 2007) & site %in% c('Farafenni', 'Sibanor', 'Kuntaur',  'Basse', 'Brikama', 'Sibanor', 'Essau'))]
anc <- anc[!(country=='CIV' & year==2008)]
anc <- anc[!(country=='BFA' & year==2008 & site == 'Sindou')]
anc <- anc[!(country=='BFA' & year==2009 & site == 'Ouahigouya')]
anc <- anc[!(country=='CMR' & year %in% c(2008, 2009) & site == 'Sangmelima')]
anc <- anc[!(country=='CMR' & year==2012 & site %in% c('HR Maroua', 'HR Ebolowa', 'Doual (HD Nylon)', 'Douala (HD Congo II)', 'Pouma', 'Saa', 'Ngaoundal', 'Molokol (Mokong)'))]
anc <- anc[!(country=='GHA' & year==2008 & site == 'Nadowli')]
anc <- anc[!(country=='NGA' & year==2005 & site %in% c('Abuja (Garki/Wuse)', 'Mubi, Lagos (Badagry)', 'Lagos (Ikeja)'))]
anc <- anc[!(country=='NGA' & year==2008 & site %in% c('Taraba (Zing)', 'Sokoto (Tambawal)', 'Borno (Maiduguri)', 'Kaduna (Kafanchan)', 'Plateau (Jos)', 'Abuja Bwari'))]
anc <- anc[!(country=='MDG' & year==2018 & site == 'Sainte Marie')]

anc <- anc[!(country=='COD' & year==2007 & site %in% c('KINGASANI', 'BUKAVU', 'MBANDAKA', 'LUBUMBASHI', 'BINZA METEO'))]
anc <- anc[!(country=='COD' & year==2011 & site %in% c('LUKULA', 'LODJA'))]
anc <- anc[!(country=='GAB' & year==2002 & site %in% c('SMI CS Glass 015', 'SMI CS Louis 011'))]
anc <- anc[!(country=='GIN' & year==2004 & site %in% c('LAFOU', 'KOULE', 'LEY SARE'))]
anc <- anc[!(country=='GMB' & year==2005 & site %in% c('Basse', 'Essau'))]
anc <- anc[!(country=='KEN' & year==2001 & site %in% c('Mt. Elgon', 'Kitui', 'Fatima'))]
anc <- anc[!(country=='KEN' & year==2011 & site %in% c('Lodwar', 'Tabaka'))]
anc <- anc[!(country=='SEN' & year==2004 & site =='DAKAR')]
anc <- anc[!(country=='TCD' & year==2011 & site %in% c('Sarh', 'Pala', 'Roi Fayçal'))]
anc <- anc[!(country=='UGA' & year==2005 & site =='Kagadi Hosp')]
anc <- anc[!(country=='UGA' & year==2010 & site %in% c('Kagadi Hosp', 'Mbarara Hosp'))]
anc <- anc[!(country=='UGA' & year==2003 & site %in% c('Mbale Hosp', 'Kilembe Hosp'))]
anc <- anc[!(country=='UGA' & year==2006 & site %in% c('Tororo Hosp', 'Nyakibaale Hosp', 'Nsambya Hospital'))]
anc <- anc[!(country=='UGA' & year==2012 & site =='Jinja Hosp')]

print(paste0('after dropping outliers, nrow ANC==' , nrow(anc)))

#-----------------------------------------------------------
#Assign more relevant values
agg_level ='age-sex-disaggregated'

    data<-copy(data_as_disagg)
  
  # replace year (the year the survey started in) with int_year (the median interview year for each location)
  data[, year := NULL]
  setnames(data, "int_year", "year")
  
  # note the data type
  data[, type := "Survey microdata"]
  
  #set agg data columns as NA for microdata
  data[, agebin_ag := NA]
  
  #-----------------------------------------------------------------
  #Add age-aggregated data
    # approximate the effective sample size using the design effects from the polygon microdata
    # Note: this is *very* rough; the purpose is to make sure the report/lit data doesn't have an unfair
    # advantage compared to microdata; the calculation below excludes surveys where there is no apparent
    # design effect -- this is not super plausible, and we don't want these to skew the distribution.
    de <- data[point == 0, list(de = sum(N)/sum(N_obs)), by = 'nid,shapefile,location_code'][de < 0.99999, median(de)]
    ag_data[, N_obs := N]
    ag_data[, N := N_obs * de]
    ag_data[, sum_of_sample_weights := N] # this is a very poor approximation...
    
    
    #Assign appropriate agebin values or age aggregated values
     ag_data <- ag_data[sex_id < 3]
     
     #Set up flexible agged bins for agebin_ag
     new_s_ages=data.table(c(15,20,25,30,35,40,45,50,55), 1:9)
     
     for(i in 1:nrow(new_s_ages)){
     ag_data[start_age == new_s_ages[i, V1],sag := new_s_ages[i, V2]]
     }
     
     new_e_ages=data.table(c(19,24,29,34,39,44,49,54,59), 1:9)
     
     for(i in 1:nrow(new_e_ages)){
       ag_data[end_age == new_e_ages[i, V1],eag := new_e_ages[i, V2]]
     }
     
    ag_data[sag==eag, agebin:=sag]
    ag_data[sag!=eag, agebin_ag := paste0(sag, ':', eag)]
    ag_data[sag!=eag, agebin :=NA]
    
    setnames(ag_data, 'survey_series', 'source', skip_absent=TRUE)
    
    data <- rbind(data, ag_data[, names(data), with = F])
    
    #Exclude situations where NAIIS is being included as both report and microdata
    data<- data[!(nid==399235 & is.na(agebin)), ]
    
    
    anc[, agebin_ag := '1:8']
    data <- rbind(data, anc[, c(names(data), "site"), with = F], fill = T)

  # repeat: drop observations with missing hiv test result or hiv weight (due to missingness of age group)
  data <- data[!is.na(hiv_test),]
  
  # now that all data is compiled, assign cluster_id (uniquely identifies a location-source-age-sex group)
  data[, cluster_id := 1:.N]
  
  print(names(data))
  

  ## Convert to counts and resample the polygon data -------------------------------------------------
  data[, hiv_test := hiv_test * N]

  #Add columns needed for modeling
  data[, weight:=1]
  data[, pseudocluster:=FALSE]
  
  ##Checks-----------------------------------------------------------------------------------------------
  # Check for missing lat/long or weights
  if (sum(is.na(data$hiv_test)) > 0 | sum(is.na(data$nid)) > 0 | sum(is.na(data$cluster_id)) > 0 | sum(is.na(data$year)) > 0) stop("There is missingness in hiv_test, nid, cluster_id, or year.")
  
  ## Format and check for changes --------------------------------------------------------------------
  # subset and rename variables
  data <- data[, list(nid, country, source, site, year, latitude, longitude, N, N_obs, hiv_test, weight, type, point, shapefile, cluster_id, sum_of_sample_weights, agebin, sex_id, agebin_ag, location_code, pseudocluster)]
  setkey(data, nid, country, source, year, site, latitude, longitude, N, hiv_test, agebin, sex_id)

  
  #Change the 'hiv_test' variable name to match the indicator
  names(data)[names(data) == 'hiv_test'] <- indicator
  
  ## Save collapsed data -----------------------------------------------------------------------------
  
  # add dates of geomatched data (geomatched_version) & collapsed data (collapse_version)
  data$geomatch_date <- geomatch_date
  data$report_commit <- report_commit
  data$anc_date <- anc_date
  data$collapse_date <- collapse_date

  
  #Set naming values and limit again to correct data types
   age_range=''
   file_tag <- '_age-poly-ag_data_disaggregated'
   all_collapsed <- data[(pseudocluster == FALSE & point == FALSE) | point == TRUE,]
   all_collapsed <- all_collapsed[sex_id != 3, ] 

        
    write.csv(all_collapsed, file = paste0(<<<< FILEPATH REDACTED >>>>, indicator, file_tag, age_range, ".csv"), row.names = FALSE)
    saveRDS(all_collapsed, file = paste0(<<<< FILEPATH REDACTED >>>>, indicator, file_tag, age_range, "_", format(Sys.Date(), "%Y_%m_%d"), ".rds"))



# write dates to log file
log_file <- paste0("<<<< FILEPATH REDACTED >>>>", indicator, "_date_log.csv")
if (!file.exists(log_file)) cat("geomatch_date,report_commit,anc_date, collapse_date", file = log_file, append = T) # if the log file doesn't exist, create it with the column names
cat("\r\n", file = log_file, append = T) # make sure you are on a new line before appending dates
cat(paste(geomatch_date, report_commit, anc_date,
          collapse_date, sep = ","), file = log_file, append = TRUE)
