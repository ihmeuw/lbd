####################################################################################################
## Description: Collapse geomatched point and polygon data, add on report polygon data, and resample
##              polygons.
####################################################################################################

##############################################
##                                          ##
## Setting up and loading packages          ##
##                                          ##
##############################################

message("Setting up environment...")
rm(list=ls())

#set arguments 
input_version        <- format(Sys.time(), "%Y_%m_%d")
input_version_custom <- format(Sys.time(), "%Y_%m_%d")
prev_version         <- format(Sys.time(), "%Y_%m_%d")
sing                 <- T 
cores                <- 5 
core_repo <- paste0("<<<< FILEPATH REDACTED >>>>")

# Specify an indicator
indicator <- "ebf"

if (!dir.exists(core_repo)) core_repo <- "<<<< FILEPATH REDACTED >>>>"
commondir       <- "<<<< FILEPATH REDACTED >>>>"
indic_repo      <- "<<<< FILEPATH REDACTED >>>>"
stg             <- "1" 
stage           <- ifelse("3" %in% stg, "stage3", ifelse((("2a" %in% stg) | ("2b" %in% stg)), "stage2", ""))

# Load libraries and  MBG project functions.
message("Loading packages...")
package_list <- readLines("<<<< FILEPATH REDACTED >>>>")

source("<<<< FILEPATH REDACTED >>>>")
mbg_setup(package_list = package_list, repos = c(indic_repo, core_repo))

##############################################
##                                          ##
## Loading data                             ##
##                                          ##
##############################################

# Load data function
load_data <- function(ind, in_date, in_date_custom = NULL){
  message(paste0("Loading input dataset for ",ind,"..."))
  data <- fread(
    paste0("<<<< FILEPATH REDACTED >>>>", ind, 
           "/", ind, "_micro_", ifelse(stage == "", "", paste0(stage, "_")),
           format.Date(in_date, "%Y-%m-%d"), ".csv")
  )
  if(ind == "ebf"){
    in_date_custom <- input_version_custom
    data_custom <- fread(paste0("<<<< FILEPATH REDACTED >>>>", ind, "_micro_", in_date_custom, ".csv"))
    data <- rbindlist(list(data, data_custom), use.names = T, fill=TRUE)
  }
}

# Vector of indicators
ind_list <- as.list(indicator)
names(ind_list) <- indicator

# Create names list of all datasets
data_list <- lapply(ind_list, function(i) load_data(ind = i, in_date = input_version))

## Collapse to points & polygons --------------------------------------------------------------------------------------------
collapse <- function(i) {
  message(paste0("Collapsing ", i))
  dt <- get(i)
  # Drop observations with missing value or weight #
  message(paste0("Dropping ", nrow(dt[is.na(value) | is.na(pweight),]), " of ",  nrow(dt)," rows with either indicator or weight values = NA"))
  dt <- dt[!is.na(value) & !is.na(pweight),]
  byvars <- paste0("nid, country, source, year, point, shapefile, location_code, latitude, longitude", ifelse(i == "bfnow", ", age_month", ""))
  dt <- dt[, list(int_year = floor(median(int_year, na.rm = T)), 
                  mean     = weighted.mean(value, pweight),
                  N        = sum(pweight)^2/sum(pweight^2), 
                  sum_of_sample_weights = sum(pweight),
                  N_obs    = .N), 
           by=byvars]
  # Replace year (the year the survey started in) with int_year (the median interview year for each location) #
  dt[, year := NULL]
  setnames(dt, "int_year", "year")
  dt[, type := "Survey microdata"]
  return(dt)
}
alldata <- lapply(ind_list, collapse)

## Load and subset survey report data ------------------------------------------------------------------------------------------
# Load report data
for (i in indicator){
  if (i == "ebf"){
    report <- fread(paste0("<<<< FILEPATH REDACTED >>>>", i, "/", i, "_report_extraction.csv"))
    report <- report[, list(nid, iso3, source, point, location_name, shapefile, location_code, 
                            latitude, longitude, year, mean, N, ebf_lb, ebf_ub)]
    setnames(report, c("iso3"), c("country"))
    
    # Set negative ebf_lb to zero
    report[ebf_lb < 0, ebf_lb := 0]
    
    # Subset to point == 0
    report <- report[point == 0]
    
    # Drop countries outside Africa
    locs <- fread("<<<< FILEPATH REDACTED >>>>")
    stg_list <- locs[Stage %in% as.character(stg), iso3]
    report[, country := substr(country, 1, 3)] 
    message(paste0("Dropping "), nrow(report[!(country %in% stg_list),]), " rows outside of ", stage)
    report <- report[country %in% stg_list,]
    
    # Drop NIDs we have microdata for
    report <- report[!nid %in% alldata$ebf$nid]
    
    # Change estimates from a percentage to a proportion
    report[, `:=` (mean = mean/100, ebf_lb = ebf_lb/100, ebf_ub = ebf_ub/100)]
    
    # Estimate N if only CIs are available (using a normal approximation, assuming 95% CIs)
    report[is.na(N) & !is.na(ebf_lb) & !is.na(ebf_ub),
           N := (mean * (1 - mean)) / (((ebf_ub - ebf_lb) / (2 * 1.96))^2)]
    ci0 <- report[is.na(N) & ebf_lb == 0 & ebf_ub == 0, .N]
    if (ci0 > 0) warning(paste("Report data dropped because no sample size was provided and lower and upper bound are both 0 for", ci0, "source-country-year-locations"))
    report <- report[!(is.na(N) & ebf_lb == 0 & ebf_ub == 0),]
    report[, c("ebf_lb", "ebf_ub") := NULL]
    report[, type := "Survey report"]
    
    ## Combine report with collapsed microdata ------------------------------------------------------------------------------------
    de <- alldata$ebf[point == 0, list(de = sum(N)/sum(N_obs)), by = "nid,shapefile,location_code"][de < 0.99999, median(de)]
    report[, N_obs := N]
    report[, N := N_obs * de]
    
    # Combine with collapsed microdata
    report <- report[, list(nid, country, source, type, point, shapefile, location_code, latitude, longitude, year, mean, N, N_obs)]
    report[, sum_of_sample_weights := N] 
    alldata$ebf <- rbind(alldata$ebf, report, fill = T)
    
  }
  else next
}

## Resample polygons ------------------------------------------------------------------------------------------------------------

countdata <- lapply(indicator, function(i) {dt <- alldata[[i]][, mean := mean*N]; return(dt)})
names(countdata) <- indicator
resamp_stg <- "africa-yem"
resamp_data <- lapply(indicator, function(i) {resample_polygons(data=countdata[[i]], indic = "mean", cores = cores, 
                                                                gaul_list = get_gaul_codes(resamp_stg), 
                                                                shapefile_version = shapefile_version)})

unique_nids_before_resample <- unique(rbindlist(countdata, fill = TRUE)[, nid])
unique_nids_after_resample <- unique(rbindlist(resamp_data, fill = TRUE)[, nid])

message(paste0("Number of NIDs before resample: ", length(unique_nids_before_resample)))
message(paste0("Number NIDs after resample:  ", length(unique_nids_after_resample)))

names(resamp_data) <- indicator

## Format data & save ------------------------------------------------------------------------------------------------------------
message("Final formatting & saving data...")
for(i in indicator) {
  message(i)
  save_data <- resamp_data[[i]]
  if("svy_id" %in% colnames(save_data)) setnames(save_data, old = c("svy_id"), new = c("nid"))

  save_data[, date_geo_match := input_version]
  save_data[, collapse_date := format.Date(Sys.Date(), "%Y-%m-%d")]
  save_data <- save_data[, list(nid, country, source, year, latitude, longitude, N, N_obs, mean, weight, 
                                                 sum_of_sample_weights, point, shapefile, location_code, date_geo_match, collapse_date)] 
  
  save_data <- save_data[, cluster_id := .GRP, keyby = "nid, country, source, year, latitude, longitude, N, mean"] 
  
  setnames(save_data, "mean", i)
  
  ## Save input data--------------------------------------------------------------------------------------------------
  
  all_collapsed <- save_data
  
  write.csv(all_collapsed, file = paste0("<<<< FILEPATH REDACTED >>>>", i, "/", i, "_", 
                                         ifelse(stage == "", "", paste0(stage, "_")),
                                         format.Date(input_version, "%Y-%m-%d"), ".csv"), row.names=F)
  
  saveRDS(all_collapsed, file = paste0("<<<< FILEPATH REDACTED >>>>", i, "/", i, "_", 
                                       ifelse(stage == "", "", paste0(stage, "_")), input_version, ".rds"))
}

