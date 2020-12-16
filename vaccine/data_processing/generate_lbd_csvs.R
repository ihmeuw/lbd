# HEADER -------------------------------------------------------------------------
# Project: Vaccines: Data Prep
# Purpose: Take poly-resampled data and create the necessary data sets 
#          for ordinal regression
# Details: 
# source("FILEPATH/generate_lbd_csvs.R")
# Outputs: pointpoly-positioned csvs (in out_dir) ready for mbg
#********************************************************************************* 

# LOAD PACKAGES ---------------------------------------------------------------
username <- Sys.info()[["user"]]

core_repo      <- sprintf('FILEPATH/%s/FILEPATH/', username); setwd(core_repo)
indic_repo     <- sprintf('/FILEPATH/%s/FILEPATH/', username)
commondir      <- sprintf('/FILEPATH')
vaccines_repo  <- paste0("/FILEPATH/", username, "/FILEPATH/")


source("/FILEPATH/load_packages.R")
invisible(load_packages(c("proto", "findpython", "getopt", "argparse", "data.table", "magrittr", "survey", "binom", "parallel", "plyr", "dplyr", "haven",
                          "rgeos", "raster", "rgdal", "dismo", "gbm", "foreign", "doParallel", "grid", "gridExtra", "gtools", "ggplot2",
                          "sf", "fasterize", "assertthat")))

source(paste0(indic_repo, 'FILEPATH/misc_vaccine_functions.R'))
source(paste0(indic_repo, "FILEPATH/data_prep_functions.R"))
source(paste0(core_repo, "FILEPATH/graph_data_coverage.R"))

package_list   <- c(t(read.csv(sprintf('%s/package_list.csv', commondir), header=FALSE)))
str_match <- stringr::str_match


# GLOBAL VARIABLES ------------------------------------------------------------

outdir             <- "FILEPATH"
outlier_path       <- "FILEPATH/vaccination.csv"
lbd_regions_path   <- "FILEPATH/lbd_regions.csv"
resampled_lit_path <- "FILEPATH/lit_extraction_output_lbd/"

### Load MICS decider results
mics_decision_path <- "FILEPATH/mics_comparison_decisions_by_me_nid.csv"
mics_decisions     <- fread(mics_decision_path)


### Filter Table
record_data_drop       <- FALSE
filter_table_path      <- "FILEPATH"
save_filter_table_path <- "FILEPATH/final"
filter_table_columns   <- c("survey_id", "total", "missing_indicator", "n_with_age_data", "n_without_age_data", "n_12_59_months", "n_outside_age_range", 
                            "n_with_mcv_data", "n_without_mcv_data", "survey_in_year_range", "n_resampled", "n_not_resampled", "n_after_2000", "n_before_2000", 
                            "outlier", "n_after_outlier", "gps_clusters", "polygon_clusters", "n_gps_located", "n_polygon_located", "record_filter_successful")

if(record_data_drop){
  filter_table_path <- "FILEPATH"
  filter_table <- data.table()
  message("Loading all filter tables: beginning of generate csvs")
  files <- list.files(filter_table_path)
  files <- files[grepl(".csv", files)]
  for(file in files){
    temp <- fread(paste0(filter_table_path, file), nrows = 1)
    temp <- temp[, ..filter_table_columns]
    filter_table <- rbind(filter_table, temp)
  }
  filter_table <- filter_table[survey_id != 0]
}


# HELPER FUNCTIONS ------------------------------------------------------------

# Given vaccine prefix, load and rbind all data from resampled folder
load_resampled_data <- function(vax_prefix){
  df_vax <- data.table()
  for(file in list.files(paste0(extraction_root, "/FILEPATH/", vax_prefix, "/"))){
    temp <- data.table(read.csv(paste0(extraction_root, "/FILEPATH/", vax_prefix, "/", file))) 
    if(nrow(temp) > 0){
      df_vax <- data.table(rbind(df_vax, temp, fill=T))
    }
  }  
  return(df_vax)
}


# Data table, Data table, String -> Data table
# Given the vaccine data and the decider results, remove the rows from the vaccine data which will be replaced by the report data
merge_with_decisions <- function(dt_vaccine_data, mics_decisions, decider_vax, is_ratio){
  
  # Subset decider data to rows relevant to current vaccine in loop
  vaccines_of_interest <- unique(mics_decisions$me_name)[grepl(paste(paste0("vacc_", decider_vax), collapse = "|"), unique(mics_decisions$me_name))]
  mics_decisions_vax_specific <- mics_decisions[me_name %in% vaccines_of_interest, ]
  
  # Multi-dose vaccines must have same swap decision for all doses. Set swap decision to decision from first dose. 
  # For ratios between vaccines, set the decision for each vaccine based on the first dose. If either vaccine is swapped, ratio should be swapped
  mics_decisions_vax_specific <- mics_decisions_vax_specific[order(nid, ihme_loc_id, me_name), ]
  if (!is_ratio) {
    # Decider results should be standardized across all doses of multi-dose vaccines (except mcv and rcv)
    if (!(decider_vax %in% c("mcv", "rcv"))) {
      multi_dose_decisions_standardized <-
        ifelse((mics_decisions_vax_specific[, .("N_decisions" = length(unique(decision))), by=.(ihme_loc_id, nid)][N_decisions > 1, .N] > 0), 
               FALSE, 
               TRUE)
      if (!multi_dose_decisions_standardized) {
        stop(paste0("Conflicting decider results for different doses of ", decider_vax, ". \nStop and fix (should be accounted for in decider code)\n "))
      } else {
        # Reformat decider data to merge with vaccine data
        setnames(mics_decisions_vax_specific, c("ihme_loc_id", "nid"), c("country", "svy_id"))
        mics_decisions_vax_specific[, me_name := NULL]
        mics_decisions_vax_specific <- unique(mics_decisions_vax_specific)
        
        # Merge vaccine data with decider data
        dt_vaccine_data <- merge(dt_vaccine_data, mics_decisions_vax_specific, by=c("country", "svy_id"), all.x=TRUE)
        return(dt_vaccine_data)
      }
    } else {
      # message("")
      # my.name <- readline(prompt="Enter name: ")
    }
  } else {
    stop("merge_with_decisions not configured to merge ratio decisions - ratio estimates are no longer produced per Emily H's request 10/01/20")
  }
  
}

# Add lit extractions except for nids that already exist
add_lit_extractions <- function(dt_vaccine_data, dt_lit_extraction, decider_vax){
  
  if ("ihme_loc_id" %in% names(dt_lit_extraction)) {
    setnames(dt_lit_extraction, old = "ihme_loc_id", new = "country")
  }
  
  dt_lit_extraction  <- subset(dt_lit_extraction, select = names(dt_lit_extraction) %in% names(dt_vaccine_data))
  redundant_surveys  <- unique(dt_lit_extraction$svy_id)[unique(dt_lit_extraction$svy_id) %in% unique(dt_vaccine_data$svy_id)]
  dt_lit_extraction  <- dt_lit_extraction[!svy_id %in% redundant_surveys, ]
  dt_vaccine_data    <- rbind(dt_vaccine_data, dt_lit_extraction, fill=T)
  return(dt_vaccine_data)
}

# Load resampled lit extraction based on vaccine and date provided. If no date is provided, use most recent
load_literature_extraction <- function(path, vaccine_to_run, date = NULL) {
  
  # If date is provided, look for folder with appropriate date
  resample_folder_dates <- list.files(path)
  if (!is.null(date)) {
    date <- gsub("-|/", "_", date)
    folder_found <- date %in% resample_folder_dates
    if (!folder_found) {
      message(paste0("Date provided to load_literature_extraction() not found in resampled literature dates.\n\nDate provided: ", date, 
              "\nFolder dates: ", paste(resample_folder_dates, collapse = " ")))
      stop()
    } else {
      date_path <- paste0(path, date, "/")
    }
    
  # Otherwise, use most recent date
  } else {
    folder_dates         <- as.Date(resample_folder_dates, format = "%Y_%m_%d")
    ordered_folder_dates <- resample_folder_dates[rev(order(folder_dates))]
    most_recent_folder   <- ordered_folder_dates[[1]]
    date_path            <- paste0(path, most_recent_folder, "/")
    message("load_literature_extraction() function not given date argument. Defaulting to most recent lit resampling: \n",
            date_path)
  }
  
  # Identify vaccine-specific resample file at file path
  file <- list.files(date_path)[grepl(paste0(vaccine_to_run, "_lit_resample"), list.files(date_path))]
  if (length(file) == 0) {
    stop("Vaccine not found in resampled file path\n", "Vaccine: ", vaccine_to_run, "\nPath: ", date_path)
  } else {
    lit_resample_path <- paste0(date_path, file)
    lit_resample <- fread(lit_resample_path)
    return(lit_resample)
  }
}


# BEGIN LOOP ------------------------------------------------------------------

for (vaccine_to_run in c("mcv", "bcg", "hepb", "pcv", "polio", "yfv", "rota")) {
  is_ratio <- grepl("ratio", vaccine_to_run)
  
  # Set up the vaccine depending on the "vaccine_to_run" parameter
  message(vaccine_to_run)
  
  vaccine_list 	    	  <- NULL
  graph_vaccines    	  <- NULL
  final_model_vaccines  <- NULL
  graph_vaccine_titles  <- NULL
  
  if ("hepb3_dpt3_ratio" %in% vaccine_to_run) {
    outlier_names = c("vacc_hepb3", "vacc_dpt3")
    decider_vax   = c("hepb", "dpt")
  }
  
  if ("hib3_dpt3_ratio" %in% vaccine_to_run) {
    outlier_names = c("vacc_hib3", "vacc_dpt3")
    decider_vax   = c("hib", "dpt")
  }
  
  if ("pcv3_dpt3_ratio" %in% vaccine_to_run) {
    outlier_names = c("vacc_pcv3", "vacc_dpt3")
    decider_vax   = c("pcv", "dpt")
  }
  
  if ("rotac_dpt3_ratio" %in% vaccine_to_run) {
    outlier_names = c("vacc_rotac", "vacc_dpt3")
    decider_vax   = c("rota", "dpt")
  }
  
  if ("mcv2_mcv1_ratio" %in% vaccine_to_run) {
    outlier_names = c("vacc_mcv1", "vacc_mcv2")
    decider_vax   = c("mcv")
  }
  
  if("dpt3_timeliness_ratio" %in% vaccine_to_run) {
    outlier_names = c("vacc_dpt3")
    decider_vax   = c("dpt")
  }
  
  if ("polio" %in% vaccine_to_run) {
    if(is.null(vaccine_list)) vaccine_list <- list()
    vaccine_list$polio <- add_vaccine(prefix  = "polio", 
                                      title   = "Polio",
                                      doses   = 3,
                                      age_min = 12,
                                      age_max = 59)$polio
    
    vaccines <- names(vaccine_list)
    if(is.null(graph_vaccines))       graph_vaccines        <- character()
    if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
    if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
    graph_vaccines <- c(graph_vaccines, "polio3_cov")
    final_model_vaccines <- c(final_model_vaccines, "polio3_cov", "polio1_cov", "polio1_3_abs_dropout", "polio1_3_rel_dropout")
    graph_vaccine_titles <- c(graph_vaccine_titles, "Polio3")
    outlier_names = c("vacc_polio1", "vacc_polio2", "vacc_polio3")
  }
  
  if ("bcg" %in% vaccine_to_run) {
    
    vaccine_list <- add_vaccine(prefix = "bcg", 
                                title = "BCG",
                                doses = 1,
                                age_min = 12,
                                age_max = 59)
    
    vaccines <- names(vaccine_list)
    
    graph_vaccines <- c("bcg_cov")
    final_model_vaccines <- c("bcg_cov")
    graph_vaccine_titles <- c("BCG")
    outlier_names = c("vacc_bcg")
    
  }
  
  if ("rcv" %in% vaccine_to_run) {
    
    vaccine_list <- add_vaccine(prefix = "rcv", 
                                title = "RCV",
                                doses = 1,
                                age_min = 12,
                                age_max = 59)
    
    vaccines <- names(vaccine_list)
    
    graph_vaccines <- c("rcv_cov")
    final_model_vaccines <- c("rcv_cov")
    graph_vaccine_titles <- c("RCV")
    outlier_names = c("vacc_rcv")
    
  }
  if ("yfv" %in% vaccine_to_run) {
    
    vaccine_list <- add_vaccine(prefix = "yfv", 
                                title = "YFV",
                                doses = 1,
                                age_min = 12,
                                age_max = 59)
    
    vaccines <- names(vaccine_list)
    
    graph_vaccines <- c("yfv_cov")
    final_model_vaccines <- c("yfv_cov")
    graph_vaccine_titles <- c("YFV")
    outlier_names = c("vacc_yfv")
    
  }
  if ("mcv" %in% vaccine_to_run) {
    
    vaccine_list <- add_vaccine(prefix = "mcv", 
                                title = "MCV",
                                doses = 2,
                                age_min = 12,
                                age_max = 59)
    
    vaccines <- names(vaccine_list)
    
    graph_vaccines <- c("mcv1_cov", "mcv2_cov")
    final_model_vaccines <- c("mcv12_cond", "mcv1_cov", "mcv2_cov")
    graph_vaccine_titles <- c("MCV1", "MCV2")
    outlier_names = c("vacc_mcv1")
    vaccine_list$mcv$all_doses <- c("1", "2")
    
  }  
  if ("dpt" %in% vaccine_to_run) {
    if(is.null(vaccine_list)) vaccine_list <- list()
    vaccine_list$dpt <- add_vaccine(prefix = "dpt", 
                                    title = "DPT",
                                    doses = 3,
                                    age_min = 12,
                                    age_max = 59)$dpt
    
    vaccine_list$dpt["cond_doses"] <- c(12)
    vaccine_list$dpt["cond_vaccines"] <- "dpt12_cond"
    
    vaccines <- names(vaccine_list)
    if(is.null(graph_vaccines))       graph_vaccines        <- character()
    if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
    if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
    graph_vaccines <- c(graph_vaccines, "dpt3_cov")
    final_model_vaccines <- c(final_model_vaccines, "dpt12_cond", "dpt3_cov", "dpt1_cov", "dpt1_3_abs_dropout", "dpt1_3_rel_dropout")
    graph_vaccine_titles <- c(graph_vaccine_titles, "DPT3")
    outlier_names = c("vacc_dpt1", "vacc_dpt2", "vacc_dpt3")
  }
  if ("pcv" %in% vaccine_to_run) {
    if(is.null(vaccine_list)) vaccine_list <- list()
    vaccine_list$pcv <- add_vaccine(prefix = "pcv", 
                                    title = "PCV",
                                    doses = 3,
                                    age_min = 12,
                                    age_max = 59)$pcv
    
    vaccines <- names(vaccine_list)
    if(is.null(graph_vaccines))       graph_vaccines        <- character()
    if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
    if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
    graph_vaccines <- c(graph_vaccines, "pcv3_cov")
    final_model_vaccines <- c(final_model_vaccines, "pcv3_cov", "pcv1_cov", "pcv1_3_abs_dropout", "pcv1_3_rel_dropout")
    graph_vaccine_titles <- c(graph_vaccine_titles, "PCV3")
  }
  if ("hib" %in% vaccine_to_run) {
    if(is.null(vaccine_list)) vaccine_list <- list()
    vaccine_list$hib <- add_vaccine(prefix = "hib", 
                                    title = "hib",
                                    doses = 3,
                                    age_min = 12,
                                    age_max = 59)$hib
    
    vaccines <- names(vaccine_list)
    if(is.null(graph_vaccines))       graph_vaccines        <- character()
    if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
    if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
    graph_vaccines <- c(graph_vaccines, "hib3_cov")
    final_model_vaccines <- c(final_model_vaccines, "hib3_cov", "hib1_cov", "hib1_3_abs_dropout", "hib1_3_rel_dropout")
    graph_vaccine_titles <- c(graph_vaccine_titles, "hib3")
    outlier_names = c("vacc_hib", "vacc_hib1", "vacc_hib2", "vacc_hib3")
    
  }
  if ("hepb" %in% vaccine_to_run) {
    if(is.null(vaccine_list)) vaccine_list <- list()
    vaccine_list$hepb <- add_vaccine(prefix = "hepb", 
                                     title = "hepb",
                                     doses = 3,
                                     age_min = 12,
                                     age_max = 59)$hepb
    
    vaccines <- names(vaccine_list)
    if(is.null(graph_vaccines))       graph_vaccines        <- character()
    if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
    if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
    graph_vaccines <- c(graph_vaccines, "hepb3_cov")
    final_model_vaccines <- c(final_model_vaccines, "hepb3_cov", "hepb1_cov", "hepb1_3_abs_dropout", "hepb1_3_rel_dropout")
    graph_vaccine_titles <- c(graph_vaccine_titles, "hepb3")
  }
  if ("rota" %in% vaccine_to_run) {
    if(is.null(vaccine_list)) vaccine_list <- list()
    vaccine_list$rota <- add_vaccine(prefix = "rota", 
                                     title = "ROTA",
                                     doses = 1,
                                     age_min = 12,
                                     age_max = 59)$rota
    vaccine_list$rota$all_titles <- c("rota0", "rotac")
    vaccine_list$rota$cond_vaccines <- c("rotac_cond", "rota0_cond")
    vaccine_list$rota$cond_doses <- c("c", "0")
    vaccine_list$rota$all_doses <- c("0", "c")
    
    vaccines <- names(vaccine_list)
    if(is.null(graph_vaccines))       graph_vaccines        <- character()
    if(is.null(final_model_vaccines)) final_model_vaccines  <- character()
    if(is.null(graph_vaccine_titles)) graph_vaccine_titles  <- character()
    graph_vaccines <- c(graph_vaccines, "rotac_cov")
    final_model_vaccines <- c(final_model_vaccines, "rotac_cov")
    graph_vaccine_titles <- c(graph_vaccine_titles, "rotac")
    outlier_names = c("vacc_rota", "vacc_rotac", "vacc_rota1", "vacc_rota2", "vacc_rota3")
  }
  
  
  vax = vaccine_to_run
  extraction_root = "FILEPATH"
  
  # Define objects
  if (!grepl("ratio", vaccine_to_run)) {
    vax_prefix  <- vaccine_list[[vax]][["prefix"]]
    vax_title   <- vaccine_list[[vax]][["title"]]
    vax_doses   <- vaccine_list[[vax]][["all_doses"]]
    cond_vax    <- vaccine_list[[vax]][["cond_vaccines"]]
    all_vax     <- vaccine_list[[vax]][["all_vaccines"]]
    decider_vax <- vaccine_to_run
  } else {
    vax_prefix <- vaccine_to_run
  }
  
  # LOAD DATA -----------------------------------------------------------------------
  
  # load resampled microdata
  df_vax <- load_resampled_data(vax_prefix)
  raw_data_copy <- copy(df_vax)
  message(vaccine_to_run, " microdata loaded")
  
  # Merge microdata with decider results
  df_vax <- merge_with_decisions(df_vax, mics_decisions, decider_vax, is_ratio)
  
  # Drop rows which will be swapped out for resampled literature extractions
  removal_decisions <- c("use_report_extractions", "drop")
  df_vax <- df_vax[!decision %in% removal_decisions, ]  
  df_vax[, decision := NULL]
  
  # Load literature extractions
  if (!is_ratio) {
    dt_lit_extraction <- load_literature_extraction(resampled_lit_path, vaccine_to_run)
    # Merge microdata with literature extractions
    df_vax <- add_lit_extractions(df_vax, dt_lit_extraction, vaccine_to_run)
  }
  
  # APPLY OUTLIERING ----------------------------------------------------------------
  
  parse_nid <- function(arow) {
    
    if (is.na(arow$year_id)) {
      # NAs: just return as NA    
      return(arow)
    } else if (grepl("-", arow$year_id)) {
      # dashes: range inclusive
      year_start <- as.numeric(stringr::str_match(arow$year_id, "(.*)-.*")[,2])
      year_end   <- as.numeric(stringr::str_match(arow$year_id, ".*-(.*)")[,2])
      return(data.table(nid     = arow$nid,
                        year_id = year_start:year_end))
    } else if (grepl("/", arow$year_id)) {
      # slashes: just the specified elements
      years <- as.numeric(unlist(strsplit(arow$year_id, "/")))
      return(data.table(nid = arow$nid,
                        year_id =years))
    } else {
      return(arow)
    }
  }
  
  outlier <- fread(outlier_path)
  # Don't include GBD subnational outliers in outliering process
  subnational_rows <- grepl(outlier$ihme_loc_id, pattern = "_")
  outlier <- outlier[!subnational_rows, ]
  #
  outlier <- outlier[lbd == 1 & (me_name == "" | me_name %in% outlier_names),] #keep only lbd outliers specific to the vax or all vaccines
  outlier$year_id <- as.character(outlier$year_id)
  outlier_yrs <- outlier[!is.na(as.numeric(year_id)),] #identify year ranges we need to split
  outlier_nid <- outlier[is.na(year_id) | year_id == "" | batch_outlier==1,]  #outlier the entire nid if batch = 1 or year_id is blank
  outlier_y2  <- outlier_yrs[year_id != "",] #ignore missing years
  outlier_yf  <- lapply(1:nrow(outlier_y2), function(i)parse_nid(outlier_y2[i,])) %>% rbindlist #parse years 
  outlier_yrs <- rbind(outlier_yrs, outlier_yf, fill=T) #bind split years back into year level data
  outlier_ids <- paste(outlier_yrs$nid, outlier_yrs$year_id, sep = "_")
  df_vax[ , outlier := 0]
  df_vax[paste(svy_id, year_id, sep="_") %in% outlier_ids , outlier := 1]
  df_vax[(as.character(svy_id)) %in% unique(as.character(outlier_nid$nid)) , outlier := 1]
  
  
  ####  

  
  # ADD MBG REGIONS -----------------------------------------------------------------
  
  lbd_regions <- fread(lbd_regions_path)
  lbd_regions <- lbd_regions[, .(country, region)]
  df_vax      <- merge(df_vax, lbd_regions, all.x = TRUE, by = "country")
  
  # RECORD DATA DROPS ---------------------------------------------------------------
  
  
  if (record_data_drop & vaccine_to_run == "mcv") {
    
    lbd_nids <- unique(df_vax$svy_id)
    
    # No data should be dropped here. Any data that is dropped here will be recorded as a resample drop
    df_vax <- df_vax[!is.na(longitude) & !is.na(latitude) & !is.na(point), ]
    # df_vax <- df_vax[is.na(mcv_dose_1), ]
    
    # Rename Lit Filter Tables
    lit_filter_tables <- list.files(filter_table_path)[grepl("lit_", list.files(filter_table_path))]
    for (table in lit_filter_tables) {
      file.rename(from = paste0(filter_table_path, table), 
                  to = paste0(filter_table_path, sub(pattern = "lit_", replacement = "", table)))
    }
    
    # Record data drops for all surveys
    aggregated_filter_table <- data.table()
    lbd_nids <- unique(df_vax$svy_id)
    
    for (nid in lbd_nids) {
      
      filter_table <- fread(paste0(filter_table_path, nid, ".csv"))
      
      filter_table[survey_id == nid, n_resampled       := df_vax[svy_id == nid, sum(N_obs * weight)]]
      filter_table[survey_id == nid, n_not_resampled   := n_with_mcv_data - n_resampled]
      filter_table[survey_id == nid, n_after_2000      := df_vax[svy_id == nid & year_id > 1998, sum(N_obs * weight)]]
      filter_table[survey_id == nid, n_before_2000     := df_vax[svy_id == nid & year_id < 1999, sum(N_obs * weight)]]
      filter_table[survey_id == nid, outlier           := df_vax[svy_id == nid & outlier == 1, .N] > 0]
      filter_table[survey_id == nid, n_after_outlier   := df_vax[svy_id == nid & outlier == 0, sum(weight * N_obs)]]
      filter_table[survey_id == nid, gps_clusters      := df_vax[svy_id == nid & outlier == 0 & point == 1, .(latitude, longitude), by=.(latitude, longitude)][, .N]]
      filter_table[survey_id == nid, polygon_clusters  := df_vax[svy_id == nid & outlier == 0 & point == 0, .(location_code, shapefile), by=.(location_code, shapefile)][, .N]]
      filter_table[survey_id == nid, n_gps_located     := df_vax[svy_id == nid & outlier == 0 & point == 1, sum(N_obs * weight)]]
      filter_table[survey_id == nid, n_polygon_located := df_vax[svy_id == nid & outlier == 0 & point == 0, sum(N_obs * weight)]]
      
      filter_table[survey_id == nid, record_filter_successful := TRUE]
      write.csv(filter_table, file = paste0(filter_table_path, nid, ".csv"))
      
      aggregated_filter_table <- rbind(aggregated_filter_table, 
                                       filter_table)
    }
    
    # Make sure surveys with missing indicators or outside year range have 
    # removed numbers to avoid skewing totals. Redundant Check
    aggregated_filter_table[missing_indicator == TRUE | survey_in_year_range == FALSE,
                            setdiff(names(filter_table), c("survey_id", "missing_indicator", "survey_in_year_range", "outlier", "record_filter_successful")) := NA]
    write.csv(aggregated_filter_table, file = paste0("FILEPATH/data_cleaning_aggregated", nid, ".csv"))
    
  } 
  
  # Drop Rows with Missing Data
  df_vax <- df_vax[!is.na(longitude) & !is.na(latitude) & !is.na(point), ]
  
  # Add a rowID common between all derivative data sets
  # (for use in creation of holdouts)
  df_vax$row_id <- seq.int(nrow(df_vax))
  
  # Define list of vaccines to save csvs for
  # Need _cov for last dose, _cond for intermediate doses, and _cov for
  # any doses that are going to be the ultimate models (for validation)
  
  if(!grepl("ratio", vaccine_to_run)){
    last_vax  <- all_vax[length(all_vax)]
    last_dose <- max(vax_doses)
    
    csv_vaccines <- unique(c(cond_vax, 
                             paste0(vax_prefix, last_dose, "_cov"), 
                             final_model_vaccines))
    
    if(vaccine_to_run == "rota") { csv_vaccines <- "rotac_cov" }
    if(vax_prefix=="rota") { df_vax[, rota_dose_0 := N - rota_dose_1 - rota_dose_2 - rota_dose_3]}
    if(vax_prefix=="polio"){ df_vax[, polio_dose_0 := N - polio_dose_1 - polio_dose_2 - polio_dose_3]}
    if(vax_prefix=="dpt"){ df_vax[, dpt_dose_0 := N - dpt_dose_1 - dpt_dose_2 - dpt_dose_3]}
    
    if(vax_prefix=="bcg"){ df_vax[, bcg_dose_0 := N - bcg_dose_1]}
    if(vax_prefix=="rcv"){ df_vax[, rcv_dose_0 := N - rcv_dose_1]}
    if(vax_prefix=="yfv"){ df_vax[, yfv_dose_0 := N - yfv_dose_1]}
    if(vax_prefix=="dpt"){ df_vax[, dpt_dose_0 := N - (dpt_dose_1 + dpt_dose_2 + dpt_dose_3)]}
    if(vax_prefix=="pcv"){ df_vax[, pcv_dose_0 := N - (pcv_dose_1 + pcv_dose_2 + pcv_dose_3)]}
    if(vax_prefix=="hepb"){ df_vax[, hepb_dose_0 := N - (hepb_dose_1 + hepb_dose_2 + hepb_dose_3)]}
    if(vax_prefix=="hib"){ df_vax[, hib_dose_0 := N - (hib_dose_1 + hib_dose_2 + hib_dose_3)]}
    if(vax_prefix=="mcv"){ 
      df_vax[, mcv_dose_0 := N - mcv_dose_1]
      df_vax[!is.na(mcv_dose_2), mcv_dose_0 := N - (mcv_dose_1 + mcv_dose_2)]
    }
  }
  
  
  # SAVE CSVS ----------------------------------------------------------------------
  
  names(df_vax)[names(df_vax) == "survey_name"] <- "source"
  names(df_vax)[names(df_vax) == "year_id"]     <- "year"
  
  if (!grepl("ratio", vaccine_to_run)){
    # Catch if single dose
    if (last_dose == 1) csv_vaccines <- paste0(vax_prefix, last_dose, "_cov")
    
    message(paste0("\nSaving .csv files for ", vax_title, "...\n\n"))
    
    for (ind in csv_vaccines) {
      
      # Create csv for the last dose -----------------------------------------
      
      df_temp <- copy(df_vax)
      
      # Extract dose number
      if(ind == "rotac_cov") {
        dose <- 1
        last_dose <- 1
      } else {
        dose <- str_match(ind, "^[a-z]*([0-9]*)_.*")[2] %>% as.numeric 
      }
      
      if (grepl("_cov", ind)) {
        
        if (vax_prefix == "mcv" & dose == 1) {
          df_temp[, outcome := mcv_dose_1]
          df_temp[!is.na(mcv_dose_2), outcome := rowSums(.SD),
                  .SDcols = paste0(vax, "_dose_", dose:last_dose)] 
        } else {
          df_temp[, outcome := rowSums(.SD),
                  .SDcols = paste0(vax, "_dose_", dose:last_dose)]
        }
        # Drop extra columns
        drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
        df_temp <- subset(df_temp, select = !names(df_temp) %in% drop_vars)
      } else if (grepl("_cond", ind)) {
        
        if (dose == 12){
          # In the event that no data exists for dpt2, dpt1 is actually dose-specific dpt1/2
          # rather than dose-specific dpt1. Check whether dpt2 exists and calculate accordingly
          if (ind == "dpt12_cond") dose <- 2
          if (ind == "mcv12_cond") {
            # MCV 1/2 conditional is the proportion of people in regions where mcv2 has been introduced who receive
            # MCV1 given that they didn't receive MCV2. Subset to survey with MCV2 data
            df_temp <- df_temp[!is.na(mcv_dose_2), ]
            dose <- 1
          }
          
          # Identify doses up to conditional dose and rename
          setnames(df_temp, paste0(vax_prefix, "_dose_", 0:dose), paste0("d", 0:dose))
          
          # Remove doses greater than conditional dose
          drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
          df_temp   <- subset(df_temp, select = !(names(df_temp) %in% drop_vars))
          
          # Overwrite N with sum of remaining variables
          df_temp[, N := rowSums(.SD), .SDcols = paste0("d", 0:dose)]
          
          # Add dose-specific 1 to dose-specific 2 for 1-2 conditional
          df_temp[, paste0("d", dose) := rowSums(.SD), .SDcols = paste0("d", 1:dose)]
          # df_temp[, d2 := d2 + d1]
          
        } else {
          # Identify doses up to conditional dose and rename
          setnames(df_temp, paste0(vax_prefix, "_dose_", 0:dose), paste0("d", 0:dose))
          
          # Drop other dose columns (ignoring: conditional)
          drop_vars <- names(df_temp)[grepl(paste0(vax_prefix, "_dose"), names(df_temp))]
          
          # Subset to just the columns of interest
          df_temp <- subset(df_temp, select = !(names(df_temp) %in% drop_vars))
          
          # Overwrite N with sum of remaining variables
          df_temp[, N := rowSums(.SD), .SDcols = paste0("d", 0:dose)]
        }
        
        # Drop extra columns
        df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:(dose-1))))
        
        # Rename to "outcome" for standardization
        setnames(df_temp, paste0("d", dose), "outcome")
        
        # Drop if N = 0 (none in that row have [dose] or fewer doses)
        df_temp <- subset(df_temp, N > 0)
        
      } else if (grepl("1_3_abs_dropout", ind)) {
        
        setnames(df_temp, paste0(vax_prefix, "_dose_", 0:last_dose), paste0("d", 0:last_dose))
        df_temp[, outcome := d1 + d2]  
        
        # Drop extra columns
        df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:last_dose)))
        
      } else if (grepl("1_3_rel_dropout", ind)) {
        
        # Get the denominator (P doses > 1)
        df_temp[, denom := rowSums(.SD),
                .SDcols = paste0(vax, "_dose_", 1:last_dose)]
        
        setnames(df_temp, paste0(vax_prefix, "_dose_", 0:last_dose), paste0("d", 0:last_dose))
        
        # Calculate relative dropout
        df_temp[, outcome := (d1 + d2) / denom]  
        
        # Drop extra columns
        df_temp <- subset(df_temp, select = !(names(df_temp) %in% paste0("d", 0:last_dose)))
        df_temp[, denom := NULL]
        
      }
      
      #Deal with 100% plus coverage
      df_over <- df_temp[outcome > N,]
      df_over[ , large_diff := as.numeric(outcome/N > 1.02)]
      df_over <- df_over[large_diff == 1, ]
      
      if(nrow(df_over) > 0){
        x <- unique(df_over$svy_id)
        n <- length(unique(df_over$svy_id))
        message("WARNING for ", ind, ": several nids have observations with > 100% coverage and a difference of more than 2%.")
        message("NIDS: ", paste(x, collapse = ", "))
      } 
      df_temp[(outcome - N) > .01, outcome := N]
      
      # Drop Negatives
      df_temp <- df_temp[outcome >= 0, ]
      df_temp <- df_temp[!is.na(outcome), ]
      
      # Rename outcome to indicator
      setnames(df_temp, "outcome", ind)
      
      
      # Write
      
      if(ind %in% c("mcv1_cov", "bcg_cov", "polio3_cov", "dpt1_cov", "dpt3_cov", "hib3_cov", "hepb_cov", "pcv_cov", "rotac_cov", "bcg1_cov", "hepb3_cov", "pcv3_cov")){
        write.csv(df_temp, paste0(outdir, ind, "_outliers_inc.csv"))
      }
      
      df_temp <- df_temp[outlier==0, ]
      if(df_temp[is.na(point), .N] > 0) message("data is being saved for which point == NA. Consider removing")
      write.csv(df_temp, paste0(outdir, ind, ".csv"))
    }
  } else {
    ### For Ratios, drop outliers and "infinite" ratios (!0/0 or 0/0)
    df_vax <- df_vax[outlier==0, ]
    df_vax <- df_vax[get(vaccine_to_run) != Inf & get(vaccine_to_run) != -Inf, ]
    write.csv(df_vax, paste0(outdir, vaccine_to_run, ".csv"))
  }
}
