#----HEADER-------------------------------------------------------------------------------------------------------------
# Purpose: Process extracted unit record data from UbCov by NID
#          
#               Creates indicators of interest
#               Apply restrictions (date formatting, age restrictions, vaccine definitions)
#               Explore card versus recall
#               Address missingness
#   
# Run:     source("FILEPATH/process.r", echo=TRUE)
# Inputs:  data  ->  NID-specific dataset from FILEPATH to process
#***********************************************************************************************************************


#----MISCELLANEOUS FUNCTIONS--------------------------------------------------------------------------------------------

job_hold <- function(job_name, file_list=NULL, obj=NULL) {
  
  ## Give it a sec to launch
  Sys.sleep(5)
  
  ## Start timer
  start.time <- proc.time()
  
  ## Wait for job to finish
  flag <-  0
  while (flag == 0) {
    ## Check if job is done
    if (system(paste0("qstat -r | grep ", job_name, "|wc -l"), intern=T) == 0) {
      ## If so, set flag to 1
      flag <- 1
    } else {
      Sys.sleep(5)
    }
  }
  
  ## End Timer
  job.runtime <- proc.time() - start.time
  job.runtime <- job.runtime[3]
  
  ## Give it another sec
  Sys.sleep(10)
  
  
  ## Check for the file list
  if (!is.null(file_list)) {
    missing_list <- NULL
    for (file in file_list) {
      ## Ensure that all files are there
      if (!file.exists(file)) {
        missing_list <- rbind(missing_list, file)
        ## Check obj if hdf
      } else {
        if (grepl(".h5", file_list[1])) {
          if (!(obj %in% h5ls(file_list)$name)) {
            missing_list <- rbind(missing_list, file)
          }
        }
      }
    }
    
    ## If missing_list > 0, break
    if (length(missing_list) > 0) {
      stop(paste0("Job failed: ", job_name,
                  "\nTime elapsed: ", job.runtime,
                  "\nYou are missing the following files: ", toString(missing_list)))
    }
  }
  
  ## Complete
  print(paste0("Job ", job_name, " has completed. Time elapsed: ", job.runtime))
  
}

### vetting code
vet <- function(data, vaccine_stems_in_dataset, vaccines_all_in_dataset, age_drops) {
  
  ### card fraction
  # antigen-specific
  cf_dataset <- data.table(ihme_loc_id=unique(data$ihme_loc_id))
  for (vax in vaccine_stems_in_dataset) {
    if (all(c(paste0(vax, "_card"), paste0(vax, "_dose")) %in% names(data))) {
      assign(paste0(vax, "_cf"), nrow(data[get(paste0(vax, "_card"))==1 & !is.na(get(paste0(vax, "_dose")))]) / nrow(data[!is.na(get(paste0(vax, "_dose")))]))
      cf_dataset[, (paste0(vax, "_card_fraction")) := get(paste0(vax, "_cf"))]
    }
  }
  # overall
  card_stems_in_dataset <- paste0(vaccine_stems_in_dataset, "_card")[paste0(vaccine_stems_in_dataset, "_card") %in% names(data)]
  data[, cf_numerator   := rowSums(.SD == 1, na.rm=TRUE), .SDcols=card_stems_in_dataset]
  data[, cf_denominator := rowSums(.SD == "0" | .SD == "At least 1" | .SD == "1" | .SD == "2" | .SD == "3" | .SD == "4" | .SD == "5" |
                                     .SD == 0 | .SD == 1 | .SD == 2 | .SD == 3 | .SD == 4 | .SD == 5, na.rm=TRUE),
       .SDcols=paste0(vaccine_stems_in_dataset, "_dose")[paste0(vaccine_stems_in_dataset, "_dose") %in% names(data)]]
  cf_overall <- nrow(data[!is.na(cf_numerator) & cf_numerator > 0]) / nrow(data[!is.na(cf_denominator) & cf_denominator > 0])
  data[, c("cf_numerator", "cf_denominator") := NULL]
  
  ### proportion from card / recall
  # antigen-specific
  cr_dataset <- data.table(ihme_loc_id=unique(data$ihme_loc_id))
  for (vax in vaccine_stems_in_dataset) {
    if (all(c(paste0(vax, "_dose_from_card"), paste0(vax, "_dose_from_recall")) %in% names(data))) {
      assign(paste0(vax, "_cr"),
             nrow(data[get(paste0(vax, "_dose_from_card"))==get(paste0(vax, "_dose")) & !is.na(get(paste0(vax, "_dose")))]) /
               nrow(data[!is.na(get(paste0(vax, "_dose")))]))
      cr_dataset[, (paste0(vax, "_proportion_from_card")) := get(paste0(vax, "_cr"))]
    }
  }
  
  ### missingness
  # antigen-specific
  missingness_dataset <- data.table(ihme_loc_id = unique(data$ihme_loc_id))
  for (vax in vaccine_stems_in_dataset) {
    if (paste0(vax, "_dose") %in% names(data)) {
      assign(paste0("missingness_", vax), nrow(data[is.na(get(paste0(vax, "_dose")))]) / nrow(data))
      missingness_dataset[, (paste0("missingness_", vax)) := get(paste0("missingness_", vax))]
    }
  }
  # overall vax
  data[, missingness := rowMeans(.SD == "0" | .SD == "At least 1" | .SD == "1" | .SD == "2" | .SD == "3" | .SD == "4" | .SD == "5" |
                                   .SD == 0 | .SD == 1 | .SD == 2 | .SD == 3 | .SD == 4 | .SD == 5, na.rm=TRUE),
       .SDcols=paste0(vaccine_stems_in_dataset, "_dose")[paste0(vaccine_stems_in_dataset, "_dose") %in% names(data)]]
  missingness_overall <- nrow(data[is.na(missingness) | is.nan(missingness)]) / nrow(data)
  data[, missingness := NULL]
  # overall vax in kids with non-missing age info
  data[!is.na(age_year), missingness_noage := rowMeans(.SD == "0" | .SD == "At least 1" | .SD == "1" | .SD == "2" | .SD == "3" | .SD == "4" | .SD == "5" |
                                                         .SD == 0 | .SD == 1 | .SD == 2 | .SD == 3 | .SD == 4 | .SD == 5, na.rm=TRUE),
       .SDcols=paste0(vaccine_stems_in_dataset, "_dose")[paste0(vaccine_stems_in_dataset, "_dose") %in% names(data)]]
  missingness_overall_noage <- nrow(data[is.na(missingness_noage) | is.nan(missingness_noage)]) / nrow(data)
  data[, missingness_noage := NULL]
  # lat / long
  if ("latitude" %in% names(data)) missingness_lat <- nrow(data[is.na(latitude)]) / nrow(data) else missingness_lat <- 1
  if ("longitude" %in% names(data)) missingness_long <- nrow(data[is.na(longitude)]) / nrow(data) else missingness_long <- 1
  if ("geospatial_id" %in% names(data)) missingness_geo <- nrow(data[is.na(geospatial_id)]) / nrow(data)
  if ("geo_id" %in% names(data)) missingness_geo <- nrow(data[is.na(geo_id)]) / nrow(data)
  if ((!"geo_id" %in% names(data) & (!"geospatial_id" %in% names(data)))) missingness_geo <- 1
  # survey design
  if ("psu" %in% names(data)) missingness_psu <- nrow(data[is.na(psu)]) / nrow(data) else missingness_psu <- 1
  if ("strata" %in% names(data)) missingness_strata <- nrow(data[is.na(strata)]) / nrow(data) else missingness_strata <- 1
  if ("pweight" %in% names(data)) missingness_pweight <- nrow(data[is.na(pweight)]) / nrow(data) else missingness_pweight <- 1
  if ("pweight_admin_1" %in% names(data)) missingness_pweight_admin_1 <- nrow(data[is.na(pweight_admin_1)]) / nrow(data) else missingness_pweight_admin_1 <- NA_integer_
  if ("pweight_admin_2" %in% names(data)) missingness_pweight_admin_2 <- nrow(data[is.na(pweight_admin_2)]) / nrow(data) else missingness_pweight_admin_2 <- NA_integer_
  if ("pweight_admin_3" %in% names(data)) missingness_pweight_admin_3 <- nrow(data[is.na(pweight_admin_3)]) / nrow(data) else missingness_pweight_admin_3 <- NA_integer_
  
  ### save single-row data table of information
  vet_log   <- data.table("nid"                     = unique(data$nid),
                          "datestamp"               = format(lubridate::with_tz(Sys.time(), tzone="America/Los_Angeles"), "%Y-%m-%d %H:%M"),
                          "ihme_loc_id"             = unique(data$ihme_loc_id),
                          "survey_name"             = unique(data$survey_name),
                          "year_start"              = unique(data$year_start),
                          "year_end"                = unique(data$year_end),
                          #"years"                   = ifelse(!"year_id" %in% names(data), "None", paste(sort(unique(data$year_id)), collapse=", ")), # now we are making cohorts in prep_for_tabulation()
                          "sample_size"             = nrow(data),
                          "included_vaccines"       = paste(vaccines_all_in_dataset, collapse=", "),
                          "latitude"                = ifelse("latitude" %in% names(data), 1, 0),
                          "longitude"               = ifelse("longitude" %in% names(data), 1, 0),
                          "card_fraction"           = cf_overall,
                          "missingness_vax"         = missingness_overall,
                          "missingness_vax_noage"   = missingness_overall_noage,
                          "missingness_age"         = nrow(data[is.na(age_year)]) / nrow(data),
                          "exclusions_age"          = age_drops,
                          "missingness_lat"         = missingness_lat,
                          "missingness_long"        = missingness_long,
                          "missingness_strata"      = missingness_strata,
                          "missingness_psu"         = missingness_psu,
                          "missingness_pweight"     = missingness_pweight,
                          "missingness_pweight_ad1" = missingness_pweight_admin_1,
                          "missingness_pweight_ad2" = missingness_pweight_admin_2,
                          "missingness_pweight_ad3" = missingness_pweight_admin_3
  )
  vet_log   <- merge(vet_log, cf_dataset, by="ihme_loc_id")
  vet_log   <- merge(vet_log, cr_dataset, by="ihme_loc_id")
  vet_log   <- merge(vet_log, missingness_dataset, by="ihme_loc_id")
  
  ### save
  write.csv(vet_log, file.path(extraction_root, "log/details", paste0(nid, ".csv")), row.names=FALSE)
  
}

### save datestamp log
datestamp <- function(date) {
  
  ### save NID to date file to see progress
  date_extraction_root <- file.path(extraction_root, "log/datestamp", date)
  if (!dir.exists(date_extraction_root)) dir.create(date_extraction_root, recursive=TRUE)
  cat(paste0("\n", nid, "\n"), file=file.path(date_extraction_root, paste0(nid, ".txt")), append=TRUE)
  
}

### load dataset by NID
load_nid <- function(nid, folder) {
  
  if (folder=="raw") {
    all_datasets  <- list.files(file.path(extraction_root, folder))
    dataset       <- all_datasets[grep(paste0("_", nid, ".csv"), all_datasets)]
  } else {
    dataset       <- paste0(nid, ".csv")
  }
  if (length(dataset) < 1 & paste0(nid, ".txt") %in% list.files(file.path(extraction_root, "log/datestamp", date))) {
    cat("\nnidnotexist", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
  }
  if (length(dataset) < 1){
    stop(paste0("No files for NID '", nid, "' exist in '", file.path(extraction_root, folder), "'"))
  }
  # read and return
  data            <- lapply(file.path(extraction_root, folder, dataset), fread) %>% rbindlist(., fill=TRUE)
  message(paste0("Reading in ", paste(file.path(extraction_root, folder, dataset), collapse=", ")))
  return(data)
  
}

### make sure all needed columns exist, and at least one vaccine column
check_completeness <- function(data, cols, date, nid, team) {
  
  ### check
  if (all(cols %in% names(data))) {
    if (length(vaccines[vaccines %in% names(data)]) > 0) {
      check <- TRUE
    } else {
      cat("\nnovaxind", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
      message("Missing at least one vaccine indicator; not moving on to tabulation")
      check <- FALSE
    }
  } else {
    cat("\nmissingind", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
    if (team=="lbd" & !"geo_id" %in% names(data) & !"geospatial_id" %in% names(data)) cat("\nmissinggeoind", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
    if (!"year_id" %in% names(data)) cat("\nmissingyearid", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
    message(paste0("Missing ", paste(cols[!cols %in% names(data)], collapse=", "), " from dataset; not moving on to tabulation"))
    check <- FALSE
  }
  
  ### age indicators
  if ("age_year" %in% names(data)) { if (nrow(data[!is.na(age_year)]) == 0) {
    cat("\nnoageind", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
    message(paste0("Missing age data from all observations in dataset; not moving on to tabulation"))
    check <- FALSE
  } }
  
  ### return
  return(check)
  
}

### function to rename columns if they are a different name conditional on the wrong name existing
check_and_rename_col <- function(data, old_name, new_name, quiet=TRUE) {
  
  if (old_name %in% names(data)) {
    setnames(data, old_name, new_name)
    if (!quiet) message(paste0("Renamed ", old_name, " to ", new_name))
  }
  else {
    if (!quiet) message(paste0("No column named ", old_name, " to rename"))
  }
  
}

### make full coverage indicator
make_full_coverage <- function(data, targets) {
  
  ### create full coverage binary if dataset includes all indicators
  if (all(indicators %in% names(data))) {
    data[eval(parse(text=paste(paste0("!is.na(", indicators, ")"), collapse=" & "))),
         full := ifelse(eval(parse(text=paste(paste0(indicators, "==1"), collapse=" & "))), 1, 0)]
    data[eval(parse(text=paste(paste0("!is.na(", indicators, ")"), collapse=" & "))),
         paste(indicators, collapse="_") := ifelse(eval(parse(text=paste(paste0(indicators, "==1"), collapse=" & "))), 1, 0)]
    min_age <- max(targets[vaccine %in% indicators, age_cohort])
    data[age_year < min_age, full := NA_integer_]
  }
  
  ### create full coverage binary if dataset includes all indicators
  if (all(indicators_no_mcv2 %in% names(data))) {
    data[eval(parse(text=paste(paste0("!is.na(", indicators_no_mcv2, ")"), collapse=" & "))),
         paste(indicators_no_mcv2, collapse="_") := ifelse(eval(parse(text=paste(paste0(indicators_no_mcv2, "==1"), collapse=" & "))), 1, 0)]
    min_age <- max(targets[vaccine %in% indicators_no_mcv2, age_cohort])
    data[age_year < min_age, full := NA_integer_]
  }
  
  ### create binary if dataset includes all following vaccines
  list_indicators <- c("three_indicators", "four_indicators_1", "four_indicators_2", "four_indicators_3", "four_indicators_4")
  for (list_ in list_indicators) {
    if (all(get(list_) %in% names(data))) {
      data[eval(parse(text=paste(paste0("!is.na(", get(list_), ")"), collapse=" & "))),
           paste(get(list_), collapse="_") := ifelse(eval(parse(text=paste(paste0(get(list_), "==1"), collapse=" & "))), 1, 0)]
      min_age <- max(targets[vaccine %in% get(list_), age_cohort])
      data[age_year < min_age, paste(get(list_), collapse="_") := NA_integer_]
    }
  }
  
  ### compute conditional probabilities
  for (probability in probabilities_to_model) {
    # list indicators
    all    <- strsplit(gsub("_", " ", probability), split=" +")[[1]]
    root   <- gsub("\\s*|_.*", "", probability)
    others <- all[all != root]
    # compute conditional probability of the root antigen (conditional on the other anitgens)
    if (all(all %in% names(data))) {
      data[eval(parse(text=paste(paste0("!is.na(", all, ")"), collapse=" & "))) &
             eval(parse(text=paste(paste0(others, "==1"), collapse=" & "))),
           (paste0("correlation_", probability)) := ifelse(get(root)==1, 1, 0)]
      # if child is younger than the minimum target age, exclude
      min_age <- max(targets[vaccine %in% all, age_cohort])
      data[age_year < min_age, (paste0("correlation_", probability)) := NA_integer_]
    }
  }
  
  ### return
  return(data)
  
}

### make rotavirus c coverage from rota doses
make_rotac <- function(data) {
  
  ### read in dosage dataset
 doses <- readRDS("FILEPATH/vaccine_schedule.rds")[me_name=="vacc_rotac" & ihme_loc_id==unique(data$ihme_loc_id), .(doses)]
 intro <- readRDS("FILEPATH/vaccine_intro.rds")[me_name=="vacc_rotac" & ihme_loc_id==unique(data$ihme_loc_id), cv_intro] %>% unique
  
  ### cap doses at 3; if doses == 0, set to 2
  if (nrow(doses)==0) doses <- 0
  doses <- ifelse(doses >= 3, 3, 2)
  rota_dose <- paste0("rota", doses)
  
  ### set rota c dosess if rotavirus has been introduced
  if (length(intro)==0) stop(paste0("BREAK | Missing introduction year for ", unique(data$ihme_loc_id), "; need to prep introduction frame for this geography before continuing"))
  if (intro < unique(data$int_year) & rota_dose %in% names(data)) data$rotac <- data[, rota_dose, with=FALSE]
  
  ### return
  return(data)
  
}

### set special target populations
set_target_ages <- function(data, vaccines.=vaccines) {
  
  ### read in schedule dataset
  schedule <- readRDS("FILEPATH/vaccine_target.rds")
  schedule_antigens <- unique(schedule$me_name)
  # schedule_antigens <- "vacc_rcv1"
  schedule <- schedule[ihme_loc_id %in% unique(data$ihme_loc_id)]
  
  ### set antigen target ages
  vax_targets <- data.table()
  for (vax_target in paste0("vacc_", vaccines.)) {
    assign(paste0(vax_target, "_target"),
           ifelse(vax_target %in% schedule_antigens,
                  ifelse(vax_target %in% unique(schedule$me_name), unique(schedule[me_name==vax_target, age_cohort]), NA),
                  1))
    vax_targets <- rbind(vax_targets, data.table(me_name=vax_target, age_cohort=get(paste0(vax_target, "_target"))))
  }
  
  ### return
  #data <- merge(data, vax_targets, by="me_name", all.x=TRUE)
  #return(data)
  return(vax_targets)
  
}

### reshape long for tabulation
prep_for_tabulation <- function(nid, team="gbd", vaccines.=vaccines, filter_table=NA) {
  
  ### read in data
  data <- load_nid(nid, folder="processed")
  
  ### create full coverage
  vax_targets <- set_target_ages(data=data)[, vaccine := gsub("vacc_", "", me_name)]
  if (team=="gbd") data <- make_full_coverage(data, targets=vax_targets)
  
  ### create card/recall only indicators
  if (team=="gbd") {
    for (ea in c(vaccines, "rotac")) {
      if (ea %in% names(data)) data[exists(paste0(gsub("[0-9]", "", ea), "_card"))==1, (paste0(ea, "_CARD")) := get(ea)]
      if (ea %in% names(data)) data[exists(paste0(gsub("[0-9]", "", ea), "_card"))==0 & !is.na(get(ea)), (paste0(ea, "_RECALL")) := get(ea)]
    }
  }
  
  ### save datestamp if LBD only and datestamp doesn't already exist
  if (lbd_only) datestamp(date=date)
  
  ### check completeness (will stop function if fails)
  tab_cols <- c("nid", "survey_name", "survey_module", "file_path", "ihme_loc_id", "year_start", "year_end", "age_year")
  gbd_cols <- c("strata", "psu", paste0("pweight", c("_admin_1", "_admin_2", "_admin_3", "_admin_4", "_admin_5")), paste0("admin_", 1:5, "_id"))
  lbd_cols <- c("geospatial_id", "strata", "psu", "pweight", paste0("pweight", c("_admin_1", "_admin_2", "_admin_3", "_admin_4", "_admin_5")))
  if (team=="lbd") tab_cols <- c(tab_cols, "geospatial_id") %>% unique
  if (team=="gbd") tab_cols <- c(tab_cols, "pweight") %>% unique
  check <- check_completeness(data, cols=tab_cols, date=date, nid=nid, team=team)
  
  if (check) {
    
    ### add fix for U/R splits without pweights
    parent <- unique(data$ihme_loc_id)
    if (team=="gbd" & "admin_1_urban_id" %in% names(data) & "pweight_admin_1" %in% names(data)) {
      data_ur <- copy(data)[!is.na(admin_1_urban_id)][, ihme_loc_id := admin_1_urban_id]
      data_ur[, pweight := pweight_admin_1]
      data_ur[, c("admin_1_urban_id", "pweight_admin_1", "admin_1_id", "admin_2_id", "pweight_admin_2", "admin_3_id", "pweight_admin_3") := NULL]
      data_ur[ihme_loc_id==parent, ihme_loc_id := NA_character_]
      data <- rbind(data, data_ur[!is.na(ihme_loc_id)], fill=TRUE)
    }
    
    ### add fix for admin2
    if (team=="gbd" & "admin_2_id" %in% names(data) & "pweight_admin_2" %in% names(data)) {
      data_a2 <- copy(data)[!is.na(admin_2_id)][, ihme_loc_id := admin_2_id]
      data_a2[, pweight := pweight_admin_2]
      data_a2[, c("admin_1_urban_id", "pweight_admin_1", "admin_1_id", "admin_2_id", "pweight_admin_2", "admin_3_id", "pweight_admin_3") := NULL]
      data_a2[ihme_loc_id==parent, ihme_loc_id := NA_character_]
      data <- rbind(data, data_a2[!is.na(ihme_loc_id)], fill=TRUE)
    }
    
    ### add fix for admin2
    if (team=="gbd" & "admin_3_id" %in% names(data) & "pweight_admin_3" %in% names(data)) {
      data_a3 <- copy(data)[!is.na(admin_3_id)][, ihme_loc_id := admin_3_id]
      data_a3[, pweight := pweight_admin_3]
      data_a3[, c("admin_1_urban_id", "pweight_admin_1", "admin_1_id", "admin_2_id", "pweight_admin_2", "admin_3_id", "pweight_admin_3") := NULL]
      data_a3[ihme_loc_id==parent, ihme_loc_id := NA_character_]
      data <- rbind(data, data_a3[!is.na(ihme_loc_id)], fill=TRUE)
    }
    
    ### for LBD ordinal regression, drop children where we don't know exact number of doses
    if (team=="lbd") { 
      for (vax_stem_col in unique(gsub("[0-9]", "", vaccines))) {
        if (paste0(vax_stem_col, "_dose") %in% names(data) & !vax_stem_col %in% c("mcv", "mmr", "rcv")) {
          vax_dose_cols <- vaccines[grep(vax_stem_col, vaccines)]
          if (all(vax_dose_cols %in% names(data))) {
            data[eval(parse(text=paste(paste0("is.na(", vax_dose_cols, ")"), collapse=" | "))), (vax_dose_cols) := NA]
          } else {
            data[, vax_dose_cols[vax_dose_cols %in% names(data)] := NA]
          }
        }
      }
    }
    
    ### drop extra columns
    cols       <- c(tab_cols, get(paste0(team, "_cols")))[c(tab_cols, get(paste0(team, "_cols"))) %in% names(data)] %>% unique
    vax_cols   <- vaccines.[vaccines. %in% names(data)]
    ready_data <- data[, c(cols, vax_cols), with=FALSE]
    ready_data[, individual := row_number(nid)]
    ready_data <- melt.data.table(ready_data, id.vars=cols, measure.vars=vax_cols, variable.name="me_name")
    ready_data[, me_name := paste0("vacc_", me_name)]
    
    ### Record data to be dropped
    ### To Do: make compatible with vaccine_to_record variable
    if(record_data_drop & team == "lbd"){
      filter_table[survey_id == nid, n_with_mcv_data    := ready_data[me_name == "vacc_mcv1" & !is.na(value), .N]]
      filter_table[survey_id == nid, n_without_mcv_data := n_12_59_months - n_with_mcv_data]
    }
    
    ### Data dropped here
    ready_data <- ready_data[!is.na(value)]
    
    
    ### birth cohort assignment
    vax_targets <- set_target_ages(ready_data, vaccines.=vaccines.)
    ready <- merge(ready_data, vax_targets, by="me_name", all.x=TRUE)
   
    ready <- ready[age_year >= age_cohort, ]
    ready[, year_id := floor((year_start + year_end) / 2) - age_year]
    ready[, c("age_cohort") := NULL]
    ready <- ready[!is.na(year_id)]
    
    message("Birth cohort assignment complete")
    
    ### return
    if(record_data_drop){
      return(list(check, ready, filter_table))
    } else {
      return(list(check, ready))
    }
  } else {
    ### return
    return(list(check))
  }
}

### check for success
check_nids <- function(nid) {
  
  # check files to see if they exist (if so, deemed a success)
  if (paste0(nid, ".csv") %in% list.files(file.path(extraction_root, "processed")))     check_process <- "Success" else check_process <- "Failure"
  if (paste0(nid, ".csv") %in% list.files(file.path(extraction_root, "tabulated/gbd"))) check_gbd     <- "Success" else check_gbd     <- "Failure"
  if (paste0(nid, ".csv") %in% list.files(file.path(extraction_root, "tabulated/lbd"))) check_lbd     <- "Success" else check_lbd     <- "Failure"
  if (lbd_only) { check_process <- "Not launched (LBD only)"; check_gbd <- "Not launched (LBD only)" }
  if (gbd_only) { check_lbd <- "Not launched (GBD only)" }
  
  # read log messaging
  note <- readLines(file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")))
  
  # save reason
  reason <- ""
  if (check_process=="Failure") {
    if (any(grepl("abc", note)))           reason <- paste(reason, sep=ifelse(nchar(reason) >= 1, "; ", ""))
  }
  if (check_gbd=="Failure" | check_lbd=="Failure") {
    if (any(grepl("nidnotexist", note)))   reason <- paste(reason, "This NID doesn't exist yet in <FILEPATH>", sep=ifelse(nchar(reason) >= 1, "; ", ""))
    if (any(grepl("novaxind", note)))      reason <- paste(reason, "No vaccine indicators found in dataset", sep=ifelse(nchar(reason) >= 1, "; ", ""))
    if (any(grepl("noageind", note)))      reason <- paste(reason, "All observations in dataset missing age information", sep=ifelse(nchar(reason) >= 1, "; ", ""))
    if (any(grepl("missingind", note)))    reason <- paste(reason, "Dataset entirely missing at least one required survey design variable", sep=ifelse(nchar(reason) >= 1, "; ", ""))
    if (any(grepl("missinggeoind", note))) reason <- paste(reason, "Dataset entirely missing geospatial_id", sep=ifelse(nchar(reason) >= 1, "; ", ""))
    if (any(grepl("missingyearid", note))) reason <- paste(reason, "Dataset entirely missing year_id (likely because no children with non-missing age older than target age range)", sep=ifelse(nchar(reason) >= 1, "; ", ""))
    if (any(grepl("DataPre2000", note)))   reason <- paste(reason, "Dataset pre-2000 (skipping LBD steps)", sep=ifelse(nchar(reason) >= 1, "; ", ""))
  }
  
  # write log
  log <- fread(file.path(extraction_root, "log/errors/log.csv"))
  nid_log <- data.table(nid=nid, checked=format(lubridate::with_tz(Sys.time(), tzone="America/Los_Angeles"), "%Y-%m-%d %H:%M"),
                        check_process=check_process, check_gbd=check_gbd, check_lbd=check_lbd, notes=reason)
  write.csv(rbind(log, nid_log, fill=TRUE), file.path(extraction_root, "log/errors/log.csv"), row.names=FALSE)
  
  # message
  message(paste0("|| NID ", nid, " ||"))
  message(paste0("|| Processing:     ", check_process, " ||"))
  message(paste0("|| GBD tabulation: ", check_gbd, " ||"))
  message(paste0("|| LBD tabulation: ", check_lbd, " ||"))
  
}

### clear out old NID files in case re-running old NID
clean_up <- function(nid, lbd_only=FALSE, gbd_only=FALSE) {
  
  # check files to see if they exist (if so, clear out for new round of processing)
  if (!lbd_only) if (paste0(nid, ".csv") %in% list.files(file.path(extraction_root, "processed")))     system(paste0("rm ", file.path(extraction_root, "processed", paste0(nid, ".csv"))))
  if (!lbd_only) if (paste0(nid, ".csv") %in% list.files(file.path(extraction_root, "tabulated/gbd"))) system(paste0("rm ", file.path(extraction_root, "tabulated/gbd", paste0(nid, ".csv"))))
  if (!gbd_only) if (paste0(nid, ".csv") %in% list.files(file.path(extraction_root, "tabulated/lbd"))) system(paste0("rm ", file.path(extraction_root, "tabulated/lbd", paste0(nid, ".csv"))))
  if (!lbd_only) if (paste0(nid, ".csv") %in% list.files(file.path(extraction_root, "log/details")))   system(paste0("rm ", file.path(extraction_root, "log/details", paste0(nid, ".csv"))))
  
}

## Record data dropped at each stage of filtering. To be used in publications data flow flowcharts. 
create_filter_table <- function(){
  filter_table <- data.table(survey_id = 0, total = 0, missing_indicator = FALSE, #If missing indicator, should not have been extracted
                             n_with_age_data = 0, n_without_age_data  = 0,
                             n_12_59_months  = 0, n_outside_age_range = 0,
                             n_with_mcv_data = 0, n_without_mcv_data  = 0,
                             survey_in_year_range = TRUE,
                             n_resampled   = 0, n_not_resampled   = 0,
                             n_after_2000  = 0, n_before_2000     = 0,
                             outlier   = FALSE, n_after_outlier   = 0,
                             gps_clusters  = 0, polygon_clusters  = 0,
                             n_gps_located = 0, n_polygon_located = 0, 
                             record_filter_successful = FALSE)
  return(filter_table)
}

load_filter_table <- function(){
  all_files = list.files(filter_table_path)
  if(paste0(nid, ".csv") %in% all_files){
    filter_table <- fread(paste0(filter_table_path, nid, ".csv"))
  } else {
    message('Filter table not present in file path. Make sure the process script was run with "record_data_drops" as TRUE')
    record_data_drop <= FALSE
    return(NULL)
  }
  return(filter_table)
}

# Where birth dates are missing, determine birth dates from interview date and age_month. Create
# birth date indicators if not present in data.
prepare_birth_dates <- function(data) {

  # Determine what birth data is present
  contains_birth_month <- "birth_month" %in% names(data)
  contains_birth_year  <- "birth_year"  %in% names(data)
  
  # Add columns if missing
  if (!contains_birth_month) {
    data <- cbind(data, data.table("birth_month" = NA_integer_))
  }
  if (!contains_birth_year) {
    data <- cbind(data, data.table("birth_year" = NA_integer_))
  }
  
  # Add missing birth data from age and interview date - requires age_month, int_month, and int_year
  date_columns <- c("age_month", "int_year", "int_month")
  if (all(date_columns %in% names(data))) {
    
    # Floor all fractional age months (precautionary)
    data[, age_month := floor(age_month)]
    
    # Set missing birth_years
    data[is.na(birth_year) & !is.na(age_month) & !is.na(int_year) & !is.na(int_month), 
         birth_year := ifelse(age_month < int_month, 
                              int_year, 
                              ifelse((age_month %% 12) == (int_month %% 12),
                                     int_year + (((int_month - age_month) %/% 12) - 1),
                                     int_year + ((int_month - age_month) %/% 12)))]
    
    # Set missing birth_months
    data[is.na(birth_month) & !is.na(age_month) & !is.na(int_year) & !is.na(int_month), 
         birth_month := ifelse((age_month %% 12) == (int_month %% 12), 
                               12, 
                               (int_month - age_month) %% 12)]
  }
  return(data)
}

# Compare birth date with DPT3 card date to get DPT3 timeliness ratio
get_dpt3_timeliness_ratio <- function(data) {
  
  # Set default to NA
  data[, dpt3_timeliness_ratio := NA]
  
  # If card year is equal to birth year, DPT3 timeliness is TRUE 
  data[birth_year == card_dpt3_year, dpt3_timeliness_ratio := TRUE]
  
  # If card year is more than a year greater than birth year, DPT3 timeliness is FALSE
  data[(birth_year + 1) < card_dpt3_year, dpt3_timeliness_ratio := FALSE]
  
  # If card year is 1 year greater than birth year, compare months. If difference is exactly  12 months, compare days. If no days, drop
  data[((birth_year + 1) == card_dpt3_year) & (((12 - birth_month) + card_dpt3_month) < 12), dpt3_timeliness_ratio := TRUE]  # Month difference < 12
  data[((birth_year + 1) == card_dpt3_year) & (((12 - birth_month) + card_dpt3_month) > 12), dpt3_timeliness_ratio := FALSE] # Month difference > 12
  if (all(c("card_dpt3_day", "birth_day") %in% names(data))) {                                                         # Month difference = 12
    data[((birth_year + 1) == card_dpt3_year) & ((12 - birth_month) + card_dpt3_month == 12) & (card_dpt3_day >= birth_day), dpt3_timeliness_ratio := FALSE]
    data[((birth_year + 1) == card_dpt3_year) & ((12 - birth_month) + card_dpt3_month == 12) & (card_dpt3_day < birth_day),  dpt3_timeliness_ratio := TRUE]
  }
  
  # Return data
  return(data)
}

#***********************************************************************************************************************


#----START FUNCTION-----------------------------------------------------------------------------------------------------

### list out vaccines (order of doses matter - ascending!)
vaccines <- c(paste0("dpt", 1:3),
              paste0("pent", 1:3),
              paste0("tetra", 1:3),
              paste0("polio", 1:3),
              paste0("hepb", 1:3),
              paste0("hib", 1:3),
              paste0("pcv", 1:3),
              paste0("rota", 1:3),
              paste0("mcv", 1:2),
              paste0("mmr", 1:2),
              paste0("rcv", 1:2),
              "bcg",
              "yfv")

### list indicators of interest for full coverage models
# full coverage components
indicators         <- c("dpt3", "mcv1", "polio3", "mcv2", "hib3", "hepb3", "pcv3", "rotac")
indicators_no_mcv2 <- c("dpt3", "mcv1", "polio3", "hib3", "hepb3", "pcv3", "rotac")
three_indicators   <- c("dpt3", "polio3", "mcv1")
four_indicators_1  <- c("hepb3", three_indicators)
four_indicators_2  <- c("hib3", three_indicators)
four_indicators_3  <- c("pcv3", three_indicators)
four_indicators_4  <- c("rotac", three_indicators)

# conditional probabilities
double <- c("mcv1_dpt3", "mcv2_mcv1", "hepb3_dpt3", "hib3_dpt3", "pcv3_dpt3", "rotac_dpt3")
triple <- c("polio3_mcv1_dpt3")
new    <- c("mcv2", "hepb3", "rotac", "pcv3", "hib3")
quad   <- c(paste0(new, "_polio3_mcv1_dpt3"))
probabilities_to_model <- c(double, triple, quad,
                            "rotac_pcv3_hib3_hepb3_polio3_mcv1_dpt3",
                            "rotac_hib3_hepb3_polio3_mcv1_dpt3",
                            "pcv3_rotac_hib3_hepb3_polio3_mcv1_dpt3",
                            "hib3_rotac_pcv3_hepb3_polio3_mcv1_dpt3",
                            "hib3_hepb3_polio3_mcv1_dpt3",
                            "pcv3_hib3_hepb3_polio3_mcv1_dpt3",
                            "hepb3_rotac_pcv3_hepb3_polio3_mcv1_dpt3",
                            "mcv2_rotac_pcv3_hib3_hepb3_polio3_mcv1_dpt3",
                            "rotac_mcv2_pcv3_hib3_hepb3_polio3_mcv1_dpt3",
                            "mcv2_pcv3_hib3_hepb3_polio3_mcv1_dpt3")


#----PROCESS------------------------------------------------------------------------------------------------------------

### open function
process <- function(nid) {
  
  ## save datestamp log
  datestamp(date=date)
  
  ### get NID dataset
  data <- load_nid(nid, folder="raw")
  
  ## Load filter table, assign nid and total number of starting children
  if(record_data_drop) {
    filter_table <- create_filter_table()
    filter_table[nrow(filter_table), c("survey_id", "total") := list(nid, data[, .N])]
  }
  #***********************************************************************************************************************
  
  
  #----SET VARIABLES------------------------------------------------------------------------------------------------------
  
  ### tell me you're starting!
  message("Starting vaccination dose-specific custom processing")
  
  ### extract two versions of each indicator
  # one based on card and recall
  # count as vaccinated if: recall confirmation, date is present on card, or marked on card but no date or recall
  # one based only on card
  # count as vaccinated if: date present, or marked on card but no date
  
  ### current year
  current_year <- as.integer(format(lubridate::with_tz(Sys.time(), tzone="America/Los_Angeles"), "%Y"))
  
  ### Account for Buddhist calendar
  buddhist_nids <- c(12732,
                     148649,
                     296646,
                     331377)
  
  ### Get birth dates from age and interview dates where birth dates are missing
  data <- prepare_birth_dates(data)
  
  ### list out vaccines (order of doses matter - ascending!)
  mult <- c("dpt", "pent", "tetra", "polio", "mcv", "mmr", "hepb", "hib", "pcv", "rota", "rcv")
  one  <- c("bcg", "yfv")
  
  ### load file for vaccine dosage information
  vacc_dose <- fread("FILEPATH/vacc_dose.csv")
  all_stems <- unique(vacc_dose$var)
  
  ### Identify vaccine stems and all vaccine/doses available in dataset.
  ### If no vaccine-related indicators are present in dataset, stop processing
  vaccine_stems_in_dataset <- c()
  for (stem in all_stems) {
    max <- vacc_dose[var == stem, refdose]
    if (max > 1) possible_doses <- 1:max else possible_doses <- ""
    possible_vaccine_columns <- c(paste0("recall_", stem, "_ever"),
                                  paste0("recall_", stem, "_times"),
                                  paste0("card_", stem, possible_doses, "_dmy"),
                                  paste0("card_", stem, possible_doses, "_day"),
                                  paste0("card_", stem, possible_doses, "_month"),
                                  paste0("card_", stem, possible_doses, "_year"),
                                  paste0("card_", stem, possible_doses, "_date"),
                                  paste0("response_", stem, possible_doses))
    if (length(possible_vaccine_columns[possible_vaccine_columns %in% names(data)]) > 0) {
      vaccine_stems_in_dataset <- c(vaccine_stems_in_dataset, stem)
    }
  }
  if (length(vaccine_stems_in_dataset) == 0) {
    cat("\nnovaxind", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
    stop("No vaccine indicators in dataset; not moving on to data processing")
  }
  vaccines_all_in_dataset <- vaccines[grep(paste(vaccine_stems_in_dataset, collapse="|"), vaccines)]
  
  # For MICS, drop antigens without both card and recall indicators
  if (unique(data$survey_name) %in% c("UNICEF_MICS")){
    
    stem_list <- fread("FILEPATH/vacc_stem_lookup.csv") 
    vacc_list <- c()
    for (VACC in vaccines_all_in_dataset) {
      
      stem <- subset(stem_list, vacc== VACC)
      stem <- stem$stem
      
      recall_columns <- c(paste0("recall_", stem, "_ever"),
                          paste0("recall_", stem, "_times"))
      
      card_columns <- c(paste0("card_", VACC, "_day"),
                        paste0("card_", VACC, "_date"))
      
      if (length(recall_columns[recall_columns %in% names(data)]) == 2 & length(card_columns[card_columns %in% names(data)]) > 0) {
        vacc_list <- c(vacc_list, VACC)
      }
      if(VACC %in% c("dpt1", "pcv1", "polio1", "mcv1", "bcg", "hepb1", "hib1", "rota1", "yfv", "pent1", "tetra1", "mmr1", "rcv1")){
        if (paste0("recall_", stem, "_ever") %in% names(data) & length(card_columns[card_columns %in% names(data)]) > 0) {
          vacc_list <- c(vacc_list, VACC)
        }  
      }
    }
    vaccines_all_in_dataset <- unique(vacc_list)
    
    # Subset vaccine stems to those represented in dataset
    vaccine_stems_in_dataset <- lapply(vaccine_stems_in_dataset, function(vs) {
      if(any(grepl(vs, vaccines_all_in_dataset)) == TRUE) {
        return(vs)
      } else {
        return(NULL)
      }
    })  %>% unlist
    
    if (length(vaccine_stems_in_dataset) == 0) {
      cat("\nnovaxind", file=file.path(extraction_root, "log/datestamp", date, paste0(nid, ".txt")), append=TRUE)
      stop("No vaccine indicators in dataset; not moving on to data processing")
    }
  }
  

  ### Drop Measles SIA doses
  measles_campaign_columns <- names(data)[grepl("recall_measles_campaign", names(data))]
  if (length(measles_campaign_columns) > 0) {
    measles_campaign_equals_1_expression <- paste(paste0(measles_campaign_columns, " == 1 "), collapse = "| ")
    if("card_mcv1_day" %in% names(data))  data[card_mcv1_day  == 66 & eval(parse(text = measles_campaign_equals_1_expression)), c("card_mcv1_day", "card_mcv1_month", "card_mcv1_year") := NA]
    if("card_mcv2_day" %in% names(data))  data[card_mcv2_day  == 66 & eval(parse(text = measles_campaign_equals_1_expression)), c("card_mcv2_day", "card_mcv2_month", "card_mcv2_year") := NA]
    
    if("card_mcv1_date" %in% names(data)) data[card_mcv1_date == 66 & eval(parse(text = measles_campaign_equals_1_expression)), card_mcv1_date := NA]
    if("card_mcv2_date" %in% names(data)) data[card_mcv2_date == 66 & eval(parse(text = measles_campaign_equals_1_expression)), card_mcv2_date := NA]
  } 
    
  ### Correct for alternative iso3 names
  data[ihme_loc_id == "ALG", ihme_loc_id := "DZA"]
  data[ihme_loc_id == "PAL", ihme_loc_id := "PSE"]
  
  ### assign true, false, missing, etc. variable labels
  code_variables <- c(paste0("has_vacc_card_", c("seen", "not_seen", "no")),
                      paste0("recall_source_", c("card", "recall")),
                      paste0("recall_ever_", c("true", "false")),
                      "recall_times_missing",
                      paste0("card_", c("nodate", "recall", "missing", "no", "date_format")),
                      paste0("response_", c("yes_card", "yes_recall", "no", "missing")))
  
  for (code_variable in code_variables) {
    if (code_variable %in% names(data)) {
      code_variable_new_name <- paste0(code_variable, "_code")
      assign(code_variable_new_name, data[, get(code_variable)] %>% unique %>% as.character)
      if (grepl(",", as.character(data[, get(code_variable)] %>% unique)) & !grepl(", ", as.character(data[, get(code_variable)] %>% unique))) {
        assign(code_variable_new_name, strsplit(get(code_variable_new_name), split=",")[[1]]) }
      if (grepl(", ", as.character(data[, get(code_variable)] %>% unique))) assign(code_variable_new_name, strsplit(get(code_variable_new_name), split=", ")[[1]])
      if (grepl(" ", as.character(data[, get(code_variable)] %>% unique)) & !grepl(", ", as.character(data[, get(code_variable)] %>% unique))) {
        assign(code_variable_new_name, gsub(" ", "", get(code_variable_new_name))) }
      if (length(grep(">=", as.character(get(code_variable_new_name)))) > 0) {
        for (each in grep(">=", as.character(get(code_variable_new_name)))) {
          current_code <- get(code_variable_new_name)[each]
          new_code     <- gsub(">=", "", current_code) %>% gsub(" ", "", .) %>% as.integer
          new_codes    <- (new_code):999999
          assign(code_variable_new_name, get(code_variable_new_name)[-each])
          assign(code_variable_new_name, c(get(code_variable_new_name), new_codes))
        }
      }
      if (length(grep("<=", as.character(get(code_variable_new_name)))) > 0) {
        for (each in grep("<=", as.character(get(code_variable_new_name)))) {
          current_code <- get(code_variable_new_name)[each]
          new_code     <- gsub("<=", "", current_code) %>% gsub(" ", "", .) %>% as.integer
          new_codes    <- -999:(new_code)
          assign(code_variable_new_name, get(code_variable_new_name)[-each])
          assign(code_variable_new_name, c(get(code_variable_new_name), new_codes))
        }
      }
      if (length(grep(">", as.character(get(code_variable_new_name)))) > 0) {
        for (each in grep(">", as.character(get(code_variable_new_name)))) {
          current_code <- get(code_variable_new_name)[each]
          new_code     <- gsub(">", "", current_code) %>% gsub(" ", "", .) %>% as.integer
          new_codes    <- (new_code + 1):999999
          assign(code_variable_new_name, get(code_variable_new_name)[-each])
          assign(code_variable_new_name, c(get(code_variable_new_name), new_codes))
        }
      }
      if (length(grep("<", as.character(get(code_variable_new_name)))) > 0) {
        for (each in grep("<", as.character(get(code_variable_new_name)))) {
          current_code <- get(code_variable_new_name)[each]
          new_code     <- gsub("<", "", current_code) %>% gsub(" ", "", .) %>% as.integer
          new_codes    <- -999:(new_code - 1)
          assign(code_variable_new_name, get(code_variable_new_name)[-each])
          assign(code_variable_new_name, c(get(code_variable_new_name), new_codes))
        }
      }
      if (grepl("missing", code_variable)) assign(code_variable_new_name, c(get(code_variable_new_name), "NA"))
      # make character string
      assign(code_variable_new_name, as.character(get(code_variable_new_name)))
    }
  }
  
  ### cleanup ages and sexes
  if ("child_age_month" %in% names(data)) {
    if ("age_month" %in% names(data)) data[is.na(age_month) & !is.na(child_age_month), age_month := child_age_month]
    data[!is.na(child_age_month), age_month := child_age_month]
  }
  if ("child_age_year" %in% names(data)) {
    if ("age_year" %in% names(data)) data[is.na(age_year) & !is.na(child_age_year), age_year := child_age_year]
    data[!is.na(child_age_year), age_year := child_age_year]
  }
  if ("child_sex_id" %in% names(data)) {
    if ("sex_id" %in% names(data)) data[is.na(sex_id) & !is.na(child_sex_id), sex_id := child_sex_id]
    data[!is.na(child_sex_id), sex_id := child_sex_id]
  }
  
  #***********************************************************************************************************************
  
  
  #----EXTRACT------------------------------------------------------------------------------------------------------------
  
  ###########################
  ### 1. PREP             ###
  ###########################
  
  ### vaccine card variable
  ### "decode" card information
  # 1 = seen
  # 2 = owns card but not seen
  # 3 = no card
  # first check to see if variable exists in dataset
  if ("has_vacc_card" %in% names(data)) {
    # convert to string
    data[, has_vacc_card := as.character(has_vacc_card)]
    for (condition in c("seen", "not_seen", "no")) {
      var <- paste0("has_vacc_card_", condition)
      if (var %in% names(data)) {
        # assign has_vacc_card variable to the value label tagged in the dataset
        data[has_vacc_card %in% get(paste0(var, "_code")), has_vacc_card := condition]
        data[, (var) := NULL]
      }
    }
    data[!has_vacc_card %in% c("seen", "not_seen", "no"), has_vacc_card := NA_character_]
  }
  
  ### vaccination-specific extraction loop
  for (vacc in vaccines_all_in_dataset) {
    if (!unique(data$survey_name) %in% c("UNICEF_MICS")){
      ###
      
      ### print out progress status
      message(paste0("|| ", vacc, " ||"))
      
      ### identify whether or not vaccine has multiple doses
      if (grepl(gsub("[0-9]", "", vacc), paste0(mult, collapse="|"))) {
        # save stem name
        vacc_stem <- mult[grep(gsub("[0-9]", "", vacc), mult)]
        # save dose number
        if (grepl("[1-9]", vacc)) {
          dose_num  <- substr(vacc, nchar(vacc), nchar(vacc)) %>% as.integer
          dose_stem <- dose_num
        }
        ### else if vaccine only has one dose
      } else if (grepl(vacc, paste0(one, collapse="|"))) {
        vacc_stem <- one[grep(vacc, one)]
        dose_num  <- 1
        dose_stem <- ""
      }
      
      ### generate variables of interest
      for (var in paste0(vacc_stem, c("_dose", "_card", "_ever", "_dose_from_card", "_dose_from_recall"))) {
        if(!var %in% names(data)) {
          data[, (var) := NA_integer_]
        }
      }
      
      ###########################
      ### 2. MATERNAL RECALL  ###
      ###########################
      
      ### remember ever having child receive this vaccine?
      # if recall_ever vars exist, decode
      if (paste0("recall_", vacc_stem, "_ever") %in% names(data) &
          "recall_ever_true" %in% names(data) &
          "recall_ever_false" %in% names(data)) {
        
        # convert to string
        data[, paste0("recall_", vacc_stem, "_ever") := as.character(get(paste0("recall_", vacc_stem, "_ever")))]
        
        # if recall ever receiving this vacc_stem, code vacc_stem_ever to 1 (else to 0)
        data[get(paste0("recall_", vacc_stem, "_ever")) %in% recall_ever_true_code,  paste0(vacc_stem, "_ever") := 1]
        data[get(paste0("recall_", vacc_stem, "_ever")) %in% recall_ever_false_code, paste0(vacc_stem, "_ever") := 0]
      }
      
      ### how many times have you received this vaccine?
      # single numeric variable with doses numbers
      recall_doses_variable <- paste0("recall_", vacc_stem, "_times")
      if (recall_doses_variable %in% names(data)) {
        
        # convert to string
        data[, (recall_doses_variable) := as.numeric(get(recall_doses_variable))]
        
        # if recall times is missing, replace vacc_stem_doses with NA from recall doses
        if ("recall_times_missing" %in% names(data)) {
          data[(!is.na(get(recall_doses_variable))) & (get(recall_doses_variable) %in% as.numeric(recall_times_missing_code)),
               (recall_doses_variable) := NA_integer_]
        }
        
        # if you remember receiving this vaccine a certain number of times, replace vacc_stem_doses (first defined from ever var) with recall doses
        data[( is.na(get(paste0(vacc_stem, "_dose"))) | get(recall_doses_variable) > get(paste0(vacc_stem, "_dose")) ) & !is.na(get(recall_doses_variable)),
             c(paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_recall")) := as.integer(get(recall_doses_variable))]
        
      }
      
      ###########################
      ### 3. CARD COVERAGE    ###
      ###########################
      
      ### splitting up of DMY and date variables, depending on what was entered into UbCov codebook,
      ### is still done within UbCov framework; see topic-specific custom code to adjust this if necessary
      
      ### clean up the card dates with missing values
      # convert date to day, month, and year
      card_vacc_root <- paste0("card_", vacc_stem, dose_stem)
      for (date_var in c("date")) {
        # set antigen-dose-specific variable
        variable <- paste0(card_vacc_root, "_", date_var)
        # make numeric
        if (variable %in% names(data)) {
          data[, (variable):= as.character(get(variable))]
          # convert to missing if value is in card_missing, card_nodate, or card_recall value labels
          codes <- c("NA", ".", "")
          for (card_val in c("missing", "nodate", "recall", "no")) {
            if (paste0("card_", card_val) %in% names(data)) codes <- c(codes, get(paste0("card_", card_val, "_code")))
          }
          data[get(variable) %in% unique(codes), (variable) := NA_character_]

          ### make dates into day/month/year variabes
          # set the format of the date
          if ("card_date_format" %in% names(data)) {
            if (card_date_format_code %in% c("dd/mm/yyyy", "DD/MM/YYYY")) date_format_string <- "%d/%m/%Y"
            if (card_date_format_code %in% c("ddmmyy")) date_format_string <- "%d%m%Y"
          } else {
            date_format_string <- ifelse(grepl("/", data[, variable, with=FALSE]), "%d/%m/%Y", "%d%B%Y")
          }
          
          # split out
          if (paste0(card_vacc_root, "_day") %in% names(data)) {
            data[is.na(get(paste0(card_vacc_root, "_day"))) & !is.na(get(variable)), paste0(card_vacc_root, "_day") := format(as.Date(get(variable), date_format_string), "%d")]
          } else {
            data[!is.na(get(variable)), paste0(card_vacc_root, "_day") := format(as.Date(get(variable), date_format_string), "%d")]
          }
          if (paste0(card_vacc_root, "_month") %in% names(data)) {
            data[is.na(get(paste0(card_vacc_root, "_month"))) & !is.na(get(variable)), paste0(card_vacc_root, "_month") := format(as.Date(get(variable), date_format_string), "%m")]
          } else {
            data[!is.na(get(variable)), paste0(card_vacc_root, "_month") := format(as.Date(get(variable), date_format_string), "%m")]
          }
          if (paste0(card_vacc_root, "_year") %in% names(data)) {
            data[is.na(get(paste0(card_vacc_root, "_year"))) & !is.na(get(variable)), paste0(card_vacc_root, "_year") := format(as.Date(get(variable), date_format_string), "%Y")]
          } else {
            data[!is.na(get(variable)), paste0(card_vacc_root, "_year") := format(as.Date(get(variable), date_format_string), "%Y")]
          }
          # clear out date if subset to day, month, and year
          data[!is.na(get(paste0(card_vacc_root, "_day"))) & !is.na(get(paste0(card_vacc_root, "_month"))) & !is.na(get(paste0(card_vacc_root, "_year"))) &
                 !is.na(get(paste0(card_vacc_root, "_date"))), paste0(card_vacc_root, "_date") := NA_character_]
        }
      }
      # day, month, and year
      card_vacc_root <- paste0("card_", vacc_stem, dose_stem)
      for (date_var in c("day", "month", "year")) {
        # set antigen-dose-specific variable
        variable <- paste0(card_vacc_root, "_", date_var)
        # make numeric
        if (variable %in% names(data)) {
          data[, (variable) := as.character(get(variable))]
          # convert to dates to missing if date value is in card_missing, card_nodate, or card_recall value labels
          codes <- c("NA", ".", "")
          for (card_val in c("missing", "nodate", "recall", "no")) {
            if (paste0("card_", card_val) %in% names(data)) codes <- c(codes, get(paste0("card_", card_val, "_code")))
          }
          data[get(variable) %in% unique(codes), (variable) := NA_character_]

          # interpret dates as numbers
          data[, (variable) := as.integer(get(variable))]
          # Buddhist calendar
          if (unique(data$nid) %in% buddhist_nids & date_var=="year") data[, (variable) := get(variable) - 543]
          # replace crazy date values with missingness
          if (date_var=="day")   data[get(variable) > 31, (variable) := NA_integer_]
          if (date_var=="month") data[get(variable) > 12, (variable) := NA_integer_]
          if (date_var=="year")  data[get(variable) > unique(as.integer(data$year_end)), (variable) := NA_integer_]
        }
      }
      
      ### replace doses with the card dates if month and year of vaccination not missing (!!)
      # month and year vars
      if (paste0(card_vacc_root, "_month") %in% names(data) & paste0(card_vacc_root, "_year") %in% names(data)) {
        
        ### if month and year vars not missing, replace doses with doses with month/year
        data[!is.na(get(paste0(card_vacc_root, "_month"))) & get(paste0(card_vacc_root, "_month")) != 0 &
               !is.na(get(paste0(card_vacc_root, "_year")))  & get(paste0(card_vacc_root, "_year")) != 0,
             c(paste0(vacc_stem, "_card"), paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card")) := list(1, (dose_num), (dose_num))]
        if (dose_num == 1) {
          data[(get(paste0(card_vacc_root, "_month")) == 0 | get(paste0(card_vacc_root, "_year")) == 0) & get(paste0(vacc_stem, "_card")) == 1, paste0(vacc_stem, "_dose_from_card") := 0]
          data[(get(paste0(card_vacc_root, "_month")) == 0 | get(paste0(card_vacc_root, "_year")) == 0) & is.na(get(paste0(vacc_stem, "_dose"))), paste0(vacc_stem, "_dose") := 0]
        }
                
      }
      
      ### replace doses with the card dates if date of vaccination not missing
      # applies to cases we are unable to parse out into day, month, and year
      if (paste0(card_vacc_root, "_date") %in% names(data)) {
        
        if (all(c(paste0(card_vacc_root, "_month"), paste0(card_vacc_root, "_year")) %in% names(data))) {
          ### if month and year vars missing but date is not missing, replace doses with doses with date
          data[is.na(get(paste0(card_vacc_root, "_month"))) & is.na(get(paste0(card_vacc_root, "_year"))) & !is.na(get(paste0(card_vacc_root, "_date"))),
               paste0(vacc_stem, "_card") := 1]
          data[is.na(get(paste0(card_vacc_root, "_month"))) & is.na(get(paste0(card_vacc_root, "_year"))) & get(paste0(vacc_stem, "_card")) == 1 &
                 !is.na(get(paste0(card_vacc_root, "_date"))) & get(paste0(card_vacc_root, "_date")) != 0,
               c(paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card")) := list((dose_num), (dose_num))]
          if (dose_num == 1) {
            data[is.na(get(paste0(card_vacc_root, "_month"))) & is.na(get(paste0(card_vacc_root, "_year"))) & get(paste0(card_vacc_root, "_date")) == 0 & get(paste0(vacc_stem, "_card")) == 1,
                 paste0(vacc_stem, "_dose_from_card") := 0]
            data[is.na(get(paste0(vacc_stem, "_dose"))) & is.na(get(paste0(card_vacc_root, "_month"))) & is.na(get(paste0(card_vacc_root, "_year"))) & get(paste0(card_vacc_root, "_date")) == 0,
                 paste0(vacc_stem, "_dose") := 0]
          }
          
        } else {
          ### if date is not missing, replace doses with doses with date
          data[!is.na(get(paste0(card_vacc_root, "_date"))), paste0(vacc_stem, "_card") := 1]
          data[!is.na(get(paste0(card_vacc_root, "_date"))) & get(paste0(card_vacc_root, "_date")) != 0 & get(paste0(vacc_stem, "_card")) == 1,
               c(paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card")) := list((dose_num), (dose_num))]
          if (dose_num == 1) {
            data[get(paste0(card_vacc_root, "_date")) == 0 & get(paste0(vacc_stem, "_card")) == 1, paste0(vacc_stem, "_dose_from_card") := 0]
            data[get(paste0(card_vacc_root, "_date")) == 0 & is.na(get(paste0(vacc_stem, "_dose"))), paste0(vacc_stem, "_dose") := 0]
          }
          
        }
        
        ### if the child has no response, no date, but has card --> infer that not vaccinated
        if ("has_vacc_card" %in% names(data) & dose_num == 1) {
          data[is.na(get(paste0(vacc_stem, "_dose"))) &               # hasn't responded yet
                 has_vacc_card == "seen",                               # child has vaccine card that was seen
               c(paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card"), paste0(vacc_stem, "_card")) := list(0, 0, 1)]
        }
      }
      
      ###########################
      ### 4. AGE AT VAX       ###
      ###########################
      
      ### split out if single date variable
      if (all(c(paste0("card_", vacc, "_dmy"), "card_date_format") %in% names(data))) {
        ### ...
      }
      
      ### calculate age at vaccination
      age_vars <- c(paste0(card_vacc_root, "_month"), paste0(card_vacc_root, "_year"), "birth_month", "birth_year")
      if (all(age_vars %in% names(data))) {
        
        # make numeric
        for (age_var in age_vars) if (!class(data[, get(age_var)]) %in% c("integer", "numeric")) data[, (age_var) := get(age_var) %>% as.integer]
        
        # check that card dates arent entirely missing
        if (sum(!is.na(data[, get(paste0(card_vacc_root, "_month"))])) > 0 &
            sum(!is.na(data[, get(paste0(card_vacc_root, "_year")) ])) > 0) {
          
          # clean up dates
          # fixing two-digit years
          max_year <- max(data[!is.na(get(paste0(card_vacc_root, "_year"))), (paste0(card_vacc_root, "_year")), with=FALSE], na.rm=TRUE)
          if (max_year < 100) {
            # if year less than 2000, always add 1900 to two-digit years
            data[year_end <  2000, paste0(card_vacc_root, "_year") := 1900 + get(paste0(card_vacc_root, "_year"))]
            # if year_end greater than 2000, allow for years that cross over 2000 by allowing years up to the current year
            data[year_end >= 2000, paste0(card_vacc_root, "_year") := ifelse( (2000 + get(paste0(card_vacc_root, "_year"))) <= current_year,
                                                                              2000 + get(paste0(card_vacc_root, "_year")),
                                                                              1900 + get(paste0(card_vacc_root, "_year")))]
          }
          
          # put logical range on month and year
          month_range <- 1:12
          year_range  <- 1900:current_year
          data[!get(paste0(card_vacc_root, "_month")) %in% month_range, paste0(card_vacc_root, "_month") := NA_integer_]
          data[!get(paste0(card_vacc_root, "_year")) %in% year_range, paste0(card_vacc_root, "_year") := NA_integer_]
          
          # calculate age in months at vaccination
          data[, paste0(vacc_stem, "_", dose_stem, "_age_month") := (12 * (get(paste0(card_vacc_root, "_year")) - 1900) + get(paste0(card_vacc_root, "_month")))
               - (12 * (birth_year - 1900) + birth_month)]
        }
      }
      
      ###########################
      ### 5. CARD & RECALL    ###
      ###########################
      
      ### address variable with responses for both card and recall (i.e. DHS)
      response_vacc_root <- paste0("response_", vacc_stem, dose_stem)
      if (response_vacc_root %in% names(data)) {
        # convert to string
        data[, (response_vacc_root) := as.character(get(response_vacc_root))]
        # if have card or recall dose
        if ("response_yes_card" %in% names(data)) {
          # change dose number to match card/recall in response column
          data[get(response_vacc_root) %in% response_yes_card_code, paste0(vacc_stem, "_dose") := (dose_num)]
          # apply card response to card doses
          data[get(response_vacc_root) %in% response_yes_card_code, paste0(vacc_stem, "_card") := 1]
          data[get(response_vacc_root) %in% response_yes_card_code, paste0(vacc_stem, "_dose_from_card") := (dose_num)]
        }
        if ("response_yes_recall" %in% names(data)) {
          # change dose number to match card/recall in response column
          data[get(response_vacc_root) %in% response_yes_recall_code & ( is.na(get(paste0(vacc_stem, "_dose"))) | (dose_num) > get(paste0(vacc_stem, "_dose")) ), paste0(vacc_stem, "_dose") := (dose_num)]
          data[get(response_vacc_root) %in% response_yes_recall_code & ( is.na(get(paste0(vacc_stem, "_dose_from_recall"))) | (dose_num) > get(paste0(vacc_stem, "_dose_from_recall")) ), paste0(vacc_stem, "_dose_from_recall") := (dose_num)]
        }
        # apply no's if even first dose isn't received
        if ("response_no" %in% names(data) & dose_num == 1) {
          # if 0 dose, assume don't remember a dose!
          data[get(response_vacc_root) %in% response_no_code & is.na(get(paste0(vacc_stem, "_dose"))), paste0(vacc_stem, "_dose") := 0]
          data[get(response_vacc_root) %in% response_no_code & is.na(get(paste0(vacc_stem, "_dose_from_recall"))), paste0(vacc_stem, "_dose_from_recall") := 0]
          # if has vaccine card but reports no, assume 0 doses is from vaccine card
          if ("has_vacc_card" %in% names(data)) {
            data[get(response_vacc_root) %in% response_no_code & has_vacc_card == "seen", paste0(vacc_stem, "_card") := 1]
            data[get(response_vacc_root) %in% response_no_code & has_vacc_card == "seen" & is.na(get(paste0(vacc_stem, "_dose_from_card"))), paste0(vacc_stem, "_dose_from_card") := 0]
          }
        }
      }
      
      ###########################
      ### 6. CLEAN            ###
      ###########################
      
      ### replace (vacc)_card = 0 if not filled out
      data[(is.na(get(paste0(vacc_stem, "_card"))) | get(paste0(vacc_stem, "_card")) != 1) & !is.na(get(paste0(vacc_stem, "_dose"))),
           paste0(vacc_stem, "_card") := 0]
      
      ### recall gateway (never had any vaccinations)
      if ("ever_vaccinated" %in% names(data)) {
        data[, ever_vaccinated := as.integer(ever_vaccinated)]
        data[is.na(get(paste0(vacc_stem, "_dose"))) & ever_vaccinated == 0, paste0(vacc_stem, "_dose") := 0]
        data[is.na(get(paste0(vacc_stem, "_dose_from_recall"))) & ever_vaccinated == 0, paste0(vacc_stem, "_dose_from_recall") := 0]
      }
      
      ### close vaccine-specific extraction loop
    }
    
    
    
    ########################################################################################################
    ########################################################################################################
    ###########################################   MICS!!!    ###############################################
    ########################################################################################################
    ########################################################################################################

    if(unique(data$survey_name) %in% c("UNICEF_MICS")){
      ### print out progress status
      message(paste0("|| ", vacc, " ||"))
      
      ### identify whether or not vaccine has multiple doses
      if (grepl(gsub("[0-9]", "", vacc), paste0(mult, collapse="|"))) {
        # save stem name
        vacc_stem <- mult[grep(gsub("[0-9]", "", vacc), mult)]
        # save dose number
        if (grepl("[1-9]", vacc)) {
          dose_num  <- substr(vacc, nchar(vacc), nchar(vacc)) %>% as.integer
          dose_stem <- dose_num
        }
        ### else if vaccine only has one dose
      } else if (grepl(vacc, paste0(one, collapse="|"))) {
        vacc_stem <- one[grep(vacc, one)]
        dose_num  <- 1
        dose_stem <- ""
      }
      
      ### generate variables of interest
      
      # if(dose_num==1){
      for (var in paste0(vacc_stem, c("_dose", "_card", "_ever", "_dose_from_card", "_dose_from_recall"))) {
        if(!var %in% names(data)) {
          data[, (var) := NA_integer_]
        }
      }
      # }
      
      ###########################
      ### 2. MATERNAL RECALL  ###
      ###########################
      
      card_vacc_root <- paste0("card_", vacc_stem, dose_stem)
      
      if(dose_num==1){
        # Set up initial values
        data[, paste0(vacc_stem, "_dose_from_recall") := NA_integer_]
        
        # Set zeros if recall ever is in the list of false codes
        data[as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code,
             paste0(vacc_stem, "_dose_from_recall") := 0]
      }
      
      if (paste0(card_vacc_root, "_day") %in% names(data)) {
        if (dose_num == 1) {
          ### If card day variable is 66 or recall_ever is 1, then we can say that there was a first dose recalled
          data[(get(paste0(card_vacc_root, "_day")) == 66 | get(paste0("recall_", vacc_stem, "_ever")) ==1),	
               c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]	
          # Flag those eligible for first dose recall	
          data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 	
                  get(paste0(card_vacc_root, "_day")) == 66 | get(paste0("recall_", vacc_stem, "_ever")) ==1),	
               paste0(vacc,"_recall_eligible") := TRUE]	
        } else if (dose_num == 2) {	
          ### If card day variable is 66 or recall_ever is 1 AND recall is exactly 2, then we can say that the second dose was recalled	
          data[get(paste0(card_vacc_root, "_day")) == 66 | (get(paste0("recall_",vacc_stem, "_ever")) == 1 & get(paste0("recall_",vacc_stem, "_times"))) == dose_num,	
               c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]	
          # Flag those eligible for second dose recall	
          data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 	
                  get(paste0(card_vacc_root, "_day")) == 66 | get(paste0("recall_", vacc_stem, "_ever")) ==1) & 	
                 !(as.character(get(paste0("recall_", vacc_stem, "_times"))) %in% recall_times_missing_code), 	
               paste0(vacc, "_recall_eligible") := TRUE]	
        } else if (dose_num > 2) {	
          ### If card day variable is 66 or recall_ever is 1 AND recall is between dose_num and 7, give them the dose_num as recall	
          data[get(paste0(card_vacc_root, "_day")) == 66 | (get(paste0("recall_",vacc_stem, "_ever")) == 1 & get(paste0("recall_",vacc_stem, "_times")) >= dose_num & get(paste0("recall_",vacc_stem, "_times")) <= 7),	
               c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]	
          # Flag those eligible for next dose recall	
          data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 	
                  get(paste0(card_vacc_root, "_day")) == 66 | get(paste0("recall_", vacc_stem, "_ever")) ==1) & 	
                 !(as.character(get(paste0("recall_", vacc_stem, "_times"))) %in% recall_times_missing_code), 	
               paste0(vacc, "_recall_eligible") := TRUE]
        }
      }
      
      ### replace doses with the card dates if date of vaccination not missing
      # applies to cases we are unable to parse out into day, month, and year
      if (paste0(card_vacc_root, "_date") %in% names(data)) {
        
        if (all(c(paste0(card_vacc_root, "_day")) %in% names(data))) {
          
          if (dose_num == 1) {
            data[is.na(get(paste0(card_vacc_root, "_day"))) & (get(paste0("recall_",vacc_stem, "_ever")) == 1),
                 c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]
            # Flag those eligible for first dose recall
            data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 
                    get(paste0("recall_", vacc_stem, "_ever")) ==1),
                 paste0(vacc,"_recall_eligible") := TRUE]
          } else if (dose_num ==2 ) {
            data[is.na(get(paste0(card_vacc_root, "_day"))) & (get(paste0("recall_",vacc_stem, "_ever")) == 1 & paste0("recall_",vacc_stem, "_times") == dose_num),
                 c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]
            # Flag those eligible for second dose recall
            data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 
                    get(paste0("recall_", vacc_stem, "_ever")) ==1) & 
                   !(as.character(get(paste0("recall_", vacc_stem, "_times"))) %in% recall_times_missing_code), 
                 paste0(vacc, "_recall_eligible") := TRUE]
            
          } else if (dose_num > 2) {
            data[is.na(get(paste0(card_vacc_root, "_day"))) & (get(paste0("recall_",vacc_stem, "_ever")) == 1 & get(paste0("recall_",vacc_stem, "_times")) >= dose_num & get(paste0("recall_",vacc_stem, "_times")) <= 7),
                 c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]
            # Flag those eligible for next dose recall
            data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 
                    get(paste0("recall_", vacc_stem, "_ever")) ==1) & 
                   !(as.character(get(paste0("recall_", vacc_stem, "_times"))) %in% recall_times_missing_code), 
                 paste0(vacc, "_recall_eligible") := TRUE]
            
          } else {
            ### if date is not missing, replace doses with doses with date
            if (dose_num == 1) {
              data[(get(paste0("recall_",vacc_stem, "_ever")) == 1),
                   c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]
              # Flag those eligible for first dose recall
              data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 
                      get(paste0("recall_", vacc_stem, "_ever")) ==1),
                   paste0(vacc,"_recall_eligible") := TRUE]
            } else if (dose_num ==2) {
              data[(get(paste0("recall_",vacc_stem, "_ever")) == 1 & paste0("recall_",vacc_stem, "_times") == dose_num),
                   c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]
              # Flag those eligible for second dose recall
              data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 
                      get(paste0("recall_", vacc_stem, "_ever")) ==1) & 
                     !(as.character(get(paste0("recall_", vacc_stem, "_times"))) %in% recall_times_missing_code), 
                   paste0(vacc, "_recall_eligible") := TRUE]
            } else if (dose_num > 2) {
              data[(get(paste0("recall_",vacc_stem, "_ever")) == 1 & get(paste0("recall_",vacc_stem, "_times")) >= dose_num & get(paste0("recall_",vacc_stem, "_times")) <= 7),
                   c(paste0(vacc_stem, "_dose_from_recall")) := list((dose_num))]
              # Flag those eligible for next dose recall
              data[(as.character(get(paste0("recall_", vacc_stem, "_ever"))) %in% recall_ever_false_code | 
                      get(paste0("recall_", vacc_stem, "_ever")) ==1) & 
                     !(as.character(get(paste0("recall_", vacc_stem, "_times"))) %in% recall_times_missing_code), 
                   paste0(vacc, "_recall_eligible") := TRUE]
            }
          }
        }
      }
      
      
      ###########################
      ### 3. CARD COVERAGE    ###
      ###########################
      
      
      ### splitting up of DMY and date variables, depending on what was entered into UbCov codebook,
      ### is still done within UbCov framework; see topic-specific custom code to adjust this if necessary
      
      ### clean up the card dates with missing values
      # convert date to day, month, and year
      card_vacc_root <- paste0("card_", vacc_stem, dose_stem)
      if(dose_num == 1){data[,paste0(vacc_stem, "_dose_from_card") := NA_integer_]}
      if(dose_num == 1){data[,paste0(vacc_stem, "_true_zero") := 0]}
      
      
      for (date_var in c("date")) {
        # set antigen-dose-specific variable
        variable <- paste0(card_vacc_root, "_", date_var)
        # make numeric
        if (variable %in% names(data)) {
          data[, (variable):= as.character(get(variable))]
          # convert to missing if value is in card_missing, card_nodate, or card_recall value labels
          codes <- c("NA", ".", "")
          for (card_val in c("missing")) {
            # for (card_val in c("missing", "nodate", "recall", "no")) {
            if (paste0("card_", card_val) %in% names(data)) codes <- c(codes, get(paste0("card_", card_val, "_code")))
          }
          data[get(variable) %in% unique(codes), (variable) := NA_character_]
          # # convert to no (0) if value is in card_no value labels
          for (card_val in c("no")) {
            if (paste0("card_", card_val) %in% names(data)) {
              codes <- get(paste0("card_", card_val, "_code"))
              data[get(variable) %in% c(codes), (variable) := "0"]
            }
          }
          ### make dates into day/month/year variabes
          # set the format of the date
          if ("card_date_format" %in% names(data)) {
            if (card_date_format_code %in% c("dd/mm/yyyy", "DD/MM/YYYY")) date_format_string <- "%d/%m/%Y"
            if (card_date_format_code %in% c("ddmmyy")) date_format_string <- "%d%m%Y"
          } else {
            date_format_string <- ifelse(grepl("/", data[, variable, with=FALSE]), "%d/%m/%Y", "%d%B%Y")
          }
          
          # split out
          if (paste0(card_vacc_root, "_day") %in% names(data)) {
            data[is.na(get(paste0(card_vacc_root, "_day"))) & !is.na(get(variable)), paste0(card_vacc_root, "_day") := format(as.Date(get(variable), date_format_string), "%d")]
          } else {
            data[!is.na(get(variable)), paste0(card_vacc_root, "_day") := format(as.Date(get(variable), date_format_string), "%d")]
          }
          if (paste0(card_vacc_root, "_month") %in% names(data)) {
            data[is.na(get(paste0(card_vacc_root, "_month"))) & !is.na(get(variable)), paste0(card_vacc_root, "_month") := format(as.Date(get(variable), date_format_string), "%m")]
          } else {
            data[!is.na(get(variable)), paste0(card_vacc_root, "_month") := format(as.Date(get(variable), date_format_string), "%m")]
          }
          
          if (paste0(card_vacc_root, "_year") %in% names(data)) {
            data[is.na(get(paste0(card_vacc_root, "_year"))) & !is.na(get(variable)), paste0(card_vacc_root, "_year") := format(as.Date(get(variable), date_format_string), "%Y")]
          } else {
            data[!is.na(get(variable)), paste0(card_vacc_root, "_year") := format(as.Date(get(variable), date_format_string), "%Y")]
          }
          # clear out date if subset to day, month, and year
          data[!is.na(get(paste0(card_vacc_root, "_day"))) & !is.na(get(paste0(card_vacc_root, "_month"))) & !is.na(get(paste0(card_vacc_root, "_year"))) &
                 !is.na(get(paste0(card_vacc_root, "_date"))), paste0(card_vacc_root, "_date") := NA_character_]
        }
      }
      # day, month, and year
      card_vacc_root <- paste0("card_", vacc_stem, dose_stem)
      for (date_var in c("day")) {
        # for (date_var in c("day", "month", "year")) {
        # set antigen-dose-specific variable
        variable <- paste0(card_vacc_root, "_", date_var)
        # make numeric
        if (variable %in% names(data)) {
          data[, (variable) := as.character(get(variable))]
          # convert to dates to missing if date value is in card_missing, card_nodate, or card_recall value labels
          codes <- c("NA", ".", "")
          for (card_val in c("missing")) {
            # for (card_val in c("missing", "nodate", "recall", "no")) {
            if (paste0("card_", card_val) %in% names(data)) codes <- c(codes, get(paste0("card_", card_val, "_code")))
          }
          data[get(variable) %in% unique(codes), (variable) := NA_character_]
          # # convert to no (0) if value is in card_no value labels
          for (card_val in c("no")) {
            if (paste0("card_", card_val) %in% names(data)) {
              codes <- get(paste0("card_", card_val, "_code"))
              data[get(variable) %in% c(codes), (variable) := "0"]
            }
          }
          # interpret dates as numbers
          data[, (variable) := as.integer(get(variable))]
          # Buddhist calendar
          if (unique(data$nid) %in% buddhist_nids & date_var=="year") data[, (variable) := get(variable) - 543]
          # Replace incorrect date values with missingness
          if (date_var=="day")   data[get(variable) > 31 & get(variable) != 44, (variable) := NA_integer_] 
        }
      }
      
      ### using card_day, assign doses from card
      if (paste0(card_vacc_root, "_day") %in% names(data)) {
        
        if (dose_num == 1) {
          # Set zeroes with the first dose
          data[get(paste0(card_vacc_root, "_day")) == 0 & is.na(get(paste0(vacc_stem, "_dose"))), 
               c(paste0(vacc_stem, "_card"), paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card")) := list(1, 0, 0)]
          data[get(paste0(card_vacc_root, "_day")) == 0 , paste0(vacc_stem, "_true_zero") := 1]
        }
        
        ### update using card_day only
        data[(!is.na(get(paste0(card_vacc_root, "_day"))) & get(paste0(card_vacc_root, "_day")) != 0) & ((get(paste0(card_vacc_root, "_day")) >= 1 & get(paste0(card_vacc_root, "_day")) <= 31) | get(paste0(card_vacc_root, "_day")) == 44),
             c(paste0(vacc_stem, "_card"), paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card")) := list(1, (dose_num), (dose_num))]
      }
      
      ### replace doses with the card dates if date of vaccination not missing
      # applies to cases we are unable to parse out into day, month, and year
      if (paste0(card_vacc_root, "_date") %in% names(data)) {
        
        if (all(c(paste0(card_vacc_root, "_day")) %in% names(data))) {
          ### if month and year vars missing but date is not missing, replace doses with doses with date
          data[is.na(get(paste0(card_vacc_root, "_day"))) & !is.na(get(paste0(card_vacc_root, "_date"))),
               paste0(vacc_stem, "_card") := 1]
          data[is.na(get(paste0(card_vacc_root, "_day"))) & get(paste0(vacc_stem, "_card")) == 1 &
                 !is.na(get(paste0(card_vacc_root, "_date"))) & get(paste0(card_vacc_root, "_date")) != 0,
               c(paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card")) := list((dose_num), (dose_num))]
          if (dose_num == 1) {
            # data[is.na(get(paste0(card_vacc_root, "_day"))) & get(paste0(card_vacc_root, "_date")) == 0 & get(paste0(vacc_stem, "_card")) == 1,
            #      paste0(vacc_stem, "_dose_from_card") := 0]
            data[is.na(get(paste0(vacc_stem, "_dose"))) & is.na(get(paste0(card_vacc_root, "_day"))) & get(paste0(card_vacc_root, "_date")) == 0,
                 paste0(vacc_stem, "_dose") := 0]
            data[get(paste0(card_vacc_root, "_date")) == 0, paste0(vacc_stem, "_true_zero") := 1]
          }
          
        } else {
          ### if date is not missing, replace doses with doses with date
          data[!is.na(get(paste0(card_vacc_root, "_date"))), paste0(vacc_stem, "_card") := 1]
          data[!is.na(get(paste0(card_vacc_root, "_date"))) & get(paste0(card_vacc_root, "_date")) != 0 & get(paste0(vacc_stem, "_card")) == 1,
               c(paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card")) := list((dose_num), (dose_num))]
          if (dose_num == 1) {
            # data[get(paste0(card_vacc_root, "_date")) == 0 & get(paste0(vacc_stem, "_card")) == 1, paste0(vacc_stem, "_dose_from_card") := 0]
            data[get(paste0(card_vacc_root, "_date")) == 0 & is.na(get(paste0(vacc_stem, "_dose"))), paste0(vacc_stem, "_dose") := 0]
            data[get(paste0(card_vacc_root, "_date")) == 0, paste0(vacc_stem, "_true_zero") := 1]
          }
        }
      }
      
      
      
      
      ###########################
      ### 4. AGE AT VAX       ###
      ###########################
      
      ### split out if single date variable
      if (all(c(paste0("card_", vacc, "_dmy"), "card_date_format") %in% names(data))) {
        ### ...
      }
      
      ### calculate age at vaccination
      age_vars <- c(paste0(card_vacc_root, "_month"), paste0(card_vacc_root, "_year"), "birth_month", "birth_year")
      if (all(age_vars %in% names(data))) {
        
        # make numeric
        for (age_var in age_vars) if (!class(data[, get(age_var)]) %in% c("integer", "numeric")) data[, (age_var) := get(age_var) %>% as.integer]
        
        # check that card dates arent entirely missing
        if (sum(!is.na(data[, get(paste0(card_vacc_root, "_month"))])) > 0 &
            sum(!is.na(data[, get(paste0(card_vacc_root, "_year")) ])) > 0) {
          
          # clean up dates
          # fixing two-digit years
          max_year <- max(data[!is.na(get(paste0(card_vacc_root, "_year"))), (paste0(card_vacc_root, "_year")), with=FALSE], na.rm=TRUE)
          if (max_year < 100) {
            # if year less than 2000, always add 1900 to two-digit years
            data[year_end <  2000, paste0(card_vacc_root, "_year") := 1900 + get(paste0(card_vacc_root, "_year"))]
            # if year_end greater than 2000, allow for years that cross over 2000 by allowing years up to the current year
            data[year_end >= 2000, paste0(card_vacc_root, "_year") := ifelse( (2000 + get(paste0(card_vacc_root, "_year"))) <= current_year,
                                                                              2000 + get(paste0(card_vacc_root, "_year")),
                                                                              1900 + get(paste0(card_vacc_root, "_year")))]
          }
          
          # put logical range on month and year
          month_range <- 1:12
          year_range  <- 1900:current_year
          data[!get(paste0(card_vacc_root, "_month")) %in% month_range, paste0(card_vacc_root, "_month") := NA_integer_]
          data[!get(paste0(card_vacc_root, "_year")) %in% year_range, paste0(card_vacc_root, "_year") := NA_integer_]
          
          # calculate age in months at vaccination
          data[, paste0(vacc_stem, "_", dose_stem, "_age_month") := (12 * (get(paste0(card_vacc_root, "_year")) - 1900) + get(paste0(card_vacc_root, "_month")))
               - (12 * (birth_year - 1900) + birth_month)]
        }
      }
      
      
      ###########################
      ### 6. CLEAN            ###
      ###########################
      
      ### recall gateway (never had any vaccinations)
      if ("ever_vaccinated" %in% names(data)) {
        data[, ever_vaccinated := as.integer(ever_vaccinated)]
        data[is.na(get(paste0(vacc_stem, "_dose"))) & ever_vaccinated == 0, paste0(vacc_stem, "_dose") := 0]
        data[is.na(get(paste0(vacc_stem, "_dose_from_recall"))) & ever_vaccinated == 0, paste0(vacc_stem, "_dose_from_recall") := 0]
      }
    }
  }
  
  
  #########################################################
  ##### combine card and recall doses, for MICS!
  #########################################################
  
  if(unique(data$survey_name) %in% c("UNICEF_MICS")){
    for (vacc_stem in all_stems) {
      
      if (paste0(vacc_stem, "_dose") %in% names(data)) {
        
        #######################################
        #### COMBINE CARD AND RECALL DOSES ####
        #######################################
        
        data[,paste0(vacc_stem, "_dose") := NA_integer_]
        data[!is.na(get(paste0(vacc_stem, "_dose_from_card"))) | !is.na(get(paste0(vacc_stem, "_dose_from_recall"))), paste0(vacc_stem, "_dose") := 0]
        # If just one of card or recall present, use that one
        data[is.na(get(paste0(vacc_stem, "_dose_from_card"))) & !is.na(get(paste0(vacc_stem, "_dose_from_recall"))), paste0(vacc_stem, "_dose") := get(paste0(vacc_stem, "_dose_from_recall"))]
        data[!is.na(get(paste0(vacc_stem, "_dose_from_card"))) & is.na(get(paste0(vacc_stem, "_dose_from_recall"))), paste0(vacc_stem, "_dose") := get(paste0(vacc_stem, "_dose_from_card"))]
        # If both present, use the larger value
        data[(get(paste0(vacc_stem, "_dose_from_recall")) >= get(paste0(vacc_stem, "_dose_from_card"))) & get(paste0(vacc_stem, "_dose")) == 0, paste0(vacc_stem, "_dose") := get(paste0(vacc_stem, "_dose_from_recall"))]
        data[(get(paste0(vacc_stem, "_dose_from_recall")) < get(paste0(vacc_stem, "_dose_from_card"))) & get(paste0(vacc_stem, "_dose")) == 0, paste0(vacc_stem, "_dose") := get(paste0(vacc_stem, "_dose_from_card"))]
      }
      
    }
  }
  
  ### once have all vaccines filled out, fill in card info for missing antigens (if child has card for other vaccines)
  for (vacc_stem in vaccine_stems_in_dataset) {
    
    ###########################
    ### 7. ASSUME CARD      ###
    ###########################
    
    ## if have reported yes with card to other antigens, assume has a card, replace with 0 (FOR FIRST DOSE ONLY)
    for (other_card_hints in vaccine_stems_in_dataset[vaccine_stems_in_dataset != vacc_stem]) {
      if (all(c(paste0(other_card_hints, "_dose"), paste0(other_card_hints, "_card"), paste0(other_card_hints, "_dose_from_card")) %in% names(data))) {
        # assume has card
        data[!is.na(get(paste0(other_card_hints, "_dose"))) & !is.na(get(paste0(other_card_hints, "_dose_from_card"))) & get(paste0(other_card_hints, "_card")) == 1,
             paste0(vacc_stem, "_card") := 1]
        # assume doses 0
        data[!is.na(get(paste0(other_card_hints, "_dose"))) & !is.na(get(paste0(other_card_hints, "_dose_from_card"))) & get(paste0(other_card_hints, "_card")) == 1 &
               is.na(get(paste0(vacc_stem, "_dose"))), paste0(vacc_stem, "_dose") := 0]
        # assume card doses 0
        data[!is.na(get(paste0(other_card_hints, "_dose"))) & !is.na(get(paste0(other_card_hints, "_dose_from_card"))) & get(paste0(other_card_hints, "_card")) == 1 &
               is.na(get(paste0(vacc_stem, "_dose_from_card"))), paste0(vacc_stem, "_dose_from_card") := 0]
      }
    }
    
    ## if the child has no response, no date, but has card --> infer that not vaccinated
    if ("has_vacc_card" %in% names(data) & !unique(data$survey_name) %in% c("UNICEF_MICS")) {
      data[is.na(get(paste0(vacc_stem, "_dose"))) &                # hasn't responded yet (date is missing)
             has_vacc_card == "seen",                              # child has vaccine card that was seen
           c(paste0(vacc_stem, "_dose"), paste0(vacc_stem, "_dose_from_card"), paste0(vacc_stem, "_card")) := list(0, 0, 1)]
    }
    
  }
  
  #***********************************************************************************************************************
  
  
  #----MULTI-DOSE VACCINES------------------------------------------------------------------------------------------------
  
  ### mark multi doses
  # mark DPT, Hib, and HepB if received PENT
  
  if("pent" %in% vaccine_stems_in_dataset){
    unconventional_penta_filepath <- "FILEPATH/penta_unconventional_component.csv"
    
    unconventional_penta_countries <- fread(unconventional_penta_filepath)
    if(unique(data$nid) %in% unconventional_penta_countries$survey_id){
      contains_hepb <- unconventional_penta_countries[survey_id == nid, HepB] == 1
      contains_hib  <- unconventional_penta_countries[survey_id == nid, HiB]  == 1
      contains_ipv  <- unconventional_penta_countries[survey_id == nid, IPV]  == 1
      
      if(contains_hepb & !contains_hib){
        pent_vaccines <- c("dpt", "polio", "hepb")
      } else if (contains_hib & ! contains_hepb){
        pent_vaccines <- c("dpt", "polio", "hib")
      } else {
        mesage(paste0('Dataset contains pentavalent vaccine comprised of unconventional component vaccines. Code is not currently designed
                    to account for this configuration of component vaccines. Please make sure document containing data on unconventional
                    pentavalent vaccines is up-to-date at: ', unconventional_penta_filepath, '. Stopping...'))
        stop()
      }
    } else{
      pent_vaccines <- c("dpt", "hib", "hepb")
    }
  }

  # mark DPT and HepB if received Tetra
  tetra_vaccines <- c("dpt", "hepb")
  # mark MCV and RCV if received MMR
  mmr_vaccines <- c("mcv", "rcv")
  
  ### loop through multi-component vaccines
  for (multi_vacc in c("pent", "tetra", "mmr")) {
    
    ###########################
    ### 1. ASSUME DOSES     ###
    ###########################
    
    # apply to all three dose variables
    for (col_in_loop in c("_dose", "_dose_from_card", "_dose_from_recall")) {
      
      # set MULTI_DOSE_VACC vaccine object
      multi_vacc_column <- paste0(multi_vacc, col_in_loop)
      
      if (multi_vacc %in% vaccine_stems_in_dataset & multi_vacc_column %in% names(data)) {
        
        # loop through antigens included in MULTI_DOSE_VACC
        for (vacc in get(paste0(multi_vacc, "_vaccines"))) {
          
          # set component object
          component_vacc_column <- paste0(vacc, col_in_loop)
          
          if (component_vacc_column %in% names(data)) {
            # if already have some non-missing doses of that component vaccine: if received MULTI_DOSE_VACC and MULTI_DOSE_VACC > vacc, replace each vacc component with MULTI_DOSE_VACC doses
            data[!is.na(get(multi_vacc_column)) & ( (get(multi_vacc_column) > get(component_vacc_column)) | is.na(get(component_vacc_column)) ), (component_vacc_column) := get(multi_vacc_column)]
          } else {
            # if haven't reported receiving doses of that component vaccine, replace each vacc component with MULTI_DOSE_VACC doses
            data[!is.na(get(multi_vacc_column)), (component_vacc_column) := get(multi_vacc_column)]
          }
        }
      }
    }
    
    
    ###########################
    ### 2. APPLY CARD       ###
    ###########################
    
    # apply card indicator
    if (multi_vacc %in% vaccine_stems_in_dataset & paste0(multi_vacc, "_card") %in% names(data)) {
      # loop through antigens included in MULTI_DOSE_VACC
      for (vacc in get(paste0(multi_vacc, "_vaccines"))) {
        # if know card status of combination vaccine, apply that to missing component vaccine card status
        if (paste0(vacc, "_card") %in% names(data)) {
          data[!is.na(paste0(multi_vacc, "_card")) & is.na(get(paste0(vacc, "_card"))), (paste0(vacc, "_card")) := get(paste0(multi_vacc, "_card"))]
        } else {
          data[!is.na(paste0(multi_vacc, "_card")), (paste0(vacc, "_card")) := get(paste0(multi_vacc, "_card"))]
        }
      }
    }
    
    ###########################
    ### 3. AGE AT VAX       ###
    ###########################
    
    # use MULTI_DOSE_VACC age at vaccination if available
    if (multi_vacc %in% vaccine_stems_in_dataset) { 
      for (vacc in get(paste0(multi_vacc, "_vaccines"))) {
        for (dose in c("", "1", "2", "3")) { 
          if (paste0(vacc, dose) %in% vaccines) {
            if (!paste0(vacc, dose, "_age_month") %in% names(data)) data[, paste0(vacc, dose, "_age_month") := NA_integer_]
            if (paste0(multi_vacc, dose, "_age_month") %in% names(data)) {
              data[!is.na(get(paste0(multi_vacc, dose, "_age_month"))) & is.na(get(paste0(vacc, dose, "_age_month"))), paste0(vacc, dose, "_age_month") := get(paste0(multi_vacc, dose, "_age_month"))]
            }
          } 
        }
      }
      
      
      # print
      if (multi_vacc %in% vaccine_stems_in_dataset){
        message(paste0("Combination vaccine components tagged for || ", multi_vacc, " (", paste(get(paste0(multi_vacc, "_vaccines")), collapse=", "), ") ||"))
      }
    }
  }
  
  #***********************************************************************************************************************
  
  
  #----SPECIAL INDICATORS-------------------------------------------------------------------------------------------------
  
  ### DPT3 timeliness ratio
  if (all(c("card_dpt3_day", "card_dpt3_month", "card_dpt3_year") %in% names(data))) {
    data <- get_dpt3_timeliness_ratio(data)
  }
  
  #***********************************************************************************************************************
  
  
  #----CLEANUP------------------------------------------------------------------------------------------------------------
  
  ### go through each antigen (not dose-specific)
  # list all antigens extracted
  outvacc <- c()
  for (vacc_stem in all_stems) {
    if (paste0(vacc_stem, "_dose") %in% names(data)) {
      
      ### replace doses with max dose if beyond
      max <- vacc_dose[var == vacc_stem, refdose]
      data[as.integer(get(paste0(vacc_stem, "_dose"))) > max & !is.na(get(paste0(vacc_stem, "_dose"))), paste0(vacc_stem, "_dose") := max]
      
      ### drop column if no related variables filled out
      levels <- unique(data[, get(paste0(vacc_stem, "_dose"))]) %>% as.integer
      if (all(is.na(levels) | levels=="")) {
        data[, paste0(vacc_stem, "_dose") := NULL]
        data[, paste0(vacc_stem, "_card") := NULL]
        data[, paste0(vacc_stem, "_dose_from_card") := NULL]
        data[, paste0(vacc_stem, "_dose_from_recall") := NULL]
      } else {
        outvacc <- c(outvacc, vacc_stem)
      }
      
      ### create binaries for each dose-specific vaccine
      if (paste0(vacc_stem, "_dose") %in% names(data)) {
        for (dose in 1:max) {
          if (max == 1) newvar <- vacc_stem else newvar <- paste0(vacc_stem, dose)
          data[, (newvar) := ifelse(get(paste0(vacc_stem, "_dose")) >= dose, 1, 0)]
          # drop if nothing but 0s
          if (max(data[, get(newvar)], na.rm=TRUE) == 0) data[, (newvar) := NULL]
          # handle eligibility for recall doses -- if not marked as eligible for this dose and not getting the dose from card, set the binary indicator variable to NA
          if (paste0(vacc_stem, dose, "_recall_eligible") %in% names(data)) {
            data[is.na(get(paste0(vacc_stem, dose, "_recall_eligible"))) & get(newvar) == 1 & (get(paste0(vacc_stem, "_dose_from_card")) < dose | is.na(get(paste0(vacc_stem, "_dose_from_card")))),
                 (newvar) := NA]  
          }
        }
      }
    }
  }
  
  
  
  if (!unique(data$survey_name) %in% c("UNICEF_MICS")){
    
    ### apply recall ever gateway to just first dose (leaving other doses as NA if already missing)
    ### allows us to boost first-dose coverage as appropriate without forcing consecutive doses to 0
    for (vacc_stem in all_stems) {
      
      ### set objects
      # identify whether or not vaccine has multiple doses
      if (grepl(vacc_stem, paste0(one, collapse="|"))) {
        first_dose_stem <- ""
      } else {
        first_dose_stem <- 1
      }
      
      if (paste0(vacc_stem, "_ever") %in% names(data)) {
        
        ### generate variables of interest
        for (var in paste0(vacc_stem, c("_dose", "_card", "_ever", "_dose_from_card", "_dose_from_recall"))) {
          if(!var %in% names(data)) {
            data[, (var) := NA_integer_]
          }
        }
        
        ### if recall ever receiving a dose of a vaccine, code vacc_stem_dose to 1 (else to 0)
        # when none
        data[is.na(get(paste0(vacc_stem, "_dose"))) & get(paste0(vacc_stem, "_ever")) == 0, c(paste0(vacc_stem, first_dose_stem), paste0(vacc_stem, "_dose")) := list(0, 0)]
        data[is.na(get(paste0(vacc_stem, "_dose_from_recall"))) & get(paste0(vacc_stem, "_ever")) == 0, paste0(vacc_stem, "_dose_from_recall") := 0]
        # when at least one
        if (nrow(data[(is.na(get(paste0(vacc_stem, "_dose"))) | get(paste0(vacc_stem, "_dose")) == 0) & get(paste0(vacc_stem, "_ever")) == 1]) > 0) data[, paste0(vacc_stem, "_dose") := as.character(get(paste0(vacc_stem, "_dose")))]
        data[(is.na(get(paste0(vacc_stem, "_dose"))) | get(paste0(vacc_stem, "_dose")) == "0") & get(paste0(vacc_stem, "_ever")) == 1, c(paste0(vacc_stem, first_dose_stem), paste0(vacc_stem, "_dose")) := list(1, "At least 1")]
        if (nrow(data[(is.na(get(paste0(vacc_stem, "_dose_from_recall"))) | get(paste0(vacc_stem, "_dose_from_recall")) == 0) & get(paste0(vacc_stem, "_ever")) == 1]) > 0) data[, paste0(vacc_stem, "_dose_from_recall") := as.character(get(paste0(vacc_stem, "_dose_from_recall")))]
        data[(is.na(get(paste0(vacc_stem, "_dose_from_recall"))) | get(paste0(vacc_stem, "_dose_from_recall")) == "0") & get(paste0(vacc_stem, "_ever")) == 1, paste0(vacc_stem, "_dose_from_recall") := "At least 1"]
        
        ### replace (vacc)_card = 0 if not filled out
        data[(is.na(get(paste0(vacc_stem, "_card"))) | get(paste0(vacc_stem, "_card")) != 1) & !is.na(get(paste0(vacc_stem, "_dose"))),
             paste0(vacc_stem, "_card") := 0]
        
      }
    }
  }
  
  ### print output created by the extraction
  message(paste0("The following vaccines have been processed: ", paste(outvacc, collapse=", ")))
  
  #***********************************************************************************************************************
  
  
  #----RESTRICTIONS-------------------------------------------------------------------------------------------------------
  ### need to check if any records are being dropped here for "filtering" documentation for publications
  ### make complete rotavirus indicator, copying the appropriate column over based on the country's dose schedule
  if ("rota" %in% vaccine_stems_in_dataset) data <- make_rotac(data)
  

  ### age restrictions (drop: missing, younger than 1 year, older than 5 years)
  if (!"age_year" %in% names(data)) data[, age_year := NA_integer_]
  if (!class(data$age_year) %in% c("numeric", "integer")) data[, age_year := as.numeric(age_year)]
  if ("age_month" %in% names(data)) {
    if (!class(data$age_month) %in% c("numeric", "integer")) data[, age_month := as.numeric(age_month)]
    data[is.na(age_year) & !is.na(age_month), age_year := age_month / 12]
  }
  
  ### No data should be dropped prior to this point

  # Calculate dropped children from age missingness and age range
  n_age_data_drop  <- data[is.na(age_year), .N]
  n_age_range_drop <- data[age_year >= 5 | age_year < 1, .N]
  non_def_ages     <- (n_age_data_drop + n_age_range_drop) / nrow(data)
  
  
  # Filter Table: Record data drops from age missingness and age range, save table
  if(record_data_drop) {
    filter_table[, n_with_age_data     := total - n_age_data_drop]
    filter_table[, n_without_age_data  := n_age_data_drop]
    filter_table[, n_12_59_months      := n_with_age_data - n_age_range_drop]
    filter_table[, n_outside_age_range := n_age_range_drop]
    
    if(filter_table[, n_12_59_months] == 0 & filter_table[, missing_indicator] == FALSE){
      message("Filter Table: All data is outside range. Setting remaining filter points to 0")
      filter_table[, c("n_with_mcv_data", "n_without_mcv_data", "n_resampled", "n_not_resampled",
                       "n_after_2000", "n_before_2000", "n_after_outlier",
                       "gps_clusters", "polygon_clusters", "n_gps_located", "n_polygon_located") := list(0)]
      write.csv(filter_table, file = paste0(filter_table_path, nid, ".csv"), row.names = FALSE)
    } else {
      write.csv(filter_table, file = paste0(filter_table_path, nid, ".csv"), row.names = FALSE)
    }
  }
  
  
  # Drop Point 1: Drop data from age missingness and age range
  data <- data[!(age_year >= 5 | age_year < 1 | is.na(age_year))]
  data[, age_year := floor(age_year)]
  
  ### message
  message("Age restrictions complete")
  
  #***********************************************************************************************************************
  
  
  #----END FUNCTION-------------------------------------------------------------------------------------------------------
  
  ### drop unneeded vars
  # keep just needed columns
  metadata_cols <- c("survey_name", "nid", "ihme_loc_id", "year_start", "year_end", "survey_module", "file_path",
                     "smaller_site_unit", "strata", "strata_recode", "psu", "psu_recode", "hh_id", "line_id",
                     "hhweight", paste0("pweight", c("", "_admin_1", "_admin_2", "_admin_3")),
                     paste0("admin_", 1:5), paste0("admin_", 1:5, "_mapped"), paste0("admin_", 1:5, "_id"), paste0("admin_", 1:5, "_urban_id"),
                     "geospatial_id", "latitude", "longitude", "urban",
                     "sex_id", "age_categorical", "age_day", "age_month", "age_year", "int_month", "int_year")
  vaccine_cols  <- c("ever_vaccinated", "has_vacc_card",
                     paste0(all_stems, "_ever"), paste0(all_stems, "_dose"), paste0(all_stems, "_card"),
                     paste0(all_stems, "_dose_from_card"), paste0(all_stems, "_dose_from_recall"),
                     vaccines, paste0(vaccines, "_age_month"),
                     "dpt3_timeliness_ratio",
                     "rotac")
  keep_cols     <- c(metadata_cols, vaccine_cols)[c(metadata_cols, vaccine_cols) %in% names(data)]
  data <- data[, keep_cols, with=FALSE]
   
  # drop columns entirely missing
  for (ea_keep_col in keep_cols) {
    if (all(is.na(data[, get(ea_keep_col)]))) data[, (ea_keep_col) := NULL]
  }
  
  ### logs
  # save vetting log
  vet(data, vaccine_stems_in_dataset=vaccine_stems_in_dataset, vaccines_all_in_dataset=vaccines_all_in_dataset, age_drops=non_def_ages)
  # message
  message("Logs saved")
  
  ### write file
  write.csv(data, file.path(extraction_root, "processed", paste0(nid, ".csv")), row.names=FALSE)
  message(paste0("Processed data for NID ", nid, " saved"))
  
  ### return T/F if function worked
  return(nrow(data) > 1)
  
  ### end of function
}