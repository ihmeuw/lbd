##########################################################
###### LBD age crosswalk              
###### This code crosswalks survey data to a specified       
###### age range when reported ages are narrower than the         
###### target (e.g., 0-3 rather than 0-5).                                     
##########################################################

###### Load required packages
age_crosswalk_setup <- function() {
  if (!require("pacman")) install.packages("pacman")
  pacman::p_load(MASS, rio, ggplot2, data.table, plyr, dplyr, gtools)
  
  path <- '<<<< FILEPATH REDACTED >>>>'
  
  ## Install googlesheets, if needed, and then load package; ditto for ggalt, geepack, geeM, MESS
  install.packages("googlesheets", lib=path) ### NOTE: Change path as necessary
  install.packages('ggalt', path)
  library(googlesheets, lib.loc=path)
  library(ggalt, lib.loc=path)
  library(geepack, lib.loc=path)
  library(geeM, lib.loc=path)
  library(MESS, lib.loc=path)
}

###### Get tracking data from Google Sheet. Note: The function currently assumes that a tracking sheet ("sheet_name")
###### with a specific structure is used to identify nids requiring adjustment.
get_tracking_sheet <- function(sheet_name) {
  ### View Google Sheets to which the user has access
  
  ### Retrieve data from tracking sheet
  tr <- gs_title(sheet_name) ### Grab information about tracking sheet, by name
  tracking <- drive_get(id = '<<<< ID REDACTED >>>>') %>%
    dplyr::select(id) %>%
    combine() %>%
    gs_key(lookup = FALSE,
           visibility = "private") %>%
    gs_read_csv(skip=1)
  tracking <- as.data.table(gs_read(ss=tr, ws="Tracking", skip=1)) ### Extract data from tracking sheet; skip the first line (unnecessary headers)
  tracking <- tracking[, c("NID", "Youngest Age", "Oldest Age")]
  colnames(tracking) <- c("nid", "Youngest_age", "Oldest_age")
  
  ### Convert months to years
  reported_as_months <- tracking[grepl("months", Youngest_age)]
  reported_as_months$Youngest_age <- as.numeric(gsub("months", " ", reported_as_months$Youngest_age))/12
  tracking[nid %in% reported_as_months$nid, "Youngest_age"] <- as.character(reported_as_months$Youngest_age)
  
  reported_as_months <- tracking[grepl("months", Oldest_age)]
  reported_as_months$Oldest_age <- as.numeric(gsub("months", " ", reported_as_months$Oldest_age))/12
  tracking[nid %in% reported_as_months$nid, "Oldest_age"] <- as.character(reported_as_months$Oldest_age)
  
  return(tracking)
}

###### Estimate prevalence at age category midpoints so that average prevalence under
###### linear interpolation matches the prevalence reported by GBD as closely as possible.
fit_linear_interpolation_function <- function(intervals_fit) {
  ### Grab intervals with GBD prevalence estimates
  prev_intervals <- intervals_fit[!is.na(midpoint_prevalence)]
  prev_intervals$GBD_cases <- prev_intervals$midpoint_prevalence*prev_intervals$population
  
  ### Grab smallest intervals with GBD population estimates
  pop_intervals <- intervals_fit[(age_group_id != 28) & !((age_group_years_start >= 1) & !(is.na(midpoint_prevalence)))]
  
  ### Optimization function for linear interpolation inflexion points
  optimizeMeasure <- function(vec) {
    ### Create temporary data tables
    temp <- prev_intervals
    temp_pop <- pop_intervals
    
    vec[which(temp$age_group_id == 164)] <- logit(0.00001)
    vec[which(temp$age_group_id == 2)] <- logit(temp[age_group_id == 2, midpoint_prevalence])
    vec[which(temp$age_group_id == 3)] <- logit(temp[age_group_id == 3, midpoint_prevalence])
    
    ### Loop through population intervals and calculate predicted cases
    for (a in 2:nrow(temp_pop)) {
      temp_pop[a, "pred_cases"] <- MESS::auc(prev_intervals$midpoint, inv.logit(vec), from=as.numeric(temp_pop[a, "age_group_years_start"]), to=as.numeric(temp_pop[a, "age_group_years_end"]), type="linear")*temp_pop[a, "population"]/(temp_pop[a, "age_group_years_end"] - temp_pop[a, "age_group_years_start"])
    }
    
    ### Collapse case estimates for GBD categories
    for (b in 1:nrow(temp)) {
      temp[b, "pred_cases"] <- sum(temp_pop[age_group_years_start >= temp[b, age_group_years_start] & age_group_years_end <= temp[b, age_group_years_end], pred_cases])
    }

    ### Return sum of absolute prediction errors, the optimization metric.
    return(sum((abs(temp$pred_cases - temp$GBD_cases)), na.rm=TRUE))
  }
  
  ### Set initial parameters at GBD category prevalence estimates
  InitialGuess <- logit(prev_intervals$midpoint_prevalence)
  InitialGuess[InitialGuess == -Inf] <- logit(0.00001) 
  
  ### Optimize midpoint prevalence estimates to minimize overall prediction error
  Optim <- optim(InitialGuess, optimizeMeasure, control=list(maxit=2000))
  
  ### Adjust prevalence estimates as in optimization function
  vec <- Optim$par

  ### Set prevalence for birth, Early and Late Neonatal stages to their reported levels, to avoid strange fits
  vec[which(prev_intervals$age_group_id == 164)] <- logit(0.00001)
  vec[which(prev_intervals$age_group_id == 2)] <- logit(prev_intervals[age_group_id == 2, midpoint_prevalence])
  vec[which(prev_intervals$age_group_id == 3)] <- logit(prev_intervals[age_group_id == 3, midpoint_prevalence])
  
  prev_intervals$pred <- inv.logit(vec)
  return(prev_intervals)
}

###### Estimate the population in custom age intervals.
calculate_age_interval_populations <- function(intervals) {
  drop_ids <- NA
  
  ### Loop through age intervals and calculate populations for custom age intervals
  for (a in 1:nrow(intervals)) {
    if (!((intervals[a, type] %in% c("GBD Interval", "GBD Single Year")) | (intervals[a, age_group_years_start] == intervals[a, age_group_years_end]) | (!is.na(intervals[a, population])))) {
      message("Calculating population for custom interval...")
      
      ### Identify the smallest enclosing parent population.
      parent_interval <- intervals[!is.na(population) & age_group_years_start <= intervals[a, age_group_years_start] & age_group_years_end >= intervals[a, age_group_years_end]]
      parent_interval <- parent_interval[which(population == min(parent_interval$population))]
      
      ### Calculate population as a fraction of the parent interval, assuming an equal population distribution across the parent interval (as a simplifying assumption).
      intervals[a, "population"] <- parent_interval$population * ((intervals[a, age_group_years_end - age_group_years_start])/(parent_interval$age_group_years_end - parent_interval$age_group_years_start))
      drop_ids <- c(drop_ids, parent_interval$age_group_id)
    }
  }
  
  ### Drop parent intervals (those that are composed of child intervals)
  drop_ids <- unique(drop_ids[!is.na(drop_ids)])
  
  ### Return the revised intervals
  return(intervals[age_group_years_start != age_group_years_end & !(age_group_id %in% drop_ids)])
}

#### Perform age crosswalk for a single nid-year. "fractional_prev_input" indicates whether estimates in data_unadjusted
#### represent prevalence (fractional; =FALSE) or number of positive individuals (=TRUE).
age_crosswalk_by_nid_year <- function(data_unadjusted, min_age_inclusive=0, max_age_exclusive=125, fractional_prev_input=FALSE, adjust_logit_space=FALSE, return_all_GBD_age_categories=FALSE, aggregate_by_nid_year=FALSE, gbd_round_id=5) {
  #### Setup
  temp_nid <- unique(data_unadjusted[, svy_id])
  temp_year <- unique(floor(data_unadjusted[, int_year]))
  location_ID <- location_hierarchy[ihme_loc_id == data_unadjusted[1, country], location_id]
  start_age <- as.numeric(as.character(surveys_non_standard[nid == temp_nid, Youngest_age]))
  end_age <- as.numeric(as.character(surveys_non_standard[nid == temp_nid, Oldest_age]))
  
  message(paste0("Starting age: ", start_age))
  message(paste0("Ending age: ", end_age))
  
  #### Identify nearest GBD years before and after survey year
  year_floor <- unique(floor((floor(temp_year) - 1990)/5)*5 + 1990)
  year_ceiling <- unique(ceiling((floor(temp_year) - 1990)/5)*5 + 1990)
  if (year_floor == 2015) year_floor <- 2010
  if (year_ceiling >= 2015) year_ceiling <- 2017
  
  #### Pull GBD results for location and years
  GBD_results <- sex_collapsed[location_id == location_ID & year_id %in% c(temp_year, year_floor, year_ceiling) & age_group_id %in% age_ids,]
  
  #### Perform linear interpolation of prevalence estimates between bounding GBD year estimates
  years_elapsed <- temp_year - year_floor
  interpolated_by_year <- GBD_results[year_id == temp_year]
  
  for (b in age_ids) {
    interpolated_by_year[age_group_id == b, "mean"] <- ((GBD_results[age_group_id == b & year_id == year_ceiling, mean] - GBD_results[age_group_id == b & year_id == year_floor, mean])/max(c(1, year_ceiling - year_floor))) * years_elapsed + GBD_results[age_group_id == b & year_id == year_floor, mean]
    interpolated_by_year[age_group_id == b, "lower"] <- ((GBD_results[age_group_id == b & year_id == year_ceiling, lower] - GBD_results[age_group_id == b & year_id == year_floor, lower])/max(c(1, year_ceiling - year_floor))) * years_elapsed + GBD_results[age_group_id == b & year_id == year_floor, lower]
    interpolated_by_year[age_group_id == b, "upper"] <- ((GBD_results[age_group_id == b & year_id == year_ceiling, upper] - GBD_results[age_group_id == b & year_id == year_floor, upper])/max(c(1, year_ceiling - year_floor))) * years_elapsed + GBD_results[age_group_id == b & year_id == year_floor, upper]
  }
  
  #### Add birth age group, if needed
  if (min(interpolated_by_year$age_group_years_start) == 0) {
    interpolated_by_year <- rbind(interpolated_by_year[1,], interpolated_by_year)
    interpolated_by_year[1, c("age_group_id", "population", "age_group_years_end", "mean") := .(164, 0, 0, 0)]
  }
  
  #### Define intervals for which population and prevalence need to be calculated
  intervals <- interpolated_by_year[order(age_group_years_start)]
  intervals[is.na(mean), "type"] <- "GBD Single Year"
  intervals[!is.na(mean), "type"] <- "GBD Interval"
  
  ### Identify midpoint age of age interval
  intervals$midpoint <- (intervals$age_group_years_start + intervals$age_group_years_end)/2 
  colnames(intervals)[colnames(intervals) == "mean"] <- "midpoint_prevalence"
  
  ### Estimate prevalence at age category midpoints using an optimization algorithm
  interpolation_points <- fit_linear_interpolation_function(intervals)
  
  #### Add single-year cutpoints
  intervals <- intervals[, c("age_group_id", "midpoint_prevalence", "age_group_years_start", "age_group_years_end", "type", "population", "midpoint")]
  intervals <- rbind(intervals, data.table(age_group_id=NA, midpoint_prevalence=NA, age_group_years_start=(max(1, start_age)):end_age, age_group_years_end=(max(1, start_age)):end_age, type="Single-Year Cutpoints", population=NA, midpoint=(max(1, start_age)):end_age))
  
  #### Add additional cut points based on non-standard age ranges
  if (nrow(intervals[age_group_years_start == start_age & age_group_years_end == start_age]) == 0) {
    previous_pt <- max(intervals[age_group_years_end < start_age, age_group_years_end])
    next_pt <- min(intervals[age_group_years_start > start_age, age_group_years_start])
    intervals <- rbind(intervals, data.table(age_group_id=NA, midpoint_prevalence=NA, age_group_years_start=start_age, age_group_years_end=start_age, type="Study-level Boundary", population=NA, midpoint=start_age))
    intervals <- rbind(intervals, data.table(age_group_id=NA, midpoint_prevalence=NA, age_group_years_start=previous_pt, age_group_years_end=start_age, type="Study-level Interval", population=NA, midpoint=(previous_pt + start_age)/2))
    intervals <- rbind(intervals, data.table(age_group_id=NA, midpoint_prevalence=NA, age_group_years_start=start_age, age_group_years_end=next_pt, type="Study-level Interval", population=NA, midpoint=(next_pt + start_age)/2))
  }
  if (nrow(intervals[age_group_years_start == end_age & age_group_years_end == end_age]) == 0) {
    previous_pt <- max(intervals[age_group_years_end < end_age, age_group_years_end])
    next_pt <- min(intervals[age_group_years_start > end_age, age_group_years_start])
    intervals <- rbind(intervals, data.table(age_group_id=NA, midpoint_prevalence=NA, age_group_years_start=end_age, age_group_years_end=end_age, type="Study-level Boundary", population=NA, midpoint=end_age))
    intervals <- rbind(intervals, data.table(age_group_id=NA, midpoint_prevalence=NA, age_group_years_start=previous_pt, age_group_years_end=end_age, type="Study-level Interval", population=NA, midpoint=(previous_pt + end_age)/2))
    intervals <- rbind(intervals, data.table(age_group_id=NA, midpoint_prevalence=NA, age_group_years_start=end_age, age_group_years_end=next_pt, type="Study-level Interval", population=NA, midpoint=(next_pt + end_age)/2))
  }
  
  #### Sort intervals and add a unique ID
  intervals <- intervals[order(midpoint)]
  intervals$unique_id <- 1:nrow(intervals)
  
  #### Cycle through rows and calculate populations for custom intervals
  intervals <- calculate_age_interval_populations(intervals)
  
  #### Drop intervals that are no longer needed
  intervals <- intervals[!(type == "GBD Single Year" & age_group_years_start < 1) & !(type == "GBD Single Year" & age_group_years_start >= max_age_exclusive) & !(type == "GBD Interval" & age_group_years_start >= 1)]
  
  #### Determine baseline prevalence for each remaining age interval
  intervals$interpolated_prevalence <- intervals$midpoint_prevalence ### Use GBD estimates directly where available
  for (a in 1:nrow(intervals)) { ### Otherwise estimate prevalence using inferred interpolation points and AUC
    if (is.na(intervals[a, "interpolated_prevalence"])) {
      intervals[a, "interpolated_prevalence"] <- ((MESS::auc(interpolation_points[, midpoint], (interpolation_points[, pred]), from=as.numeric(intervals[a, age_group_years_start]), to=as.numeric(intervals[a, age_group_years_end]), type="linear"))/(intervals[a, age_group_years_end] - intervals[a, age_group_years_start]))
    }
  }
  
  #### Calculate baseline population-weighted prevalence for study interval and for Under-5
  start_category <- intervals[age_group_years_start == start_age, unique_id]
  end_category <- intervals[age_group_years_end == end_age, unique_id]
  baseline_prev_study_interval <- sum(intervals[age_group_years_start >= start_age & age_group_years_end <= end_age, interpolated_prevalence*population])/sum(intervals[age_group_years_start >= start_age & age_group_years_end <= end_age, population])
  baseline_prev_target <- sum(intervals[age_group_years_start >= min_age_inclusive & age_group_years_end <= max_age_exclusive, interpolated_prevalence*population])/sum(intervals[age_group_years_start >= min_age_inclusive & age_group_years_end <= max_age_exclusive, population])
  
  #### Determine scaling factor and calculate adjusted prevalence
  #temp_data <- non_standard_data[svy_id == temp_nid & floor(int_year) == temp_year,]
  temp_data <- data_unadjusted
  
  if (fractional_prev_input) { # had_indicator represents prevalence
    if (adjust_logit_space) { # Calculate adjustment in logit space
      temp_data$scaling <- logit(temp_data$had_indicator) - logit(baseline_prev_study_interval)
      temp_data$adjusted_had_indicator <- inv.logit(temp_data$scaling + logit(baseline_prev_target))
    } else {  # Calculate adjustment in probability space
      temp_data$scaling <- temp_data$had_indicator/baseline_prev_study_interval
      temp_data$adjusted_had_indicator <- temp_data$scaling*baseline_prev_target
    }
    temp_data[adjusted_had_indicator > 1, adjusted_had_indicator := 1]
  } else { # had_indicator represents positive individuals
    if (adjust_logit_space) { # Calculate adjustment in logit space
      temp_data$scaling <- logit(temp_data$had_indicator/temp_data$N) - logit(baseline_prev_study_interval)
      temp_data$adjusted_had_indicator <- inv.logit(temp_data$scaling + logit(baseline_prev_target))
    } else { # Calculate adjustment in probability space
      temp_data$scaling <- (temp_data$had_indicator/temp_data$N)/baseline_prev_study_interval
      temp_data$adjusted_had_indicator <- temp_data$scaling*baseline_prev_target
    }
    temp_data[adjusted_had_indicator > N, adjusted_had_indicator := N]
  }
  
  all_GBD_categories <- NULL
  GBD_adjusted <- NULL
  print(temp_data[, c("had_indicator", "scaling", "adjusted_had_indicator")])
  
  if (!aggregate_by_nid_year) {
    final_scaling <- NULL
  } else {
    final_scaling <- logit(baseline_prev_target) - logit(baseline_prev_study_interval)
    print(paste0("Final scaling (baseline_prev_target - baseline_prev_study_interval): ", final_scaling))
  }
  
  if (return_all_GBD_age_categories) { # Optionally retrieve baseline prevalence and calculate adjusted prevalence for each GBD age category
    age_metadata <- get_age_metadata(age_group_set_id=12, gbd_round_id=gbd_round_id)
     
    #### Retrieve GBD prevalence estimates
    prev_by_age <- as.data.table(get_model_results("epi", gbd_id=modelable_entity, measure_id=5, location_id=location_ID, location_set_id=35, year_id=c(temp_year, year_floor, year_ceiling), age_group_id=age_metadata$age_group_id, sex_id='all', status='best', gbd_round_id=gbd_round_id))
    
    #### Perform linear interpolation of prevalence estimates between bounding GBD year estimates
    years_elapsed <- temp_year - year_floor
    interpolated_by_year <- prev_by_age[year_id == year_floor]
    
    for (b in age_ids) {
      interpolated_by_year[age_group_id == b, "mean"] <- ((prev_by_age[age_group_id == b & year_id == year_ceiling, mean] - prev_by_age[age_group_id == b & year_id == year_floor, mean])/max(c(1, year_ceiling - year_floor))) * years_elapsed + prev_by_age[age_group_id == b & year_id == year_floor, mean]
      interpolated_by_year[age_group_id == b, "lower"] <- ((prev_by_age[age_group_id == b & year_id == year_ceiling, lower] - prev_by_age[age_group_id == b & year_id == year_floor, lower])/max(c(1, year_ceiling - year_floor))) * years_elapsed + prev_by_age[age_group_id == b & year_id == year_floor, lower]
      interpolated_by_year[age_group_id == b, "upper"] <- ((prev_by_age[age_group_id == b & year_id == year_ceiling, upper] - prev_by_age[age_group_id == b & year_id == year_floor, upper])/max(c(1, year_ceiling - year_floor))) * years_elapsed + prev_by_age[age_group_id == b & year_id == year_floor, upper]
    }
    
    interpolated_by_year$year_id <- temp_year
    
    #### Retrieve GBD population estimates
    popNumbers <- merge(data.table(get_population(age_group_id=prev_by_age$age_group_id, location_id=location_ID, year_id=temp_year, sex_id=1:3, single_year_age=F)), age_metadata)
    
    #### Sex-collapse GBD prevalence estimates
    pop_prev <- as.data.table(merge(popNumbers, interpolated_by_year, by=c("age_group_id", "sex_id", "year_id"), all=TRUE))
    sex_collapsed <- pop_prev[sex_id == 3]
    for (a in 1:nrow(sex_collapsed)) {
      temp <- pop_prev[sex_id %in% 1:2 & year_id == sex_collapsed[a, year_id] & age_group_id == sex_collapsed[a, age_group_id]]
      sex_collapsed[a, mean := sum(temp[, population * mean])/sum(temp$population)]
    }
    
    all_GBD_categories <- sex_collapsed
    
    #### Determine scaling factor and calculate adjusted prevalence
    temp_GBD_data <- non_standard_data[svy_id == temp_nid & floor(int_year) == temp_year,]
    GBD_adjusted <- vector("list", nrow(temp_GBD_data))
      
    for (c in 1:length(GBD_adjusted)) {
      GBD_adjusted[[c]] <- sex_collapsed
      if (adjust_logit_space) { # Calculate adjustment in logit space
          scaling <- logit(temp_GBD_data[c, had_indicator]) - logit(baseline_prev_study_interval)
          GBD_adjusted[[c]]$adjusted_mean <- inv.logit(logit(GBD_adjusted[[c]]$mean) + scaling)
       } else {  # Calculate adjustment in probability space
          scaling <- temp_GBD_data[c, had_indicator]/baseline_prev_study_interval
          GBD_adjusted[[c]]$adjusted_mean <- GBD_adjusted[[c]]$mean*scaling
      }
    }
  }
  
  return(list("adjusted_had_indicator"=temp_data$adjusted_had_indicator, "all_GBD_categories"=all_GBD_categories, "GBD_adjusted"=GBD_adjusted, "scaling"=final_scaling))
}

#### Retrieve population and prevalence data from GBD
retrieve_GBD_pop_prev <- function(non_standard_data, min_age_inclusive=0, max_age_exclusive=125, gbd_round_id=5, modelable_entity=1181) {
  #### Get location hierarchy and age metadata from central database
  location_hierarchy <- get_location_metadata(location_set_id=35, gbd_round_id=gbd_round_id)
  age_metadata <- get_age_metadata(age_group_set_id=12, gbd_round_id=gbd_round_id)
  non_standard_data <- merge(non_standard_data, location_hierarchy[, c("ihme_loc_id", "location_id")], by.x="country", by.y="ihme_loc_id", all.x=TRUE)
  
  #### Identify location IDs
  location_IDs <- location_hierarchy[ihme_loc_id %in% unique(non_standard_data$country), location_id]
  
  #### Identify nearest GBD years before and after survey years
  year_floor <- unique(floor((non_standard_data$int_year - 1990)/5)*5 + 1990)
  year_ceiling <- unique(ceiling((non_standard_data$int_year - 1990)/5)*5 + 1990)
  year_floor[year_floor == 2015] <- 2010
  year_ceiling[year_ceiling == 2015] <- 2017
  year_IDs <- unique(c(year_floor, year_ceiling))
  
  #### Detect needed age categories; retrieve up to 3 categories below start age and 3 categories above end age to improve interpolation
  age_metadata$unique_id <- as.integer(rownames(age_metadata))
  age_groups_row_start <- max(1, max(age_metadata[age_group_years_start <= min_age_inclusive, unique_id]) - 3)
  age_groups_row_end <- min(nrow(age_metadata), min(age_metadata[age_group_years_end >= max_age_exclusive, unique_id]) + 3)
  age_groups_row_start2 <- max(1, max(age_metadata[age_group_years_start <= min(as.numeric(as.character(surveys_non_standard$Youngest_age))), unique_id]) - 3)
  age_groups_row_end2 <- min(nrow(age_metadata), min(age_metadata[age_group_years_end >= max(as.numeric(as.character(surveys_non_standard$Oldest_age))), unique_id]) + 3)
  age_ids1 <- age_metadata[age_groups_row_start:age_groups_row_end, age_group_id]
  age_ids2 <- age_metadata[age_groups_row_start2:age_groups_row_end2, age_group_id]
  age_ids <- sort(unique(c(age_ids1, age_ids2)))
  
  #### Retrieve population data for all non-standard surveys
  non_standard_data$int_year <- floor(non_standard_data$int_year)
  year_IDs_exact <- unique(c(year_IDs, non_standard_data$int_year))
  popNumbers1 <- merge(data.table(get_population(age_group_id=age_ids, location_id=location_IDs, year_id=year_IDs_exact, sex_id=1:3, single_year_age=F, gbd_round_id=gbd_round_id)), age_metadata[age_groups_row_start:age_groups_row_end,])
  popNumbers2 <- merge(data.table(get_population(age_group_id=age_ids, location_id=location_IDs, year_id=year_IDs_exact, sex_id=1:3, single_year_age=F, gbd_round_id=gbd_round_id)), age_metadata[age_groups_row_start2:age_groups_row_end2,])
  popNumbers <- unique(rbind(popNumbers1, popNumbers2))
  
  #### Add single-year age categories
  age_ids_single_years <- data.table(age_group_id=c(28, 49:147, 295:300, 302:305, 346:355, 366:370), age_group_years_start=c(0:124), age_group_years_end=c(1:125))
  age_ids_single_years1 <- age_ids_single_years[age_group_years_start >= min(popNumbers1$age_group_years_start) & age_group_years_end <= max(popNumbers1$age_group_years_end)]
  age_ids_single_years2 <- age_ids_single_years[age_group_years_start >= min(popNumbers2$age_group_years_start) & age_group_years_end <= max(popNumbers2$age_group_years_end)]
  popNumbers_single_years1 <- merge(data.table(get_population(age_group_id=age_ids_single_years1$age_group_id, location_id=location_IDs, year_id=year_IDs_exact, sex_id=1:3, single_year_age=T, gbd_round_id=gbd_round_id)), age_ids_single_years1)
  popNumbers_single_years2 <- merge(data.table(get_population(age_group_id=age_ids_single_years2$age_group_id, location_id=location_IDs, year_id=year_IDs_exact, sex_id=1:3, single_year_age=T, gbd_round_id=gbd_round_id)), age_ids_single_years2)
  popNumbers_single_years <- unique(rbind(popNumbers_single_years1, popNumbers_single_years2))
  popNumbers <- rbind(popNumbers[, c("age_group_id", "location_id", "year_id", "sex_id", "population", "age_group_years_start", "age_group_years_end")], popNumbers_single_years[, c("age_group_id", "location_id", "year_id", "sex_id", "population", "age_group_years_start", "age_group_years_end")])
  
  age_ids <- sort(unique(popNumbers$age_group_id))
  
  #### Retrieve GBD prevalence estimates
  prev_by_age <- get_model_results("epi", gbd_id=modelable_entity, measure_id=5, location_id=location_IDs, location_set_id=35, year_id=year_IDs_exact, age_group_id=age_ids, sex_id='all', status='best', gbd_round_id=gbd_round_id)
  
  #### Sex-collapse GBD prevalence estimates
  pop_prev <- merge(popNumbers, prev_by_age, all=TRUE)
  sex_collapsed <- pop_prev[sex_id == 3]
  for (a in 1:nrow(sex_collapsed)) {
    temp <- pop_prev[sex_id %in% 1:2 & location_id == sex_collapsed[a, location_id] & year_id == sex_collapsed[a, year_id] & age_group_id == sex_collapsed[a, age_group_id]]
    sex_collapsed[a, mean := sum(temp[, population * mean])/sum(temp$population)]
  }
  
  return(list("non_standard_data"=non_standard_data, "sex_collapsed"=sex_collapsed, "popNumbers"=popNumbers, "location_hierarchy"=location_hierarchy, "age_metadata"=age_metadata, "age_ids"=age_ids))
}
