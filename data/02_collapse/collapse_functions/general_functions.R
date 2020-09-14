#Grabs the most recently created/updated file in a given directory by a given pattern
get_latest_file <- function(folder, use_pattern){
  files <- file.info(list.files(folder, pattern = use_pattern, full.names=TRUE))
  files <- files[with(files, order(as.POSIXct(ctime), decreasing = TRUE)), ]
  latest_postextraction <- unlist(strsplit(rownames(files)[1], '/'))
  latest_postextraction <- latest_postextraction[length(latest_postextraction)]
  message(paste0('You are using file: ', latest_postextraction))
  return(latest_postextraction)
}

#Initial cleaning
#Sets up diagnostic sheet, renames variables, cleans int_month,
#Merges on GBD scalars, removes surveys missing DB,
#Subsetting to relevant years, ages, etc
#Fixes geography info
#Initialized case def variables
initial_cleaning <- function(data){
  message('Starting Initial Cleaning')
  
  message('Setting Up Diagnostic Sheet')
  nids_diag <- unique(lri_data[,c('nid','iso3','year_start')])
  
  message('Renaming and fixing varialbes')
  
  data <- data[,-c('latitude','longitude')]
  setnames(data, old = c('iso3','year_start','int_month','lat','long'),
           new = c('country','start_year','month','latitude','longitude'))
  
  data[,country := substr(country,1,3)]
  
  data[,month := floor(month)]
  
  num_cols <- c('start_year','latitude','longitude')
  data[, (num_cols) := lapply(.SD, function(x){as.numeric(x)}), .SDcols=num_cols]
  
  message('Merge locs & seasonality scalars onto dataset')
  data <- join(data, locs, by = 'country')
  scalar <- scalar[,c('scalar','month','region_name')]
  data <- join(data, scalar, by=c('region_name','month'))
  
  ### CLEAN UP DATA ###
  #Change NA's to 9's in symptom columns
  data[is.na(had_cough), had_cough := 9]
  data[is.na(had_fever), had_fever := 9]
  data[is.na(diff_breathing), diff_breathing := 9]
  data[is.na(chest_symptoms), chest_symptoms := 9]
  
  message('Dropping surveys completely missing Difficulty Breathing responses')
  data$tabulate <- ave(data$diff_breathing, data$nid, FUN= function(x) min(x))
  data <- subset(data, tabulate != 9)
  nids_diag <- nids_diag[(!nid %in% data$nid),reason_dropped := 'No DB']
  
  
  message('Removing outliers and other excluded surveys')
  data <- data[survey_name != "CVD_GEMS",]
  outls <- fread('<<<< FILEPATH REDACTED >>>>')
  data <- data[!nid %in% outls$nid]
  nids_diag[is.na(reason_dropped) & (!nid %in% data$nid),reason_dropped := 'Intentionally excluded']
  
  #create new point column
  message('Create new point column')
  data[, point := NULL]
  data <- data[(!is.na(latitude) & !is.na(longitude)) & latitude != '' & longitude != '', point := 1]
  data <- data[is.na(latitude) & is.na(longitude) & shapefile != '' & location_code != '' & !is.na(shapefile) & !is.na(location_code), point := 0]
  
  message('Removing data without geography info')
  data <- data[!is.na(point),]
  nids_diag[is.na(reason_dropped) & (!nid %in% data$nid),reason_dropped := 'No geodata']
  
  message('Fixing weights, removing invalid weights later')
  data[is.na(pweight) & !is.na(hhweight), pweight := hhweight]
  data[point == 1, pweight := 1]
  data <- data[!is.na(pweight),]
  data <- data[pweight > 0,]
  nids_diag[is.na(reason_dropped) & (!nid %in% data$nid),reason_dropped := 'No or invalid weights']
  
  message('Dropping surveys before 2000')
  data[,true_year := round(weighted.mean(x = int_year, w = pweight)), by = c('nid','country')]
  data <- data[true_year > 1999,]
  nids_diag[is.na(reason_dropped) & (!nid %in% data$nid),reason_dropped := 'Pre 2000 survey (weighted survey year)']
  
  data[,age_year := floor(age_year)]
  data <- data[age_year <= 4,]
  nids_diag[is.na(reason_dropped) & (!nid %in% data$nid),reason_dropped := 'No ages under 5']
  
  setnames(data, old = 'sex_id', new = 'child_sex')
  data$child_sex[is.na(data$child_sex)] = 3
  
  message('Assigning Case Definitions')
  data <- data[, cv_fever := numeric()]
  data <- data[, cv_good := numeric()]
  data <- data[, case := 0]
  
  assign('nids_diag',nids_diag, envir=.GlobalEnv)
  
  return(data)
}

#Assigned case defintions based off of present indicators
assign_case_defs <- function(n){
  temp <- lri_data[nid == n,]
  exist_chest <- ifelse(min(temp$chest_symptoms, na.rm=T) == 0, 1, 0)
  exist_fever <- ifelse(min(temp$had_fever, na.rm=T) == 0, 1, 0)
  if(exist_fever==1 & exist_chest==1) {
    temp <- temp[, cv_fever := 1]
    temp <- temp[, cv_good := 1]
    temp <- temp[had_fever == 1 & chest_symptoms == 1, case := 1]
  } else if (exist_fever==1 & exist_chest==0) {
    temp <- temp[, cv_fever := 1]
    temp <- temp[, cv_good := 0]
    temp <- temp[had_fever == 1 & diff_breathing == 1, case := 1]
  } else if (exist_fever==0 & exist_chest==1) {
    temp <- temp[, cv_fever := 0]
    temp <- temp[, cv_good := 1]
    temp <- temp[chest_symptoms == 1, case := 1]
  } else if (exist_fever==0 & exist_chest==0) {
    temp <- temp[, cv_fever := 0]
    temp <- temp[, cv_good := 0]
    temp <- temp[diff_breathing == 1, case := 1]
  }
  return(temp)
}

#Set up for case definition scalars
lri_scalars <- function(n, data = data){
  message(n)
  temp <- subset(data, nid.new == n)
  #Set indicator = 1 if survey has chest/fever, set 0 if not
  exist_chest <- ifelse(min(temp$chest_symptoms, na.rm=T) == 0, 1, 0)
  exist_fever <- ifelse(min(temp$had_fever, na.rm=T) == 0, 1, 0)
  recall_period <- max(temp$recall_period, na.rm=T)

  #Refine case defs
  temp$good_no_fever <- 0
  temp$poor_no_fever <- 0
  if(exist_fever==1 & exist_chest==1) {
    temp[,overall_lri := ifelse(had_fever==1 & chest_symptoms==1, 1, 0)]
    temp[,good_no_fever := ifelse(chest_symptoms==1,1,0)]
  } else if (exist_fever==1 & exist_chest==0) {
    temp[,overall_lri := ifelse(had_fever==1 & diff_breathing==1, 1, 0)]
    temp[,poor_no_fever := ifelse(diff_breathing==1,1,0)]
  } else if (exist_fever==0 & exist_chest==1) {
    temp[,overall_lri := ifelse(chest_symptoms==1, 1, 0)]
    temp[,good_no_fever := 0]
  } else if (exist_fever==0 & exist_chest==0) {
    temp[,overall_lri := ifelse(diff_breathing==1, 1, 0)]
    temp[,poor_no_fever := 0]
  }
  #Add missing month seasonality scalars
  if(n %in% mm_scalar$nid) {
    mm_scalar = as.data.table(mm_scalar)
    mm_short = mm_scalar[, c('nid', 'country', 'avg_scalar')]
    mm_short[,scalar := as.numeric(NA)] #add on NA column so it only merges on rows missing scalars
    #merge in missing month scalars
    temp = merge(temp, mm_short, by = c('nid', 'country', 'scalar'), all.x = TRUE)
    temp[,scalar := ifelse(!is.na(avg_scalar), avg_scalar, scalar)]
    temp = within(temp, rm(avg_scalar))
  }
  #Apply seasonality scalars
  scalar.dummy <- max(temp$scalar)
  if(is.na(scalar.dummy) | scalar.dummy == -Inf) {
    temp[,scalar_lri := overall_lri]
  } else {
    temp[,scalar_lri := overall_lri*scalar]
    temp[,good_no_fever := good_no_fever*scalar]
    temp[,poor_no_fever := poor_no_fever*scalar]
  }

  #Apply survey design & collapse
  dclus <- svydesign(id=~geospatial_id, weights=~pweight, data=temp)
  prev <- data.table(svyby(~scalar_lri, ~child_sex + age_year, dclus, svymean, na.rm=T))
  prev$base_lri <- svyby(~overall_lri, ~child_sex + age_year, dclus, svymean, na.rm=T)$overall_lri
  prev$good_no_fever <- svyby(~good_no_fever, ~child_sex + age_year, dclus, svymean, na.rm=T)$good_no_fever
  prev$poor_no_fever <- svyby(~poor_no_fever, ~child_sex + age_year, dclus, svymean, na.rm=T)$poor_no_fever
  prev$sample_size <- svyby(~scalar_lri, ~child_sex + age_year, dclus, unwtd.count, na.rm=T)$count

  prev[,cases := scalar_lri * sample_size]
  prev[,country := unique(temp$country)]
  #prev[,location := unique(temp$subname)]
  prev[,start_year := unique(temp$start_year)]
  prev[,end_year := unique(temp$end_year)[1]]
  prev[,cv_diag_valid_good := exist_chest]
  prev[,cv_had_fever := exist_fever]
  prev[,nid := unique(temp$nid)]
  prev[,sex := ifelse(child_sex==2,'Female','Male')]
  prev[,age_start := age_year]
  prev[,age_end := age_year + 1]
  prev[,recall_period := recall_period]
  prev[,survey := unique(temp$survey_name)]
  prev[,notes := paste0('LRI prevalence adjusted for seasonality. The original value was ', round(base_lri,4))]
  
  return(prev)
}

#Generate case definition scalars
cross_walk_scalars <- function(df){
  df[,country := as.character(country)]
  df[,recall_period := as.numeric(recall_period)]
  
  #Replace missing values so they can be used
  df[,missing := ifelse(scalar_lri == 0,1,0)]
  df[,cv_had_fever := ifelse(missing == 1, ifelse(good_no_fever + poor_no_fever != 0, 0, cv_had_fever), cv_had_fever)]
  df[,scalar_lri := ifelse(missing == 1, ifelse(good_no_fever!=0, good_no_fever, poor_no_fever), scalar_lri)]
  
  prevs <- df[,c('scalar_lri','recall_period')]
  prevs[scalar_lri == 0, scalar_lri := 0.0001]
  prevs[scalar_lri > 1, scalar_lri := 1]
  
  #Convert to logit space
  prevs[,ln_lri := logit(scalar_lri)]
  
  prevs[,duration := duration$mean]
  prevs$draw <- inv.logit(prevs$ln_lri)
  prevs[,point := (draw * prevs$duration) / (duration + (recall_period*7) - 1)]
  #Append to original data
  df$mean <- prevs$point
  
  ### GENERATE CASE DEFINITION CROSSWALK COEFFICIENTS ###
  df[,ln_mean := log(mean)]
  
  df[,cv_no_fever := (1-cv_had_fever)]
  
  #Tack on geography information
  locs <- locs[, c('country', 'region_name', 'super_region_name', 'location_id')]
  df[,country := ifelse(country == 'KOSOVO', 'SRB', country)]
  df[,country := ifelse(country == 'KEN_44798', 'KEN', country)]
  df <- join(df, locs, by='country')
  
  
  ##Regress cv_no_fever WITH chest_symptoms##
  mod <- lmer(ln_mean ~ cv_no_fever + (1|age_end) + (1|region_name), data=subset(df, cv_diag_valid_good==1))
  cv_no_fever_good = fixef(mod)[2]
  df$mean <- with(df, ifelse(cv_diag_valid_good==1, ifelse(cv_no_fever==1, exp(ln_mean - cv_no_fever_good),mean), mean))
  
  ##Regress cv_no_fever WITHOUT chest_symptoms##
  mod <- lmer(ln_mean ~ cv_no_fever + (1|age_end) + (1|region_name), data=subset(df, cv_diag_valid_good==0))
  cv_no_fever_poor <- fixef(mod)[2]
  xwalks = c(cv_no_fever_good, cv_no_fever_poor)
  return(xwalks)
}

#actually apply the crosswalking scalars
apply_crosswalks <- function(data = lri_data){
  message('Find and Apply Survey Scalars')
  data$urban_name <- ifelse(data$urban==1,"Urban","Rural")
  data$nid.new <- ifelse(data$country=="ZAF", paste0(data$country,"_",data$nid), 
                        ifelse(data$country=="IND", paste0(data$country,"_",data$urban_name,"_",data$nid),data$nid))
  df <- mclapply(unique(data$nid.new), lri_scalars, data = data, mc.cores = cores)
  df <- rbindlist(df)
  xwalks <- cross_walk_scalars(df)
  message('Fever good: ',xwalks[1])
  message('Fever poor: ',xwalks[2])
  
  #Applying 'fever - good' and 'fever - poor' crosswalks from GBD data
  data <- data[cv_fever == 0 & cv_good == 1 & case != 0, case := case/exp(xwalks[1])]
  data <- data[cv_fever == 0 & cv_good == 0 & case != 0, case := case/exp(xwalks[2])]
  
  #Applying DisMod 'self-report' crosswalk to all values
  data <- data[, case := case/dm.coeffs$`self-report`]
  
  #Applying DisMod 'poor validity' crosswalk to those without chest symptoms
  data <- data[cv_good == 0, case := case/dm.coeffs$poor_validity]
  
  #Applying seasonality scalars
  mm_short = mm_scalar[, c("nid", "country", "avg_scalar")]
  mm_short = mm_short[, scalar := as.numeric(NA)] #add on NA column so it only merges on rows missing scalars
  #merge in missing month scalars
  data = merge(data, mm_short, by = c("nid", "country", "scalar"), all.x = TRUE)
  data = data[!is.na(avg_scalar), scalar := avg_scalar]
  #after pulling in missing month scalars, now drop rows still missing scalars
  no_seas <- unique(data[is.na(scalar),c('nid','country','start_year')])
  write.csv(no_seas,'<<<< FILEPATH REDACTED >>>>', row.names = F)
  data[is.na(scalar), scalar := 1]
  #apply seasonality scalars
  data[,case_pre_seas := case]
  data <- data[, case := case*scalar]
  return(data)
}

#Aggregate lri data
#only thing of note is how year is handled
#multiple years per survey messes with the models a bit
#changed to finding one year per survey based off of the weighted.mean of years
collapse <- function(dt) {
  byvars <- 'start_year,country,location_code,shapefile,survey_series,nid,latitude,longitude,point'
  dt <- unique(dt[, list(int_year = floor(median(int_year, na.rm=T)), #find the median interview year for each location
                         has_lri = weighted.mean(case, pweight), #prevalence
                         N_obs = .N, #orignal number of observations
                         N = sum(pweight)^2/sum(pweight^2), #weighted number of observations
                         sum_of_sample_weights = sum(pweight), #used for Mosser's time series plots
                         recall = unique(recall_period_weeks)), #recall period used per survey
                  by = byvars])
  # replace year with single value per survey
  dt[, year := round(weighted.mean(int_year, N)), by = 'nid']
  # create cluster id
  dt[, cluster_id := .I]
  return(dt)
}

#Bind all CSV files in a given resampling directory
rbind_files <- function(file_loc){
  #grab all resmapled csvs
  resampled <- list.files(file_loc, full.names = T, pattern = '.csv', recursive = F)
  all_poly_data <- mclapply(resampled, fread, mc.cores = cores)
  all_poly_data <- rbindlist(all_poly_data, fill=T, use.names=T)
  
  #change names of columns from resampling
  all_poly_data <- all_poly_data[,-c('latitude','longitude')]
  names(all_poly_data)[names(all_poly_data)=='lat'] <- 'latitude'
  names(all_poly_data)[names(all_poly_data)=='long'] <- 'longitude'
  return(all_poly_data)
}

#Make sure that all surveys made it through resampling, but don't stop code
check_resample <- function(dt, dt_no_resamp){
  dt_no_resamp <- dt_no_resamp[start_year > 1999,]
  dt_no_resamp <- dt_no_resamp[!nid %in% unique(dt$nid),]
  
  resamp_errors <- unique(dt_no_resamp[,c('nid','country','shapefile')])
  if(nrow(resamp_errors) != 0){
    message('THESE SURVEYS DID NOT MAKE IT THROUGH RESAMPLING')
    print(resamp_errors)
  }
  write.csv(resamp_errors, '<<<< FILEPATH REDACTED >>>>')
}

#Final cleanup
#Recall adjustments done here
final_cleanup <- function(){
  message('Binding Resampled Polygons')
  all_poly_data <- rbind_files(file_loc = '<<<< FILEPATH REDACTED >>>>')
  
  #read in coverage data so script can run without starting over
  coverage_data <- get_latest_file(folder_out, 'coverage_data')
  coverage_data <- fread(paste0(folder_out, coverage_data))
  
  #grab all point data
  all_point_data <- coverage_data[point == 1, ]
  all_point_data <- all_point_data[, weight := 1]
  
  message('Binding Points and Polys Back Together')
  all_collapsed <- as.data.table(rbind.fill(all_point_data, all_poly_data))
  
  dur_days <- duration$mean
  all_collapsed$recall <- as.numeric(all_collapsed$recall)
  
  message('Adjusting for Recall Period')
  all_collapsed <- all_collapsed[, rate := has_lri * dur_days / (dur_days + (recall * 7) - 1)]
  all_collapsed <- all_collapsed[, has_lri := rate * N]
  
  message('Checking for dropped surveys in resampling')
  check_resample(all_collapsed, lri_data)
  
  return(all_collapsed)
}