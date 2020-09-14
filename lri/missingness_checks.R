rm(list = ls())
library(plyr)
library(dplyr)
library(data.table)

#set arguments
by_cluster <- TRUE #plot missingness by cluster
new_dir <- TRUE #make a new directory with today's date?
missing_lim <- 0.25 #set limit for cluster-level missingness

#make directory
if (new_dir == TRUE){
  today <- Sys.Date()
  today <- gsub("-", "_", today)
  dir.create('<<<< FILEPATH REDACTED >>>>')
}

################################################################
# (1) load the collapse
################################################################
##############################################
# General setup
##############################################
cores <- 4
repo <- '<<<< FILEPATH REDACTED >>>>'
setwd(repo)

module_date <- Sys.Date()
module_date <- gsub("-", "_", module_date)
folder_in <- '<<<< FILEPATH REDACTED >>>>'
folder_out <- '<<<< FILEPATH REDACTED >>>>'

source('collapse_functions/general_functions.R')
source('<<<< FILEPATH REDACTED >>>>/mbg_central/setup.R')
package_list <- c('mgcv', 'nlme', 'lme4', 'survey', 'foreign', 'rgeos', 'data.table','raster',
                  'rgdal','seegSDM','seegMBG','plyr','dplyr', 'foreach', 'doParallel', 'feather','boot')
mbg_setup(package_list = package_list, repos='<<<< FILEPATH REDACTED >>>>')

#shared functions
source('<<<< FILEPATH REDACTED >>>>/get_outputs.R')
source('<<<< FILEPATH REDACTED >>>>/get_draws.R')
source('<<<< FILEPATH REDACTED >>>>/get_ids.R')
source('<<<< FILEPATH REDACTED >>>>/get_population.R')
source('<<<< FILEPATH REDACTED >>>>/get_location_metadata.R')

#libraries
library(ggplot2)
library(ggrepel)

##############################################
# Load and Clean Data
##############################################
latest_postextraction <- get_latest_file(folder_in, '*.feather')

#load most recent geomatched data
lri_data <- data.table(read_feather(paste0(folder_in, latest_postextraction)))

#Read in supplementary information
locs <- fread('<<<< FILEPATH REDACTED >>>>') #location metadata
locs <- locs[,c('location_name','location_id','super_region_name','super_region_id',
                'region_name','region_id','parent_id', 'ihme_loc_id')]
setnames(locs, old = c('ihme_loc_id'), new = c('country'))
scalar <- fread('<<<< FILEPATH REDACTED >>>>') #seasonality scalars
mm_scalar <- data.table(read.csv('<<<< FILEPATH REDACTED >>>>')) #seasonality scalars & info from surveys without month microdata
dm.coeffs <- fread('<<<< FILEPATH REDACTED >>>>') #DisMod crosswalk coefficients, for 2016 (not 2015)
duration <- fread('<<<< FILEPATH REDACTED >>>>') #LRI duration draws (FOR GBD 2016)

#Initial cleaning: Fixing/renaming columns, subsetting ages, prep case definitions
lri_data <- initial_cleaning(lri_data)

#assign cluster_id to lri_data if plotting missingness by cluster
if (by_cluster == TRUE){
  #collapse to obtain index of cluster_ids
  lri_collapse <- mclapply(unique(lri_data$nid), assign_case_defs, mc.cores=cores)
  lri_collapse <- rbindlist(lri_collapse)
  lri_collapse <- apply_crosswalks(lri_collapse)
  lri_collapse <- collapse(lri_collapse)
  
  #find index for the polygons by cluster_id
  lri_data_poly <- filter(lri_data, point == 0)
  index <- dplyr::select(lri_collapse, c('location_code','nid', 'cluster_id','point','survey_series','shapefile')) %>%
    filter(point == 0)
  lri_data_poly_cluster <- merge(lri_data_poly, index, by = c('location_code','nid','survey_series','shapefile', 'point'), all.x = T)
  
  #find index for points by cluster_id
  lri_data_point <- filter(lri_data, point == 1)
  index <- dplyr::select(lri_collapse, c('latitude','longitude', 'cluster_id','point','survey_series', 'nid','location_code')) %>%
    filter(point == 1)
  lri_data_point_cluster <- merge(lri_data_point, index, by = c('latitude','longitude','point','survey_series', 'nid','location_code'), all.x = T)
  
  #bind points and polygons
  lri_data_cluster <- rbind(lri_data_point_cluster, lri_data_poly_cluster)
  if (nrow(lri_data_cluster) != nrow(lri_data)){
    message('lri_data_cluster does not have same number of rows as lri_data--look into this')
  }
  lri_data <- lri_data_cluster
}

#subset to surveys in input data
input_data <- fread('<<<< FILEPATH REDACTED >>>>')
nids_to_keep <- unique(input_data$nid)
lri_data <- filter(lri_data, nid %in% nids_to_keep)
lri_data_save <- lri_data
if (by_cluster == TRUE){
  missing_cluster_id <- dplyr::filter(lri_data, is.na(cluster_id))
  ids <- unique(missing_cluster_id$nid)
  if (length(ids != 0)){
    ids <- paste(ids, sep="", collapse=", ") 
    message(paste0('nids ', ids, ' are getting dropped in collapse currently and are not assigned a cluster_id'))
  }
}

#figure out if cough, db, chest, and fever asked at all
lri_data$asked_db <- NA
lri_data$asked_cough <- NA
lri_data$asked_chest <- NA
lri_data$asked_fever <- NA

nids <- unique(lri_data$nid)
for (survey in nids){
  print(survey)
  nid_data <- filter(lri_data, nid == survey)
  
  if (all(nid_data$diff_breathing == 9)){
    lri_data <- mutate(lri_data, asked_db = replace(asked_db, nid == survey, 0))
  } else {
    lri_data <- mutate(lri_data, asked_db = replace(asked_db, nid == survey, 1))
  }
  
  if (all(nid_data$had_cough == 9)){
    lri_data <- mutate(lri_data, asked_cough = replace(asked_cough, nid == survey, 0))
  } else {
    lri_data <- mutate(lri_data, asked_cough = replace(asked_cough, nid == survey, 1))
  }
  
  if (all(nid_data$had_fever == 9)){
    lri_data <- mutate(lri_data, asked_fever = replace(asked_fever, nid == survey, 0))
  } else {
    lri_data <- mutate(lri_data, asked_fever = replace(asked_fever, nid == survey, 1))
  }
  
  if (all(nid_data$chest_symptoms == 9)){
    lri_data <- mutate(lri_data, asked_chest = replace(asked_chest, nid == survey, 0))
  } else {
    lri_data <- mutate(lri_data, asked_chest = replace(asked_chest, nid == survey, 1))
  }
  
}

#make missingness indicator: 0 = asked, answered | 1 = asked, not answered | 9 = not asked, not answered
lri_data$db_missing <- NA
lri_data$cough_missing <- NA
lri_data$fever_missing <- NA
lri_data$chest_missing <- NA

#if sx_asked == 0, set missing indicator to 9
lri_data <- mutate(lri_data, db_missing = replace(db_missing, asked_db == 0, 9))
lri_data <- mutate(lri_data, cough_missing = replace(cough_missing, asked_cough == 0, 9))
lri_data <- mutate(lri_data, fever_missing = replace(fever_missing, asked_fever == 0, 9))
lri_data <- mutate(lri_data, chest_missing = replace(chest_missing, asked_chest == 0, 9))

#if has_sx == 9 & sx_asked == 1, indicator to 1
lri_data <- mutate(lri_data, db_missing = replace(db_missing, diff_breathing == 9 & asked_db == 1, 1))
lri_data <- mutate(lri_data, cough_missing = replace(cough_missing, had_cough == 9 & asked_cough == 1, 1))
lri_data <- mutate(lri_data, fever_missing = replace(fever_missing, had_fever == 9 & asked_fever == 1, 1))
lri_data <- mutate(lri_data, chest_missing = replace(chest_missing, chest_symptoms == 9 & asked_chest == 1, 1))

#if has_sx == 1 or 0 & sx_asked == 1, indicator to 0
lri_data <- mutate(lri_data, db_missing = replace(db_missing, diff_breathing == 1 | diff_breathing == 0 & asked_db == 1, 0))
lri_data <- mutate(lri_data, cough_missing = replace(cough_missing, had_cough == 1 | had_cough == 0 & asked_cough == 1, 0))
lri_data <- mutate(lri_data, fever_missing = replace(fever_missing, had_fever == 1 | had_fever == 0 & asked_fever == 1, 0))
lri_data <- mutate(lri_data, chest_missing = replace(chest_missing, chest_symptoms == 1 | chest_symptoms == 0 & asked_chest == 1, 0))

#this data took a long time to produce: save a backup csv just in case
write.csv(lri_data, '<<<< FILEPATH REDACTED >>>>')

############################################################################
# Calculate missingness by NID and plot
############################################################################
if (by_cluster == F){
  #make proprotion missing table by nid
  #start an empty table for results by nid
  missing_table <- list()
  missing_table$nid <- unique(lri_data$nid)
  missing_table$prop_cough_missing <- NA
  missing_table$prop_db_missing <- NA
  missing_table$prop_chest_missing <- NA
  missing_table$prop_fever_missing <- NA
  missing_table <- as.data.frame(missing_table)
  
  nids <- unique(lri_data$nid)
  for (survey in nids){
    print(survey)
    nid_data <- filter(lri_data, nid == survey)
    
    cough_data <- filter(nid_data, diff_breathing == 1 | diff_breathing == 9)
    missing <- filter(cough_data, cough_missing == 1)
    answered <- filter(cough_data, cough_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$nid == survey,]$prop_cough_missing <- m_prop
    
    db_data <- filter(nid_data, had_cough == 1 | had_cough == 9)
    missing <- filter(db_data, db_missing == 1)
    answered <- filter(db_data, db_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$nid == survey,]$prop_db_missing <- m_prop
    
    case_data <- filter(nid_data, (had_cough == 1 | had_cough == 9) & (diff_breathing == 1 | diff_breathing == 9))
    missing <- filter(case_data, chest_missing == 1)
    answered <- filter(case_data, chest_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$nid == survey,]$prop_chest_missing <- m_prop
    
    case_data <- filter(nid_data, (had_cough == 1 | had_cough == 9) & (diff_breathing == 1 | diff_breathing == 9))
    missing <- filter(case_data, fever_missing == 1)
    answered <- filter(case_data, fever_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$nid == survey,]$prop_fever_missing <- m_prop
  }
  
  #save in case of future errors
  write.csv(missing_table, '<<<< FILEPATH REDACTED >>>>')
  
  ######################################################################
  # (2) Make plots
  ######################################################################
  
  ## Plot prevalence vs. missing for each survey, case defintion, year w/ GBD 95% CI for reference
  #collapse input_data by nid (analogous to TS plots)
  nid_table <- input_data %>%
    group_by(nid, start_year, country, point, survey_series) %>%
    dplyr::summarize(prev = weighted.mean(x = has_lri/N, w = N, na.rm = T),
                     N = sum(N*weight, na.rm = T))
  
  #merge prevalences onto the input table
  nid_table <- merge(nid_table, missing_table, by = 'nid') %>%
    as.data.frame()
  
  #get GBD results
  #location metadata
  locs = get_location_metadata(location_set_id = 9, gbd_round_id = 5)
  
  #gbd data for LRI
  gbd_data <- get_outputs('cause',
                          cause_id = 322, 
                          measure_id = 5, 
                          metric_id = 3, 
                          year_id = c(2000:2017), 
                          location_id = locs[,location_id], 
                          age_group_id = 1, 
                          sex_id = 3, 
                          gbd_round_id = 5, 
                          compare_version_id = 377)
  
  #merge on locs info to nid_table
  locs <- select(locs, c('ihme_loc_id', 'location_name', 'location_id'))
  nid_table <- merge(nid_table, locs, by.x = 'country', by.y = 'ihme_loc_id')
  
  loc_ids <- unique(nid_table$location_id)
  
  #open pdf
  pdf(('<<<< FILEPATH REDACTED >>>>'))
  
  #loop over countries and plot
  for (loc in loc_ids){
    loc_data <- filter(nid_table, location_id == loc)
    country_name <- unique(loc_data$location_name)
    print(country_name)
    plot <- ggplot(data = loc_data, aes(x = prop_cough_missing, y = prev)) + geom_point(aes(color = start_year, size = N, shape = survey_series)) + 
      ggtitle(paste0(country_name, ' | cough')) + ylab('LRI prevalence') + 
      geom_text_repel(aes(x = prop_cough_missing, y = prev, label = nid), size = 2)
    plot(plot)
    
    loc_data <- filter(nid_table, location_id == loc)
    country_name <- unique(loc_data$location_name)
    plot <- ggplot(data = loc_data, aes(x = prop_db_missing, y = prev)) + geom_point(aes(color = start_year, size = N, shape = survey_series)) + 
      ggtitle(paste0(country_name, ' | difficulty breathing')) + ylab('LRI prevalence') + 
      geom_text_repel(aes(x = prop_db_missing, y = prev, label = nid), size = 2)
    plot(plot)
    
    loc_data <- filter(nid_table, location_id == loc)
    country_name <- unique(loc_data$location_name)
    plot <- ggplot(data = loc_data, aes(x = prop_fever_missing, y = prev)) + geom_point(aes(color = start_year, size = N, shape = survey_series)) + 
      ggtitle(paste0(country_name, ' | fever')) + ylab('LRI prevalence') + 
      geom_text_repel(aes(x = prop_fever_missing, y = prev, label = nid), size = 2)
    plot(plot)
    
    loc_data <- filter(nid_table, location_id == loc)
    country_name <- unique(loc_data$location_name)
    plot <- ggplot(data = loc_data, aes(x = prop_chest_missing, y = prev)) + geom_point(aes(color = start_year, size = N, shape = survey_series)) + 
      ggtitle(paste0(country_name, ' | chest')) + ylab('LRI prevalence') +
      geom_text_repel(aes(x = prop_chest_missing, y = prev, label = nid), size = 2)
    
    plot(plot)
  }
  dev.off()
  
}

###########################################################################
# calculate and plot missingness by cluster
###########################################################################
if (by_cluster == TRUE){
  #make proprotion missing table by cluster_id
  #start an empty table for results by cluster_id
  missing_table <- list()
  missing_table$cluster_id <- unique(lri_data$cluster_id)
  missing_table$prop_cough_missing <- NA
  missing_table$prop_db_missing <- NA
  missing_table$prop_chest_missing <- NA
  missing_table$prop_fever_missing <- NA
  missing_table <- as.data.frame(missing_table) %>%
    filter(!is.na(cluster_id))
  
  cluster_ids <- unique(missing_table$cluster_id)
  for (cluster in cluster_ids){
    print(cluster)
    cluster_data <- filter(lri_data, cluster_id == cluster)
    
    #cough missing: 1) filter for yes or missing DB 2) missing / (missing + answered)
    cough_data <- dplyr::filter(cluster_data, diff_breathing == 1 | diff_breathing == 9)
    missing <- dplyr::filter(cough_data, cough_missing == 1)
    answered <- dplyr::filter(cough_data, cough_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$cluster_id == cluster,]$prop_cough_missing <- m_prop
    
    #db missing: 1) dplyr::filter for yes or missing cough 2) missing / (missing + answered)
    db_data <- dplyr::filter(cluster_data, had_cough == 1 | had_cough == 9)
    missing <- dplyr::filter(db_data, db_missing == 1)
    answered <- dplyr::filter(db_data, db_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$cluster_id == cluster,]$prop_db_missing <- m_prop
    
    #chest missing: 1) dplyr::filter for cough and db answered or missing 2) missing / (missing + answered)
    case_data <- dplyr::filter(cluster_data, (had_cough == 1 | had_cough == 9) & (diff_breathing == 1 | diff_breathing == 9))
    missing <- dplyr::filter(case_data, chest_missing == 1)
    answered <- dplyr::filter(case_data, chest_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$cluster_id == cluster,]$prop_chest_missing <- m_prop
    
    #fever missing: 1) dplyr::filter for cough and db answered or missing 2) missing / (missing + answered)
    case_data <- dplyr::filter(cluster_data, (had_cough == 1 | had_cough == 9) & (diff_breathing == 1 | diff_breathing == 9))
    missing <- dplyr::filter(case_data, fever_missing == 1)
    answered <- dplyr::filter(case_data, fever_missing == 0)
    n_missing <- nrow(missing)
    n_answered <- nrow(answered)
    m_prop <- n_missing/(n_missing + n_answered)
    
    missing_table[missing_table$cluster_id == cluster,]$prop_fever_missing <- m_prop
  }
  
  #bind on cluster size, nid
  info <- dplyr::select(input_data, c('cluster_id','nid','N','survey_series','country','year')) %>%
    unique.data.frame()
  missing_table <- merge(missing_table, info, by = 'cluster_id')
  
  #create indicator for missingness greater than a certain cutoff
  # 1 if over limit, 0 if missingness below limit
  missing_table$cough_lim <- NA
  missing_table$db_lim <- NA
  missing_table$chest_lim <- NA
  missing_table$fever_lim <- NA

  missing_table <- mutate(missing_table, cough_lim = replace(cough_lim, prop_cough_missing > missing_lim, 1))
  missing_table <- mutate(missing_table, cough_lim = replace(cough_lim, prop_cough_missing <= missing_lim | is.nan(prop_cough_missing), 0))
  
  missing_table <- mutate(missing_table, db_lim = replace(db_lim, prop_db_missing > missing_lim, 1))
  missing_table <- mutate(missing_table, db_lim = replace(db_lim, prop_db_missing <= missing_lim | is.nan(prop_db_missing), 0))
  
  missing_table <- mutate(missing_table, chest_lim = replace(chest_lim, prop_chest_missing > missing_lim, 1))
  missing_table <- mutate(missing_table, chest_lim = replace(chest_lim, prop_chest_missing <= missing_lim | is.nan(prop_chest_missing), 0))
  
  missing_table <- mutate(missing_table, fever_lim = replace(fever_lim, prop_fever_missing > missing_lim, 1))
  missing_table <- mutate(missing_table, fever_lim = replace(fever_lim, prop_fever_missing <= missing_lim | is.nan(prop_fever_missing), 0))

  #save in case of future errors
  write.csv(missing_table, '<<<< FILEPATH REDACTED >>>>')
  
  #for each nid, what proportion of clusters in each indicator are over the threshold?
  nid_col <- dplyr::select(missing_table, c('nid','cough_lim','db_lim','chest_lim','fever_lim','survey_series','country','year')) %>%
    melt(measure.vars = c('cough_lim','db_lim','chest_lim','fever_lim'), variable.name = 'variable', value.name = 'over') %>%
    dcast(survey_series + nid + country + year ~ variable, value.var = 'over', fun.agg = function(x) mean(x))
  
  #bind on full country names to nid table
  locs <- get_location_metadata(location_set_id = 35) %>%
    dplyr::select(c('ihme_loc_id','location_name'))
  nid_col <- merge(nid_col, locs, by.x = 'country', by.y = 'ihme_loc_id')
  sample_sz <- dplyr::select(input_data, c('N','nid')) %>%
    dcast(formula = nid ~ ., sep = '_', value.var = 'N', fun.agg = function(x) sum(x)) %>%
    rename(c('.' = 'sum_N'))
  
  nid_col <- merge(nid_col, sample_sz, by = 'nid')
  nid_col$nid <- as.character(nid_col$nid)
  
  #save in case of future errors
  write.csv(nid_col, '<<<< FILEPATH REDACTED >>>>')
  
  ######################################################################
  # (2) Make plots
  ######################################################################
  #bar plot for each country-sx with x = nid, y = prop clusters missing greater than threshold, width of bar is N
  #open pdf
  pdf('<<<< FILEPATH REDACTED >>>>')

  #loop over countries and plot
  for (loc in unique(nid_col$country)){
    loc_data <- dplyr::filter(nid_col, country == loc)
    country_name <- unique(loc_data$location_name)
    print(country_name)
    sample_sum <- sum(loc_data$sum_N)
    loc_data$sample_prop <- loc_data$sum_N/sample_sum
    plot <- ggplot(data = loc_data, aes(x = nid, y = cough_lim)) + geom_col(aes(fill = survey_series) ) +
      ggtitle(paste0(country_name, ' | cough')) + ylab(paste0('proportion clusters with missingness > ', missing_lim)) + xlab('NID') + coord_cartesian(ylim = c(0,1))
    plot(plot)

    plot <- ggplot(data = loc_data, aes(x = nid, y = db_lim)) + geom_col(aes(fill = survey_series)) +
      ggtitle(paste0(country_name, ' | difficulty breathing')) + ylab(paste0('proportion clusters with missingness > ', missing_lim)) + xlab('NID') + coord_cartesian(ylim = c(0,1))
    plot(plot)

    plot <- ggplot(data = loc_data, aes(x = nid, y = chest_lim)) + geom_col(aes(fill = survey_series)) +
      ggtitle(paste0(country_name, ' | chest')) + ylab(paste0('proportion clusters with missingness > ', missing_lim)) + xlab('NID') + coord_cartesian(ylim = c(0,1))
    plot(plot)

    plot <- ggplot(data = loc_data, aes(x = nid, y = fever_lim)) + geom_col(aes(fill = survey_series)) +
      ggtitle(paste0(country_name, ' | fever')) + ylab(paste0('proportion clusters with missingness > ', missing_lim)) + xlab('NID') + coord_cartesian(ylim = c(0,1))
    plot(plot)
  }
  dev.off()
}

melt <- melt(missing_table, id.vars = c('cluster_id', 'V1', 'nid','N','survey_series','country','year'), 
             measure.vars = c('prop_cough_missing','prop_db_missing','prop_chest_missing','prop_fever_missing'))
locs <- get_location_metadata(location_set_id = 35) %>%
  dplyr::select(c('ihme_loc_id','location_name'))
melt <- merge(melt, locs, by.x = 'country', by.y = 'ihme_loc_id')
nid_list <- unique(melt$nid)


pdf(width = 20, '<<<< FILEPATH REDACTED >>>>')
for (survey in nid_list){
  melt_filter <- filter(melt, nid == survey)
  country_name <- unique(melt_filter$location_name)
  year <- unique(melt_filter$year)
  survey_series <- unique(melt_filter$survey_series)
  print(paste(country_name, survey_series, year))
  plot <- ggplot(data = melt_filter, aes(x = cluster_id, y = variable)) + geom_tile(aes(fill = value)) + scale_fill_gradient(low = "white", high = "steelblue") + ggtitle(paste(country_name, survey_series, year))
  plot(plot)
}
dev.off()
