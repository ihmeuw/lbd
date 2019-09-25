######################################################################################
#### Draft of script to syncronize LRI data processing between GBD and Geospatial ####
####                  Scott Swartz | 6/29/2017 | sswartz@uw.edu                   ####
####                 ADJUSTED BACK TO GBD 2016 PROCESS Aug. 2017                  ####
####             missing month scalars incorporated to Geosp. 8/24/17             ####
####                      also incorporated to GBD 9/28/17                        ####
####                    included Kish approximation 11/6/17                       ####
######################################################################################

rm(list=ls())

sing_image = TRUE #whether running on RStudio singularity image

### PREP ###
#Set filepaths
root <- '<<<< FILEPATH REDACTED >>>>' #filepath root
currentdata <- '<<<< FILEPATH REDACTED >>>>' #filepath to most recent geomatched dataset

library(boot)
#Geospatial setup for the cluster
#Set number of cores
numcores = 2
#Set repo
repo <- '<<<< FILEPATH REDACTED >>>>' #path to MBG central code
setwd(repo)
package_lib <- '<<<< FILEPATH REDACTED >>>>' #path to package library
.libPaths(package_lib)
#Source MBG central code
source('collapse_functions.R')
source('covariate_functions.R')
source('holdout_functions.R')
source('mbg_functions.R')
source('misc_functions.R')
source('polygon_functions.R')
source('post_estimation_functions.R')
source('prep_functions.R')
source('seegMBG_transform_functions.R')
source('shiny_functions.R')
source('shapefile_functions.R')

#load packages
package_list <- c('survey', 'foreign', 'rgeos', 'data.table','raster','rgdal',
                  'plyr', 'foreach', 'doParallel', 'ggplot2', 'lme4', 'dplyr',
                  'matrixStats', 'gridExtra', 'feather')
library(pacman)
for(package in package_list) {
  p_load(package, character.only=TRUE)
}



#Read in supplementary information
locs <- '<<<< FILEPATH REDACTED >>>>' #ihme location metadata
locs <- locs[,c("location_name","location_id","super_region_name","super_region_id",
                "region_name","region_id","parent_id", "ihme_loc_id")]
months <- '<<<< FILEPATH REDACTED >>>>' #seasonality scalars
mm_scalar <- '<<<< FILEPATH REDACTED >>>>' #seasonality scalars & info from surveys without month microdata
dm.coeffs <- '<<<< FILEPATH REDACTED >>>>' #DisMod crosswalk coefficients
duration <- '<<<< FILEPATH REDACTED >>>>' #LRI duration draws

#Read in dataset
lri <- as.data.table(read_feather(currentdata))
setnames(lri, old = c('iso3', 'year_start'), new = c('ihme_loc_id', 'start_year'))


#Merge locs & seasonality scalars onto dataset
lri <- join(lri, locs, by = "ihme_loc_id")
scalar <- months[,c("scalar","month","region_name")]
setnames(lri, old = 'int_month', new = 'month')
lri <- join(lri, scalar, by=c("region_name","month"))


### CLEAN UP DATA ###
#Change NA's to 9's in symptom columns
lri$had_fever[is.na(lri$had_fever)] = 9
lri$had_cough[is.na(lri$had_cough)] = 9
lri$diff_breathing[is.na(lri$diff_breathing)] = 9
lri$chest_symptoms[is.na(lri$chest_symptoms)] = 9

#Drop any NIDs completely missing difficulty breathing responses
lri$tabulate <- ave(lri$diff_breathing, lri$nid, FUN= function(x) min(x))
lri <- subset(lri, tabulate != 9)

#Drop rows w/o sex, round ages, drop over 5s
setnames(lri, old = 'sex_id', new = 'child_sex')
lri$child_sex[is.na(lri$child_sex)] = 3
lri$age_year <- floor(lri$age_year)
lri <- subset(lri, age_year <= 4)

#Reassign pweight = 1 if missing or 0
#Include pwt_adj indicator if reassigned to 1
lri <- lri[, pwt_adj := 0]
lri <- lri[is.na(pweight) | pweight == 0, pwt_adj := 1]
lri <- lri[is.na(pweight), pweight := 1]
lri <- lri[pweight == 0, pweight := 1]

#Copy data for Geospatial use later
lri_data <- copy(lri)

lri <- lri[, start_n := .N, by = c("nid", "ihme_loc_id", "start_year")]


#Drop if missing PSU
lri <- subset(lri, !is.na(psu))


lri <- lri[,psudrop_n := .N, by = c("nid", "ihme_loc_id", "start_year")]


#Create a dummy NID for each GBD subnational
lri$urban_name <- ifelse(lri$urban==1,"Urban","Rural")
lri$nid.new <- ifelse(lri$ihme_loc_id=="ZAF", paste0(lri$location_name,"_",lri$nid),
                      ifelse(lri$ihme_loc_id=="IND", paste0(lri$location_name,"_",lri$urban_name,"_",lri$nid),lri$nid))
lri$subname <- ifelse(lri$ihme_loc_id=="ZAF", lri$location_name, ifelse(lri$ihme_loc_id=="IND", paste0(lri$location_name,", ",lri$urban_name), "none"))
nid.new <- unique(lri$nid.new)

### COLLAPSE TO COUNTRY-YEAR-AGE-SEX ###
#Data frame for outputs
df <- data.frame()

#Loop over surveys to 1) Assign case definitions, 2) Incorporate seasonality, 3) Apply surveydesign & collapse to country-year-age-sex
for(n in nid.new) {
  message(n)
  temp <- subset(lri, nid.new == n)

  #Set indicator = 1 if survey has chest/fever, set 0 if not
  exist_chest <- ifelse(min(temp$chest_symptoms, na.rm=T) == 0, 1, 0)
  exist_fever <- ifelse(min(temp$had_fever, na.rm=T) == 0, 1, 0)

  #Pull in recall period (if bsswlank, assume it's 2 weeks)
  temp$recall_period <- ifelse(is.na(temp$recall_period_weeks), 2, temp$recall_period_weeks)
  recall_period <- max(temp$recall_period, na.rm=T)

  #Assign case definitions
  #We want prevalence w/o fever for crosswalk
  temp$good_no_fever <- 0
  temp$poor_no_fever <- 0
  if(exist_fever==1 & exist_chest==1) {
    temp$overall_lri = ifelse(temp$had_fever==1 & temp$chest_symptoms==1, 1, 0)
    temp$good_no_fever = ifelse(temp$chest_symptoms==1,1,0)
  } else if (exist_fever==1 & exist_chest==0) {
    temp$overall_lri = ifelse(temp$had_fever==1 & temp$diff_breathing==1, 1, 0)
    temp$poor_no_fever = ifelse(temp$diff_breathing==1,1,0)
  } else if (exist_fever==0 & exist_chest==1) {
    temp$overall_lri = ifelse(temp$chest_symptoms==1, 1, 0)
    temp$good_no_fever = 0
  } else if (exist_fever==0 & exist_chest==0) {
    temp$overall_lri = ifelse(temp$diff_breathing==1, 1, 0)
    temp$poor_no_fever = 0
  }

  #Add missing month seasonality scalars
  if(n %in% mm_scalar$nid) {
    mm_scalar = as.data.table(mm_scalar)
    mm_short = mm_scalar[, c("nid", "country", "avg_scalar")]
    mm_short$scalar = as.numeric(NA) #add on NA column so it only merges on rows missing scalars
    #merge in missing month scalars
    temp = merge(temp, mm_short, by.x = c("nid", "ihme_loc_id", "scalar"), by.y = c("nid", "country", "scalar"), all.x = TRUE)
    temp$scalar = ifelse(!is.na(temp$avg_scalar), temp$avg_scalar, temp$scalar)
    temp = within(temp, rm(avg_scalar))
  }
  #Apply seasonality scalars
  scalar.dummy <- max(temp$scalar)
  if(is.na(scalar.dummy)) {
    temp$scalar_lri <- temp$overall_lri
  } else {
    temp$scalar_lri <- temp$overall_lri*temp$scalar
    temp$good_no_fever <- temp$good_no_fever*temp$scalar
    temp$poor_no_fever <- temp$poor_no_fever*temp$scalar
  }

  #Apply survey design & collapse
  dclus <- svydesign(id=~psu, weights=~pweight, data=temp)
  prev <- svyby(~scalar_lri, ~child_sex + age_year, dclus, svymean, na.rm=T)
  prev$base_lri <- svyby(~overall_lri, ~child_sex + age_year, dclus, svymean, na.rm=T)$overall_lri
  prev$good_no_fever <- svyby(~good_no_fever, ~child_sex + age_year, dclus, svymean, na.rm=T)$good_no_fever
  prev$poor_no_fever <- svyby(~poor_no_fever, ~child_sex + age_year, dclus, svymean, na.rm=T)$poor_no_fever
  prev$sample_size <- svyby(~scalar_lri, ~child_sex + age_year, dclus, unwtd.count, na.rm=T)$count
  prev$cases <- prev$scalar_lri * prev$sample_size
  prev$ihme_loc_id <- unique(temp$ihme_loc_id)
  prev$location <- unique(temp$subname)
  prev$start_year <- unique(temp$start_year)
  prev$end_year <- unique(temp$end_year)[1]
  prev$cv_diag_valid_good <- exist_chest
  prev$cv_had_fever <- exist_fever
  prev$nid <- unique(temp$nid)
  prev$sex <- ifelse(prev$child_sex==2,"Female","Male")
  prev$age_start <- prev$age_year
  prev$age_end <- prev$age_year + 1
  prev$recall_period <- recall_period
  prev$survey <- unique(temp$survey_name)
  prev$notes <- paste0("LRI prevalence adjusted for seasonality. The original value was ", round(prev$base_lri,4))

  message(nrow(prev))
  df <- rbind.data.frame(df, prev)
}


uniqlri = unique(lri[, c("nid", "ihme_loc_id", "start_year", "start_n", "psudrop_n")])
dt_df = as.data.table(df)
dt_df = dt_df[, collapsed_n := sum(sample_size), by = c("nid", "ihme_loc_id", "start_year")]
dt_df = merge(dt_df, uniqlri, by = c("nid", "ihme_loc_id", "start_year"))


#Export for reference
write.csv(df, '<<<< FILEPATH REDACTED >>>>', row.names = F)

#Replace missing values so they can be used
df$missing <- ifelse(df$scalar_lri == 0,1,0)
df$cv_had_fever <- ifelse(df$missing == 1, ifelse(df$good_no_fever + df$poor_no_fever != 0, 0, df$cv_had_fever), df$cv_had_fever)
df$scalar_lri <- ifelse(df$missing == 1, ifelse(df$good_no_fever!=0, df$good_no_fever, df$poor_no_fever), df$scalar_lri)

prevs <- df

### CONVERT FROM PERIOD --> POINT PREVALENCE, GENERATE DRAWS ###
#Use (period*duration)/(duration + recall - 1)
prevs <- prevs[,c("scalar_lri","recall_period")]

prevs$recall_period <- as.numeric(prevs$recall_period)

prevs$scalar_lri[prevs$scalar_lri==0] <- 0.0001

#Convert to logit space
prevs$ln_lri <- logit(prevs$scalar_lri)
prevs$duration <- duration$mean
prevs$draw <- inv.logit(prevs$ln_lri)
prevs$point <- (prevs$draw * prevs$duration) / (prevs$duration + (prevs$recall_period*7) - 1)
#Append to original data
df$period_mean <- prevs$draw
df$mean <- prevs$point

#Convert to cases
df$cases <- df$sample_size * df$mean
df$notes <- paste("Converted to point prevalence using duration draws. Period prevalence was", round(df$scalar_lri,4),".",df$notes)

#Add duration w/ uncertainty
df$duration <- rowMeans(duration[,7:1006])
df$duration_lower <- as.numeric(quantile(duration[,7:1006], 0.025))
df$duration_upper <- as.numeric(quantile(duration[,7:1006], 0.975))

dt_df2 <- as.data.table(df)
dt_df2 <- dt_df2[, point_n := sum(sample_size), by = c("nid", "ihme_loc_id", "start_year")]
dt_df2 <- unique(dt_df2[, c("nid", "ihme_loc_id", "start_year", "point_n")])
dt_df <- merge(dt_df, dt_df2, by = c("nid", "ihme_loc_id", "start_year"))


#Export for reference
write.csv(df, '<<<< FILEPATH REDACTED >>>>', row.names = F)

### GENERATE CASE DEFINITION CROSSWALK COEFFICIENTS ###
#Keep record of previous values
df$ln_mean <- log(df$mean)

df$cv_no_fever <- (1-df$cv_had_fever)

#Tack on geography information
locs <- locs[, c("ihme_loc_id", "region_name", "super_region_name", "location_id")]
df$ihme_loc_id <- ifelse(df$ihme_loc_id == "KOSOVO", "SRB", df$ihme_loc_id)
df$ihme_loc_id <- ifelse(df$ihme_loc_id == "KEN_44798", "KEN", df$ihme_loc_id)
df <- join(df, locs, by="ihme_loc_id")


##Regress cv_no_fever WITH chest_symptoms##
#ihme_loc_id has a lot of variation, Chris uses region_name for random effect because it is more stable
mod <- lmer(ln_mean ~ cv_no_fever + (1|age_end) + (1|region_name), data=subset(df, cv_diag_valid_good==1))
cv_no_fever_good = fixef(mod)[2]
df$mean <- with(df, ifelse(cv_diag_valid_good==1, ifelse(cv_no_fever==1, exp(ln_mean - cv_no_fever_good),mean), mean))

mod2 <- lmer(ln_mean ~ cv_no_fever + factor(age_end) + (1|region_name), data=subset(df, cv_diag_valid_good==1))

#Bootstrap for variance of adjustment
val <- fixef(mod)[2]
var <- vcov(mod)[2,2]

#Take exponentiated draws for beta_coefficient and variance from those
variance_draws <- exp(rnorm(n=1000, mean=val, sd=sqrt(var)))
variance <- rowSds(matrix(variance_draws, nrow=1))^2

# #Combine variances using equation: mean(A)^2 * variance(B) + mean(B)^2 * variance(A) + variance(A)*variance(B)
temp_variance <- df$mean^2*variance + (exp(val))^2*df$standard_error^2 + df$standard_error^2*variance

##Regress cv_no_fever WITHOUT chest_symptoms##
mod <- lmer(ln_mean ~ cv_no_fever + (1|age_end) + (1|region_name), data=subset(df, cv_diag_valid_good==0))
cv_no_fever_poor <- fixef(mod)[2]
df$mean <- with(df, ifelse(cv_diag_valid_good==0, ifelse(cv_no_fever==1, exp(ln_mean - fixef(mod)[2]),mean), mean))

#Bootstrap for variance of adjustment
val <- fixef(mod)[2]
var <- vcov(mod)[2,2]

#Take exponentiated draws for beta_coefficient and variance from those
variance_draws <- exp(rnorm(n=1000, mean=val, sd=sqrt(var)))
variance <- rowSds(matrix(variance_draws, nrow=1))^2

#Combine variances using equation mean(A)^2 * variance(B) + mean(B)^2 * variance(A) + variance(A)*variance(B)
temp_variance <- df$mean^2*variance + (exp(val))^2*df$standard_error^2 + df$standard_error^2*variance

#Clean up standard errors
if(dropping_zeros == TRUE) {
  df$standard_error <- with(df, ifelse(cv_diag_valid_good==0, ifelse(cv_no_fever==1, sqrt(temp_variance),standard_error), standard_error))
}

#Reassign case count using new mean
df$cases <- df$sample_size * df$mean

dt_df3 <- as.data.table(df)
dt_df3 <- dt_df3[, xwalk_n := sum(sample_size), by = c("nid", "ihme_loc_id", "start_year")]
dt_df3 <- unique(dt_df3[, c("nid", "ihme_loc_id", "start_year", "xwalk_n")])
dt_df <- merge(dt_df, dt_df3, by = c("nid", "ihme_loc_id", "start_year"))


#Export
write.csv(df, '<<<< FILEPATH REDACTED >>>>', row.names = F)

#Final crosswalk coefficients (was previously log, now exponentiated):
xwalks = c(cv_no_fever_good, cv_no_fever_poor)



##################################################################################################################
###                                   GEOSPATIAL COLLAPSE (MUST RUN ON CLUSTER)                               ####
##################################################################################################################

## INITIAL CLEANUP ##
setnames(lri_data, old = c('ihme_loc_id'), new = c('country'))
lri_data$country <- substr(lri_data$country, 1, 3)
lri_data <- lri_data[, N := 1]
lri_data$country <- as.character(lri_data$country)
#fix select multi-year surveys so start_year is the actual interview year
lri_data <- lri_data[nid %in% multi_srv, year := int_year]
lri_data <- lri_data[nid %in% multi_srv, start_year := int_year]


## ASSIGN CASE DEFINITIONS ##
lri_data <- lri_data[, cv_fever := numeric()]
lri_data <- lri_data[, cv_good := numeric()]
lri_data <- lri_data[, case := 0]
nidlist <- unique(lri_data$nid)
for(n in nidlist) {
  temp <- lri_data[nid == n,]
  exist_chest <- ifelse(min(temp$chest_symptoms, na.rm=T) == 0, 1, 0)
  exist_fever <- ifelse(min(temp$had_fever, na.rm=T) == 0, 1, 0)
  if(exist_fever==1 & exist_chest==1) {
    lri_data <- lri_data[nid == n, cv_fever := 1]
    lri_data <- lri_data[nid == n, cv_good := 1]
    lri_data <- lri_data[nid == n & had_fever == 1 & chest_symptoms == 1, case := 1]
  } else if (exist_fever==1 & exist_chest==0) {
    lri_data <- lri_data[nid == n, cv_fever := 1]
    lri_data <- lri_data[nid == n, cv_good := 0]
    lri_data <- lri_data[nid == n & had_fever == 1 & diff_breathing == 1, case := 1]
  } else if (exist_fever==0 & exist_chest==1) {
    lri_data <- lri_data[nid == n, cv_fever := 0]
    lri_data <- lri_data[nid == n, cv_good := 1]
    lri_data <- lri_data[nid == n & chest_symptoms == 1, case := 1]
  } else if (exist_fever==0 & exist_chest==0) {
    lri_data <- lri_data[nid == n, cv_fever := 0]
    lri_data <- lri_data[nid == n, cv_good := 0]
    lri_data <- lri_data[nid == n & diff_breathing == 1, case := 1]
  }
}


## CROSSWALK CASE DEFINITIONS & APPLY SEASONALITY SCALARS ##
#Apply 'fever - good' and 'fever - poor' crosswalks from GBD data (only to LRI case rows)
lri_data <- lri_data[cv_fever == 0 & cv_good == 1 & case != 0, case := case/exp(xwalks[1])]
lri_data <- lri_data[cv_fever == 0 & cv_good == 0 & case != 0, case := case/exp(xwalks[2])]

#Apply DisMod 'self-report' crosswalk to all values
lri_data <- lri_data[, case := case/dm.coeffs$`self-report`]

#Apply DisMod 'poor validity' crosswalk to those without chest symptoms
lri_data <- lri_data[cv_good == 0, case := case/dm.coeffs$poor_validity]

##Apply seasonality scalars
#prep missing month scalars
mm_scalar = as.data.table(mm_scalar)
mm_short = mm_scalar[, c("nid", "country", "avg_scalar")]
mm_short = mm_short[, scalar := as.numeric(NA)] #add on NA column so it only merges on rows missing scalars
#merge in missing month scalars
lri_data = lri_data[, -c(11, 15)]
lri_data = merge(lri_data, mm_short, by = c("nid", "country", "scalar"), all.x = TRUE)
lri_data = lri_data[!is.na(avg_scalar), scalar := avg_scalar]
lri_data = lri_data[, -c("avg_scalar")]
#after pulling in missing month scalars, now drop rows still missing scalars
lri_data <- lri_data[!is.na(scalar),]
#apply seasonality scalars
lri_data <- lri_data[, case := case*scalar]


## PULL OUT RECALL PERIOD INFO FOR LATER ##
recall <- lri_data[, c("nid", "country", "survey_series", "start_year", "recall_period_weeks")]
recall <- unique(recall)
write.csv(recall, '<<<< FILEPATH REDACTED >>>>')

## SPLIT POINTS & POLYGONS ##
lri_data_point <- lri_data[point == 1,]
lri_data_point <- lri_data_point[!is.na(lat)]
lri_data_point <- lri_data_point[!is.na(long)]
lri_data_poly <- lri_data[point != 1,]

## COLLAPSE TO POLYGONS ##
#Collapse each NID separately to the most granular shapefiles
collapse_each_nid <- function(this_nid) {

  message(paste0('Collapsing NID: ', this_nid))
  test_poly_data <- lri_data_poly[nid==this_nid,]
  test_poly_data$strata <- 0
  names(test_poly_data)[names(test_poly_data)=='cluster_number'] <- 'psu'
  names(test_poly_data)[names(test_poly_data)=='weight'] <- 'pweight'

  # Check for missings
  if(length(test_poly_data$pweight[is.na(test_poly_data$pweight)])>0) {
    message(paste0(length(test_poly_data$pweight[is.na(test_poly_data$pweight)]), ' / ', length(test_poly_data$pweight), ' are missing pweight'))
    return(NULL)
  } else {
    collapse_polys <- function(x) {
      by_vars <- c('start_year', 'country', 'location_code', 'shapefile', 'survey_series', 'nid')
      #calculate raw, unweighted sample sizes & # of cases
      test_poly_data <- test_poly_data[, N_orig := .N, by = by_vars]
      test_poly_data <- test_poly_data[, cases_orig := sum(case), by = by_vars]
      test_poly_data <- test_poly_data[, origmean := cases_orig/N_orig]
      #generate weighted Ns with Kish
      test_poly_data <- test_poly_data[, weights2 := (pweight^2)]
      test_poly_data <- test_poly_data[, n_eff := (sum(pweight)^2)/(sum(weights2)), by = by_vars]
      #generate weighted means with Kish
      test_poly_data <- test_poly_data[, wtmean := weighted.mean(case, pweight), by = by_vars]
      test_poly_data <- test_poly_data[, sum_of_sample_weights := sum(pweight), by = by_vars]
      collapsed <- test_poly_data[, c('start_year', 'country', 'location_code', 'shapefile', 'survey_series','nid',
                                      'wtmean', 'n_eff', 'origmean', 'N_orig','sum_of_sample_weights','point')]
      collapsed <- unique(collapsed) #shorten to one row per polygon
      names(collapsed)[names(collapsed)=='wtmean'] <- x
      names(collapsed)[names(collapsed)=='n_eff'] <- 'N'
      collapsed <- as.data.frame(collapsed)
      return(collapsed)
    }
    polys <- c('case')
    polys <- lapply(polys, collapse_polys)
    merged_polys <- Reduce(function(...) merge(..., all=T), polys)
    return(merged_polys)
  }
}

# Subset out problematic NIDs: 1) UGA LSMS_ISA 2013 throwing an error in collapse_by, 2) PSE Health Survey 1996 currently has no PSU
# New problematic NIDS: 1) BGD HH Svy, unknown issue; 2/3/4) NGA HH survey and 2x PSE surveys missing NIDs, 5) NIC RHS, incomplete PSUs
# Actually point data...dropping SUSENAS 2011 SEPT 85262
lri_data_poly <- lri_data_poly[!is.na(lri_data_poly$case),]
poly_nids <- mclapply(poly_nids, collapse_each_nid, mc.cores=numcores)
poly_nids <- poly_nids[sapply(poly_nids, function(x) class(x) != "try-error")]


##Collapse point data (w/ Kish approximation)
lri_data_point <- lri_data_point[, c('nid', 'start_year', 'country', 'geospatial_id', 'lat', 'long', 'survey_series', 'case', 'pweight', 'psu'), with=F]
by_vars = c('survey_series', 'start_year', 'country', 'nid', 'geospatial_id', 'psu')
lri_data_point <- lri_data_point[!is.na(pweight),] #drop if no pweights
#calculate weighted N with Kish
all_point_data <- lri_data_point[, weights2 := (pweight^2)]
all_point_data <- all_point_data[, n_eff := (sum(pweight)^2)/(sum(weights2)), by = by_vars]
#calculate weighted mean with Kish
all_point_data <- all_point_data[, wtmean := weighted.mean(case, pweight), by = by_vars]
#calculate raw, unweighted sample sizes & # of cases
all_point_data <- all_point_data[, N_orig := .N, by = by_vars]
all_point_data <- all_point_data[, cases_orig := sum(case), by = by_vars]
all_point_data <- all_point_data[, origmean := cases_orig/N_orig]
all_point_data <- all_point_data[, sum_of_sample_weights := sum(pweight), by = by_vars]
#drop too-unique columns
all_point_data <- all_point_data[, -c('geospatial_id', 'psu', 'case', 'pweight', 'int_month', 'int_year', 'weights2', 'cases_orig')]
all_point_data <- unique(all_point_data) #unique to shorten to one row per psu

names(all_point_data)[names(all_point_data)=='wtmean'] <- 'case'
names(all_point_data)[names(all_point_data)=='n_eff'] <- 'N'
all_point_data$point <- 1

#clean up poly data
all_poly_data <- rbindlist(poly_nids)
all_poly_data$point <- 0

# Append polygon data to point data
collapsed <- rbind(all_poly_data, all_point_data, fill=TRUE)
setnames(collapsed, old = c('case'), new = c('has_lri'))

# Export dataset as collapse checkpoint (has preserved original mean & SS without Kish)
write.csv(collapsed, '<<<< FILEPATH REDACTED >>>>')

#launch resample polygons
source('resample_parent.R')

#read in all resampled CSVs
resampled <- list.files('<<<< FILEPATH REDACTED >>>>', full.names = T, pattern = ".csv", recursive = F)
cluster <- TRUE
if(cluster == TRUE) {
  message("Make second cluster")
  cl <- makeCluster(numcores)
  clusterEvalQ(cl, .libPaths('<<<< FILEPATH REDACTED >>>>'))
  message("Register cluster")
  registerDoParallel(cl)
  message("Start foreach")
  #Read in each .dta file in parallel - returns a list of data frames
  top <- foreach(i=1:length(resampled), .packages = c('haven')) %dopar% {
    dta <- read.csv(resampled[i])
    return(dta)
  }
  message("Foreach finished")
  message("Closing cluster")
  stopCluster(cl)
} else if(cluster == FALSE) {
  top <- foreach(i=1:length(resampled)) %do% {
    message(paste0("Reading in: ", resampled[i]))
    dta <- read.csv(resampled[i])
    return(dta)
  }
}

#rbind all resamples together
all_poly_data <- rbindlist(top, fill=T, use.names=T)
rm(top)
names(all_poly_data)[names(all_poly_data)=="lat"] <- "latitude"
names(all_poly_data)[names(all_poly_data)=="long"] <- "longitude"
all_poly_data <- all_poly_data[,c('longitude', 'latitude', 'weight', 'shapefile', 'start_year', 'country',
                                  'source', 'location_code', 'svy_id', 'lriprevalence', 'N', 'N_orig',
                                  'point','sum_of_sample_weights')]

#set weight of all not resampled points to 1
all_point_data <- collapsed[point == 1, ]
all_point_data <- all_point_data[, weight := 1]
all_point_data <- all_point_data[, pseudocluster := FALSE]
names(all_point_data)[names(all_point_data)=="lat"] <- "latitude"
names(all_point_data)[names(all_point_data)=="long"] <- "longitude"
#assign count values
all_point_data <- all_point_data[, has_lri_count := lriprevalence * N]
all_point_data <- all_point_data[, names(all_poly_data), with=FALSE]
#all_poly_data <- all_poly_data[, has_lri := has_lri * N]
all_collapsed <- rbind(all_point_data, all_poly_data)


## FORMAT, DURATION, MAL-ED & SAVE ##
#Format and save
setnames(all_collapsed, c('start_year', 'svy_id', 'source'),  c('year', 'nid', 'survey_series'))
all_collapsed <- all_collapsed[, c('nid', 'survey_series', 'year','country','latitude','longitude','N','lriprevalence','weight','sum_of_sample_weights','point'), with = FALSE]
all_collapsed <- all_collapsed[, cluster_id := .GRP, by = c('year','latitude','longitude','country','survey_series')]
names(all_collapsed)[names(all_collapsed)=='survey_series'] <- 'source'
all_collapsed <- all_collapsed[order(cluster_id,latitude,longitude,lriprevalence,year,source)]

#subset to relevant years
all_collapsed <- subset(all_collapsed, year >= 2000)

#keep only acceptable points
all_collapsed <- all_collapsed[!is.na(latitude)]
all_collapsed <- all_collapsed[!is.na(longitude)]
all_collapsed <- all_collapsed[latitude>=-90 & latitude<=90]
all_collapsed <- all_collapsed[longitude>=-180 & longitude<=180]

#Merge recall periods, written earlier, by survey series, country, start year
recall <- read.csv('<<<< FILEPATH REDACTED >>>>')
recall_un = unique(recall[, c("country", "survey_series", "start_year", "recall_period_weeks")])
all_collapsed <- merge(all_collapsed, recall_un, by.x = c("source", "year", "country"), by.y = c("survey_series", "start_year", "country"))
all_collapsed$recall_period_weeks <- as.numeric(all_collapsed$recall_period_weeks)
all_collapsed <- all_collapsed[, recall_days := recall_period_weeks*7]
all_collapsed <- all_collapsed[, -c("recall_period_weeks")]

#Convert to point prevalence [Period Prev * duration / (duration + recall - 1)]
dur_days = duration$mean
all_collapsed <- all_collapsed[, rate := lriprevalence * dur_days / (dur_days + recall_days - 1)]
all_collapsed <- all_collapsed[, lriprevalence := rate * N] #convert back to counts

#In clusters where LRI > N (due to tiny samples and every child having LRI), cap at N
all_collapsed <- all_collapsed[lriprevalence > N, lriprevalence := N]

write.csv(all_collapsed, file = '<<<< FILEPATH REDACTED >>>>')
