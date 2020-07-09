#########################################################################
###### RHF definition crosswalk              
###### This code crosswalks survey data to a standard       
###### RHF definition ('recommended home fluids' or 'acceptable fluids')
#########################################################################

###### Load required packages
definition_crosswalk_setup <- function() {
  if (!require('pacman')) install.packages('pacman')
  pacman::p_load(MASS, rio, ggplot2, data.table, plyr, dplyr, gtools, splines, lme4)

  ### Install googlesheets, if needed, and then load package; ditto for ggalt, geepack, geeM, MESS
  library(googlesheets)
}

###### Get tracking data from Google Sheet. Note: The function currently assumes that a tracking sheet ('sheet_name')
###### with a specific structure is used to identify nids requiring adjustment.
get_tracking_sheet_definition <- function(sheet_name) {
  ### Activate google sheet token
  load(paste0(indicator_repo, '<<<< FILEPATH REDACTED >>>>'))
  gs_auth(token = ttt)
  
  ### Retrieve data from tracking sheet
  tr <- gs_title('Diarrhea ORZ Tracking Sheet') ### Grab information about tracking sheet, by name
  tracking <- as.data.table(gs_read(ss=tr, ws='Tracking', skip=1)) ### Extract data from tracking sheet; skip the first line (unnecessary headers)
  tracking <- tracking[, c('nid', 'ORT Extraction Status', 'Sugar and salt solutions', 
                           'Other home fluids', 'Other liquid foods', 'Recommended or acceptable home fluids')]
  colnames(tracking) <- c('nid', 'extraction_status', 'sss', 'ofluid', 'lfood', 'rahf')
  
  return(tracking)
}

#### Adjust non-standard definitions
fit_definition_crosswalk <- function(data, use_country_fixed_effects=TRUE, use_temporal_effect=TRUE) {
  if (!use_country_fixed_effects & !use_temporal_effect) { # Apply a single model to all locations and times
    m <- glm(mean ~ 1 + as.factor(def), weights=N, data=data)
    coefs <- coef(m)[2]
  } else if (use_country_fixed_effects & !use_temporal_effect) { # Use country fixed effects but no temporal effect
    m <- glm(mean ~ 1 + as.factor(def) + as.factor(country_new), weights=N, data=data)
    coefs <- coef(m)[2]
  } else if (use_country_fixed_effects & use_temporal_effect) { # Use country fixed effects and a temporal effect
    m <- glm(mean ~ 1 + as.factor(def) + as.factor(country_new) + ns(year, df=3), weights=N, data=data)
    coefs <- coef(m)[2]
  } else if (!use_country_fixed_effects & use_temporal_effect) { # Use temporal effect but no country fixed effects
    m <- glm(mean ~ 1 + as.factor(def) + ns(year, df=3), weights=N, data=data)
  } 
  return(m)
}

###############################################################################################
###### Perform definition crosswalk

#### Initial setup
definition_crosswalk_setup()
tracking <- get_tracking_sheet_definition(sheet_name='Diarrhea ORZ Tracking Sheet')

### Remove excluded surveys
tracking <- tracking[extraction_status == 'Completed' & !is.na(rahf)]
tracking[, extraction_status := NULL]

#### Create definition categories
tracking[, def := lapply(.SD, paste0, ofluid, lfood, rahf),
         .SDcols = 'sss']

### Merge with data
dat <- merge(coverage_data, tracking, by = 'nid', all.x = T, allow.cartesian = T)
dat[, country_new := substr(country, 1, 3)]

### Perform some data checks
if (nrow(dat[is.na(mean)]) > 0) {
  message('\n', nrow(dat[is.na(mean)]), ' rows in data have NA mean values. They are being removed now.')
  message('This will affect nid(s): ', paste0(unique(dat[is.na(mean), nid]), sep = ' '), '\n')
  dat <- dat[!is.na(mean)]
}
if (nrow(dat[is.na(def)]) > 0) {
  message('\n', nrow(dat[is.na(def)]), ' rows in data are missing definition groups. They are being removed now.')
  message('This will affect nid(s): ', paste0(unique(dat[is.na(def), nid]), sep = ' '), '\n')
  dat <- dat[!is.na(def)]
}

### Separate data to adjust and data not to adjust
message('\n', nrow(dat[def == '0000']), ' rows in data do not ask RHF questions according to the tracking sheet. They will not be crosswalked.\n')
dat_unadjust <- dat[def == '0000']
dat <- dat[def != '0000']

### Fit crosswalk model and grab coefficients
def_model <- fit_definition_crosswalk(dat, use_country_fixed_effects = TRUE, use_temporal_effect = TRUE)
c <- data.table(def = levels(as.factor(dat$def)), coef = c(0, coef(def_model)[2:15]))
dat <- merge(dat, c, by = 'def', allow.cartesian = T)

### Apply adjustments
dat[, mean_mod := inv.logit(logit(mean) - coef)]

### Clean adjusted data
coverage_data <- copy(dat)
coverage_data[!is.na(mean_mod), mean := mean_mod]
coverage_data[, c('mean_mod', 'country_new', 'def', 'coef', 'sss', 'ofluid', 'lfood', 'rahf') := NULL]

### Bind with unajusted data
dat_unadjust[, c('country_new', 'def', 'sss', 'ofluid', 'lfood', 'rahf') := NULL]
coverage_data <- rbind(coverage_data, dat_unadjust)

##### END definition crosswalk
###############################################################################################