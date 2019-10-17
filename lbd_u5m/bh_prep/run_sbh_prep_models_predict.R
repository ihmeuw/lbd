## #######################################################
## PURPOSE: Apply SBH indirect prediction for a survey
## INPUTS: model object, cs object
## #######################################################

## SGE LOG HEADER
message(' ... ... ... ... ... ... ... ... ... ... ')
message('Run notes:')
print(R.Version())
print(Sys.info())
message(' ... ... ... ... ... ... ... ... ... ... ')



## #######################################################
## SETUP 
message('\nSETUP ... ')

# Load in args from qsub call
run_date    <- as.character(commandArgs()[4])
snid        <- as.numeric(commandArgs()[5])
agg_to_nid  <- as.logical(commandArgs()[6])

user       <- Sys.info()['user']

# print out run info for the log
message(sprintf('\nUSER: %s',user))
message(sprintf('\nNID: %s',snid))
message(sprintf('\nRUN_DATE: %s',run_date))
message(sprintf('\nagg_to_nid: %s',agg_to_nid))


# Get output directory
root   <- '<<<< FILEPATH REDACTED >>>>'

core_repo <- sprintf("<<<< FILEPATH REDACTED >>>>", user)
source(paste0(core_repo, '<<<< FILEPATH REDACTED >>>>'))

# Load packages
packlist <- c('data.table','parallel','seegMBG','survival','magrittr','scales','lme4','Biograph','ggplot2','pryr')
load_R_packages(packlist)

message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))


# load custom functions
setwd(sprintf('<<<< FILEPATH REDACTED >>>>', user))
source('sbh_utils.R')


# get best run date of fitted model from the model fitting folder
best_rd <- as.character(read.table(sprintf('<<<< FILEPATH REDACTED >>>>',root))[1,1])
message(sprintf('BEST RUN_DATE OF FITTED MODEL TO BE USED: %s', best_rd))

# references to some key directories needed here
mod_dir     <- sprintf('<<<< FILEPATH REDACTED >>>>', root, best_rd)
rd_dir      <- sprintf('<<<< FILEPATH REDACTED >>>>', root, run_date)
rd_nid_dir  <- sprintf('%s/%s', rd_dir, snid)
dat_dir     <- '<<<< FILEPATH REDACTED >>>>'

# make output directory for this nid within this run_date
suppressWarnings(dir.create(rd_nid_dir))

## load configs from the best run date
# ab times
ab_times <- fread(sprintf('<<<< FILEPATH REDACTED >>>>',mod_dir))

# config values into environment
load_config(mod_dir)

# hard code for now, should always be 10 or less
cores <- 10
## #######################################################




## #######################################################
## LOAD AND SUBSET DATA
message('\nLoading Global SBH datataset and subsetting... ')

data <- setDT(readRDS(sprintf('<<<< FILEPATH REDACTED >>>>',dat_dir)))


# subset to survey
d <- subset(data, nid == snid)
rm(data)

# turn subnational ids into country codes
d[, country := substr(country,1,3)]

# load in metadata for the countries to be able to get region
info <- fread(sprintf('<<<< FILEPATH REDACTED >>>>',root)) 
d    <- merge(d, info, by = 'country', all.x = TRUE)

# identify the region for this country, and thus the model run to pull from
region <- d$region[1]
MAP_region <- d$MAP_gbdregion[1] # Full name to match on MAP/POB distribution later on

# print out info about this survey and region
message(sprintf('RAW SURVEY INFO:\n  COUNTRY: %s \n  NID: %i \n  SURVEY: %s \n  YEAR: %i \n  NROWS: %s',
                d$country[1],snid,d$source[1],d$year[1], prettyNum(nrow(d),big.mark=',')))

# write raw dataset for this survey
saveRDS(d, sprintf('<<<< FILEPATH REDACTED >>>>',rd_nid_dir))

# load in and properly label sdi info, important covariate for later
sdi <-fread(sprintf('<<<< FILEPATH REDACTED >>>>',root)) [,c('location_id','year_id','sdi'),with=FALSE]
setnames(sdi,c('location_id','year_id'),c('loc_id','year_entering'))
sdi[, year_entering := year_entering - base_year]

# load in the centre scale info for this model run
cs <- fread(sprintf('<<<< FILEPATH REDACTED >>>>', mod_dir, region))

message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))

## #######################################################





## #######################################################
## CLEAN DATA AND MESSAGE OUT QUALITY INFORMATION
message('\nVetting Data Quality for this survey... ')


# check ceb and ced
message(sprintf('WARNING: There are %s rows (%f%%) missing CED.', 
        prettyNum(sum(is.na(d$ced)),big.mark = ','), sum(is.na(d$ced))/nrow(d)*100))
message(sprintf('WARNING: There are %s rows (%f%%) missing CEB.', 
        prettyNum(sum(is.na(d$ceb)),big.mark = ','), sum(is.na(d$ceb))/nrow(d)*100))
message(sprintf('WARNING: There are %s rows (%f%%) where ced>ceb.', 
                prettyNum(sum(d$ced>d$ceb,na.rm=TRUE),big.mark = ','), 
                sum(d$ced>d$ceb,na.rm=TRUE)/nrow(d)*100))

print(summary(d$ceb))
print(summary(d$ced))


# check for missing geographic identifiers, make geo unit identifier for later on collapse
d[, location_code := as.character(location_code)]
d$location_code[d$location_code=='NA']=''
d[, latnum        := as.numeric(latnum)]
d[, longnum       := as.numeric(longnum)]
d$point[is.na(d$point) & (d$location_code != "" & d$shapefile != "")] = 0
d$point[is.na(d$point) & (!is.na(d$latnum) & !is.na(d$longnum))] = 1
d[point==0, geo_unit := paste0(nid,location_code,shapefile)]
d[point==1, geo_unit := paste0(nid,cluster_number)]

message(sprintf('These data contain %s rows of point data and %s (%f%%) rows of polygon data.\n',
                prettyNum(sum(d$point==1),big.mark=','),
                prettyNum(sum(d$point==0),big.mark=','),
                sum(d$point==0)/nrow(d)*100))
d[, no_geog := (is.na(latnum) | is.na(longnum)) & (shapefile == '' | location_code == '')]
d$no_geog[is.na(d$no_geog)] <- TRUE
if(sum(d$no_geog) > 0)
  message(sprintf('WARNING: There are %s rows (%f%%) with NO GEOGRAPHIC IDENTIFIERS.\n',
                  prettyNum(sum(d$no_geog),big.mark=','),
                  sum(d$no_geog)/nrow(d)*100))

if(sum(d$no_geog)==nrow(d))
  stop('THIS SURVEY HAS NO GEOGRAPHIC IDENTIFIERS')
  
  
# check for missing maternal age
setnames(d, 'age_year', 'mothage_atint')
message(sprintf('WARNING: There are %s rows (%f%%) missing maternal age.', 
          prettyNum(sum(is.na(d$mothage_atint)),big.mark = ','), sum(is.na(d$mothage_atint))/nrow(d)*100))

print(summary(d$mothage_atint))



## if an entire survey is missing weights, then we give them a 1 weight.
#   This follows on the explicit choice to include surveys without weights as unweighted
#   The only non-weighted rows remaining after this will be polygon rows in surveys with
#   weights, theyll be dropped later on
if(sum(is.na(d$weight))/nrow(d)==1){
  message('WARNING: There are no survey weights in this survey and they will be set to 1.')
  d$weight <- 1
} else if(sum(is.na(d$weight))==0){
  message('All rows have survey weights')
} else {
  message(sprintf('WARNING: There are %s rows (%f%%) with missing survey weights that will be dropped.',
            prettyNum(sum(is.na(d$weight)),big.mark = ','), sum(is.na(d$weight))/nrow(d)*100))
}

# check for missing survey weight
## Point data will not be weighted, so we convert at point row weights to 1
d$weight[d$point == 1] <- 1

print(summary(d$weight))
d$svywgt <- d$weight

# check for missing interview date
if(all(is.na(d$interview_date_cmc))){
  message('All interview dates were missing so using the survey year for this survey')
  d[, yrint := as.double(year)]
} else{
  message('All interview dates were NOT missing so using the cmc variable for interview year.')
  d[, yrint := cmc_as_year(interview_date_cmc)]
}
if(any(is.na(d$yrint))){
  message('Some int dates still seem to be missing, so filling those in with the survey year')
  d$yrint[is.na(d$yrint)] <- as.double(d$year[1])
}


if(sum(is.na(d$yrint)) > 0)
  message(sprintf('WARNING: There are %s rows (%f%%) with a missing yrint variable.\n',
                  prettyNum(sum(is.na(d$yrint)),big.mark=','),
                  sum(is.na(d$yrint))/nrow(d)*100))
print(summary(d$yrint))

# scale yrint to base year to make comparable with other timining variables used here
d[, yrint := yrint - base_year]

# start a variable indicating rows to drop
if(agg_to_nid == TRUE){
  d[, must_drop := is.na(ceb) | is.na(ced) | ced > ceb                   | is.na(mothage_atint) | is.na(weight)]
} else {
  d[, must_drop := is.na(ceb) | is.na(ced) | ced > ceb | no_geog == TRUE | is.na(mothage_atint) | is.na(weight)]
}
if(sum(d$must_drop) > 0)
  message(sprintf('WARNING: Dropping %s rows (%f%%) from the final dataset due to missingess in key variables.',
                  prettyNum(sum(d$must_drop),big.mark=','),sum(d$must_drop)/nrow(d)*100))

# mothage needs to be int
d[, mothage_atint := round(mothage_atint, 0)]

d <- subset(d, must_drop == FALSE)


# keep only variables of interest
d <- d[, c('region','nid','source','country','loc_id','year','strata','svywgt',
           'point','latnum','longnum','location_code','shapefile',
           'geo_unit','yrint','ceb','ced','mothage_atint'),
       with = FALSE]


# drop any women who reported 0 ceb. Not a quality issue just need to be dropped
message(sprintf('Finally, dropping %s women (%f%%) who reported 0 CEB.',
                prettyNum(sum(d$ceb==0),big.mark=','),sum(d$ceb==0)/nrow(d)*100))
d <- subset(d, ceb > 0)



# make a mother_id variable
d$mother_id <- 1:nrow(d)


# Print new row count
message(sprintf('THERE ARE %s MOTHERS WITH COMPLETE INFORMATION IN THIS SURVEY
                 WITH %s CEB, AND A RAW CED/CEB RATIO OF %f',
                prettyNum(nrow(d),big.mark=','),
                prettyNum(sum(d$ceb),big.mark=','),
                sum(d$ced)/sum(d$ceb)))

message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))


## #######################################################





## #######################################################
## loading model object
message('\nLOAD MODEL OBJECT ... ')

# load in the fitted mgcv::bam model object for the model fit for this region
mod <- readRDS(sprintf('<<<< FILEPATH REDACTED >>>>', mod_dir, region))

message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))
## #######################################################






## #######################################################
## IF THE SURVEY IS LARGE THE REST OF THE SCRIPT MUST HAPPEN IN CHUNKS
## CONSERVATIVELY WE WILL SUBSET OUT CHUNKS OF ~100K OBSERVATIONS BUT WITHOUT CUTTING ACROSS GEO_UNIT TO MAINTAIN REPRESENTATIVENESS

# set mxrows to max nrow in a geo_unit or 10ok, which ever is higher
mxrows <- max( max(d[,.(n=.N),by=geo_unit]$n) , 1e5)
message(sprintf('Max number of rows for a subset set at %s rows.', prettyNum(mxrows,big.mark=',')))

# assign subsets to geo_units, get close to mxrows as possibles
subs <- d[, .(gu_nrows = .N), by = geo_unit]
subs[, subset := ceiling( cumsum(gu_nrows)/mxrows ) ]
d <- merge(d,subs,by='geo_unit')

# INFO PRINT
message(sprintf('There were %i subsets assigned, each with the following number of rows:',max(d$subset)))
print(table(d$subset))

# print a small msg explaining number of subsets expected
writeLines(as.character(max(d$subset)),sprintf('%s/expected_subsets_to_save.txt',rd_nid_dir))


#### #######################################################
## For aggregation purposes we will typically use geo_unit. 
## in some cases, we will can to just aggregate to the nid. 
## That is what is being set here. 
if(agg_to_nid == TRUE) {
  message('WARNING: agg_to_nid was set so only retrieving survey level estimates.')
  d[, geo_unit := nid]
}
## #######################################################

# save an original version of d
d_orig <- copy(d)

ndrawz <- 250

## #######################################################
# loop over subsets and complete the data setup (hyp childre) and prediction steps
for(sub in 1:max(d_orig$subset)){ 
  
  message(sprintf("\n\n\n\n~~~~~~~~~~ ON SUBSET %i of %i ~~~~~~~~~\n\n\n\n", sub, max(d_orig$subset)))
  message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))
  
  d <- subset(d_orig, subset == sub)
  
  
  ## #######################################################
  ## CLEAN DATA TO PREPARE FOR PREDICTION
  message('\nPREPARING DATA FOR PREDICTION ... ')
  
  # cedceb ratio
  d[, cdceb100 := ced/ceb * 100]
  
  # Expanding the frame to be by child-age
  nd  <- expand_prediction_frame(d, mothageatintvar = 'mothage_atint', abt = ab_times)
  message(sprintf('Prediction frame with %s rows (hyp child-age) created',prettyNum(nrow(nd),big.mark=',')))
  rm(d) # clear from mem
  message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))
  
  
  # add SDI
  message(' ... ... Adding year-varying covariates to child-age prediction frame')
  
  # calculate the year entering each age bin
  nd[, year_entering := floor(yrborn + tstart/12)]
  
  # get sdi by yrborn
  nd <- merge(nd,sdi,by=c('year_entering','loc_id'),all.x=TRUE)
  nd <- subset(nd,year_entering < (2017-base_year)) # limit of SDI data
  # Backfill SDI
  sdib <- nd[,.(sdib = min(sdi,na.rm=TRUE)),by=loc_id]  
  nd <- merge(nd,sdib,by='loc_id')
  nd$sdi[is.na(nd$sdi)] <- nd$sdib[is.na(nd$sdi)]
  
  
  # Re-scale variables that were in unit scale in model fitting
  message('Center scaling variables in new prediction frame')
  for(v in cs$name[cs$mean!=0]) {
    message(sprintf('center scaling %s',v))
    nd[[paste0(v,'_orig')]] <- nd[[v]]
    nd[[v]] <- centscale(nd[[v]],v)
  }
  
  
  # load in MAP/POB distribution for this region
  message('Imputing pct_ceb and birthorder variables for each hypothetical child')
  load(sprintf('<<<< FILEPATH REDACTED >>>>',root))
  map <- subset(setDT(distributions$ceb),gbdregion == MAP_region)
  map[,agegroup:=as.numeric(substr(agegroup, 1, 2))]
  
  
  # get birthorder imputed from the map
  nd <- get_pct_ceb(newdata = nd, mapp = map, mothageatintvar = 'mothage_atint', cebvar = 'ceb')
  message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))
  
  
  # centrescale pct_ceb and birthorder if needed
  nd$pct_ceb <- centscale(var=nd$pct_ceb_orig,varname='pct_ceb')
  nd$birthorder_orig <- nd$pct_ceb_orig*nd$ceb # make sure ceb hasnt been center scaled, shouldnt be if wasnt used in model
  nd$birthorder <- centscale(var=nd$birthorder_orig,varname='birthorder')
  
  
  # make sure that we didnt have a bug in center scaling and that ranges of model data and prediction data are similar
  message('Just a check to see that model data and prediction frame vars in the same range, what with all the center scaling')
  vars <- get_fe_from_gamformula(formula)
  for(v in vars){
    message(paste0('\n',v))
    message('Model Frame:'); print(summary(mod$model[[v]]))
    message('Prediction Frame:'); print(summary(nd[[v]]))
  }
  
  # make sure we also have country ab
  nd[, country_ab := as.factor(paste0(country,"_",ab))]
  
  message(' Completed setting up prediction frame, proceeding to predicting hazards')
  
  ##################################################################
  
  
  
  
  ##################################################################
  ## PREDICTION, get draw hazard estimates for each hyp-child, age bin in the data
  
  # save nd info by geo_unit to append later on
  if(agg_to_nid){
    geo_unit_info <- unique(nd[,c('nid','year','source','loc_id','geo_unit'),with=FALSE])  
  } else {
    geo_unit_info <- unique(nd[,c('nid','year','source','loc_id','point',
                                  'longnum','latnum','location_code','shapefile','geo_unit'),with=FALSE])
  }
  setnames(geo_unit_info, 'year', 'svyyr')
  
  
  # Making individual level predictions on newdata
  message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))
  
  res <- predict_Dsurv(model     = mod,
                       new_data  = nd,
                       yrbornvar = 'yrborn_orig',
                       ndraws    = ndrawz, 
                       summarize = FALSE,
                       ncores    = cores)
  rm(nd) # clear from mem
  
  message(' Prediction is complete, aggregating period estimates')
  message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))
  
  # agregate over period using MAP to get exp number entering each bin
  aggresl <- aggregate_period_estimates(est            = res,
                                        map             = map,
                                        use_svy_weight  = TRUE,
                                        yrbornvar       = 'yrborn_orig',
                                        cebvar          = 'ceb',
                                        mothageatintvar = 'mothage_atint',
                                        yrs_to_agg      = 1)
  rm(res) # clear from mem
  message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))
  
  # sep out summary and draw level estimates
  draws_agg <- aggresl$draws
  aggres    <- aggresl$res
  rm(aggresl)
  
  # set up in normal year space
  aggres[,     year_entering := year_entering + base_year - 0.5]
  draws_agg [, year_entering := year_entering + base_year - 0.5] 
  
  # merge back on the geo unit info
  aggres     <- merge(aggres,    geo_unit_info, by = 'geo_unit', all.x = TRUE)
  draws_agg  <- merge(draws_agg, geo_unit_info, by = 'geo_unit', all.x = TRUE)
  
  # SAVE FINAL AGGREGATED YEARLY SUMMARY DATA
  saveRDS(aggres, sprintf('<<<< FILEPATH REDACTED >>>>',rd_nid_dir,sub))
  
  # SAVE FINAL AGGREGATED YEARLY DRAW-LEVEL DATA
  saveRDS(draws_agg, sprintf('<<<< FILEPATH REDACTED >>>>',rd_nid_dir,sub))
  rm(draws_agg) # clear from mem
  
  # plot out the results just for fun
  pdf(sprintf('<<<< FILEPATH REDACTED >>>>',rd_nid_dir,sub))
  for(g in unique(aggres$geo_unit)){
    g <- ggplot(subset(aggres, geo_unit == g), 
                aes(x = year_entering, y = haz, group = ab, color = ab, fill = ab)) + 
                theme_bw() + ggtitle(sprintf('geo_unit: %s',g)) +
                geom_point(aes(size=EEB),alpha=0.5)+
                geom_ribbon(aes(ymax = haz_upper, ymin = haz_lower), alpha = 0.5)
    
    plot(g)
  }
  dev.off()
  rm(aggres) # clear from mem
  
}

# 
message(sprintf('\nCURRENT MEMORY USAGE: %.3f GB\n', round(mem_used()/1e9,3)))

message('DONE')
  
##################################################################
  
# write a note to detect later that this run has finished
writeLines('WOOHOO', sprintf('%s/fin.txt',rd_nid_dir))
  
  
